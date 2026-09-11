#!/usr/bin/env python3
"""
Adult-weight predictor from breed composition, validated by cross-validation on
customer-reported weights.

Why this exists
---------------
The stage-11 weight (dense LMM breed-level prediction blended with the Darwin's
Ark size PRS) was calibrated on 94 cohort dogs and claimed r=0.92 / MAE 5.1 kg.
On 878 held-out ProsperK9 customers (2026-09-10) it gave r=0.63, MAE 8.2 kg,
bias +4.1 kg and slope(pred~actual)=0.57: small dogs were over-predicted by up
to +114%, giants under-predicted by ~30%.

What this script does
---------------------
1. Loads customer weights (kit -> lb), ages (from the microbiome CSV,
   sample_id = DW<kit>), and each kit's breed_result.json + prs_result.json.
2. Drops dogs under MIN_AGE years (puppies are weighed at swab time, so their
   label is not an adult weight).
3. Compares, by repeated 5-fold CV (never in-sample):
     M0  current pipeline pred_kg
     M1  linear recalibration       y ~ a + b*pred
     M1L log-linear recalibration   log y ~ a + b*log pred
     M2  NNLS breed model           y ~ sum_b w_b * p_b, w >= 0, ridge toward a
                                    prior (lambda chosen by nested CV)
     M2L same in log space          log y ~ sum_b w_b * p_b
     M3  M2 stacked with genomic covariates (prs_z, DA size PRS, dense pred,
         height) via inner out-of-fold breed predictions + OLS
   and reports r, MAE, bias overall and per actual-weight bin.
4. Refits the chosen model on all adults and writes
   reference_json/weight_breed_model.json for stage 11.

Usage
-----
  fit_weight_model.py <kits_dir> <weights.tsv> <age_weight.tsv> <run_dog_pipeline.sh> <out_dir>

kits_dir holds pk-<kit>/{breed_result.json,prs_result.json}; weights.tsv is
kit<TAB>weight_lb; age_weight.tsv is sample_id<TAB>age<TAB>weight (age in years).
"""
import csv
import json
import os
import re
import sys

import numpy as np
from scipy.optimize import nnls, lsq_linear

LB_PER_KG = 2.20462
MIN_AGE = 1.0
BINS = [(0, 8), (8, 15), (15, 25), (25, 40), (40, 999)]
LAM_GRID = [0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
N_REPEATS, N_FOLDS = 5, 5
BOX = 1.5        # bounded model: each breed's weight stays within [prior/BOX, prior*BOX]


# ── helpers ──────────────────────────────────────────────────────────────────
def extract_dict(src, name):
    m = re.search(r'^%s\s*=\s*\{' % re.escape(name), src, re.M)
    i, depth = m.start(), 0
    while i < len(src):
        if src[i] == '{':
            depth += 1
        elif src[i] == '}':
            depth -= 1
            if depth == 0:
                break
        i += 1
    ns = {}
    exec(src[m.start():i + 1], {}, ns)
    return ns[name]


def tokens(name):
    out = set()
    for w in re.sub(r'[^a-z ]', ' ', name.lower()).split():
        w = re.sub(r'ies$', 'y', w).rstrip('s')
        if w and w not in ('dog', 'the'):
            out.add(w)
    return frozenset(out)


def join_prior(breeds, table):
    """Panel label -> table value. Exact token-set match first, then a unique
    subset/superset match (GERMAN_SHEPHERD ~ 'German Shepherd Dogs')."""
    tt = {tokens(k): v for k, v in table.items()}
    out = {}
    for b in breeds:
        t = tokens(b.replace('_', ' '))
        if t in tt:
            out[b] = tt[t]
            continue
        cands = [v for k, v in tt.items() if (t <= k or k <= t) and len(t & k) >= 2]
        if len(cands) == 1:
            out[b] = cands[0]
    return out


def metrics(y, p):
    y, p = np.asarray(y, float), np.asarray(p, float)
    out = {'n': int(len(y)),
           'r': float(np.corrcoef(y, p)[0, 1]) if len(y) > 2 else float('nan'),
           'mae': float(np.mean(np.abs(p - y))),
           'bias': float(np.mean(p - y)),
           'slope': float(np.polyfit(y, p, 1)[0]) if len(y) > 2 else float('nan')}
    for lo, hi in BINS:
        m = (y >= lo) & (y < hi)
        if m.sum() >= 3:
            out[f'bin_{lo}_{hi}'] = {'n': int(m.sum()),
                                     'mae': float(np.mean(np.abs(p[m] - y[m]))),
                                     'bias': float(np.mean(p[m] - y[m])),
                                     'bias_pct': float(100 * np.mean(p[m] - y[m]) / np.mean(y[m]))}
    return out


def fmt(name, m):
    s = f"{name:<34} n={m['n']:<4} r={m['r']:.3f} MAE={m['mae']:5.2f} bias={m['bias']:+5.2f} slope={m['slope']:.2f} |"
    for lo, hi in BINS:
        b = m.get(f'bin_{lo}_{hi}')
        s += f" {b['bias']:+5.1f}({b['bias_pct']:+4.0f}%)" if b else '    n/a     '
    return s


def folds(n, seed):
    rng = np.random.RandomState(seed)
    idx = rng.permutation(n)
    return [idx[k::N_FOLDS] for k in range(N_FOLDS)]


# ── models ───────────────────────────────────────────────────────────────────
def fit_nnls_ridge(X, y, prior, lam):
    """min ||y - Xw||^2 + lam ||w - prior||^2, w >= 0 (augmented NNLS)."""
    if lam > 0:
        s = np.sqrt(lam)
        Xa = np.vstack([X, s * np.eye(X.shape[1])])
        ya = np.concatenate([y, s * prior])
    else:
        Xa, ya = X, y
    w, _ = nnls(Xa, ya, maxiter=50 * X.shape[1])
    return w


def fit_box_ridge(X, ly, lprior, lam, box=BOX):
    """Log-space, box-constrained ridge toward the prior:
    min ||ly - Xw||^2 + lam ||w - lprior||^2,  lprior - log(box) <= w <= lprior + log(box).
    The box keeps every breed's implied purebred weight within a factor of
    `box` of its prior, so a toy breed can't be driven to ~0 kg by mixes."""
    s_ = np.sqrt(max(lam, 1e-6))
    Xa = np.vstack([X, s_ * np.eye(X.shape[1])])
    ya = np.concatenate([ly, s_ * lprior])
    lb, ub = lprior - np.log(box), lprior + np.log(box)
    return lsq_linear(Xa, ya, bounds=(lb, ub), lsmr_tol='auto', max_iter=2000).x


def box_cv_lambda(X, ly, lprior, seed):
    best, best_mae = None, np.inf
    fl = folds(len(ly), seed + 1000)
    for lam in LAM_GRID:
        pred = np.zeros(len(ly))
        for te in fl:
            tr = np.setdiff1d(np.arange(len(ly)), te)
            pred[te] = np.exp(X[te] @ fit_box_ridge(X[tr], ly[tr], lprior, lam))
        mae = np.mean(np.abs(pred - np.exp(ly)))
        if mae < best_mae:
            best, best_mae = lam, mae
    return best


def box_inner_oof(X, ly, lprior, lam, seed):
    pred = np.zeros(len(ly))
    for te in folds(len(ly), seed + 2000):
        tr = np.setdiff1d(np.arange(len(ly)), te)
        pred[te] = np.exp(X[te] @ fit_box_ridge(X[tr], ly[tr], lprior, lam))
    return pred


def marginal_prior(X, y, prior_table, min_support=2.0):
    """Per-breed prior: AKC table where it joins; otherwise the
    proportion-weighted mean weight of the dogs carrying the breed (shrunk to
    the global mean when support < min_support dogs)."""
    supp = X.sum(0)
    marg = (X * y[:, None]).sum(0) / np.maximum(supp, 1e-9)
    gm = float(y.mean())
    marg = (supp * marg + min_support * gm) / (supp + min_support)
    return np.where(np.isfinite(prior_table), prior_table, marg)


def nnls_cv_lambda(X, y, prior, seed):
    """Pick lambda by inner 5-fold CV (MAE)."""
    best, best_mae = None, np.inf
    fl = folds(len(y), seed + 1000)
    for lam in LAM_GRID:
        pred = np.zeros(len(y))
        for te in fl:
            tr = np.setdiff1d(np.arange(len(y)), te)
            w = fit_nnls_ridge(X[tr], y[tr], prior, lam)
            pred[te] = X[te] @ w
        mae = np.mean(np.abs(pred - y))
        if mae < best_mae:
            best, best_mae = lam, mae
    return best


def inner_oof(X, y, prior, lam, seed):
    pred = np.zeros(len(y))
    for te in folds(len(y), seed + 2000):
        tr = np.setdiff1d(np.arange(len(y)), te)
        pred[te] = X[te] @ fit_nnls_ridge(X[tr], y[tr], prior, lam)
    return pred


def ols(A, y):
    A1 = np.column_stack([np.ones(len(y)), A])
    b, *_ = np.linalg.lstsq(A1, y, rcond=None)
    return b


HUBER_C = 6.0   # kg; residuals beyond this get down-weighted (label noise / obesity)
LOG_HUBER_C = 0.3   # same idea in log space (~35%)


def huber_ols(A, y, c=HUBER_C, iters=20):
    A1 = np.column_stack([np.ones(len(y)), A])
    wts = np.ones(len(y))
    for _ in range(iters):
        b, *_ = np.linalg.lstsq(A1 * np.sqrt(wts)[:, None], y * np.sqrt(wts), rcond=None)
        r = np.abs(y - A1 @ b)
        wts = np.minimum(1.0, c / np.maximum(r, 1e-9))
    return b


def fit_nnls_huber(X, y, prior, lam, c=HUBER_C, iters=10):
    wts = np.ones(len(y))
    for _ in range(iters):
        sw = np.sqrt(wts)[:, None]
        w = fit_nnls_ridge(X * sw, y * sw[:, 0], prior, lam)
        r = np.abs(y - X @ w)
        wts = np.minimum(1.0, c / np.maximum(r, 1e-9))
    return w


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    kits_dir, weights_tsv, age_tsv, pipeline_sh, out_dir = sys.argv[1:6]
    os.makedirs(out_dir, exist_ok=True)

    weight_lb = {}
    for line in open(weights_tsv):
        k, v = line.split('\t')[:2]
        try:
            weight_lb[k] = float(v)
        except ValueError:
            pass
    age = {}
    for r in csv.DictReader(open(age_tsv), delimiter='\t'):
        try:
            age[r['sample_id'].replace('DW', '', 1)] = float(r['age'])
        except ValueError:
            pass

    rows = []
    for d in sorted(os.listdir(kits_dir)):
        kit = d[3:] if d.startswith('pk-') else d
        if kit not in weight_lb:
            continue
        try:
            br = json.load(open(os.path.join(kits_dir, d, 'breed_result.json')))
            pr = json.load(open(os.path.join(kits_dir, d, 'prs_result.json')))
        except (OSError, ValueError):
            continue
        wk = pr.get('physical_traits', {}).get('weight_kg')
        if not wk:
            continue
        comp = {r['breed']: float(r['proportion'])
                for r in (br.get('breed_composition_raw') or br.get('breed_composition') or [])}
        rows.append({'kit': kit, 'y': weight_lb[kit] / LB_PER_KG, 'age': age.get(kit, np.nan),
                     'comp': comp,
                     'pred': float(wk['pred_kg']),
                     'dense': float(wk.get('pred_kg_dense', wk['pred_kg'])),
                     'prs_z': float(wk.get('prs_z', 0.0)),
                     'da': float(wk.get('da_size_prs', 0.0)),
                     'height': float(pr['physical_traits'].get('height_cm', {}).get('pred_cm', np.nan))})
    print(f'{len(rows)} dogs with weight + breed_result + prs_result')

    breeds = sorted({b for r in rows for b in r['comp']})
    bi = {b: i for i, b in enumerate(breeds)}
    print(f'{len(breeds)} breed labels')

    # ── baseline on everyone, then on adults ────────────────────────────────
    y_all = np.array([r['y'] for r in rows]); p_all = np.array([r['pred'] for r in rows])
    print('\nCurrent pipeline pred_kg:')
    print(fmt('  all ages', metrics(y_all, p_all)))
    ages = np.array([r['age'] for r in rows])
    keep = ages >= MIN_AGE
    print(f'  {(~keep).sum()} dogs under {MIN_AGE} y dropped (puppies weighed at swab time)')
    rows = [r for r, k in zip(rows, keep) if k]

    n = len(rows)
    y = np.array([r['y'] for r in rows])
    X = np.zeros((n, len(breeds)))
    for i, r in enumerate(rows):
        for b, p in r['comp'].items():
            X[i, bi[b]] = p
    # NNLS on proportions: renormalise rows so a dog whose components were
    # trimmed (<5e-5 dropped in stage 9) still sums to 1.
    X = X / np.maximum(X.sum(1, keepdims=True), 1e-9)
    pred0 = np.array([r['pred'] for r in rows])
    cov_all = np.column_stack([[r['prs_z'] for r in rows], [r['da'] for r in rows],
                               [r['dense'] for r in rows], [r['height'] for r in rows]])
    cov_all = np.where(np.isnan(cov_all), np.nanmean(cov_all, 0), cov_all)
    COV_NAMES = ['prs_z', 'da_size_prs', 'dense_kg', 'height_cm']
    print(fmt(f'  adults >= {MIN_AGE} y', metrics(y, pred0)))

    # AKC breed-standard prior where the name joins, else global mean
    src = open(pipeline_sh).read()
    bw = extract_dict(src, 'BREED_WEIGHT_KG')
    joined = join_prior(breeds, bw)
    prior_akc = np.full(len(breeds), np.nan)
    for b, v in joined.items():
        prior_akc[bi[b]] = v
    n_joined = int(np.isfinite(prior_akc).sum())
    print(f'{n_joined}/{len(breeds)} breed labels join the AKC weight table (prior); rest -> global mean')

    # ── repeated k-fold CV ──────────────────────────────────────────────────
    oof = {k: np.zeros((N_REPEATS, n)) for k in
           ['M0', 'P0', 'M1', 'M1Q', 'M1H', 'M1L', 'M2_mean', 'M2_akc', 'M2H_akc', 'M2L_akc',
            'M3_prs_z', 'M3_da', 'M3_dense', 'M3_height', 'M3_all', 'M3_blend', 'M3H_blend',
            'M2B', 'M3B_blend', 'M3BH_blend', 'M3BH_da', 'M3BL_da', 'M3BL_blend', 'M3BL_both']}
    lam_chosen = {'M2_mean': [], 'M2_akc': [], 'M2L_akc': [], 'M2B': []}
    for rep in range(N_REPEATS):
        for te in folds(n, rep):
            tr = np.setdiff1d(np.arange(n), te)
            ytr, Xtr, Xte = y[tr], X[tr], X[te]
            oof['M0'][rep, te] = pred0[te]
            b = ols(pred0[tr, None], ytr)
            oof['M1'][rep, te] = b[0] + b[1] * pred0[te]
            b = ols(np.column_stack([pred0[tr], pred0[tr] ** 2]), ytr)
            oof['M1Q'][rep, te] = b[0] + b[1] * pred0[te] + b[2] * pred0[te] ** 2
            b = huber_ols(pred0[tr, None], ytr)
            oof['M1H'][rep, te] = b[0] + b[1] * pred0[te]
            b = ols(np.log(pred0[tr, None]), np.log(ytr))
            oof['M1L'][rep, te] = np.exp(b[0] + b[1] * np.log(pred0[te]))

            gm = float(ytr.mean())
            pr_mean = np.full(len(breeds), gm)
            pr_akc = np.where(np.isfinite(prior_akc), prior_akc, gm)
            lam = nnls_cv_lambda(Xtr, ytr, pr_mean, rep); lam_chosen['M2_mean'].append(lam)
            oof['M2_mean'][rep, te] = Xte @ fit_nnls_ridge(Xtr, ytr, pr_mean, lam)
            lam_a = nnls_cv_lambda(Xtr, ytr, pr_akc, rep); lam_chosen['M2_akc'].append(lam_a)
            w_akc = fit_nnls_ridge(Xtr, ytr, pr_akc, lam_a)
            oof['M2_akc'][rep, te] = Xte @ w_akc
            oof['P0'][rep, te] = Xte @ pr_akc
            w_h = fit_nnls_huber(Xtr, ytr, pr_akc, lam_a)
            oof['M2H_akc'][rep, te] = Xte @ w_h
            ly = np.log(ytr)
            pr_l = np.log(np.where(np.isfinite(prior_akc), prior_akc, np.exp(ly.mean())))
            lam_l = nnls_cv_lambda(Xtr, ly, pr_l, rep); lam_chosen['M2L_akc'].append(lam_l)
            oof['M2L_akc'][rep, te] = np.exp(Xte @ fit_nnls_ridge(Xtr, ly, pr_l, lam_l))

            # stacking: inner OOF breed pred + covariates -> OLS
            inner = inner_oof(Xtr, ytr, pr_akc, lam_a, rep)
            outer = Xte @ w_akc
            for name, cols in [('M3_prs_z', [0]), ('M3_da', [1]), ('M3_dense', [2]),
                               ('M3_height', [3]), ('M3_all', [0, 1, 2, 3])]:
                A_tr = np.column_stack([inner, cov_all[tr][:, cols]])
                A_te = np.column_stack([outer, cov_all[te][:, cols]])
                b = ols(A_tr, ytr)
                oof[name][rep, te] = b[0] + A_te @ b[1:]
            A_tr = np.column_stack([inner, pred0[tr]]); A_te = np.column_stack([outer, pred0[te]])
            b = ols(A_tr, ytr); oof['M3_blend'][rep, te] = b[0] + A_te @ b[1:]
            b = huber_ols(A_tr, ytr); oof['M3H_blend'][rep, te] = b[0] + A_te @ b[1:]

            # bounded log-space breed model (marginal prior for un-joined breeds)
            lpr = np.log(marginal_prior(Xtr, ytr, prior_akc))
            lam_b = box_cv_lambda(Xtr, ly, lpr, rep); lam_chosen['M2B'].append(lam_b)
            w_b = fit_box_ridge(Xtr, ly, lpr, lam_b)
            outer_b = np.exp(Xte @ w_b)
            oof['M2B'][rep, te] = outer_b
            inner_b = box_inner_oof(Xtr, ly, lpr, lam_b, rep)
            A_tr = np.column_stack([inner_b, pred0[tr]]); A_te = np.column_stack([outer_b, pred0[te]])
            b = ols(A_tr, ytr); oof['M3B_blend'][rep, te] = b[0] + A_te @ b[1:]
            b = huber_ols(A_tr, ytr); oof['M3BH_blend'][rep, te] = b[0] + A_te @ b[1:]
            A_tr = np.column_stack([inner_b, cov_all[tr][:, 1]]); A_te = np.column_stack([outer_b, cov_all[te][:, 1]])
            b = huber_ols(A_tr, ytr); oof['M3BH_da'][rep, te] = b[0] + A_te @ b[1:]
            # log-space stacks: multiplicative PRS effect, so a toy dog can't be dragged to ~0 kg
            for name, cols_tr, cols_te in [
                    ('M3BL_da', [np.log(inner_b), cov_all[tr][:, 1]], [np.log(outer_b), cov_all[te][:, 1]]),
                    ('M3BL_blend', [np.log(inner_b), np.log(pred0[tr])], [np.log(outer_b), np.log(pred0[te])]),
                    ('M3BL_both', [np.log(inner_b), np.log(pred0[tr]), cov_all[tr][:, 1]],
                                  [np.log(outer_b), np.log(pred0[te]), cov_all[te][:, 1]])]:
                b = huber_ols(np.column_stack(cols_tr), ly, c=LOG_HUBER_C)
                oof[name][rep, te] = np.exp(b[0] + np.column_stack(cols_te) @ b[1:])

    print(f'\n{N_REPEATS}x{N_FOLDS}-fold CV on {n} adults (mean over repeats; per-bin bias kg (%) for actual '
          + ', '.join(f'{lo}-{hi}' for lo, hi in BINS) + ' kg):')
    summary = {}
    for k, P in oof.items():
        ms = [metrics(y, P[r]) for r in range(N_REPEATS)]
        agg = {'n': n}
        for key in ['r', 'mae', 'bias', 'slope']:
            agg[key] = float(np.mean([m[key] for m in ms]))
        for lo, hi in BINS:
            bk = f'bin_{lo}_{hi}'
            if bk in ms[0]:
                agg[bk] = {kk: float(np.mean([m[bk][kk] for m in ms])) for kk in ms[0][bk]}
        summary[k] = agg
        print(fmt(k, agg))
    for k, v in lam_chosen.items():
        print(f'  {k}: lambda chosen per fold {sorted(set(v))} (median {np.median(v)})')

    print('\nCalibration: bias kg (%) binned by PREDICTED weight (a calibrated model is ~0 in every bin):')
    for k, P in oof.items():
        line = f'{k:<12}'
        for lo, hi in BINS:
            m = (P >= lo) & (P < hi)
            if m.sum() >= 15:
                d = (P - y[None, :])[m]
                line += f' {np.mean(d):+5.1f}({100 * np.mean(d) / np.mean(np.broadcast_to(y, P.shape)[m]):+4.0f}%) n={int(m.sum() / N_REPEATS):<4}'
            else:
                line += '        n/a        '
        summary[k]['calib_by_pred'] = line
        print(line)

    # paired bootstrap (dogs) of MAE differences, using per-dog mean OOF pred
    rng = np.random.RandomState(7)
    bs = rng.randint(0, n, size=(2000, n))
    mean_pred = {k: P.mean(0) for k, P in oof.items()}
    ae = {k: np.abs(mean_pred[k] - y) for k in oof}
    print('\nMAE difference vs M0 and vs M1 (kg, paired bootstrap 95% CI over dogs):')
    for k in oof:
        for ref in ['M0', 'M1']:
            if k == ref:
                continue
            d = (ae[k] - ae[ref])[bs].mean(1)
            summary[k][f'dmae_vs_{ref}'] = [float(np.mean(d)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]
        s0, s1 = summary[k].get('dmae_vs_M0'), summary[k].get('dmae_vs_M1')
        print(f"  {k:<12}" + (f" vs M0 {s0[0]:+5.2f} [{s0[1]:+5.2f},{s0[2]:+5.2f}]" if s0 else ' ' * 30)
              + (f"   vs M1 {s1[0]:+5.2f} [{s1[1]:+5.2f},{s1[2]:+5.2f}]" if s1 else ''))

    # ── choose winner (lowest CV MAE; tie -> simpler) ───────────────────────
    order = ['M1', 'M1H', 'M1Q', 'M1L', 'M2_mean', 'M2_akc', 'M2H_akc', 'M2L_akc', 'M3_prs_z', 'M3_da', 'M3_dense', 'M3_height', 'M3_all', 'M3_blend', 'M3H_blend', 'M2B', 'M3B_blend', 'M3BH_blend', 'M3BH_da', 'M3BL_da', 'M3BL_blend', 'M3BL_both']
    winner = min(order, key=lambda k: (round(summary[k]['mae'], 2), order.index(k)))
    print(f'\nLowest CV MAE: {winner}')
    if len(sys.argv) > 6:
        winner = sys.argv[6]
        print(f'Shipping (chosen on the purebred sanity table as well as MAE): {winner}')

    # ── refit winner on all adults, write model ─────────────────────────────
    gm = float(y.mean())
    pr_akc = np.where(np.isfinite(prior_akc), prior_akc, gm)
    lam_a = float(np.median(lam_chosen['M2_akc']))
    model = {'model': winner, 'n_train': n, 'min_age_years': MIN_AGE,
             'truth': 'ProsperK9 customer-reported weights (lb), adults only',
             'cv': {'design': f'{N_REPEATS}x{N_FOLDS}-fold, lambda by nested CV',
                    'summary': summary, 'winner': winner},
             'breeds': breeds}
    if winner in ('M2B', 'M3B_blend', 'M3BH_blend', 'M3BH_da', 'M3BL_da', 'M3BL_blend', 'M3BL_both'):
        ly = np.log(y)
        lpr = np.log(marginal_prior(X, y, prior_akc))
        lam_b = float(np.median(lam_chosen['M2B']))
        w_b = fit_box_ridge(X, ly, lpr, lam_b)
        model.update({'space': 'log', 'lambda': lam_b, 'box': BOX,
                      'prior_source': f'AKC breed-standard table for {n_joined} labels, proportion-weighted marginal mean for the rest',
                      'breed_weight': {b: round(float(np.exp(w_b[i])), 3) for i, b in enumerate(breeds)},
                      'breed_prior': {b: round(float(np.exp(lpr[i])), 3) for i, b in enumerate(breeds)},
                      'default_kg': round(float(np.exp(ly.mean())), 3)})
        if winner.startswith('M3BL'):
            inner_b = box_inner_oof(X, ly, lpr, lam_b, 0)
            cols = {'M3BL_da': [cov_all[:, 1]], 'M3BL_blend': [np.log(pred0)],
                    'M3BL_both': [np.log(pred0), cov_all[:, 1]]}[winner]
            names = {'M3BL_da': ['da_size_prs'], 'M3BL_blend': ['log_pred_kg_current'],
                     'M3BL_both': ['log_pred_kg_current', 'da_size_prs']}[winner]
            b = huber_ols(np.column_stack([np.log(inner_b)] + cols), ly, c=LOG_HUBER_C)
            model['stack'] = {'space': 'log', 'intercept': float(b[0]), 'coef_log_breed_pred': float(b[1]),
                              'covariates': {nm: float(b[2 + j]) for j, nm in enumerate(names)},
                              'robust': f'Huber c={LOG_HUBER_C} in log space',
                              'formula': 'pred_kg = exp(intercept + coef_log_breed_pred*log(breed_pred_kg) + sum(coef*cov))'}
        elif winner != 'M2B':
            inner_b = box_inner_oof(X, ly, lpr, lam_b, 0)
            col = pred0 if winner != 'M3BH_da' else cov_all[:, 1]
            fitf = huber_ols if winner.startswith('M3BH') else ols
            b = fitf(np.column_stack([inner_b, col]), y)
            model['stack'] = {'intercept': float(b[0]), 'coef_breed_pred': float(b[1]),
                              'covariate': 'pred_kg_current' if winner != 'M3BH_da' else 'da_size_prs',
                              'coef_covariate': float(b[2]), 'robust': winner.startswith('M3BH')}
    elif winner.startswith('M2L'):
        ly = np.log(y)
        pr_l = np.log(np.where(np.isfinite(prior_akc), prior_akc, np.exp(ly.mean())))
        lam_l = float(np.median(lam_chosen['M2L_akc']))
        w = fit_nnls_ridge(X, ly, pr_l, lam_l)
        model.update({'space': 'log', 'lambda': lam_l, 'prior_default_log_kg': float(ly.mean()),
                      'breed_weight': {b: round(float(np.exp(w[i])), 3) for i, b in enumerate(breeds)}})
    elif winner.startswith('M2') or winner.startswith('M3'):
        w = fit_nnls_ridge(X, y, pr_akc, lam_a)
        model.update({'space': 'linear', 'lambda': lam_a, 'prior_default_kg': gm,
                      'breed_weight': {b: round(float(w[i]), 3) for i, b in enumerate(breeds)}})
        if winner.startswith('M3'):
            cols = {'M3_prs_z': [0], 'M3_da': [1], 'M3_dense': [2], 'M3_height': [3],
                    'M3_all': [0, 1, 2, 3]}[winner]
            inner = inner_oof(X, y, pr_akc, lam_a, 0)
            b = ols(np.column_stack([inner, cov_all[:, cols]]), y)
            model['stack'] = {'intercept': float(b[0]), 'coef_breed_pred': float(b[1]),
                              'covariates': {COV_NAMES[c]: float(b[2 + j]) for j, c in enumerate(cols)}}
    else:
        if winner == 'M1':
            b = ols(pred0[:, None], y)
            model.update({'space': 'linear', 'recal': {'intercept': float(b[0]), 'slope': float(b[1])}})
        else:
            b = ols(np.log(pred0[:, None]), np.log(y))
            model.update({'space': 'log', 'recal': {'intercept': float(b[0]), 'slope': float(b[1])}})

    out = os.path.join(out_dir, 'weight_breed_model.json')
    json.dump(model, open(out, 'w'), indent=1)
    print(f'wrote {out}')
    with open(os.path.join(out_dir, 'cv_predictions.tsv'), 'w') as f:
        f.write('kit\tage\tactual_kg\tM0_kg\t' + winner + '_oof_kg\ttop_breed\ttop_prop\n')
        for i, r in enumerate(rows):
            tb = max(r['comp'].items(), key=lambda kv: kv[1]) if r['comp'] else ('', 0)
            f.write(f"{r['kit']}\t{r['age']:.2f}\t{y[i]:.2f}\t{pred0[i]:.2f}\t{oof[winner].mean(0)[i]:.2f}\t{tb[0]}\t{tb[1]:.3f}\n")
    if 'breed_weight' in model:
        bwk = model['breed_weight']
        supp = X.sum(0)
        print('\nFitted per-breed adult weights (kg) = what a 100% purebred would be predicted, 30 best-supported breeds:')
        for i in np.argsort(-supp)[:30]:
            b = breeds[i]
            pa = prior_akc[i]
            print(f"  {b:<34} {bwk[b]:6.1f}  support={supp[i]:6.1f} dogs  akc_prior={'%.1f' % pa if np.isfinite(pa) else '   - '}"
                  + (f"  prior_used={model['breed_prior'][b]:.1f}" if 'breed_prior' in model else ''))


if __name__ == '__main__':
    main()
