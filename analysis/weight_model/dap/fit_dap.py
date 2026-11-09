#!/usr/bin/env python3
"""
Refit the stage-11 adult-weight model with the Dog Aging Project (DAP) as extra
training data, and test every candidate on ProsperK9 customers, the population
the model serves.

    fit_dap.py BUILD_DIR run_dog_pipeline.sh CURRENT_MODEL_JSON OUT_DIR

BUILD_DIR comes from build_dap.py (dap.npz, dap_pheno.tsv, customers.json).
Same model family as the shipped M3BL_da_sex (fit_weight_model.py): per-breed
log weights by box-constrained ridge toward the AKC prior, stacked in log space
with the Darwin's Ark size PRS and sex (Huber). Candidates:

  NOW    shipped weight_breed_model.json applied to today's (v2.3) compositions
  CUST   customer-only refit on v2.3 labels, 5x5-fold CV
  DAP    fit on DAP adults only, applied to every customer (fully external)
  DAPR   DAP breed weights, stack (calibration) refit on customer training folds
  POOL   breed weights on DAP + customer training folds, stack on customers
Plus DAP-internal 5-fold CV of the covariates customers don't have (neuter, age).
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from fit_weight_model import (extract_dict, join_prior, metrics, fmt, folds, fit_box_ridge,
                              box_cv_lambda, box_inner_oof, marginal_prior, huber_ols,
                              LOG_HUBER_C, BINS)

MIN_AGE = 1.0
N_REP = 5
build, pipeline_sh, cur_model, out_dir = sys.argv[1:5]
os.makedirs(out_dir, exist_ok=True)

# ── data ─────────────────────────────────────────────────────────────────────
z = np.load(f'{build}/dap.npz')
breeds = [str(b) for b in z['breeds']]
bi = {b: i for i, b in enumerate(breeds)}
qi = {str(d): i for i, d in enumerate(z['ids'])}
ph = [l.rstrip('\n').split('\t') for l in open(f'{build}/dap_pheno.tsv')][1:]
D = []
for dog, w, age, male, neut, status, breed in ph:
    if dog not in qi or not age or not male:
        continue
    D.append((qi[dog], float(w), float(age), float(male), float(neut) if neut else np.nan, status, breed))
age_d = np.array([r[2] for r in D])
keep = age_d >= MIN_AGE
D = [r for r, k in zip(D, keep) if k]
Xd = z['Q'][[r[0] for r in D]].astype(np.float64)
Xd /= Xd.sum(1, keepdims=True)
yd = np.array([r[1] for r in D]); lyd = np.log(yd)
dad = z['da_prs'][[r[0] for r in D]]
maled = np.array([r[3] for r in D]); neutd = np.array([r[4] for r in D]); aged = np.array([r[2] for r in D])
statd = np.array([r[5] for r in D]); breedd = np.array([r[6] for r in D])
print(f'DAP: {len(D)} genotyped adults (>= {MIN_AGE} y) with weight, sex, age '
      f'({(statd == "Mixed breed").sum()} mixed)')

C = [r for r in json.load(open(f'{build}/customers.json'))
     if r['age'] is not None and r['age'] >= MIN_AGE and r['comp'] and r['da'] is not None]
Xc = np.zeros((len(C), len(breeds)))
unknown = 0.0
for i, r in enumerate(C):
    for b, p in r['comp'].items():
        if b in bi: Xc[i, bi[b]] = p
        else: unknown += p
Xc /= np.maximum(Xc.sum(1, keepdims=True), 1e-9)
yc = np.array([r['y_kg'] for r in C]); lyc = np.log(yc)
dac = np.array([float(r['da']) for r in C])
malec = np.array([0.5 if r['male'] is None else r['male'] for r in C])
nowc = np.array([float(r['pred_now']) for r in C])
print(f'Customers: {len(C)} adults with weight, v2.3 composition and DA PRS; '
      f'{int((malec == 0.5).sum())} without a sex call; label mass outside the panel {unknown:.3f}')
print(f'DA size PRS: DAP mean {dad.mean():.1f} sd {dad.std():.1f}, r(log kg) {np.corrcoef(dad, lyd)[0,1]:.3f}; '
      f'customers mean {dac.mean():.1f} sd {dac.std():.1f}, r(log kg) {np.corrcoef(dac, lyc)[0,1]:.3f}')

src = open(pipeline_sh).read()
akc = join_prior(breeds, extract_dict(src, 'BREED_WEIGHT_KG'))
prior_akc = np.array([akc.get(b, np.nan) for b in breeds])
print(f'{np.isfinite(prior_akc).sum()}/{len(breeds)} labels join the AKC prior table')

def stack_fit(lbp, cov, ly):
    return huber_ols(np.column_stack([lbp] + cov), ly, c=LOG_HUBER_C)

def stack_pred(b, lbp, cov):
    return np.exp(b[0] + np.column_stack([lbp] + cov) @ b[1:])

def summarise(name, P):
    """P: (reps, n) OOF predictions on customers"""
    ms = [metrics(yc, P[r]) for r in range(P.shape[0])]
    agg = {'n': len(yc)}
    for k in ('r', 'mae', 'bias', 'slope'):
        agg[k] = float(np.mean([m[k] for m in ms]))
    for lo, hi in BINS:
        bk = f'bin_{lo}_{hi}'
        if bk in ms[0]:
            agg[bk] = {kk: float(np.mean([m[bk][kk] for m in ms])) for kk in ms[0][bk]}
    print(fmt(name, agg))
    return agg

summary, oof = {}, {}

# ── NOW: shipped model on today's compositions ──────────────────────────────
cm = json.load(open(cur_model))
lbp_now = np.array([sum(p * np.log(max(cm['breed_weight'].get(b, cm['default_kg']), 0.5))
                        for b, p in r['comp'].items()) / max(sum(r['comp'].values()), 1e-9) for r in C])
st = cm['stack']
now = np.exp(st['intercept'] + st['coef_log_breed_pred'] * lbp_now
             + st['covariates']['da_size_prs'] * dac + st['covariates']['male'] * malec)
miss_mass = np.array([sum(p for b, p in r['comp'].items() if b not in cm['breed_weight']) for r in C])
print(f'\nShipped model: {(miss_mass > 0.05).sum()} of {len(C)} customers have >5% ancestry in labels '
      f'it has no weight for (mean {miss_mass.mean():.3f})')
print('Customer adults (bias kg / % by actual weight bin):')
oof['NOW'] = now[None, :]; summary['NOW'] = summarise('NOW shipped, v2.3 comps (in-sample)', oof['NOW'])
oof['REPORTED'] = nowc[None, :]; summary['REPORTED'] = summarise('value in the report today', oof['REPORTED'])

# ── DAP fit (all DAP adults) ────────────────────────────────────────────────
def fit_breed(X, ly, seed=0, box=1.5, lam=None):
    lpr = np.log(marginal_prior(X, np.exp(ly), prior_akc))
    if lam is None:
        lam = box_cv_lambda(X, ly, lpr, seed)
    return fit_box_ridge(X, ly, lpr, lam, box=box), lpr, lam

res = {box: fit_breed(Xd, lyd, box=box) for box in (1.5, 2.0)}
w_d, lpr_d, lam_d = res[1.5]
print(f'\nDAP breed weights: lambda {lam_d} (box 1.5), {res[2.0][2]} (box 2.0)')
# out-of-fold breed prediction on DAP for the stack (same recipe as the shipped fit)
inner_d = np.log(box_inner_oof(Xd, lyd, lpr_d, lam_d, 0))
b_dap = stack_fit(inner_d, [dad, maled], lyd)
print(f'DAP stack: log kg = {b_dap[0]:.3f} + {b_dap[1]:.3f} log breed + {b_dap[2]:.6f} DA + {b_dap[3]:.3f} male '
      f'(male x{np.exp(b_dap[3]):.3f})')
lbp_c_dap = Xc @ w_d
oof['DAP'] = stack_pred(b_dap, lbp_c_dap, [dac, malec])[None, :]
print('\nCustomer adults:')
summary['DAP'] = summarise('DAP only -> customers (external)', oof['DAP'])
w_d2 = res[2.0][0]
b_dap2 = stack_fit(Xd @ w_d2, [dad, maled], lyd)    # in-sample stack: box 2.0 is a sensitivity check only
oof['DAP_box2'] = stack_pred(b_dap2, Xc @ w_d2, [dac, malec])[None, :]
summary['DAP_box2'] = summarise('DAP only, box 2.0 (sensitivity)', oof['DAP_box2'])

# ── customer CV: CUST, DAPR, POOL ───────────────────────────────────────────
n = len(yc)
for k in ('CUST', 'DAPR', 'POOL', 'POOL_w'):
    oof[k] = np.zeros((N_REP, n))
for rep in range(N_REP):
    for te in folds(n, rep):
        tr = np.setdiff1d(np.arange(n), te)
        # CUST: shipped recipe on customers alone
        w_c, lpr_c, lam_c = fit_breed(Xc[tr], lyc[tr], seed=rep)
        inn = np.log(box_inner_oof(Xc[tr], lyc[tr], lpr_c, lam_c, rep))
        b = stack_fit(inn, [dac[tr], malec[tr]], lyc[tr])
        oof['CUST'][rep, te] = stack_pred(b, Xc[te] @ w_c, [dac[te], malec[te]])
        # DAPR: DAP breed weights, calibration from customers
        b = stack_fit(Xc[tr] @ w_d, [dac[tr], malec[tr]], lyc[tr])
        oof['DAPR'][rep, te] = stack_pred(b, Xc[te] @ w_d, [dac[te], malec[te]])
        # POOL: breed weights on DAP + customer training dogs (lambda from the DAP fit)
        for name, wt in (('POOL', 1), ('POOL_w', 4)):
            Xp = np.vstack([Xd] + [Xc[tr]] * wt); lyp = np.concatenate([lyd] + [lyc[tr]] * wt)
            lpr_p = np.log(marginal_prior(Xp, np.exp(lyp), prior_akc))
            w_p = fit_box_ridge(Xp, lyp, lpr_p, lam_d)
            # stack on customers, using DAP-only weights for the in-fold breed term
            # (customer rows are in the pooled fit, so their in-sample term would be optimistic)
            b = stack_fit(Xc[tr] @ w_d, [dac[tr], malec[tr]], lyc[tr])
            oof[name][rep, te] = stack_pred(b, Xc[te] @ w_p, [dac[te], malec[te]])
    print(f'  CV repeat {rep + 1}/{N_REP}', flush=True)
for k, label in (('CUST', 'customer-only refit (v2.3), CV'), ('DAPR', 'DAP breeds + customer calibration, CV'),
                 ('POOL', 'pooled breeds + customer calibration, CV'), ('POOL_w', 'pooled, customers x4, CV')):
    summary[k] = summarise(label, oof[k])

# paired bootstrap MAE differences vs CUST and NOW
rng = np.random.RandomState(7); bs = rng.randint(0, n, size=(2000, n))
ae = {k: np.abs(P.mean(0) - yc) for k, P in oof.items()}
print('\nMAE difference (kg, paired bootstrap 95% CI over customers):')
for k in oof:
    line = f'  {k:<10}'
    for ref in ('NOW', 'CUST'):
        if k != ref:
            d = (ae[k] - ae[ref])[bs].mean(1)
            summary[k][f'dmae_vs_{ref}'] = [float(d.mean()), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]
            line += f'  vs {ref} {d.mean():+5.2f} [{np.percentile(d, 2.5):+5.2f},{np.percentile(d, 97.5):+5.2f}]'
    print(line)

# ── DAP-internal CV: what else explains weight ──────────────────────────────
print('\nDAP 5-fold CV (DAP adults; log-space stack on out-of-fold breed term):')
nd = len(yd)
variants = {'breed only': [], '+DA PRS': [dad], '+DA +male': [dad, maled],
            '+DA +male +neutered': [dad, maled, np.nan_to_num(neutd, nan=np.nanmean(neutd))],
            '+DA +male +neut +age<2': [dad, maled, np.nan_to_num(neutd, nan=np.nanmean(neutd)), (aged < 2).astype(float)]}
pd_ = {k: np.zeros(nd) for k in variants}
coefs = {}
for te in folds(nd, 11):
    tr = np.setdiff1d(np.arange(nd), te)
    w_f, lpr_f, lam_f = fit_breed(Xd[tr], lyd[tr], lam=lam_d)
    inn = np.log(box_inner_oof(Xd[tr], lyd[tr], lpr_f, lam_f, 3))
    for k, cov in variants.items():
        b = stack_fit(inn, [c[tr] for c in cov], lyd[tr])
        pd_[k][te] = stack_pred(b, Xd[te] @ w_f, [c[te] for c in cov])
        coefs.setdefault(k, []).append(b)
mixed = statd == 'Mixed breed'
for k in variants:
    m_all, m_mix, m_pure = metrics(yd, pd_[k]), metrics(yd[mixed], pd_[k][mixed]), metrics(yd[~mixed], pd_[k][~mixed])
    b = np.mean(coefs[k], 0)
    extra = ' '.join(f'x{np.exp(c):.3f}' for c in b[2:]) if len(b) > 2 else ''
    print(f'  {k:<24} all r={m_all["r"]:.3f} MAE={m_all["mae"]:.2f} | mixed r={m_mix["r"]:.3f} MAE={m_mix["mae"]:.2f} '
          f'| purebred r={m_pure["r"]:.3f} MAE={m_pure["mae"]:.2f}   multipliers {extra}')
    summary[f'DAPCV {k}'] = {'all': m_all, 'mixed': m_mix, 'purebred': m_pure}

# ── purebred sanity: DAP owner-labelled purebreds vs fitted breed weight ────
print('\nImplied 100% purebred weight (kg): DAP fit vs shipped vs AKC prior vs DAP purebred median (n)')
pure = {}
for i, r in enumerate(D):
    if r[5] == 'Purebred':
        top = int(np.argmax(Xd[i]))
        if Xd[i, top] >= 0.85:
            pure.setdefault(breeds[top], []).append(yd[i])
for b in sorted(pure, key=lambda b: -len(pure[b]))[:40]:
    j = bi[b]
    print(f'  {b:<32} DAP {np.exp(w_d[j]):5.1f}  shipped {cm["breed_weight"].get(b, float("nan")):5.1f}  '
          f'AKC {prior_akc[j]:5.1f}  median {np.median(pure[b]):5.1f} (n={len(pure[b])})')

json.dump({'summary': summary, 'dap_lambda': lam_d, 'dap_stack': [float(x) for x in b_dap],
           'n_dap': nd, 'n_customers': n}, open(f'{out_dir}/dap_fit_summary.json', 'w'), indent=1)
np.savez(f'{out_dir}/dap_fit.npz', breeds=np.array(breeds), w_dap=w_d, lpr_dap=lpr_d, b_dap=b_dap,
         **{f'oof_{k}': v for k, v in oof.items()}, yc=yc)
print(f'\nwrote {out_dir}/dap_fit_summary.json')
