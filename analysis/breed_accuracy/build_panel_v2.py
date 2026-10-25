#!/usr/bin/env python3
"""
Breed panel v2: Parker + Dog10K + Hayward (+ DAP later) on the shipped 131,353-site grid.

    python3 build_panel_v2.py --out DIR [--min-n 6] [--holdout 0.2] [--seed 1]
        [--no-hayward] [--dap DAP_NPZ] [--drop IDS.txt ...] [--lasso 0.3]

Inputs are fixed paths relative to the repo (see CONFIG). Writes DIR/{phat.npy,
breeds.txt,sites.tsv,reference_counts.tsv} from ALL dogs, and benchmarks a
stratified hold-out (per breed, per source) scored with Stage 9's estimator
(non-negative lasso, lambda 0.3) against frequencies from the training dogs
only. Reports top-1 by reference size and by source, mean self-proportion, and
the spurious mass on synthetic F1 crosses of held-out purebreds.
"""
import json, sys, os, re
from collections import Counter, defaultdict
import numpy as np
from scipy.optimize import nnls
from scipy.linalg import solve_triangular
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from merged_panel import read_parker, read_dog10k

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
S = os.environ.get('SCRATCH', '/private/tmp/claude-501/-Users-matteopellegrini-Downloads-dogs/0fb26461-2afd-4142-936c-a58295d35959/scratchpad')
CONFIG = dict(
    parker=f'{ROOT}/COSMO/analysis/cosmo_parker_full', clust=f'{ROOT}/COSMO/analysis/parker_clust.txt',
    dog10k_gt=f'{S}/dog10k_gt.txt', dog10k_samples=f'{S}/dog10k_samples.txt',
    harmonized=f'{S}/harmonized.json', hayward=f'{S}/hayward_dos.npz',
    hayward_map=f'{ROOT}/analysis/breed_accuracy/hayward_label_map.json',
    sites=f'{ROOT}/breed_panel/sites.tsv', dups=f'{ROOT}/analysis/breed_accuracy/parker_dog10k_duplicates.tsv',
)
args = sys.argv[1:]
opt = lambda k, d: args[args.index(k) + 1] if k in args else d
OUT = opt('--out', f'{S}/panel_v2'); MIN_N = int(opt('--min-n', 6)); HOLD = float(opt('--holdout', 0.2))
SEED = int(opt('--seed', 1)); LASSO = float(opt('--lasso', 0.3)); USE_HW = '--no-hayward' not in args
DAP = opt('--dap', None)
drop = set()
for i, a in enumerate(args):
    if a == '--drop': drop |= set(l.split()[0] for l in open(args[i + 1]) if l.strip() and not l.startswith('#'))
os.makedirs(OUT, exist_ok=True)

# --- load everything on the Parker grid -------------------------------------
P, p_ids, p_truth, keys, a1 = read_parker(CONFIG['parker'], CONFIG['clust'])
key_index = {k: i for i, k in enumerate(keys)}
label_of_code = json.load(open(CONFIG['harmonized']))
for l in open(CONFIG['dups']):
    f = l.split('\t')
    if len(f) >= 5 and re.match(r'^0\.9', f[4]): drop.add(f[0])          # Parker copy of a Dog10K dog
cols, ids, labels, source = [], [], [], []
def add(dos, dids, labs, src):
    keep = [j for j, d in enumerate(dids) if d not in drop and labs[j]]
    cols.append(dos[:, keep]); ids.extend(dids[j] for j in keep); labels.extend(labs[j] for j in keep); source.extend([src] * len(keep))
    print(f'{src}: {len(keep)} dogs kept of {len(dids)}', file=sys.stderr)
add(P, p_ids, [label_of_code.get(p_truth.get(i, ''), None) for i in p_ids], 'parker')
D, d_ids, d_truth = read_dog10k(CONFIG['dog10k_gt'], CONFIG['dog10k_samples'], key_index, a1)
add(D, d_ids, [label_of_code.get(d_truth.get(i, ''), None) for i in d_ids], 'dog10k')
if USE_HW:
    z = np.load(CONFIG['hayward']); hmap = json.load(open(CONFIG['hayward_map']))
    add(z['dos'], list(z['ids']), [hmap.get(str(l)) for l in z['labels']], 'hayward')
if DAP:
    z = np.load(DAP); add(z['dos'], list(z['ids']), [str(l) or None for l in z['labels']], 'dap')
X = np.concatenate(cols, axis=1); del cols
labels = np.array(labels); source = np.array(source); ids = np.array(ids)
# restrict to the shipped site grid
site_rows = []
with open(CONFIG['sites']) as fh:
    next(fh)
    for l in fh:
        c, pos, s1, s2 = l.rstrip('\n').split('\t'); site_rows.append(key_index[(c, pos)])
site_rows = np.array(site_rows); X = X[site_rows]
print(f'matrix: {X.shape[0]} sites x {X.shape[1]} dogs; labels {len(set(labels))}', file=sys.stderr)

# --- min-n and the full panel -------------------------------------------------
cnt = Counter(labels)
breeds = sorted(b for b, n in cnt.items() if n >= MIN_N and b != 'COYOTE')
bidx = {b: j for j, b in enumerate(breeds)}
def freqs(cols_sel):
    Pm = np.full((X.shape[0], len(breeds)), 0.5, dtype=np.float32)
    for j, b in enumerate(breeds):
        sel = [i for i in cols_sel if labels[i] == b]
        if not sel: continue
        sub = X[:, sel].astype(np.float32); sub[sub < 0] = np.nan
        with np.errstate(invalid='ignore'): m = np.nanmean(sub, axis=1) / 2.0
        Pm[:, j] = np.nan_to_num(m, nan=0.5)
    return Pm
Pall = freqs(range(X.shape[1]))
np.save(f'{OUT}/phat.npy', Pall)
open(f'{OUT}/breeds.txt', 'w').write('\n'.join(breeds) + '\n')
os.system(f'cp {CONFIG["sites"]} {OUT}/sites.tsv')
with open(f'{OUT}/reference_counts.tsv', 'w') as fh:
    fh.write('breed\tn\tparker\tdog10k\thayward\tdap\n')
    for b in sorted(cnt, key=lambda b: -cnt[b]):
        sel = labels == b
        fh.write(f'{b}\t{cnt[b]}\t' + '\t'.join(str(int((source[sel] == s).sum())) for s in ('parker', 'dog10k', 'hayward', 'dap')) + ('' if b in bidx else '\t(below min-n)') + '\n')
print(f'panel: {len(breeds)} populations with n>={MIN_N} ({sum(cnt[b] for b in breeds)} dogs); {len(cnt) - len(breeds)} labels below min-n', file=sys.stderr)

# --- held-out benchmark -------------------------------------------------------
rng = np.random.default_rng(SEED)
held = []
for b in breeds:
    for s in ('parker', 'dog10k', 'hayward', 'dap'):
        idx = np.where((labels == b) & (source == s))[0]
        if len(idx) >= 4: held.extend(rng.choice(idx, size=max(1, int(round(HOLD * len(idx)))), replace=False))
held = np.array(sorted(held)); train = np.setdiff1d(np.arange(X.shape[1]), held)
Ptr = freqs(train)
# Stage 9 estimator on full-grid (missing sites mean-imputed to the panel-average dosage).
def score(x):
    valid = x >= 0
    Pv = Ptr[valid]; xv = x[valid].astype(np.float32)
    A = np.vstack([Pv, np.ones((1, Pv.shape[1]), dtype=np.float32)]); b = np.hstack([xv, [1.0]]).astype(np.float32)
    G = (A.T @ A).astype(np.float64); rhs = (A.T @ b).astype(np.float64)
    if LASSO: rhs -= 0.5 * LASSO * float(np.mean(np.diag(Pv.T @ Pv)))
    G[np.diag_indices_from(G)] += 1e-6
    R = np.linalg.cholesky(G).T
    q, _ = nnls(R, solve_triangular(R.T, rhs, lower=True))
    return q / (q.sum() + 1e-12)
res = []
for i in held:
    q = score(X[:, i]); top = int(q.argmax())
    res.append((labels[i], source[i], breeds[top], q[bidx[labels[i]]]))
res_by = defaultdict(list)
for lab, src, pred, selfp in res:
    n = cnt[lab]; bucket = 'n=6-7' if n < 8 else 'n=8-19' if n < 20 else 'n>=20'
    for k in ('all', f'src:{src}', bucket): res_by[k].append((pred == lab, selfp))
print(f'\nHELD-OUT ({len(held)} dogs, {HOLD:.0%} per breed per source, lasso {LASSO}, min-n {MIN_N}, {len(breeds)} populations):')
for k in sorted(res_by): v = res_by[k]; print(f'  {k:14s} top-1 {100*np.mean([a for a,_ in v]):5.1f}%  mean self {100*np.mean([s for _,s in v]):5.1f}%  (n={len(v)})')
wrong = [(lab, src, pred) for lab, src, pred, _ in res if pred != lab]
print('  errors by true breed:', Counter(l for l, _, _ in wrong).most_common(15))
print('  errors (true -> predicted):', Counter(f'{l}->{p}' for l, _, p in wrong).most_common(12))
# synthetic F1 crosses from held-out purebreds: spurious mass = share assigned to non-parents
pure = [i for i in held if score(X[:, i]).argmax() == bidx[labels[i]]][:400]
sp, err, found = [], [], 0
for t in range(200):
    i, j = rng.choice(pure, 2, replace=False)
    if labels[i] == labels[j]: continue
    xi, xj = X[:, i], X[:, j]; ok = (xi >= 0) & (xj >= 0)
    # one allele from each parent per site
    child = np.full_like(xi, -1); child[ok] = rng.binomial(1, xi[ok] / 2) + rng.binomial(1, xj[ok] / 2)
    q = score(child); a, b_ = bidx[labels[i]], bidx[labels[j]]
    sp.append(1 - q[a] - q[b_]); err.append(abs(q[a] - .5) + abs(q[b_] - .5)); found += (q[a] > .1 and q[b_] > .1)
print(f'  synthetic F1 ({len(sp)}): spurious mass {100*np.mean(sp):.1f}%  proportion error {100*np.mean(err)/2:.1f}pp  both parents found {found}/{len(sp)}')
