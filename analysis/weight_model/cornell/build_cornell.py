#!/usr/bin/env python3
"""
Breed compositions for the weighed Cornell dogs (Hayward 2016, CanineHD array), so
the stage-11 weight model can be trained on them the way build_dap.py does for DAP.

    build_cornell.py HAYWARD_PREFIX HAYWARD_SITES.tsv PHENOTYPES.txt BREEDS.txt PANEL_DIR OUT_DIR

HAYWARD_SITES.tsv is hayward_sites.py's output (array SNP -> canFam4 panel site and
allele relation). Dosage of the panel a1 is checked against PANEL_DIR/sites.tsv.
weight_individual in the Dryad phenotype file is sqrt(kg). Writes:
  cornell.npz        ids, breeds, Q (dogs x populations), sites_used
  cornell_pheno.tsv  dog_id, weight_kg, male, breed
"""
import sys, os
import numpy as np
from scipy.optimize import nnls
from scipy.linalg import solve_triangular

LASSO = 0.3
hw, sites_tsv, pheno, breeds_txt, panel_dir, out_dir = sys.argv[1:7]
os.makedirs(out_dir, exist_ok=True)

# weighed dogs
hdr, *rows = [l.rstrip('\n').split('\t') for l in open(pheno)]
wi = hdr.index('weight_individual')
wt = {r[0]: float(r[wi]) ** 2 for r in rows if r[wi] not in ('NA', '-9', '')}
txt = open(breeds_txt, newline='').read().replace('\r\n', '\n').replace('\r', '\n')
label = dict(l.split('\t')[:2] for l in txt.split('\n')[1:] if '\t' in l)

fam = [l.split() for l in open(hw + '.fam')]
col = [i for i, f in enumerate(fam) if f[1] in wt]
ids = [fam[i][1] for i in col]
male = [1 if fam[i][4] == '1' else 0 if fam[i][4] == '2' else None for i in col]
n_ind, n = len(fam), len(col)

# panel sites
psites = [l.split('\t') for l in open(f'{panel_dir}/sites.tsv')][1:]
pidx = {(c, int(p)): (j, a1) for j, (c, p, a1, _) in enumerate(psites)}
P = np.load(f'{panel_dir}/phat.npy').astype(np.float32)
breeds = [l.strip() for l in open(f'{panel_dir}/breeds.txt') if l.strip()]

# genotypes (SNP-major bed) -> panel a1 dosage
n_snp = sum(1 for _ in open(hw + '.bim'))
raw = np.fromfile(hw + '.bed', dtype=np.uint8)
assert tuple(raw[:3]) == (0x6c, 0x1b, 1), 'not a SNP-major bed'
data = raw[3:].reshape(n_snp, (n_ind + 3) // 4)
cb, sh = np.array(col) // 4, 2 * (np.array(col) % 4)
X = np.full((P.shape[0], n), 255, dtype=np.uint8)
lut = np.array([2, 255, 1, 0], dtype=np.uint8)          # bed code -> dosage of the bim A1
stats = {'mapped': 0, 'a1_swap': 0, 'a1_mismatch': 0}
for l in open(sites_tsv):
    if l.startswith('hw_index'): continue
    hi, _, c, pos, a1, a2, rel = l.rstrip('\n').split('\t')
    k = pidx.get((c, int(pos)))
    if k is None: continue
    j, pa1 = k
    d = lut[(data[int(hi), cb] >> sh) & 3]
    swap = rel in ('swap', 'flipswap')                  # array A1 is the Parker a2
    if pa1 != a1:                                       # v2.3 a1 vs the Parker grid a1
        if pa1 == a2: swap = not swap; stats['a1_swap'] += 1
        else: stats['a1_mismatch'] += 1; continue
    if swap: d = np.where(d == 255, 255, 2 - d).astype(np.uint8)
    X[j] = np.where(d == 255, 255, d * 100).astype(np.uint8)
    stats['mapped'] += 1
print(f'{n} weighed dogs; panel sites {stats}')

miss = X == 255
site_ok = ~miss.all(1)
print(f'panel sites present {site_ok.sum()}/{P.shape[0]}; per-dog extra missing max {int(miss[site_ok].sum(0).max())}')
Pv = P[site_ok]
pen = 0.5 * LASSO * float(np.mean(np.einsum('ij,ij->j', Pv, Pv)))
G = Pv.T.astype(np.float64) @ Pv + 1.0
G[np.diag_indices_from(G)] += 1e-6
R = np.linalg.cholesky(G).T
Q = np.zeros((n, len(breeds)), dtype=np.float32)
used = np.full(n, int(site_ok.sum()))
xb_all = X[site_ok]
for j in range(n):
    bad = xb_all[:, j] == 255
    x = xb_all[:, j].astype(np.float64) / 100.0
    if bad.any():                                       # refit this dog on its own sites
        v = ~bad; Pd = Pv[v]
        Gd = Pd.T.astype(np.float64) @ Pd + 1.0
        Gd[np.diag_indices_from(Gd)] += 1e-6
        Rd = np.linalg.cholesky(Gd).T
        rhs = Pd.T.astype(np.float64) @ x[v] + 1.0 - 0.5 * LASSO * float(np.mean(np.einsum('ij,ij->j', Pd, Pd)))
        q, _ = nnls(Rd, solve_triangular(Rd.T, rhs, lower=True)); used[j] = int(v.sum())
    else:
        rhs = Pv.T.astype(np.float64) @ x + 1.0 - pen
        q, _ = nnls(R, solve_triangular(R.T, rhs, lower=True))
    Q[j] = q / (q.sum() + 1e-12)
    if (j + 1) % 500 == 0: print(f'  scored {j + 1}/{n}', flush=True)
np.savez_compressed(f'{out_dir}/cornell.npz', ids=np.array(ids), breeds=np.array(breeds), Q=Q, sites_used=used)
with open(f'{out_dir}/cornell_pheno.tsv', 'w') as fo:
    fo.write('dog_id\tweight_kg\tmale\tbreed\n')
    for i, m in zip(ids, male):
        fo.write(f'{i}\t{wt[i]:.3f}\t{"" if m is None else m}\t{label.get(i, "")}\n')
print(f'wrote {out_dir}/cornell.npz')
