#!/usr/bin/env python3
"""
Hayward 2016 dosages on the Parker site grid (dosage of panel a1), plus the
same-animal check against Parker and Dog10K.

    python3 hayward_dosage.py HAYWARD_PREFIX SITES.tsv PARKER_PREFIX CLUST.txt DOG10K_GT.txt DOG10K_SAMPLES.txt OUT_NPZ

SITES.tsv is hayward_sites.py's output. Saves OUT_NPZ with dos (sites x dogs,
-1 missing), ids, labels (breeds.txt), and prints cross-panel IBS matches.
"""
import sys
import numpy as np
from merged_panel import read_parker, read_dog10k

hw, sites, pprefix, clust, gt, samples, out = sys.argv[1:8]
fam = [l.split() for l in open(hw + '.fam')]
ids = [f[1] for f in fam]
n_ind = len(ids)
n_snp = sum(1 for _ in open(hw + '.bim'))
raw = np.fromfile(hw + '.bed', dtype=np.uint8)
assert tuple(raw[:3]) == (0x6c, 0x1b, 1), 'not a SNP-major bed'
bps = (n_ind + 3) // 4
data = raw[3:].reshape(n_snp, bps)
codes = np.empty((n_snp, bps * 4), dtype=np.uint8)
for k in range(4):
    codes[:, k::4] = (data >> (2 * k)) & 3
dosA1 = np.array([2, -1, 1, 0], dtype=np.int8)[codes[:, :n_ind]]   # dosage of Hayward A1

# Parker grid
P, p_ids, p_truth, keys, a1 = read_parker(pprefix, clust)
key_index = {k: i for i, k in enumerate(keys)}
H = np.full((len(keys), n_ind), -1, dtype=np.int8)
n = 0
for l in open(sites):
    if l.startswith('hw_index'): continue
    hi, _, c, pos, pa1, pa2, rel = l.rstrip('\n').split('\t')
    j = key_index[(c, pos)]
    d = dosA1[int(hi)]
    if rel in ('swap', 'flipswap'):      # Hayward A1 is the panel a2
        d = np.where(d >= 0, 2 - d, -1).astype(np.int8)
    H[j] = d; n += 1
print(f'Hayward: {n_ind} dogs, {n} panel sites filled of {len(keys)}', file=sys.stderr)

labels = {}
for l in open(hw.rsplit('/', 1)[0] + '/breeds.txt', newline=''):
    p = l.rstrip('\r\n').split('\t')
    if len(p) >= 2 and p[0] != 'dogID': labels[p[0]] = p[1]
np.savez_compressed(out, dos=H, ids=np.array(ids), labels=np.array([labels.get(i, '') for i in ids]))

def ibs_matrix(A, B):
    okA = (A >= 0).mean(1) >= 0.95; okB = (B >= 0).mean(1) >= 0.95
    use = np.where(okA & okB)[0]
    A = A[use].astype(np.float32); B = B[use].astype(np.float32)
    mA = (A >= 0); mB = (B >= 0)
    A[~mA] = 0; B[~mB] = 0
    n_both = mA.T.astype(np.float32) @ mB.astype(np.float32)
    same = sum(((A == g) & mA).astype(np.float32).T @ ((B == g) & mB).astype(np.float32) for g in (0, 1, 2))
    off2 = ((A == 0) & mA).astype(np.float32).T @ ((B == 2) & mB).astype(np.float32) + ((A == 2) & mA).astype(np.float32).T @ ((B == 0) & mB).astype(np.float32)
    off1 = n_both - same - off2
    return (2 * same + off1) / (2 * np.maximum(n_both, 1)), len(use)

def report(name, ibs, a_ids, a_lab, b_ids, b_lab, thr=0.97):
    best = ibs.max(1)
    print(f'\n== {name}: best-match IBS percentiles 50/90/99/max = {np.round(np.percentile(best, [50, 90, 99, 100]), 4)}')
    hits = np.argwhere(ibs >= thr)
    print(f'pairs with IBS >= {thr}: {len(hits)}')
    for i, j in hits[:80]:
        print(f'  {a_ids[i]}\t{a_lab(a_ids[i])}\t{b_ids[j]}\t{b_lab(b_ids[j])}\t{ibs[i, j]:.4f}')

ibs, nu = ibs_matrix(H, P); print(f'Hayward vs Parker at {nu} sites', file=sys.stderr)
report('Hayward vs Parker', ibs, ids, lambda i: labels.get(i, ''), p_ids, lambda i: p_truth.get(i, ''))
# Hayward vs Hayward (internal duplicates), upper triangle only
ibs_hh, nu = ibs_matrix(H, H); np.fill_diagonal(ibs_hh, 0); ibs_hh = np.triu(ibs_hh)
report('Hayward internal', ibs_hh, ids, lambda i: labels.get(i, ''), ids, lambda i: labels.get(i, ''))
D, d_ids, d_truth = read_dog10k(gt, samples, key_index, a1)
ibs, nu = ibs_matrix(H, D); print(f'Hayward vs Dog10K at {nu} sites', file=sys.stderr)
report('Hayward vs Dog10K', ibs, ids, lambda i: labels.get(i, ''), d_ids, lambda i: d_truth.get(i, ''))
