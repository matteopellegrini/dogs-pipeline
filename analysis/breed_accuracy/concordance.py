#!/usr/bin/env python3
"""
Same-animal check across reference panels.

    python3 concordance.py PARKER_PREFIX CLUST.txt DOG10K_GT.txt DOG10K_SAMPLES.txt [--min-ibs 0.97]

A dog genotyped on the Parker array AND sequenced in Dog10K would be counted
twice in the merged reference and would leak into the "held-out" Dog10K
benchmark. Two records of the same animal agree at ~all sites (array error +
imputation error only); unrelated dogs of one breed agree at ~70-80% of sites,
first-degree relatives ~85%. We report every cross-panel pair with identity-
by-state above --min-ibs, plus the distribution so the gap is visible.
IBS here = mean over sites typed in both of (2 - |dosage difference|) / 2.
"""
import sys
import numpy as np
from merged_panel import read_parker, read_dog10k

prefix, clust, gt, samples = sys.argv[1:5]
min_ibs = float(sys.argv[sys.argv.index('--min-ibs') + 1]) if '--min-ibs' in sys.argv else 0.97

A, a_ids, a_truth, keys, a1 = read_parker(prefix, clust)
key_index = {k: i for i, k in enumerate(keys)}
B, b_ids, b_truth = read_dog10k(gt, samples, key_index, a1)
print(f'Parker {A.shape[1]} dogs, Dog10K {B.shape[1]} dogs, {A.shape[0]} sites', file=sys.stderr)

# Use sites typed in >=95% of both panels; mean-impute the rest so matrices are dense.
okA = (A >= 0).mean(1) >= 0.95
okB = (B >= 0).mean(1) >= 0.95
use = np.where(okA & okB)[0]
A = A[use].astype(np.float32); B = B[use].astype(np.float32)
A[A < 0] = np.nan; B[B < 0] = np.nan
print(f'{len(use)} sites used', file=sys.stderr)

# |a-b| summed over sites via (a-b)^2 = a^2 - 2ab + b^2 on 0/1/2 dosages: equal iff squared diff 0.
# Exact IBS needs |a-b|, so compute with one-hot encodings: IBS2 = matches of identical dosage;
# IBS = (2*n_same + n_offby1) / (2*n_sites). One-hot: for g in 0,1,2 count agreement.
mA = ~np.isnan(A); mB = ~np.isnan(B)
A0 = np.nan_to_num(A); B0 = np.nan_to_num(B)
n_both = mA.T.astype(np.float32) @ mB.astype(np.float32)               # sites typed in both
same = np.zeros((A.shape[1], B.shape[1]), dtype=np.float32)
for g in (0.0, 1.0, 2.0):
    ag = ((A0 == g) & mA).astype(np.float32); bg = ((B0 == g) & mB).astype(np.float32)
    same += ag.T @ bg
# off-by-one = typed-both - same - off-by-two; off-by-two = (A==0&B==2)+(A==2&B==0)
a0 = ((A0 == 0) & mA).astype(np.float32); a2 = ((A0 == 2) & mA).astype(np.float32)
b0 = ((B0 == 0) & mB).astype(np.float32); b2 = ((B0 == 2) & mB).astype(np.float32)
off2 = a0.T @ b2 + a2.T @ b0
off1 = n_both - same - off2
ibs = (2 * same + off1) / (2 * np.maximum(n_both, 1))

best = ibs.max(1)
print('cross-panel best-match IBS per Parker dog: percentiles 50/90/99/max =',
      np.round(np.percentile(best, [50, 90, 99, 100]), 4), file=sys.stderr)
hits = np.argwhere(ibs >= min_ibs)
print(f'pairs with IBS >= {min_ibs}: {len(hits)}')
print('parker_id\tparker_breed\tdog10k_id\tdog10k_breed\tibs\tn_sites')
for i, j in hits:
    print(f'{a_ids[i]}\t{a_truth.get(a_ids[i], "?")}\t{b_ids[j]}\t{b_truth.get(b_ids[j], "?")}\t{ibs[i, j]:.4f}\t{int(n_both[i, j])}')
# Also show the top 10 near-misses so the gap between "same animal" and "close relative" is visible.
flat = np.argsort(best)[::-1][:12]
print('\ntop cross-panel matches (incl. below threshold):')
for i in flat:
    j = int(ibs[i].argmax())
    print(f'{a_ids[i]}\t{a_truth.get(a_ids[i], "?")}\t{b_ids[j]}\t{b_truth.get(b_ids[j], "?")}\t{ibs[i, j]:.4f}')
