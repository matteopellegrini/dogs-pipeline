#!/usr/bin/env python3
"""
Self-check the Parker reference panel: does each reference dog predict as its
own breed?

    python3 check_reference.py PLINK_PREFIX PHAT.txt CLUST.txt

Runs exactly the projection Stage 9 runs for a customer dog — NNLS of
dosage-of-a1 against the SCOPE Phat matrix with a sum-to-1 row appended — but
on the reference individuals themselves.

This is an IN-SAMPLE check: every dog contributed to its own breed's allele
frequencies, so it is the easiest possible test and accuracy here is an upper
bound. That is what makes failures informative. A reference dog that cannot
recover its own breed, when that breed's frequencies were estimated partly from
it, marks a reference population that is too small, mislabelled, or not
genetically distinct from a neighbour.

Speed note: solving 1356 separate NNLS problems on a 143933 x 177 design is
prohibitive, so the normal equations are used instead. For A = [P; 1...1],
min ||Aq - b|| is equivalent to min ||Rq - R^-T A^T b|| where A^T A = R^T R
(Cholesky). Each solve becomes 177x177 and the answer is identical.
"""
import sys

import numpy as np
from scipy.linalg import cho_factor, solve_triangular
from scipy.optimize import nnls

CHUNK = 200


def read_bed(prefix, n_ind, n_snp):
    raw = np.fromfile(prefix + '.bed', dtype=np.uint8)
    assert raw[0] == 108 and raw[1] == 27 and raw[2] == 1, 'not a PLINK1 SNP-major bed'
    bps = (n_ind + 3) // 4
    data = raw[3:].reshape(n_snp, bps)
    codes = np.empty((n_snp, bps * 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k::4] = (data >> (2 * k)) & 3
    codes = codes[:, :n_ind]
    # PLINK1: 0=hom a1, 1=missing, 2=het, 3=hom a2  ->  dosage of a1
    lut = np.array([2, -1, 1, 0], dtype=np.int8)
    return lut[codes]


def main():
    prefix, phat_path, clust_path = sys.argv[1:4]

    fam = [l.split() for l in open(prefix + '.fam')]
    ids = [f[1] for f in fam]
    n_ind = len(fam)
    n_snp = sum(1 for _ in open(prefix + '.bim'))

    truth = {}
    order = []
    for line in open(clust_path):
        p = line.split()
        if len(p) >= 3:
            truth[p[1]] = p[2]
            if p[2] not in order:
                order.append(p[2])
    labels = order
    K = len(labels)

    P = np.loadtxt(phat_path, dtype=np.float32)
    assert P.shape == (n_snp, K), f'Phat {P.shape} vs {(n_snp, K)}'

    dos = read_bed(prefix, n_ind, n_snp)
    miss = dos < 0
    print(f"individuals {n_ind}  snps {n_snp}  breeds {K}")
    print(f"missing genotypes: {100*miss.mean():.3f}%")

    exp = (2.0 * P.mean(axis=1)).astype(np.float32)   # expected dosage under panel average

    # Normal equations for A = [P ; ones(1,K)]
    G = (P.T @ P) + np.ones((K, K), dtype=np.float32)
    G = G.astype(np.float64)
    G[np.diag_indices_from(G)] += 1e-6
    R = np.linalg.cholesky(G).T                      # G = R^T R

    # Leave-one-out. Phat is the per-breed empirical allele frequency (verified:
    # corr 0.9997+ against the genotypes), so a dog's own contribution can be
    # removed analytically instead of re-running SCOPE 1355 times:
    #     p_loo = (n*p - x/2)/(n-1)   =>   delta = p_loo - p = (p - x/2)/(n-1)
    # Only column b changes, so the Gram matrix takes a rank-2 update rather
    # than a full recomputation.
    #
    # Without this the test is circular. A breed with n=1 has Phat equal to that
    # single dog's genotype and recovers itself perfectly by construction, which
    # says nothing about whether the panel can identify that breed in a new dog.
    loo = '--loo' in sys.argv
    idx_of = {b: i for i, b in enumerate(labels)}
    counts = {}
    for d, b in truth.items():
        counts[b] = counts.get(b, 0) + 1

    top1 = []
    untestable = 0
    for s in range(0, n_ind, CHUNK):
        e = min(s + CHUNK, n_ind)
        X = dos[:, s:e].astype(np.float32)
        m = miss[:, s:e]
        if m.any():
            X[m] = np.broadcast_to(exp[:, None], X.shape)[m]
        B = (P.T @ X) + 1.0                           # (K, chunk)
        for j in range(B.shape[1]):
            i = s + j
            bb = B[:, j].astype(np.float64)
            Gd = G
            b_lab = truth.get(ids[i])
            if loo and b_lab is not None and counts.get(b_lab, 0) > 1:
                col = idx_of[b_lab]
                x = X[:, j]
                delta = ((P[:, col] - x / 2.0) / (counts[b_lab] - 1)).astype(np.float32)
                u = (P.T @ delta).astype(np.float64)
                Gd = G.copy()
                Gd[col, :] += u
                Gd[:, col] += u
                Gd[col, col] += float(delta @ delta)
                bb[col] += float(delta @ x)
            elif loo and b_lab is not None:
                untestable += 1
                top1.append(None)
                continue
            try:
                Rd = np.linalg.cholesky(Gd).T
            except np.linalg.LinAlgError:
                top1.append(None)
                continue
            y = solve_triangular(Rd.T, bb, lower=True)
            q, _ = nnls(Rd, y)
            top1.append(labels[int(np.argmax(q))])
        print(f"  scored {e}/{n_ind}", end='\r', flush=True)
    print(" " * 30, end='\r')
    if loo:
        print(f"mode: LEAVE-ONE-OUT   (untestable, n=1 breeds: {untestable})")
    else:
        print("mode: IN-SAMPLE (upper bound; circular for small breeds)")

    # ── results ──────────────────────────────────────────────────
    have = [(i, ids[i]) for i in range(n_ind) if ids[i] in truth and top1[i] is not None]
    correct = [i for i, d in have if top1[i] == truth[d]]
    print(f"\nreference dogs scored: {len(have)}")
    print(f"top-1 == own breed   : {len(correct)}/{len(have)} "
          f"({100*len(correct)/len(have):.1f}%)\n")

    from collections import defaultdict
    per = defaultdict(lambda: [0, 0])
    wrong = []
    for i, d in have:
        b = truth[d]
        per[b][1] += 1
        if top1[i] == b:
            per[b][0] += 1
        else:
            wrong.append((d, b, top1[i]))

    bad = sorted(((c / n, n, b) for b, (c, n) in per.items()), key=lambda x: (x[0], -x[1]))
    print("worst-performing reference populations:")
    print(f"  {'breed':<34} {'n':>3}  {'self-recovered':>14}")
    for frac, n, b in bad[:18]:
        print(f"  {b:<34} {n:>3}  {int(frac*n):>6}/{n} ({100*frac:>3.0f}%)")

    sizes = {b: per[b][1] for b in per}
    small = [b for b in per if sizes[b] <= 3]
    big = [b for b in per if sizes[b] >= 8]
    def acc(bs):
        c = sum(per[b][0] for b in bs); t = sum(per[b][1] for b in bs)
        return 100 * c / t if t else float('nan'), t
    a_s, n_s = acc(small)
    a_b, n_b = acc(big)
    print(f"\nreference size vs self-recovery:")
    print(f"  populations with n<=3 : {a_s:.0f}% correct over {n_s} dogs ({len(small)} breeds)")
    print(f"  populations with n>=8 : {a_b:.0f}% correct over {n_b} dogs ({len(big)} breeds)")

    print(f"\nall {len(wrong)} misassigned reference dogs (first 30):")
    for d, b, got in wrong[:30]:
        print(f"  {d:<18} labelled {b:<16} -> predicted {got}")


if __name__ == '__main__':
    main()
