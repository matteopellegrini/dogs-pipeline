#!/usr/bin/env python3
"""
Head-to-head benchmark: NNLS vs SCOPE Qhat on the same held-out dogs.

    python3 benchmark_eval.py BENCH_PREFIX [SCOPE_QHAT.txt]

BENCH_PREFIX is the output of make_merged_plink.py --holdout, which writes a
.bed containing every dog but a .frq.strat computed from the TRAINING dogs only.
The held-out IDs are in BENCH_PREFIX.holdout.

This is the comparison that was missing. The numbers quoted so far came from
different protocols — SCOPE Qhat in-sample (96.83%) against NNLS leave-one-out
(97.20%) — which cannot settle which estimator is better. Here both score the
same dogs against the same training frequencies, so the difference is the
estimator and nothing else.

Why the held-out dogs stay in the .bed: SCOPE estimates its latent subspace over
the whole cohort, so it cannot score a dog that is absent. Deployment would
likewise run the reference panel plus the customer dog. Their genotypes never
enter the frequencies, which is what makes the test honest.
"""
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import nnls


def main():
    pre = sys.argv[1]
    qhat_path = sys.argv[2] if len(sys.argv) > 2 else None

    breeds = [l.strip() for l in open(pre + '.breeds')]
    fam = [l.split() for l in open(pre + '.fam')]
    truth = [f[0] for f in fam]
    ids = [f[1] for f in fam]
    held = set(l.strip() for l in open(pre + '.holdout') if l.strip())
    n_ind = len(ids)
    n_snp = sum(1 for _ in open(pre + '.bim'))
    bidx = {b: j for j, b in enumerate(breeds)}

    # training frequencies straight out of the file SCOPE was given
    P = np.zeros((n_snp, len(breeds)), dtype=np.float32)
    with open(pre + '.frq.strat') as fh:
        fh.readline()
        i = -1
        seen = 0
        for line in fh:
            f = line.split()
            j = bidx.get(f[2])
            if j is None:
                continue
            if seen % len(breeds) == 0:
                i += 1
            P[i, j] = float(f[5])
            seen += 1
    print(f"panel {n_snp} SNPs x {len(breeds)} breeds; held out {len(held)} dogs")

    raw = np.fromfile(pre + '.bed', dtype=np.uint8)
    bps = (n_ind + 3) // 4
    data = raw[3:].reshape(n_snp, bps)
    codes = np.empty((n_snp, bps * 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k::4] = (data >> (2 * k)) & 3
    dos = np.array([2, -1, 1, 0], dtype=np.int8)[codes[:, :n_ind]]

    test = [j for j in range(n_ind) if ids[j] in held]
    K = len(breeds)
    G = ((P.T @ P) + np.ones((K, K), dtype=np.float32)).astype(np.float64)
    G[np.diag_indices_from(G)] += 1e-6
    R = np.linalg.cholesky(G).T
    exp = (2.0 * P.mean(axis=1)).astype(np.float32)

    nnls_top = {}
    for n, j in enumerate(test):
        X = dos[:, j].astype(np.float32)
        m = X < 0
        if m.any():
            X[m] = exp[m]
        y = solve_triangular(R.T, ((P.T @ X) + 1.0).astype(np.float64), lower=True)
        q, _ = nnls(R, y)
        nnls_top[j] = breeds[int(np.argmax(q))]
        if n % 100 == 0:
            print(f"   NNLS {n}/{len(test)}", end='\r', flush=True)
    print(" " * 30, end='\r')

    scope_top = {}
    if qhat_path:
        Q = np.loadtxt(qhat_path)
        if Q.shape[0] != len(breeds):
            Q = Q.T
        for j in test:
            scope_top[j] = breeds[int(np.argmax(Q[:, j]))]

    cnt = Counter(truth)
    def report(name, pred):
        per = defaultdict(lambda: [0, 0])
        conf = defaultdict(Counter)
        for j in test:
            t = truth[j]
            per[t][1] += 1
            if pred[j] == t:
                per[t][0] += 1
            else:
                conf[t][pred[j]] += 1
        ok = sum(v[0] for v in per.values())
        tot = sum(v[1] for v in per.values())
        print(f"\n=== {name} ===")
        print(f"held-out dogs {tot}   correctly assigned {ok}  ({100*ok/tot:.2f}%)")
        for lo, hi, lab in ((5, 7, 'n=5-7'), (8, 19, 'n=8-19'), (20, 10**9, 'n>=20')):
            grp = [b for b in per if lo <= cnt[b] <= hi]
            c = sum(per[b][0] for b in grp); t = sum(per[b][1] for b in grp)
            if t:
                print(f"   {lab:<7} {100*c/t:5.1f}%  over {t:>4} dogs, {len(grp):>3} breeds")
        vill = sum(k for b in conf for g, k in conf[b].items()
                   if b.startswith('VILL') and g.startswith('VILL'))
        clup = sum(k for b in conf for g, k in conf[b].items()
                   if b.startswith('CLUP') and g.startswith('CLUP'))
        terr = sum(k for b in conf for g, k in conf[b].items())
        if terr:
            print(f"   of {terr} errors: {vill} village->village, {clup} wolf->wolf, "
                  f"{terr-vill-clup} other")
        return per, conf

    pn, cn = report("NNLS projection (what Stage 9 does today)", nnls_top)
    if scope_top:
        ps, cs = report("SCOPE Qhat (supervised)", scope_top)
        agree = sum(nnls_top[j] == scope_top[j] for j in test)
        both = sum(nnls_top[j] == truth[j] and scope_top[j] == truth[j] for j in test)
        only_n = sum(nnls_top[j] == truth[j] and scope_top[j] != truth[j] for j in test)
        only_s = sum(scope_top[j] == truth[j] and nnls_top[j] != truth[j] for j in test)
        print(f"\n=== head to head on {len(test)} held-out dogs ===")
        print(f"   the two agree      : {agree} ({100*agree/len(test):.1f}%)")
        print(f"   both correct       : {both}")
        print(f"   only NNLS correct  : {only_n}")
        print(f"   only SCOPE correct : {only_s}")
        print(f"   neither correct    : {len(test)-both-only_n-only_s}")


if __name__ == '__main__':
    main()
