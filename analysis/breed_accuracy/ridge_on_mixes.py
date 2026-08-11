#!/usr/bin/env python3
"""
Does ridge regularization help on MIXED dogs?

    python3 ridge_on_mixes.py PROD_PREFIX BENCH_PREFIX BCF [BCF ...]

The earlier sweep (regularization_sweep.py) measured top-1 accuracy on
purebreds and found ridge exactly neutral. That metric is blind to the tail:
a dog can be called correctly while carrying a fringe of spurious 2-3%
components, which is precisely what a mixed-breed customer report shows.

So this asks two things instead:

  1. what happens to real mixed dogs (Cosmo, Kiki) as lambda increases, and
  2. whether proportion error against KNOWN truth improves, using the synthetic
     crosses where the true composition is 50/50 or 50/25/25.

(2) is the one that can actually decide it. (1) alone would only show that the
numbers move, not that they move the right way.

A prediction worth stating before looking: L2 shrinkage spreads mass across
correlated predictors rather than concentrating it, so with 260 partly
collinear breeds ridge may well make the tail worse, not better. Sparsity is
what would clean a tail, and NNLS's non-negativity already supplies some.
"""
import os
import sys

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import nnls

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_real_dog import load_panel, dosages_from_bcf          # noqa: E402
from synthetic_hybrids import load as load_bench, PAIRS, TRIPLES  # noqa: E402

LAMBDAS = [0.0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]


def solver(P, lam):
    K = P.shape[1]
    G = (P.T @ P).astype(np.float64) + np.ones((K, K))
    if lam:
        G = G + lam * float(np.mean(np.diag(P.T @ P))) * np.eye(K)
    G = G + 1e-6 * np.eye(K)
    R = np.linalg.cholesky(G).T

    def solve(x):
        y = solve_triangular(R.T, ((P.T @ x) + 1.0).astype(np.float64), lower=True)
        q, _ = nnls(R, y)
        return q / (q.sum() + 1e-12)
    return solve


def main():
    prod, bench = sys.argv[1], sys.argv[2]
    bcfs = [a for a in sys.argv[3:] if a.endswith('.bcf')]

    # ── part 1: real mixed dogs ──────────────────────────────────
    breeds, bim, P = load_panel(prod)
    exp = (2.0 * P.mean(axis=1)).astype(np.float32)
    import tempfile
    fd, bed = tempfile.mkstemp(suffix='.bed')
    with os.fdopen(fd, 'w') as fh:
        for c, pos in sorted(((b[0] if b[0].startswith('chr') else 'chr' + b[0]),
                              int(b[3])) for b in bim):
            fh.write(f"{c}\t{pos-1}\t{pos}\n")
    dogs = {}
    for b in bcfs:
        x = dosages_from_bcf(b, bim, bed)
        m = ~np.isfinite(x)
        x[m] = exp[m]
        dogs[os.path.basename(b).replace('_imputed_dog10k.bcf', '')] = x
    os.unlink(bed)

    for lam in LAMBDAS:
        solve = solver(P, lam)
        print(f"\n{'='*58}\nlambda = {lam:g}\n{'='*58}")
        for name, x in dogs.items():
            q = solve(x)
            o = np.argsort(-q)
            pood = sum(q[i] for i in range(len(breeds)) if 'POODLE' in breeds[i])
            tail = sum(q[i] for i in o[3:] if q[i] < 0.05)
            print(f"  {name}:  poodle-total {100*pood:.1f}%   "
                  f"top: " + ", ".join(f"{breeds[i].replace('_',' ').title()} "
                                       f"{100*q[i]:.1f}%" for i in o[:4]))
            print(f"        components >=1%: {int((q >= 0.01).sum())}   "
                  f"mass in sub-5% tail: {100*tail:.1f}%")

    # ── part 2: synthetic crosses, where truth is known ──────────
    print(f"\n\n{'='*58}\nproportion error on synthetic crosses (ground truth)\n{'='*58}")
    b2, truth2, ids2, held2, P2, dos2 = load_bench(bench)
    from collections import defaultdict
    by = defaultdict(list)
    for j, s in enumerate(ids2):
        if s in held2:
            by[truth2[j]].append(j)
    bidx = {b: j for j, b in enumerate(b2)}
    exp2 = (2.0 * P2.mean(axis=1)).astype(np.float32)
    rng = np.random.default_rng(15052011)

    def clean(j):
        v = dos2[:, j].astype(np.float32)
        m = v < 0
        if m.any():
            v[m] = exp2[m]
        return v

    def gamete(g):
        return (rng.random(g.shape) < (g / 2.0)).astype(np.float32)

    crosses = []
    for a, b in PAIRS:
        if by.get(a) and by.get(b):
            crosses.append(([a, b], [0.5, 0.5]))
    for a, b, c in TRIPLES:
        if by.get(a) and by.get(b) and by.get(c):
            crosses.append(([a, b, c], [0.5, 0.25, 0.25]))

    REPS = 10
    samples = []
    for comp, expct in crosses:
        for _ in range(REPS):
            d = [clean(rng.choice(by[b])) for b in comp]
            if len(comp) == 2:
                x = gamete(d[0]) + gamete(d[1])
            else:
                x = gamete(d[0]) + gamete(gamete(d[1]) + gamete(d[2]))
            samples.append((comp, expct, x))

    print(f"  {len(samples)} simulated dogs from {len(crosses)} crosses\n")
    print(f"  {'lambda':>8}  {'mean abs err':>13}  {'all comps in top-k':>19}  "
          f"{'spurious >=5%':>14}")
    for lam in LAMBDAS:
        solve = solver(P2, lam)
        errs, hits, spur = [], 0, []
        for comp, expct, x in samples:
            q = solve(x)
            errs += [abs(q[bidx[b]] - e) for b, e in zip(comp, expct)]
            names = [b2[i] for i in np.argsort(-q)[:len(comp)]]
            hits += all(b in names for b in comp)
            idx = {bidx[b] for b in comp}
            spur.append(sum(q[i] for i in range(len(b2))
                            if i not in idx and q[i] >= 0.05))
        print(f"  {lam:>8g}  {100*np.mean(errs):>12.2f}pp  "
              f"{hits:>10}/{len(samples):<8}  {100*np.mean(spur):>13.2f}%")


if __name__ == '__main__':
    main()
