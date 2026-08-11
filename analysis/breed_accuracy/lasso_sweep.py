#!/usr/bin/env python3
"""
Non-negative lasso for breed proportions — does L1 clean the tail?

    python3 lasso_sweep.py PROD_PREFIX BENCH_PREFIX BCF [BCF ...]

Ridge failed on mixes because L2 spreads mass across correlated predictors, and
the tail we want gone IS correlated-breed leakage. L1 is the opposite: under
q >= 0 it drives small coefficients to exactly zero.

An earlier note in this repo claimed L1 was degenerate here because the simplex
fixes sum(q). That was wrong. The sum-to-1 row is inert against ~131k data rows,
so the fit picks its own scale (~2, since P is a frequency and x a dosage) and
an L1 penalty genuinely bites.

It costs nothing to add, because with q >= 0 the penalty is linear:

    ||Aq - b||^2 + lam * 1'q  =  q'Gq - 2(A'b - lam/2 * 1)'q + const

so it is the same Gram matrix with a shifted right-hand side. lam is quoted
relative to mean(diag(G)) so it means the same thing at any SNP count.

Three metrics, because they can disagree:
  * proportion error on synthetic crosses, where truth is known  <- decisive
  * top-1 accuracy on held-out purebreds, to catch collateral damage
  * tail composition of two real mixed dogs, which is what a customer sees
"""
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import nnls

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_real_dog import load_panel, dosages_from_bcf          # noqa: E402
from synthetic_hybrids import load as load_bench, PAIRS, TRIPLES  # noqa: E402

LAMBDAS = [0.0, 1e-2, 3e-2, 1e-1, 3e-1, 1.0, 3.0, 10.0, 30.0]


def make_solver(P):
    K = P.shape[1]
    PtP = (P.T @ P).astype(np.float64)
    scale = float(np.mean(np.diag(PtP)))
    G = PtP + np.ones((K, K)) + 1e-6 * np.eye(K)
    R = np.linalg.cholesky(G).T

    def solve(x, lam):
        rhs = (P.T @ x).astype(np.float64) + 1.0
        if lam:
            rhs = rhs - 0.5 * lam * scale          # the L1 penalty, exactly
        y = solve_triangular(R.T, rhs, lower=True)
        q, _ = nnls(R, y)
        s = q.sum()
        return q / (s + 1e-12), s
    return solve


def main():
    prod, bench = sys.argv[1], sys.argv[2]
    bcfs = [a for a in sys.argv[3:] if a.endswith('.bcf')]

    # ── synthetic crosses: known truth ───────────────────────────
    b2, truth2, ids2, held2, P2, dos2 = load_bench(bench)
    solve2 = make_solver(P2)
    bidx = {b: j for j, b in enumerate(b2)}
    exp2 = (2.0 * P2.mean(axis=1)).astype(np.float32)
    by = defaultdict(list)
    for j, s in enumerate(ids2):
        if s in held2:
            by[truth2[j]].append(j)
    rng = np.random.default_rng(15052011)

    def clean(j):
        v = dos2[:, j].astype(np.float32)
        m = v < 0
        if m.any():
            v[m] = exp2[m]
        return v

    def gam(g):
        return (rng.random(g.shape) < (g / 2.0)).astype(np.float32)

    crosses = [([a, b], [.5, .5]) for a, b in PAIRS if by.get(a) and by.get(b)]
    crosses += [([a, b, c], [.5, .25, .25]) for a, b, c in TRIPLES
                if by.get(a) and by.get(b) and by.get(c)]
    samples = []
    for comp, e in crosses:
        for _ in range(10):
            d = [clean(rng.choice(by[b])) for b in comp]
            x = (gam(d[0]) + gam(d[1])) if len(comp) == 2 else \
                (gam(d[0]) + gam(gam(d[1]) + gam(d[2])))
            samples.append((comp, e, x))

    # ── held-out purebreds: collateral damage check ──────────────
    pure = [(j, truth2[j]) for j, s in enumerate(ids2) if s in held2]

    print(f"synthetic crosses: {len(samples)} dogs from {len(crosses)} crosses")
    print(f"held-out purebreds: {len(pure)}\n")
    print(f"  {'lambda':>7} | {'mix err':>8} {'comps found':>12} | "
          f"{'purebred top-1':>14} | {'mean #comp>=1%':>14}")
    print("  " + "-" * 66)
    for lam in LAMBDAS:
        errs, hits = [], 0
        ncomp = []
        for comp, e, x in samples:
            q, _ = solve2(x, lam)
            errs += [abs(q[bidx[b]] - t) for b, t in zip(comp, e)]
            names = [b2[i] for i in np.argsort(-q)[:len(comp)]]
            hits += all(b in names for b in comp)
            ncomp.append(int((q >= 0.01).sum()))
        ok = 0
        for j, t in pure:
            v = dos2[:, j].astype(np.float32)
            m = v < 0
            if m.any():
                v[m] = exp2[m]
            q, _ = solve2(v, lam)
            ok += (b2[int(np.argmax(q))] == t)
        print(f"  {lam:>7g} | {100*np.mean(errs):>7.2f}pp {hits:>6}/{len(samples):<5} | "
              f"{100*ok/len(pure):>13.2f}% | {np.mean(ncomp):>14.1f}")

    # ── real mixed dogs ──────────────────────────────────────────
    if not bcfs:
        return
    breeds, bim, P = load_panel(prod)
    solve1 = make_solver(P)
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

    for name, x in dogs.items():
        print(f"\n{'='*66}\n{name}\n{'='*66}")
        print(f"  {'lambda':>7} | {'poodle':>7} {'labrador':>8} | "
              f"{'#>=1%':>6} {'tail<5%':>8} | top non-poodle")
        for lam in LAMBDAS:
            q, _ = solve1(x, lam)
            pood = sum(q[i] for i in range(len(breeds)) if 'POODLE' in breeds[i])
            lab = q[breeds.index('LABRADOR_RETRIEVER')] if 'LABRADOR_RETRIEVER' in breeds else 0
            tail = sum(v for v in q if v < 0.05)
            o = [i for i in np.argsort(-q) if 'POODLE' not in breeds[i]][:2]
            top = ", ".join(f"{breeds[i].replace('_',' ').title()} {100*q[i]:.1f}%" for i in o)
            print(f"  {lam:>7g} | {100*pood:>6.1f}% {100*lab:>7.1f}% | "
                  f"{int((q>=0.01).sum()):>6} {100*tail:>7.1f}% | {top}")


if __name__ == '__main__':
    main()
