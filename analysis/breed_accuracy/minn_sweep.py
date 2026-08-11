#!/usr/bin/env python3
"""
How small can a reference population be before it does more harm than good?

    python3 minn_sweep.py BENCH_PREFIX [--lasso 0.3]

A thin reference has noisy allele frequencies, which by chance align with
residual variation in any dog — so it absorbs unexplained variance from
everyone. Matteo spotted this in the wild: Kiki, a Labradoodle, was picking up
3.8% Catalan Sheepdog (n=5), while every genuine component of her profile came
from a reference with n>=10.

The earlier held-out benchmark could not see this. It asked whether a thin breed
recognises ITSELF — n=5-7 scored 90.1%, which looked acceptable — and never
asked whether it wrongly absorbs ancestry from OTHER dogs. A reference can pass
the first test while failing the second, and the second is what a customer sees.

So this measures SPURIOUS MASS: on simulated crosses of known composition, the
total proportion assigned to breeds that are not parents. That number should be
zero and is not, and it is the quantity a minimum-n threshold trades against
breed coverage.

Caveat when reading the purebred column: it scores only breeds still in the
panel, so part of any rise is the denominator shedding hard cases rather than
the model improving. The 'breeds lost' column is the honest cost — dogs of those
breeds can no longer be called at all.
"""
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import nnls

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from synthetic_hybrids import load, PAIRS, TRIPLES   # noqa: E402

MIN_NS = [5, 6, 7, 8, 9, 10, 11, 12]
REPS = 10


def main():
    pre = sys.argv[1]
    args = sys.argv[2:]
    lam = float(args[args.index('--lasso') + 1]) if '--lasso' in args else 0.3

    breeds, truth, ids, held, P, dos = load(pre)
    n = Counter(truth)
    exp = (2.0 * P.mean(axis=1)).astype(np.float32)
    by = defaultdict(list)
    for j, s in enumerate(ids):
        if s in held:
            by[truth[j]].append(j)

    def clean(j):
        v = dos[:, j].astype(np.float32)
        m = v < 0
        if m.any():
            v[m] = exp[m]
        return v

    rng = np.random.default_rng(7)

    def gam(g):
        return (rng.random(g.shape) < (g / 2.0)).astype(np.float32)

    crosses = [([a, b], [.5, .5]) for a, b in PAIRS if by.get(a) and by.get(b)]
    crosses += [([a, b, c], [.5, .25, .25]) for a, b, c in TRIPLES
                if all(by.get(z) for z in (a, b, c))]
    samples = []
    for comp, e in crosses:
        for _ in range(REPS):
            d = [clean(rng.choice(by[b])) for b in comp]
            x = (gam(d[0]) + gam(d[1])) if len(comp) == 2 else \
                (gam(d[0]) + gam(gam(d[1]) + gam(d[2])))
            samples.append((comp, e, x))
    pure = [(j, truth[j]) for j, s in enumerate(ids) if s in held]

    print(f"{len(samples)} simulated crosses, {len(pure)} held-out purebreds, "
          f"lasso lambda={lam}\n")
    print(f"  {'min-n':>6} {'breeds':>7} | {'mix err':>8} {'SPURIOUS':>9} | "
          f"{'purebred':>9} {'lost':>6}")
    print("  " + "-" * 56)
    for minn in MIN_NS:
        keep = [j for j, b in enumerate(breeds) if n[b] >= minn]
        Pk = P[:, keep]
        bk = [breeds[j] for j in keep]
        K = len(keep)
        bidx = {b: i for i, b in enumerate(bk)}
        scale = float(np.mean(np.diag(Pk.T @ Pk)))
        G = (Pk.T @ Pk).astype(np.float64) + np.ones((K, K)) + 1e-6 * np.eye(K)
        R = np.linalg.cholesky(G).T

        def solve(x):
            rhs = (Pk.T @ x).astype(np.float64) + 1.0 - 0.5 * lam * scale
            q, _ = nnls(R, solve_triangular(R.T, rhs, lower=True))
            return q / (q.sum() + 1e-12)

        errs, spur = [], []
        for comp, e, x in samples:
            if not all(c in bidx for c in comp):
                continue
            q = solve(x)
            errs += [abs(q[bidx[c]] - t) for c, t in zip(comp, e)]
            idx = {bidx[c] for c in comp}
            spur.append(sum(q[i] for i in range(K) if i not in idx))
        ok = tot = 0
        for j, t in pure:
            if t not in bidx:
                continue
            ok += (bk[int(np.argmax(solve(clean(j))))] == t)
            tot += 1
        print(f"  {minn:>6} {K:>7} | {100*np.mean(errs):>7.2f}pp "
              f"{100*np.mean(spur):>8.1f}% | {100*ok/tot:>8.2f}% "
              f"{len(breeds)-K:>6}")


if __name__ == '__main__':
    main()
