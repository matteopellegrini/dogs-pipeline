#!/usr/bin/env python3
"""
Should village dogs be dropped from the panel, since they are not breeds?

    python3 village_ablation.py BENCH_PREFIX [--min-n 7] [--lasso 0.3]

The case for dropping them is that "Village Dog East Asia 2.9%" means nothing to
a pet owner, and that they act as a residual sink. The case against is that the
residual has to go somewhere: remove the sink and that mass lands on whichever
BREED fits least badly, which is a confident wrong answer rather than a vague
one. Imported rescues are also a large share of the US pet population, so some
customers genuinely carry this ancestry.

Measured three ways — effect on breed dogs (synthetic crosses and held-out
purebreds), and effect on village dogs themselves.

Result: dropping them changes breed-dog accuracy by nothing at all, and every
held-out village dog is then mis-called as a basal breed (Shar-Pei, Basenji,
Thai Ridgeback, Chow Chow) with no signal that anything is wrong. So the label
is a DISPLAY problem, not a reason to remove the population.
"""
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import nnls

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from synthetic_hybrids import load, PAIRS, TRIPLES   # noqa: E402


def main():
    pre = sys.argv[1]
    a = sys.argv[2:]
    minn = int(a[a.index('--min-n') + 1]) if '--min-n' in a else 7
    lam = float(a[a.index('--lasso') + 1]) if '--lasso' in a else 0.3

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

    crosses = [([x, y], [.5, .5]) for x, y in PAIRS if by.get(x) and by.get(y)]
    crosses += [([x, y, z], [.5, .25, .25]) for x, y, z in TRIPLES
                if all(by.get(w) for w in (x, y, z))]
    samples = []
    for comp, e in crosses:
        for _ in range(10):
            d = [clean(rng.choice(by[b])) for b in comp]
            samples.append((comp, e,
                            (gam(d[0]) + gam(d[1])) if len(comp) == 2
                            else (gam(d[0]) + gam(gam(d[1]) + gam(d[2])))))
    pure = [(j, truth[j]) for j, s in enumerate(ids)
            if s in held and not truth[j].startswith('VILLAGE')]
    vill = [(j, truth[j]) for j, s in enumerate(ids)
            if s in held and truth[j].startswith('VILLAGE')]

    for label, drop in (("village dogs KEPT", False), ("village dogs REMOVED", True)):
        keep = [j for j, b in enumerate(breeds)
                if n[b] >= minn and not (drop and b.startswith('VILLAGE'))]
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
            q = solve(x)
            errs += [abs(q[bidx[c]] - t) if c in bidx else t
                     for c, t in zip(comp, e)]
            idx = {bidx[c] for c in comp if c in bidx}
            spur.append(sum(q[i] for i in range(K) if i not in idx))
        ok = sum(bk[int(np.argmax(solve(clean(j))))] == t for j, t in pure)
        vok, vcalls = 0, Counter()
        for j, t in vill:
            top = bk[int(np.argmax(solve(clean(j))))]
            vok += (top == t)
            vcalls[top] += 1
        print(f"\n=== {label} ({K} breeds) ===")
        print(f"  synthetic crosses : err {100*np.mean(errs):.2f}pp   "
              f"spurious {100*np.mean(spur):.1f}%")
        print(f"  held-out purebreds: {100*ok/len(pure):.2f}%  ({ok}/{len(pure)})")
        print(f"  held-out VILLAGE dogs ({len(vill)}): correctly called {vok}")
        print("     most common calls: " + ", ".join(
            f"{b.replace('_',' ').title()} ({k})" for b, k in vcalls.most_common(4)))


if __name__ == '__main__':
    main()
