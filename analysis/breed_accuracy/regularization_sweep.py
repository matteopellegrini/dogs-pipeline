#!/usr/bin/env python3
"""
Does a regularized projection beat plain NNLS for breed assignment?

    python3 regularization_sweep.py BENCH_PREFIX

Scores the same held-out dogs under several penalties. All of them are cheap
because the fit goes through the normal equations: for A = [P ; w*1...1],
min ||Aq - b|| has G = P'P + w^2*J and rhs = P'x + w^2, so a penalty is just an
edit to G and the right-hand side, and each solve stays K x K.

Variants
--------
plain      what Stage 9 does today: NNLS with a single sum-to-1 row appended,
           then normalise. With 131,353 data rows against that one constraint
           row, the simplex is imposed after the fact rather than during it.
simplex-w  the same, but the constraint row weighted so it actually binds.
           Worth separating from regularization proper: it changes the feasible
           set, not the penalty.
ridge      + lambda*||q||^2, i.e. G + lambda*I. Shrinks toward zero and spreads
           mass across correlated breeds — which may help stability on
           near-identical pairs or may simply blur the very distinction that
           top-1 depends on.
shrink     + lambda*||q - uniform||^2, shrinking toward the panel average
           instead of toward zero.

lambda is quoted relative to mean(diag(G)) so the numbers mean the same thing
regardless of SNP count.
"""
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import nnls


def load(pre):
    breeds = [l.strip() for l in open(pre + '.breeds')]
    fam = [l.split() for l in open(pre + '.fam')]
    truth = [f[0] for f in fam]
    ids = [f[1] for f in fam]
    held = set(l.strip() for l in open(pre + '.holdout') if l.strip())
    n_ind = len(ids)
    n_snp = sum(1 for _ in open(pre + '.bim'))
    bidx = {b: j for j, b in enumerate(breeds)}

    P = np.zeros((n_snp, len(breeds)), dtype=np.float32)
    with open(pre + '.frq.strat') as fh:
        fh.readline()
        i, seen = -1, 0
        for line in fh:
            f = line.split()
            j = bidx.get(f[2])
            if j is None:
                continue
            if seen % len(breeds) == 0:
                i += 1
            P[i, j] = float(f[5])
            seen += 1

    raw = np.fromfile(pre + '.bed', dtype=np.uint8)
    bps = (n_ind + 3) // 4
    data = raw[3:].reshape(n_snp, bps)
    codes = np.empty((n_snp, bps * 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k::4] = (data >> (2 * k)) & 3
    dos = np.array([2, -1, 1, 0], dtype=np.int8)[codes[:, :n_ind]]
    return breeds, truth, ids, held, P, dos


def main():
    pre = sys.argv[1]
    breeds, truth, ids, held, P, dos = load(pre)
    K = len(breeds)
    test = [j for j in range(len(ids)) if ids[j] in held]
    exp = (2.0 * P.mean(axis=1)).astype(np.float32)
    PtP = (P.T @ P).astype(np.float64)
    scale = float(np.mean(np.diag(PtP)))
    print(f"{len(test)} held-out dogs, {K} breeds, mean diag(G) = {scale:.0f}\n")

    X = {}
    for j in test:
        v = dos[:, j].astype(np.float32)
        m = v < 0
        if m.any():
            v[m] = exp[m]
        X[j] = v

    def run(name, w=1.0, lam=0.0, toward=None):
        G = PtP + (w ** 2) * np.ones((K, K))
        if lam:
            G = G + lam * scale * np.eye(K)
        G = G + 1e-6 * np.eye(K)
        try:
            R = np.linalg.cholesky(G).T
        except np.linalg.LinAlgError:
            print(f"   {name:<28} not positive definite")
            return None
        ok = 0
        for j in test:
            b = (P.T @ X[j]).astype(np.float64) + (w ** 2)
            if lam and toward is not None:
                b = b + lam * scale * toward
            y = solve_triangular(R.T, b, lower=True)
            q, _ = nnls(R, y)
            if breeds[int(np.argmax(q))] == truth[j]:
                ok += 1
        print(f"   {name:<28} {ok:>4}/{len(test)}   {100*ok/len(test):6.2f}%")
        return ok

    uniform = np.full(K, 1.0 / K)
    print("baseline")
    base = run("plain (Stage 9 today)", w=1.0)

    print("\nsum-to-1 constraint weight")
    for w in (10.0, 100.0, 1000.0, 10000.0):
        run(f"simplex w={w:g}", w=w)

    print("\nridge, shrinking toward zero")
    for lam in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2):
        run(f"ridge lambda={lam:g}", w=1.0, lam=lam)

    print("\nshrinkage toward the panel average")
    for lam in (1e-5, 1e-4, 1e-3):
        run(f"shrink lambda={lam:g}", w=1.0, lam=lam, toward=uniform)

    print(f"\nbaseline was {base}/{len(test)} — anything below this is a regression")


if __name__ == '__main__':
    main()
