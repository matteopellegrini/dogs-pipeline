#!/usr/bin/env python3
"""
Synthesise two- and three-breed crosses from held-out dogs and check whether the
panel recovers their true composition.

    python3 synthetic_hybrids.py BENCH_PREFIX [--reps 20] [--list]

Every benchmark so far has scored purebreds, but most customer dogs are mixes,
and a mix is a harder and different question: not "is the top call right" but
"are both parent breeds reported, in roughly the right proportions".

How the crosses are made
------------------------
Mendelian transmission, one allele drawn from each parent:

    gamete(g) ~ Bernoulli(g / 2)          g in {0,1,2} = dosage of a1
    F1(A, B)  = gamete(A) + gamete(B)                    50/50
    3-way     = gamete(A) + gamete(F1(B, C))             50/25/25

Drawing alleles independently per SNP ignores linkage, which would matter for
an LD-aware method but not here: the projection matches per-SNP allele
frequencies and never looks at haplotype structure, so the estimate depends only
on the genome-wide allele mixture. The 3-way is a real cross (A x an F1), which
is why it is 50/25/25 rather than a third each — an equal three-way mix is not
something a dog can be.

Both parents are drawn from the HELD-OUT dogs, whose genotypes were excluded
when the panel frequencies were computed. Otherwise a cross of two panel
members would partly be scoring the panel against itself.
"""
import re
import sys
from collections import defaultdict

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import nnls

# Crosses chosen to cover what customers actually send in, plus deliberate
# hard cases. Names are canonical labels from harmonize.py.
PAIRS = [
    ('LABRADOR_RETRIEVER', 'STANDARD_POODLE'),      # Labradoodle
    ('GOLDEN_RETRIEVER', 'STANDARD_POODLE'),        # Goldendoodle
    ('GERMAN_SHEPHERD', 'SIBERIAN_HUSKY'),          # like DOGS-Gen-6
    ('CHIHUAHUA', 'DACHSHUND'),
    ('BEAGLE', 'BOXER'),
    ('LABRADOR_RETRIEVER', 'GERMAN_SHEPHERD'),
    ('COCKER_SPANIEL', 'CAVALIER_KING_CHARLES_SPANIEL'),
    ('BORDER_COLLIE', 'AUSTRALIAN_SHEPHERD'),       # hard: closely related herders
    ('PARSON_RUSSELL_TERRIER', 'JACK_RUSSELL_TERRIER'),  # hard: near-identical
]
TRIPLES = [
    ('LABRADOR_RETRIEVER', 'STANDARD_POODLE', 'GOLDEN_RETRIEVER'),
    ('GERMAN_SHEPHERD', 'SIBERIAN_HUSKY', 'BEAGLE'),
    ('CHIHUAHUA', 'DACHSHUND', 'POMERANIAN'),
]


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
    args = sys.argv[2:]
    reps = int(args[args.index('--reps') + 1]) if '--reps' in args else 20
    breeds, truth, ids, held, P, dos = load(pre)
    K = len(breeds)
    rng = np.random.default_rng(15052011)

    by_breed = defaultdict(list)
    for j, s in enumerate(ids):
        if s in held:
            by_breed[truth[j]].append(j)

    if '--list' in args:
        for b, v in sorted(by_breed.items(), key=lambda kv: -len(kv[1]))[:40]:
            print(f"   {b:<34} {len(v)} held-out donors")
        return

    exp = (2.0 * P.mean(axis=1)).astype(np.float32)
    G = ((P.T @ P) + np.ones((K, K), dtype=np.float32)).astype(np.float64)
    G[np.diag_indices_from(G)] += 1e-6
    R = np.linalg.cholesky(G).T
    bidx = {b: j for j, b in enumerate(breeds)}

    def clean(j):
        v = dos[:, j].astype(np.float32)
        m = v < 0
        if m.any():
            v[m] = exp[m]
        return v

    def gamete(g):
        """one transmitted allele per SNP: Bernoulli(dosage/2)"""
        return (rng.random(g.shape) < (g / 2.0)).astype(np.float32)

    def project(x):
        y = solve_triangular(R.T, ((P.T @ x) + 1.0).astype(np.float64), lower=True)
        q, _ = nnls(R, y)
        return q / (q.sum() + 1e-12)

    def run(components, expected):
        """components: list of breed names; expected: list of true proportions"""
        missing = [b for b in components if b not in by_breed or not by_breed[b]]
        if missing:
            return None
        hits_all = hits_any = 0
        err = []
        top = defaultdict(int)
        for _ in range(reps):
            donors = [clean(rng.choice(by_breed[b])) for b in components]
            if len(components) == 2:
                x = gamete(donors[0]) + gamete(donors[1])
            else:
                inner = gamete(donors[1]) + gamete(donors[2])
                x = gamete(donors[0]) + gamete(inner)
            q = project(x)
            order = np.argsort(-q)
            names = [breeds[i] for i in order[:len(components)]]
            top[breeds[order[0]]] += 1
            if all(b in names for b in components):
                hits_all += 1
            if any(b in names for b in components):
                hits_any += 1
            err.append([abs(q[bidx[b]] - e) for b, e in zip(components, expected)])
        err = np.array(err)
        label = ' x '.join(b.replace('_', ' ').title() for b in components)
        exp_s = '/'.join(f"{100*e:.0f}" for e in expected)
        print(f"\n   {label}   (true {exp_s})")
        print(f"      all components in top-{len(components)} : "
              f"{hits_all}/{reps}    any: {hits_any}/{reps}")
        print(f"      mean abs error per component  : "
              + ", ".join(f"{100*m:.1f}pp" for m in err.mean(axis=0)))
        common = sorted(top.items(), key=lambda kv: -kv[1])[:2]
        print(f"      most frequent top-1           : "
              + ", ".join(f"{b.replace('_',' ').title()} ({k})" for b, k in common))
        return hits_all

    print(f"panel {K} breeds, {len(held)} held-out donors, {reps} replicates per cross")
    print("\n" + "=" * 62)
    print("TWO-BREED CROSSES (F1, true 50/50)")
    print("=" * 62)
    ok2 = tot2 = 0
    for a, b in PAIRS:
        r = run([a, b], [0.5, 0.5])
        if r is None:
            print(f"\n   {a} x {b}: skipped, no held-out donor")
        else:
            ok2 += r; tot2 += reps
    print("\n" + "=" * 62)
    print("THREE-BREED CROSSES (A x F1(B,C), true 50/25/25)")
    print("=" * 62)
    ok3 = tot3 = 0
    for a, b, c in TRIPLES:
        r = run([a, b, c], [0.5, 0.25, 0.25])
        if r is None:
            print(f"\n   {a} x {b} x {c}: skipped, no held-out donor")
        else:
            ok3 += r; tot3 += reps
    print("\n" + "=" * 62)
    if tot2:
        print(f"two-breed  : all components recovered in {ok2}/{tot2} ({100*ok2/tot2:.0f}%)")
    if tot3:
        print(f"three-breed: all components recovered in {ok3}/{tot3} ({100*ok3/tot3:.0f}%)")


if __name__ == '__main__':
    main()
