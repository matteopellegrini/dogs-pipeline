#!/usr/bin/env python3
"""
Project every Dog10K sample onto the Parker panel and check its breed call.

    # extract Parker sites from the Dog10K panel first
    awk 'BEGIN{OFS="\\t"} {c=$1; if(c !~ /^chr/) c="chr"c; print c,$4-1,$4}' parker.bim \\
      | sort -k1,1 -k2,2n > sites.bed
    bcftools query -R sites.bed -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' dog10k.bcf > gt.txt

    python3 project_dog10k.py parker.bim PHAT.txt CLUST.txt NAMES.json gt.txt SAMPLES.txt

Why this is a better test than leave-one-out
--------------------------------------------
Dog10K samples were never used to build Parker's Phat, so this is genuinely
held out rather than merely leave-one-out. Sample IDs carry the breed code as a
prefix (ACKR000001), so labels come free.

It also settles code-identity questions that leave-one-out cannot. If Dog10K's
SSHP dogs project onto Parker's SSHP, the code means the same thing in both
panels; if they land somewhere else, one of the two labellings is wrong. That
matters because PARKER_NAMES holds 236 codes for a 177-population panel and at
least one panel code was given the wrong name from that superset (NELK).

Applies the same panel modifications as Stage 9 so the test matches production.
"""
import json
import re
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import nnls

SALUKI_MERGE = {'SALU_ArabPen': 'SALU', 'SALU_CentAsia': 'SALU', 'SALU_Tribal': 'SALU'}
WOLF_MERGE = {f'WOLF-{r}': 'WOLF' for r in
              ('China', 'Croatia', 'India', 'Israel', 'Italy', 'Portugal', 'Yellowstone')}


def main():
    bim_path, phat_path, clust_path, names_path, gt_path, samples_path = sys.argv[1:7]
    names = json.load(open(names_path))

    # ── Parker SNPs, in Phat row order ───────────────────────────
    pos_index, a1, a2 = {}, [], []
    for i, line in enumerate(open(bim_path)):
        f = line.split()
        c = f[0] if f[0].startswith('chr') else 'chr' + f[0]
        pos_index[(c, f[3])] = i
        a1.append(f[4])
        a2.append(f[5])
    n_snp = len(a1)

    labels = []
    counts = Counter()
    for line in open(clust_path):
        p = line.split()
        if len(p) >= 3:
            counts[p[2]] += 1
            if p[2] not in labels:
                labels.append(p[2])
    P = np.loadtxt(phat_path, dtype=np.float32)
    assert P.shape == (n_snp, len(labels)), f'{P.shape} vs {(n_snp, len(labels))}'

    # ── same merges Stage 9 applies ──────────────────────────────
    MERGE = {}
    if '--no-merge' not in sys.argv:
        MERGE = dict(SALUKI_MERGE)
        MERGE.update(WOLF_MERGE)
    tgt = [MERGE.get(b, b) for b in labels] if MERGE else list(labels)
    new = []
    for t in tgt:
        if t not in new:
            new.append(t)
    M = np.zeros((n_snp, len(new)), dtype=np.float32)
    for j, t in enumerate(new):
        srcs = [i for i, x in enumerate(tgt) if x == t]
        w = np.array([float(counts[labels[i]]) for i in srcs])
        M[:, j] = (P[:, srcs] * w).sum(axis=1) / w.sum()
    P, labels = M, new
    K = len(labels)
    print(f"panel: {n_snp} SNPs x {K} populations (after Saluki and wolf merges)")

    samples = [s.strip() for s in open(samples_path) if s.strip()]
    n_ind = len(samples)

    # ── read genotypes into a dosage-of-a1 matrix ────────────────
    dos = np.full((n_snp, n_ind), -1, dtype=np.int8)
    seen = skipped = 0
    for line in open(gt_path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 4 + n_ind:
            continue
        idx = pos_index.get((f[0], f[1]))
        if idx is None:
            continue
        ref, alt = f[2], f[3]
        if ref not in (a1[idx], a2[idx]) or alt not in (a1[idx], a2[idx]):
            skipped += 1
            continue
        ref_is_a1 = (ref == a1[idx])
        row = dos[idx]
        for j, g in enumerate(f[4:4 + n_ind]):
            if len(g) < 3 or g[0] not in '01' or g[2] not in '01':
                continue
            n_alt = (g[0] == '1') + (g[2] == '1')
            row[j] = (2 - n_alt) if ref_is_a1 else n_alt
        seen += 1
    print(f"sites usable: {seen}  allele-mismatched: {skipped}  "
          f"missing genotypes: {100*(dos < 0).mean():.2f}%")

    exp = (2.0 * P.mean(axis=1)).astype(np.float32)
    G = ((P.T @ P) + np.ones((K, K), dtype=np.float32)).astype(np.float64)
    G[np.diag_indices_from(G)] += 1e-6
    R = np.linalg.cholesky(G).T

    top = []
    for s in range(0, n_ind, 200):
        e = min(s + 200, n_ind)
        X = dos[:, s:e].astype(np.float32)
        m = dos[:, s:e] < 0
        if m.any():
            X[m] = np.broadcast_to(exp[:, None], X.shape)[m]
        B = (P.T @ X) + 1.0
        for j in range(B.shape[1]):
            y = solve_triangular(R.T, B[:, j].astype(np.float64), lower=True)
            q, _ = nnls(R, y)
            q = q / (q.sum() + 1e-12)
            o = np.argsort(-q)[:3]
            top.append([(labels[k], float(q[k])) for k in o])
        print(f"  projected {e}/{n_ind}", end='\r', flush=True)
    print(" " * 30, end='\r')

    # ── score against the breed code in the sample ID ────────────
    code_of = []
    for s in samples:
        m = re.match(r'^([A-Za-z\-]+?)\d+$', s)
        code_of.append(m.group(1) if m else s)
    inpanel = set(labels)
    per = defaultdict(lambda: [0, 0])
    conf = defaultdict(Counter)
    shared = 0
    for s, c, t in zip(samples, code_of, top):
        if c not in inpanel:
            continue
        shared += 1
        per[c][1] += 1
        if t[0][0] == c:
            per[c][0] += 1
        else:
            conf[c][t[0][0]] += 1
    ok = sum(v[0] for v in per.values())
    print(f"\nDog10K samples whose breed code exists in the Parker panel: {shared} "
          f"({len(per)} breeds)")
    print(f"top-1 matches the sample's own code: {ok}/{shared} ({100*ok/shared:.1f}%)\n")

    print("breeds where Dog10K disagrees most (n>=3):")
    bad = sorted(((v[0] / v[1], v[1], b) for b, v in per.items() if v[1] >= 3))
    for frac, n, b in bad[:20]:
        got = ', '.join(f"{g}({k})" for g, k in conf[b].most_common(2))
        print(f"   {b:<8} {names.get(b,b)[:26]:<27} {int(frac*n)}/{n} -> {got}")

    for q in ('SSHP', 'NELK', 'COLL', 'SMCD', 'GDJK'):
        rows = [(s, t) for s, c, t in zip(samples, code_of, top) if c == q]
        if rows:
            print(f"\n--- Dog10K {q} samples (n={len(rows)}) ---")
            for s, t in rows[:8]:
                print("   " + s + "  ->  " +
                      ", ".join(f"{b} {v:.2f} [{names.get(b,b)[:22]}]" for b, v in t))


if __name__ == '__main__':
    main()
