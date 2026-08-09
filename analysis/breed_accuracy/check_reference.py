#!/usr/bin/env python3
"""
Leave-one-out self-check of the Parker reference panel.

    python3 check_reference.py PREFIX PHAT.txt CLUST.txt NAMES.json [flags]

    --loo            leave-one-out (otherwise in-sample, which is meaningless)
    --no-wolves      drop the 7 n=1 wolf populations
    --merge-wolves   pool the 7 wolf populations into one WOLF (n=7)
    --merge-saluki   pool SALU_ArabPen/CentAsia/Tribal into SALU
    --akc-display    score on the AKC display name, not the raw population

Runs the same projection Stage 9 runs for a customer dog — NNLS of dosage-of-a1
against Phat with a sum-to-1 row — on the reference individuals themselves.

In-sample is circular and always ~100%: Phat is the per-breed empirical allele
frequency (verified, corr 0.9997+ against the genotypes), so every dog helped
define its own target and an n=1 population recovers itself by construction.

Leave-one-out is the real test, and Phat being an empirical mean makes it cheap:
a dog's own contribution comes out analytically as
    p_loo = (n*p - x/2)/(n-1)   =>   delta = (p - x/2)/(n-1)
Only that one column changes, so the Gram matrix takes a rank-2 update instead
of re-running SCOPE once per dog.

The flags separate effects that are easy to conflate. --no-wolves is a pure
model change. --merge-saluki also relabels 9 dogs, so some of its gain is the
taxonomy becoming correct rather than the model discriminating better.
--akc-display changes only the scoring, not the model.
"""
import json
import sys
from collections import defaultdict

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import nnls

CHUNK = 200
SALUKI_MERGE = {'SALU_ArabPen': 'SALU', 'SALU_CentAsia': 'SALU', 'SALU_Tribal': 'SALU'}
WOLF_MERGE = {f'WOLF-{r}': 'WOLF' for r in
              ('China', 'Croatia', 'India', 'Israel', 'Italy', 'Portugal', 'Yellowstone')}
AKC_DISPLAY = {'SPOO': 'Poodle', 'MPOO': 'Poodle', 'TPOO': 'Poodle',
               'COLL': 'Collie', 'SSHP': 'Collie',
               'XOLO': 'Xoloitzcuintli', 'MXOL': 'Xoloitzcuintli'}


def read_bed(prefix, n_ind, n_snp):
    raw = np.fromfile(prefix + '.bed', dtype=np.uint8)
    assert raw[0] == 108 and raw[1] == 27 and raw[2] == 1
    bps = (n_ind + 3) // 4
    data = raw[3:].reshape(n_snp, bps)
    codes = np.empty((n_snp, bps * 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k::4] = (data >> (2 * k)) & 3
    return np.array([2, -1, 1, 0], dtype=np.int8)[codes[:, :n_ind]]


def main():
    prefix, phat_path, clust_path, names_path = sys.argv[1:5]
    flags = set(a for a in sys.argv[5:])
    loo = '--loo' in flags
    names = json.load(open(names_path))

    fam = [l.split() for l in open(prefix + '.fam')]
    ids = [f[1] for f in fam]
    n_ind, n_snp = len(fam), sum(1 for _ in open(prefix + '.bim'))

    truth, labels = {}, []
    for line in open(clust_path):
        p = line.split()
        if len(p) >= 3:
            truth[p[1]] = p[2]
            if p[2] not in labels:
                labels.append(p[2])

    P = np.loadtxt(phat_path, dtype=np.float32)
    assert P.shape == (n_snp, len(labels))
    dos = read_bed(prefix, n_ind, n_snp)
    miss = dos < 0
    counts = defaultdict(int)
    for b in truth.values():
        counts[b] += 1

    # ── panel modifications, applied exactly as the pipeline does ────
    if '--no-wolves' in flags:
        keep = [i for i, b in enumerate(labels) if not b.upper().startswith('WOLF')]
        dropped = len(labels) - len(keep)
        P = P[:, keep]
        labels = [labels[i] for i in keep]
        print(f"dropped {dropped} wolf populations")

    MERGE = {}
    if '--merge-saluki' in flags:
        MERGE.update(SALUKI_MERGE)
    if '--merge-wolves' in flags:
        MERGE.update(WOLF_MERGE)
    if MERGE:
        tgt = [MERGE.get(b, b) for b in labels]
        new = []
        for t in tgt:
            if t not in new:
                new.append(t)
        M = np.zeros((P.shape[0], len(new)), dtype=P.dtype)
        for j, t in enumerate(new):
            srcs = [i for i, x in enumerate(tgt) if x == t]
            w = np.array([float(counts[labels[i]]) for i in srcs])
            M[:, j] = (P[:, srcs] * w).sum(axis=1) / w.sum()
        P, labels = M, new
        truth = {d: MERGE.get(b, b) for d, b in truth.items()}
        counts = defaultdict(int)
        for b in truth.values():
            counts[b] += 1
        for t in sorted(set(MERGE.values())):
            print(f"merged -> {t} (n={counts[t]})")

    idx_of = {b: i for i, b in enumerate(labels)}
    K = len(labels)
    exp = (2.0 * P.mean(axis=1)).astype(np.float32)
    G = ((P.T @ P) + np.ones((K, K), dtype=np.float32)).astype(np.float64)
    G[np.diag_indices_from(G)] += 1e-6

    top1 = [None] * n_ind
    wolf_frac = {}
    for s in range(0, n_ind, CHUNK):
        e = min(s + CHUNK, n_ind)
        X = dos[:, s:e].astype(np.float32)
        m = miss[:, s:e]
        if m.any():
            X[m] = np.broadcast_to(exp[:, None], X.shape)[m]
        B = (P.T @ X) + 1.0
        for j in range(B.shape[1]):
            i = s + j
            lab = truth.get(ids[i])
            if lab is None or lab not in idx_of:
                continue                       # e.g. a wolf after --no-wolves
            bb = B[:, j].astype(np.float64)
            Gd = G
            if loo:
                if counts[lab] <= 1:
                    continue                   # n=1: LOO removes the population
                col = idx_of[lab]
                x = X[:, j]
                delta = ((P[:, col] - x / 2.0) / (counts[lab] - 1)).astype(np.float32)
                u = (P.T @ delta).astype(np.float64)
                Gd = G.copy()
                Gd[col, :] += u
                Gd[:, col] += u
                Gd[col, col] += float(delta @ delta)
                bb[col] += float(delta @ x)
            try:
                Rd = np.linalg.cholesky(Gd).T
            except np.linalg.LinAlgError:
                continue
            y = solve_triangular(Rd.T, bb, lower=True)
            q, _ = nnls(Rd, y)
            top1[i] = labels[int(np.argmax(q))]
            qn = q / (q.sum() + 1e-12)
            if 'WOLF' in idx_of:
                wolf_frac[i] = float(qn[idx_of['WOLF']])

    def disp(code):
        if '--akc-display' in flags:
            return AKC_DISPLAY.get(code) or names.get(code, code)
        return code

    scored = [(i, ids[i]) for i in range(n_ind) if top1[i] is not None and ids[i] in truth]
    per = defaultdict(lambda: [0, 0])
    wrong = []
    for i, d in scored:
        t = truth[d]
        per[t][1] += 1
        if disp(top1[i]) == disp(t):
            per[t][0] += 1
        else:
            wrong.append((d, t, top1[i]))
    ok = sum(v[0] for v in per.values())
    tot = sum(v[1] for v in per.values())

    tag = ' '.join(sorted(flags - {'--loo'})) or '(baseline panel)'
    print(f"\n=== {tag} ===")
    print(f"scored {tot} dogs   correct {ok}  ({100*ok/tot:.2f}%)   errors {len(wrong)}")

    small = [b for b in per if per[b][1] <= 3]
    big = [b for b in per if per[b][1] >= 8]
    for name, grp in (('n<=3', small), ('n>=8', big)):
        c = sum(per[b][0] for b in grp); t = sum(per[b][1] for b in grp)
        if t:
            print(f"   {name}: {100*c/t:.0f}% over {t} dogs ({len(grp)} populations)")
    if wolf_frac:
        dogs = [v for i, v in wolf_frac.items() if truth.get(ids[i]) != 'WOLF']
        wolves = [v for i, v in wolf_frac.items() if truth.get(ids[i]) == 'WOLF']
        dogs.sort()
        if dogs:
            print(f"   wolf ancestry assigned to the {len(dogs)} NON-wolf dogs: "
                  f"median {100*dogs[len(dogs)//2]:.2f}%  "
                  f"95th {100*dogs[int(len(dogs)*0.95)]:.2f}%  max {100*dogs[-1]:.2f}%")
        top = sorted(((v, ids[i], truth.get(ids[i])) for i, v in wolf_frac.items()
                      if truth.get(ids[i]) != 'WOLF'), reverse=True)[:10]
        print("   highest wolf fraction among non-wolf dogs:")
        for v, d, t in top:
            print(f"     {100*v:5.1f}%  {d:<18} {t} ({names.get(t, t)})")
        if wolves:
            print(f"   wolf ancestry recovered in the {len(wolves)} wolves: "
                  + ", ".join(f"{100*v:.0f}%" for v in sorted(wolves, reverse=True)))
    if wrong:
        print("   remaining errors:")
        for d, t, g in wrong[:40]:
            print(f"     {d:<18} {t:<16} -> {g}")


if __name__ == '__main__':
    main()
