#!/usr/bin/env python3
"""
Build and evaluate a merged Parker + Dog10K breed reference.

    python3 merged_panel.py PARKER_PREFIX CLUST.txt DOG10K_GT.txt DOG10K_SAMPLES.txt \\
        [--panel parker|dog10k|merged] [--min-n 2] [--cross]

SCOPE itself is not needed. Its Phat is the per-breed empirical allele frequency
(verified against the genotypes at corr 0.9997+), so a panel for any set of
individuals is just their per-breed allele frequencies at the shared SNPs.

Evaluation is leave-one-out: a dog's own contribution is removed analytically,
p_loo = (n*p - x/2)/(n-1), with a rank-2 Gram update rather than refitting.

--cross runs the fully held-out test instead: build the panel from ONE source
and predict the OTHER source's dogs. Nothing is shared between fit and test, so
it cannot flatter itself the way leave-one-out can — leave-one-out is what hid
the Saluki regression earlier.

The platform question this is meant to answer: Parker is array genotypes at SNPs
ascertained to discriminate breeds, Dog10K is unascertained WGS. If per-breed
allele frequencies from the two sources disagree systematically, pooling them
injects batch structure that the model would learn as "breed". The
--panel merged run reports that agreement directly for the 119 shared breeds.
"""
import re
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import nnls


def read_parker(prefix, clust_path):
    fam = [l.split() for l in open(prefix + '.fam')]
    ids = [f[1] for f in fam]
    n_ind = len(fam)
    bim = [l.split() for l in open(prefix + '.bim')]
    n_snp = len(bim)
    raw = np.fromfile(prefix + '.bed', dtype=np.uint8)
    bps = (n_ind + 3) // 4
    data = raw[3:].reshape(n_snp, bps)
    codes = np.empty((n_snp, bps * 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k::4] = (data >> (2 * k)) & 3
    dos = np.array([2, -1, 1, 0], dtype=np.int8)[codes[:, :n_ind]]
    truth = {}
    for line in open(clust_path):
        p = line.split()
        if len(p) >= 3:
            truth[p[1]] = p[2]
    keys = [((b[0] if b[0].startswith('chr') else 'chr' + b[0]), b[3]) for b in bim]
    a1 = [b[4] for b in bim]
    return dos, ids, truth, keys, a1


def read_dog10k(gt_path, samples_path, key_index, a1):
    samples = [s.strip() for s in open(samples_path) if s.strip()]
    n_ind = len(samples)
    dos = np.full((len(key_index), n_ind), -1, dtype=np.int8)
    for line in open(gt_path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 4 + n_ind:
            continue
        idx = key_index.get((f[0], f[1]))
        if idx is None:
            continue
        ref, alt = f[2], f[3]
        if ref != a1[idx] and alt != a1[idx]:
            continue
        ref_is_a1 = (ref == a1[idx])
        row = dos[idx]
        for j, g in enumerate(f[4:4 + n_ind]):
            if len(g) < 3 or g[0] not in '01' or g[2] not in '01':
                continue
            n_alt = (g[0] == '1') + (g[2] == '1')
            row[j] = (2 - n_alt) if ref_is_a1 else n_alt
    truth = {}
    for s in samples:
        m = re.match(r'^([A-Za-z\-]+?)\d+$', s)
        truth[s] = m.group(1) if m else s
    return dos, samples, truth


def freqs(dos, cols, labels_of, breeds):
    """per-breed allele frequency of a1, ignoring missing calls"""
    P = np.zeros((dos.shape[0], len(breeds)), dtype=np.float32)
    for j, b in enumerate(breeds):
        sel = [i for i in cols if labels_of[i] == b]
        sub = dos[:, sel].astype(np.float32)
        sub[sub < 0] = np.nan
        with np.errstate(invalid='ignore'):
            m = np.nanmean(sub, axis=1) / 2.0
        P[:, j] = np.nan_to_num(m, nan=0.5)
    return P


def project(P, X, loo_col=None, loo_n=0, G=None):
    K = P.shape[1]
    bb = (P.T @ X) + 1.0
    Gd = G
    if loo_col is not None and loo_n > 1:
        delta = ((P[:, loo_col] - X / 2.0) / (loo_n - 1)).astype(np.float32)
        u = (P.T @ delta).astype(np.float64)
        Gd = G.copy()
        Gd[loo_col, :] += u
        Gd[:, loo_col] += u
        Gd[loo_col, loo_col] += float(delta @ delta)
        bb[loo_col] += float(delta @ X)
    try:
        R = np.linalg.cholesky(Gd).T
    except np.linalg.LinAlgError:
        return None
    y = solve_triangular(R.T, bb.astype(np.float64), lower=True)
    q, _ = nnls(R, y)
    return q


def main():
    prefix, clust, gt, samp = sys.argv[1:5]
    args = sys.argv[5:]
    which = args[args.index('--panel') + 1] if '--panel' in args else 'merged'
    min_n = int(args[args.index('--min-n') + 1]) if '--min-n' in args else 2
    cross = '--cross' in args

    pdos, pids, ptruth, keys, a1 = read_parker(prefix, clust)
    key_index = {k: i for i, k in enumerate(keys)}
    ddos, dids, dtruth = read_dog10k(gt, samp, key_index, a1)

    # keep SNPs typed in both sources
    ok = (pdos >= 0).any(axis=1) & (ddos >= 0).any(axis=1)
    print(f"SNPs usable in both sources: {ok.sum()} of {len(keys)}")
    pdos, ddos = pdos[ok], ddos[ok]

    dos = np.concatenate([pdos, ddos], axis=1)
    ids = pids + dids
    src = np.array([0] * len(pids) + [1] * len(dids))
    labels_of = {}
    for i, s in enumerate(ids):
        labels_of[i] = (ptruth if i < len(pids) else dtruth).get(s)
    del pdos, ddos

    use = {'parker': [i for i in range(len(ids)) if src[i] == 0 and labels_of[i]],
           'dog10k': [i for i in range(len(ids)) if src[i] == 1 and labels_of[i]],
           'merged': [i for i in range(len(ids)) if labels_of[i]]}[which]
    cnt = Counter(labels_of[i] for i in use)
    breeds = sorted(b for b, n in cnt.items() if n >= min_n)
    cols = [i for i in use if labels_of[i] in breeds]
    print(f"panel={which}  individuals={len(cols)}  breeds(n>={min_n})={len(breeds)}")

    P = freqs(dos, cols, labels_of, breeds)
    bidx = {b: j for j, b in enumerate(breeds)}
    K = len(breeds)
    G = ((P.T @ P) + np.ones((K, K), dtype=np.float32)).astype(np.float64)
    G[np.diag_indices_from(G)] += 1e-6
    exp = (2.0 * P.mean(axis=1)).astype(np.float32)

    # platform agreement for breeds present in BOTH sources
    if which == 'merged':
        shared = [b for b in breeds
                  if any(src[i] == 0 and labels_of[i] == b for i in cols)
                  and any(src[i] == 1 and labels_of[i] == b for i in cols)]
        if shared:
            Pp = freqs(dos, [i for i in cols if src[i] == 0], labels_of, shared)
            Pd = freqs(dos, [i for i in cols if src[i] == 1], labels_of, shared)
            same = np.mean((Pp - Pd) ** 2, axis=0)
            other = []
            for j in range(min(len(shared), 40)):
                k = (j + 7) % len(shared)
                other.append(np.mean((Pp[:, j] - Pd[:, k]) ** 2))
            print(f"\nplatform check on {len(shared)} breeds present in both sources:")
            print(f"   same breed, Parker vs Dog10K frequencies : {same.mean():.4f}")
            print(f"   DIFFERENT breeds across sources          : {np.mean(other):.4f}")
            print("   (same-breed distance must be much smaller, or pooling"
                  " would inject platform structure)")

    # ── evaluate ─────────────────────────────────────────────────
    if cross:
        test = [i for i in range(len(ids))
                if labels_of[i] in bidx and src[i] != (0 if which == 'parker' else 1)]
        mode = f"HELD OUT: panel={which}, testing the other source"
    else:
        test = cols
        mode = "leave-one-out"

    per = defaultdict(lambda: [0, 0])
    conf = defaultdict(Counter)
    for n, i in enumerate(test):
        b = labels_of[i]
        X = dos[:, i].astype(np.float32)
        m = X < 0
        if m.any():
            X[m] = exp[m]
        q = project(P, X, None if cross else bidx[b],
                    cnt[b] if not cross else 0, G)
        if q is None:
            continue
        got = breeds[int(np.argmax(q))]
        per[b][1] += 1
        if got == b:
            per[b][0] += 1
        else:
            conf[b][got] += 1
        if n % 200 == 0:
            print(f"   scored {n}/{len(test)}", end='\r', flush=True)
    print(" " * 40, end='\r')

    ok_n = sum(v[0] for v in per.values())
    tot = sum(v[1] for v in per.values())
    print(f"\n=== {which} panel, {mode} ===")
    print(f"individuals scored {tot}   top-1 correct {ok_n}  ({100*ok_n/tot:.2f}%)")
    for lo, hi, lab in ((2, 3, 'n=2-3'), (4, 7, 'n=4-7'), (8, 10**9, 'n>=8')):
        grp = [b for b in per if lo <= cnt[b] <= hi]
        c = sum(per[b][0] for b in grp); t = sum(per[b][1] for b in grp)
        if t:
            print(f"   {lab:<7} {100*c/t:5.1f}%  over {t:>5} dogs, {len(grp):>3} breeds")
    allerr = [(b, g, k) for b in conf for g, k in conf[b].items()]
    tot_err = sum(k for _, _, k in allerr)
    vill = sum(k for b, g, k in allerr if b.startswith('VILL') and g.startswith('VILL'))
    def stem(x):
        return re.sub(r'_(Italy|China|ArabPen|CentAsia|Tribal|NAfrica|Mali)$', '', x)
    samestem = sum(k for b, g, k in allerr
                   if stem(b) != b or stem(g) != g) and sum(
                   k for b, g, k in allerr if stem(b) == stem(g))
    print(f"\n   error breakdown of {tot_err} misassignments:")
    print(f"     village dog -> other village dog : {vill:>4}  ({100*vill/tot_err:.0f}%)")
    print(f"     regional variant of same code    : {samestem:>4}  ({100*samestem/tot_err:.0f}%)")
    print(f"     everything else                  : {tot_err-vill-samestem:>4}")
    worst = sorted(((v[0]/v[1], v[1], b) for b, v in per.items() if v[1] >= 4))[:15]
    print("\n   weakest breeds (n>=4):")
    for frac, n, b in worst:
        got = ', '.join(f"{g}({k})" for g, k in conf[b].most_common(2))
        print(f"     {b:<12} {int(frac*n):>3}/{n:<3} -> {got}")


if __name__ == '__main__':
    main()
