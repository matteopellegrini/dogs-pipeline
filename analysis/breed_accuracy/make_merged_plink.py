#!/usr/bin/env python3
"""
Write a merged Parker + Dog10K PLINK fileset plus a matching .frq.strat, ready
for supervised SCOPE.

    python3 make_merged_plink.py PARKER_PREFIX CLUST.txt DOG10K_GT.txt \\
        DOG10K_SAMPLES.txt OUT_PREFIX [--step 4] [--min-n 2]

The .bed is written directly from the dosage matrix rather than by merging two
PLINK filesets. SCOPE's documentation warns that alleles must be coded
consistently between the frequency file and the target data, and a PLINK merge
of an array dataset with a BCF invites exactly that class of error (REF/ALT vs
A1/A2, strand flips, ID mismatches). Both sources are already oriented to
dosage-of-a1 against the Parker .bim here, so writing the pack directly keeps
one source of truth for allele order.

PLINK1 bed codes, SNP-major: 00 = hom A1, 01 = missing, 10 = het, 11 = hom A2.
"""
import re
import sys
from collections import Counter

import numpy as np


def read_parker(prefix):
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
    return dos, ids, bim


def main():
    prefix, clust, gt, samp, out = sys.argv[1:6]
    args = sys.argv[6:]
    step = int(args[args.index('--step') + 1]) if '--step' in args else 4
    min_n = int(args[args.index('--min-n') + 1]) if '--min-n' in args else 2

    pdos, pids, bim = read_parker(prefix)
    keys = {}
    for i, b in enumerate(bim):
        c = b[0] if b[0].startswith('chr') else 'chr' + b[0]
        keys[(c, b[3])] = i
    a1 = [b[4] for b in bim]

    dsamples = [s.strip() for s in open(samp) if s.strip()]
    ddos = np.full((len(bim), len(dsamples)), -1, dtype=np.int8)
    for line in open(gt):
        f = line.rstrip('\n').split('\t')
        if len(f) < 4 + len(dsamples):
            continue
        idx = keys.get((f[0], f[1]))
        if idx is None:
            continue
        ref, alt = f[2], f[3]
        if ref != a1[idx] and alt != a1[idx]:
            continue
        ref_is_a1 = (ref == a1[idx])
        row = ddos[idx]
        for j, g in enumerate(f[4:4 + len(dsamples)]):
            if len(g) < 3 or g[0] not in '01' or g[2] not in '01':
                continue
            na = (g[0] == '1') + (g[2] == '1')
            row[j] = (2 - na) if ref_is_a1 else na

    ok = np.where(((pdos >= 0).any(axis=1)) & ((ddos >= 0).any(axis=1)))[0][::step]
    print(f"SNPs kept: {len(ok)} (every {step}th of those typed in both sources)")

    dos = np.concatenate([pdos[ok], ddos[ok]], axis=1)
    ids = pids + dsamples
    del pdos, ddos

    truth = {}
    for line in open(clust):
        p = line.split()
        if len(p) >= 3:
            truth[p[1]] = p[2]
    for s in dsamples:
        m = re.match(r'^([A-Za-z\-]+?)\d+$', s)
        truth[s] = m.group(1) if m else s

    cnt = Counter(truth[s] for s in ids if s in truth)
    breeds = sorted(b for b, n in cnt.items() if n >= min_n)
    keep_ind = [j for j, s in enumerate(ids) if truth.get(s) in set(breeds)]
    dos = dos[:, keep_ind]
    ids = [ids[j] for j in keep_ind]
    n_ind, n_snp = len(ids), len(ok)
    print(f"individuals: {n_ind}   breeds (n>={min_n}): {len(breeds)}")

    with open(out + '.fam', 'w') as fh:
        for s in ids:
            fh.write(f"{truth[s]} {s} 0 0 0 -9\n")
    with open(out + '.bim', 'w') as fh:
        for i in ok:
            b = bim[i]
            fh.write(f"{b[0]}\t{b[1]}\t0\t{b[3]}\t{b[4]}\t{b[5]}\n")

    # 2 -> 00 hom A1, 1 -> 10 het, 0 -> 11 hom A2, missing -> 01
    lut = np.array([0b11, 0b10, 0b00, 0b01], dtype=np.uint8)   # index by dosage, -1 -> last
    bps = (n_ind + 3) // 4
    with open(out + '.bed', 'wb') as fh:
        fh.write(bytes([0x6c, 0x1b, 0x01]))
        pad = np.full((n_snp, bps * 4 - n_ind), -1, dtype=np.int8)
        d2 = np.concatenate([dos, pad], axis=1)
        codes = lut[d2]
        packed = np.zeros((n_snp, bps), dtype=np.uint8)
        for k in range(4):
            packed |= (codes[:, k::4] << (2 * k)).astype(np.uint8)
        packed.tofile(fh)
    print(f"wrote {out}.bed/.bim/.fam")

    # frequency file, cluster-major within SNP, matching PLINK --freq --within
    sub = {}
    for b in breeds:
        sel = [j for j, s in enumerate(ids) if truth[s] == b]
        x = dos[:, sel].astype(np.float32)
        x[x < 0] = np.nan
        sub[b] = x
    freq, nobs = {}, {}
    for b in breeds:
        with np.errstate(invalid='ignore'):
            freq[b] = np.nan_to_num(np.nanmean(sub[b], axis=1) / 2.0, nan=0.0)
        nobs[b] = 2 * np.sum(~np.isnan(sub[b]), axis=1)
    with open(out + '.frq.strat', 'w') as fh:
        fh.write(" CHR                  SNP     CLST   A1   A2      MAF    MAC  NCHROBS\n")
        buf = []
        for i, si in enumerate(ok):
            b_ = bim[si]
            for b in breeds:
                f = float(freq[b][i]); nc = int(nobs[b][i])
                buf.append(f"  {b_[0]:>2} {b_[1]:>20} {b:>8} {b_[4]:>4} {b_[5]:>4} "
                           f"{f:>8.4g} {int(round(f*nc)):>6} {nc:>8}\n")
            if len(buf) > 400000:
                fh.writelines(buf); buf = []
        fh.writelines(buf)
    with open(out + '.breeds', 'w') as fh:
        fh.write("\n".join(breeds) + "\n")
    print(f"wrote {out}.frq.strat  ({n_snp} SNPs x {len(breeds)} clusters)")


if __name__ == '__main__':
    main()
