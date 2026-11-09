#!/usr/bin/env python3
"""
One chromosome of the Dog Aging Project VCF -> breed-panel dosages + partial
Darwin's Ark size PRS, computed the way stage 9 / stage 11 compute them for a
customer (expected dosage from the imputation, allele-checked, missing DA sites
filled with 2*af).

    bcftools query -r CHR -T SITES -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' dap.vcf.gz \
      | extract_chr.py CHR panel_sites.tsv darwins_ark_size_prs.tsv.gz OUT.npz

OUT.npz: panel_idx (row in panel sites.tsv), panel_dos (uint8, 100 x dosage of
a1, 255 = missing; rows x dogs), da_prs (per dog, already includes the 2*af
baseline of every DA site on this chromosome), da_matched, da_total.
"""
import gzip, sys
import numpy as np

chrom, panel_path, da_path, out = sys.argv[1:5]

panel = {}
with open(panel_path) as f:
    next(f)
    for i, line in enumerate(f):
        c, pos, a1, a2 = line.rstrip('\n').split('\t')
        if c == chrom:
            panel[int(pos)] = (i, a1, a2)

da, da_base, da_total = {}, 0.0, 0
with gzip.open(da_path, 'rt') as f:
    for line in f:
        c, pos, ea, oa, beta, af = line.split()
        if c != chrom:
            continue
        beta, af = float(beta), float(af)
        da.setdefault(int(pos), []).append((ea, oa, beta, af))
        da_base += beta * 2.0 * af
        da_total += 1

idx_out, dos_out = [], []
da_prs, da_matched, n = None, 0, None
for line in sys.stdin:
    f = line.rstrip('\n').split('\t')
    pos, ref, alt = int(f[1]), f[2], f[3]
    p, d = panel.get(pos), da.get(pos)
    if p is None and d is None:
        continue
    try:
        ds = np.array(f[4:], dtype=np.float64)
    except ValueError:
        ds = np.array([float(v) if v not in ('.', '') else np.nan for v in f[4:]])
    if n is None:
        n = len(ds)
        da_prs = np.full(n, da_base)
    if p is not None:
        i, a1, a2 = p
        if {ref, alt} == {a1, a2}:
            x = 2.0 - ds if ref == a1 else ds          # dosage of a1, as stage 9
            q = np.where(np.isnan(x), 255, np.clip(np.rint(x * 100), 0, 200)).astype(np.uint8)
            idx_out.append(i)
            dos_out.append(q)
    if d is not None:
        for ea, oa, beta, af in d:
            if {ea, oa} != {ref, alt}:
                continue
            de = ds if ea == alt else 2.0 - ds          # dosage of the effect allele
            de = np.where(np.isnan(de), 2.0 * af, de)
            da_prs += beta * (de - 2.0 * af)
            da_matched += 1

if n is None:
    sys.exit(f'{chrom}: no records')
np.savez_compressed(out, panel_idx=np.array(idx_out, dtype=np.int32),
                    panel_dos=np.array(dos_out, dtype=np.uint8).reshape(len(dos_out), n),
                    da_prs=da_prs, da_matched=da_matched, da_total=da_total)
print(f'{chrom}: {len(idx_out)}/{len(panel)} panel sites, DA {da_matched}/{da_total}, {n} dogs')
