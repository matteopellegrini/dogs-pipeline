#!/usr/bin/env python3
"""
One chromosome of the Dog Aging Project VCF -> partial Darwin's Ark behaviour-factor
PRS for every DAP dog, computed as stage 11 does for a customer (expected dosage of
the effect allele, allele-checked; sites the dog lacks contribute beta * 2af).

    bcftools query -r CHR -T SITES -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' dap.vcf.gz \
      | extract_factors.py CHR DA_DIR OUT.npz trait1 trait2 ...

OUT.npz: traits, prs (traits x dogs, includes the 2af baseline of every site on
this chromosome), matched, total (per trait).
"""
import gzip, sys
import numpy as np

chrom, da_dir, out, *traits = sys.argv[1:]
sites, base, total = {}, np.zeros(len(traits)), np.zeros(len(traits), dtype=int)
for k, t in enumerate(traits):
    with gzip.open(f'{da_dir}/wts_{t}.tsv.gz', 'rt') as f:
        for line in f:
            c, pos, ea, oa, beta, af = line.split()
            if c != chrom:
                continue
            beta, af = float(beta), float(af)
            sites.setdefault(int(pos), []).append((k, ea, oa, beta, af))
            base[k] += beta * 2.0 * af
            total[k] += 1

prs, matched, n = None, np.zeros(len(traits), dtype=int), None
for line in sys.stdin:
    f = line.rstrip('\n').split('\t')
    hits = sites.get(int(f[1]))
    if hits is None:
        continue
    ref, alt = f[2], f[3]
    try:
        ds = np.array(f[4:], dtype=np.float64)
    except ValueError:
        ds = np.array([float(v) if v not in ('.', '') else np.nan for v in f[4:]])
    if n is None:
        n = len(ds)
        prs = np.repeat(base[:, None], n, axis=1)
    for k, ea, oa, beta, af in hits:
        if {ea, oa} != {ref, alt}:
            continue
        de = ds if ea == alt else 2.0 - ds
        prs[k] += beta * (np.where(np.isnan(de), 2.0 * af, de) - 2.0 * af)
        matched[k] += 1
if n is None:
    sys.exit(f'{chrom}: no records')
np.savez_compressed(out, traits=np.array(traits), prs=prs, matched=matched, total=total)
print(f'{chrom}: ' + ', '.join(f'{t} {m}/{tt}' for t, m, tt in zip(traits, matched, total)) + f'; {n} dogs')
