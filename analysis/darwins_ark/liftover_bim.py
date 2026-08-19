#!/usr/bin/env python3
"""Map every Darwin's Ark bim site canFam3.1 -> canFam4.
Output: snpid chr3 pos3 chr4 pos4 strand  (only same-chromosome, unique maps).
Strand '-' means the canFam4 forward strand is the complement — scorers must
flip alleles there."""
import sys, gzip
from pyliftover import LiftOver
lo = LiftOver('canFam3ToCanFam4.over.chain.gz')
bim = 'GeneticData/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155.bim'
out = gzip.open('bim_canFam4_map.tsv.gz','wt')
n = ok = cross = multi = un = 0
for line in open(bim):
    f = line.split()
    n += 1
    chrom = 'chr'+f[0] if not f[0].startswith('chr') else f[0]
    pos = int(f[3])
    r = lo.convert_coordinate(chrom, pos-1)
    if not r:
        un += 1; continue
    if len(r) > 1:
        multi += 1; continue
    c4, p4, strand, _ = r[0]
    if c4 != chrom:
        cross += 1; continue
    out.write(f'{f[1]}\t{chrom}\t{pos}\t{c4}\t{p4+1}\t{strand}\n')
    ok += 1
out.close()
print(f'{n} sites: {ok} mapped, {un} unmapped, {cross} cross-chrom, {multi} multi')
