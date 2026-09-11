#!/usr/bin/env python3
# Constant-rate genetic map (1 cM/Mb) in PLINK format for FLARE/Beagle at the
# 131k panel sites. No canine recombination map is on the cluster; for
# multi-megabase ancestry segments a uniform rate is an acceptable pilot
# approximation (the dog genome-wide average is ~1 cM/Mb). One file per
# chromosome is not needed: FLARE accepts a genome-wide PLINK map.
out = open('constant_1cM_per_Mb.map', 'w')
n = 0
with open('../breed_panel/sites.tsv') as f:
    next(f)
    for l in f:
        chrom, pos = l.split('\t')[:2]
        out.write(f"{chrom}\t{chrom}_{pos}\t{int(pos)/1e6:.6f}\t{pos}\n")
        n += 1
out.close()
print('map rows:', n)
