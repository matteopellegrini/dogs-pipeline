#!/usr/bin/env python3
"""assoc + liftover map -> canFam4 scoring weights.

    make_weights.py ASSOC P_CUT OUT.tsv

OUT columns: chr4 pos4 effect_allele oth_allele beta af
Alleles are complemented where the chain says the canFam4 locus is on the
minus strand — the chain's strand call is authoritative, so palindromic
SNPs stay usable."""
import gzip, sys
COMP=str.maketrans('ACGT','TGCA')
assoc, pcut, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
m={}
for l in gzip.open('bim_canFam4_map.tsv.gz','rt'):
    f=l.split()
    m[f[0]]=(f[3],int(f[4]),f[5])
n=w=0
with open(out,'w') as o:
    with open(assoc) as fh:
        next(fh)
        for l in fh:
            f=l.split('\t')
            n+=1
            if float(f[11])>pcut: continue
            hit=m.get(f[1])
            if not hit: continue
            c4,p4,strand=hit
            ea,oa=f[4],f[5]
            if strand=='-':
                ea,oa=ea.translate(COMP),oa.translate(COMP)
            o.write(f'{c4}\t{p4}\t{ea}\t{oa}\t{f[7]}\t{f[6]}\n')
            w+=1
print(f'{w} weights (of {n} tested) at p<={pcut}', file=sys.stderr)
