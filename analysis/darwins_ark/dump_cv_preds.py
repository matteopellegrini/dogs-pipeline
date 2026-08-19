#!/usr/bin/env python3
"""Held-out CV predictions per dog at p<=0.1, written per trait for plotting."""
import numpy as np, os, sys
B='GeneticData/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155'
CUT=0.1
fam=[l.split()[0] for l in open(B+'.fam')]
n=len(fam)
snp_row={l.split()[1]:i for i,l in enumerate(open(B+'.bim'))}
m=len(snp_row)
cols=[l.rstrip('\n').split('\t') for l in open('gemma/columns.tsv')]
pheno=[l.rstrip('\n').split('\t') for l in open('gemma/pheno.tsv')]
LUT=np.zeros((256,4),dtype=np.float32)
for byte in range(256):
    for k in range(4):
        c=(byte>>(2*k))&0b11
        LUT[byte,k]=2 if c==0 else (1 if c==2 else (0 if c==3 else np.nan))
bpf=(n+3)//4
def score(bvec):
    rows_used=np.where(bvec!=0)[0]
    sc=np.zeros(n)
    with open(B+'.bed','rb') as fh:
        fh.read(3)
        CH=40000
        for s in range(0,len(rows_used),CH):
            rows=rows_used[s:s+CH]
            fh.seek(3+int(rows[0])*bpf)
            span=int(rows[-1])-int(rows[0])+1
            buf=np.frombuffer(fh.read(span*bpf),dtype=np.uint8).reshape(span,bpf)
            dos=LUT[buf[rows-rows[0]]].reshape(len(rows),-1)[:,:n]
            mu=np.nanmean(dos,axis=1)
            bad=np.isnan(dos)
            if bad.any(): dos=np.where(bad,mu[:,None],dos)
            sc+=bvec[rows]@dos
    return sc
traits=sys.argv[1:]
for trait in traits:
    if os.path.exists(f'cv_preds_{trait}.tsv'): continue
    full_col=next(int(c) for t,f,c in cols if t==trait and f=='full')
    yfull=[r[full_col-1] for r in pheno]
    outrows=[]
    for fold in '12345':
        fn=f'output/{trait}_{fold}.assoc.txt'
        if not os.path.exists(fn): continue
        col=next(int(c) for t,f,c in cols if t==trait and f==fold)
        ytr=[r[col-1] for r in pheno]
        te=[i for i in range(n) if yfull[i]!='NA' and ytr[i]=='NA']
        beta=np.zeros(m,dtype=np.float32); pv=np.ones(m,dtype=np.float32)
        with open(fn) as fh:
            next(fh)
            for line in fh:
                f2=line.split('\t')
                i=snp_row.get(f2[1])
                if i is not None:
                    beta[i]=float(f2[7]); pv[i]=float(f2[11])
        sc=score(np.where(pv<=CUT,beta,0.0))
        # z within fold model's scored distribution for cross-fold comparability
        mu,sd=sc.mean(),sc.std()
        for i in te:
            outrows.append((fam[i],fold,(sc[i]-mu)/sd,float(yfull[i])))
    with open(f'cv_preds_{trait}.tsv','w') as fh:
        fh.write('dog\tfold\tprs_z\ttruth\n')
        for r in outrows: fh.write(f'{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]}\n')
    print(trait, len(outrows), 'held-out predictions', flush=True)
