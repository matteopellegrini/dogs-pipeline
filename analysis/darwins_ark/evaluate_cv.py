#!/usr/bin/env python3
"""Dog-level CV: held-out r for PRS at several p-cut densities.
For each (trait, fold): betas from the fold's GEMMA run (held-out dogs' NA'd),
score ALL cuts in one streaming pass over the bed, correlate held-out dogs'
scores with their true phenotype."""
import numpy as np, sys

CUTS=[1.0, 0.1, 0.01, 1e-4, 1e-6, 1e-8]
TRAITS=['size','human_sociability','biddability']
B='GeneticData/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155'

fam=[l.split()[0] for l in open(B+'.fam')]
n=len(fam)
snp_row={}
for i,l in enumerate(open(B+'.bim')):
    snp_row[l.split()[1]]=i
m=len(snp_row)
cols=[l.rstrip('\n').split('\t') for l in open('gemma/columns.tsv')]
pheno=[l.rstrip('\n').split('\t') for l in open('gemma/pheno.tsv')]

LUT=np.zeros((256,4),dtype=np.float32)
for byte in range(256):
    for k in range(4):
        c=(byte>>(2*k))&0b11
        LUT[byte,k]=2 if c==0 else (1 if c==2 else (0 if c==3 else np.nan))
bpf=(n+3)//4

def score_multi(bmat, rows_used):
    """bmat: (ncuts, m) float32; returns (ncuts, n) scores; streams bed."""
    scores=np.zeros((bmat.shape[0],n))
    order=np.sort(rows_used)
    with open(B+'.bed','rb') as fh:
        assert fh.read(3)==b'\x6c\x1b\x01'
        CH=40000
        for s in range(0,len(order),CH):
            rows=order[s:s+CH]
            fh.seek(3+int(rows[0])*bpf)
            span=int(rows[-1])-int(rows[0])+1
            buf=np.frombuffer(fh.read(span*bpf),dtype=np.uint8).reshape(span,bpf)
            dos=LUT[buf[rows-rows[0]]].reshape(len(rows),-1)[:,:n]
            mu=np.nanmean(dos,axis=1)
            nanmask=np.isnan(dos)
            if nanmask.any():
                dos=np.where(nanmask, mu[:,None], dos)
            scores+=bmat[:,rows]@dos
    return scores

print(f"{'trait':<20}{'fold':>5}"+''.join(f"{c:>10g}" for c in CUTS))
agg={ (t,c):([],[]) for t in TRAITS for c in CUTS }
for trait in TRAITS:
    full_col=next(int(c) for t,f,c in cols if t==trait and f=='full')
    yfull=[r[full_col-1] for r in pheno]
    for fold in '12345':
        col=next(int(c) for t,f,c in cols if t==trait and f==fold)
        ytr=[r[col-1] for r in pheno]
        te=[i for i in range(n) if yfull[i]!='NA' and ytr[i]=='NA']
        beta=np.zeros(m,dtype=np.float32); pv=np.ones(m,dtype=np.float32)
        with open(f'output/{trait}_{fold}.assoc.txt') as fh:
            next(fh)
            for line in fh:
                f2=line.split('\t')
                i=snp_row.get(f2[1])
                if i is not None:
                    beta[i]=float(f2[7]); pv[i]=float(f2[11])
        bmat=np.stack([np.where(pv<=c,beta,0.0) for c in CUTS])
        rows_used=np.where(pv<1.5)[0]   # every tested SNP
        sc=score_multi(bmat,rows_used)
        truth=np.array([float(yfull[i]) for i in te])
        line=f"{trait:<20}{fold:>5}"
        for ci,c in enumerate(CUTS):
            pred=sc[ci,te]
            r=np.corrcoef(pred,truth)[0,1]
            agg[(trait,c)][0].extend(pred.tolist()); agg[(trait,c)][1].extend(truth.tolist())
            line+=f"{r:>10.3f}"
        print(line,flush=True)
print('\npooled held-out r:')
print(f"{'trait':<20}     "+''.join(f"{c:>10g}" for c in CUTS))
for trait in TRAITS:
    line=f"{trait:<20}     "
    for c in CUTS:
        p,t=agg[(trait,c)]
        line+=f"{np.corrcoef(p,t)[0,1]:>10.3f}"
    print(line)
