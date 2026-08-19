#!/usr/bin/env python3
"""Platform refinement: keep only weight sites that exist (allele-matched)
in our Dog10K imputed BCFs, then rebuild each trait's cohort reference on
exactly that subset — removing the beta*2af fallback dilution that
compressed external dogs' percentiles toward the middle."""
import gzip, subprocess, sys, numpy as np

TRAITS=['human_sociability','arousal_level','toy_motor_patterns','biddability',
        'agonistic_threshold','dog_sociability','environmental_engagement',
        'proximity_seeking','size','white_fur','curly_tail','ear_shape',
        'fur_length','fur_texture']
COMP=str.maketrans('ACGT','TGCA')
B='GeneticData/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155'

print('loading platform sites...', flush=True)
plat={}
for l in gzip.open('platform_sites.tsv.gz','rt'):
    f=l.split()
    plat[f[0]+':'+f[1]]=(f[2],f[3])
print(f'{len(plat)} platform sites', flush=True)

print('loading liftover map...', flush=True)
lift={}
for l in gzip.open('bim_canFam4_map.tsv.gz','rt'):
    f=l.split()
    lift[f[0]]=(f[3],f[4],f[5])   # snpid -> chr4,pos4,strand
fam=[l.split()[0] for l in open(B+'.fam')]
n=len(fam)
snp_row={l.split()[1]:i for i,l in enumerate(open(B+'.bim'))}

LUT=np.zeros((256,4),dtype=np.float32)
for byte in range(256):
    for k in range(4):
        c=(byte>>(2*k))&0b11
        LUT[byte,k]=2 if c==0 else (1 if c==2 else (0 if c==3 else np.nan))
bpf=(n+3)//4

for trait in TRAITS:
    kept=[]; rows=[]; betas=[]
    with open(f'output/{trait}_full.assoc.txt') as fh:
        next(fh)
        for line in fh:
            f=line.split('\t')
            if float(f[11])>0.1: continue
            hit=lift.get(f[1])
            if not hit: continue
            c4,p4,strand=hit
            ea,oa=f[4],f[5]
            if strand=='-': ea,oa=ea.translate(COMP),oa.translate(COMP)
            pl=plat.get(c4+':'+p4)
            if not pl or {ea,oa}!={pl[0],pl[1]}: continue
            kept.append((c4,p4,ea,oa,f[7],f[6]))
            r=snp_row.get(f[1])
            rows.append(r); betas.append(float(f[7]))
    with gzip.open(f'artifact_wts_{trait}.tsv.gz','wt') as o:
        for k in kept: o.write('\t'.join(k)+'\n')
    # cohort reference on the SAME subset
    rows=np.array(rows); betas=np.array(betas)
    order=np.argsort(rows)
    rows_s=rows[order]; b_s=betas[order]
    sc=np.zeros(n)
    with open(B+'.bed','rb') as fh:
        fh.read(3)
        CH=40000
        for s in range(0,len(rows_s),CH):
            rr=rows_s[s:s+CH]
            fh.seek(3+int(rr[0])*bpf)
            span=int(rr[-1])-int(rr[0])+1
            buf=np.frombuffer(fh.read(span*bpf),dtype=np.uint8).reshape(span,bpf)
            dos=LUT[buf[rr-rr[0]]].reshape(len(rr),-1)[:,:n]
            mu=np.nanmean(dos,axis=1)
            bad=np.isnan(dos)
            if bad.any(): dos=np.where(bad,mu[:,None],dos)
            sc+=b_s[s:s+CH]@dos
    np.save(f'artifact_ref_{trait}.npy', np.sort(sc))
    print(f'{trait:<26} {len(kept):>7} platform sites kept', flush=True)
