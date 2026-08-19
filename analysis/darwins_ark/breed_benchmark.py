#!/usr/bin/env python3
"""Benchmark our 230-breed composition panel against Darwin's Ark published
supervised-ADMIXTURE ancestry, on their genotyped dogs.

Pipeline-verbatim math: NNLS with the non-negative-lasso RHS shift
(BREED_LASSO=0.3), dosage = E[a1] oriented to the panel's a1/a2 with chain
strand flips, per-SNP mean imputation of missing genotypes.
"""
import csv, gzip, collections, numpy as np
from scipy.optimize import nnls
from scipy.linalg import solve_triangular

PANEL='../breed_panel'
B='GeneticData/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155'
LASSO=0.3
COMP=str.maketrans('ACGT','TGCA')

# panel sites in Phat row order
sites=[]
with open(f'{PANEL}/sites.tsv') as f:
    next(f)
    for line in f:
        c,p,a1,a2=line.rstrip('\n').split('\t')
        sites.append((c,int(p),a1,a2))
pos_index={(c,p):i for i,(c,p,_,_) in enumerate(sites)}
n_sites=len(sites)
P=np.load(f'{PANEL}/phat.npy')
breed_labels=[l.strip() for l in open(f'{PANEL}/breeds.txt') if l.strip()]
print(f'panel: {n_sites} sites x {len(breed_labels)} breeds')

# liftover map: canFam4 (chr,pos) -> (DA snpid, strand)
lift={}
for l in gzip.open('bim_canFam4_map.tsv.gz','rt'):
    f=l.split()
    lift[(f[3],int(f[4]))]=(f[0],f[5])

# DA bim: snpid -> (row, a1_bim, a2_bim); dosage in bed = copies of a1_bim
fam=[l.split()[0] for l in open(B+'.fam')]
n=len(fam)
da_row={}
for i,l in enumerate(open(B+'.bim')):
    f=l.split()
    da_row[f[1]]=(i,f[4],f[5])

# panel site -> (da bed row, flip01) where dosage_of_panel_a1 = d or 2-d
rows=[]; flips=[]; site_idx=[]
for i,(c,p,a1,a2) in enumerate(sites):
    hit=lift.get((c,p))
    if not hit: continue
    snpid,strand=hit
    da=da_row.get(snpid)
    if not da: continue
    r,b1,b2=da
    if strand=='-':
        b1,b2=b1.translate(COMP),b2.translate(COMP)
    if {a1,a2}!={b1,b2}: continue
    rows.append(r); flips.append(a1!=b1); site_idx.append(i)
rows=np.array(rows); flips=np.array(flips); site_idx=np.array(site_idx)
print(f'{len(rows)} panel sites matched in Darwin’s Ark data '
      f'({100*len(rows)/n_sites:.1f}%), {flips.sum()} orientation flips')

# extract dosages for ALL DA dogs at those rows (chunked)
LUT=np.zeros((256,4),dtype=np.float32)
for byte in range(256):
    for k in range(4):
        cc=(byte>>(2*k))&0b11
        LUT[byte,k]=2 if cc==0 else (1 if cc==2 else (0 if cc==3 else np.nan))
bpf=(n+3)//4
order=np.argsort(rows)
D=np.empty((len(rows),n),dtype=np.float32)
with open(B+'.bed','rb') as fh:
    fh.read(3)
    CH=40000
    for s in range(0,len(order),CH):
        sel=order[s:s+CH]; rr=rows[sel]
        fh.seek(3+int(rr[0])*bpf)
        span=int(rr[-1])-int(rr[0])+1
        buf=np.frombuffer(fh.read(span*bpf),dtype=np.uint8).reshape(span,bpf)
        D[sel]=LUT[buf[rr-rr[0]]].reshape(len(rr),-1)[:,:n]
# orient to panel a1; mean-impute
mu=np.nanmean(D,axis=1)
nanm=np.isnan(D)
D=np.where(nanm,mu[:,None],D)
D[flips]=2.0-D[flips]
print('dosage matrix ready', D.shape)

# shared Gram (all dogs share the matched-site mask)
P_v=P[site_idx,:].astype(np.float64)
A=np.vstack([P_v,np.ones((1,P_v.shape[1]))])
G=A.T@A
lam_shift=0.5*LASSO*float(np.mean(np.diag(P_v.T@P_v)))
G[np.diag_indices_from(G)]+=1e-6
R=np.linalg.cholesky(G).T

# truth ancestry
anc=collections.defaultdict(dict)
for r in csv.DictReader(open('DarwinsArk/DarwinsArk/DarwinsArk_20191115_breedcalls.csv')):
    anc[r['dog']][r['breed']]=float(r['pct'])

def norm(b): return b.strip().lower().replace('_',' ').replace('-',' ').replace("'",'')
ours_by_norm={norm(b):j for j,b in enumerate(breed_labels)}
# their breed name -> our column (token fallback)
import re
def tokens(name): return frozenset(w.rstrip('s') for w in re.sub(r'[^a-z ]',' ',norm(name)).split() if w)
ours_tokens={tokens(b):j for j,b in enumerate(breed_labels)}
def map_breed(b):
    j=ours_by_norm.get(norm(b))
    if j is None: j=ours_tokens.get(tokens(b))
    return j

results=[]
for di,dog in enumerate(fam):
    truth=anc.get(dog)
    if not truth: continue
    x=D[:,di].astype(np.float64)
    b=np.hstack([x,[1.0]])
    rhs=A.T@b-lam_shift
    q,_=nnls(R,solve_triangular(R.T,rhs,lower=True))
    q=q/(q.sum()+1e-12)
    results.append((dog,q,truth))
print(f'{len(results)} dogs scored')
np.save('breed_benchmark_results.npy',
        {'fam':[r[0] for r in results],'q':np.array([r[1] for r in results]),
         'truth':[r[2] for r in results],'breed_labels':breed_labels},allow_pickle=True)

# metrics
top1=0; n_eval=0; tvs=[]; rs=[]
mixed_top1=0; mixed_n=0
for dog,q,truth in results:
    tmax=max(truth.values())
    tt_sorted=sorted(truth.items(),key=lambda kv:-kv[1])
    our_top=breed_labels[int(np.argmax(q))]
    their_top=tt_sorted[0][0]
    jt=map_breed(their_top)
    match=(jt is not None and jt==int(np.argmax(q)))
    n_eval+=1; top1+=match
    if tmax<0.45: mixed_n+=1; mixed_top1+=match
    # proportion vector comparison on their breed space
    ours_on_theirs=[]; theirs=[]
    for tb,tp in truth.items():
        j=map_breed(tb)
        if j is None: continue
        ours_on_theirs.append(q[j]); theirs.append(tp)
    if len(theirs)>=3:
        rs.append(np.corrcoef(ours_on_theirs,theirs)[0,1])
        tvs.append(0.5*sum(abs(a-b) for a,b in zip(ours_on_theirs,theirs)))
print(f'\ntop-1 breed agreement: {top1}/{n_eval} = {100*top1/n_eval:.1f}%')
print(f'  highly admixed (<45% top): {mixed_top1}/{mixed_n} = {100*mixed_top1/max(mixed_n,1):.1f}%')
print(f'per-dog proportion correlation (median): {np.median(rs):.3f}')
print(f'per-dog total-variation distance (median): {np.median(tvs):.3f}')
