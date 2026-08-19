#!/usr/bin/env python3
"""GEMMA inputs for individual-level Darwin's Ark PRS.

pheno.tsv     one column per (trait x fold): 'full' + 5 dog-level CV folds
              (real individual phenotypes, so dog-level folds are legitimate —
              no breed pseudo-replication as in the AKC harness)
covar_<t>.txt per-trait covariates, paper-matched: intercept, sex, data type
              (low-pass vs Axiom), age at survey (trait-specific), and size
              for the behaviour factors (not for size itself)
columns.tsv   trait, fold, 1-based column for gemma -n
"""
import csv, hashlib

FAM='GeneticData/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155.fam'
D='DarwinsArk/DarwinsArk'
N_FOLDS=5

fam=[l.split()[0] for l in open(FAM)]
famset=set(fam)

meta={}
for r in csv.DictReader(open(f'{D}/DarwinsArk_20191115_dogs.csv')):
    if r['id'] in famset:
        meta[r['id']]=r

scores={}   # dog -> factor -> (score.norm, age)
for r in csv.DictReader(open(f'{D}/DarwinsArk_20191115_factor_scores.csv')):
    if r['factor'] and r['dog'] in famset and r['score.norm'] not in ('','NA'):
        f=int(r['factor'])
        if f<=8:
            age=None
            try: age=float(r['age'])
            except (ValueError,TypeError): pass
            scores.setdefault(r['dog'],{})[f]=(float(r['score.norm']), age)

FACTOR_NAMES={1:'human_sociability',2:'arousal_level',3:'toy_motor_patterns',
              4:'biddability',5:'agonistic_threshold',6:'dog_sociability',
              7:'environmental_engagement',8:'proximity_seeking'}

def sizeval(d):
    v=meta[d].get('size','NA')
    try: return float(v)
    except ValueError: return None

fold_of={d:int(hashlib.sha1(d.encode()).hexdigest(),16)%N_FOLDS+1 for d in fam}
open('gemma/folds.tsv','w').write(''.join(f'{d}\t{fold_of[d]}\n' for d in fam))

traits=[(FACTOR_NAMES[i], i) for i in range(1,9)]+[('size', None)]
cols, header = [], []
for tname, fnum in traits:
    vals, ages = [], []
    for d in fam:
        if fnum is None:
            v=sizeval(d); a=None
        else:
            v,a=scores.get(d,{}).get(fnum,(None,None))
        vals.append(v); ages.append(a)
    n_ok=sum(v is not None for v in vals)
    for fold in ['full']+[str(k) for k in range(1,N_FOLDS+1)]:
        col=['NA' if (v is None or (fold!='full' and fold_of[d]==int(fold))) else f'{v}'
             for d,v in zip(fam,vals)]
        cols.append(col); header.append((tname,fold))

    # per-trait covariates (GEMMA drops NA-pheno dogs itself; covars must be
    # complete for everyone, so mean-impute)
    def mi(xs):
        have=[x for x in xs if x is not None]
        mu=sum(have)/len(have) if have else 0.0
        return [x if x is not None else mu for x in xs]
    sexv=[1.0 if meta[d]['sex']=='M' else 0.0 for d in fam]
    dtv=[1.0 if meta[d]['geno_lowpass']=='TRUE' else 0.0 for d in fam]
    agev=mi(ages)
    szv=mi([sizeval(d) for d in fam])
    with open(f'gemma/covar_{tname}.txt','w') as fh:
        for i in range(len(fam)):
            row=[1.0,sexv[i],dtv[i],agev[i]]
            if fnum is not None: row.append(szv[i])
            fh.write('\t'.join(f'{x:.4f}' for x in row)+'\n')
    print(f'{tname:<26} n={n_ok}')

with open('gemma/pheno.tsv','w') as fh:
    for row in zip(*cols): fh.write('\t'.join(row)+'\n')
with open('gemma/columns.tsv','w') as fh:
    for i,(t,f) in enumerate(header,1): fh.write(f'{t}\t{f}\t{i}\n')
print(f'{len(header)} phenotype columns written')
