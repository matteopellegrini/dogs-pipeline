"""Does the mitochondrial share of reads change with age?

mtDNA copy-number proxy = chrM reads / autosomal reads (samtools idxstats on the
archived whole-genome BAMs; logs/mito_idxstats.tsv from logs/mito_idx.sh), so
the bacterial share of a swab does not enter. Ages are collection-time ages
from the microbiome age-clock training table. Covariates: sex and predicted
adult weight (prs_result.json), top breed, depth, unmapped (bacterial) share.

  envs/genomics/bin/python3 analysis/mito_age.py   # run on Hoffman from $D"""
import gzip, json, os, sys
import numpy as np
from scipy import stats
import statsmodels.formula.api as smf
import pandas as pd

D = '/u/project/pellegrini/gkislik/dogs'
idx = {}
for line in open(f'{D}/logs/mito_idxstats.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) == 5: idx[f[0]] = tuple(int(x) for x in f[1:])
ages = {r['sample']: r['age'] for r in json.load(gzip.open(sys.argv[1] if len(sys.argv) > 1 else f'{D}/logs/microbiome_training_table.json.gz'))
        if r.get('age') is not None}
rows = []
for s, age in ages.items():
    if s not in idx: continue
    m, auto, other, unm = idx[s]
    if auto < 1_000_000: continue
    res = f'{D}/results_prosper/{s}'
    try:
        w = json.load(open(f'{res}/prs_result.json'))['physical_traits']['weight_kg']
        top = json.load(open(f'{res}/breed_result.json'))['breed_composition'][0]
    except Exception:
        continue
    rows.append(dict(sample=s, age=age, mt=m / auto * 1e3, log_mt=np.log10(m / auto), sex=w.get('sex_used'),
                     kg=w.get('pred_kg_breed') or w['pred_kg'], breed=top['breed_name'],
                     pure=top['proportion'] >= 0.85, log_auto=np.log10(auto), bact=unm / (unm + m + auto + other)))
df = pd.DataFrame(rows)
df['batch'] = df['sample'].str[3:9]  # kit-number prefix = production/sequencing batch
df.to_csv(f'{D}/logs/mito_age_table.tsv', sep='\t', index=False)
print('dogs', len(df), 'age range', df.age.min(), df.age.max(), 'median', df.age.median())
print('chrM reads per 1,000 autosomal: median %.2f IQR %.2f-%.2f' % tuple(df.mt.quantile([.5, .25, .75])))
print('Spearman age vs mt: rho %.3f p %.2g' % stats.spearmanr(df.age, df.mt))
for lo, hi in [(0, 1), (1, 3), (3, 6), (6, 9), (9, 12), (12, 25)]:
    x = df[(df.age >= lo) & (df.age < hi)]
    print(f'  age {lo}-{hi}: n={len(x):4d} median {x.mt.median():.2f}')
for f in ['log_mt ~ age', 'log_mt ~ age + sex + np.log(kg) + log_auto + bact',
          'log_mt ~ age + I(age**2) + sex + np.log(kg) + log_auto + bact']:
    fit = smf.ols(f, df).fit()
    b, p = fit.params['age'], fit.pvalues['age']
    print(f'{f}\n   age beta {b:+.4f} log10/yr ({(10**b - 1) * 100:+.1f}%/yr) p {p:.2g}  R2 {fit.rsquared:.3f}')
# within breed: purebreds of breeds with >= 15 dogs, breed fixed effect
pb = df[df.pure]; big = pb.breed.value_counts(); pb = pb[pb.breed.isin(big[big >= 8].index)]
fit = smf.ols('log_mt ~ age + C(breed) + sex + log_auto + bact', pb).fit()
print(f'within-breed ({pb.breed.nunique()} breeds, n={len(pb)}): age beta {fit.params["age"]:+.4f} '
      f'({(10**fit.params["age"] - 1) * 100:+.1f}%/yr) p {fit.pvalues["age"]:.2g}')
print('Spearman bact vs age: rho %.3f p %.2g' % stats.spearmanr(df.age, df.bact))
print('Spearman bact vs mt:  rho %.3f p %.2g' % stats.spearmanr(df.bact, df.mt))
big = df.batch.value_counts(); db = df[df.batch.isin(big[big >= 5].index)]
fit = smf.ols('log_mt ~ age + C(batch) + sex + np.log(kg) + log_auto + bact', db).fit()
lo, hi = fit.conf_int().loc['age']
print(f'with batch ({db.batch.nunique()} batches, n={len(db)}): age {(10**fit.params["age"] - 1) * 100:+.1f}%/yr '
      f'(95% CI {(10**lo - 1) * 100:+.1f} to {(10**hi - 1) * 100:+.1f}) p {fit.pvalues["age"]:.2g}')
for bt, x in db.groupby('batch'):
    if len(x) >= 30: print(f'  batch {bt}: n={len(x)} median age {x.age.median()} rho %.2f p %.2g' % stats.spearmanr(x.age, x.mt))
print('Spearman kg vs mt: rho %.2f p %.2g' % stats.spearmanr(df.kg, df.mt), ' sex medians', df.groupby('sex').mt.median().round(3).to_dict())

