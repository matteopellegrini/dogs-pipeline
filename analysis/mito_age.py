"""Does mitochondrial DNA in the swab change with age?

mtDNA copies per cell = 2 x (chrM reads / 16,728 bp) / (autosomal reads / 2,228,550,668 bp):
mitochondrial depth over nuclear (autosomal) depth, from samtools idxstats on
the archived whole-genome BAMs (logs/mito_idxstats.tsv, logs/mito_idx.sh).
Normalising by nuclear reads keeps the bacterial share of the swab out of it.
Ages: customer-entered age at order time (= age at swab) from
prosper_sample_metadata.tsv; optionally Healthbook birth dates
(logs/healthbook_dobs.json from /api/admin/dobs) for kits without one.
Covariates: DNA sex, predicted adult weight (prs_result.json), top breed,
depth, unmapped (mostly bacterial) share, kit batch (kit-number prefix).

  envs/genomics/bin/python3 analysis/mito_age.py   # on Hoffman from $D"""
import json, os
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

D = '/u/project/pellegrini/gkislik/dogs'
CHRM, AUTO = 16_728, 2_228_550_668
idx = {}
for line in open(f'{D}/logs/mito_idxstats.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) == 5: idx[f[0]] = tuple(int(x) for x in f[1:])
meta = pd.read_csv(f'{D}/prosper_sample_metadata.tsv', sep='\t', dtype={'kit': str})
ages = {r.sample_id: (r.age_years, 'order form') for r in meta.itertuples() if pd.notna(r.age_years)}
hb = f'{D}/logs/healthbook_dobs.json'
if os.path.exists(hb):  # Healthbook birth dates, age at swab already computed
    for k, a in json.load(open(hb)).items():
        ages.setdefault(f'pk-{k}', (a, 'healthbook'))
rows = []
for s, (age, src) in ages.items():
    if s not in idx: continue
    m, auto, other, unm = idx[s]
    if auto < 1_000_000: continue
    res = f'{D}/results_prosper/{s}'
    try:
        w = json.load(open(f'{res}/prs_result.json'))['physical_traits']['weight_kg']
        top = json.load(open(f'{res}/breed_result.json'))['breed_composition'][0]
    except Exception:
        continue
    cn = 2 * (m / CHRM) / (auto / AUTO)
    rows.append(dict(sample=s, age=age, age_src=src, cn=cn, log_cn=np.log10(cn), sex=w.get('sex_used'),
                     kg=w.get('pred_kg_breed') or w['pred_kg'], breed=top['breed_name'], pure=top['proportion'] >= 0.85,
                     log_auto=np.log10(auto), bact=unm / (unm + m + auto + other), batch=s[3:9]))
df = pd.DataFrame(rows)
df.to_csv(f'{D}/logs/mito_age_table.tsv', sep='\t', index=False)
print('dogs', len(df), df.age_src.value_counts().to_dict(), 'age range %.1f-%.0f median %.1f' % (df.age.min(), df.age.max(), df.age.median()))
print('mtDNA copies per cell: median %.0f IQR %.0f-%.0f' % tuple(df.cn.quantile([.5, .25, .75])))
print('Spearman age vs copies: rho %.3f p %.2g' % stats.spearmanr(df.age, df.cn))
for lo, hi in [(0, 1), (1, 3), (3, 6), (6, 9), (9, 12), (12, 30)]:
    x = df[(df.age >= lo) & (df.age < hi)]
    if len(x): print(f'  age {lo}-{hi}: n={len(x):4d} median {x.cn.median():.0f} copies')
pct = lambda b: (10 ** b - 1) * 100
def report(name, f, d):
    fit = smf.ols(f, d).fit(); lo, hi = fit.conf_int().loc['age']
    print(f'{name} (n={len(d)}): age {pct(fit.params["age"]):+.1f}%/yr (95% CI {pct(lo):+.1f} to {pct(hi):+.1f}) p {fit.pvalues["age"]:.2g}  R2 {fit.rsquared:.2f}')
big = df.batch.value_counts(); db = df[df.batch.isin(big[big >= 5].index)]
report('age only', 'log_cn ~ age', df)
report('+ sex, size, depth, bacteria', 'log_cn ~ age + sex + np.log(kg) + log_auto + bact', df)
report('+ kit batch', 'log_cn ~ age + C(batch) + sex + np.log(kg) + log_auto + bact', db)
report('without bacteria term', 'log_cn ~ age + C(batch) + sex + np.log(kg) + log_auto', db)
for bt, x in db.groupby('batch'):
    if len(x) >= 30: print(f'  batch {bt}: n={len(x)} rho %.2f p %.2g' % stats.spearmanr(x.age, x.cn))
pb = df[df.pure]; vc = pb.breed.value_counts(); pb = pb[pb.breed.isin(vc[vc >= 8].index)]
if pb.breed.nunique() >= 2:
    print('purebred breeds with >= 8 aged dogs:', dict(vc[vc >= 8]))
    report('within breed', 'log_cn ~ age + C(breed) + sex + log_auto + bact', pb)
print('bacteria vs age rho %.2f; bacteria vs copies rho %.2f; size vs copies rho %.2f' % (
    stats.spearmanr(df.age, df.bact)[0], stats.spearmanr(df.bact, df.cn)[0], stats.spearmanr(df.kg, df.cn)[0]))
print('sex medians', df.groupby('sex').cn.median().round(0).to_dict())
