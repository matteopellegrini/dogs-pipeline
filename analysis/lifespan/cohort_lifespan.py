"""Cohort check of the lifespan card's estimate (mirror of dogs-app/lib/lifespan.ts).
  python3 cohort_lifespan.py lifespanTable.json RESULTS_DIR > table.tsv"""
import glob, json, math, os, sys
T = json.load(open(sys.argv[1]))
out = ['kit\tmedian\tcoverage\tpred_kg\tsex\ttop']
for d in sorted(glob.glob(f'{sys.argv[2]}/*/')):
    try: comp = json.load(open(d + 'breed_result.json'))['breed_composition']
    except Exception: continue
    try: w = json.load(open(d + 'prs_result.json'))['physical_traits']['weight_kg']
    except Exception: w = {}
    parts = []
    for b in comp:
        parts += [(c['code'], c['proportion']) for c in b['components']] if b.get('components') else [(b['breed'], b['proportion'])]
    kg = w.get('pred_kg')
    size = T['size_curve']['a'] + T['size_curve']['b'] * math.log(kg) if kg else T['overall_median']
    tot = cov = acc = 0
    for code, p in parts:
        if p <= 0: continue
        tot += p; row = T['breeds'].get(code)
        if row: cov += p; acc += p * math.log(row['median'])
        else: acc += p * math.log(size)
    if not tot: continue
    sex = w.get('sex_used')
    m = math.exp(acc / tot) * T['sex'].get(sex, 1)
    top = max(comp, key=lambda b: b['proportion'])
    out.append(f"{os.path.basename(d.rstrip('/'))}\t{m:.2f}\t{cov/tot:.2f}\t{kg}\t{sex}\t{top.get('breed_name')} {top['proportion']:.2f}")
print('\n'.join(out))
