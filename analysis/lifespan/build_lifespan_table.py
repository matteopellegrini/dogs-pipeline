"""Breed lifespan table for the Genetic Health lifespan card (2026-11-10).

Source: McMillan et al. 2024 Sci Rep 14:531 (doi 10.1038/s41598-023-50458-w),
584,734 UK dogs (284,734 deceased), Kaplan-Meier median survival for 155 breeds
(Supplementary Table 3) and survival by decimal year (Supplementary Table 1).
Supplements downloaded to analysis/lifespan/src (untracked) and dumped with
docx_tables.py.

Writes dogs-app/lib/lifespanTable.json:
  breeds      our panel code -> median / 95% CI / n from the paper
  size_curve  median lifespan ~ a + b*log(kg) fitted over the matched breeds
              (kg = DAP-trained breed weights, reference_json/weight_breed_model.json);
              used for ancestry with no published figure (village dogs, rare breeds)
  quantiles   the crossbred survival curve's quartiles as multiples of its median;
              the card scales them to the dog's median (accelerated failure time,
              the paper's own 'x-fold faster time to death' framing)
  sex         female / male multipliers from the paper's sex medians (12.7 / 12.4)
"""
import json, math, re
import numpy as np

D = '/Users/matteopellegrini/Downloads/dogs'
S = f'{D}/analysis/lifespan/src'
norm = lambda s: re.sub(r'[^a-z0-9]+', '', s.replace('_', ' ').lower())

W = json.load(open(f'{D}/reference_json/weight_breed_model.json'))
codes = {norm(c): c for c in W['breeds']}
# Paper name -> our panel code where the names differ.
ALIAS = {
    'Akita': 'AMERICAN_AKITA',              # the paper collapses American + Japanese Akita
    'Bulldog': 'ENGLISH_BULLDOG', 'American Cocker Spaniel': 'COCKER_SPANIEL',
    'St Bernard': 'SAINT_BERNARD', 'Belgian Malenois': 'BELGIAN_MALINOIS',
    'Caucasian Shepherd Dog': 'CAUCASIAN_OVCHARKA', 'Chinese Shar Pei': 'SHAR_PEI',
    'Kelpie': 'AUSTRALIAN_KELPIE', 'Petite Basset Griffon Vendeen': 'PETIT_BASSET_GRIFFON_VENDEEN',
    'Presa Canario': 'DOGO_CANARIO', 'Swedish Valhund': 'SWEDISH_VALLHUND',
    'Toy Manchester Terrier': 'MANCHESTER_TERRIER',
}
# One panel population covers several of the paper's breeds: pool them,
# deaths-weighted (the panel cannot tell a miniature from a standard Dachshund;
# the paper reports a single Poodle figure for all three sizes).
POOL = {'Dachshund': 'DACHSHUND', 'Miniature Dachshund': 'DACHSHUND'}
SPLIT = {'Poodle': ['STANDARD_POODLE', 'MINIATURE_POODLE', 'TOY_POODLE']}
MIN_N = 100   # rare breeds (e.g. Caucasian Shepherd, 52 dogs, 5.4 y) are too noisy to apply
assert all(v in W['breed_weight'] for v in list(ALIAS.values()) + list(POOL.values()) + sum(SPLIT.values(), []))
rows, pooled = {}, {}
for line in open(f'{S}/moesm8.tsv'):
    c = line.rstrip('\n').split('\t')
    if len(c) < 6 or c[0] in ('', 'Breed', 'Crossbreed') or c[0].startswith('##'): continue
    name, n, dead, med = c[0], int(c[1]), int(c[2]), float(c[3])
    f = lambda x: float(x) if re.match(r'^[0-9.]+$', x) else None
    entry = dict(paper=name, median=med, lo=f(c[4]), hi=f(c[5]), n=n, deaths=dead)
    if name in POOL: pooled.setdefault(POOL[name], []).append(entry); continue
    targets = SPLIT.get(name) or [ALIAS.get(name) or codes.get(norm(name)) or codes.get(norm(re.sub(r' Dog$', '', name)))]
    if targets[0] is None: print('not in panel:', name, n, med); continue
    if n < MIN_N: print('too few dogs:', name, n, med); continue
    for t in targets: rows[t] = entry
for code, es in pooled.items():
    dd = sum(e['deaths'] for e in es)
    rows[code] = dict(paper=' + '.join(e['paper'] for e in es), median=round(sum(e['median'] * e['deaths'] for e in es) / dd, 1),
                      lo=None, hi=None, n=sum(e['n'] for e in es), deaths=dd)
print('panel populations with a figure:', len(rows))

# Size curve over the matched breeds (one point per paper breed, sqrt-deaths weighted).
pts = {}
for code, e in rows.items():
    pts.setdefault(e['paper'], (math.log(W['breed_weight'][code]), e['median'], math.sqrt(e['deaths'])))
x, y, w = map(np.array, zip(*pts.values()))
A = np.vstack([np.ones_like(x), x]).T * w[:, None]
a, b = np.linalg.lstsq(A, y * w, rcond=None)[0]
res = y - (a + b * x)
r = np.corrcoef(x, y)[0, 1]
print(f'size curve: median = {a:.2f} + {b:.2f}*ln(kg)  r={r:.2f}  resid sd {res.std():.2f} y  '
      f'(5 kg {a+b*math.log(5):.1f}, 20 kg {a+b*math.log(20):.1f}, 45 kg {a+b*math.log(45):.1f})')

# Survival-curve shape (Supplementary Table 1, purebreds: 473,681 dogs).
srows = [l.rstrip('\n').split('\t') for l in open(f'{S}/moesm6.tsv')]
srows = [q for q in srows if q and re.match(r'^[0-9.]+$', q[0])]
yrs, surv = [float(q[0]) for q in srows], [float(q[3]) for q in srows]
age_at = lambda p: next(t for t, v in zip(yrs, surv) if v <= 1 - p)
med0 = age_at(0.5)
quant = {k: round(age_at(p) / med0, 3) for k, p in (('p25', .25), ('p75', .75), ('p90', .9))}
print('quantile ratios', quant)

out = dict(
    source='McMillan KM et al. 2024, Longevity of companion dog breeds: those at risk from early death. Sci Rep 14:531',
    doi='10.1038/s41598-023-50458-w',
    population='584,734 UK dogs (284,734 deceased) from breed registries, vet groups, insurers and rehoming charities',
    overall_median=12.5, crossbred_median=12.0, purebred_median=12.7,
    size_curve=dict(a=round(float(a), 3), b=round(float(b), 3), r=round(float(r), 2), resid_sd=round(float(res.std()), 2),
                    note='median lifespan of the paper\'s breeds against ln(DAP adult weight); used for ancestry without a published figure'),
    quantiles=quant,
    sex=dict(female=round(12.7 / 12.55, 4), male=round(12.4 / 12.55, 4)),
    breeds={k: rows[k] for k in sorted(rows)},
)
json.dump(out, open(f'{D}/dogs-app/lib/lifespanTable.json', 'w'), indent=1)
print('wrote dogs-app/lib/lifespanTable.json')
