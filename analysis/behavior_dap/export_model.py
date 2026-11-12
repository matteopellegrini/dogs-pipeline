#!/usr/bin/env python3
"""
Write reference_json/behavior_dap_model.json: the shipped Dog Aging Project behaviour models.

    export_model.py TABLE.npz OUT_JSON

Per trait (fit_behavior.py composites): ridge on breed composition + sex, with the three
poodle sizes pooled into one POODLE column (doodles split arbitrarily across sizes, and the
Toy Poodle column carried toy-breed effects into them). Lambda by 5-fold CV. Calibration:
logistic regression of 'owner score above the DAP mean' on the out-of-fold prediction, so a
report shows the chance a dog with this ancestry and sex is above the average dog (most dogs
land near 50%: ancestry explains a small part of behaviour). A dog's score:
    s = intercept + male*m + sum_g p_g * group_dev[g]  (composition pooled to groups, normalised)
    chance_above = 100 / (1 + exp(-(cal_a + cal_b * s)))
"""
import datetime, json, sys
import numpy as np

src = open(__file__.replace('export_model.py', 'fit_behavior.py')).read().split('\nLAMBDAS')[0]
ns = {}; exec(src.split('\ntab = sys.argv')[0], ns)
COMPOSITES = ns['COMPOSITES']
tab, out = sys.argv[1:3]

POOL = {'STANDARD_POODLE': 'POODLE', 'MINIATURE_POODLE': 'POODLE', 'TOY_POODLE': 'POODLE'}
SHIP = {   # trait: (display, blurb, comparison phrase, reports instead of)
    'stranger_friendliness': ('Friendliness with strangers', 'Relaxed and friendly with unfamiliar people, rather than fearful or reactive.', 'are friendlier with strangers', ['Openness To Strangers', 'human_sociability']),
    'dog_friendliness': ('Friendliness with other dogs', 'Calm around unfamiliar dogs, rather than fearful or reactive.', 'are calmer around unfamiliar dogs', ['Good With Other Dogs', 'dog_sociability']),
    'affection': ('Closeness to people', 'Follows their people around and likes to sit close.', 'are more people-focused', ['Affectionate With Family', 'proximity_seeking']),
    'trainability': ('Trainability', 'Responds to sit and stay, stays focused, walks without pulling.', 'are more trainable', ['Trainability Level', 'biddability']),
    'energy': ('Energy', 'How energetic and active, day to day.', 'are more energetic', ['Energy Level', 'arousal_level']),
    'adaptability': ('Calm with noise and new things', 'Unbothered by loud noises, new objects and new situations.', 'are calmer with noise and new things', ['Adaptability Level', 'environmental_engagement']),
    'watchdog': ('Watchdog alertness', 'Reacts to visitors at the door and people near the yard.', 'are more watchful', ['Watchdog/Protective Nature', 'agonistic_threshold']),
    'barking': ('Barking', 'How often they bark, including when left alone.', 'bark more', ['Barking Level']),
    'appetite': ('Appetite', 'How food-motivated they are.', 'are more food-motivated', []),
    'house_soiling': ('House-training accidents', 'Urinating or defecating indoors, especially when alone.', 'have more house-training accidents', []),
}
NOT_SHIPPED = {'playfulness': 'held-out r 0.09', 'excitability': 'r 0.15', 'chewing': 'r 0.07',
               'prey_chase': 'calibrated chance stays 38-62% for 90% of customers',
               'escaping': 'stays 41-63%', 'handling_sensitivity': 'stays 41-62%',
               'separation_distress': 'stays 41-56%', 'resource_guarding': 'stays 44-60%'}
LAMBDAS, K = [0.3, 1, 3, 10, 30, 100], 5

t = np.load(tab, allow_pickle=True)
br = [str(b) for b in t['breeds']]
groups = sorted(set(POOL.get(b, b) for b in br))
gi = {g: i for i, g in enumerate(groups)}
Q0 = t['Q'].astype(np.float64); Q0 /= Q0.sum(1, keepdims=True)
Q = np.zeros((len(Q0), len(groups)))
for j, b in enumerate(br):
    Q[:, gi[POOL.get(b, b)]] += Q0[:, j]
items = [str(i) for i in t['items']]; Y = t['Y']
male, neut, age, mixed = t['male'], t['neutered'], t['age'], t['mixed'].astype(bool)
base = np.isfinite(male) & np.isfinite(neut) & np.isfinite(age) & (age >= 1)

def zcol(c):
    v = Y[:, items.index(c)].astype(float)
    return (v - np.nanmean(v)) / np.nanstd(v)

def fit(A, y, lam, p):
    R = np.zeros(A.shape[1]); R[1:1 + p] = lam
    return np.linalg.solve(A.T @ A + np.diag(R), A.T @ y)

def folds(n, seed):
    idx = np.random.RandomState(seed).permutation(n)
    return [idx[k::K] for k in range(K)]

def oof(A, y, lam, p, seed):
    pr = np.zeros(len(y))
    for te in folds(len(y), seed):
        tr = np.setdiff1d(np.arange(len(y)), te)
        pr[te] = A[te] @ fit(A[tr], y[tr], lam, p)
    return pr

def logistic(x, z, iters=50):
    a, b = 0.0, 1.0
    X = np.column_stack([np.ones_like(x), x])
    w = np.array([a, b])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-X @ w))
        g = X.T @ (z - p); H = X.T @ (X * (p * (1 - p))[:, None])
        w += np.linalg.solve(H + 1e-9 * np.eye(2), g)
    return w

model = {'model': 'DAP_behavior_breed_sex_v2', 'built': datetime.date.today().isoformat(),
         'groups': groups, 'pool': POOL,
         'source': 'US pet-dog owner survey (behaviour items scored 0-4), genotyped dogs >= 1 y',
         'method': ('Per trait, owner answers are z-scored, signed and averaged; a ridge regression on breed '
                    'composition (poodle sizes pooled) plus sex is fit on all US pet-dog adults, and a '
                    'logistic calibration on held-out predictions turns it into the chance a dog with this '
                    'ancestry and sex is above the average dog.'),
         'not_shipped': NOT_SHIPPED, 'traits': {}}
print(f"{'trait':<22} {'lam':>4} {'r':>5} {'r_mix':>6} | calibration by OOF-pred quintile: predicted% / observed%")
for name, (display, blurb, phrase, replaces) in SHIP.items():
    spec = COMPOSITES[name][0]
    M = np.column_stack([s * zcol(c) for c, s in spec])
    k = np.isfinite(M).sum(1)
    y = np.where(k >= max(1, (M.shape[1] + 1) // 2), np.nanmean(np.where(np.isfinite(M), M, np.nan), 1), np.nan)
    y = (y - np.nanmean(y)) / np.nanstd(y)
    m = base & np.isfinite(y)
    A = np.column_stack([np.ones(m.sum()), Q[m], male[m]]); yy = y[m]; p = Q.shape[1]
    best = max(LAMBDAS, key=lambda l: np.corrcoef(oof(A, yy, l, p, 100), yy)[0, 1])
    pr = oof(A, yy, best, p, 0)
    r_all, r_mix = np.corrcoef(pr, yy)[0, 1], np.corrcoef(pr[mixed[m]], yy[mixed[m]])[0, 1]
    above = (yy > 0).astype(float)
    ca, cb = logistic(pr, above)
    pc = 1 / (1 + np.exp(-(ca + cb * pr)))
    qs = np.quantile(pr, [0, .2, .4, .6, .8, 1])
    cal = []
    for lo, hi in zip(qs[:-1], qs[1:]):
        s_ = (pr >= lo) & (pr <= hi)
        cal.append(f'{100 * pc[s_].mean():.0f}/{100 * above[s_].mean():.0f}')
    co = fit(A, yy, best, p)
    model['traits'][name] = {
        'display': display, 'blurb': blurb, 'phrase': phrase, 'replaces': replaces,
        'items': [[c, s] for c, s in spec], 'lambda': best, 'n_dogs': int(m.sum()),
        'intercept': round(float(co[0]), 6), 'male': round(float(co[-1]), 6),
        'group_dev': [round(float(x), 5) for x in co[1:-1]],
        'cal_a': round(float(ca), 5), 'cal_b': round(float(cb), 5),
        'base_rate_above': round(float(above.mean()), 4),
        'cv_r': round(float(r_all), 3), 'cv_r_mixed': round(float(r_mix), 3),
    }
    print(f'{name:<22} {best:>4} {r_all:>5.2f} {r_mix:>6.2f} | ' + '  '.join(cal))
json.dump(model, open(out, 'w'))
print(f'wrote {out}: {len(model["traits"])} traits, {len(groups)} groups')
