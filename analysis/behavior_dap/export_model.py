#!/usr/bin/env python3
"""
Write reference_json/behavior_dap_model.json: the shipped Dog Aging Project behaviour
models (breed composition + sex), refit on every eligible DAP adult with the lambda
fit_behavior.py chose by CV, plus each trait's held-out accuracy and the reference
distribution (percentiles of DAP dogs' own predictions) a report percentile is read from.

    export_model.py TABLE.npz CV.json OUT_JSON
A dog's score: s = intercept + sum_b p_b * breed_dev[b] + male_coef * male (male 1/0, 0.5 if
unknown), over breed_composition_raw normalised to the labels in the model; percentile =
position of s in ref_quantiles (101 DAP quantiles, 0..100).
"""
import datetime, json, sys
import numpy as np
sys.argv, _argv = sys.argv[:1] + ['x', 'x'], sys.argv            # import composites only
import importlib.util
spec = importlib.util.spec_from_file_location('fb', __file__.replace('export_model.py', 'fit_behavior.py'))
src = open(spec.origin).read().split('\nLAMBDAS')[0]
ns = {}; exec(src, ns)
COMPOSITES = ns['COMPOSITES']
sys.argv = _argv
tab, cv_json, out = sys.argv[1:4]

SHIP = {   # trait: (display, blurb, reports instead of)
    'stranger_friendliness': ('Friendliness with strangers', 'Relaxed and friendly with unfamiliar people, rather than fearful or reactive.', ['Openness To Strangers', 'human_sociability']),
    'dog_friendliness': ('Friendliness with other dogs', 'Calm around unfamiliar dogs, rather than fearful or reactive.', ['Good With Other Dogs', 'dog_sociability']),
    'trainability': ('Trainability', 'Responds to sit and stay, stays focused, walks without pulling.', ['Trainability Level', 'biddability']),
    'energy': ('Energy', 'How energetic and active, day to day.', ['Energy Level', 'arousal_level']),
    'affection': ('Closeness to people', 'Follows their people around and likes to sit close.', ['Affectionate With Family', 'proximity_seeking']),
    'watchdog': ('Watchdog alertness', 'Reacts to visitors at the door and people near the yard.', ['Watchdog/Protective Nature', 'agonistic_threshold']),
    'barking': ('Barking', 'How often they bark, including when left alone.', ['Barking Level']),
    'adaptability': ('Calm with noise and new things', 'Unbothered by loud noises, new objects and new situations.', ['Adaptability Level', 'environmental_engagement']),
    'prey_chase': ('Chasing instinct', 'Urge to chase squirrels, birds and other small animals.', []),
    'appetite': ('Appetite', 'How food-motivated they are.', []),
    'escaping': ('Escape attempts', 'Tendency to slip out of the house or yard.', []),
    'handling_sensitivity': ('Sensitivity to grooming', 'Nervousness about nail trims and baths.', []),
    'separation_distress': ('Distress when left alone', 'Restlessness, barking or scratching when alone.', []),
    'resource_guarding': ('Guarding food and toys', 'Growling or snapping when food or toys are taken away.', []),
    'house_soiling': ('House-training accidents', 'Urinating or defecating indoors, especially when alone.', []),
}
NOT_SHIPPED = {'playfulness': 'held-out r 0.09: ancestry barely predicts it', 'excitability': 'r 0.15, not chosen',
               'chewing': 'r 0.07: ancestry barely predicts it'}

t = np.load(tab, allow_pickle=True)
Q = t['Q'].astype(np.float64); Q /= Q.sum(1, keepdims=True)
items = [str(i) for i in t['items']]; Y = t['Y']
male, neut, age = t['male'], t['neutered'], t['age']
base = np.isfinite(male) & np.isfinite(neut) & np.isfinite(age) & (age >= 1)
cv = {r['trait']: r for r in json.load(open(cv_json))}

def zcol(c):
    v = Y[:, items.index(c)].astype(float)
    return (v - np.nanmean(v)) / np.nanstd(v)

model = {'model': 'DAP_behavior_breed_sex', 'built': datetime.date.today().isoformat(),
         'breeds': [str(b) for b in t['breeds']],
         'source': 'Dog Aging Project HLES owner survey (42 behaviour items, 0-4), genotyped dogs >= 1 y',
         'method': ('Per trait, owner answers are z-scored, signed and averaged into one score; a ridge regression '
                    'on breed composition (stage-9 estimator, v2.3 panel) plus sex is fit on all Dog Aging Project '
                    'adults. Accuracy is 5-fold held-out correlation with the owner score.'),
         'not_shipped': NOT_SHIPPED, 'traits': {}}
for name, (display, blurb, replaces) in SHIP.items():
    spec_items = COMPOSITES[name][0]
    M = np.column_stack([s * zcol(c) for c, s in spec_items])
    k = np.isfinite(M).sum(1)
    y = np.where(k >= max(1, (M.shape[1] + 1) // 2), np.nanmean(np.where(np.isfinite(M), M, np.nan), 1), np.nan)
    y = (y - np.nanmean(y)) / np.nanstd(y)
    m = base & np.isfinite(y)
    lam = cv[name]['lambda']
    A = np.column_stack([np.ones(m.sum()), Q[m], male[m]])
    R = np.zeros(A.shape[1]); R[1:1 + Q.shape[1]] = lam
    co = np.linalg.solve(A.T @ A + np.diag(R), A.T @ y[m])
    pred = A @ co
    model['traits'][name] = {
        'display': display, 'blurb': blurb, 'replaces': replaces,
        'items': [[c, s] for c, s in spec_items], 'lambda': lam, 'n_dogs': int(m.sum()),
        'intercept': round(float(co[0]), 6), 'male': round(float(co[-1]), 6),
        'breed_dev': [round(float(x), 5) for x in co[1:-1]],
        'ref_quantiles': [round(float(x), 5) for x in np.percentile(pred, np.arange(101))],
        'cv_r': round(cv[name]['r_dep'], 3), 'cv_r_mixed': round(cv[name]['r_dep_mixed'], 3),
        'cv_r_purebred': round(cv[name]['r_dep_pure'], 3),
        'today_akc_r': None if cv[name].get('r_akc') is None else round(cv[name]['r_akc'], 3),
    }
    print(f'{name:<22} n={m.sum()} lam={lam} cv_r={cv[name]["r_dep"]:.2f} pred sd {pred.std():.3f}')
json.dump(model, open(out, 'w'))
print(f'wrote {out}: {len(model["traits"])} traits')
