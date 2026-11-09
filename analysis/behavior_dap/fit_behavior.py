#!/usr/bin/env python3
"""
Breed-composition behaviour models trained on Dog Aging Project owner answers, and
how well today's report scores (AKC-trained LMM PRS; Darwin's Ark factors blended
with ancestry) predict the same owner answers.

    fit_behavior.py TABLE.npz [FACTORS.npz] OUT_DIR

Composite traits: items z-scored (DAP-wide), signed so higher = more of the named
trait, averaged over answered items (>= half required). Model per composite:
    y = a + Q b + covariates, ridge on b only (b = breed deviations from the mean),
lambda by 5-fold CV. 'deployable' = breed + sex (what a customer report has);
'full' adds neuter and age. Accuracy = held-out r (5-fold, dogs >= 1 y), overall and
within mixed-breed and purebred dogs. Today's scores are fixed predictors (never fit on
DAP), so their r with the owner composite is a fair comparison.
"""
import json, os, sys
import numpy as np

COMPOSITES = {   # name: ([(item, sign)], report AKC trait, AKC sign, DA factor)
    'stranger_friendliness': ([('db_fear_level_unknown_human_away_from_home', -1), ('db_fear_level_unknown_human_touch', -1),
                               ('db_aggression_level_on_leash_unknown_human', -1)], 'Openness To Strangers', 1, 'human_sociability'),
    'dog_friendliness': ([('db_aggression_level_on_leash_unknown_dog', -1), ('db_aggression_level_unknown_aggressive_dog', -1),
                          ('db_fear_level_unknown_dogs', -1), ('db_fear_level_unknown_aggressive_dog', -1)], 'Good With Other Dogs', 1, 'dog_sociability'),
    'trainability': ([('db_training_obeys_sit_command_frequency', 1), ('db_training_obeys_stay_command_frequency', 1),
                      ('db_training_distraction_frequency', -1), ('db_pulls_leash_frequency', -1)], 'Trainability Level', 1, 'biddability'),
    'energy': ([('db_energetic_frequency', 1), ('db_hyperactive_frequency', 1), ('pa_activity_level', 1)], 'Energy Level', 1, 'arousal_level'),
    'playfulness': ([('db_playful_frequency', 1)], 'Playfulness Level', 1, 'toy_motor_patterns'),
    'affection': ([('db_attention_seeking_follows_humans_frequency', 1), ('db_attention_seeking_sits_close_to_humans_frequency', 1)],
                  'Affectionate With Family', 1, 'proximity_seeking'),
    'watchdog': ([('db_aggression_level_delivery_workers_at_home', 1), ('db_aggression_level_unknown_human_near_yard', 1),
                  ('db_barks_frequency', 1)], 'Watchdog/Protective Nature', 1, 'agonistic_threshold'),
    'barking': ([('db_barks_frequency', 1), ('db_left_alone_barking_frequency', 1)], 'Barking Level', 1, None),
    'adaptability': ([('db_fear_level_loud_noises', -1), ('db_fear_level_unknown_objects_outside', -1),
                      ('db_fear_level_unknown_situations', -1)], 'Adaptability Level', 1, 'environmental_engagement'),
    'separation_distress': ([('db_left_alone_restlessness_frequency', 1), ('db_left_alone_barking_frequency', 1),
                             ('db_left_alone_scratching_frequency', 1)], None, 1, None),
    'house_soiling': ([('db_urinates_in_home_frequency', 1), ('db_urinates_alone_frequency', 1), ('db_defecates_alone_frequency', 1)], None, 1, None),
    'prey_chase': ([('db_chases_birds_frequency', 1), ('db_chases_squirrels_frequency', 1)], None, 1, None),
    'resource_guarding': ([('db_aggression_level_toys_taken_away', 1), ('db_aggression_level_food_taken_away', 1),
                           ('db_aggression_level_approached_while_eating', 1)], None, 1, None),
    'handling_sensitivity': ([('db_fear_level_nails_clipped_at_home', 1), ('db_fear_level_bathed_at_home', 1)], None, 1, None),
    'escaping': ([('db_escapes_home_or_property_frequency', 1)], None, 1, None),
    'excitability': ([('db_excitement_level_before_walk', 1), ('db_excitement_level_before_car_ride', 1)], None, 1, None),
    'chewing': ([('db_chews_inappropriate_objects_frequency', 1)], None, 1, None),
    'appetite': ([('df_appetite', 1)], None, 1, None),
}
LAMBDAS = [0.3, 1, 3, 10, 30, 100]
K = 5

tab = sys.argv[1]; fac = sys.argv[2] if len(sys.argv) > 3 else None; out = sys.argv[-1]
os.makedirs(out, exist_ok=True)
t = np.load(tab, allow_pickle=True)
Q = t['Q'].astype(np.float64); Q /= Q.sum(1, keepdims=True)
items = [str(i) for i in t['items']]; Y = t['Y']
male, neut, age, mixed = t['male'], t['neutered'], t['age'], t['mixed'].astype(bool)
akc = {str(n): t['akc_score'][:, j] for j, n in enumerate(t['akc_traits'])}
blend = {}
if fac:
    f = np.load(fac, allow_pickle=True)
    blend = {str(n): f['blended'][j] for j, n in enumerate(f['traits'])}
base = np.isfinite(male) & np.isfinite(neut) & np.isfinite(age) & (age >= 1)

def zcol(c):
    v = Y[:, items.index(c)].astype(float)
    return (v - np.nanmean(v)) / np.nanstd(v)

def composite(spec):
    M = np.column_stack([s * zcol(c) for c, s in spec])
    k = np.isfinite(M).sum(1)
    y = np.where(k >= max(1, (M.shape[1] + 1) // 2), np.nanmean(np.where(np.isfinite(M), M, np.nan), 1), np.nan)
    return (y - np.nanmean(y)) / np.nanstd(y)

def fit(Qt, Ct, yt, lam):
    """y = a + Q b + C g, ridge lam on b only (Q rows sum to 1, so b are breed deviations)."""
    n, p = Qt.shape
    A = np.column_stack([np.ones(n), Qt, Ct])
    R = np.zeros(A.shape[1]); R[1:1 + p] = lam
    return np.linalg.solve(A.T @ A + np.diag(R), A.T @ yt)

def pred(coef, Qx, Cx):
    return np.column_stack([np.ones(len(Qx)), Qx, Cx]) @ coef

def folds(n, seed):
    idx = np.random.RandomState(seed).permutation(n)
    return [idx[k::K] for k in range(K)]

def cv(Qd, Cd, y, seed=0):
    best = None
    for lam in LAMBDAS:                                   # lambda by inner CV on the first split
        p = np.zeros(len(y))
        for te in folds(len(y), seed + 100):
            tr = np.setdiff1d(np.arange(len(y)), te)
            p[te] = pred(fit(Qd[tr], Cd[tr], y[tr], lam), Qd[te], Cd[te])
        r = np.corrcoef(p, y)[0, 1]
        if best is None or r > best[1]: best = (lam, r)
    lam = best[0]
    p = np.zeros(len(y))
    for te in folds(len(y), seed):
        tr = np.setdiff1d(np.arange(len(y)), te)
        p[te] = pred(fit(Qd[tr], Cd[tr], y[tr], lam), Qd[te], Cd[te])
    return p, lam

def r(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[m], b[m])[0, 1]) if m.sum() > 30 else float('nan')

rows, models = [], {}
hdr = (f"{'trait':<22} {'n':>5} | DAP breed+sex: {'all':>5} {'mixed':>5} {'pure':>5} | full {'all':>5} | "
       f"today AKC {'all':>5} {'mixed':>5} | today DA {'all':>5} {'mixed':>5} | sex x neut age")
print(hdr)
for name, (spec, akc_t, akc_sign, da_t) in COMPOSITES.items():
    y_all = composite(spec)
    m = base & np.isfinite(y_all)
    y = y_all[m]; Qd = Q[m]; mx = mixed[m]
    Cs = np.column_stack([male[m]])
    Cf = np.column_stack([male[m], neut[m], age[m], age[m] ** 2, (age[m] < 2).astype(float)])
    p_dep, lam = cv(Qd, Cs, y)
    p_full, _ = cv(Qd, Cf, y)
    co = fit(Qd, Cf, y, lam)
    row = {'trait': name, 'n': int(m.sum()), 'lambda': lam,
           'r_dep': r(p_dep, y), 'r_dep_mixed': r(p_dep[mx], y[mx]), 'r_dep_pure': r(p_dep[~mx], y[~mx]),
           'r_full': r(p_full, y), 'coef_male': float(co[-5]), 'coef_neut': float(co[-4]),
           'akc_trait': akc_t, 'da_factor': da_t}
    if akc_t:
        a = akc_sign * akc[akc_t][m]
        row.update(r_akc=r(a, y), r_akc_mixed=r(a[mx], y[mx]))
    if da_t and da_t in blend:
        b = blend[da_t][m]
        row.update(r_da=r(b, y), r_da_mixed=r(b[mx], y[mx]))
    rows.append(row)
    g = lambda k: f"{row[k]:+.2f}" if k in row and np.isfinite(row[k]) else '   - '
    print(f"{name:<22} {row['n']:>5} | {'':14} {g('r_dep'):>5} {g('r_dep_mixed'):>5} {g('r_dep_pure'):>5} | "
          f"{'':4} {g('r_full'):>5} | {'':9} {g('r_akc'):>5} {g('r_akc_mixed'):>5} | {'':8} {g('r_da'):>5} {g('r_da_mixed'):>5} "
          f"| {row['coef_male']:+.2f} {row['coef_neut']:+.2f}")
    # deployable model on all eligible dogs (breed + sex)
    cd = fit(Qd, Cs, y, lam)
    models[name] = {'lambda': lam, 'intercept': float(cd[0]), 'breed_dev': cd[1:-1].tolist(), 'male': float(cd[-1]),
                    'n': int(m.sum()), 'items': spec, 'y_sd_note': 'composite is z-scored over DAP dogs'}

json.dump(rows, open(f'{out}/behavior_cv.json', 'w'), indent=1)
json.dump({'breeds': [str(b) for b in t['breeds']], 'models': models}, open(f'{out}/behavior_models.json', 'w'))
print(f'wrote {out}/behavior_cv.json, behavior_models.json')
