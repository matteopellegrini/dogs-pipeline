#!/usr/bin/env python3
"""
Write reference_json/weight_breed_model.json from the Dog Aging Project fit
(fit_dap.py chose the recipe): per-breed log weights for all panel labels by
box-constrained ridge (box 2.0 around the AKC prior or the marginal prior),
lambda by 5-fold CV, then a log-space Huber stack on the out-of-fold breed term
and sex. No Darwin's Ark PRS (DAP CV: MAE 4.11 -> 4.10 kg with it). Validation
is external: ProsperK9 customer adults, none used in fitting.

    export_model.py BUILD_DIR run_dog_pipeline.sh OUT_JSON
"""
import json, os, sys, datetime
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from fit_weight_model import (extract_dict, join_prior, metrics, fmt, folds, fit_box_ridge,
                              marginal_prior, huber_ols, LOG_HUBER_C, LAM_GRID)

BOX, MIN_AGE = 2.0, 1.0
build, pipeline_sh, out_json = sys.argv[1:4]
z = np.load(f'{build}/dap.npz')
breeds = [str(b) for b in z['breeds']]; bi = {b: i for i, b in enumerate(breeds)}
qi = {str(d): i for i, d in enumerate(z['ids'])}
rows = [l.rstrip('\n').split('\t') for l in open(f'{build}/dap_pheno.tsv')][1:]
D = [(qi[d], float(w), float(m)) for d, w, a, m, n, s, b in rows if d in qi and a and m and float(a) >= MIN_AGE]
X = z['Q'][[r[0] for r in D]].astype(np.float64); X /= X.sum(1, keepdims=True)
y = np.array([r[1] for r in D]); ly = np.log(y); male = np.array([r[2] for r in D])

akc = join_prior(breeds, extract_dict(open(pipeline_sh).read(), 'BREED_WEIGHT_KG'))
prior_akc = np.array([akc.get(b, np.nan) for b in breeds])
lpr = np.log(marginal_prior(X, y, prior_akc))

def oof(lam, seed):
    p = np.zeros(len(ly))
    for te in folds(len(ly), seed):
        tr = np.setdiff1d(np.arange(len(ly)), te)
        p[te] = X[te] @ fit_box_ridge(X[tr], ly[tr], lpr, lam, box=BOX)
    return p
lam = min(LAM_GRID, key=lambda l: np.mean(np.abs(np.exp(oof(l, 1000)) - y)))
w = fit_box_ridge(X, ly, lpr, lam, box=BOX)
inner = oof(lam, 2000)
b = huber_ols(np.column_stack([inner, male]), ly, c=LOG_HUBER_C)
dap_cv = metrics(y, np.exp(b[0] + b[1] * inner + b[2] * male))
print(f'DAP {len(y)} adults; lambda {lam}; stack log kg = {b[0]:.3f} + {b[1]:.3f} log breed + {b[2]:.3f} male')
print(fmt('DAP 5-fold (OOF breed term)', dap_cv))

C = [r for r in json.load(open(f'{build}/customers.json')) if r['age'] is not None and r['age'] >= MIN_AGE and r['comp']]
lbp = np.array([sum(p * w[bi[k]] for k, p in r['comp'].items() if k in bi) / sum(p for k, p in r['comp'].items() if k in bi) for r in C])
mc = np.array([0.5 if r['male'] is None else r['male'] for r in C]); yc = np.array([r['y_kg'] for r in C])
cust = metrics(yc, np.exp(b[0] + b[1] * lbp + b[2] * mc))
print(fmt('customers (external)', cust))

name = 'DAP_breed_sex'
model = {
    'model': name, 'built': datetime.date.today().isoformat(),
    'n_train': len(y), 'min_age_years': MIN_AGE,
    'truth': 'Dog Aging Project owner-reported weights (HLES dd_weight_lbs), genotyped dogs >= 1 y',
    'training': {'source': 'Dog Aging Project', 'n_dogs': len(y),
                 'composition': 'stage-9 estimator (v2.3 panel, lasso 0.3) on GLIMPSE-imputed DAP genotypes',
                 'dap_internal_cv': {k: dap_cv[k] for k in ('r', 'mae', 'bias', 'slope')}},
    'cv': {'design': f'external: trained on {len(y)} Dog Aging Project adults, scored on {len(yc)} ProsperK9 '
                     f'customer adults with owner-reported weights, none used in fitting',
           'summary': {name: cust}},
    'space': 'log', 'lambda': lam, 'box': BOX,
    'prior_source': f'AKC breed-standard table for {int(np.isfinite(prior_akc).sum())} labels, '
                    'proportion-weighted marginal mean for the rest',
    'breeds': breeds,
    'breed_weight': {k: round(float(np.exp(w[i])), 3) for i, k in enumerate(breeds)},
    'breed_prior': {k: round(float(np.exp(lpr[i])), 3) for i, k in enumerate(breeds)},
    'default_kg': round(float(np.exp(np.median(ly))), 3),
    'stack': {'space': 'log', 'intercept': float(b[0]), 'coef_log_breed_pred': float(b[1]),
              'covariates': {'male': float(b[2])}, 'robust': f'Huber c={LOG_HUBER_C} in log space',
              'formula': 'pred_kg = exp(intercept + coef_log_breed_pred*log(breed_pred_kg) + sum(coef*cov)); '
                         'male is 1/0 from the X-coverage sex call, 0.5 if unknown',
              'male_multiplier': round(float(np.exp(b[2])), 4),
              'n_male': int((male == 1).sum()), 'n_female': int((male == 0).sum())},
    'previous': {'model': 'M3BL_da_sex', 'n_train': 577, 'note': 'customer-trained, 231 labels, Darwin\'s Ark size PRS; '
                 'on the same customers: r 0.754, MAE 5.43 kg in-sample, 5.56 kg honest CV'},
}
json.dump(model, open(out_json, 'w'), indent=1)
print(f'wrote {out_json}: {len(breeds)} labels')
