#!/usr/bin/env python3
"""
Does Cornell (Hayward 2016) train a better weight model than the Dog Aging Project?
Same recipe as fit_dap.py (box-constrained ridge breed weights, log-space Huber stack),
breed + sex only (the DA PRS can't be computed on the array and adds nothing in DAP).
Every candidate is fit without customers and scored on the same customer adults.

    fit_cornell.py CORNELL_BUILD DAP_BUILD run_dog_pipeline.sh
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from fit_weight_model import (extract_dict, join_prior, metrics, fmt, fit_box_ridge, box_cv_lambda,
                              box_inner_oof, marginal_prior, huber_ols, LOG_HUBER_C)

cb, db, pipeline_sh = sys.argv[1:4]
zc, zd = np.load(f'{cb}/cornell.npz'), np.load(f'{db}/dap.npz')
breeds = [str(b) for b in zd['breeds']]
assert breeds == [str(b) for b in zc['breeds']]
bi = {b: i for i, b in enumerate(breeds)}

def load(z, pheno, wcol, mcol, extra):
    qi = {str(d): i for i, d in enumerate(z['ids'])}
    rows = [l.rstrip('\n').split('\t') for l in open(pheno)]
    h = rows[0]; out = []
    for r in rows[1:]:
        d = dict(zip(h, r))
        if d['dog_id'] in qi and d[mcol] != '' and extra(d):
            out.append((qi[d['dog_id']], float(d[wcol]), float(d[mcol]), d.get('breed', ''), d.get('breed_status', '')))
    X = z['Q'][[o[0] for o in out]].astype(np.float64); X /= X.sum(1, keepdims=True)
    return X, np.log([o[1] for o in out]), np.array([o[2] for o in out]), [o[3] for o in out], [o[4] for o in out]

Xk, lyk, mk, labk, _ = load(zc, f'{cb}/cornell_pheno.tsv', 'weight_kg', 'male', lambda d: True)
Xd, lyd, md, labd, statd = load(zd, f'{db}/dap_pheno.tsv', 'weight_kg', 'male', lambda d: d['age'] != '' and float(d['age']) >= 1)
print(f'Cornell {len(lyk)} weighed dogs (median {np.exp(np.median(lyk)):.1f} kg); DAP {len(lyd)} adults (median {np.exp(np.median(lyd)):.1f} kg)')

# Cornell breed-call sanity: top population vs the study label
top = [breeds[i] for i in Xk.argmax(1)]
mx = Xk.max(1)
norm = lambda s: s.upper().replace(' ', '_')
lab_ok = [(norm(l) in norm(t) or norm(t) in norm(l)) for l, t in zip(labk, top)]
village = np.array(['village' in l for l in labk])
print(f'Cornell: {np.mean(mx >= 0.85):.0%} have one population >= 85%; top call matches the study label '
      f'(substring) for {np.mean(lab_ok):.0%}; village dogs {village.sum()}; '
      f'"mixed"-labelled {sum("mix" in l for l in labk)}')

C = [r for r in json.load(open(f'{db}/customers.json')) if r['age'] is not None and r['age'] >= 1 and r['comp'] and r['da'] is not None]
Xc = np.zeros((len(C), len(breeds)))
for i, r in enumerate(C):
    for b, p in r['comp'].items():
        if b in bi: Xc[i, bi[b]] = p
Xc /= Xc.sum(1, keepdims=True)
yc = np.array([r['y_kg'] for r in C]); mc = np.array([0.5 if r['male'] is None else r['male'] for r in C])
print(f'Customers: {len(yc)} adults\n')

akc = join_prior(breeds, extract_dict(open(pipeline_sh).read(), 'BREED_WEIGHT_KG'))
prior_akc = np.array([akc.get(b, np.nan) for b in breeds])

def train(X, ly, male, box=1.5, lam=None):
    lpr = np.log(marginal_prior(X, np.exp(ly), prior_akc))
    lam = box_cv_lambda(X, ly, lpr, 0) if lam is None else lam
    w = fit_box_ridge(X, ly, lpr, lam, box=box)
    inner = np.log(box_inner_oof(X, ly, lpr, lam, 0, box=box) if 'box' in box_inner_oof.__code__.co_varnames
                   else box_inner_oof(X, ly, lpr, lam, 0))
    b = huber_ols(np.column_stack([inner, male]), ly, c=LOG_HUBER_C)
    return w, b, lam

def score(name, w, b, X=Xc, y=yc, male=mc):
    pred = np.exp(b[0] + b[1] * (X @ w) + b[2] * male)
    m = metrics(y, pred); print(fmt(name, m)); return pred, m

res, preds = {}, {}
for box in (1.5, 2.0):
    res[('DAP', box)] = train(Xd, lyd, md, box)
    res[('CORNELL', box)] = train(Xk, lyk, mk, box)
    res[('DAP+CORNELL', box)] = train(np.vstack([Xd, Xk]), np.concatenate([lyd, lyk]), np.concatenate([md, mk]), box)
print('Customer adults, models trained without any customer (breed + sex):')
for (src, box), (w, b, lam) in res.items():
    preds[(src, box)], _ = score(f'{src:<12} box {box} (lam {lam}, male x{np.exp(b[2]):.2f})', w, b)

# paired bootstrap, Cornell and pooled vs DAP at the same box
rng = np.random.RandomState(7); bs = rng.randint(0, len(yc), size=(2000, len(yc)))
print('\nMAE difference vs DAP at the same box (kg, paired bootstrap 95% CI):')
for box in (1.5, 2.0):
    ref = np.abs(preds[('DAP', box)] - yc)
    for src in ('CORNELL', 'DAP+CORNELL'):
        d = (np.abs(preds[(src, box)] - yc) - ref)[bs].mean(1)
        print(f'  {src:<12} box {box}: {d.mean():+.2f} [{np.percentile(d, 2.5):+.2f}, {np.percentile(d, 97.5):+.2f}]')

# cross-population: each source predicting the other
print('\nCross-checks:')
w, b, _ = res[('CORNELL', 1.5)]; score('Cornell model -> DAP adults', w, b, Xd, np.exp(lyd), md)
w, b, _ = res[('DAP', 1.5)]; score('DAP model -> Cornell dogs', w, b, Xk, np.exp(lyk), mk)

# purebred medians in each population
print('\nPurebred median weight (kg), dogs with one population >= 85%: Cornell (n) vs DAP (n), fitted 100% weight Cornell / DAP')
def pure(X, ly):
    d = {}
    for i in np.where(X.max(1) >= 0.85)[0]:
        d.setdefault(breeds[int(X[i].argmax())], []).append(np.exp(ly[i]))
    return d
pk, pdp = pure(Xk, lyk), pure(Xd, lyd)
wk, wd = res[('CORNELL', 1.5)][0], res[('DAP', 1.5)][0]
for br in sorted(pk, key=lambda b: -len(pk[b]))[:30]:
    dm = f'{np.median(pdp[br]):5.1f} ({len(pdp[br]):>3})' if br in pdp else '    -      '
    print(f'  {br:<32} {np.median(pk[br]):5.1f} ({len(pk[br]):>3})  {dm}   fit {np.exp(wk[bi[br]]):5.1f} / {np.exp(wd[bi[br]]):5.1f}')
nk = sum(1 for j in range(len(breeds)) if (Xk[:, j] > 0.25).sum() >= 3)
nd = sum(1 for j in range(len(breeds)) if (Xd[:, j] > 0.25).sum() >= 3)
print(f'\nLabels with >= 3 dogs at > 25% ancestry: Cornell {nk}, DAP {nd} (of {len(breeds)})')
