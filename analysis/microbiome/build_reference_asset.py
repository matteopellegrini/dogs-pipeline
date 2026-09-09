#!/usr/bin/env python3
"""Build the frontend microbiome reference asset from the panel-of-normals.

For every species seen in >=2% of the 1,254 panel dogs, ship:
  - q: 101 abundance quantiles (%) across the panel dogs that CARRY the
       species (zeros excluded). Zeros-included quantiles made every
       deeply-characterized dog score >=50th everywhere: detecting a
       species at all already beats the whole zero block, so percentile
       mostly measured detection, which tracks sequencing depth. Against
       carriers only, the percentile measures LEVEL given presence;
       absence renders as its own state in the report, not a number;
  - rho: Spearman correlation of abundance vs age across aged panel dogs,
       the display's red/green "increases/decreases with age" direction —
       model-independent, so it never disagrees with a per-report retrain;
  - prev: fraction of panel dogs carrying the species at all.

Output: dogs-app/lib/microbiomeReference.json (static asset; the report
frontend keys into it by full clade string).
"""
import json
import numpy as np

PANEL = 'reference_panel/microbiome_panel.json'
OUT = 'dogs-app/lib/microbiomeReference.json'
MIN_PREV = 0.02

p = json.load(open(PANEL))
dogs = p['dogs']
n = len(dogs)
ages = np.array([d.get('age') if d.get('age') is not None else np.nan for d in dogs])
aged = ~np.isnan(ages)
print(f'panel dogs: {n} (aged: {int(aged.sum())})')

clades = {}
for d in dogs:
    for c in d['species']:
        clades[c] = clades.get(c, 0) + 1

def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    # midranks for ties (zeros dominate) — without this the sign is unstable
    for arr, r in ((x, rx), (y, ry)):
        vals, inv, cnt = np.unique(arr, return_inverse=True, return_counts=True)
        sums = np.zeros(len(vals)); np.add.at(sums, inv, r)
        r[:] = (sums / cnt)[inv]
    rx -= rx.mean(); ry -= ry.mean()
    den = np.sqrt((rx**2).sum() * (ry**2).sum())
    return float((rx*ry).sum() / den) if den > 0 else 0.0

asset = {}
for c, cnt in clades.items():
    prev = cnt / n
    if prev < MIN_PREV:
        continue
    v = np.array([d['species'].get(c, 0.0) for d in dogs])
    q = np.percentile(v[v > 0], np.arange(101))
    rho = spearman(v[aged], ages[aged])
    asset[c] = {'q': [round(float(x), 4) for x in q],
                'rho': round(rho, 3),
                'prev': round(prev, 3)}

# Prevalence for EVERY species ever seen in the panel (no floor) — the
# rare-species card needs "seen in 0.3% of dogs" for taxa far below the
# 2% quantile floor; species absent even from this map were never seen in
# any panel dog. Tiny: one rounded float per clade.
prev_all = {c: round(cnt / n, 4) for c, cnt in clades.items()}

json.dump({'meta': {'n_dogs': n, 'n_aged': int(aged.sum()),
                    'n_species': len(asset), 'min_prevalence': MIN_PREV,
                    'built': p['meta'].get('built'), 'source_version': p['meta'].get('version')},
           'species': asset, 'prev_all': prev_all}, open(OUT, 'w'))
print(f'{len(asset)} species (+{len(prev_all)} prevalence-only) -> {OUT}')
