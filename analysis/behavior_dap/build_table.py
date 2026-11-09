#!/usr/bin/env python3
"""
Dog Aging Project behaviour table for the breed-composition behaviour models.

    build_table.py DAP_WEIGHT_DIR PANEL_DIR LMM_BETAS DOGS_TSV OUT.npz

DAP_WEIGHT_DIR is the weight build dir (build/dap.npz with Q from the stage-9 estimator, and
parts/chr*.npz with panel-site dosages of a1). Adds, for every genotyped DAP dog:
  - the 42 HLES behaviour items (db_*, 0-4; codes >= 98 missing), activity level,
    activity intensity and appetite;
  - covariates: male, neutered, age at HLES, mixed-breed flag, owner breed label;
  - today's report scores for the 14 AKC traits (dense GEMMA LMM PRS on the panel
    sites, z against the stored reference, predicted_score = clip(mean + z*sd, 1, 5)),
    exactly as stage 11 computes them.
"""
import glob, gzip, json, sys
import numpy as np
import pandas as pd

wdir, panel_dir, lmm_path, dogs_tsv, out = sys.argv[1:6]
z = np.load(f'{wdir}/build/dap.npz')
ids = [str(i) for i in z['ids']]
Q = z['Q'].astype(np.float32)
breeds = [str(b) for b in z['breeds']]

# panel dosages (a1 x 100, 255 missing)
sites = pd.read_csv(f'{panel_dir}/sites.tsv', sep='\t')
X = np.full((len(sites), len(ids)), 255, dtype=np.uint8)
for f in sorted(glob.glob(f'{wdir}/parts/chr*.npz')):
    p = np.load(f)
    X[p['panel_idx']] = p['panel_dos']
key = {(c, int(pos)): (i, a1) for i, (c, pos, a1) in enumerate(zip(sites.chrom, sites.pos, sites.a1))}

L = json.load(gzip.open(lmm_path, 'rt'))
row, flip, af = [], [], []
for c, pos, ea, oa, a in L['sites']:
    k = key.get((c, int(pos)))
    row.append(-1 if k is None else k[0]); flip.append(k is not None and k[1] != ea); af.append(a)
row, flip, af = np.array(row), np.array(flip), np.array(af)
ok = row >= 0
D = np.repeat((2 * af)[:, None], len(ids), axis=1).astype(np.float32)   # missing -> 2af
xs = X[row[ok]].astype(np.float32)
miss = xs == 255
xs = xs / 100.0
xs[flip[ok]] = 2.0 - xs[flip[ok]]
xs[miss] = np.broadcast_to((2 * af[ok])[:, None], xs.shape)[miss]
D[ok] = xs
print(f'LMM sites {len(row)}, on panel {ok.sum()}, flipped {flip.sum()}, missing dosages {int(miss.sum())}')
akc_traits = [t for t in L['traits'] if t not in ('height_cm', 'weight_kg')]
akc_z, akc_score = np.zeros((len(ids), len(akc_traits))), np.zeros((len(ids), len(akc_traits)))
for j, t in enumerate(akc_traits):
    ti = L['traits'][t]
    beta = np.array(json.loads(ti['beta']) if isinstance(ti['beta'], str) else ti['beta'], dtype=np.float64)
    ref = np.array(json.loads(ti['ref_prs_sorted']) if isinstance(ti['ref_prs_sorted'], str) else ti['ref_prs_sorted'])
    prs = beta @ D
    zz = (prs - ref.mean()) / (ref.std() + 1e-8)
    akc_z[:, j] = zz
    akc_score[:, j] = np.clip(float(ti['truth_mean']) + zz * float(ti['truth_sd']), 1, 5)
    print(f'  {t}: score mean {akc_score[:, j].mean():.2f} sd {akc_score[:, j].std():.2f}')

dogs = pd.read_csv(dogs_tsv, sep='\t', low_memory=False, dtype={'dog_id': str}).set_index('dog_id').reindex(ids)
num = lambda c: pd.to_numeric(dogs[c], errors='coerce')
items = [c for c in dogs.columns if c.startswith('db_')] + ['pa_activity_level', 'pa_avg_activity_intensity', 'df_appetite']
Y = np.column_stack([num(c).where(num(c) < 98).values for c in items]).astype(np.float32)
sex = dogs['Sex_Class_at_HLES'].fillna('')
male = np.where(sex == '', np.nan, sex.str.startswith('Male').astype(float))
neut = np.where(sex == '', np.nan, sex.str.contains('spayed|neutered').astype(float))
age = num('Estimated_Age_Years_at_HLES').values
mixed = dogs['Breed_Status'].eq('Mixed breed').values
label = dogs['Breed'].fillna('').astype(str).values if 'Breed' in dogs else np.array([''] * len(ids))
print(f'{len(ids)} dogs; {np.isfinite(Y).any(1).sum()} with any behaviour answer; mixed {mixed.sum()}')
np.savez_compressed(out, ids=np.array(ids), breeds=np.array(breeds), Q=Q, items=np.array(items), Y=Y,
                    male=male, neutered=neut, age=age, mixed=mixed, label=label,
                    akc_traits=np.array(akc_traits), akc_z=akc_z, akc_score=akc_score)
print(f'wrote {out}')
