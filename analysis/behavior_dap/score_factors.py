#!/usr/bin/env python3
"""
Darwin's Ark behaviour factors for DAP dogs, as stage 11 reports them: honest z from
the per-chromosome PRS parts (extract_factors.py), then the ancestry blend
(w0 + w_prs_z * z + w_ancestry * breed-composition expectation), NaN where the
expectation covers < 50% of the dog's ancestry (stage 11 then shows the raw percentile).

    score_factors.py PARTS_DIR MANIFEST.json.gz TABLE.npz OUT.npz
"""
import ast, glob, gzip, json, sys
import numpy as np

parts, manifest, tab, out = sys.argv[1:5]
M = json.load(gzip.open(manifest, 'rt'))
t = np.load(tab, allow_pickle=True)
Q = t['Q'].astype(np.float64); breeds = [str(b) for b in t['breeds']]
prs, traits = None, None
for f in sorted(glob.glob(f'{parts}/chr*.npz')):
    p = np.load(f)
    traits = [str(x) for x in p['traits']]
    prs = p['prs'] if prs is None else prs + p['prs']
print(f'{len(glob.glob(parts + "/chr*.npz"))} chromosomes, {prs.shape[1]} dogs')
means = M['breed_trait_means']
num = lambda v: ast.literal_eval(v) if isinstance(v, str) else v
Z, B, A = [], [], []
for k, tr in enumerate(traits):
    info = M['traits'][tr]
    ref = np.array(num(info['ref_prs_sorted']), dtype=float)
    z = (prs[k] - ref.mean()) / ((ref.std() + 1e-9) * float(info['z_scale']))
    v = np.array([means.get(b, {}).get(tr, np.nan) for b in breeds])
    cov = Q[:, np.isfinite(v)].sum(1)
    anc = np.where(cov >= 0.5, (Q[:, np.isfinite(v)] @ v[np.isfinite(v)]) / np.maximum(cov, 1e-9), np.nan)
    bl = num(info['ancestry_blend'])
    blended = bl['w0'] + bl['w_prs_z'] * z + bl['w_ancestry'] * anc
    Z.append(z); A.append(anc); B.append(np.where(np.isfinite(blended), blended, np.nan))
    print(f'  {tr}: z sd {z.std():.2f}, ancestry covered {np.isfinite(anc).mean():.0%}')
np.savez_compressed(out, traits=np.array(traits), z=np.array(Z), ancestry=np.array(A), blended=np.array(B))
print(f'wrote {out}')
