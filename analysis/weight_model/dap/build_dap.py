#!/usr/bin/env python3
"""
Assemble the weight-model training table from the Dog Aging Project extract
(extract.sbatch) plus the ProsperK9 customer table, so the fit can run anywhere.

    build_dap.py W_DIR PANEL_DIR DAP_DOGS_TSV KITS_LIST CUST_WEIGHTS CUST_AGE OUT_DIR

DAP dogs get a breed composition from their panel-site dosages with the stage-9
estimator exactly as a customer does (NNLS on P with the sum-to-one row and the
non-negative lasso, lambda 0.3), and the Darwin's Ark size PRS from
extract_chr.py. Writes:
  dap.npz        ids, Q (dogs x populations), da_prs, sites_used
  dap_pheno.tsv  dog_id, weight_kg, age, male, neutered, breed_status, breed
  customers.json one row per weighed kit: kit, y_kg, age, male, da, pred_now, comp
"""
import csv, glob, json, os, sys
import numpy as np
from scipy.optimize import nnls
from scipy.linalg import solve_triangular

LASSO = 0.3
w_dir, panel_dir, dap_tsv, kits_list, cust_w, cust_age, out_dir = sys.argv[1:8]
os.makedirs(out_dir, exist_ok=True)

ids = [l.strip() for l in open(f'{w_dir}/samples.txt') if l.strip()]
P = np.load(f'{panel_dir}/phat.npy').astype(np.float32)
breeds = [l.strip() for l in open(f'{panel_dir}/breeds.txt') if l.strip()]
n_sites, n = P.shape[0], len(ids)
X = np.full((n_sites, n), 255, dtype=np.uint8)
da = np.zeros(n)
matched = total = 0
for f in sorted(glob.glob(f'{w_dir}/parts/chr*.npz')):
    z = np.load(f)
    if z['panel_dos'].shape[1] != n:
        sys.exit(f'{f}: {z["panel_dos"].shape[1]} dogs, expected {n}')
    X[z['panel_idx']] = z['panel_dos']
    da += z['da_prs']; matched += int(z['da_matched']); total += int(z['da_total'])
print(f'{len(glob.glob(w_dir + "/parts/chr*.npz"))} chromosomes; DA sites matched {matched}/{total}')

miss = X == 255
site_ok = ~miss.all(1)
print(f'panel sites present {site_ok.sum()}/{n_sites}; per-dog extra missing: '
      f'max {int((miss[site_ok]).sum(0).max())}')

# stage-9 estimator on the shared site set; dogs with their own missing sites refit alone
Pv = P[site_ok]
A_ones = np.ones((1, Pv.shape[1]), dtype=np.float32)
pen = 0.5 * LASSO * float(np.mean(np.einsum('ij,ij->j', Pv, Pv)))
G = (Pv.T.astype(np.float64) @ Pv) + 1.0      # A'A with the appended ones row
G[np.diag_indices_from(G)] += 1e-6
R = np.linalg.cholesky(G).T
Q = np.zeros((n, len(breeds)), dtype=np.float32)
used = np.full(n, int(site_ok.sum()))
for s in range(0, n, 500):
    blk = X[site_ok, s:s + 500]
    xb = blk.astype(np.float32) / 100.0
    bad = blk == 255
    for j in range(blk.shape[1]):
        if bad[:, j].any():                        # refit this dog on its own sites
            v = ~bad[:, j]
            Pd = Pv[v]
            Gd = Pd.T.astype(np.float64) @ Pd + 1.0
            Gd[np.diag_indices_from(Gd)] += 1e-6
            Rd = np.linalg.cholesky(Gd).T
            rhs = Pd.T.astype(np.float64) @ xb[v, j] + 1.0 \
                - 0.5 * LASSO * float(np.mean(np.einsum('ij,ij->j', Pd, Pd)))
            q, _ = nnls(Rd, solve_triangular(Rd.T, rhs, lower=True))
            used[s + j] = int(v.sum())
        else:
            rhs = Pv.T.astype(np.float64) @ xb[:, j] + 1.0 - pen
            q, _ = nnls(R, solve_triangular(R.T, rhs, lower=True))
        Q[s + j] = q / (q.sum() + 1e-12)
    print(f'  scored {min(s + 500, n)}/{n}', flush=True)
np.savez_compressed(f'{out_dir}/dap.npz', ids=np.array(ids), breeds=np.array(breeds), Q=Q, da_prs=da,
                    sites_used=used)

# DAP phenotypes
with open(dap_tsv) as fh, open(f'{out_dir}/dap_pheno.tsv', 'w') as out:
    out.write('dog_id\tweight_kg\tage\tmale\tneutered\tbreed_status\tbreed\n')
    for r in csv.DictReader(fh, delimiter='\t'):
        try:
            w = float(r['dd_weight_lbs']) / 2.20462
        except ValueError:
            continue
        out.write('\t'.join([r['dog_id'], f'{w:.3f}', r['Estimated_Age_Years_at_HLES'],
                             {'1.0': '1', '2.0': '0'}.get(r['dd_sex'], ''),
                             {'True': '1', 'False': '0'}.get(r['dd_spayed_or_neutered'], ''),
                             r['Breed_Status'], r['Breed']]) + '\n')

# customers (kit dirs from the republish list: pk, barcode, pub, work)
wlb = dict(l.split('\t')[:2] for l in open(cust_w) if '\t' in l)
age = dict(l.split('\t')[:2] for l in open(cust_age) if '\t' in l)
rows = []
for line in open(kits_list):
    f = line.rstrip('\n').split('\t')
    kit, pub = f[1], f[2]
    if kit not in wlb:
        continue
    try:
        br = json.load(open(f'{pub}/breed_result.json'))
        wk = json.load(open(f'{pub}/prs_result.json'))['physical_traits']['weight_kg']
    except (OSError, ValueError, KeyError):
        continue
    sex = None
    try:
        sex = json.load(open(f'{pub}/coverage_1mb.json')).get('_meta', {}).get('predicted_sex')
    except (OSError, ValueError):
        pass
    rows.append({'kit': kit, 'y_kg': float(wlb[kit]) / 2.20462,
                 'age': float(age[kit]) if kit in age else None,
                 'male': {'male': 1.0, 'female': 0.0}.get(sex),
                 'da': wk.get('da_size_prs'), 'pred_now': wk.get('pred_kg'),
                 'comp': {r['breed']: float(r['proportion']) for r in br.get('breed_composition_raw') or []}})
json.dump(rows, open(f'{out_dir}/customers.json', 'w'))
print(f'{len(rows)} customers with weight + results')
