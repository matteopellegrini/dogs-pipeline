"""Which DAP dogs produce the customer 'relatives': hit counts per DAP dog vs
its heterozygosity, and the customer background after dropping outliers."""
import gzip, json
import numpy as np
REF = '/u/project/pellegrini/gkislik/dogs/relatives_ref'
rng = np.random.default_rng(11)
B = np.load(f'{REF}/geno_dap.npy'); dm = json.load(open(f'{REF}/meta_dap.json'))['samples']
G = np.load(f'{REF}/geno.npy', mmap_mode='r'); gm = json.load(gzip.open(f'{REF}/meta.json.gz', 'rt'))['samples']
src = np.array([s['source'] for s in gm])
def king(q, R):
    use = q >= 0; Rs, qs = R[use], q[use]; valid = Rs >= 0
    het, hA, hB = qs == 1, qs == 0, qs == 2
    Nhh = (Rs[het] == 1).sum(0); Nopp = (Rs[hA] == 2).sum(0) + (Rs[hB] == 0).sum(0)
    den = valid[het].sum(0) + (Rs == 1).sum(0)
    return np.where(den > 0, (Nhh - 2.0 * Nopp) / np.maximum(den, 1), -1.0)
het = (B == 1).sum(0) / (B >= 0).sum(0)
print('DAP het rate percentiles 1/25/50/75/99/99.9:', np.percentile(het, [1, 25, 50, 75, 99, 99.9]).round(3).tolist())
cust = rng.choice(np.where(src == 'prosperkits')[0], 150, replace=False)
hits = np.zeros(B.shape[1], int); PHI = []
for j in cust:
    phi = king(np.asarray(G[:, j]), B); hits += phi >= 0.0884; PHI.append(phi)
PHI = np.array(PHI)
order = np.argsort(hits)[::-1]
print('DAP dogs hit by >=10% of customers:', int((hits >= 15).sum()), ' >=50%:', int((hits >= 75).sum()))
for j in order[:12]: print(f'  {dm[j]["id"]:10s} {dm[j]["label"][:28]:28s} het {het[j]:.3f} hits {hits[j]}/150')
for cut in (0.40, 0.38, 0.36, 0.35):
    keep = het < cut
    mx = PHI[:, keep].max(1)
    print(f'drop het >= {cut} ({int((~keep).sum())} dogs): customers with any >=0.0884: {int((PHI[:, keep] >= 0.0884).any(1).sum())}/150, max phi p50 {np.median(mx):.3f}, total pairs {int((PHI[:, keep] >= 0.0884).sum())}')
# same comparison against the existing customer block, for scale
C = np.asarray(G[:, src == 'prosperkits'])
cc = [ (lambda p: (p[(p < 0.45)] >= 0.0884).sum())(king(np.asarray(G[:, j]), C)) for j in cust[:60]]
print('existing reference: 60 customers vs 1,773 customer dogs, pairs >=0.0884 (excluding self):', int(np.sum(cc)))
