"""Why DAP dogs look related to everyone: dosage distribution and how the
customer-vs-DAP kinship background moves with a stricter hard-call rule
(dosages farther than `tol` from 0/1/2 treated as missing)."""
import gzip, json, os
import numpy as np
D = '/u/project/pellegrini/gkislik/dogs'; REF, W = f'{D}/relatives_ref', f'{D}/dap_weight'
rng = np.random.default_rng(3)
keys = {}
for i, line in enumerate(open(f'{REF}/keys.tsv')):
    c, p, ref, alt = line.rstrip('\n').split('\t'); keys[(c, int(p))] = (i, ref, alt)
panel = [l.rstrip('\n').split('\t') for l in open(f'{D}/breed_panel_v23/sites.tsv')][1:]
X = np.full((len(keys), 7828), np.nan, dtype=np.float32)
for ch in sorted({k[0] for k in keys}):
    z = np.load(f'{W}/parts/{ch}.npz'); pidx, pdos = z['panel_idx'], z['panel_dos']
    for r, pi in enumerate(pidx):
        c, pos, a1, a2 = panel[pi]; k = keys.get((c, int(pos)))
        if k is None: continue
        i, ref, alt = k; d = pdos[r].astype(np.float32); m = d == 255; d /= 100
        x = d if (a1, a2) == (alt, ref) else (2 - d if (a1, a2) == (ref, alt) else None)
        if x is None: continue
        x[m] = np.nan; X[i] = x
frac = lambda lo, hi: float(np.nanmean((X > lo) & (X < hi)))
print('DAP dosage: within 0.1 of 0/1/2:', round(float(np.nanmean(np.minimum(np.minimum(np.abs(X), np.abs(X - 1)), np.abs(X - 2)) <= 0.1)), 3),
      ' in (0.3,0.7):', round(frac(0.3, 0.7), 3), ' in (1.3,1.7):', round(frac(1.3, 1.7), 3))

G = np.load(f'{REF}/geno.npy', mmap_mode='r'); gm = json.load(gzip.open(f'{REF}/meta.json.gz', 'rt'))['samples']
src = np.array([s['source'] for s in gm])
def king(q, R):
    use = q >= 0; Rs, qs = R[use], q[use]; valid = Rs >= 0
    het, hA, hB = qs == 1, qs == 0, qs == 2
    Nhh = (Rs[het] == 1).sum(0); Nopp = (Rs[hA] == 2).sum(0) + (Rs[hB] == 0).sum(0)
    den = valid[het].sum(0) + (Rs == 1).sum(0)
    return np.where(den > 0, (Nhh - 2.0 * Nopp) / np.maximum(den, 1), -1.0)
for name, idx in (('dog10k', np.where(src == 'dog10k')[0]), ('prosperkits', np.where(src == 'prosperkits')[0])):
    S = np.asarray(G[:, rng.choice(idx, 200, replace=False)]); print(name, 'het rate median', round(float(np.median((S == 1).sum(0) / (S >= 0).sum(0))), 3))
cust = rng.choice(np.where(src == 'prosperkits')[0], 100, replace=False)
sub = rng.choice(7828, 2000, replace=False); Xs = X[:, sub]
for tol in (0.5, 0.3, 0.2, 0.1):
    r = np.rint(Xs); B = np.where(np.isnan(Xs) | (np.abs(Xs - r) > tol), -1, r).astype(np.int8)
    het = float(np.median((B == 1).sum(0) / np.maximum((B >= 0).sum(0), 1))); miss = float((B < 0).mean())
    mx = []; n2 = 0
    for j in cust:
        phi = king(np.asarray(G[:, j]), B); mx.append(phi.max()); n2 += int((phi >= 0.0884).sum())
    print(f'tol {tol}: missing {miss:.2f} het {het:.3f} | 100 customers vs 2000 DAP: max phi p50 {np.median(mx):.3f} p95 {np.percentile(mx,95):.3f}, pairs >=0.0884: {n2}')
