"""Validate the DAP relatives block before it joins the reference.
Same-breed background must stay under the second-degree threshold (0.0884) the
way it does for Dog10K + cohort dogs; DAP dogs are low-pass imputed, so the
check that matters is whether imputation inflates kinship between unrelated
dogs of one breed. KING-robust exactly as stage 9b."""
import gzip, json, collections
import numpy as np
REF = '/u/project/pellegrini/gkislik/dogs/relatives_ref'
rng = np.random.default_rng(7)

def king(q, R):
    use = q >= 0; Rs, qs = R[use], q[use]
    valid = Rs >= 0
    het, hA, hB = qs == 1, qs == 0, qs == 2
    Nhh = (Rs[het] == 1).sum(0); Nopp = (Rs[hA] == 2).sum(0) + (Rs[hB] == 0).sum(0)
    den = valid[het].sum(0) + (Rs == 1).sum(0)
    return np.where(den > 0, (Nhh - 2.0 * Nopp) / np.maximum(den, 1), -1.0)

B = np.load(f'{REF}/geno_dap.npy'); dm = json.load(open(f'{REF}/meta_dap.json'))['samples']
G = np.load(f'{REF}/geno.npy', mmap_mode='r'); gm = json.load(gzip.open(f'{REF}/meta.json.gz', 'rt'))['samples']
dl = np.array([s['label'] for s in dm])
print('DAP missing per dog: median', float(np.median((B < 0).mean(0))), 'max', float((B < 0).mean(0).max()))
print('DAP het rate median', float(np.median((B == 1).sum(0) / (B >= 0).sum(0))))

pure = [j for j, s in enumerate(dm) if not s['label'].endswith(' mix')]
top = [b for b, n in collections.Counter(dl[pure]).most_common(25)]
qs = [j for b in top for j in rng.choice([j for j in pure if dl[j] == b], size=min(12, sum(dl[pure] == b)), replace=False)]
same, cross, dup, hits = [], [], 0, collections.Counter()
for j in qs:
    phi = king(B[:, j], B); phi[j] = -1
    sb = dl == dl[j]
    same += list(phi[sb]); cross += list(phi[~sb][rng.choice((~sb).sum(), 300, replace=False)])
    dup += int((phi >= 0.45).sum())
    hits['first'] += int(((phi >= 0.177) & (phi < 0.45)).sum()); hits['second'] += int(((phi >= 0.0884) & (phi < 0.177)).sum())
same, cross = np.array(same), np.array(cross)
print(f'DAP->DAP purebred queries {len(qs)} ({len(top)} breeds)')
print(f'  same-breed pairs n={len(same)} p50 {np.percentile(same,50):.3f} p99 {np.percentile(same,99):.3f} p99.9 {np.percentile(same,99.9):.3f}  >=0.0884 {np.mean(same>=0.0884):.4f}')
print(f'  cross-breed p99 {np.percentile(cross,99):.3f}; same-dog pairs {dup}; first-degree {hits["first"]}; second-degree {hits["second"]}')

# Reference background for comparison: Dog10K purebreds vs Dog10K, same breeds
gl = np.array([s['label'] for s in gm]); src = np.array([s['source'] for s in gm])
d10 = [j for j in range(len(gm)) if src[j] == 'dog10k' and gl[j] in top]
G10 = np.asarray(G[:, src == 'dog10k']); l10 = gl[src == 'dog10k']
s10 = []
for j in rng.choice(d10, size=min(150, len(d10)), replace=False):
    phi = king(np.asarray(G[:, j]), G10); s10 += list(phi[(l10 == gl[j])][phi[(l10 == gl[j])] < 0.45])
s10 = np.array(s10)
print(f'Dog10K->Dog10K same-breed (for reference) n={len(s10)} p99 {np.percentile(s10,99):.3f}  >=0.0884 {np.mean(s10>=0.0884):.4f}')

# Customer dogs (already in the reference) as queries against the DAP block
cust = [j for j in range(len(gm)) if src[j] == 'prosperkits']
cnt = collections.Counter(); top_phi = []
for j in rng.choice(cust, size=400, replace=False):
    phi = king(np.asarray(G[:, j]), B)
    m = phi.max(); top_phi.append(m)
    cnt['same_dog'] += int((phi >= 0.45).sum()); cnt['first'] += int(((phi >= 0.177) & (phi < 0.45)).sum())
    cnt['second'] += int(((phi >= 0.0884) & (phi < 0.177)).sum())
    cnt['dogs_with_any'] += int((phi >= 0.0884).any())
print('400 customers vs DAP:', dict(cnt), 'max phi p50/p99', np.percentile(top_phi, [50, 99]).round(3).tolist())
