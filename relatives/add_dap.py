"""Relatives reference v4: append the Dog Aging Project dogs (2026-11-09).

DAP genotypes are the per-chromosome panel dosages already extracted for the
weight model (dap_weight/parts/chr*.npz: 100 x dosage of breed-panel allele a1,
255 = missing; dap_weight/samples.txt = column order). Every relatives site is
a v2.3 panel site, so each is re-expressed as the ALT count of keys.tsv and
rounded exactly as stage 9b rounds a query (<0.5, <1.5, else 2).

Writes STAGING files next to the live reference (relatives_ref/geno_dap.npy,
meta_dap.json); assemble with --commit after validate_dap.py.
  python3 relatives/add_dap.py            # build the DAP block
  python3 relatives/add_dap.py --commit   # geno.npy/meta.json.gz := v3 + DAP (v3 kept as *_v3_backup)"""
import gzip, json, os, sys
import numpy as np

D = '/u/project/pellegrini/gkislik/dogs'
REF, W = f'{D}/relatives_ref', f'{D}/dap_weight'

if '--commit' in sys.argv:
    G = np.load(f'{REF}/geno.npy'); meta = json.load(gzip.open(f'{REF}/meta.json.gz', 'rt'))
    B = np.load(f'{REF}/geno_dap.npy'); dm = json.load(open(f'{REF}/meta_dap.json'))
    assert G.shape[0] == B.shape[0] and not any(s['source'] == 'dap' for s in meta['samples'])
    os.replace(f'{REF}/geno.npy', f'{REF}/geno_v3_backup.npy')
    os.replace(f'{REF}/meta.json.gz', f'{REF}/meta_v3_backup.json.gz')
    np.save(f'{REF}/geno.npy', np.hstack([G, B]))
    meta['samples'] += dm['samples']
    meta['built_v4'] = '2026-11-09'
    meta['note_v4'] = f"+{len(dm['samples'])} Dog Aging Project dogs (rounded imputed dosages; owner-reported adult weights; 47 likely-mixed samples with heterozygosity >= 0.38 excluded)"
    json.dump(meta, gzip.open(f'{REF}/meta.json.gz', 'wt'))
    print('v4:', G.shape[1], '+', B.shape[1], '=', G.shape[1] + B.shape[1], 'reference dogs')
    raise SystemExit

keys = {}
for i, line in enumerate(open(f'{REF}/keys.tsv')):
    c, p, ref, alt = line.rstrip('\n').split('\t'); keys[(c, int(p))] = (i, ref, alt)
panel = [l.rstrip('\n').split('\t') for l in open(f'{D}/breed_panel_v23/sites.tsv')][1:]
ids = [l.strip() for l in open(f'{W}/samples.txt') if l.strip()]
B = np.full((len(keys), len(ids)), -1, dtype=np.int8)
filled = flipped = 0
for ch in sorted({k[0] for k in keys}):
    f = f'{W}/parts/{ch}.npz'
    if not os.path.exists(f): continue
    z = np.load(f)
    pidx, pdos = z['panel_idx'], z['panel_dos']       # NpzFile re-reads on every [] access
    for r, pi in enumerate(pidx):
        c, pos, a1, a2 = panel[pi]
        k = keys.get((c, int(pos)))
        if k is None: continue
        i, ref, alt = k
        d = pdos[r].astype(np.float32)
        miss = d == 255
        d = d / 100.0
        if (a1, a2) == (alt, ref): x = d
        elif (a1, a2) == (ref, alt): x = 2.0 - d; flipped += 1
        else: continue
        g = np.where(x < 0.5, 0, np.where(x < 1.5, 1, 2)).astype(np.int8)
        g[miss] = -1
        B[i] = g; filled += 1
print(f'sites filled {filled}/{len(keys)} ({flipped} a1=REF)', 'missing frac', round(float((B < 0).mean()), 4))

# Labels and weights: owner-reported breed for purebreds, else our ancestry
# estimate ("<top breed> mix" below 85%); adult (>= 1 y) owner-reported weight.
z = np.load(f'{W}/build/dap.npz')
assert list(z['ids']) == ids
pretty = lambda code: code.replace('_', ' ').title()
ph = {}
for line in list(open(f'{W}/build/dap_pheno.tsv'))[1:]:
    did, wkg, age, male, neut, status, breed = line.rstrip('\n').split('\t')
    ph[did] = dict(w=float(wkg) if wkg else None, age=float(age) if age else None, status=status, breed=breed)
samples = []
for j, did in enumerate(ids):
    q = z['Q'][j]; t = int(np.argmax(q))
    p = ph.get(did, {})
    if p.get('status') == 'Purebred' and p.get('breed'):
        label = p['breed']
    else:
        label = pretty(str(z['breeds'][t])) + ('' if q[t] >= 0.85 else ' mix')
    s = {'id': f'DAP-{did}', 'label': label, 'source': 'dap'}
    if p.get('w') and (p.get('age') or 0) >= 1: s['weight_kg'] = round(p['w'], 1)
    samples.append(s)
# Drop likely mixed (contaminated) samples: heterozygosity far above the DAP
# range (median 0.30, 99th percentile 0.37). In validation a 0.55-het Labrador
# came out as a "relative" of 148 of 150 customers; removing the 47 dogs at
# >= 0.38 brought customers with any match from 150/150 to 8/150, the same
# rate as the existing customer reference (2026-11-09).
MAX_HET = 0.38
het = (B == 1).sum(0) / np.maximum((B >= 0).sum(0), 1)
keep = het < MAX_HET
B = B[:, keep]; samples = [s for s, k in zip(samples, keep) if k]
print(f'dropped {int((~keep).sum())} DAP dogs with heterozygosity >= {MAX_HET}')
np.save(f'{REF}/geno_dap.npy', B)
json.dump({'samples': samples}, open(f'{REF}/meta_dap.json', 'w'))
print('DAP dogs', len(samples), 'with weight', sum('weight_kg' in s for s in samples))
