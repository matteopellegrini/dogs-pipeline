#!/usr/bin/env python3
"""Build the 50kb-window coverage panel-of-normals from the full cohort.

Inputs: work*/<sample>/analysis/coverage_50kb.tsv.gz (identical fixed grid,
47,091 windows, chr1-38 + chrX; samtools bedcov base counts).

Per sample: divide by the sample's autosomal median -> relative coverage.
Sex from chrX/autosome ratio (females ~1.0, males ~0.5).
Per window: median + MAD across panel dogs (chrX separately by sex), plus a
bad-window mask (median < 0.25 = effectively unmappable at 50kb;
MAD/median > 0.5 = unstable).

Excluded from the panel: dogs with aneuploidy-scale anomalies (any autosome
whose chromosome-median deviates >20% from 1.0 — catches trisomies and big
mosaic events; those dogs are validation cases, not normals) and dogs with
grossly unstable profiles (>5% of good windows beyond +-50%).

Outputs (Hoffman):
  reference_panel/coverage_50kb_panel.npz   (medians, MADs, mask, coords)
  reference_panel/coverage_50kb_panel_meta.json
  /u/project/pellegrini_archive/data/dogs_cnv/cov50k_matrix.npz  (full
  normalized matrix, float16, for segmentation development)
"""
import gzip, json, os, sys, glob
import numpy as np

D = '/u/project/pellegrini/gkislik/dogs'
os.chdir(D)

def load_grid(path):
    chroms, counts = [], []
    with gzip.open(path, 'rt') as f:
        for l in f:
            c = l.rstrip('\n').split('\t')
            chroms.append(c[0]); counts.append(float(c[3]))
    return chroms, np.array(counts, dtype=np.float64)

samples = []
for pat, root in (('work_prosper/pk-*/analysis/coverage_50kb.tsv.gz', 'prosper'),
                  ('work/DOGS-Gen-*/analysis/coverage_50kb.tsv.gz', 'gen')):
    for p in sorted(glob.glob(pat)):
        samples.append((p.split('/')[1], root, p))
print(f'input samples: {len(samples)}')

chroms = None
cols, names, groups = [], [], []
for s, root, p in samples:
    ch, v = load_grid(p)
    if chroms is None:
        chroms = np.array(ch)
        auto = chroms != 'chrX'
        xmask = ~auto
    if len(v) != len(chroms):
        print('SKIP bad grid:', s, len(v)); continue
    med = np.median(v[auto][v[auto] > 0])
    if not np.isfinite(med) or med <= 0:
        print('SKIP zero-coverage:', s); continue
    cols.append((v / med).astype(np.float32)); names.append(s); groups.append(root)
M = np.stack(cols, axis=1)   # windows x dogs
print('matrix:', M.shape)

# Sex from chrX median
xratio = np.median(M[xmask], axis=0)
sex = np.where(xratio > 0.75, 'F', 'M')
print('sex inferred: F =', int((sex == 'F').sum()), '| M =', int((sex == 'M').sum()))

# Aneuploidy / instability screen on autosomes
keep = np.ones(M.shape[1], dtype=bool)
chrom_ids = sorted(set(chroms[auto]))
excl = []
for j in range(M.shape[1]):
    prof = M[:, j]
    bad = False
    for c in chrom_ids:
        cm = np.median(prof[chroms == c])
        if abs(cm - 1.0) > 0.20:
            bad = True; excl.append((names[j], c, round(float(cm), 3))); break
    if not bad:
        frac_off = np.mean(np.abs(prof[auto] - 1.0) > 0.5)
        if frac_off > 0.05:
            bad = True; excl.append((names[j], 'unstable', round(float(frac_off), 3)))
    keep[j] = not bad
print(f'panel dogs kept: {int(keep.sum())} | excluded: {len(excl)}')
for e in excl[:15]: print('  excluded:', e)

K = M[:, keep]
ksex = sex[keep]
med_auto = np.median(K, axis=1)
mad = 1.4826 * np.median(np.abs(K - med_auto[:, None]), axis=1)
# chrX by sex
medF = np.median(K[:, ksex == 'F'], axis=1) if (ksex == 'F').any() else med_auto
medM_ = np.median(K[:, ksex == 'M'], axis=1) if (ksex == 'M').any() else med_auto
madF = 1.4826 * np.median(np.abs(K[:, ksex == 'F'] - medF[:, None]), axis=1)
madM = 1.4826 * np.median(np.abs(K[:, ksex == 'M'] - medM_[:, None]), axis=1)

bad_window = (med_auto < 0.25) | (mad / np.maximum(med_auto, 1e-6) > 0.5)
print(f'bad windows masked: {int(bad_window.sum())} of {len(bad_window)} ({100*bad_window.mean():.1f}%)')

# Read-length pooling check: 151bp vs 101bp per-window medians
rl = {}
for s in names:
    q = (f'results_prosper/{s}/qc_result.json' if s.startswith('pk-')
         else f'results/{s.lower().replace("dogs-gen","dogs-gen")}/qc_result.json')
    try: rl[s] = json.load(open(q)).get('read_length_bp')
    except Exception: rl[s] = None
kept_names = [n for n, k in zip(names, keep) if k]
i151 = [i for i, n in enumerate(kept_names) if rl.get(n) == 151]
i101 = [i for i, n in enumerate(kept_names) if rl.get(n) == 101]
if len(i151) > 30 and len(i101) > 30:
    m151 = np.median(K[:, i151], axis=1)[auto]
    m101 = np.median(K[:, i101], axis=1)[auto]
    ok = (m151 > 0.25) & (m101 > 0.25)
    r = np.corrcoef(np.log2(m151[ok]), np.log2(m101[ok]))[0, 1]
    frac2 = float(np.mean(np.abs(np.log2(m151[ok] / m101[ok])) > 0.32))
    print(f'read-length check (n151={len(i151)}, n101={len(i101)}): per-window median corr r={r:.4f}, {100*frac2:.2f}% windows differ >25%')

np.savez_compressed(f'{D}/reference_panel/coverage_50kb_panel.npz',
    chrom=chroms, med=med_auto.astype(np.float32), mad=mad.astype(np.float32),
    medF=medF.astype(np.float32), madF=madF.astype(np.float32),
    medM=medM_.astype(np.float32), madM=madM.astype(np.float32),
    bad=bad_window)
meta = {'n_input': len(names), 'n_panel': int(keep.sum()),
        'excluded': [list(e) for e in excl],
        'n_windows': int(len(chroms)), 'bad_windows': int(bad_window.sum()),
        'built': '2026-09-03', 'window_bp': 50000,
        'normalization': 'per-sample autosomal median; chrX stats by inferred sex',
        'read_lengths': {str(k): sum(1 for n in kept_names if rl.get(n) == k)
                         for k in sorted({v for v in rl.values() if v}, key=int)}}
json.dump(meta, open(f'{D}/reference_panel/coverage_50kb_panel_meta.json', 'w'), indent=1)
A = '/u/project/pellegrini_archive/data/dogs_cnv'
os.makedirs(A, exist_ok=True)
np.savez_compressed(f'{A}/cov50k_matrix.npz',
    M=M.astype(np.float16), names=np.array(names), sex=sex, chrom=chroms, keep=keep)
print('panel + matrix written')
