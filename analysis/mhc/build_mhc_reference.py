#!/usr/bin/env python3
"""
Build reference_json/mhc_reference.json for stage 9e (MHC / DLA diversity).

For each window (classII, classI) it stores:
  window, core          coordinates (canFam4 chr12)
  n_sites               Dog10K panel sites in the window
  cohort_het_sorted     sorted het_rate of the 2,053 processed dogs -> percentiles
  cohort_homozygous_frac fraction of dogs whose two core haplotypes differ at <= HOM_THR
  het_vs_froh           (intercept, slope) of het_rate ~ F_ROH in the cohort, so a
                        dog's MHC diversity can be compared with what its genome-wide
                        inbreeding predicts
  breeds                per panel breed (Dog10K label): n, mean het, homozygous frac
Inputs: mhc/{WIN}_sites.tsv, mhc/panel_{WIN}_gt.tsv, mhc/panel_samples.txt,
        mhc/cohort_mhc_{WIN}.tsv (from mhc_diversity.py), relatives_ref/meta.json.gz
"""
import gzip, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

D = os.environ.get('D', '/u/project/pellegrini/gkislik/dogs')
HOM_THR = 0.02
WINDOWS = {'classII': ((2300000, 2660000), (2330000, 2600000)),
           'classI':  ((1000000, 1200000), (1030000, 1175000))}


def gt_to_haps(gt_fields):
    n = len(gt_fields); h1 = np.full(n, -1, np.int8); h2 = np.full(n, -1, np.int8)
    for i, g in enumerate(gt_fields):
        if len(g) >= 3 and g[0] in '01' and g[2] in '01':
            h1[i] = int(g[0]); h2[i] = int(g[2])
    return h1, h2


def main():
    out = {'built': '2026-09-13', 'hom_thr': HOM_THR, 'chrom': 'chr12', 'windows': {}}
    meta = json.load(gzip.open(f'{D}/relatives_ref/meta.json.gz', 'rt'))
    label = {s['id']: s['label'] for s in meta['samples']}
    ids = [l.strip() for l in open(f'{D}/mhc/panel_samples.txt')]
    for win, (window, core) in WINDOWS.items():
        pos = np.array([int(l.split('\t')[1]) for l in open(f'{D}/mhc/{win}_sites.tsv')])
        rows = [l.rstrip('\n').split('\t') for l in open(f'{D}/mhc/panel_{win}_gt.tsv')]
        gts = np.array([r[4:] for r in rows]); assert gts.shape == (len(pos), len(ids))
        cm = (pos >= core[0]) & (pos <= core[1])
        breeds = {}
        for j, sid in enumerate(ids):
            h1, h2 = gt_to_haps(gts[:, j]); ok = (h1 >= 0) & (h2 >= 0)
            het = float(((h1 != h2) & ok).sum() / max(1, ok.sum()))
            dist = float((h1[cm] != h2[cm]).mean())
            b = breeds.setdefault(label.get(sid, '?'), {'n': 0, 'het': [], 'hom': 0})
            b['n'] += 1; b['het'].append(het); b['hom'] += int(dist <= HOM_THR)
        breed_out = {k: {'n': v['n'], 'mean_het': round(float(np.mean(v['het'])), 4), 'homozygous_frac': round(v['hom'] / v['n'], 3)}
                     for k, v in breeds.items() if v['n'] >= 3}
        crows = [l.rstrip('\n').split('\t') for l in open(f'{D}/mhc/cohort_mhc_{win}.tsv')][1:]
        het = np.array([float(r[1]) for r in crows]); dist = np.array([float(r[3]) for r in crows])
        froh = np.array([float(r[8]) if r[8] else np.nan for r in crows]); ok = ~np.isnan(froh)
        slope, intercept = np.polyfit(froh[ok], het[ok], 1)
        # Breed context from OUR cohort too: the Dog10K subset in the panel is
        # skewed to rare breeds (Poodle 1, German Shepherd 2, no Husky/Rottweiler),
        # while the cohort has dozens of each common breed. Cohort dogs whose top
        # breed is >= 50% count for that breed (n >= 5); cohort entries take
        # precedence over panel entries of the same name.
        cb = {}
        for r in crows:
            if r[9] and r[10] and float(r[10]) >= 0.5:
                e = cb.setdefault(r[9], {'n': 0, 'het': [], 'hom': 0}); e['n'] += 1; e['het'].append(float(r[1])); e['hom'] += int(float(r[3]) <= HOM_THR)
        for k, v in cb.items():
            if v['n'] >= 5:
                breed_out[k] = {'n': v['n'], 'mean_het': round(float(np.mean(v['het'])), 4), 'homozygous_frac': round(v['hom'] / v['n'], 3), 'source': 'prosperk9 cohort (top breed >= 50%)'}
        for k, v in breed_out.items(): v.setdefault('source', 'Dog10K panel')
        out['windows'][win] = {
            'window': list(window), 'core': list(core), 'n_sites': int(len(pos)),
            'cohort_n': int(len(het)), 'cohort_het_sorted': [round(float(x), 4) for x in np.sort(het)],
            'cohort_mean_het': round(float(het.mean()), 4), 'cohort_homozygous_frac': round(float((dist <= HOM_THR).mean()), 3),
            'het_vs_froh': {'intercept': round(float(intercept), 4), 'slope': round(float(slope), 4),
                            'resid_sd': round(float(np.std(het[ok] - (intercept + slope * froh[ok]))), 4)},
            'panel_mean_het': round(float(np.mean([np.mean(v['het']) for v in breeds.values()])), 4),
            'breeds': breed_out,
        }
        print(f'{win}: sites {len(pos)}, cohort {len(het)} (mean het {het.mean():.3f}, homozygous {(dist<=HOM_THR).mean()*100:.0f}%), '
              f'{len(breed_out)} breeds with n>=3, het ~ {intercept:.3f} + {slope:.3f}*F_ROH')
    json.dump(out, open(f'{D}/reference_json/mhc_reference.json', 'w'))
    print('wrote reference_json/mhc_reference.json')


if __name__ == '__main__':
    main()
