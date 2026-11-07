#!/usr/bin/env python3
"""Merge per-kit disease-indel imputations (cluster/indel_window_cohort.sh: WORK/indel_imp/<NAME>/<kit>.tsv, one
CHROM POS REF ALT GT:GP line) into WORK/indel_imputed.tsv.gz, the cohort matrix panel_v3_breed_prior.py fit reads
(same layout as cohort_imputed.tsv.gz: CHROM POS REF ALT then one GT:GP column per kit, '.' where a kit has none).
Only windows listed in reference_panel/indel_windows/windows.tsv are merged.
  python3 analysis/omia/indel_window_merge.py WORK_DIR
"""
import gzip, os, pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
work = pathlib.Path(sys.argv[1])
kits = [l.split('\t')[0] for l in open(work / 'manifest.tsv') if l.strip()]
rows = []
for line in open(ROOT / 'reference_panel/indel_windows/windows.tsv'):
    name, c, p, ref, alt = line.rstrip('\n').split('\t')
    d = work / 'indel_imp' / name
    if not d.is_dir():
        continue
    vals, n = [], 0
    for k in kits:
        v = '.'
        f = d / f'{k}.tsv'
        if f.exists():
            for l in open(f):
                x = l.rstrip('\n').split('\t')
                if x[:4] == [c, p, ref, alt] and ':' in x[4]:
                    v = x[4]; n += 1
        vals.append(v)
    rows.append([c, p, ref, alt] + vals)
    print(f'{name}: {n} of {len(kits)} kits imputed')
with gzip.open(work / 'indel_imputed.tsv.gz', 'wt') as o:
    o.write('\t'.join(['CHROM', 'POS', 'REF', 'ALT'] + kits) + '\n')
    for r in rows:
        o.write('\t'.join(r) + '\n')
print(f'{len(rows)} indel windows -> {work}/indel_imputed.tsv.gz')
