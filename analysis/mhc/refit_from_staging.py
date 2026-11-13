#!/usr/bin/env python3
"""Refit reference_json/mhc_reference.json on the panel v2 re-imputed cohort (results_imp2 staging).

The reference's cohort fields (cohort_het_sorted percentiles, homozygous share, het ~ F_ROH fit, cohort breed
context) were fitted on the old imputation and the old stage 12 F_ROH, which with panel v2 collapses to ~0.
This writes the builder's cohort tables from the staged reports instead of the September extraction:
het_rate / core_hap_distance from each kit's mhc_result.json (they depend only on its imputed BCF), f_roh from
its inbreeding_result.json (depth-corrected stage 12), top breed from its mhc_result.json. Then it runs the
unchanged build_mhc_reference.py on them in a scratch D, so panel and breed logic stay identical.

usage: refit_from_staging.py LIST OUT_JSON
LIST rows: name, barcode, pub dir, work dir; staging dir = results_imp2/<basename of pub dir>.
Downsampled replicates (name contains -DS) are left out of the cohort.
"""
import json, os, subprocess, sys, tempfile

D = '/u/project/pellegrini/gkislik/dogs'
lst, out_json = sys.argv[1], sys.argv[2]

rows = {'classII': [], 'classI': []}
skipped = 0
for line in open(lst):
    f = line.rstrip('\n').split('\t')
    if len(f) < 3 or '-DS' in f[0].upper():
        continue
    st = os.path.join(D, 'results_imp2', os.path.basename(f[2].rstrip('/')))
    if not os.path.exists(f'{st}/.done'):
        skipped += 1
        continue
    try:
        mhc = json.load(open(f'{st}/mhc_result.json'))
        ib = json.load(open(f'{st}/inbreeding_result.json'))
    except (OSError, ValueError):
        skipped += 1
        continue
    if 'bcftools roh' not in ib.get('method', ''):
        sys.exit(f'{st}: inbreeding_result.json is not from the depth-corrected stage 12; run cluster/inbreeding_pass.sh first')
    froh = ib.get('f_roh')
    for win in rows:
        w = (mhc.get('windows') or {}).get(win) or {}
        if w.get('status') != 'ok' or w.get('core_hap_distance') is None:
            continue
        # builder columns: 1 het_rate, 3 core_hap_distance, 8 f_roh, 9 top_breed, 10 top_prop
        rows[win].append([f[0], w['het_rate'], '', w['core_hap_distance'], '', '', '', '',
                          '' if froh is None else froh, mhc.get('top_breed') or '', mhc.get('top_breed_proportion') or ''])

with tempfile.TemporaryDirectory(dir=f'{D}/scratch') as tmp:
    os.makedirs(f'{tmp}/mhc'); os.makedirs(f'{tmp}/reference_json')
    os.symlink(f'{D}/relatives_ref', f'{tmp}/relatives_ref')
    for fn in os.listdir(f'{D}/mhc'):
        if fn.startswith('panel_') or fn.endswith('_sites.tsv'):
            os.symlink(f'{D}/mhc/{fn}', f'{tmp}/mhc/{fn}')
    for win, rs in rows.items():
        with open(f'{tmp}/mhc/cohort_mhc_{win}.tsv', 'w') as fh:
            fh.write('sample\thet_rate\thet_pct_cohort\tcore_hap_distance\tgroup1\tgroup2\thomozygous\tdepth\tf_roh\ttop_breed\ttop_prop\n')
            for r in rs:
                fh.write('\t'.join(map(str, r)) + '\n')
    subprocess.run([sys.executable, f'{D}/analysis/mhc/build_mhc_reference.py'], env={**os.environ, 'D': tmp}, check=True)
    ref = json.load(open(f'{tmp}/reference_json/mhc_reference.json'))

ref['built'] = '2026-11-13'
ref['cohort_source'] = 'panel v2 re-imputed cohort (results_imp2), depth-corrected stage 12 F_ROH'
json.dump(ref, open(out_json, 'w'))
print(f'cohort rows: classII {len(rows["classII"])}, classI {len(rows["classI"])}; kits without staged results {skipped}; wrote {out_json}')
