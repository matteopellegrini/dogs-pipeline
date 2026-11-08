#!/usr/bin/env python3
"""Rewrite only the CBD103 KB deletion entry (cat218) of published omia_result.json files from each kit's
coat_color.json K locus read call (panel_v3_apply.k_locus_genotype). Each touched report is backed up as
omia_result.pre-<TAG>.json (timestamp kept); only reports whose entry changes are rewritten, so
cluster/republish_files.sh <TAG> republishes exactly those.
  python3 kb_backfill.py TAG LIST [LIST ...]     (LIST rows: pk barcode pub work)
"""
import json, os, shutil, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panel_v3_apply import k_locus_genotype, is_gt, KB_SITE

def main(tag, *lists):
    n = collections.Counter()
    for lst in lists:
        for line in open(lst):
            f = line.rstrip('\n').split('\t')
            if len(f) < 3 or not f[0]:
                continue
            pub = f[2]
            try:
                o = json.load(open(f'{pub}/omia_result.json'))
                coat = json.load(open(f'{pub}/coat_color.json'))
            except (OSError, ValueError):
                n['missing'] += 1
                continue
            e = next((v for v in o['variants'] if v.get('site_id') == KB_SITE), None)
            key = e and next((k for k, x in e.items() if is_gt(x) and x.get('source') == 'panel_v3'), None)
            if not key:
                n['no_panel_v3_call'] += 1
                continue
            old = e[key]
            new = k_locus_genotype(coat, old)
            if new is None or new == old:
                n['unchanged'] += 1
                continue
            bak = f'{pub}/omia_result.pre-{tag}.json'
            if not os.path.exists(bak):
                shutil.copy2(f'{pub}/omia_result.json', bak)
            e[key] = new
            json.dump(o, open(f'{pub}/omia_result.json', 'w'), indent=1)
            n[f"{old['zygosity']} -> {new['zygosity']}"] += 1
    for k, v in sorted(n.items(), key=lambda x: -x[1]):
        print(v, k)

if __name__ == '__main__':
    main(*sys.argv[1:])
