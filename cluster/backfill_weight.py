#!/usr/bin/env python3
"""
Recompute physical_traits.weight_kg in published prs_result.json files with the
current reference_json/weight_breed_model.json, by running stage 11's own weight
block (extracted from run_dog_pipeline.sh) on each kit's stored files: breed_result.json,
coverage_1mb.json (sex) and the stored da_size_prs / blend / dense values. No BCF needed.
Backs up prs_result.json as prs_result.pre-<TAG>.json (once) before rewriting.

    backfill_weight.py run_dog_pipeline.sh REF_JSON_DIR TAG LIST [--dry-run] > changes.tsv
LIST: pk, barcode, pub dir, work dir (tab-separated), as for cluster/republish_files.sh.
"""
import json, os, sys, shutil, io, contextlib
import numpy as np

pipe, ref_json, tag, lst = sys.argv[1:5]
dry = '--dry-run' in sys.argv
src = open(pipe).read()
start = src.index("    _wbm, _wcomp = None, []\n")
end = src.index("\n    else:\n        phys_traits['weight_kg'].update({\n            'pred_kg': round(_blend_kg", start)
block = '\n'.join(l[4:] if l.startswith('    ') else l for l in src[start:end].split('\n'))
block = block.replace('$REF_JSON', ref_json)
code = compile(block, 'stage11_weight_block', 'exec')

print('pk\told_kg\tnew_kg\tbreed_kg\tsex\tstatus')
n = {'updated': 0, 'unchanged': 0, 'skipped': 0}
for line in open(lst):
    f = line.rstrip('\n').split('\t')
    if len(f) < 3: continue
    pk, pub = f[0], f[2]
    path = f'{pub}/prs_result.json'
    try:
        d = json.load(open(path))
        w = d['physical_traits']['weight_kg']
    except Exception as e:
        print(f'{pk}\t\t\t\t\tskip:{type(e).__name__}'); n['skipped'] += 1; continue
    if 'da_size_prs' not in w:
        print(f'{pk}\t{w.get("pred_kg")}\t\t\t\tskip:no_da_prs'); n['skipped'] += 1; continue
    old = w.get('pred_kg')
    m, t = (w.get('da_sites_matched') or '0/0').split('/')
    g = {'json': json, 'np': np, 'PUB': pub, 'phys_traits': d['physical_traits'],
         '_da_prs': float(w['da_size_prs']), '_da_matched': int(m), '_da_rows': range(int(t)),
         '_blend_kg': float(w.get('pred_kg_blend', old)), '_dense_kg': float(w.get('pred_kg_dense', old))}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(code, g)
    if not (g.get('_wbm') and g.get('_wcomp')):
        print(f'{pk}\t{old}\t\t\t\tskip:no_model_or_comp'); n['skipped'] += 1; continue
    w = d['physical_traits']['weight_kg']
    status = 'updated' if w.get('method_note', '').startswith('Breed-composition adult-weight model: per-breed weights fit to') else 'unexpected'
    print(f'{pk}\t{old}\t{w["pred_kg"]}\t{w.get("pred_kg_breed")}\t{w.get("sex_used")}\t{status}')
    if status != 'updated':
        n['skipped'] += 1; continue
    n['updated'] += 1
    if dry: continue
    bak = f'{pub}/prs_result.pre-{tag}.json'
    if not os.path.exists(bak):
        shutil.copy2(path, bak)
    tmp = path + '.tmp'
    json.dump(d, open(tmp, 'w'), indent=2)
    os.replace(tmp, path)
print(f'# {n}', file=sys.stderr)
