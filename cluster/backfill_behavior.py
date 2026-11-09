#!/usr/bin/env python3
"""
Add / refresh dap_behavior in published prs_result.json files by running stage 11's own
DAP behaviour block (between the '# ── DAP behaviour (begin)' and '(end)' markers in
run_dog_pipeline.sh) on each kit's stored breed_result.json and coverage_1mb.json.
Backs up prs_result.json as prs_result.pre-<TAG>.json (once) before rewriting.

    backfill_behavior.py run_dog_pipeline.sh REF_JSON_DIR TAG LIST [--dry-run] > changes.tsv
LIST: pk, barcode, pub dir, work dir (tab-separated), as for cluster/republish_files.sh.
"""
import io, json, os, shutil, sys, contextlib
import numpy as np

pipe, ref_json, tag, lst = sys.argv[1:5]
dry = '--dry-run' in sys.argv
src = open(pipe).read()
start = src.index('# ── DAP behaviour (begin)')
end = src.index('# ── DAP behaviour (end)', start)
code = compile(src[start:end].replace('$REF_JSON', ref_json), 'stage11_dap_behavior', 'exec')

traits = None
n = {'updated': 0, 'skipped': 0}
for line in open(lst):
    f = line.rstrip('\n').split('\t')
    if len(f) < 3: continue
    pk, pub = f[0], f[2]
    path = f'{pub}/prs_result.json'
    try:
        d = json.load(open(path))
    except Exception as e:
        print(f'{pk}\tskip:{type(e).__name__}'); n['skipped'] += 1; continue
    g = {'json': json, 'np': np, 'PUB': pub}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(code, g)
    b = g.get('dap_behavior')
    if not b:
        print(f'{pk}\tskip:no_behavior'); n['skipped'] += 1; continue
    if traits is None:
        traits = list(b['traits'])
        print('pk\tsex\t' + '\t'.join(traits))
    print(f'{pk}\t{b["sex_used"]}\t' + '\t'.join(f'{b["traits"][t]["percentile"]:.0f}' for t in traits))
    n['updated'] += 1
    if dry: continue
    bak = f'{pub}/prs_result.pre-{tag}.json'
    if not os.path.exists(bak):
        shutil.copy2(path, bak)
    d['dap_behavior'] = b
    tmp = path + '.tmp'
    json.dump(d, open(tmp, 'w'), indent=2)
    os.replace(tmp, path)
print(f'# {n}', file=sys.stderr)
