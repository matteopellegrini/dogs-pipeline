#!/usr/bin/env python3
"""Compare re-imputed staging reports (results_imp2/<kit>) against the published ones, kit by kit.

usage: compare_staging.py LIST OUT_DIR
LIST rows: name, barcode, pub dir, work dir (cluster/reimpute_list.tsv); kits without results_imp2/<kit>/.done are skipped.
Writes OUT_DIR/kits.tsv (one row per kit), OUT_DIR/omia_flips.tsv (every condition call that changed),
OUT_DIR/coat_traits.tsv (every coat locus / trait call that changed) and prints a cohort summary. A condition call counts
as changed only when what the app shows changes (affected, or carrier probability crossing 5%), not threshold jitter."""
import json, os, statistics, sys
from collections import Counter

D = '/u/project/pellegrini/gkislik/dogs'
lst, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)

def load(d, f):
    try:
        with open(os.path.join(d, f)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None

def breed(b):
    return {x['breed']: x['proportion'] for x in (b or {}).get('breed_composition', [])}

def l1(a, b):
    return sum(abs(a.get(k, 0) - b.get(k, 0)) for k in set(a) | set(b)) / 2

def omia_calls(o, name):
    calls = {}
    for v in (o or {}).get('variants', []):
        s = v.get(name) or next((x for k, x in v.items() if isinstance(x, dict) and 'zygosity' in x), None)
        if s:
            calls[v['site_id']] = (s.get('zygosity'), bool(s.get('affected')), s.get('call_confidence'), v.get('trait', ''), v.get('category', ''),
                                   float(s.get('p_carrier') or 0))
    return calls

def coat(c):
    return {k: '/'.join(map(str, l.get('predicted_alleles') or [])) + ' ' + str(l.get('confidence'))
            for k, l in ((c or {}).get('loci') or {}).items()}

def traits(t):
    return {x['id']: (x.get('status'), x.get('result')) for x in (t or {}).get('traits', [])}

def weight(p):
    w = ((p or {}).get('physical_traits') or {}).get('weight_kg') or {}
    return w.get('pred_kg')

def num(d, *ks):
    for k in ks:
        d = (d or {}).get(k) if isinstance(d, dict) else None
    return d

rows, flips, ct = [], [], []
for line in open(lst):
    name, bc, pub, _ = line.rstrip('\n').split('\t')[:4]
    st = f'{D}/results_imp2/{os.path.basename(pub)}'
    if not os.path.exists(f'{st}/.done'):
        continue
    r = {'kit': name}
    q = load(pub, 'qc_result.json') or {}
    r['depth'] = q.get('mean_depth') or q.get('depth') or ''
    bo, bn = breed(load(pub, 'breed_result.json')), breed(load(st, 'breed_result.json'))
    r['breed_top_old'] = max(bo, key=bo.get) if bo else ''
    r['breed_top_new'] = max(bn, key=bn.get) if bn else ''
    r['breed_l1'] = round(l1(bo, bn), 4)
    lo, ln = num(load(pub, 'local_ancestry.json'), 'proportions') or {}, num(load(st, 'local_ancestry.json'), 'proportions') or {}
    r['lai_l1'] = round(l1(lo, ln), 4)
    oo, on = omia_calls(load(pub, 'omia_result.json'), name), omia_calls(load(st, 'omia_result.json'), name)
    r['affected_old'] = sum(1 for c in oo.values() if c[1] and c[4] != 'coat')
    r['affected_new'] = sum(1 for c in on.values() if c[1] and c[4] != 'coat')
    nflip = 0
    for sid in sorted(set(oo) | set(on)):
        a, b = oo.get(sid), on.get(sid)
        vis = lambda c: (c[1], c[5] >= 0.05) if c else None   # what the app shows: affected, or an uncertain call at >= 5%
        if vis(a) != vis(b):
            nflip += 1
            flips.append((name, sid, (a or b)[3], (a or b)[4], a[0] if a else '-', a[2] if a else '-', b[0] if b else '-', b[2] if b else '-',
                          'gain' if b and b[1] and not (a and a[1]) else 'loss' if a and a[1] and not (b and b[1]) else 'zyg'))
    r['omia_changed'] = nflip
    co, cn = coat(load(pub, 'coat_color.json')), coat(load(st, 'coat_color.json'))
    r['coat_changed'] = sum(co.get(k) != cn.get(k) for k in set(co) | set(cn))
    ct += [(name, 'coat', k, co.get(k, '-'), cn.get(k, '-')) for k in sorted(set(co) | set(cn)) if co.get(k) != cn.get(k)]
    to, tn = traits(load(pub, 'traits_result.json')), traits(load(st, 'traits_result.json'))
    r['traits_changed'] = sum(to.get(k) != tn.get(k) for k in set(to) | set(tn))
    r['traits_called_old'] = sum(1 for v in to.values() if v[0] == 'called')
    r['traits_called_new'] = sum(1 for v in tn.values() if v[0] == 'called')
    ct += [(name, 'trait', k, ' | '.join(map(str, to.get(k, '-'))), ' | '.join(map(str, tn.get(k, '-')))) for k in sorted(set(to) | set(tn)) if to.get(k) != tn.get(k)]
    po, pn = load(pub, 'prs_result.json'), load(st, 'prs_result.json')
    r['weight_old'], r['weight_new'] = weight(po), weight(pn)
    r['froh_old'] = num(load(pub, 'inbreeding_result.json'), 'f_roh')
    r['froh_new'] = num(load(st, 'inbreeding_result.json'), 'f_roh')
    r['gF_old'] = num(load(pub, 'inbreeding_froh_dog10k_result.json'), 'sample_froh')
    r['gF_new'] = num(load(st, 'inbreeding_froh_dog10k_result.json'), 'sample_froh')
    ro, rn = load(pub, 'relatives_result.json') or {}, load(st, 'relatives_result.json') or {}
    mo = {m.get('id') or m.get('sample') or json.dumps(m, sort_keys=True)[:60] for m in ro.get('matches', [])}
    mn = {m.get('id') or m.get('sample') or json.dumps(m, sort_keys=True)[:60] for m in rn.get('matches', [])}
    r['relatives_old'], r['relatives_new'], r['relatives_same'] = len(mo), len(mn), len(mo & mn)
    r['mhc_old'] = json.dumps(num(load(pub, 'mhc_result.json'), 'summary'), sort_keys=True)
    r['mhc_new'] = json.dumps(num(load(st, 'mhc_result.json'), 'summary'), sort_keys=True)
    r['high_fx_old'] = len((load(pub, 'functional_variants.json') or {}).get('high_variants', []))
    r['high_fx_new'] = len((load(st, 'functional_variants.json') or {}).get('high_variants', []))
    rows.append(r)

cols = list(rows[0]) if rows else []
with open(f'{out}/kits.tsv', 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        f.write('\t'.join('' if r[c] is None else str(r[c]) for c in cols) + '\n')
with open(f'{out}/omia_flips.tsv', 'w') as f:
    f.write('kit\tsite_id\ttrait\tcategory\told_zyg\told_conf\tnew_zyg\tnew_conf\tkind\n')
    for x in flips:
        f.write('\t'.join(map(str, x)) + '\n')
with open(f'{out}/coat_traits.tsv', 'w') as f:
    f.write('kit\tkind\tid\told\tnew\n')
    for x in ct:
        f.write('\t'.join(map(str, x)) + '\n')

n = len(rows)
print(f'kits compared: {n}')
if not n:
    sys.exit()
med = lambda xs: statistics.median(xs) if xs else float('nan')
print(f"top breed changed: {sum(r['breed_top_old'] != r['breed_top_new'] for r in rows)}")
print(f"breed composition L1: median {med([r['breed_l1'] for r in rows]):.3f}, >0.10 in {sum(r['breed_l1'] > 0.10 for r in rows)}, max {max(r['breed_l1'] for r in rows):.3f}")
print(f"local ancestry L1: median {med([r['lai_l1'] for r in rows]):.3f}, >0.10 in {sum(r['lai_l1'] > 0.10 for r in rows)}")
k = Counter(x[8] for x in flips)
print(f"condition calls changed: {len(flips)} in {sum(r['omia_changed'] > 0 for r in rows)} kits (gain {k['gain']}, loss {k['loss']}, zygosity only {k['zyg']})")
print(f"affected (non-coat) per kit: old {sum(r['affected_old'] for r in rows)}, new {sum(r['affected_new'] for r in rows)}")
for (sid, trait, kind), c in Counter((x[1], x[2], x[8]) for x in flips if x[8] != 'zyg').most_common(15):
    print(f'   {kind:4s} {c:4d}  {sid} {trait[:60]}')
print(f"coat loci changed: {sum(r['coat_changed'] for r in rows)} in {sum(r['coat_changed'] > 0 for r in rows)} kits; by locus {dict(Counter(x[2] for x in ct if x[1] == 'coat'))}")
print(f"traits changed: {sum(r['traits_changed'] for r in rows)} in {sum(r['traits_changed'] > 0 for r in rows)} kits; called old {sum(r['traits_called_old'] for r in rows)} new {sum(r['traits_called_new'] for r in rows)}")
dw = [abs(r['weight_new'] - r['weight_old']) for r in rows if r['weight_old'] is not None and r['weight_new'] is not None]
print(f"PRS weight |change| kg: median {med(dw):.2f}, >3 kg in {sum(d > 3 for d in dw)}, max {max(dw) if dw else 0:.1f}")
df = [r['froh_new'] - r['froh_old'] for r in rows if r['froh_old'] is not None and r['froh_new'] is not None]
print(f"F_ROH change: median {med(df):+.4f}, mean {statistics.mean(df) if df else 0:+.4f}, |>0.03| in {sum(abs(d) > 0.03 for d in df)}")
dg = [r['gF_new'] - r['gF_old'] for r in rows if r['gF_old'] is not None and r['gF_new'] is not None]
print(f"genotype F change: median {med(dg):+.4f}")
print(f"relatives: old {sum(r['relatives_old'] for r in rows)}, new {sum(r['relatives_new'] for r in rows)}, kept {sum(r['relatives_same'] for r in rows)}")
print(f"MHC summary changed: {sum(r['mhc_old'] != r['mhc_new'] for r in rows)}")
print(f"high-effect variants per kit: old median {med([r['high_fx_old'] for r in rows])}, new {med([r['high_fx_new'] for r in rows])}")
