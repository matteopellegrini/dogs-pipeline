#!/usr/bin/env python3
"""Compare report-level calls (breed, coat, OMIA conditions, traits) of each Cosmo run with COSMO-SUPER (~18x).
usage: compare_calls.py NAME=RESULTS_DIR ...   (first argument is the reference)"""
import json, os, sys

runs = [a.split('=', 1) for a in sys.argv[1:]]
def load(d, f):
    p = os.path.join(d, f)
    return json.load(open(p)) if os.path.exists(p) else None

def breed(d):
    b = load(d, 'breed_result.json')
    return [(x['breed_name'], round(100 * x['proportion'])) for x in b['breed_composition'] if x['proportion'] >= 0.05] if b else None

def coat(d):
    c = load(d, 'coat_color.json')
    return {k: '/'.join(v.get('predicted_alleles') or ['?']) for k, v in c['loci'].items()} if c else {}

def omia(d):
    o = load(d, 'omia_result.json')
    out = {}
    for v in (o or {}).get('variants', []):
        key = f"{v['gene']} {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}"
        call = next((x for k, x in v.items() if isinstance(x, dict) and 'zygosity' in x), None)
        out[key] = (call or {}).get('zygosity', 'none'), v.get('trait', ''), v.get('category', '')
    return out

def traits(d):
    t = load(d, 'traits_result.json')
    return {x['id']: ';'.join(str(l.get('n_alt')) for l in x.get('loci', [])) for x in (t or {}).get('traits', [])}

ref_name, ref_dir = runs[0]
R = dict(breed=breed(ref_dir), coat=coat(ref_dir), omia=omia(ref_dir), traits=traits(ref_dir))
print(f'REFERENCE {ref_name}: breed {R["breed"]}')
print('  coat ' + ' '.join(f'{k}={v}' for k, v in R['coat'].items()))
print('  OMIA non-ref: ' + '; '.join(f'{k.split()[0]} {z} ({t[:40]})' for k, (z, t, c) in R['omia'].items() if z not in ('ref', 'none', 'hom_ref', 'not_detected', 'no_call')))
for name, d in runs[1:]:
    b, c, o, t = breed(d), coat(d), omia(d), traits(d)
    cdiff = [f'{k} {R["coat"].get(k)}->{v}' for k, v in c.items() if R['coat'].get(k) != v]
    shared = [k for k in R['omia'] if k in o]
    odiff = [f'{k.split()[0]} {R["omia"][k][0]}->{o[k][0]} ({R["omia"][k][1][:30]})' for k in shared if R['omia'][k][0] != o[k][0]]
    tdiff = [f'{k} {R["traits"].get(k)}->{v}' for k, v in t.items() if R['traits'].get(k) != v]
    print(f'\n== {name}: breed {b}')
    print(f'  coat diffs ({len(cdiff)}): ' + '; '.join(cdiff))
    print(f'  OMIA diffs ({len(odiff)}/{len(shared)}): ' + '; '.join(odiff))
    print(f'  trait diffs ({len(tdiff)}/{len(t)}): ' + '; '.join(tdiff))
