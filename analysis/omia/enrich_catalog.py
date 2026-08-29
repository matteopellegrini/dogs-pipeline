#!/usr/bin/env python3
"""Enrich reference_json/omia_variants.json: dedupe, inherit traits, describe.

1. DEDUPE. Merge rows that describe the same variant: identical
   (chrom,pos,ref,alt), or same gene + same normalized trait with positions
   within 50bp (coordinate drift from different sources — the SLC13A1 pair).
   The surviving row keeps the richest fields; 'panel' strings are unioned.
2. TRAIT INHERITANCE. A row with no trait whose gene has exactly one distinct
   named trait elsewhere in the catalogue inherits that trait.
3. DESCRIPTIONS. Every row with a trait gets a 1-2 sentence owner-facing
   'description' from descriptions.py (condition-level, family rules).

Usage: python3 analysis/omia/enrich_catalog.py
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from descriptions import describe

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DB = os.path.join(ROOT, 'reference_json', 'omia_variants.json')

def tnorm(t):
    return " ".join(str(t or "").split()).strip().rstrip(",")

def richness(v):
    return sum(1 for f in ('trait','ref','alt','hgvs_c','hgvs_p','variant_breed','deleterious') if v.get(f))

def merge(keep, drop):
    for f, val in drop.items():
        if f == 'panel':
            panels = {p.strip() for p in (str(keep.get('panel') or '') + ',' + str(val or '')).split(',') if p.strip() and p.strip() != 'None'}
            keep['panel'] = ', '.join(sorted(panels))
        elif not keep.get(f) and val:
            keep[f] = val
    return keep

def main():
    db = json.load(open(DB))
    vs = db['variants']

    # pass 1: exact-position dedupe
    seen = {}
    out = []
    dropped = 0
    for v in vs:
        k = (v.get('chrom'), v.get('pos'), str(v.get('ref')), str(v.get('alt')))
        if k in seen and v.get('pos'):
            merge(seen[k], v); dropped += 1
            print('MERGED exact:', v.get('gene'), k[0], k[1])
        else:
            seen[k] = v; out.append(v)
    vs = out

    # pass 2: same gene+trait within 50bp
    out = []
    for v in vs:
        hit = None
        for u in out:
            if (u.get('gene') == v.get('gene') and tnorm(u.get('trait')) and
                    tnorm(u.get('trait')) == tnorm(v.get('trait')) and
                    u.get('chrom') == v.get('chrom') and u.get('pos') and v.get('pos') and
                    abs(int(u['pos']) - int(v['pos'])) <= 50 and int(u['pos']) != int(v['pos'])):
                hit = u; break
        if hit:
            keep, drop = (hit, v) if richness(hit) >= richness(v) else (v, hit)
            if keep is v:
                out.remove(hit); out.append(v)
            merge(keep, drop); dropped += 1
            print('MERGED near:', v.get('gene'), tnorm(v.get('trait'))[:40], hit.get('pos'), '/', v.get('pos'))
        else:
            out.append(v)
    vs = out

    # pass 3: trait inheritance by gene
    from collections import defaultdict
    g2t = defaultdict(set)
    for v in vs:
        if tnorm(v.get('trait')):
            g2t[v.get('gene')].add(tnorm(v.get('trait')))
    inherited = 0
    for v in vs:
        if not tnorm(v.get('trait')) and len(g2t.get(v.get('gene'), set())) == 1:
            v['trait'] = next(iter(g2t[v['gene']]))
            v['trait_inherited_from_gene'] = True
            inherited += 1

    # pass 4: descriptions
    described = 0
    missing = set()
    for v in vs:
        d = describe(v.get('trait'))
        if d:
            v['description'] = d; described += 1
        elif tnorm(v.get('trait')):
            missing.add(tnorm(v.get('trait')))

    db['variants'] = vs
    json.dump(db, open(DB, 'w'), indent=1)
    print(f'\nrows: {len(vs)} (merged away {dropped}) · traits inherited: {inherited} · described: {described}')
    if missing:
        print('traits still without description:', len(missing))
        for t in sorted(missing)[:20]: print('  -', t[:70])

if __name__ == '__main__':
    sys.exit(main())
