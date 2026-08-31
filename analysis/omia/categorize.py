#!/usr/bin/env python3
"""Assign each catalogue variant a display category: coat | trait | disease.

Drives the dashboard's tab split: coat variants render on the Coat Color tab,
single-gene trait variants on Trait Scores, and everything else (the default)
on Genetic Health Variants. Rule for dual-nature variants: if it can hurt the
dog, the health tab owns it (CDDY/IVDD is disease; CDPA short legs is trait).

Writes `category` into reference_json/omia_variants.json and emits
dogs-app/lib/variantCategories.json ({"chrom:pos": category}) so the frontend
can categorize already-published staged files without any pipeline rerun.

Usage: python3 analysis/omia/categorize.py
"""
import json, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DB = os.path.join(ROOT, 'reference_json', 'omia_variants.json')
OUT_MAP = os.path.join(ROOT, 'dogs-app', 'lib', 'variantCategories.json')

# Substring rules on the trait name, case-insensitive. First match wins;
# anything unmatched (including the 76 unnamed rows) defaults to disease.
COAT = [
    'brown color', 'cocoa color', 'color dilution', 'recessive red', 'ancient red',
    'grizzle', 'recessive black', 'hidden patterning', 'harlequin',
    'red pigment intensity', 'eye color', 'albinism',
]
TRAIT = [
    'coat length', 'curly coat', 'hair ridge', 'shedding', 'taillessness',
    'muscular hypertrophy',
]
# Exact overrides for dual-nature calls that the substring rules would misfile.
# (ITGA10 "Chondrodystrophy" is the Norwegian Elkhound dwarfism DISEASE and
# stays in the disease default; the benign short-leg CDPA is not catalogued.)
EXACT = {
    'chondrodystrophy and intervertebral disc disease (cddy/ivdd, type i ivdd)': 'disease',
}

def categorize(trait):
    t = ' '.join((trait or '').lower().split())
    if t in EXACT:
        return EXACT[t]
    for kw in COAT:
        if kw in t:
            return 'coat'
    for kw in TRAIT:
        if kw in t:
            return 'trait'
    return 'disease'

def main():
    db = json.load(open(DB))
    counts = {'coat': 0, 'trait': 0, 'disease': 0}
    cmap = {}
    for v in db['variants']:
        c = categorize(v.get('trait') or v.get('phene_name'))
        v['category'] = c
        counts[c] += 1
        if v.get('chrom') and v.get('pos'):
            cmap[f"{v['chrom']}:{v['pos']}"] = c
    json.dump(db, open(DB, 'w'), indent=1)
    json.dump(cmap, open(OUT_MAP, 'w'), indent=0, sort_keys=True)
    print('categories:', counts, '| map entries:', len(cmap))
    for c in ('coat', 'trait'):
        print(f'-- {c} --')
        for v in db['variants']:
            if v['category'] == c:
                print('  ', v.get('gene'), '|', (v.get('trait') or '')[:60])

if __name__ == '__main__':
    main()
