#!/usr/bin/env python3
"""
Build a canonical breed label for every Parker and Dog10K population code.

    python3 harmonize.py PIPELINE.sh DOG10K_TABLE_S1.xlsx CLUST.txt OUT.json
        [--village per-country|region|pooled]

Two panels, two code vocabularies, one biology. Parker calls Belgian Tervuren
TURV, Dog10K calls it TERV; Parker CIRN_Italy is Dog10K CIRN; Dog10K uses both
FXTE and SMFX for Smooth Fox Terrier. Left alone these count as errors while
being pure bookkeeping, and worse, a customer could see the same breed listed
twice.

Names come from the two authoritative sources rather than from guesswork:
PARKER_NAMES inside the pipeline, and Table S1 of the Dog10K supplement
(Sample Name -> Breed/Type, plus a Category column separating breed dogs from
village dogs, wolves and coyotes).

Non-breed categories need a decision, not a default:

  wolves   pooled into one GRAY_WOLF. Dog10K has 68 across 12 countries; the
           regional populations are not separable from each other (CLUPCN dogs
           predict as CLUPRU), and "does this dog have wolf ancestry" is the
           question a customer actually asks.
  coyotes  kept separate as COYOTE. Only 4, so they fall below any sensible
           min-n and drop out — which is right, they are not dogs.
  village  --village selects the grouping, because this is exactly the case
           where pooling can backfire. Merging the regional Salukis made a
           broad attractor that swallowed Anatolian Shepherds, so the choice
           between per-country, continental region and one global pool is
           settled by held-out benchmark, not by assertion.
"""
import json
import re
import sys
import zipfile
from collections import Counter, defaultdict
from xml.etree import ElementTree as ET

NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

REGION = {
    'China': 'EastAsia', 'Vietnam': 'EastAsia', 'Thailand': 'EastAsia',
    'Philippines': 'EastAsia', 'Taiwan': 'EastAsia', 'Indonesia': 'EastAsia',
    'Papua New Guinea': 'Oceania', 'French Polynesia': 'Oceania', 'Fiji': 'Oceania',
    'Iran': 'WestAsia', 'Afghanistan': 'WestAsia', 'Qatar': 'WestAsia',
    'Lebanon': 'WestAsia', 'Turkey': 'WestAsia', 'Israel': 'WestAsia',
    'Uzbekistan': 'CentralAsia', 'Tajikistan': 'CentralAsia',
    'Kazakhstan': 'CentralAsia', 'Kyrgyzstan': 'CentralAsia', 'Mongolia': 'CentralAsia',
    'India': 'SouthAsia', 'Nepal': 'SouthAsia', 'Sri Lanka': 'SouthAsia',
    'Pakistan': 'SouthAsia', 'Bangladesh': 'SouthAsia',
    'Kenya': 'Africa', 'Liberia': 'Africa', 'Congo': 'Africa', 'Nigeria': 'Africa',
    'Namibia': 'Africa', 'Uganda': 'Africa', 'Egypt': 'Africa', 'Morocco': 'Africa',
    'Mali': 'Africa', 'Botswana': 'Africa',
    'Mexico': 'Americas', 'Peru': 'Americas', 'Brazil': 'Americas',
    'Colombia': 'Americas', 'Bahamas': 'Americas', 'Puerto Rico': 'Americas',
    'USA': 'Americas', 'Belize': 'Americas', 'Costa Rica': 'Americas',
}


def read_table_s1(path):
    z = zipfile.ZipFile(path)
    shared = [''.join(t.text or '' for t in si.iter(f'{NS}t'))
              for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall(f'{NS}si')]
    root = ET.fromstring(z.read('xl/worksheets/sheet2.xml'))
    out = []
    for row in root.iter(f'{NS}row'):
        cells = {}
        for c in row.findall(f'{NS}c'):
            col = re.match(r'[A-Z]+', c.get('r')).group()
            v = c.find(f'{NS}v'); t = c.get('t')
            cells[col] = ((shared[int(v.text)] if t == 's' else v.text)
                          if v is not None else '').strip()
        out.append(cells)
    return [r for r in out[3:] if r.get('A')]


def parker_names(pipeline_path):
    s = open(pipeline_path).read()
    m = re.search(r'PARKER_NAMES\s*=\s*\{(.*?)\n\}', s, re.S)
    # names may contain escaped apostrophes ("Cirneco dell\'Etna"), so the
    # value pattern has to allow backslash escapes or it truncates mid-name.
    pairs = re.findall(r"'((?:[^'\\]|\\.)*)'\s*:\s*'((?:[^'\\]|\\.)*)'", m.group(1))
    return {k.replace("\\'", "'"): v.replace("\\'", "'") for k, v in pairs}


def canon(name):
    """normalise a breed name so the two vocabularies line up"""
    n = re.sub(r'\([^)]*\)', ' ', name or '')
    n = re.sub(r'[^A-Za-z]+', ' ', n).upper().strip()
    n = re.sub(r'\s+', ' ', n)
    # the two sources differ only by these qualifiers on otherwise identical breeds
    n = re.sub(r'^CHINESE (SHAR PEI)$', r'\1', n)
    n = re.sub(r'^AMERICAN (COCKER SPANIEL)$', r'\1', n)
    n = re.sub(r'\bDOG$', '', n).strip()
    # PLINK cluster IDs are whitespace-delimited in .frq.strat, so a label with
    # a space silently shifts every field after it — for SCOPE as well as for
    # our own readers. Join with underscores.
    return n.replace(' ', '_')


def main():
    pipeline, xlsx, clust, out = sys.argv[1:5]
    args = sys.argv[5:]
    village = args[args.index('--village') + 1] if '--village' in args else 'region'

    pn = parker_names(pipeline)
    rows = read_table_s1(xlsx)

    # Dog10K: code -> (name, category)
    d_name, d_cat = {}, {}
    for r in rows:
        m = re.match(r'^([A-Za-z\-]+?)\d+$', r['A'])
        if not m:
            continue
        c = m.group(1)
        d_name.setdefault(c, Counter())[r.get('B', '')] += 1
        d_cat.setdefault(c, Counter())[r.get('C', '')] += 1
    d_name = {c: v.most_common(1)[0][0] for c, v in d_name.items()}
    d_cat = {c: v.most_common(1)[0][0] for c, v in d_cat.items()}

    parker_codes = set()
    for line in open(clust):
        p = line.split()
        if len(p) >= 3:
            parker_codes.add(p[2])

    label = {}
    for c in parker_codes:
        nm = pn.get(c, c)
        if c.upper().startswith('WOLF'):
            label[c] = 'GRAY_WOLF'
        else:
            label[c] = canon(nm) or c
    for c, nm in d_name.items():
        cat = d_cat.get(c, '')
        if cat == 'Wolf':
            label[c] = 'GRAY_WOLF'
        elif cat == 'Coyote':
            label[c] = 'COYOTE'
        elif cat == 'Village_Dogs':
            if village == 'pooled':
                label[c] = 'VILLAGE_DOG'
            elif village == 'region':
                label[c] = 'VILLAGE_' + REGION.get(nm, 'Other')
            else:
                label[c] = 'VILLAGE_' + re.sub(r'[^A-Za-z]', '', nm)
        elif cat == 'Mixed/Other':
            label[c] = 'MIXED_OTHER'
        else:
            label[c] = canon(nm) or c

    json.dump(label, open(out, 'w'), indent=0, sort_keys=True)

    groups = defaultdict(list)
    for c, l in label.items():
        groups[l].append(c)
    multi = {l: cs for l, cs in groups.items() if len(cs) > 1}
    print(f"codes mapped: {len(label)}   canonical labels: {len(groups)}")
    print(f"labels drawing on more than one code: {len(multi)}\n")
    for l, cs in sorted(multi.items())[:24]:
        src = ', '.join(sorted(cs))
        print(f"   {l:<34} <- {src}")
    if village != 'pooled':
        vg = {l: len(cs) for l, cs in groups.items() if l.startswith('VILLAGE')}
        print(f"\nvillage groupings ({village}): {len(vg)} -> {sorted(vg)[:8]}")
    print(f"\nwritten: {out}")


if __name__ == '__main__':
    main()
