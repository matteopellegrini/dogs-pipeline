"""Hide owner-reported Dog Aging Project weights that cannot be right.

Relatives cards show a matched DAP dog's adult weight, and a few survey entries
are clearly wrong (a 2.7 kg "German Shepherd mix", a 52 kg "Chihuahua mix").
A weight is dropped when it is outside 1.2-110 kg, or far from the median
weight of DAP dogs whose owners report that breed (n >= 5): purebreds below
0.2x or above 2.5x, mixes below 0.2x or above 5x (and over 20 kg), since the
other side of a mix can be much bigger. Overweight dogs and size varieties
(standard Dachshunds, toy Poodles) stay. On the v4 reference this hides 6 of
6,775 weights (2026-11-09).

  python3 relatives/dap_weights.py   # patch relatives_ref/meta.json.gz + meta_dap.json in place"""
import collections, gzip, json, os, shutil, statistics


def implausible(samples):
    """Indices of DAP samples whose weight_kg should not be shown."""
    dap = [(i, s) for i, s in enumerate(samples) if s.get('source') == 'dap' and s.get('weight_kg')]
    base = lambda b: b[:-4] if b.endswith(' mix') else b
    by = collections.defaultdict(list)
    for _, s in dap:
        if not s['label'].endswith(' mix'): by[s['label']].append(s['weight_kg'])
    med = {b: statistics.median(v) for b, v in by.items() if len(v) >= 5}
    bad = set()
    for i, s in dap:
        x, md = s['weight_kg'], med.get(base(s['label']))
        if x < 1.2 or x > 110: bad.add(i)
        elif md is None: continue
        elif s['label'].endswith(' mix'):
            if x / md < 0.2 or (x / md > 5 and x > 20): bad.add(i)
        elif x / md < 0.2 or x / md > 2.5: bad.add(i)
    return bad


if __name__ == '__main__':
    REF = '/u/project/pellegrini/gkislik/dogs/relatives_ref'
    meta = json.load(gzip.open(f'{REF}/meta.json.gz', 'rt'))
    bad = implausible(meta['samples'])
    for i in bad:
        s = meta['samples'][i]
        print('hide', s['label'], s.pop('weight_kg'))
    shutil.copy(f'{REF}/meta.json.gz', f'{REF}/meta.pre-dapweight.json.gz')
    with gzip.open(f'{REF}/meta.json.gz.tmp', 'wt') as f: json.dump(meta, f)
    os.replace(f'{REF}/meta.json.gz.tmp', f'{REF}/meta.json.gz')  # running 9b jobs never see a partial file
    dm = json.load(open(f'{REF}/meta_dap.json'))
    for i in implausible(dm['samples']): dm['samples'][i].pop('weight_kg')
    json.dump(dm, open(f'{REF}/meta_dap.json', 'w'))
    print('hidden', len(bad))
