#!/usr/bin/env python3
# Score a lai_hmm.py painting (<prefix>.lai.json) against Embark's painting of
# Cosmo (embark_cosmo_segments.json): positions as chromosome fractions,
# unordered haplotype pairs, at 5-class (Embark's labels) and family level.
#   python3 compare_hmm_embark.py <embark_segments.json> <prefix.lai.json> [...]
import json, sys
from collections import defaultdict
emb_segs = json.load(open(sys.argv[1]))
emb = defaultdict(lambda: defaultdict(list))
for s in emb_segs: emb[s['chrom']][s['hap']].append((s['start_frac'], s['end_frac'], s['breed']))
def cls5(b):
    b = b.upper()
    if b in ('STANDARD_POODLE', 'POODLE_STANDARD'): return 'POODLE_STANDARD'
    if 'POODLE' in b: return 'POODLE_SMALL'
    if b in ('COCKER_SPANIEL',): return 'COCKER_SPANIEL'
    if 'ENGLISH_COCKER' in b: return 'ENGLISH_COCKER_SPANIEL'
    if 'LABRADOR' in b: return 'LABRADOR_RETRIEVER'
    return 'OTHER'
def fam(b):
    b = b.upper(); return 'POODLE' if 'POODLE' in b else 'COCKER' if 'COCKER' in b else 'LAB' if 'LABRADOR' in b else 'OTHER'
def emb_at(c, frac):
    out = []
    for hap in (1, 2):
        b = None
        for f0, f1, br in emb[c][hap]:
            if f0 <= frac <= f1: b = br; break
        out.append(b)
    return out
for path in sys.argv[2:]:
    res = json.load(open(path)); segs = res['segments']
    bychr = defaultdict(list)
    for s in segs: bychr[int(s['chrom'].replace('chr', ''))].append(s)
    n = f5 = ff = 0; low = {}
    for c, ss in bychr.items():
        if c not in emb: continue
        L = max(s['end'] for s in ss); okc = nc = 0
        for s in ss:
            for frac in (x / 20.0 for x in range(21)):           # 21 points per segment span, dense enough
                pos = s['start'] + frac * (s['end'] - s['start'])
                e = emb_at(c, pos / L)
                if None in e: continue
                a = s['anc']; n += 1; nc += 1
                ok5 = sorted(map(cls5, e)) == sorted(map(cls5, a)); okf = sorted(map(fam, e)) == sorted(map(fam, a))
                f5 += ok5; ff += okf; okc += okf
        low[c] = round(100 * okc / nc) if nc else None
    worst = sorted((v, k) for k, v in low.items() if v is not None)[:4]
    print(f'{path}: gen={res["gen"]} temper={res.get("temper")} segments={res["n_segments"]} | vs Embark: 5-class {100*f5/n:.1f}%  family {100*ff/n:.1f}%  (n={n}); lowest chromosomes {worst}')
    print('   proportions:', {k: round(v, 3) for k, v in sorted(res['proportions'].items(), key=lambda kv: -kv[1])})
