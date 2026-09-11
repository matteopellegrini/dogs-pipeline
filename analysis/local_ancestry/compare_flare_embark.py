#!/usr/bin/env python3
# Compare a FLARE painting of cosmo3 with Embark's (embark_cosmo_segments.json,
# extracted exactly from Embark's PDF vector painting).
#
#   python3 compare_flare_embark.py <flare_out_prefix> [<flare_out_prefix> ...]
#
# FLARE .anc.vcf.gz carries AN1/AN2 (most probable ancestry per haplotype)
# per site; ancestry indices map to names via the .model file's first line.
# Positions are compared as FRACTIONS of chromosome length (Embark's bars are
# scaled to chromosome length; assemblies differ), and haplotype phase is
# arbitrary between the two, so agreement is scored on the UNORDERED pair of
# ancestries at each evaluation point. Our 230 breeds are collapsed to
# Embark's five classes + OTHER.
import gzip, json, sys
from collections import defaultdict

EMB = {'STANDARD_POODLE': 'POODLE_STANDARD', 'MINIATURE_POODLE': 'POODLE_SMALL', 'TOY_POODLE': 'POODLE_SMALL',
       'LABRADOR_RETRIEVER': 'LABRADOR_RETRIEVER', 'COCKER_SPANIEL': 'COCKER_SPANIEL',
       'ENGLISH_COCKER_SPANIEL': 'ENGLISH_COCKER_SPANIEL'}
CLASSES = ['POODLE_STANDARD', 'POODLE_SMALL', 'LABRADOR_RETRIEVER', 'COCKER_SPANIEL', 'ENGLISH_COCKER_SPANIEL', 'OTHER']

emb = json.load(open('embark_cosmo_segments.json'))
emb_by = defaultdict(list)
for s in emb:
    emb_by[(s['chrom'], s['hap'])].append((s['start_frac'], s['end_frac'], s['breed']))

def emb_pair(chrom, frac):
    pair = []
    for h in (1, 2):
        b = next((br for a, e, br in emb_by[(chrom, h)] if a <= frac < e), None)
        pair.append(b or 'OTHER')
    return tuple(sorted(pair))

def load_flare(prefix):
    names = open(prefix + '.model').readline().split()
    lens = {}
    try:   # canFam4 chromosome lengths (from the reference .fai), if staged
        for l in open('chrom_lengths.tsv'):
            c, n = l.split()[:2]; lens[c] = int(n)
    except FileNotFoundError:
        pass
    calls = defaultdict(list)   # chrom -> [(pos, a1, a2)]
    with gzip.open(prefix + '.anc.vcf.gz', 'rt') as f:
        for l in f:
            if l.startswith('##ANCESTRY'):
                names = {}
                for a in l.strip().split('=', 1)[1].strip('<>').split(','):
                    k, v = a.strip().split('='); names[int(v)] = k
                continue
            if l.startswith('##contig'):
                cid = l.split('ID=')[1].split(',')[0].split('>')[0]
                if 'length=' in l: lens[cid] = int(l.split('length=')[1].split('>')[0].split(',')[0])
                continue
            if l.startswith('#'):
                fmt_i = None; continue
            p = l.rstrip('\n').split('\t')
            fmt = p[8].split(':'); i1 = fmt.index('AN1'); i2 = fmt.index('AN2')
            g = p[9].split(':')
            calls[p[0]].append((int(p[1]), names[int(g[i1])], names[int(g[i2])]))
    return names, lens, calls

def summarize(prefix):
    names, lens, calls = load_flare(prefix)
    # global proportions from FLARE's own global file
    glob = {}
    with gzip.open(prefix + '.global.anc.gz', 'rt') as f:
        hdr = f.readline().split(); vals = f.readline().split()
        glob = dict(zip(hdr[1:], map(float, vals[1:])))
    coll = defaultdict(float)
    for b, v in glob.items(): coll[EMB.get(b, 'OTHER')] += v
    print(f'\n=== {prefix} ===')
    print('FLARE global (collapsed to Embark classes):', {k: round(100*coll[k], 1) for k in CLASSES})
    top = sorted(glob.items(), key=lambda kv: -kv[1])[:8]
    print('FLARE global top breeds:', [(b, round(100*v, 1)) for b, v in top])
    # position-level agreement on a grid of evaluation points per chromosome
    agree = both = exact = 0; n = 0
    per_chrom = {}
    for c in range(1, 39):
        cid = f'chr{c}'
        pts = calls.get(cid, [])
        if not pts: continue
        L = lens.get(cid) or max(p for p, _, _ in pts)
        ca = 0; cn = 0
        for k in range(50):
            frac = (k + 0.5) / 50
            target = frac * L
            # nearest FLARE site
            pos, a1, a2 = min(pts, key=lambda t: abs(t[0] - target))
            fp = tuple(sorted((EMB.get(a1, 'OTHER'), EMB.get(a2, 'OTHER'))))
            ep = emb_pair(c, frac)
            n += 1; cn += 1
            if fp == ep: exact += 1; ca += 1
            if fp[0] in ep or fp[1] in ep: agree += 1
        per_chrom[c] = ca / cn if cn else None
    print(f'position-level agreement vs Embark at {n} points: exact unordered pair {100*exact/n:.1f}% | at least one haplotype {100*agree/n:.1f}%')
    worst = sorted(per_chrom.items(), key=lambda kv: kv[1])[:5]
    print('lowest-agreement chromosomes:', [(c, round(100*v)) for c, v in worst])
    # segment count (switch density) per haplotype as a fragmentation measure
    nseg = 0
    for cid, pts in calls.items():
        for hap in (1, 2):
            prev = None
            for pos, a1, a2 in pts:
                a = a1 if hap == 1 else a2
                if a != prev: nseg += 1; prev = a
    print(f'FLARE segments (both haplotypes, 38 chr): {nseg}   [Embark painting: 526]')

for pre in sys.argv[1:]:
    summarize(pre)
