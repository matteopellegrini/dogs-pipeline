# Extract Embark's "breed mix by chromosome" painting as numeric segments.
# Layout: 4-column grid of chromosome labels ("1".."38"); right of each label,
# two stacked bars (haplotypes) made of coloured rects. Segment coordinates are
# expressed as fractions of the chromosome bar (bars are scaled to chromosome
# length), so they can be compared to any assembly by relative position.
import pdfplumber, json
from collections import defaultdict
COL = {(1.0,0.8078,0.2039):'POODLE_STANDARD', (0.3412,0.7569,0.9137):'POODLE_SMALL',
       (0.6549,0.4627,0.651):'COCKER_SPANIEL', (0.9451,0.3569,0.3608):'LABRADOR_RETRIEVER',
       (0.9765,0.6235,0.1176):'ENGLISH_COCKER_SPANIEL'}
pdf = pdfplumber.open('/Users/matteopellegrini/Downloads/MDN78-XKDPP-3EWKZ-4R25P (1) (4).pdf')
p = pdf.pages[0]
rects = [r for r in p.rects if tuple(r.get('non_stroking_color') or ()) in COL]
# chromosome labels: numeric words below the legend (y > 400)
labels = [w for w in p.extract_words() if w['text'].isdigit() and 1 <= int(w['text']) <= 38 and w['top'] > 400]
labels = sorted(labels, key=lambda w: (round(w['top']), w['x0']))
# assign each rect to the nearest label to its left in the same grid row band (within ~22pt vertically)
segs = defaultdict(list)   # (chrom, hap) -> list of (x0,x1,breed)
for r in rects:
    cands = [w for w in labels if w['x0'] < r['x0'] and (w['top'] - 2 <= r['top'] <= w['top'] + 9)]
    if not cands: continue
    w = max(cands, key=lambda w: w['x0'])
    chrom = int(w['text'])
    segs[(chrom, r['top'])].append((r['x0'], r['x1'], COL[tuple(r['non_stroking_color'])]))
# per chromosome: two distinct 'top' values = haplotype 1 (upper) and 2 (lower)
out = []; hap_of = {}
for chrom in range(1, 39):
    tops = sorted({k[1] for k in segs if k[0] == chrom})
    if len(tops) != 2: print('chrom', chrom, 'has', len(tops), 'rows'); 
    for h, t in enumerate(tops[:2], start=1):
        ss = sorted(segs[(chrom, t)])
        x0 = min(s[0] for s in ss); x1 = max(s[1] for s in ss); L = x1 - x0
        # merge adjacent same-breed rects
        merged = []
        for a, b, br in ss:
            if merged and merged[-1][2] == br and a - merged[-1][1] < 0.6: merged[-1] = (merged[-1][0], b, br)
            else: merged.append((a, b, br))
        for a, b, br in merged:
            out.append({'chrom': chrom, 'hap': h, 'start_frac': round((a - x0) / L, 4), 'end_frac': round((b - x0) / L, 4), 'breed': br, 'bar_pt': round(L, 1)})
json.dump(out, open('embark_cosmo_segments.json', 'w'), indent=0)
# self-check: genome-wide proportions weighted by bar length (∝ chromosome length)
tot = defaultdict(float); allL = 0
for s in out:
    w = (s['end_frac'] - s['start_frac']) * s['bar_pt']; tot[s['breed']] += w; allL += w
print('segments:', len(out), '| chromosomes:', len({s['chrom'] for s in out}))
for br, w in sorted(tot.items(), key=lambda kv: -kv[1]): print(f'  {br:24s} {100*w/allL:5.1f}%   (Embark states: ' + {'POODLE_STANDARD':'40.2','POODLE_SMALL':'34.9','LABRADOR_RETRIEVER':'10.4','COCKER_SPANIEL':'9.4','ENGLISH_COCKER_SPANIEL':'5.1'}[br] + '%)')
