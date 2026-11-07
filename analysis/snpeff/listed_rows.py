#!/usr/bin/env python3
"""Which functional_variants.json rows a report LISTS by default (OmiaTable PredictedVariants + lib/predictedVariants.ts),
and whether each sits inside the canFam4 span of the gene it names. Used to size a functional_variants rebuild.
  listed_rows.py <app lib dir> <gene spans tsv: symbol chrom start end> <out.tsv> <functional_variants.json> [...]
The app lib dir supplies geneDiseases.json, variantCategories.json, cohortHomCounts.json and panelHomLof.json, so the
same rule can be scored with the old or the rebuilt assets.
"""
import collections, json, os, sys

lib, spans_tsv, out, files = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
diseases = json.load(open(f'{lib}/geneDiseases.json'))['genes']
categories = json.load(open(f'{lib}/variantCategories.json'))
cohort = json.load(open(f'{lib}/cohortHomCounts.json'))['hom_counts']
panel_hom = set(json.load(open(f'{lib}/panelHomLof.json'))['hom'])
spans = collections.defaultdict(list)
for line in open(spans_tsv):
    g, c, s, e = line.rstrip('\n').split('\t')
    spans[g].append((c.replace('chr', ''), int(s), int(e)))

def less_likely(v):
    if v['af_dog10k'] is None or v['af_dog10k'] >= 0.01: return True
    if not any(a['moi'] in ('recessive', 'X-linked', '') for a in diseases.get(v['gene'], [])): return True
    return f"chr{v['chr']}:{v['pos']}:{v['alt']}" in panel_hom

rows, reports = [], 0
for f in files:
    listed = [v for v in json.load(open(f))['high_variants']
              if v['zygosity'] == 'hom_alt' and (v['af_dog10k'] is None or v['af_dog10k'] < 0.05)
              and f"chr{v['chr']}:{v['pos']}" not in categories and cohort.get(v['gene'], 0) < 5
              and v['gene'] in diseases and not less_likely(v)]
    reports += bool(listed)
    for v in listed:
        inside = any(c == v['chr'] and s <= int(v['pos']) <= e for c, s, e in spans.get(v['gene'], []))
        rows.append([os.path.basename(os.path.dirname(f)), v['gene'], v['chr'], v['pos'], v['alt'], v['effect'],
                     str(v['af_dog10k']), str(inside)])
with open(out, 'w') as fo:
    fo.write('kit\tgene\tchr\tpos\talt\teffect\taf_dog10k\tinside_gene\n')
    fo.writelines('\t'.join(r) + '\n' for r in rows)
print(f'{len(files)} reports: {len(rows)} listed rows on {reports} reports, {sum(r[-1] == "True" for r in rows)} inside the named gene')
