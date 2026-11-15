#!/usr/bin/env python3
"""Match the Panel v3 catalogue to the published frequencies of 250 disease variants in 1,054,293 dogs
(Donner et al. 2023 PLoS Genet, doi:10.1371/journal.pgen.1010651; S1 = frequencies in all/mixed/purebred dogs,
S3 = alternative allele frequency by breed). S1/S3 are converted to TSV first with analysis/omia/xlsx_rows.py.

    python3 analysis/omia/donner2023_frequencies.py known_variants/donner2023 > known_variants/donner2023/catalogue_match.tsv
Matching: OMIA variant ID where catalogue_alleles_resolved.tsv has one; otherwise the paper's CanFam3.1 position lifted to
canFam4 within 25 bp of a catalogue site of the same gene (same chromosome, any gene within 2 bp).
"""
import csv, re, sys
from collections import defaultdict
from pathlib import Path
from pyliftover import LiftOver

DOGS = Path(__file__).resolve().parents[2]
D = Path(sys.argv[1])
LO = LiftOver(str(DOGS / 'known_variants/omia_download_2026-11-01/chains/canFam3ToCanFam4.over.chain.gz'))
f = lambda x: float(x) if x not in ('', None) else None

# S1: header on line 3, data after; columns by name with the three blocks (all, mixed, purebred)
lines = [l.rstrip('\n').split('\t') for l in open(D / 'S1.tsv')]
hdr = lines[2]
def idx(name, k=0):
    return [i for i, h in enumerate(hdr) if h == name][k]
s1 = []
for r in lines[3:]:
    if len(r) < 10 or not r[0].strip():
        continue
    r += [''] * (len(hdr) - len(r))
    g = lambda name, k=0: r[idx(name, k)]
    s1.append(dict(omia=r[0].strip(), phene=g('Variant phenotype'), moi=g('Mode of inheritance'), gene=g('Gene'),
                   variant=g('Variant'), build=g('Genome build'), rank=g('Prevalence Rank'),
                   af_all=f(g('Allele frequency', 0)), af_mixed=f(g('Allele frequency', 1)), af_pure=f(g('Allele frequency', 2)),
                   n_all=sum(int(float(g(c, 0) or 0)) for c in ('Count 0 copy dogs', 'Count 1 copy dogs', 'Count 2 copy dogs')),
                   car_all=sum(f(g(c, 0)) or 0 for c in ('Frequency 1 copy dogs', 'Frequency 2 copy dogs')),
                   car_mixed=sum(f(g(c, 1)) or 0 for c in ('Frequency 1 copy dogs', 'Frequency 2 copy dogs')),
                   car_pure=sum(f(g(c, 2)) or 0 for c in ('Frequency 1 copy dogs', 'Frequency 2 copy dogs'))))

# S3: row 3 = variant id / reference per column, then pairs of rows (N dogs genotyped, alternative allele frequency) per breed
s3 = [l.rstrip('\n').split('\t') for l in open(D / 'S3.tsv')]
colid = s3[2]
by_breed = defaultdict(list)          # column -> [(breed, n, af)]
n_row = {}
for r in s3[4:]:
    if len(r) < 3:
        continue
    if r[0] == 'N dogs genotyped':
        n_row = {i: r[i] for i in range(2, len(r))}
    elif r[0] == 'Alternative allele frequency':
        for i in range(2, len(r)):
            if r[i] not in ('', 'NA', 'N/A') and n_row.get(i) not in ('', None, 'NA', 'N/A'):
                by_breed[i].append((r[1], int(float(n_row[i])), float(r[i])))
col_of = {colid[i].strip(): i for i in range(2, len(colid)) if colid[i].strip().isdigit()}

# catalogue sites (canFam4) + OMIA ids
cat = list(csv.DictReader(open(DOGS / 'known_variants/panel_v3/panel_v3_sites.tsv'), delimiter='\t'))
omia_site = defaultdict(set)
for r in csv.DictReader(open(DOGS / 'known_variants/panel_v3/catalogue_alleles_resolved.tsv'), delimiter='\t'):
    if r.get('omia_variant_id'):
        omia_site[r['omia_variant_id']].add(r['site_id'])

out = csv.writer(sys.stdout, delimiter='\t')
out.writerow(['omia_variant_id', 'paper_phenotype', 'gene', 'paper_variant_canfam3', 'canfam4_pos', 'site_ids', 'match',
              'dogs_genotyped', 'allele_freq_all', 'carrier_freq_all', 'carrier_freq_mixed', 'carrier_freq_purebred',
              'carriers_per_10k_mixed', 'top_breeds_af(n>=50)', 'breeds_with_alt'])
for v in s1:
    m = re.match(r'chr(\w+):(\d+)', v['variant'].replace(' ', ''))
    pos4 = chrom = None
    if m and 'CanFam3' in v['build']:
        hit = LO.convert_coordinate('chr' + m.group(1), int(m.group(2)) - 1)
        if hit:
            chrom, pos4 = hit[0][0], hit[0][1] + 1
    sites, how = set(omia_site.get(v['omia'], ())), 'omia_id' if omia_site.get(v['omia']) else ''
    if pos4:
        near = []
        for c in cat:
            if c['chrom'] != chrom or not c['pos'].isdigit():
                continue
            d = abs(int(c['pos']) - pos4)
            if (d <= 25 and c['gene'].upper() == v['gene'].upper()) or d <= 2:
                near.append((d, c['site_id']))
        if near and not sites:                        # nearest catalogue site(s) only: two variants of one gene can sit 20 bp apart
            dmin = min(d for d, _ in near)
            sites = {s for d, s in near if d == dmin}
            how = 'position'
    col = col_of.get(v['omia'])
    br = by_breed.get(col, [])
    top = sorted((x for x in br if x[1] >= 50 and x[2] > 0), key=lambda x: -x[2])[:6]
    out.writerow([v['omia'], v['phene'], v['gene'], v['variant'], f'{chrom}:{pos4}' if pos4 else '', ','.join(sorted(sites)),
                  how or 'none', v['n_all'], f"{v['af_all']:.3g}" if v['af_all'] is not None else '',
                  f"{v['car_all']:.3g}", f"{v['car_mixed']:.3g}", f"{v['car_pure']:.3g}", round(v['car_mixed'] * 1e4, 1),
                  '; '.join(f'{b} {a:.3f}' for b, n, a in top), sum(1 for x in br if x[2] > 0)])
