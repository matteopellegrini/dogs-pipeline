#!/usr/bin/env python3
"""Per catalogue variant: how many dogs carry it among customers, the imputation reference panel and public genomes.

Reads the outputs of analysis/omia/catalogue_public_carriers.sh (copied to DIR) plus bp_sites.tsv:
    python3 analysis/omia/catalogue_public_carriers.py DIR > DIR/catalogue_carriers.csv
A public-genome call needs DP >= 4. Carrier = at least one ALT allele. Customer counts are the live Panel v3 calls
(reads + breed prior, so a kit can be called without reads at the site); uncertain/held calls are left out.
"""
import csv, gzip, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path

DOGS = Path(__file__).resolve().parents[2]
DIR = Path(sys.argv[1])
MIN_DP = 4

CODE4 = json.load(open(DOGS / 'breed_panel_inputs' / 'harmonized_country.json'))
CODE4.update({'GESH': 'GERMAN_SHEPHERD', 'GLDR': 'GOLDEN_RETRIEVER', 'STPD': 'STANDARD_POODLE', 'MNSC': 'MINIATURE_SCHNAUZER',
              'MNPD': 'MINIATURE_POODLE', 'YRKT': 'YORKSHIRE_TERRIER', 'LBRT': 'LABRADOR_RETRIEVER', 'POMM': 'POMERANIAN',
              'TOLL': 'NOVA_SCOTIA_DUCK_TOLLING_RETRIEVER'})


def nice(label):
    return label.replace('_', ' ').title() if label else 'unlabelled'


def panel_breed(sample):
    if sample.isdigit():
        return 'US cohort dog'
    m = re.match(r'([A-Za-z]{4})', sample)
    return nice(CODE4.get(m.group(1).upper())) if m and CODE4.get(m.group(1).upper()) else 'Dog10K other'


# public-genome labels
labels = {}
for f in ('analysis/imputation_panel/pilot_v3_dogs.tsv', 'analysis/imputation_panel/panel_v3_gap_dogs.tsv'):
    for row in csv.reader(open(DOGS / f), delimiter='\t'):
        if row and not row[0].startswith('#'):
            labels[row[0]] = row[1]
for row in csv.reader(open(DOGS / 'breed_panel_inputs/sra_full_runs.txt'), delimiter='\t'):
    labels[row[0]] = row[1]
for row in csv.reader(open(DOGS / 'analysis/imputation_panel/pilot_existing_dogs.tsv'), delimiter='\t'):
    labels[row[0].removeprefix('SRA-')] = row[3]


def public_breed(dog):
    return nice(labels.get(dog.removeprefix('SRA-'), '')).replace('Pit Bull', 'Pit bull')


sites = {}      # key -> site_id
pos_site = defaultdict(list)
for sid, key, chrom, pos in csv.reader(open(DIR / 'sites.tsv'), delimiter='\t'):
    sites[key] = sid
    pos_site[(chrom, pos)].append(key)

panel_called, panel_car = Counter(), defaultdict(list)
for key, sample, gt in csv.reader(gzip.open(DIR / 'panel.tsv.gz', 'rt'), delimiter='\t'):
    if key not in sites:
        continue
    if sample == '__called__':
        panel_called[key] = int(gt)
    else:
        panel_car[key].append(sample)

pub_called, pub_car = Counter(), defaultdict(list)
for key, dog, gt, dp in csv.reader(gzip.open(DIR / 'slices.tsv.gz', 'rt'), delimiter='\t'):
    if key not in sites or '.' in gt or dp in ('.', '') or int(dp) < MIN_DP:
        continue
    pub_called[key] += 1
    if '1' in gt:
        pub_car[key].append(dog)
for chrom, pos, ref, alt, dog, gt, dp in csv.reader(gzip.open(DIR / 'gvcf.tsv.gz', 'rt'), delimiter='\t'):
    if dp in ('.', '') or int(dp) < MIN_DP or '.' in gt:
        continue
    for key in pos_site.get((chrom, pos), []):
        if alt == '.':
            pub_called[key] += 1                     # hom-ref block at the site
        elif f'{chrom}:{pos}:{ref}:{alt}' == key or alt.split(',')[0] == key.split(':')[3] and ref == key.split(':')[2]:
            pub_called[key] += 1
            if '1' in gt:
                pub_car[key].append(dog)
        else:
            pub_called[key] += 1                     # another allele at the site: not this variant

bp = {r['site_id']: r for r in csv.DictReader(open(DIR / 'bp_sites.tsv'), delimiter='\t')}
key_of = {sid: key for key, sid in sites.items()}


def top(names, n=4):
    return '; '.join(f'{b} {c}' for b, c in Counter(names).most_common(n))


w = csv.writer(sys.stdout)
w.writerow(['site_id', 'gene', 'condition', 'report_class', 'kind', 'variant', 'documented_breeds',
            'customers_called', 'customers_with_reads', 'customer_carriers',
            'reference_panel_dogs_called', 'reference_panel_carriers', 'reference_panel_carrier_breeds',
            'public_genomes_called', 'public_genome_carriers', 'public_carrier_breeds',
            'total_dogs_called', 'total_carriers'])
for sid, r in bp.items():
    key = key_of.get(sid)
    cust_car = sum(int(r[c] or 0) for c in ('carrier', 'two_copies', 'hemizygous'))
    cust_n = cust_car + int(r['clear'] or 0)       # uncertain and held calls are not counted as called
    pc = panel_called.get(key, 0) if key else 0
    pcar = panel_car.get(key, []) if key else []
    uc = pub_called.get(key, 0) if key else 0
    ucar = pub_car.get(key, []) if key else []
    variant = f"{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}"
    w.writerow([sid, r['gene'], r['condition'], r['report_class'], r['kind'], variant, r['documented_breeds'],
                cust_n, r['kits_with_reads'], cust_car,
                pc if pc else 'not on panel', len(pcar), top([panel_breed(s) for s in pcar]),
                uc, len(ucar), top([public_breed(d) for d in ucar]),
                cust_n + pc + uc, cust_car + len(pcar) + len(ucar)])
