#!/usr/bin/env python3
"""Summarise seq_allele_caller.py counts per site: dogs with reads, dogs with alt reads, and their breeds
(customers: top breed >= 25% from work_panel_v3/kit_inputs.json; public genomes: their SRA label).
  python3 analysis/omia/seqcall_summary.py COUNTS.tsv
"""
import csv, json, sys
from collections import Counter, defaultdict
from pathlib import Path

DOGS = Path(__file__).resolve().parents[2]
sys.argv += []
kits = json.load(open(DOGS / 'work_panel_v3/kit_inputs.json'))['kits']
labels = {}
for f in ('analysis/imputation_panel/pilot_v3_dogs.tsv', 'analysis/imputation_panel/panel_v3_gap_dogs.tsv'):
    for row in csv.reader(open(DOGS / f), delimiter='\t'):
        if row and not row[0].startswith('#'):
            labels[row[0]] = row[1]
for row in csv.reader(open(DOGS / 'breed_panel_inputs/sra_full_runs.txt'), delimiter='\t'):
    labels[row[0]] = row[1]
for row in csv.reader(open(DOGS / 'analysis/imputation_panel/pilot_existing_dogs.tsv'), delimiter='\t'):
    labels[row[0].removeprefix('SRA-')] = row[3]

def breed(s):
    if s.startswith('pk-') or s in kits:
        b = (kits.get(s) or {}).get('breeds') or {}
        top = max(b.items(), key=lambda x: x[1], default=(None, 0))
        return (top[0] if top[1] >= 0.25 else 'MIXED'), 'customer'
    return labels.get(s.removeprefix('SRA-'), 'unlabelled'), 'public'

stat = defaultdict(lambda: Counter())
car = defaultdict(list)
other = defaultdict(Counter)
for r in csv.DictReader(open(sys.argv[1]), delimiter='\t'):
    b, grp = breed(r['sample'])
    n_ref, n_alt = int(r['ref_reads']), int(r['alt_reads'])
    k = (r['site_id'], grp)
    if n_ref + n_alt:
        stat[k]['with_reads'] += 1
        stat[k]['reads'] += n_ref + n_alt
    if n_alt:
        stat[k]['alt_dogs'] += 1
        car[k].append(f"{r['sample']}({b} {n_alt}/{n_ref + n_alt})")
    for d in filter(None, r['detail'].split(',')):
        other[r['site_id']][d.split(':')[0]] += int(d.split(':')[1])
for k in sorted(stat):
    s = stat[k]
    print(f"{k[0]:14s} {k[1]:8s} dogs_with_reads={s['with_reads']:5d} mean_reads={s['reads'] / max(1, s['with_reads']):.1f} "
          f"dogs_with_alt={s['alt_dogs']}")
    for c in car[k][:60]:
        print('    ', c)
for sid, c in other.items():
    print(sid, 'other alleles/run lengths:', dict(c))
