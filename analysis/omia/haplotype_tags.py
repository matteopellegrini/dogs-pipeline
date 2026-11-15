#!/usr/bin/env python3
"""Do a catalogue variant's carriers share a haplotype? Correlates the variant genotype (from seq_allele_caller.py
counts in deep public genomes) with every panel SNV within +-WINDOW, using the same dogs' panel-site genotypes
(sra_panel/<BS>.panel_sites.bcf). Prints the best tag SNVs (r^2, carrier/non-carrier allele frequencies).

  python3 haplotype_tags.py COUNTS.tsv SITE_ID CHROM POS [WINDOW]
"""
import csv, os, subprocess, sys
import numpy as np

counts, sid, chrom, pos = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
W = int(sys.argv[5]) if len(sys.argv) > 5 else 1_000_000
truth = {}
for r in csv.DictReader(open(counts), delimiter='\t'):
    if r['site_id'] != sid or r['sample'].startswith(('pk-', 'dogs-gen', 'COSMO')):
        continue
    a, n = int(r['alt_reads']), int(r['ref_reads']) + int(r['alt_reads'])
    if a >= 2 and a / n >= 0.15:
        truth[r['sample']] = 2 if a / n >= 0.85 and n >= 6 else 1
    elif a == 0 and n >= 6:
        truth[r['sample']] = 0
dogs = [d for d in truth if os.path.exists(f'sra_panel/{d}.panel_sites.bcf')]
print(f'{sid}: {len(dogs)} dogs with a confident genotype + window genotypes; carriers', sum(truth[d] > 0 for d in dogs))
geno = {}
for j, d in enumerate(dogs):
    out = subprocess.run(['bcftools', 'query', '-r', f'{chrom}:{max(1, pos - W)}-{pos + W}', '-i', 'FMT/DP>=3',
                          '-f', '%POS\t%REF\t%ALT\t[%GT]\n', f'sra_panel/{d}.panel_sites.bcf'],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        p, ref, alt, gt = line.split('\t')
        if '.' in gt.replace('./.', '.') and gt in ('./.', '.'):
            continue
        g = gt.replace('|', '/')
        dose = 0 if g == '0/0' else 1 if g in ('0/1', '1/0') else 2 if g == '1/1' else None
        if dose is None:
            continue
        k = (int(p), alt if alt != '.' else None)
        geno.setdefault(int(p), {})[d] = (dose, alt)
y = np.array([truth[d] for d in dogs], float)
res = []
for p, gd in geno.items():
    alts = {a for _, a in gd.values() if a != '.'}
    if len(alts) != 1 or len(gd) < 0.8 * len(dogs):
        continue
    idx = [j for j, d in enumerate(dogs) if d in gd]
    x = np.array([gd[dogs[j]][0] for j in idx], float)
    yy = y[idx]
    if x.std() == 0:
        continue
    r = np.corrcoef(x, yy)[0, 1]
    car = yy > 0
    res.append((r * r, r, p, x[car].mean() / 2, x[~car].mean() / 2, len(idx)))
res.sort(reverse=True)
print('r2     r      pos        alt_freq_carriers alt_freq_noncarriers n')
for t in res[:15]:
    print(f'{t[0]:.3f} {t[1]:+.3f} {t[2]:>10d} {t[3]:.2f} {t[4]:.3f} {t[5]}')
