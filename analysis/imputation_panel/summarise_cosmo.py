#!/usr/bin/env python3
"""Summarise score_chrom.py outputs for the old-vs-v2 Cosmo runs (score_cosmo.sbatch, gap scoring).
Leave-own-source-out truth: COSMO3* runs vs noGEN50, admera vs noADMERA.
usage: summarise_cosmo.py OUTDIR [OUTDIR...]"""
import glob, sys
from collections import defaultdict
MAFS = ['<1%', '1-5%', '5-10%', '10-20%', '20-50%']
c = defaultdict(int); r2 = defaultdict(lambda: [0.0] * 6)
for d in sys.argv[1:]:
    for f in glob.glob(d + '/chr*.counts.tsv'):
        for line in open(f):
            run, ts, mb, conf, t, i, n = line.split('\t')
            c[(run, ts, mb, conf, int(t), int(i))] += int(n)
    for f in glob.glob(d + '/chr*.r2.tsv'):
        for line in open(f):
            x = line.split('\t'); s = r2[tuple(x[:3])]
            for k in range(6): s[k] += float(x[3 + k])
loo = lambda r: 'noADMERA' if 'admera' in r else 'noGEN50'
def rsq(run, ts, mb):
    n, sx, sy, sxy, sxx, syy = r2[(run, ts, mb)]
    if n < 2: return float('nan')
    cov = sxy / n - sx / n * sy / n; vx = sxx / n - (sx / n) ** 2; vy = syy / n - (sy / n) ** 2
    return cov * cov / (vx * vy) if vx > 0 and vy > 0 else float('nan')
runs = sorted({k[0] for k in c})
print('run\ttruth\tsites\tconcordance\tnonref_conc\terrors_per_100k\t' + '\t'.join('r2 MAF ' + m for m in MAFS))
for run in runs:
    ts = loo(run)
    sel = {k: n for k, n in c.items() if k[0] == run and k[1] == ts}
    tot = sum(sel.values()); ok = sum(n for k, n in sel.items() if k[4] == k[5])
    nr = sum(n for k, n in sel.items() if k[4] or k[5]); nrok = sum(n for k, n in sel.items() if (k[4] or k[5]) and k[4] == k[5])
    if not tot: continue
    print(f'{run}\t{ts}\t{tot}\t{ok/tot:.5f}\t{nrok/nr:.4f}\t{(tot-ok)/tot*1e5:.0f}\t' + '\t'.join(f'{rsq(run, ts, m):.3f}' for m in MAFS))
