#!/usr/bin/env python3
"""Combine score_chrom.py outputs: per run, genotype concordance vs ALL and leave-own-source-out truth.
usage: aggregate.py OUTDIR > summary.tsv"""
import glob, sys
from collections import defaultdict

LOO = {'admera-0001': 'noADMERA', 'UCLA-20001': 'noUCLA', 'UCLA-6164': 'noUCLA'}
loo = lambda r: LOO.get(r, 'noGEN50' if r.startswith('COSMO3') else 'ALL')
ORDER = ['COSMO3-DS0.1', 'COSMO3-DS0.25', 'COSMO3-DS0.5', 'COSMO3-DS1.0', 'COSMO3', 'admera-0001', 'UCLA-20001', 'UCLA-6164', 'COSMO-SUPER']
MAFS = ['<1%', '1-5%', '5-10%', '10-20%', '20-50%']

c = defaultdict(int)
r2 = defaultdict(lambda: [0.0] * 6)
for f in glob.glob(sys.argv[1] + '/chr*.counts.tsv'):
    for line in open(f):
        run, ts, mb, conf, t, i, n = line.split('\t')
        c[(run, ts, mb, conf, int(t), int(i))] += int(n)
for f in glob.glob(sys.argv[1] + '/chr*.r2.tsv'):
    for line in open(f):
        x = line.split('\t'); s = r2[tuple(x[:3])]
        for k in range(6): s[k] += float(x[3 + k])

def stats(run, ts, mbs=None, conf=None):
    tot = ok = nr = nrok = hc = 0
    cls = defaultdict(lambda: [0, 0])
    for (r, t2, mb, cf, t, i), n in c.items():
        if r != run or t2 != ts or (mbs and mb not in mbs) or (conf and cf != conf): continue
        tot += n; ok += n * (t == i); cls[t][0] += n; cls[t][1] += n * (t == i)
        if t or i: nr += n; nrok += n * (t == i)
        if cf == 'GP>=0.8': hc += n
    return tot, ok / tot, nrok / nr, {k: v[1] / v[0] for k, v in cls.items()}, hc / tot

def rsq(run, ts, mb):
    n, sx, sy, sxy, sxx, syy = r2[(run, ts, mb)]
    if n < 2: return float('nan')
    cov = sxy / n - sx / n * sy / n; vx = sxx / n - (sx / n) ** 2; vy = syy / n - (sy / n) ** 2
    return cov * cov / (vx * vy) if vx > 0 and vy > 0 else float('nan')

print('run\ttruth\tsites\tconcordance\tnonref_conc\thomref\thet\thomalt\tpct_GP>=0.8\tconc_GP>=0.8\tconc_GP<0.8\t' + '\t'.join('r2 MAF ' + m for m in MAFS))
for run in ORDER:
    for ts in dict.fromkeys(['ALL', loo(run)]):
        tot, conc, nrc, cls, hc = stats(run, ts)
        hi = stats(run, ts, conf='GP>=0.8')[1]; lo = stats(run, ts, conf='GP<0.8')[1]
        print(f'{run}\t{ts}\t{tot}\t{conc:.4f}\t{nrc:.4f}\t{cls.get(0,0):.4f}\t{cls.get(1,0):.4f}\t{cls.get(2,0):.4f}\t{hc:.3f}\t{hi:.4f}\t{lo:.4f}\t'
              + '\t'.join(f'{rsq(run, ts, m):.3f}' for m in MAFS))
