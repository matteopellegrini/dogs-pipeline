#!/usr/bin/env python3
"""Summarise score_panels.py outputs: per depth x panel, concordance, non-ref concordance, r2 by Dog10K MAF.
usage: summarise.py SCORE_TSV [...] > summary.tsv"""
import sys
from collections import defaultdict
MAFS = ['<1%', '1-5%', '5-20%', '20-50%']
c = defaultdict(int); r2 = defaultdict(lambda: [0.0] * 6)
for fn in sys.argv[1:]:
    for line in open(fn):
        x = line.rstrip('\n').split('\t')
        if x[0] == 'C':
            c[(x[1], x[2], int(x[3]), int(x[4]))] += int(x[5])
        else:
            s = r2[(x[1], x[2])]
            for k in range(6): s[k] += float(x[3 + k])
def rsq(run, mb):
    n, sx, sy, sxy, sxx, syy = r2[(run, mb)]
    cov = sxy / n - sx / n * sy / n; vx = sxx / n - (sx / n) ** 2; vy = syy / n - (sy / n) ** 2
    return cov * cov / (vx * vy) if vx > 0 and vy > 0 else float('nan')
runs = sorted({k[0] for k in c}, key=lambda r: (float(r.split('_')[0][:-1]), r))
print('depth\tpanel\tsites\tconcordance\tnonref_conc\terrors_per_100k\t' + '\t'.join('r2 MAF ' + m for m in MAFS))
for run in runs:
    tot = sum(n for (r, mb, t, i), n in c.items() if r == run)
    ok = sum(n for (r, mb, t, i), n in c.items() if r == run and t == i)
    nr = sum(n for (r, mb, t, i), n in c.items() if r == run and (t or i))
    nrok = sum(n for (r, mb, t, i), n in c.items() if r == run and (t or i) and t == i)
    d, p = run.split('_', 1)
    print(f'{d}\t{p}\t{tot}\t{ok/tot:.5f}\t{nrok/nr:.4f}\t{(tot-ok)/tot*1e5:.0f}\t' + '\t'.join(f'{rsq(run, m):.3f}' for m in MAFS))
