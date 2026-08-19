#!/usr/bin/env python3
"""Score one dog's GLIMPSE-imputed BCF against a canFam4 weights file
(Darwin's Ark individual-level PRS; py3.6-compatible for Hoffman).

    score_dog.py WEIGHTS.tsv[.gz] BCF

Weights: chr pos effect_allele other_allele beta af. Sites the dog lacks
(or with mismatched alleles) contribute beta * 2af — the training-cohort
mean — keeping scores comparable across dogs with different imputation
coverage. Prints: prs<TAB>matched/total."""
import gzip, os, subprocess, sys, tempfile

wts, bcf = sys.argv[1], sys.argv[2]
op = gzip.open if wts.endswith('.gz') else open
rows = [l.split() for l in op(wts, 'rt')]
bed = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
for r in rows:
    bed.write('%s\t%d\t%s\n' % (r[0], int(r[1]) - 1, r[1]))
bed.close()
q = subprocess.run(['bcftools', 'query', '-R', bed.name,
                    '-f', '%CHROM\t%POS\t%REF\t%ALT\t[%GP]\n', bcf],
                   stdout=subprocess.PIPE, universal_newlines=True)
os.unlink(bed.name)
dose = {}
for l in q.stdout.splitlines():
    f = l.split('\t')
    if len(f) != 5:
        continue
    gp = f[4].split(',')
    if len(gp) != 3:
        continue
    try:
        g = [float(x) for x in gp]
    except ValueError:
        continue
    dose[(f[0], f[1])] = (f[2], f[3], g[1] + 2 * g[2])
prs, matched = 0.0, 0
for c, p, ea, oa, beta, af in rows:
    beta, af = float(beta), float(af)
    hit = dose.get((c, p))
    if hit and set([ea, oa]) == set([hit[0], hit[1]]):
        d = hit[2] if ea == hit[1] else 2.0 - hit[2]
        matched += 1
    else:
        d = 2.0 * af
    prs += beta * d
print('%.4f\t%d/%d' % (prs, matched, len(rows)))
