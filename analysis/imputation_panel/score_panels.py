#!/usr/bin/env python3
"""Score the panel-test imputations of Cosmo against the COSMO-SUPER read truth (cosmo_truth/out pileups).

usage: score_panels.py CHROM OUT_TSV RUN=BCF [RUN=BCF ...]

Same truth rule as cosmo_truth/score_chrom.py (>=10 reads, GQ>=20, leave-own-source-out: COSMO3 depths are
scored without the GEN50 reads, the Admera run without the Admera reads), but every run is binned by the
Dog10K allele frequency at the site, so the panels are compared on identical sites and bins.
TRUTH_PILEUP / TRUTH_SAMPLES / TRUTH_DROP override the Cosmo truth (e.g. an SRA dog from sra_truth.py).
Writes OUT_TSV: run maf truth imputed n, plus dosage sums for r2 (n sx sy sxy sxx syy) per run x maf.
"""
import gzip, math, os, subprocess, sys
from collections import defaultdict

chrom, out = sys.argv[1], sys.argv[2]
runs = [a.split('=', 1) for a in sys.argv[3:]]
D = '/u/project/pellegrini/gkislik/dogs'
BCFTOOLS = f'{D}/envs/genomics/bin/bcftools'
SAMPLES = os.environ.get('TRUTH_SAMPLES', 'ADMERA,UCLA,GEN50,OTHER').split(',')
LE, L1E, LH = math.log(0.01), math.log(0.99), math.log(0.5)

def truth(ref, alt):
    n = ref + alt
    if n < 10:
        return None
    ll = sorted(((alt * LE + ref * L1E, 0), (n * LH, 1), (ref * LE + alt * L1E, 2)), reverse=True)
    return ll[0][1] if ll[0][0] - ll[1][0] >= math.log(100) else None

pile = {}
with gzip.open(os.environ.get('TRUTH_PILEUP', f'{D}/analysis/cosmo_truth/out/{chrom}.pileup.tsv.gz'), 'rt') as f:
    for line in f:
        c = line.rstrip('\n').split('\t')
        pile[int(c[1])] = (c[2], c[3].split(','), [[int(x) if x != '.' else 0 for x in ad.split(',')] for ad in c[4:]])

raf = {}
with gzip.open(f'{D}/scratch/panel_test/ref/raf_{chrom}.tsv.gz', 'rt') as f:
    for line in f:
        pos, ref, alt, ac, an = line.split()
        raf[(int(pos), ref, alt)] = int(ac) / int(an)

def truth_at(pos, ref, alt, drop):
    p = pile.get(pos)
    if not p or p[0] != ref:
        return None
    ia = p[1].index(alt) + 1 if alt in p[1] else None
    r = a = o = 0
    for s, ad in zip(SAMPLES, p[2]):
        if s == drop:
            continue
        rr = ad[0]; aa = ad[ia] if ia is not None and ia < len(ad) else 0
        r += rr; a += aa; o += sum(ad) - rr - aa
    return None if o > 0.1 * (r + a + o) else truth(r, a)

def mafbin(f):
    m = min(f, 1 - f)
    for lim, lab in ((0.01, '<1%'), (0.05, '1-5%'), (0.2, '5-20%')):
        if m < lim:
            return lab
    return '20-50%'

counts = defaultdict(int)
r2 = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0, 0.0])
for run, bcf in runs:
    drop = os.environ.get('TRUTH_DROP') or ('ADMERA' if run.startswith('5.3x') else 'GEN50')
    q = subprocess.Popen([BCFTOOLS, 'query', '-f', '%POS\t%REF\t%ALT[\t%GT\t%DS]\n', bcf],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    for line in q.stdout:
        pos, ref, alt, gt, ds = line.rstrip('\n').split('\t')
        pos = int(pos)
        f = raf.get((pos, ref, alt))
        if f is None or gt.startswith('.'):
            continue
        t = truth_at(pos, ref, alt, drop)
        if t is None:
            continue
        imp = int(gt[0]) + int(gt[2]); ds = float(ds); mb = mafbin(f)
        counts[(run, mb, t, imp)] += 1
        s = r2[(run, mb)]
        s[0] += 1; s[1] += ds; s[2] += t; s[3] += ds * t; s[4] += ds * ds; s[5] += t * t
    q.wait()

with open(out, 'w') as fo:
    for k, v in sorted(counts.items()):
        fo.write('C\t' + '\t'.join(map(str, k + (v,))) + '\n')
    for k, s in sorted(r2.items()):
        fo.write('R\t' + '\t'.join(map(str, k + tuple(s))) + '\n')
print(chrom, 'scored', len(runs), 'runs')
