#!/usr/bin/env python3
"""Score each Cosmo run's imputed genotypes against direct read calls from COSMO-SUPER (~18x).

usage: score_chrom.py CHROM PILEUP_TSV_GZ SAMPLES OUT_PREFIX RUN=BCF [RUN=BCF ...]

PILEUP_TSV_GZ: CHROM POS REF ALT AD... from bcftools mpileup on the COSMO-SUPER reads, one sample per
read source (SAMPLES, comma list, e.g. ADMERA,UCLA,GEN50,OTHER). Truth sets: ALL (every source) and
no<SRC> (leave one source out, so a run is never scored against its own reads).
Truth call: ref+alt >= 10, other alleles < 10%, binomial likelihoods (error 1%), GQ >= 20.
Writes OUT_PREFIX.counts.tsv (run truthset maf conf truth imputed n) and OUT_PREFIX.r2.tsv (dosage sums).
"""
import gzip, math, subprocess, sys
from collections import defaultdict

chrom, pileup, samples, out = sys.argv[1], sys.argv[2], sys.argv[3].split(','), sys.argv[4]
runs = [a.split('=', 1) for a in sys.argv[5:]]
BCFTOOLS = 'envs/genomics/bin/bcftools'
SETS = ['ALL'] + ['no' + s for s in samples if s != 'OTHER']
LE, L1E, LH = math.log(0.01), math.log(0.99), math.log(0.5)

def truth(ref, alt):
    n = ref + alt
    if n < 10:
        return None
    ll = sorted(((alt * LE + ref * L1E, 0), (n * LH, 1), (ref * LE + alt * L1E, 2)), reverse=True)
    return ll[0][1] if ll[0][0] - ll[1][0] >= math.log(100) else None

# pos -> (alts tuple, [AD list per sample])
pile = {}
with gzip.open(pileup, 'rt') as f:
    for line in f:
        c = line.rstrip('\n').split('\t')
        pile[int(c[1])] = (c[2], c[3].split(','), [[int(x) if x != '.' else 0 for x in ad.split(',')] for ad in c[4:]])

memo = {}
def truths(pos, ref, alt):
    key = (pos, alt)
    if key in memo:
        return memo[key]
    p = pile.get(pos)
    res = None
    if p and p[0] == ref:
        palts, ads = p[1], p[2]
        ia = palts.index(alt) + 1 if alt in palts else None
        per = []
        for ad in ads:
            r = ad[0]; a = ad[ia] if ia is not None and ia < len(ad) else 0
            per.append((r, a, sum(ad) - r - a))
        res = []
        for ts in SETS:
            use = [i for i, s in enumerate(samples) if ts == 'ALL' or s != ts[2:]]
            r = sum(per[i][0] for i in use); a = sum(per[i][1] for i in use); o = sum(per[i][2] for i in use)
            res.append(None if o > 0.1 * (r + a + o) else truth(r, a))
        res = tuple(res)
    memo[key] = res
    return res

def mafbin(raf):
    m = min(raf, 1 - raf)
    for lim, lab in ((0.01, '<1%'), (0.05, '1-5%'), (0.1, '5-10%'), (0.2, '10-20%')):
        if m < lim:
            return lab
    return '20-50%'

counts = defaultdict(int)
r2 = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0, 0.0])
for run, bcf in runs:
    q = subprocess.Popen([BCFTOOLS, 'query', '-r', chrom, '-f', '%POS\t%REF\t%ALT\t%INFO/RAF[\t%GT\t%DS\t%GP]\n', bcf],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    for line in q.stdout:
        pos, ref, alt, raf, gt, ds, gp = line.rstrip('\n').split('\t')
        if ',' in alt or gt.startswith('.'):
            continue
        t = truths(int(pos), ref, alt)
        if t is None:
            continue
        imp = int(gt[0]) + int(gt[2])
        ds = float(ds)
        conf = 'GP>=0.8' if max(float(x) for x in gp.split(',')) >= 0.8 else 'GP<0.8'
        mb = mafbin(float(raf))
        for ts, tg in zip(SETS, t):
            if tg is None:
                continue
            counts[(run, ts, mb, conf, tg, imp)] += 1
            s = r2[(run, ts, mb)]
            s[0] += 1; s[1] += ds; s[2] += tg; s[3] += ds * tg; s[4] += ds * ds; s[5] += tg * tg
    q.wait()

with open(out + '.counts.tsv', 'w') as f:
    for k, v in sorted(counts.items()):
        f.write('\t'.join(map(str, k + (v,))) + '\n')
with open(out + '.r2.tsv', 'w') as f:
    for k, s in sorted(r2.items()):
        f.write('\t'.join(map(str, k + tuple(s))) + '\n')
print(f'{chrom}: {len(pile)} pileup sites, {len(memo)} scored keys')
