#!/usr/bin/env python3
"""
Per-run allele counts and diploid calls at the breed-panel sites, from the full-depth SRA scan (v2.3).

    python3 sra_full_counts.py SITES.tsv RUNS.txt SRA_FULL_DIR OUT.npz

SRA_FULL_DIR holds <RUN>.panel.tsv.gz from cluster/sra-full-genotype.sh (CHROM POS REF ALT GT GQ AD DP; REF/ALT
from canFam4, hom-ref rows have no GQ). For each run and each site of SITES.tsv (chrom pos a1 a2) this keeps:
  k1, k2  reads carrying a1 / a2 (AD), for pseudo-haploid draws and for pooling split runs of one biosample;
  dip     copies of a1 from a confident diploid call (DP >= 8 and, unless hom-ref, GQ >= 20), else -1.
Turning these into one genotype per dog is sra_full_dosage.py.
"""
import gzip, sys
import numpy as np

sites_path, runs_path, d, out = sys.argv[1:5]
idx, a1, a2 = {}, [], []
with open(sites_path) as fh:
    next(fh)
    for i, l in enumerate(fh):
        c, pos, s1, s2 = l.rstrip('\n').split('\t'); idx[(c, pos)] = i; a1.append(s1); a2.append(s2)
runs, labels = [], []
for l in open(runs_path):
    f = l.rstrip('\n').split('\t')
    if f and f[0]: runs.append(f[0]); labels.append(f[1] if len(f) > 1 else '')
n = len(idx)
k1 = np.zeros((n, len(runs)), np.uint16); k2 = np.zeros((n, len(runs)), np.uint16)
dip = np.full((n, len(runs)), -1, np.int8)
for j, r in enumerate(runs):
    for l in gzip.open(f'{d}/{r}.panel.tsv.gz', 'rt'):
        c, pos, ref, alt, gt, gq, ad, dp = l.rstrip('\n').split('\t')[:8]
        i = idx.get((c, pos))
        if i is None: continue
        bases = [ref] + ([] if alt == '.' else alt.split(','))
        cnt = dict(zip(bases, (int(x) if x.isdigit() else 0 for x in ad.split(','))))
        k1[i, j] = min(cnt.get(a1[i], 0), 65535); k2[i, j] = min(cnt.get(a2[i], 0), 65535)
        g = gt.replace('|', '/').split('/')
        if len(g) != 2 or '.' in g or not dp.isdigit() or int(dp) < 8: continue
        if g != ['0', '0'] and not (gq.isdigit() and int(gq) >= 20): continue
        al = [bases[int(x)] if int(x) < len(bases) else None for x in g]
        if any(a not in (a1[i], a2[i]) for a in al): continue      # third allele: no call
        dip[i, j] = sum(a == a1[i] for a in al)
    print(f'{j + 1}/{len(runs)} {r}', file=sys.stderr, flush=True)
np.savez_compressed(out, k1=k1, k2=k2, dip=dip, runs=np.array(runs), labels=np.array(labels))
