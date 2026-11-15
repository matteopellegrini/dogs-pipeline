#!/usr/bin/env python3
"""Site overlap of the Dog10K imputation panel, the Dog Aging Project VCF and the Parker array
(site_lists.sbatch outputs in $D/scratch/panel_overlap), plus DAP imputation quality by frequency.
usage: overlap.py > overlap.txt"""
import gzip, json, os
from collections import defaultdict

D = '/u/project/pellegrini/gkislik/dogs'
O = f'{D}/scratch/panel_overlap'
CHRS = [f'chr{i}' for i in range(1, 39)] + ['chrX']

def mbin(f):
    m = min(f, 1 - f)
    return 'mono' if m == 0 else '<1%' if m < 0.01 else '1-5%' if m < 0.05 else '5-20%' if m < 0.2 else '20-50%'
BINS = ['mono', '<1%', '1-5%', '5-20%', '20-50%']

parker = defaultdict(set)
for line in open(f'{D}/COSMO/analysis/cosmo_parker_full.bim'):
    c, _, _, pos, a1, a2 = line.split()
    parker[c if c.startswith('chr') else 'chr' + c].add(int(pos))
omia = defaultdict(list)
for v in json.load(open(f'{D}/reference_json/omia_variants_v3.json'))['variants']:
    if v.get('ref') and v.get('alt') and v.get('pos'):
        omia[v['chrom']].append((int(v['pos']), v['ref'], v['alt'], v.get('category', '')))

tot = defaultdict(int)
d10k_bin, d10k_bin_in_dap = defaultdict(int), defaultdict(int)
dap_info = defaultdict(lambda: [0, 0.0, 0])   # n, sum INFO, n INFO>=0.8 by DAP allele-frequency bin
pk = defaultdict(int); om = defaultdict(int)
for c in CHRS:
    fd, fk = f'{O}/dap_{c}.tsv.gz', f'{O}/d10k_{c}.tsv.gz'
    if not (os.path.exists(fd) and os.path.exists(fk)):
        print('MISSING', c); continue
    dap, dap_pos = set(), set()
    for line in gzip.open(fd, 'rt'):
        pos, ref, alt, info, raf, af = line.rstrip('\n').split('\t')
        dap.add((int(pos), ref, alt)); dap_pos.add(int(pos))
        try:
            s = dap_info[mbin(float(af))]; i = float(info); s[0] += 1; s[1] += i; s[2] += i >= 0.8
        except ValueError:
            pass
    d10k_pos, n_d10k = set(), 0
    for line in gzip.open(fk, 'rt'):
        pos, ref, alt, ac, an = line.rstrip('\n').split('\t')
        k = (int(pos), ref, alt); n_d10k += 1; d10k_pos.add(k[0])
        b = mbin(int(ac) / int(an)) if an not in ('.', '0') else 'mono'
        d10k_bin[b] += 1
        if k in dap:
            d10k_bin_in_dap[b] += 1; tot['shared'] += 1
    tot['d10k'] += n_d10k; tot['dap'] += len(dap)
    for p in parker.get(c, ()):
        pk['parker'] += 1; a = p in d10k_pos; b = p in dap_pos
        pk['in_d10k'] += a; pk['in_dap'] += b; pk['in_both'] += a and b
    for pos, ref, alt, cat in omia.get(c, ()):
        om['omia'] += 1; om['in_d10k'] += pos in d10k_pos; om['in_dap'] += pos in dap_pos
        om['in_both'] += pos in d10k_pos and pos in dap_pos
    print('done', c, flush=True)

print(f"\nDog10K records {tot['d10k']:,}  DAP records {tot['dap']:,}  shared (same pos/REF/ALT) {tot['shared']:,}")
print(f"DAP-only {tot['dap']-tot['shared']:,}  Dog10K-only {tot['d10k']-tot['shared']:,}")
print('\nDog10K sites by Dog10K MAF: n, share present in DAP')
for b in BINS:
    if d10k_bin[b]: print(f'  {b:7s} {d10k_bin[b]:>11,}  {d10k_bin_in_dap[b]/d10k_bin[b]:.3f}')
print('\nDAP imputation INFO by DAP allele frequency: n, mean INFO, share INFO>=0.8')
for b in BINS:
    n, s, h = dap_info[b]
    if n: print(f'  {b:7s} {n:>11,}  {s/n:.3f}  {h/n:.3f}')
print(f"\nParker array positions {pk['parker']:,}: in Dog10K {pk['in_d10k']:,}, in DAP {pk['in_dap']:,}, in both {pk['in_both']:,}")
print(f"OMIA catalogue SNV/indel positions {om['omia']:,}: in Dog10K {om['in_d10k']:,}, in DAP {om['in_dap']:,}, in both {om['in_both']:,}")
