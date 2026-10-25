#!/usr/bin/env python3
"""
Map Hayward 2016 (cornell_canine, canFam3.1 array SNPs) onto the panel's canFam4 sites.

    python3 hayward_sites.py HAYWARD_PREFIX PARKER_BIM CHAIN OUT.tsv

Writes one row per Hayward SNP that lands on a panel site: hayward_index, panel
chrom, pos, and how Hayward's A1/A2 relate to the panel a1/a2 (same, swapped,
flipped strand, or both). Ambiguous A/T and C/G SNPs are kept only when the
alleles match without a strand flip, since a flip there is undetectable.
"""
import sys
from pyliftover import LiftOver

hw, pbim, chain, out = sys.argv[1:5]
comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
panel = {}
for l in open(pbim):
    c, _, _, pos, a1, a2 = l.split()[:6]
    panel[(c if c.startswith('chr') else 'chr' + c, int(pos))] = (a1, a2)
lo = LiftOver(chain)
n = n_lift = n_hit = 0
rel = {'same': 0, 'swap': 0, 'flip': 0, 'flipswap': 0, 'mismatch': 0, 'ambig_drop': 0}
with open(out, 'w') as fo:
    fo.write('hw_index\thw_id\tchrom\tpos\ta1\ta2\trelation\n')
    for i, l in enumerate(open(hw + '.bim')):
        n += 1
        c, rsid, _, pos, a1, a2 = l.split()[:6]
        if c == '39': c = 'X'
        r = lo.convert_coordinate('chr' + c, int(pos) - 1)
        if not r: continue
        n_lift += 1
        key = (r[0][0], r[0][1] + 1)
        if key not in panel: continue
        p1, p2 = panel[key]
        ambiguous = comp.get(a1) == a2
        # A/T and C/G SNPs: the array export does not say which strand, and the allele
        # letters cannot tell. Half of them came out flipped (checked against Parker/Dog10K
        # frequencies: r ~ 0), so they are dropped rather than guessed.
        if ambiguous: k = 'ambig_drop'
        elif (a1, a2) == (p1, p2): k = 'same'
        elif (a1, a2) == (p2, p1): k = 'swap'
        elif (comp.get(a1), comp.get(a2)) == (p1, p2): k = 'flip'
        elif (comp.get(a1), comp.get(a2)) == (p2, p1): k = 'flipswap'
        else: k = 'mismatch'
        rel[k] += 1
        if k in ('same', 'swap', 'flip', 'flipswap'):
            n_hit += 1
            fo.write(f'{i}\t{rsid}\t{key[0]}\t{key[1]}\t{p1}\t{p2}\t{k}\n')
print(f'hayward SNPs {n}, lifted {n_lift}, on panel sites with usable alleles {n_hit}; {rel}')
