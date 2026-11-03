#!/usr/bin/env python3
"""
Breed panel v2.3 source: full-depth SRA genomes -> one genotype per dog on the Parker grid (npz: dos, ids, labels, diploid_frac).

    python3 sra_full_dosage.py COUNTS.npz BIOSAMPLES.tsv OUT.npz --dog10k GT.txt SAMPLES.txt [--seed 1]

COUNTS.npz is sra_full_counts.py output (per run: reads for a1/a2, confident diploid call). Runs of one biosample
(ENA sample_accession) are pooled: reads summed, diploid call taken from the deepest run that made one. Per site,
a confident diploid call is used when there is one; otherwise one read is drawn among those carrying a1 or a2
(pseudo-haploid, coded 0 or 2 as in sra_ad_dosage.py). diploid_frac is the share of called sites with a diploid
call; dogs below 0.5 count at half weight toward min-n (the v2.2 rule).

v2.3 SHIPS --pseudo-only. With diploid calls the panel failed the customer-kit gate (VILLAGE_Vietnam 0 -> 3 of 300
kits at >=5%; White Swiss Shepherd, American Bulldog, Dutch Shepherd up). Against the same animals in Dog10K, the
10-15x SRR7107xxx calls lose ~7% of hets to homozygous, leaning to one allele; --min-hom-depth 15 did not fix it.
Pseudo-haploid draws from the full-depth reads pass, and every SRA dog then counts at half weight toward min-n.

Duplicates of Parker / Dog10K dogs and within the set are removed by the rate of opposite homozygotes (0 vs 2)
over sites both call, using only this dog's DIPLOID calls (a pseudo-haploid draw at a het site looks homozygous,
so plain IBS under-reads identity at low depth). A pair is the same animal when the rate is <= 0.001 over
>= 50,000 shared sites AND under a third of the dog's next-closest match: distinct dogs of very inbred breeds
(Norwegian Lundehund, Berger Picard) sit at 0.002-0.004 with each other, and tiny split runs with ~8,000 diploid
sites matched dogs of other breeds at 0.001-0.004. Written to sra_full_duplicates.tsv next to OUT.
"""
import os, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from merged_panel import read_parker, read_dog10k

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
counts, bios_path, out = sys.argv[1:4]
args = sys.argv[4:]
opt = lambda k, d: args[args.index(k) + 1] if k in args else d
rng = np.random.default_rng(int(opt('--seed', 1)))
CUT = float(opt('--opp-hom', 0.001)); MIN_SITES = int(opt('--min-sites', 50000)); GAP = 3.0
# --pseudo-only: no diploid calls (v2.2-style draws from the full-depth reads). --min-hom-depth N: a homozygous diploid
# call needs N reads at the site; at 10-15x bcftools calls ~7% of true hets homozygous (vs Dog10K, same animals).
# Within the SRA set both sides are 10-15x diploid calls that miss some hets, which hides differences: two Norwegian
# Lundehunds from different biosamples (SRR7107939/940, an extreme-bottleneck breed) sat at 0.0006. Stricter cut there.
CUT_INTERNAL = float(opt('--opp-hom-internal', 0.0005))
PSEUDO_ONLY = '--pseudo-only' in args; MIN_HOM_DP = int(opt('--min-hom-depth', 0))

z = np.load(counts); k1, k2, dip = z['k1'], z['k2'], z['dip']; runs = list(z['runs']); labs = list(z['labels'])
bios = {}
for l in list(open(bios_path))[1:]:
    f = l.rstrip('\n').split('\t'); bios[f[0]] = f[1] or f[0]
groups = defaultdict(list)
for j, r in enumerate(runs): groups[bios.get(r, r)].append(j)
dogs = sorted(groups, key=lambda b: runs[groups[b][0]])

P, p_ids, p_truth, keys, a1 = read_parker(f'{ROOT}/COSMO/analysis/cosmo_parker_full', f'{ROOT}/COSMO/analysis/parker_clust.txt')
key_index = {k: i for i, k in enumerate(keys)}
grid = []                                          # Parker-grid row of each panel site, in sites.tsv order
with open(f'{ROOT}/breed_panel/sites.tsv') as fh:
    next(fh)
    for l in fh:
        c, pos = l.split('\t')[:2]; grid.append(key_index[(c, pos)])
grid = np.array(grid)

n = len(keys)
dos = np.full((n, len(dogs)), -1, np.int8); dipo = np.full((n, len(dogs)), -1, np.int8)
ids, labels, dfrac, depth = [], [], [], []
for j, b in enumerate(dogs):
    cols = groups[b]
    a = k1[:, cols].sum(1, dtype=np.int64); t = a + k2[:, cols].sum(1, dtype=np.int64)
    d = np.full(len(a), -1, np.int8)
    for c in sorted(cols, key=lambda c: -int((k1[:, c].astype(np.int64) + k2[:, c]).sum())):
        m = (d < 0) & (dip[:, c] >= 0); d[m] = dip[m, c]
    draw = np.where(t > 0, np.where(rng.random(len(a)) < a / np.maximum(t, 1), 2, 0), -1).astype(np.int8)
    dipo[grid, j] = d                                # the duplicate check always uses the diploid calls
    if PSEUDO_ONLY: d = np.full_like(d, -1)
    elif MIN_HOM_DP: d = np.where((d != 1) & (t < MIN_HOM_DP), -1, d).astype(np.int8)
    g = np.where(d >= 0, d, draw)
    dos[grid, j] = g
    called = (g >= 0).sum()
    ids.append(runs[cols[0]] if len(cols) == 1 else '+'.join(runs[c] for c in cols))
    labels.append(labs[cols[0]]); dfrac.append((d >= 0).sum() / max(called, 1)); depth.append(float(t[t > 0].mean()) if called else 0.0)
print(f'{len(runs)} runs -> {len(dogs)} dogs; median called sites {int(np.median((dos >= 0).sum(0)))}; '
      f'diploid_frac >=0.5: {sum(f >= 0.5 for f in dfrac)}', file=sys.stderr)

def opp_hom(Adip, B):
    """rate of 0-vs-2 over sites both call; A diploid calls only."""
    mA = Adip >= 0; mB = B >= 0
    n_both = mA.T.astype(np.float32) @ mB.astype(np.float32)
    opp = ((Adip == 0) & mA).astype(np.float32).T @ ((B == 2) & mB).astype(np.float32) + \
          ((Adip == 2) & mA).astype(np.float32).T @ ((B == 0) & mB).astype(np.float32)
    return opp / np.maximum(n_both, 1), n_both

D, d_ids, d_truth = read_dog10k(opt('--dog10k', '') , args[args.index('--dog10k') + 2], key_index, a1)
dup, rows = set(), []
for name, B, b_ids, b_lab in (('Parker', P, p_ids, lambda i: p_truth.get(i, '')), ('Dog10K', D, d_ids, lambda i: d_truth.get(i, ''))):
    r_, nb = opp_hom(dipo, B); r_[nb < MIN_SITES] = 1
    for i, k in np.argwhere(r_ <= CUT):
        if r_[i, k] * GAP > np.partition(r_[i], 1)[1]: continue          # not clearly apart from the next match
        dup.add(ids[i]); rows.append((name, ids[i], labels[i], b_ids[k], b_lab(b_ids[k]), r_[i, k], nb[i, k]))
r_, nb = opp_hom(dipo, dipo); r_[nb < MIN_SITES] = 1; np.fill_diagonal(r_, 1)
r2 = np.minimum(r_, r_.T)
for i, k in np.argwhere(np.triu(r2 <= CUT_INTERNAL, 1)):
    if ids[i] in dup or ids[k] in dup: continue
    nxt = min(np.partition(np.delete(r2[i], k), 0)[0], np.partition(np.delete(r2[k], i), 0)[0])
    if r2[i, k] * GAP > nxt: continue
    lo = i if depth[i] < depth[k] else k                     # keep the deeper copy
    dup.add(ids[lo]); rows.append(('SRA internal', ids[lo], labels[lo], ids[i + k - lo], labels[i + k - lo], r2[i, k], max(nb[i, k], nb[k, i])))
with open(os.path.join(os.path.dirname(os.path.abspath(out)), 'sra_full_duplicates.tsv'), 'w') as fh:
    fh.write('comparison\tsra_dog\tsra_label\tother_id\tother_label\topp_hom\tsites\n')
    for r in rows: fh.write('\t'.join(map(str, r[:5])) + f'\t{r[5]:.4f}\t{int(r[6])}\n')
keep = [j for j, i in enumerate(ids) if i not in dup]
print(f'duplicates removed: {len(set(dup))}; kept {len(keep)}', file=sys.stderr)
np.savez_compressed(out, dos=dos[:, keep], ids=np.array([ids[j] for j in keep]), labels=np.array([labels[j] for j in keep]),
                    diploid_frac=np.array([dfrac[j] for j in keep]), mean_depth=np.array([depth[j] for j in keep]))
