#!/usr/bin/env python3
"""
Breed panel source from DIRECT read counts of the SRA top-up genomes (no imputation) -> npz on the Parker grid.

    python3 sra_ad_dosage.py AD_DIR SHEET RUNS.txt OUT.npz [--seed 1]

AD_DIR holds <RUN>.tsv.gz (CHROM POS REF ALT AD, cluster/sra-ad-extract.sh); RUNS.txt the runs to keep (the
IBS-deduplicated set from sra_dosage.py). Pseudo-haploid: at each panel site one read is drawn at random among
those carrying a1 or a2 and the dog is coded 0 or 2 copies of a1; -1 without such a read. Breed frequencies from
these calls are unbiased whatever the depth, and nothing is pulled toward the imputation reference — the GLIMPSE
dosages of the same genomes made their columns resemble every imputed customer dog.
"""
import gzip, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from merged_panel import read_parker

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
ad_dir, sheet, runs_path, out = sys.argv[1:5]
rng = np.random.default_rng(int(sys.argv[sys.argv.index('--seed') + 1]) if '--seed' in sys.argv else 1)
P, p_ids, p_truth, keys, _ = read_parker(f'{ROOT}/COSMO/analysis/cosmo_parker_full', f'{ROOT}/COSMO/analysis/parker_clust.txt')
del P
key_index = {k: i for i, k in enumerate(keys)}
alleles_at = {}                                   # Parker-grid row -> (a1, a2) of the panel site
for l in list(open(f'{ROOT}/breed_panel/sites.tsv'))[1:]:
    c, pos, s1, s2 = l.rstrip('\n').split('\t'); alleles_at[key_index[(c, pos)]] = (s1, s2)
label = {}
for l in open(sheet):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'sample_id' or len(f) < 9: continue
    lab = f[8].split()[0] if f[8].strip() else ''
    label[os.path.basename(f[1])] = lab        # SRA spellings are mapped with build_panel_v2 --alias breed_aliases.json
runs = [r.strip() for r in open(runs_path) if r.strip() and os.path.exists(f'{ad_dir}/{r.strip()}.tsv.gz')]
dos = np.full((len(keys), len(runs)), -1, dtype=np.int8)
for j, r in enumerate(runs):
    for l in gzip.open(f'{ad_dir}/{r}.tsv.gz', 'rt'):
        c, pos, ref, alt, ad = l.rstrip('\n').split('\t')[:5]
        i = key_index.get((c, pos))
        if i not in alleles_at: continue
        n = dict(zip([ref] + alt.split(','), (int(x) if x.isdigit() else 0 for x in ad.split(','))))
        s1, s2 = alleles_at[i]; k1, k2 = n.get(s1, 0), n.get(s2, 0)
        if k1 + k2 == 0: continue
        dos[i, j] = 2 if rng.random() < k1 / (k1 + k2) else 0
print(f'{len(runs)} genomes; sites called per genome median {int(np.median((dos >= 0).sum(0)))}', file=sys.stderr)
np.savez_compressed(out, dos=dos, ids=np.array(runs), labels=np.array([label[r] for r in runs]))
