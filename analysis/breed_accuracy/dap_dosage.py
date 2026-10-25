#!/usr/bin/env python3
"""
Dog Aging Project genotypes on the Parker grid + owner-reported purebred labels.

    python3 dap_dosage.py DAP_GT.txt DAP_SAMPLES.txt DOG_OVERVIEW.csv LABEL_MAP.json PARKER_PREFIX CLUST.txt OUT.npz

Keeps only dogs the owner called purebred with a breed that maps to a canonical
label. Labels are NOT trusted yet: build_panel_v2.py --dap-confirm projects each
DAP dog onto a panel built WITHOUT DAP and keeps it as reference only when it
comes back >= 85% its stated breed. That step removes mislabelled and mixed
dogs without touching anything else.
"""
import csv, json, sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from merged_panel import read_parker, read_dog10k

gt, samples, overview, lmap, pprefix, clust, out = sys.argv[1:8]
P, p_ids, p_truth, keys, a1 = read_parker(pprefix, clust)
key_index = {k: i for i, k in enumerate(keys)}
D, ids, _ = read_dog10k(gt, samples, key_index, a1)        # same text layout as the Dog10K extract
m = json.load(open(lmap))
lab = {}
for r in csv.DictReader(open(overview)):
    if r['Breed_Status'] == 'Purebred' and m.get(r['Breed']): lab[r['dog_id']] = m[r['Breed']]
keep = [j for j, i in enumerate(ids) if i in lab]
typed = (D[:, keep] >= 0).mean(0)
print(f'DAP: {len(ids)} genotyped, {len(keep)} owner-reported purebred with a usable label; median sites typed {np.median(typed):.3f}', file=sys.stderr)
np.savez_compressed(out, dos=D[:, keep], ids=np.array([ids[j] for j in keep]), labels=np.array([lab[ids[j]] for j in keep]))
