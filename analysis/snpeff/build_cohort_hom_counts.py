#!/usr/bin/env python3
"""dogs-app lib/cohortHomCounts.json: per gene, how many cohort dogs carry a two-copy predicted HIGH-impact change
that is not common (Dog10K AF < 5% or unknown) — the same rows OmiaTable PredictedVariants considers. A gene hit in
5+ dogs is treated as annotation noise and not listed. The cohort is the DOGS-Gen study dogs (results/dogs-gen-*),
the 96 dogs the threshold was set on.
  build_cohort_hom_counts.py <out.json> <functional_variants.json> [...]
"""
import collections, json, sys

out, files = sys.argv[1], sys.argv[2:]
cnt = collections.Counter()
for f in files:
    genes = {r['gene'] for r in json.load(open(f))['high_variants']
             if r['zygosity'] == 'hom_alt' and (r['af_dog10k'] is None or r['af_dog10k'] < 0.05)}
    cnt.update(genes)
json.dump({'cohort_size': len(files), 'hom_counts': dict(cnt.most_common())}, open(out, 'w'))
print(f'{out}: {len(files)} dogs, {len(cnt)} genes, {sum(v >= 5 for v in cnt.values())} hit in 5+ dogs')
