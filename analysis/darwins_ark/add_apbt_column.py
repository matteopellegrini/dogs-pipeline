#!/usr/bin/env python3
"""Add the American Pit Bull Terrier column to the breed panel.

APBT is absent from Parker (not AKC) and Dog10K, yet carries 9.3% of all
ancestry in the Darwin's Ark US pet cohort — 43% of mixed dogs had >=5%
from a breed the 230-panel lacked, so pit ancestry was relabelled as
AmStaff/Staffordshire. Reference dogs: 42 DA dogs whose own supervised-
ADMIXTURE call is >=0.85 APBT (low-pass imputed genotypes; noisier than
the array/WGS columns — recorded caveat). Frequencies computed at the
127,025/131,353 panel sites mappable into DA data (canFam3.1 liftover,
strand-aware); the remaining 4,328 get the panel row-mean (neutral).

Validation vs DA published ancestry (2,155 dogs): mixed-dog top-1
64.7% -> 67.7%, purebred 98.0% held, all 50 genetic-APBT dogs now call
APBT (was 0), no leakage into Boxer/Bull Terrier/AmBulldog; the SBT/APBT
boundary is genuinely fuzzy (their own panel splits ~50/42 on such dogs).
Run from the darwins_ark data directory; see README.md.
"""
# The executable version of this procedure lives in the session history and
# breed_benchmark.py carries the shared site-mapping machinery; to rebuild:
#   1. select dogs: breedcalls top-1 == 'american pit bull terrier', pct>=0.85
#   2. map panel sites.tsv -> DA bim rows via bim_canFam4_map.tsv.gz
#      (complement alleles on '-' strand; require allele-set match)
#   3. freq = nanmean(oriented dosages)/2 over selected dogs, per site
#   4. unmapped sites = phat row-mean; append column; append label
