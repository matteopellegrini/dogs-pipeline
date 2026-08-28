# Relatives reference artifact (relatives_ref/, gitignored — 267MB)

`relatives_ref/{geno.npy,meta.json.gz,keys.tsv}` — genotypes of 2,031 reference
dogs at the 131,353 Parker sites (breed_panel/sites.tsv), used by pipeline
Stage 9b for KING-robust kinship.

Sources:
- 1,929 Dog10K imputation-panel dogs: hard GTs extracted from
  dog10k_panel/AutoAndXPAR.Dog10K.phased.bcf (bcftools query -R sites);
  breed labels from Dog10K Table S1 (Meadows 2023 Genome Biology,
  13059_2023_3023_MOESM1_ESM.xlsx).
- 102 Prosper reference/cohort dogs: DS dosages from each
  work/<s>/analysis/glimpse2/*_imputed_dog10k.bcf, rounded to 0/1/2;
  -1 where a site is absent from the imputed BCF (~13k of 131k).

Extraction jobs: cluster-side one-offs (relatives/extract_{panel,cohort}.sh on
Hoffman). Lesson repeated: export OPENBLAS_NUM_THREADS=1 or bcftools dies on
many-core nodes under 4G vmem.

Method choice: KING-robust, NOT a GRM — global-AF GRM scores unrelated Boxers
at phi~0.9 from breed homogeneity alone. KING validation on known truth:
same-dog kit pairs phi=0.499 (theory 0.5), cross-platform same-dog 0.489,
known sib quartet 0.20-0.26, two greyhound littermate pairs 0.23-0.26,
Cosmo vs Luna -0.03 (unrelated), same-breed background p99 = 0.046 (below the
2nd-degree cutoff 0.0884). Thresholds in meta.json.gz.

To rebuild: rerun the extraction jobs, then the matrix-assembly steps recorded
in the session transcript (panel_G/cohort_C construction + Table S1 labels).
