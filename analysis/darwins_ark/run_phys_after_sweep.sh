#!/usr/bin/env bash
# Wait for the behavioural fold sweep (54 assoc files), then run the five
# physical-trait GWAS in parallel.
set -uo pipefail
B=GeneticData/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155
GEMMA=$HOME/micromamba/envs/gemma/bin/gemma
while [ "$(ls output/*.assoc.txt 2>/dev/null | wc -l)" -lt 54 ]; do sleep 300; done
echo "sweep complete, starting physical traits at $(date)"
while IFS=$'\t' read -r trait fold col; do
  out="${trait}_full"
  [ -s "output/${out}.assoc.txt" ] && continue
  "$GEMMA" -bfile "$B" -p gemma/pheno_phys.tsv -n "$col" \
    -k gemma/grm_plink.rel -c gemma/covar_phys.txt -maf 0.05 \
    -lmm 1 -o "$out" > "gemma/log_${out}.txt" 2>&1 &
done < gemma/columns_phys.tsv
wait
echo "physical traits done: $(ls output/*_full.assoc.txt | wc -l | tr -d ' ') full assoc files"
