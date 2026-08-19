#!/usr/bin/env bash
set -uo pipefail
JOBS=${1:-8}
B=GeneticData/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155
GEMMA=$HOME/micromamba/envs/gemma/bin/gemma
i=0
while IFS=$'\t' read -r trait fold col; do
  out="${trait}_${fold}"
  [ -s "output/${out}.assoc.txt" ] && continue
  "$GEMMA" -bfile "$B" -p gemma/pheno.tsv -n "$col" \
    -k gemma/grm_plink.rel -c "gemma/covar_${trait}.txt" -maf 0.05 \
    -lmm 1 -o "$out" > "gemma/log_${out}.txt" 2>&1 &
  i=$((i+1))
  if [ $((i % JOBS)) -eq 0 ]; then wait; fi
done < gemma/queue.tsv
wait
echo "DONE: $(ls output/*.assoc.txt | wc -l | tr -d ' ') assoc files"
