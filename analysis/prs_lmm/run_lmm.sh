#!/usr/bin/env bash
# Run the GEMMA sweep: one GRM, then one univariate LMM per (trait x fold)
# phenotype column prepared by build_inputs.py.
#
#   bash run_lmm.sh WORKDIR BFILE_PREFIX [JOBS]
#
# ~96 runs at n=2,895 x 131k SNPs; a few minutes each, JOBS in parallel.
set -uo pipefail
WORK=$1; BFILE=$2; JOBS=${3:-4}
GEMMA=$HOME/micromamba/envs/gemma/bin/gemma
cd "$WORK"

# GEMMA writes into ./output relative to cwd.
if [ ! -s output/grm.cXX.txt ]; then
  echo "=== GRM (centered) ==="
  # The GRM must be computed on ALL dogs regardless of phenotype missingness,
  # so every fold's run subsets the same matrix: use a dummy all-1 phenotype.
  awk '{print 1}' "$BFILE.fam" > pheno_dummy.txt
  "$GEMMA" -bfile "$BFILE" -p pheno_dummy.txt -gk 1 -o grm > grm.log 2>&1 \
    || { tail -5 output/grm.log grm.log 2>/dev/null; exit 1; }
fi

run_one() {
  local trait_slug=$1 col=$2
  local out="lmm_${trait_slug}_col${col}"
  [ -s "output/${out}.assoc.txt" ] && return 0
  "$GEMMA" -bfile "$BFILE" -p pheno.tsv -n "$col" \
           -k output/grm.cXX.txt -c covar.txt \
           -lmm 1 -o "$out" > "log_${out}.txt" 2>&1 \
    || echo "FAILED ${out} (see $WORK/log_${out}.txt)"
}

echo "=== LMM sweep ==="
i=0
while IFS=$'\t' read -r trait fold col; do
  slug=$(echo "${trait}_${fold}" | tr ' /' '__' | tr -cd 'A-Za-z0-9_')
  run_one "$slug" "$col" &
  i=$((i+1))
  if [ $((i % JOBS)) -eq 0 ]; then wait; fi
done < columns.tsv
wait
echo "done: $(ls output/*.assoc.txt 2>/dev/null | wc -l | tr -d ' ') assoc files"
