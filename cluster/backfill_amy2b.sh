#!/bin/bash
# AMY2B copy number (analysis/traits/amy2b.py) for every kit in LIST -> OUT/<pk>.json (staging; not the pub dir).
# Install into pub dirs separately, after the cohort is validated.
#   sbatch -A pi_matteop -p campus -a 1-99%40 -c 1 --mem=2G -t 1:00:00 -J amy2b \
#          cluster/backfill_amy2b.sh work_panel_v3/manifest.tsv 20 work_amy2b/out
# LIST rows: pk, barcode, pub dir, work dir. Needs the full markdup.bam (sites.bam lacks the control windows).
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
D="${SLURM_SUBMIT_DIR:-$PWD}"; cd "$D"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
LIST=$1; PER=${2:-20}; O=${3:?out dir}; mkdir -p "$O"
T=${SLURM_ARRAY_TASK_ID:-1}
start=$(( (T - 1) * PER + 1 )); end=$(( start + PER - 1 ))
ok=0; fail=0; skip=0
while IFS=$'\t' read -r pk barcode pub work; do
  [[ -n "$pk" ]] || continue
  [[ -s "$O/$pk.json" ]] && { skip=$((skip+1)); continue; }
  bam="$work/markdup.bam"
  [[ -s "$bam" ]] || { echo "AMY2B-NOBAM $pk"; fail=$((fail+1)); continue; }
  if envs/genomics/bin/python3 analysis/traits/amy2b.py "$bam" "$O/$pk.json.tmp" --fasta canFam4_plusY.fa 2>/dev/null; then
    mv "$O/$pk.json.tmp" "$O/$pk.json"; ok=$((ok+1))
  else rm -f "$O/$pk.json.tmp"; echo "AMY2B-FAIL $pk"; fail=$((fail+1)); fi
done < <(sed -n "${start},${end}p" "$LIST")
echo "AMY2B-DONE task=$T ok=$ok fail=$fail skip=$skip"
