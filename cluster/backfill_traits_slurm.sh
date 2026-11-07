#!/bin/bash
# Stage 13b backfill on Slurm: rewrite traits_result.json for published kits with the current trait catalogue.
# The previous file is kept once as traits_result.pre-$TAG.json (the republish key for cluster/republish_files.sh).
#   sbatch -A pi_matteop -p campus -a 1-99%40 -c 1 --mem=4G -t 3:00:00 -J traits25 \
#          cluster/backfill_traits_slurm.sh TAG LIST 20
# LIST rows: pk, barcode, pub dir, work dir (as logs/coat_b1101_list.tsv); kits sharing a pub dir run once.
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"; cd "$D"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONIOENCODING=utf-8
TAG=$1; LIST=$2; PER=${3:-20}
start=$(( (SLURM_ARRAY_TASK_ID - 1) * PER + 1 )); end=$(( start + PER - 1 ))
ok=0; fail=0; skip=0; seen=""
while IFS=$'\t' read -r pk barcode pub work; do
  [[ -n "$pk" && -d "$pub" ]] || continue
  case " $seen " in *" $pub "*) continue;; esac; seen="$seen $pub"
  [[ -f "$pub/traits_result.pre-$TAG.json" && "$pub/traits_result.json" -nt "$pub/traits_result.pre-$TAG.json" ]] && { skip=$((skip+1)); continue; }
  bam="$work/markdup.bam"; [[ -s "$bam" ]] || bam="$work/sites.bam"
  [[ -s "$bam" ]] || { echo "TRAITS-NOBAM $pk"; fail=$((fail+1)); continue; }
  bcf=$(ls "$work"/glimpse2/*_imputed_dog10k.bcf 2>/dev/null | head -1)
  sex=$(python3 -c "import json,sys;print(json.load(open(sys.argv[1]))['_meta'].get('predicted_sex','unknown'))" "$pub/coverage_1mb.json" 2>/dev/null || echo unknown)
  [[ -f "$pub/traits_result.json" && ! -f "$pub/traits_result.pre-$TAG.json" ]] && cp -p "$pub/traits_result.json" "$pub/traits_result.pre-$TAG.json"
  if envs/genomics/bin/python3 analysis/traits/trait_stage.py reference_json/trait_catalog.json "$pub/traits_result.json.tmp" \
       $([[ -n "$bcf" ]] && echo --bcf "$bcf") --bam "$bam" --sex "$sex" > /dev/null 2>&1; then
    mv "$pub/traits_result.json.tmp" "$pub/traits_result.json"; ok=$((ok+1))
  else rm -f "$pub/traits_result.json.tmp"; echo "TRAITS-FAIL $pk"; fail=$((fail+1)); fi
done < <(sed -n "${start},${end}p" "$LIST")
echo "TRAITS-BACKFILL-DONE task=$SLURM_ARRAY_TASK_ID ok=$ok fail=$fail skip=$skip"
