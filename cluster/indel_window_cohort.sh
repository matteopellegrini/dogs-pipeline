#!/bin/bash
# Impute one disease indel window (analysis/omia/indel_window_impute.sh) for every kit in the panel v3 manifest.
#   sbatch -A pi_matteop -p campus -a 1-99%60 -c 1 --mem=6G -t 4:00:00 -J iw_NAME cluster/indel_window_cohort.sh \
#          NAME CHROM POS REF ALT work_panel_v3/manifest.tsv 20
# Output: work_panel_v3/indel_imp/NAME/<pk>.tsv (CHROM POS REF ALT GT:GP); kits whose output exists are skipped.
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"
cd "$D"
NAME=$1; C=$2; P=$3; REF=$4; ALT=$5; MANIFEST=$6; PER=${7:-20}
R="$D/reference_panel/indel_windows/$NAME.ref.bcf"
O="$D/work_panel_v3/indel_imp/$NAME"; mkdir -p "$O" "$D/work_panel_v3/indel_imp/tmp"
start=$(( (SLURM_ARRAY_TASK_ID - 1) * PER + 1 )); end=$(( start + PER - 1 ))
ok=0; fail=0; skip=0
while IFS=$'\t' read -r pk barcode pub work; do
  [[ -n "$pk" ]] || continue
  out="$O/$pk.tsv"
  if [[ -s "$out" ]]; then skip=$((skip+1)); continue; fi
  [[ -s "$work/markdup.bam" ]] || { fail=$((fail+1)); continue; }
  if bash analysis/omia/indel_window_impute.sh "$R" "$C" "$P" "$REF" "$ALT" "$work/markdup.bam" "$D/work_panel_v3/indel_imp/tmp" > "$out.tmp" 2>/dev/null \
     && [[ -s "$out.tmp" ]]; then mv "$out.tmp" "$out"; ok=$((ok+1)); else rm -f "$out.tmp"; fail=$((fail+1)); echo "IW-FAIL $NAME $pk"; fi
done < <(sed -n "${start},${end}p" "$MANIFEST")
echo "IW-DONE $NAME task=$SLURM_ARRAY_TASK_ID ok=$ok fail=$fail skip=$skip"
