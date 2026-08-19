#!/usr/bin/env bash
# Score every dog in a sample sheet with the Darwin's Ark size PRS.
#   bash analysis/darwins_ark/score_cohort.sh sample_sheet.hoffman.tsv > size_prs_cohort.tsv
set -uo pipefail
SHEET=$1
DIR=$(dirname "$0")
echo -e "sample\tprs\tmatched"
tail -n +2 "$SHEET" | while IFS=$'\t' read -r sid _ _ out wd _rest; do
  bcf="$wd/glimpse2/${out}_imputed_dog10k.bcf"
  [ -s "$bcf" ] || { echo -e "$sid\tNA\tno_bcf" >&2; continue; }
  r=$(python3 "$DIR/score_dog.py" "$DIR/wts_size_0.1.tsv.gz" "$bcf")
  echo -e "$sid\t$r"
done
