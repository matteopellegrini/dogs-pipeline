#!/usr/bin/env bash
# Score every dog in a sample sheet with the Darwin's Ark size PRS.
#
#   bash analysis/darwins_ark/score_cohort.sh <sample_sheet.tsv> [out.tsv]
#
# qsub-friendly: UGE spools the script to /work/UGE/.../job_scripts, so $0
# does NOT point at the repo — resolve everything from SGE_O_WORKDIR (the
# submit directory) exactly as cluster/submit-array.sh does.
#   qsub -cwd -j y -o logs/size_prs.log -l h_data=4G,h_rt=6:00:00 \
#        analysis/darwins_ark/score_cohort.sh sample_sheet.hoffman.tsv size_prs_cohort.tsv
#$ -S /bin/bash
set -uo pipefail
SHEET=$1
OUT=${2:-size_prs_cohort.tsv}
ROOT="${SGE_O_WORKDIR:-$(pwd)}"
DIR="$ROOT/analysis/darwins_ark"
[ -f "$DIR/score_dog.py" ] || DIR="$(cd "$(dirname "$0")" && pwd)"
[ -f "$DIR/score_dog.py" ] || { echo "ERROR: cannot locate score_dog.py (submit from repo root)" >&2; exit 1; }

D_DEFAULT="/u/project/pellegrini/$USER/dogs"
[ -d "$D_DEFAULT/envs/genomics/bin" ] && export PATH="$D_DEFAULT/envs/genomics/bin:$PATH"
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH" >&2; exit 1; }

NROWS=$(wc -l < "$SHEET")
echo "sheet $SHEET: $((NROWS-1)) data rows" >&2
echo -e "sample\tprs\tmatched" > "$OUT"
# awk per row: tab-separated fields survive empty columns, which bash `read`
# with IFS=tab silently collapses.
for row in $(seq 2 "$NROWS"); do
  sid=$(awk -F'\t' -v r="$row" 'NR==r{print $1}' "$SHEET")
  out=$(awk -F'\t' -v r="$row" 'NR==r{print $4}' "$SHEET")
  wd=$(awk  -F'\t' -v r="$row" 'NR==r{print $5}' "$SHEET")
  [ -n "$sid" ] || continue
  bcf="$wd/glimpse2/${out}_imputed_dog10k.bcf"
  if [ ! -s "$bcf" ]; then echo -e "$sid\tNA\tno_bcf" >> "$OUT"; echo "skip $sid: no $bcf" >&2; continue; fi
  r=$(python3 "$DIR/score_dog.py" "$DIR/wts_size_0.1.tsv.gz" "$bcf" < /dev/null) || { echo -e "$sid\tNA\tscore_fail" >> "$OUT"; continue; }
  echo -e "$sid\t$r" >> "$OUT"
  echo "done $sid: $r" >&2
done
echo "wrote $OUT" >&2
