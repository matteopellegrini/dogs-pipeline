#!/usr/bin/env bash
# Score every dog in a sample sheet with the Darwin's Ark size PRS.
#
#   bash analysis/darwins_ark/score_cohort.sh <sample_sheet.tsv> [out.tsv]
#
# qsub-friendly: submit the script itself, no -b y / bash -c quoting needed:
#   qsub -cwd -j y -o logs/size_prs.log -l h_data=4G,h_rt=6:00:00 \
#        analysis/darwins_ark/score_cohort.sh sample_sheet.hoffman.tsv size_prs_cohort.tsv
#$ -S /bin/bash
set -uo pipefail
SHEET=$1
OUT=${2:-size_prs_cohort.tsv}
DIR=$(cd "$(dirname "$0")" && pwd)

# Hoffman keeps bcftools/python3 in the pipeline's genomics env, not on the
# default PATH of a compute node; the Mac has them ambient.
D_DEFAULT="/u/project/pellegrini/$USER/dogs"
[ -d "$D_DEFAULT/envs/genomics/bin" ] && export PATH="$D_DEFAULT/envs/genomics/bin:$PATH"
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH" >&2; exit 1; }

echo -e "sample\tprs\tmatched" > "$OUT"
tail -n +2 "$SHEET" | while IFS=$'\t' read -r sid _ _ out wd _rest; do
  bcf="$wd/glimpse2/${out}_imputed_dog10k.bcf"
  [ -s "$bcf" ] || { echo -e "$sid\tNA\tno_bcf" >> "$OUT"; continue; }
  r=$(python3 "$DIR/score_dog.py" "$DIR/wts_size_0.1.tsv.gz" "$bcf")
  echo -e "$sid\t$r" >> "$OUT"
  echo "done $sid: $r" >&2
done
echo "wrote $OUT" >&2
