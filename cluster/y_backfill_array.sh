#!/bin/bash
#SBATCH -o logs/%x.%A.%a.log
# Slurm array over scratch/ybackfill_list.tsv rows (name, barcode, pub dir, work dir, legacy BAM, autosomal depth,
# predicted sex): runs cluster/y_backfill.sh into <work>/y_backfill and installs y_result.json in the pub dir.
# The pub dir's previous y_result.json (or an "absent" marker when there was none) is kept once as
# y_result.pre-<TAG>.json, so FILES=y_result.json cluster/republish_files.sh <TAG> <list> picks the kit up.
#   sbatch -A pi_matteop -p campus -a 1-N%60 -c 8 --mem=48G -t 4:00:00 -J ybf --export=ALL,TAG=ybf \
#          cluster/y_backfill_array.sh scratch/ybackfill_list.tsv
set -uo pipefail
D=$SLURM_SUBMIT_DIR; cd "$D"; TAG=${TAG:-ybf}
IFS=$'\t' read -r name bc pub work bam auto sex < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")
[ -n "$name" ] || exit 0
OUTD=$work/y_backfill
if D=$D bash cluster/y_backfill.sh "$bam" "$OUTD" "$auto" "$sex"; then
  if [ ! -f "$pub/y_result.pre-$TAG.json" ]; then
    if [ -f "$pub/y_result.json" ]; then cp -p "$pub/y_result.json" "$pub/y_result.pre-$TAG.json"
    else echo '{"status": "absent_before_y_backfill"}' > "$pub/y_result.pre-$TAG.json"; fi
    touch -d '1 minute ago' "$pub/y_result.pre-$TAG.json"
  fi
  cp "$OUTD/y_result.json" "$pub/y_result.json.tmp" && mv "$pub/y_result.json.tmp" "$pub/y_result.json"
  echo "YBF-DONE $name"
else
  echo "YBF-FAIL $name"; exit 1
fi
