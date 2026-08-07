#!/bin/bash
# UGE array job: one task per sample-sheet row.
#
#   qsub -t 2-9 cluster/submit-array.sh sample_sheet.tsv
#
# Task IDs are sample-sheet row numbers — row 1 is the header, so data rows
# start at 2. The pipeline already takes the row as its second argument, so no
# translation is needed.
#
# Resources. h_data is PER SLOT on Hoffman2, so 8 x 4G = 32G total.
#
# The 40GB figure from Mac runs was misleading: that machine ran 8 GLIMPSE jobs
# AND held the bwa-mem2 index simultaneously. A real Hoffman2 run peaked at
# 1.5GB during Stage 7. Alignment (Stage 4) is the actual high-water mark, hence
# 4G/slot rather than something much smaller — confirm with:
#     qacct -j <job_id> | grep maxvmem
# and trim further if it comes in low. Over-requesting only costs queue time.
#
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=4G,h_rt=8:00:00
#$ -pe shared 8

set -euo pipefail

SHEET="${1:?usage: qsub -t 2-N cluster/submit-array.sh <sample_sheet.tsv>}"
ROW="${SGE_TASK_ID:?this script must be submitted as an array job (-t)}"

export DOGS_SITE=hoffman

echo "=== task $ROW of $JOB_ID on $(hostname) ==="
echo "    slots=${NSLOTS:-?}  sheet=$SHEET"

bash "$(dirname "$0")/../run_dog_pipeline.sh" "$SHEET" "$ROW" 1

echo "=== task $ROW finished ==="
echo "Results are staged, not published — compute nodes have no internet."
echo "From a login node, run:  bash cluster/publish-pending.sh"
