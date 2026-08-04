#!/bin/bash
# UGE array job: one task per sample-sheet row.
#
#   qsub -t 2-9 cluster/submit-array.sh sample_sheet.tsv
#
# Task IDs are sample-sheet row numbers — row 1 is the header, so data rows
# start at 2. The pipeline already takes the row as its second argument, so no
# translation is needed.
#
# Resources: peak memory observed on a full run is ~40GB. On Hoffman2 h_data is
# PER SLOT, so 8 slots x 6G = 48G total. Runtime is ~1h15m on a Mac; 8h gives
# generous headroom for slower or contended nodes.
#
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=6G,h_rt=8:00:00
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
