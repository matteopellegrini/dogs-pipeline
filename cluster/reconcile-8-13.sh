#!/bin/bash
# Reconciliation FAST LANE: rerun the fixed callers (stage 8 known variants,
# stage 13 coat/merle) for dogs whose reads are still available (sites.bam or
# kept markdup.bam) — no realignment. ~1 min per dog.
#
#   qsub -t 1-N -tc 40 cluster/reconcile-8-13.sh <sheet> <rows-file>
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=4G,h_rt=2:00:00
#$ -pe shared 4
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PIPELINE_DIR="${SGE_O_WORKDIR:?}"
SHEET="${1:?usage: qsub -t 1-N cluster/reconcile-8-13.sh <sheet> <rows-file>}"
ROWS="${2:?rows file required}"
cd "$PIPELINE_DIR"
ROW=$(awk -v i="${SGE_TASK_ID:?}" 'NR==i{print; exit}' "$ROWS")
[[ -n "$ROW" ]] || { echo "no row for task $SGE_TASK_ID"; exit 1; }
sample=$(awk -F'\t' -v r="$ROW" 'NR==r{print $4}' "$SHEET")

TO_STAGE=8  PUBLISH_RESULTS=0 bash run_dog_pipeline.sh "$SHEET" "$ROW" 8 \
  || { echo "ERROR: stage 8 failed for $sample"; exit 1; }
TO_STAGE=13 PUBLISH_RESULTS=0 bash run_dog_pipeline.sh "$SHEET" "$ROW" 13 \
  || { echo "ERROR: stage 13 failed for $sample"; exit 1; }
echo "RECONCILE-DONE $sample"
