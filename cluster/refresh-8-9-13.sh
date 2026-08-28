#!/bin/bash
# Refresh known-variants (graded confidence), relatives (enriched matches) and
# coat color (B-locus fix) on already-processed samples, without realignment.
#   qsub -t 2-97 -tc 30 cluster/refresh-8-9-13.sh sample_sheet.gen.tsv
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=4G,h_rt=2:00:00
#$ -pe shared 4
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PIPELINE_DIR="${SGE_O_WORKDIR:?}"
SHEET="${1:?usage: qsub -t 2-N cluster/refresh-8-9-13.sh <sheet>}"
ROW="${SGE_TASK_ID:?}"
cd "$PIPELINE_DIR"
sample=$(awk -F'\t' -v r="$ROW" 'NR==r{print $4}' "$SHEET")
rc=0
TO_STAGE=9  PUBLISH_RESULTS=0 bash run_dog_pipeline.sh "$SHEET" "$ROW" 8  || rc=$?
(( rc == 0 )) || { echo "ERROR: stages 8-9 failed for $sample (rc=$rc)"; exit 1; }
TO_STAGE=13 PUBLISH_RESULTS=0 bash run_dog_pipeline.sh "$SHEET" "$ROW" 13 || rc=$?
(( rc == 0 )) || { echo "ERROR: stage 13 failed for $sample (rc=$rc)"; exit 1; }
echo "REFRESH-DONE $sample"
