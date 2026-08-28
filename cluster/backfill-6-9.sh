#!/bin/bash
# Backfill staged cohort reports with the current stage-6 and stage-9 outputs
# (fixed karyotype tooltips + re-based QC + mito lineage + 231-breed call +
# relatives) without re-running alignment or imputation.
#
#   qsub -t 2-97 -tc 30 cluster/backfill-6-9.sh sample_sheet.gen.tsv
#
# Needs per-dog survivors in the work dir: coverage TSVs (stage 6 inputs),
# the imputed BCF (stage 9 input), and the raw FASTQs (stage 6d mito remap).
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=4G,h_rt=3:00:00
#$ -pe shared 4

set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

PIPELINE_DIR="${SGE_O_WORKDIR:?run via qsub from the repo root}"
SHEET="${1:?usage: qsub -t 2-N cluster/backfill-6-9.sh <sample_sheet.tsv>}"
ROW="${SGE_TASK_ID:?}"
cd "$PIPELINE_DIR"

sample=$(awk -F'\t' -v r="$ROW" 'NR==r{print $4}' "$SHEET")
[[ -n "$sample" ]] || { echo "ERROR: empty row $ROW"; exit 1; }

rc=0
TO_STAGE=6 PUBLISH_RESULTS=0 bash run_dog_pipeline.sh "$SHEET" "$ROW" 6 || rc=$?
if (( rc != 0 )); then echo "ERROR: stage 6 failed for $sample (rc=$rc)"; exit 1; fi
TO_STAGE=9 PUBLISH_RESULTS=0 bash run_dog_pipeline.sh "$SHEET" "$ROW" 9 || rc=$?
if (( rc != 0 )); then echo "ERROR: stage 9 failed for $sample (rc=$rc)"; exit 1; fi
echo "BACKFILL-DONE $sample"
