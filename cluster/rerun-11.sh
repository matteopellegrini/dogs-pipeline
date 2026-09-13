#!/bin/bash
# Rerun stage 11 (PRS / weight) ONLY on an already-processed sample. Needs the
# imputed BCF in the sample's work dir and breed_result.json + coverage_1mb.json
# in its pub dir. Nothing is uploaded (PUBLISH_RESULTS=0).
#   qsub cluster/rerun-11.sh <sheet> <row>
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.log
#$ -l h_data=4G,h_rt=1:00:00
#$ -pe shared 2
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PIPELINE_DIR="${SGE_O_WORKDIR:?}"
SHEET="${1:?usage: qsub cluster/rerun-11.sh <sheet> <row>}"
ROW="${2:?row required}"
cd "$PIPELINE_DIR"
RDP_SNAPSHOT="$PIPELINE_DIR/.rdp.$JOB_ID.11.sh"
cp run_dog_pipeline.sh "$RDP_SNAPSHOT"
trap "rm -f \"$RDP_SNAPSHOT\"" EXIT
sample=$(awk -F'\t' -v r="$ROW" 'NR==r{print $4}' "$SHEET")
TO_STAGE=11 PUBLISH_RESULTS=0 bash "$RDP_SNAPSHOT" "$SHEET" "$ROW" 11 \
  || { echo "ERROR: stage 11 failed for $sample"; exit 1; }
echo "RERUN11-DONE $sample"
