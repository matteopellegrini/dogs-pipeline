#!/bin/bash
# Array form of rerun-11.sh: stage 11 (PRS / weight) only, over the rows listed
# in <rows-file> of <sheet>. Used 2026-09-13 to roll the rebuilt weight model
# (breed composition x DA size PRS x sex) onto every published report. Needs the
# imputed BCF in the work dir; nothing is uploaded (republish_w11.sh afterwards).
#   qsub -t 1-N -tc 40 cluster/rerun-11-array.sh <sheet> <rows-file>
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=4G,h_rt=2:00:00
#$ -pe shared 2
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PIPELINE_DIR="${SGE_O_WORKDIR:?}"
SHEET="${1:?usage: qsub -t 1-N cluster/rerun-11-array.sh <sheet> <rows-file>}"
ROWS="${2:?rows file required}"
cd "$PIPELINE_DIR"
RDP_SNAPSHOT="$PIPELINE_DIR/.rdp.$JOB_ID.$SGE_TASK_ID.sh"
cp run_dog_pipeline.sh "$RDP_SNAPSHOT"
trap "rm -f \"$RDP_SNAPSHOT\"" EXIT
ROW=$(awk -v i="${SGE_TASK_ID:?}" "NR==i{print; exit}" "$ROWS")
[[ -n "$ROW" ]] || { echo "no row for task $SGE_TASK_ID"; exit 1; }
sample=$(awk -F"\t" -v r="$ROW" "NR==r{print \$4}" "$SHEET")
TO_STAGE=11 PUBLISH_RESULTS=0 bash "$RDP_SNAPSHOT" "$SHEET" "$ROW" 11 \
  || { echo "ERROR: stage 11 failed for $sample"; exit 1; }
echo "RERUN11-DONE $sample"
