#!/bin/bash
# UGE array job: one task per sample-sheet row.
#
#   qsub -t 2-9 cluster/submit-array.sh sample_sheet.tsv
#
# For large arrays, cap how many run at once with -tc. Each task needs ~15-20GB
# of scratch under $D/scratch (~7TB free), so 40 concurrent is ~720GB — ample:
#
#   qsub -t 2-97 -tc 40 cluster/submit-array.sh sample_sheet.hoffman.tsv
#
# Task IDs are sample-sheet row numbers — row 1 is the header, so data rows
# start at 2. The pipeline already takes the row as its second argument, so no
# translation is needed.
#
# Resources. h_data is PER SLOT on Hoffman2, so 8 x 6G = 48G total.
#
# Stage 15 sets the memory floor, not alignment: bowtie2 loads the MetaPhlAn
# index — forward AND mirror (.1/.2 plus .rev.1/.rev.2), ~28.5GB total — into
# RAM, and a 32GB job died with "Out of memory allocating the ebwt[] array".
# Stage 7 peaks near 1.5GB and Stage 4 holds the ~10GB bwa-mem2 index, so
# neither is the binding constraint. 48G leaves ~19GB of headroom over the
# index; do not trim below ~40G without testing Stage 15 specifically.
#
# Confirm on a real run and trim if there is headroom:
#     qacct -j <job_id> | grep maxvmem
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

# Resolve the pipeline from the SUBMIT directory, not from $0. UGE spools the
# job script into /work/UGE/.../job_scripts/, so $0 points there and
# "$(dirname $0)/.." lands in the scheduler's spool, not the repo.
# -cwd means the job starts in the submit directory; SGE_O_WORKDIR records it.
PIPELINE_ROOT="${SGE_O_WORKDIR:-$PWD}"
[[ -f "$PIPELINE_ROOT/run_dog_pipeline.sh" ]] \
  || { echo "ERROR: run_dog_pipeline.sh not found in $PIPELINE_ROOT — submit from the repo root"; exit 1; }

bash "$PIPELINE_ROOT/run_dog_pipeline.sh" "$SHEET" "$ROW" 1

echo "=== task $ROW finished ==="
echo "Results are staged, not published — compute nodes have no internet."
echo "From a login node, run:  bash cluster/publish-pending.sh"
