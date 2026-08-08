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

SHEET="${1:?usage: qsub -t 2-N cluster/submit-array.sh <sample_sheet.tsv> [from_stage]}"
# Optional resume point. Defaults to 1 (full run). Use it to pick a sample back
# up without repeating alignment, e.g. after a stage-15 failure:
#     qsub -t 3-3 cluster/submit-array.sh sample_sheet.hoffman.tsv 5
# Resume only works when the intermediates still exist, so not for runs that
# used local scratch — those are discarded with the job.
FROM_STAGE="${2:-1}"
ROW="${SGE_TASK_ID:?this script must be submitted as an array job (-t)}"

export DOGS_SITE=hoffman

echo "=== task $ROW of $JOB_ID on $(hostname) ==="
echo "    slots=${NSLOTS:-?}  sheet=$SHEET  from_stage=$FROM_STAGE"

# Resolve the pipeline from the SUBMIT directory, not from $0. UGE spools the
# job script into /work/UGE/.../job_scripts/, so $0 points there and
# "$(dirname $0)/.." lands in the scheduler's spool, not the repo.
# -cwd means the job starts in the submit directory; SGE_O_WORKDIR records it.
PIPELINE_ROOT="${SGE_O_WORKDIR:-$PWD}"
[[ -f "$PIPELINE_ROOT/run_dog_pipeline.sh" ]] \
  || { echo "ERROR: run_dog_pipeline.sh not found in $PIPELINE_ROOT — submit from the repo root"; exit 1; }

# Run the script from NODE-LOCAL disk, not from /u/project.
#
# Bash does not read a script into memory — it reads incrementally, keeping a
# file offset, and executes as it goes. Stage 3+4 blocks inside a single command
# for ~35 minutes without touching that fd; on a shared parallel filesystem a
# read resumed after that gap can come back short. Bash takes a short read as
# end-of-script and stops: EXIT trap fires, summary prints, exit status 0, and
# every remaining stage is silently skipped. That is exactly how COSMO3 (job
# 14262939) "succeeded" after Stage 4 with no error on stdout or stderr.
#
# Copying to node-local disk first costs ~100ms and removes the whole class of
# failure. Only the script and site/ move; $D still points at project storage,
# so all data and outputs are unaffected.
STAGE_DIR="${TMPDIR:-/tmp}/pipeline.${JOB_ID:-manual}.$ROW"
mkdir -p "$STAGE_DIR"
trap 'rm -rf "$STAGE_DIR"' EXIT
cp "$PIPELINE_ROOT/run_dog_pipeline.sh" "$STAGE_DIR/"
cp -r "$PIPELINE_ROOT/site" "$STAGE_DIR/"
echo "    running from node-local $STAGE_DIR"

# Column 5 of the sheet is the sample's work dir, where the pipeline writes its
# pipeline.done marker (a local-scratch run copies the marker back there too).
# Clear any marker from an earlier run FIRST — a stale one would let a truncated
# rerun pass the completion check below.
SAMPLE_OUT=$(awk -F'\t' -v r="$ROW" 'NR==r {print $5}' "$SHEET")
if [[ -n "$SAMPLE_OUT" ]]; then rm -f "$SAMPLE_OUT/pipeline.done"; fi

bash "$STAGE_DIR/run_dog_pipeline.sh" "$SHEET" "$ROW" "$FROM_STAGE"

# The pipeline writes pipeline.done as its last act. Without this check a
# truncated run reports exit 0 and the scheduler calls it a success — across 96
# tasks that is silent, undetectable data loss.
if [[ -n "$SAMPLE_OUT" && ! -e "$SAMPLE_OUT/pipeline.done" ]]; then
  echo "ERROR: task $ROW ended without $SAMPLE_OUT/pipeline.done — the run did NOT" >&2
  echo "       reach the end of the pipeline. Check which stages appear in"        >&2
  echo "       $SAMPLE_OUT/pipeline.log before treating this sample as complete."  >&2
  exit 1
fi

echo "=== task $ROW finished ==="
echo "Results are staged, not published — compute nodes have no internet."
echo "From a login node, run:  bash cluster/publish-pending.sh"
