#!/bin/bash
# Rerun stage 13 (coat color) ONLY on already-processed samples — no
# realignment, no stage 8. Used for the 2026-09-10 B-locus rule change (two
# distinct heterozygous TYRP1 brown alleles now call b/b in trans; the imputed
# phase is no longer used to call a B/b-cis carrier). Needs the imputed BCF and
# reads (kept markdup.bam or the permanent sites.bam) in the sample's work dir.
# ~1 min per dog.
#
#   qsub -t 1-N -tc 20 cluster/rerun-13.sh <sheet> <rows-file>
#
# Results land in the sheet's pub_dir (results_prosper/<kit>/coat_color.json);
# nothing is uploaded — republish from a login node afterwards.
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=4G,h_rt=1:00:00
#$ -pe shared 2
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PIPELINE_DIR="${SGE_O_WORKDIR:?}"
SHEET="${1:?usage: qsub -t 1-N cluster/rerun-13.sh <sheet> <rows-file>}"
ROWS="${2:?rows file required}"
cd "$PIPELINE_DIR"
# Run a SNAPSHOT of the pipeline script (see cluster/reconcile-8-13.sh): bash
# reads scripts incrementally, so an edit or pull mid-task poisons running
# invocations. The snapshot must live IN the repo root because the pipeline
# resolves PIPELINE_DIR from its own BASH_SOURCE dirname.
RDP_SNAPSHOT="$PIPELINE_DIR/.rdp.$JOB_ID.$SGE_TASK_ID.sh"
cp run_dog_pipeline.sh "$RDP_SNAPSHOT"
trap "rm -f \"$RDP_SNAPSHOT\"" EXIT
ROW=$(awk -v i="${SGE_TASK_ID:?}" 'NR==i{print; exit}' "$ROWS")
[[ -n "$ROW" ]] || { echo "no row for task $SGE_TASK_ID"; exit 1; }
sample=$(awk -F'\t' -v r="$ROW" 'NR==r{print $4}' "$SHEET")

TO_STAGE=13 PUBLISH_RESULTS=0 bash "$RDP_SNAPSHOT" "$SHEET" "$ROW" 13 \
  || { echo "ERROR: stage 13 failed for $sample"; exit 1; }
echo "RERUN13-DONE $sample"
