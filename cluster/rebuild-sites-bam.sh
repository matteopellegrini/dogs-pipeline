#!/bin/bash
# Recovery + hardening array: realign a dog whose markdup.bam was cleaned up,
# extract the permanent sites.bam (stage 3+4 now does this), rerun the
# read-dependent report stages (8, 13), then delete the big BAM again.
#
#   qsub -t 1-N -tc 15 cluster/rebuild-sites-bam.sh <sheet> <rows-file>
#
# <rows-file>: one sample-sheet row number per line; task i processes line i.
# Alignment (stage 4) holds the ~10GB bwa-mem2 index; 4 slots x 14G = 56G is
# the proven profile from the ProsperKits batch reshape.
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=14G,h_rt=6:00:00
#$ -pe shared 4
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PIPELINE_DIR="${SGE_O_WORKDIR:?}"
SHEET="${1:?usage: qsub -t 1-N cluster/rebuild-sites-bam.sh <sheet> <rows-file>}"
ROWS="${2:?rows file required}"
cd "$PIPELINE_DIR"
ROW=$(awk -v i="${SGE_TASK_ID:?}" 'NR==i{print; exit}' "$ROWS")
[[ -n "$ROW" ]] || { echo "no row for task $SGE_TASK_ID"; exit 1; }
sample=$(awk -F'\t' -v r="$ROW" 'NR==r{print $4}' "$SHEET")
outdir=$(awk -F'\t' -v r="$ROW" 'NR==r{print $5}' "$SHEET")

# Realign: stages 1-5 produce markdup.bam, (new) sites.bam, and the coverage
# window counts including the fixed 50kb grid for the unified CNV framework —
# the one chance to count windows before the big BAM is deleted again.
# USE_LOCAL_SCRATCH=0: the scratch profile copies back only a whitelist and
# DELETES markdup.bam with the scratch dir — stages 8/13 below then have no
# reads (this is also how the 92 dogs lost their BAMs originally). Write to
# the persistent work dir; the big BAM is removed by this script at the end.
USE_LOCAL_SCRATCH=0 TO_STAGE=5 PUBLISH_RESULTS=0 bash run_dog_pipeline.sh "$SHEET" "$ROW" 1 \
  || { echo "ERROR: stages 1-5 failed for $sample"; exit 1; }
[[ -f "$outdir/sites.bam" ]] || { echo "ERROR: sites.bam missing for $sample"; exit 1; }

# Read-dependent report stages against the fresh BAM.
TO_STAGE=8  PUBLISH_RESULTS=0 bash run_dog_pipeline.sh "$SHEET" "$ROW" 8 \
  || { echo "ERROR: stage 8 failed for $sample"; exit 1; }
TO_STAGE=13 PUBLISH_RESULTS=0 bash run_dog_pipeline.sh "$SHEET" "$ROW" 13 \
  || { echo "ERROR: stage 13 failed for $sample"; exit 1; }

# Big BAM goes away again; sites.bam (a few MB) stays forever.
rm -f "$outdir/markdup.bam" "$outdir/markdup.bam.csi" "$outdir/markdup.bam.bai"
echo "REBUILD-DONE $sample"
