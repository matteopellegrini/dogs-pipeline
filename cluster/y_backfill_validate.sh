#!/bin/bash
#SBATCH -o logs/%x.%A.%a.log
# Validation for cluster/y_backfill.sh: take a dog aligned WITH chrY, subsample, re-align the subsample to the old
# canFam4.fa (no chrY) to recreate a pre-2026-09-13 kit, back-fill it, and call the same subsample directly on its
# chrY alignment. Rows of $1: name, plusY BAM, fraction, sex, autosomal depth of the subsample.
set -euo pipefail
D=$SLURM_SUBMIT_DIR; cd "$D"; ENV=$D/envs/genomics; NT=${SLURM_CPUS_PER_TASK:-8}
export LD_LIBRARY_PATH=$ENV/lib:${LD_LIBRARY_PATH:-} OPENBLAS_NUM_THREADS=1
read -r name bam frac sex auto < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")
W=$D/scratch/y_backfill_val/$name; mkdir -p "$W"
"$ENV/bin/samtools" view -b -s "11${frac#0}" -F 0x900 "$bam" -o "$W/sub_plusY.bam"
"$ENV/bin/samtools" index "$W/sub_plusY.bam"
"$ENV/bin/samtools" collate -u -O "$W/sub_plusY.bam" "$W/col" | "$ENV/bin/samtools" fastq -n -1 "$W/r1.fq" -2 "$W/r2.fq" -0 /dev/null -s /dev/null - 2>/dev/null
"$ENV/bin/bwa-mem2" mem -t "$NT" "$D/canFam4_idx" "$W/r1.fq" "$W/r2.fq" 2>/dev/null \
  | "$ENV/bin/samtools" fixmate -m -u - - | "$ENV/bin/samtools" sort -u -@ 2 -m 1G - \
  | "$ENV/bin/samtools" markdup - "$W/legacy.bam"
"$ENV/bin/samtools" index "$W/legacy.bam"; rm -f "$W/r1.fq" "$W/r2.fq"
"$ENV/bin/python3" analysis/y_haplogroup/y_haplogroup.py "$W/sub_plusY.bam" reference_json/y_tree_oetjens2018.tsv "$W/direct.json" --auto-depth "$auto" --sex "$sex"
bash cluster/y_backfill.sh "$W/legacy.bam" "$W/backfill" "$auto" "$sex"
echo "YVAL-DONE $name auto=$auto"
