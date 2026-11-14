#!/bin/bash
# Paternal line (stage 9d) for a kit aligned before 2026-09-13 to canFam4.fa, which has no chrY.
# Re-aligns only the reads a chrY could have claimed (unmapped, MAPQ < 30, or on chrX; no duplicates or
# secondary/supplementary records) to canFam4_plusY.fa as single-end reads, then runs the stage-9d caller on that
# small BAM. Y reads have nowhere good to go on the female reference, so they sit in exactly these bins.
#   bash cluster/y_backfill.sh LEGACY_BAM OUT_DIR AUTO_DEPTH SEX     (SEX = male|female|unknown)
# Writes OUT_DIR/y_result.json (+ y_subset.bam = its chrY records, kept for audit). Slurm array wrapper: cluster/y_backfill_array.sh.
set -euo pipefail
BAM=${1:?legacy BAM}; OUTD=${2:?out dir}; AUTO=${3:-0}; SEX=${4:-unknown}
D=${D:-/u/project/pellegrini/gkislik/dogs}; ENV=$D/envs/genomics; NT=${SLURM_CPUS_PER_TASK:-8}
export LD_LIBRARY_PATH=$ENV/lib:${LD_LIBRARY_PATH:-} OPENBLAS_NUM_THREADS=1
mkdir -p "$OUTD"; tmp=$(mktemp -d "$OUTD/.ytmp.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
"$ENV/bin/samtools" view -u -F 0xD00 -e 'flag.unmap || mapq < 30 || rname == "chrX"' "$BAM" \
  | "$ENV/bin/samtools" fastq -n - > "$tmp/r.fq" 2>/dev/null
"$ENV/bin/bwa-mem2" mem -t "$NT" "$D/canFam4_plusY.fa" "$tmp/r.fq" 2>"$tmp/bwa.log" \
  | "$ENV/bin/samtools" sort -@ 2 -m 1G -o "$tmp/y.bam" -
"$ENV/bin/samtools" index "$tmp/y.bam"
# keep the chrY records only (all the caller reads); the rest of the subset is mostly oral-microbe reads
"$ENV/bin/samtools" view -b -o "$tmp/chrY.bam" "$tmp/y.bam" chrY && "$ENV/bin/samtools" index "$tmp/chrY.bam"
mv "$tmp/chrY.bam" "$OUTD/y_subset.bam"; mv "$tmp/chrY.bam.bai" "$OUTD/y_subset.bam.bai"
"$ENV/bin/python3" "$D/analysis/y_haplogroup/y_haplogroup.py" "$OUTD/y_subset.bam" "$D/reference_json/y_tree_oetjens2018.tsv" \
    "$OUTD/y_result.json" --auto-depth "$AUTO" --sex "$SEX"
echo "YBACKFILL-DONE $BAM $(grep -o '"status": "[a-z_]*"' "$OUTD/y_result.json")"
