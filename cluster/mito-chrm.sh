#!/bin/bash
# Mitochondrial prototype: per-sample chrM variant calls without the full pipeline.
#
#   qsub -t 2-97 cluster/mito-chrm.sh sample_sheet.gen.tsv
#
# The cohort's whole-genome BAMs are deleted after each pipeline run, but chrM
# is 16.7kb, so remapping the raw FASTQs against chrM alone recovers the
# mitochondrial reads for a fraction of an alignment run. All reads are mapped
# to the chrM-only index; only mapped reads (i.e. mito-like) are kept.
#
# NUMT caveat: reads from nuclear mitochondrial insertions have nowhere else to
# map with a chrM-only index, so they land here too. True mito depth is 30-500x
# on this cohort while nuclear depth is 1-7x, so NUMTs contribute at most a few
# percent of reads; haploid consensus calling is insensitive to that. Validated
# against a whole-genome-BAM chrM slice for DOGS-Gen-2 before trusting cohort-wide.
#
# Output per sample under $D/mito/:
#   vcf/<s>.vcf.gz     haploid chrM calls (AD/DP annotated, for heteroplasmy)
#   fasta/<s>.fa       consensus mitogenome
#   depth/<s>.txt      mean chrM depth
#$ -cwd
#$ -j y
#$ -r y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=3G,h_rt=4:00:00
#$ -pe shared 4

set -euo pipefail

PIPELINE_DIR="${SGE_O_WORKDIR:?run via qsub from the repo root}"
SHEET="${1:?usage: qsub -t 2-N cluster/mito-chrm.sh <sample_sheet.tsv>}"
ROW="${SGE_TASK_ID:?}"

DOGS_SITE="${DOGS_SITE:-hoffman}"
source "$PIPELINE_DIR/site/${DOGS_SITE}.sh"
S="$D/envs/genomics/bin"
REF="$D/mito/ref/chrM.fa"

fastq_dir=$(awk -F'\t' -v r="$ROW" 'NR==r{print $2}' "$SHEET")
sample=$(awk -F'\t' -v r="$ROW" 'NR==r{print $4}' "$SHEET" | tr '[:upper:]' '[:lower:]')
[[ -n "$sample" && -n "$fastq_dir" ]] || { echo "ERROR: empty row $ROW"; exit 1; }
mkdir -p "$D/mito/vcf" "$D/mito/fasta" "$D/mito/depth"

r1=$(ls "$fastq_dir"/*_R1_*.f*q.gz 2>/dev/null | sort | tr '\n' ' ')
r2=$(ls "$fastq_dir"/*_R2_*.f*q.gz 2>/dev/null | sort | tr '\n' ' ')
[[ -n "$r1" && -n "$r2" ]] || { echo "ERROR: no _R1_/_R2_ fastqs in $fastq_dir"; exit 1; }

TMP="${TMPDIR:-/tmp}/mito_$sample"
mkdir -p "$TMP"

"$S/bwa-mem2" mem -t "${NSLOTS:-4}" "$REF" <(zcat $r1) <(zcat $r2) 2>"$TMP/bwa.log" \
  | "$S/samtools" view -b -F 4 -q 20 - \
  | "$S/samtools" sort -o "$TMP/$sample.chrM.bam" -
"$S/samtools" index "$TMP/$sample.chrM.bam"

"$S/bcftools" mpileup -f "$REF" -d 4000 -q 20 -Q 20 -a AD,DP "$TMP/$sample.chrM.bam" 2>/dev/null \
  | "$S/bcftools" call -mv --ploidy 1 -Oz -o "$D/mito/vcf/$sample.vcf.gz"
"$S/bcftools" index -f "$D/mito/vcf/$sample.vcf.gz"
"$S/bcftools" consensus -f "$REF" "$D/mito/vcf/$sample.vcf.gz" \
  | sed "1s/.*/>${sample}_chrM/" > "$D/mito/fasta/$sample.fa"
"$S/samtools" depth -a "$TMP/$sample.chrM.bam" \
  | awk -v s="$sample" '{t+=$3} END {printf "%s\t%.0f\n", s, t/16728}' > "$D/mito/depth/$sample.txt"

rm -rf "$TMP"
echo "DONE $sample: $("$S/bcftools" view -H "$D/mito/vcf/$sample.vcf.gz" | wc -l) variants"
