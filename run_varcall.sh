#!/usr/bin/env bash
set -euo pipefail

cd /Users/matteopellegrini/Downloads/dogs

echo "[$(date)] Waiting for markdup to finish..."
wait $(pgrep -f "samtools markdup" 2>/dev/null) 2>/dev/null || true

# Verify markdup.bam is complete
micromamba run -n genomics samtools quickcheck markdup.bam || {
  echo "[$(date)] ERROR: markdup.bam failed quickcheck, re-running..."
  micromamba run -n genomics samtools markdup -@ 8 --write-index sorted.bam markdup.bam
}

echo "[$(date)] markdup.bam OK. Indexing if needed..."
[ -f markdup.bam.bai ] || micromamba run -n genomics samtools index -@ 8 markdup.bam

echo "[$(date)] Running bcftools mpileup | call..."
micromamba run -n genomics bcftools mpileup \
  -f canFam4.fa \
  --threads 8 \
  -q 20 -Q 20 \
  -a FORMAT/DP,FORMAT/AD \
  markdup.bam | \
micromamba run -n genomics bcftools call \
  -mv \
  --ploidy 2 \
  -o variants_raw.vcf.gz \
  -O z \
  --threads 8

echo "[$(date)] Indexing VCF..."
micromamba run -n genomics bcftools index variants_raw.vcf.gz

echo "[$(date)] Filtering variants..."
micromamba run -n genomics bcftools filter \
  -s LOWQUAL \
  -e 'QUAL<20 || DP<10' \
  variants_raw.vcf.gz | \
micromamba run -n genomics bcftools view \
  -f PASS \
  -o variants_filtered.vcf.gz \
  -O z

micromamba run -n genomics bcftools index variants_filtered.vcf.gz

echo "[$(date)] Done! Summary:"
micromamba run -n genomics bcftools stats variants_filtered.vcf.gz | grep "^SN"
