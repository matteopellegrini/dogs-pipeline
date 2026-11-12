#!/bin/bash
# Panel v3, step 1: genotype the frozen small-variant sites (known_variants/panel_v3/panel_v3_small.vcf.gz,
# 524 SNVs/indels) in every kit that still has its full alignment, and pull the existing GLIMPSE2 genotypes at
# the v3 sites that are already in the Dog10K panel. Nothing is published; output feeds the cohort-level
# allele-frequency / LD-aware calling step (analysis/omia/panel_v3_cohort.py).
#
#   per kit -> work_panel_v3/kits/<pk>.reads.bcf  (AD, DP, PL from reads; alleles forced to the catalogue's)
#              work_panel_v3/kits/<pk>.imp.bcf    (GT, DS, GP at v3 sites present in the imputed BCF)
#
#   sbatch -A pi_matteop -p campus -a 1-99%60 -c 1 --mem=4G -t 2:00:00 -J pv3geno --export=ALL,KITS=$PWD/work_panel_v3/kits2 \
#          cluster/panel_v3_genotype.sh work_panel_v3/manifest.tsv 20
# Manifest rows: pk, barcode, pub dir, work dir (same columns as the rerun lists). Task i does rows
# (i-1)*N+1 .. i*N. Kits whose two outputs already exist are skipped, so a resubmit only redoes failures.
# PV3_WORK points the site inputs and outputs at another work dir (e.g. sites added after the freeze).
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"
cd "$D"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PATH="$D/envs/genomics/bin:$PATH"
MANIFEST="${1:?manifest}"; PER="${2:-20}"
OUT="${PV3_WORK:-$D/work_panel_v3}"; KITS="${KITS:-$OUT/kits}"; mkdir -p "$KITS"
FASTA="$D/canFam4_plusY.fa"
SITES="$OUT/v3_small.regions.tsv"        # CHROM POS (for mpileup -R)
ALLELES="$OUT/v3_small.alleles.tsv.gz"   # CHROM POS REF,ALT (for call -C alleles -T)
IMPSITES="$OUT/d10k_match.pos"           # v3 sites with a matching Dog10K record
start=$(( (SLURM_ARRAY_TASK_ID - 1) * PER + 1 )); end=$(( start + PER - 1 ))
ok=0; fail=0; skip=0
while IFS=$'\t' read -r pk barcode pub work; do
  [[ -n "$pk" ]] || continue
  r="$KITS/$pk.reads.bcf"; i="$KITS/$pk.imp.bcf"
  if [[ -s "$r" && -s "$i" ]]; then skip=$((skip+1)); continue; fi
  bam="$work/markdup.bam"
  if [[ ! -s "$bam" ]]; then echo "PV3-NOBAM $pk"; fail=$((fail+1)); continue; fi
  tmp=$(mktemp "$KITS/.$pk.XXXXXX")
  # -B off (no BAQ) keeps single low-coverage reads usable. -m 1: one gapped read is enough to consider an indel
  # (the default 2 hides every 1x carrier). call -A keeps the catalogue ALT on hom-ref records (otherwise ALT=".",
  # which merges as a separate row); norm puts indels in the catalogue's left-normalised form.
  if bcftools mpileup -f "$FASTA" -R "$SITES" -q 20 -Q 20 -B -m 1 -d 1000 -a FORMAT/AD,FORMAT/DP -Ou "$bam" 2>/dev/null \
     | bcftools call -m -A -C alleles -T "$ALLELES" -Ou 2>/dev/null \
     | bcftools norm -f "$FASTA" -Ou 2>/dev/null \
     | bcftools reheader -s <(echo "$pk") 2>/dev/null \
     | bcftools view -Ob -o "$tmp" 2>/dev/null && [[ -s "$tmp" ]]; then
    mv "$tmp" "$r"; bcftools index -f "$r"
  else
    rm -f "$tmp"; echo "PV3-FAIL reads $pk"; fail=$((fail+1)); continue
  fi
  imp=$(ls "$work"/glimpse2/*_imputed_dog10k.bcf 2>/dev/null | head -1)
  if [[ -n "$imp" ]] && bcftools view -T "$IMPSITES" "$imp" 2>/dev/null \
       | bcftools reheader -s <(echo "$pk") 2>/dev/null | bcftools view -Ob -o "$i" 2>/dev/null; then
    bcftools index -f "$i"
  else
    echo "PV3-FAIL imputed $pk"; fail=$((fail+1)); continue
  fi
  ok=$((ok+1))
done < <(sed -n "${start},${end}p" "$MANIFEST")
echo "PV3-DONE task=$SLURM_ARRAY_TASK_ID rows=$start-$end ok=$ok fail=$fail skip=$skip"
