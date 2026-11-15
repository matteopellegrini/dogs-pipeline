#!/bin/bash
# Impute a disease indel into many kits at once from a window reference (indel_window_ref_reads.sh PUBLIC_ONLY=1):
# the kits' own GLIMPSE genotypes in the window (their stage-12 imputed BCF) are the target, Beagle 4.1 fills the
# indel from the reference haplotypes and reports its dosage and genotype probabilities.
#   analysis/omia/indel_window_target_impute.sh REF_BCF CHROM POS LIST OUT_TSV
# LIST: kit <tab> work_dir (the kit's <work_dir>/glimpse2/<kit>_imputed_dog10k.bcf). OUT_TSV: kit GT DS GP.
set -euo pipefail
R=$1; C=$2; P=$3; LIST=$4; OUT=$5
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export PATH="$ROOT/envs/genomics/bin:$PATH" OPENBLAS_NUM_THREADS=1
T=$(mktemp -d "$(dirname "$OUT")/.tgt.XXXX")
read -r lo hi < <(bcftools query -f '%POS\n' "$R" | awk 'NR==1{f=$1} {l=$1} END{print f, l}')
bcftools query -f '%CHROM\t%POS\n' "$R" > "$T/sites.tsv"
n=0
while IFS=$'\t' read -r kit w; do
  f="$w/glimpse2/${kit}_imputed_dog10k.bcf"
  [ -s "$f" ] || continue
  echo "$kit" > "$T/name"
  bcftools view -r "$C:$lo-$hi" -T "$T/sites.tsv" -m2 -M2 "$f" -Ou 2>/dev/null | bcftools annotate -x INFO,^FORMAT/GT -Ou \
    | bcftools reheader -s "$T/name" - -o "$T/k$n.bcf" 2>/dev/null && bcftools index -f "$T/k$n.bcf" && n=$((n+1))
done < "$LIST"
ls "$T"/k*.bcf | split -l 400 - "$T/batch."
for b in "$T"/batch.*; do bcftools merge -m none -l "$b" -Ob -o "$b.bcf" && bcftools index -f "$b.bcf"; done
ls "$T"/batch.*.bcf > "$T/batches"
bcftools merge -m none -l "$T/batches" -Oz -o "$T/target.vcf.gz"
bcftools index -t "$T/target.vcf.gz"
bcftools view -Oz -o "$T/ref.vcf.gz" "$R"
"$ROOT/envs/beagle/bin/beagle" -Xmx16g gt="$T/target.vcf.gz" ref="$T/ref.vcf.gz" out="$T/imp" impute=true gprobs=true \
  nthreads=${BEAGLE_THREADS:-4} > "$T/beagle.log" 2>&1
bcftools query -i "POS=$P" -f '[%SAMPLE\t%GT\t%DS\t%GP\n]' "$T/imp.vcf.gz" > "$OUT"
echo "$n kits imputed at $C:$P -> $OUT ($(awk '$3>=0.5' "$OUT" | wc -l) with dosage >= 0.5)"
rm -rf "$T"
