#!/bin/bash
# Reference panel window with a disease indel phased in, so GLIMPSE2 can impute the indel from the haplotype it sits on
# (the Dog10K imputation panel has SNVs only; the indel genotypes come from the Dog10K non-SNP call set for the same
# 1,929 dogs). The window's panel SNVs and the indel's genotypes are merged and re-phased together with Beagle 4.1;
# the window only has to be consistent with itself, not with the genome-wide panel phase.
#   analysis/omia/indel_window_ref.sh NAME CHROM POS REF ALT OUT_DIR [FLANK_BP]
# Output: OUT_DIR/NAME.ref.bcf (+ .csi) and OUT_DIR/NAME.carriers.tsv (Dog10K carriers of the indel).
set -euo pipefail
NAME=$1; C=$2; P=$3; REF=$4; ALT=$5; O=$6; F=${7:-1000000}
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export PATH="$ROOT/envs/genomics/bin:$PATH" OPENBLAS_NUM_THREADS=1
PANEL="$ROOT/dog10k_panel/AutoAndXPAR.Dog10K.phased_plus_disease_rh.bcf"
INDELS="$ROOT/dog10k_indels/AutoAndXPAR.nonSNPs.filter.vcf.gz"
FASTA="$ROOT/canFam4_plusY.fa"
mkdir -p "$O"; T=$(mktemp -d "$O/.$NAME.XXXX")
bcftools query -l "$PANEL" > "$T/samples.txt"
lo=$(( P > F ? P - F : 1 )); hi=$(( P + F ))
bcftools view -r "$C:$lo-$hi" -S "$T/samples.txt" -m2 -M2 -v snps -Oz -o "$T/snv.vcf.gz" "$PANEL"
bcftools index -t "$T/snv.vcf.gz"
# the indel: left-normalised, split from any multiallelic record, exactly the catalogue allele; same dogs, same order
bcftools view -r "$C:$P" -S "$T/samples.txt" "$INDELS" -Ou 2>/dev/null \
  | bcftools norm -m -any -f "$FASTA" -Ou 2>/dev/null \
  | bcftools view -i "REF=\"$REF\" && ALT=\"$ALT\"" -Ou \
  | bcftools annotate -x INFO,^FORMAT/GT -Oz -o "$T/indel.vcf.gz"
bcftools index -t "$T/indel.vcf.gz"
n=$(bcftools view -H "$T/indel.vcf.gz" | wc -l)
[[ $n -eq 1 ]] || { echo "ERROR: $n records for $C:$P $REF>$ALT in the Dog10K indel calls"; exit 1; }
bcftools query -f '[%SAMPLE\t%GT\n]' "$T/indel.vcf.gz" | awk '$2 ~ /1/' > "$O/$NAME.carriers.tsv"
bcftools concat -a -Oz -o "$T/merged.vcf.gz" "$T/snv.vcf.gz" "$T/indel.vcf.gz"
bcftools index -t "$T/merged.vcf.gz"
"$ROOT/envs/beagle/bin/beagle" -Xmx8g gt="$T/merged.vcf.gz" out="$T/phased" nthreads=${BEAGLE_THREADS:-4} > "$T/beagle.log" 2>&1
# Beagle writes no ##contig lines, which BCF needs
# and GLIMPSE2 needs AC/AN on a reference panel
bcftools reheader --fai "$FASTA.fai" "$T/phased.vcf.gz" | bcftools +fill-tags -Ob -o "$O/$NAME.ref.bcf" -- -t AC,AN
bcftools index -f "$O/$NAME.ref.bcf"
# carrier haplotypes before/after phasing must agree in number
before=$(awk '{n+=gsub(/1/,"",$2)} END{print n+0}' "$O/$NAME.carriers.tsv")
after=$(bcftools query -i "POS=$P && REF=\"$REF\"" -f '[%GT\n]' "$O/$NAME.ref.bcf" | awk '{n+=gsub(/1/,"")} END{print n+0}')
echo "$NAME: $(bcftools view -H "$O/$NAME.ref.bcf" | wc -l) sites in $C:$lo-$hi; carrier haplotypes $before before phasing, $after after"
rm -rf "$T"
