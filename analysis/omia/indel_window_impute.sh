#!/bin/bash
# Impute one disease indel for one kit from its haplotype: GLIMPSE2_phase on the kit's BAM against the indel window
# reference (analysis/omia/indel_window_ref.sh). GLIMPSE2 takes read likelihoods at the window's SNVs and fills the
# indel in from the reference haplotypes. Prints CHROM POS REF ALT GT:GP for the indel (the cohort_imputed.tsv format).
#   analysis/omia/indel_window_impute.sh REF_BCF CHROM POS REF ALT BAM [TMP_DIR]
set -euo pipefail
R=$1; C=$2; P=$3; REF=$4; ALT=$5; BAM=$6; TMPD=${7:-${TMPDIR:-/tmp}}
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export PATH="$ROOT/envs/genomics/bin:$PATH" OPENBLAS_NUM_THREADS=1
# (no `| head`: with pipefail its early exit fails the script)
read -r lo hi < <(bcftools query -f '%POS\n' "$R" | awk 'NR==1{f=$1} {l=$1} END{print f, l}')
olo=$(( P - 250000 > lo ? P - 250000 : lo )); ohi=$(( P + 250000 < hi ? P + 250000 : hi ))
T=$(mktemp -d "$TMPD/indelimp.XXXX")
"$ROOT/envs/glimpse/bin/GLIMPSE2_phase" --bam-file "$BAM" --reference "$R" --fasta "$ROOT/canFam4_plusY.fa" \
  --input-region "$C:$lo-$hi" --output-region "$C:$olo-$ohi" --threads 1 --output "$T/out.bcf" > "$T/log" 2>&1 \
  || { tail -3 "$T/log" >&2; rm -rf "$T"; exit 1; }
bcftools query -i "POS=$P && REF=\"$REF\" && ALT=\"$ALT\"" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%GP]\n' "$T/out.bcf"
rm -rf "$T"
