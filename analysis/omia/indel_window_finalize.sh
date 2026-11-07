#!/bin/bash
# Last step of analysis/omia/indel_window_ref.sh on an existing Beagle output: contig lines + AC/AN, index, carrier
# check, temp dir removed. For builds whose script was edited on disk while they ran (bash reads scripts lazily).
#   analysis/omia/indel_window_finalize.sh NAME POS REF OUT_DIR
set -euo pipefail
NAME=$1; P=$2; REF=$3; O=$4
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export PATH="$ROOT/envs/genomics/bin:$PATH" OPENBLAS_NUM_THREADS=1
T=$(ls -d "$O/.$NAME".* | head -1)
[[ -s "$T/phased.vcf.gz" ]] || { echo "ERROR: no Beagle output in $T"; exit 1; }
bcftools reheader --fai "$ROOT/canFam4_plusY.fa.fai" "$T/phased.vcf.gz" | bcftools +fill-tags -Ob -o "$O/$NAME.ref.bcf" -- -t AC,AN
bcftools index -f "$O/$NAME.ref.bcf"
before=$(awk '{n+=gsub(/1/,"",$2)} END{print n+0}' "$O/$NAME.carriers.tsv")
after=$(bcftools query -i "POS=$P && REF=\"$REF\"" -f '[%GT\n]' "$O/$NAME.ref.bcf" | awk '{n+=gsub(/1/,"")} END{print n+0}')
echo "$NAME: $(bcftools view -H "$O/$NAME.ref.bcf" | wc -l) sites; carrier haplotypes $before before phasing, $after after"
rm -rf "$T"
