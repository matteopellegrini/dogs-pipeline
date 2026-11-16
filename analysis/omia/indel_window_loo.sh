#!/bin/bash
# K-fold check of a window reference: hide each fold's indel genotype, impute it back from the fold's own window SNVs
# against the other folds, and report truth vs imputed per sample. Best case (deep SNVs); customers' GLIMPSE SNVs are noisier.
# Usage: indel_window_loo.sh REF_BCF CHROM POS OUT_TSV [K]
set -euo pipefail
R=$1 C=$2 P=$3 OUT=$4 K=${5:-5}
ROOT=/u/project/pellegrini/gkislik/dogs
export PATH="$ROOT/envs/genomics/bin:$PATH" OPENBLAS_NUM_THREADS=1
T=$(mktemp -d "$(dirname "$OUT")/.loo.XXXX")
bcftools query -i "POS=$P" -f '[%SAMPLE\t%GT\n]' "$R" > "$T/truth.tsv"
# carriers spread across folds first, then everyone else
awk -F'\t' '$2~/1/' "$T/truth.tsv" | cut -f1 > "$T/order"; awk -F'\t' '$2!~/1/' "$T/truth.tsv" | cut -f1 | shuf --random-source=<(yes) >> "$T/order"
awk -v k=$K '{print > "'"$T"'/fold" (NR-1)%k}' "$T/order"
: > "$T/all.tsv"
for f in "$T"/fold*; do
  bcftools view -S ^"$f" -Oz -o "$f.ref.vcf.gz" "$R"
  bcftools view -S "$f" -e "POS=$P" "$R" -Ou | bcftools +setGT -Ou -- -t a -n u | bcftools view -Oz -o "$f.tgt.vcf.gz"
  bcftools index -t "$f.tgt.vcf.gz"
  "$ROOT/envs/beagle/bin/beagle" -Xmx8g gt="$f.tgt.vcf.gz" ref="$f.ref.vcf.gz" out="$f.imp" impute=true gprobs=true nthreads=${BEAGLE_THREADS:-4} > "$f.log" 2>&1
  bcftools query -i "POS=$P" -f '[%SAMPLE\t%DS\t%GP\n]' "$f.imp.vcf.gz" >> "$T/all.tsv"
done
awk -F'\t' 'NR==FNR{t[$1]=$2; next} {print $1"\t"t[$1]"\t"$2"\t"$3}' "$T/truth.tsv" "$T/all.tsv" > "$OUT"
rm -rf "$T"
awk -F'\t' '{split($4,g,","); c=($2~/1/); p=1-g[1]; n[c]++; s[c]+=$3; if(p<=0.002) clr[c]++; if(p>=0.5) hi[c]++}
  END{for(c=0;c<=1;c++) printf "truth %s: n=%d meanDS=%.2f imputedClear(P<=0.2%%)=%d P>=0.5=%d\n", (c?"carrier":"non-carrier"), n[c], s[c]/n[c], clr[c], hi[c]}' "$OUT"
