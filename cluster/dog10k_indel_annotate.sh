#!/bin/bash
# Annotate the Dog10K indel call set (dog10k_indels/AutoAndXPAR.nonSNPs.filter.vcf.gz, 1,987 samples, GATK hard
# filters per callset.README.txt) like cluster/dog10k_annotate.sh does the SNV panel. Per chromosome: PASS records,
# split multiallelics and left-normalise on canFam4, allele frequency over the 1,987 dogs, genotypes dropped, snpEff
# -canon (ROS_Cfam_1.115); then a TSV of the HIGH-impact and in-frame records.
#   sbatch -A pi_matteop -p campus -a 1-39 -c 2 --mem=12G -t 6:00:00 -J d10kiann cluster/dog10k_indel_annotate.sh
# Output: work_dog10k_indel/<chrom>.ann.vcf.gz, <chrom>.damaging.tsv (CHROM POS REF ALT AF AC_Hom GENE EFFECT IMPACT HGVS_P)
#SBATCH -o logs/%x.%A.%a.log
set -euo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"
cd "$D"
export OPENBLAS_NUM_THREADS=1
E="$D/envs/genomics"
JAVA="$(find "$E/lib/jvm/bin" -name java 2>/dev/null | head -1 || true)"; [[ -n "$JAVA" ]] || JAVA="$E/bin/java"
export PATH="$E/bin:$(dirname "$JAVA"):$PATH"
I="$D/dog10k_indels"; O="$D/work_dog10k_indel"; mkdir -p "$O"
chroms=( $(for i in $(seq 1 38); do echo chr$i; done) chrX )
c=${chroms[$((SLURM_ARRAY_TASK_ID - 1))]}
# chrX: the pseudoautosomal stretch is in the autosome file, the rest in its own file
if [[ "$c" == chrX ]]; then
  bcftools concat -a -r chrX -Ou "$I/AutoAndXPAR.nonSNPs.filter.vcf.gz" "$I/chrX.NONPAR.nonSNPs.filter.vcf.gz"
else
  bcftools view -r "$c" -Ou "$I/AutoAndXPAR.nonSNPs.filter.vcf.gz"
fi \
  | bcftools view -f PASS -Ou \
  | bcftools norm -m -any -f "$D/canFam4_plusY.fa" -Ou 2>/dev/null \
  | bcftools +fill-tags -Ou -- -t AC,AN,AF,AC_Hom \
  | bcftools view -G -Ov \
  | snpEff ann -canon -noStats -noLog ROS_Cfam_1.115 2> "$O/$c.snpeff.log" \
  | bgzip -c > "$O/$c.ann.vcf.gz.tmp"
mv "$O/$c.ann.vcf.gz.tmp" "$O/$c.ann.vcf.gz"
bcftools index -t "$O/$c.ann.vcf.gz"
bcftools query -i 'INFO/ANN~"|HIGH|" || INFO/ANN~"inframe_"' -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC_Hom\t%INFO/ANN\n' "$O/$c.ann.vcf.gz" \
  | awk -F'\t' 'BEGIN{OFS="\t"} { split($7, a, ","); split(a[1], f, "|"); print $1,$2,$3,$4,$5,$6,f[4],f[2],f[3],f[11] }' \
  > "$O/$c.damaging.tsv"
echo "D10KIANN-DONE $c $(bcftools index -n "$O/$c.ann.vcf.gz") indel records, $(wc -l < "$O/$c.damaging.tsv") high/in-frame"
