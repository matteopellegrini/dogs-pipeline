#!/bin/bash
# Annotate the Dog10K imputation panel sites once (the same snpEff run stage 10 does per dog on its imputed BCF,
# whose sites are exactly these). Per chromosome: panel allele frequency (AF over the reference dogs), genotypes
# dropped, snpEff -canon with the pipeline's database; then a TSV of the HIGH-impact and missense records.
#   sbatch -A pi_matteop -p campus -a 1-39 -c 2 --mem=12G -t 6:00:00 -J d10kann cluster/dog10k_annotate.sh
# then, without an array (after all 39 tasks), the panel file stage 10 copies ANN/LOF/NMD from:
#   sbatch -A pi_matteop -p campus -c 2 --mem=8G -t 2:00:00 -J d10kannfin --dependency=afterok:<array job> cluster/dog10k_annotate.sh
# Output: work_dog10k_ann/<chrom>.ann.vcf.gz, <chrom>.damaging.tsv (CHROM POS REF ALT AF GENE EFFECT IMPACT HGVS_P AC_HOM)
#SBATCH -o logs/%x.%A.%a.log
set -euo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"
cd "$D"
export OPENBLAS_NUM_THREADS=1
E="$D/envs/genomics"
JAVA="$(find "$E/lib/jvm/bin" -name java 2>/dev/null | head -1 || true)"; [[ -n "$JAVA" ]] || JAVA="$E/bin/java"
export PATH="$E/bin:$(dirname "$JAVA"):$PATH"
DB="${SNPEFF_DB:-canFam4_NCBI_AR106}"   # built on canFam4 by cluster/build_snpeff_canfam4.sh
P="$D/dog10k_panel/AutoAndXPAR.Dog10K.phased_plus_disease_rh.bcf"
O="${D10K_ANN_WORK:-$D/work_dog10k_ann}"; mkdir -p "$O"
chroms=( $(for i in $(seq 1 38); do echo chr$i; done) chrX )
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  F="${D10K_ANN_OUT:-$D/reference_panel/dog10k_ann/dog10k_sites.ann.bcf}"; mkdir -p "$(dirname "$F")"
  for c in "${chroms[@]}"; do echo "$O/$c.ann.vcf.gz"; done > "$O/ann.list"
  bcftools concat -f "$O/ann.list" -Ou \
    | bcftools annotate -x INFO/AC,INFO/AF,INFO/NS,INFO/AN,INFO/AC_Hom -Ob -o "$F.tmp" --threads 2
  mv "$F.tmp" "$F"
  bcftools index -f "$F"
  echo "D10KANN-FINAL $(bcftools index -n "$F") sites -> $F ($DB)"
  exit 0
fi
c=${chroms[$((SLURM_ARRAY_TASK_ID - 1))]}
bcftools view -r "$c" -Ou "$P" \
  | bcftools +fill-tags -Ou -- -t AC,AN,AF,AC_Hom \
  | bcftools view -G -Ov \
  | snpEff ann -canon -noStats -noLog "$DB" 2> "$O/$c.snpeff.log" \
  | bgzip -c > "$O/$c.ann.vcf.gz.tmp"
mv "$O/$c.ann.vcf.gz.tmp" "$O/$c.ann.vcf.gz"
bcftools index -t "$O/$c.ann.vcf.gz"
# first ANN entry per allele is the canonical-transcript effect with the highest impact
bcftools query -i 'INFO/ANN~"|HIGH|" || INFO/ANN~"missense_variant"' -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/ANN\t%INFO/AC_Hom\n' "$O/$c.ann.vcf.gz" \
  | awk -F'\t' 'BEGIN{OFS="\t"} { split($6, a, ","); split(a[1], f, "|"); print $1,$2,$3,$4,$5,f[4],f[2],f[3],f[11],$7 }' \
  > "$O/$c.damaging.tsv"
echo "D10KANN-DONE $c $(bcftools index -n "$O/$c.ann.vcf.gz") sites, $(wc -l < "$O/$c.damaging.tsv") high/missense"
