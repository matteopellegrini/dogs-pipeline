#!/bin/bash
# COSMO-SUPER truth set: direct read counts at every Dog10K imputation site from Cosmo's combined ~18x BAM,
# split by read source (read-name prefix: ML15 = Admera, V350 = UCLA Zymo+Genotek lane, LH00427 = DOGS-Gen-50),
# then score each Cosmo run's imputed genotypes against ALL and leave-one-source-out truth (score_chrom.py).
#   qsub -t 1-38 cluster/cosmo-truth.sh        (one autosome per task; 1 slot schedules quickly)
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=4G,h_rt=12:00:00
set -uo pipefail
D="${SGE_O_WORKDIR:?run via qsub from the repo root}"; cd "$D"
export PATH=$D/envs/genomics/bin:$PATH
chr=chr$SGE_TASK_ID; OUT=analysis/cosmo_truth/out; mkdir -p "$OUT"
BAM=work/COSMO-SUPER/analysis/markdup.bam; G=glimpse2
SUPER=work/COSMO-SUPER/analysis/$G/COSMO-SUPER_imputed_dog10k.bcf
[ -s "$OUT/$chr.counts.tsv" ] && { echo "TRUTH-SKIP $chr"; exit 0; }
echo "[$(date +%T)] $chr on $(hostname)"
bcftools query -r "$chr" -f '%CHROM\t%POS\n' "$SUPER" 2>/dev/null | uniq > "$OUT/$chr.sites"
if [ ! -s "$OUT/$chr.pileup.tsv.gz" ]; then
  samtools view -h "$BAM" "$chr" \
  | awk 'BEGIN{OFS="\t"} /^@RG/{next}
         /^@/{print; next}
         !h{print "@RG\tID:ADMERA\tSM:ADMERA"; print "@RG\tID:UCLA\tSM:UCLA"; print "@RG\tID:GEN50\tSM:GEN50"; print "@RG\tID:OTHER\tSM:OTHER"; h=1}
         {g = $1 ~ /^ML15/ ? "ADMERA" : $1 ~ /^V350/ ? "UCLA" : $1 ~ /^LH00427/ ? "GEN50" : "OTHER"
          for (i = 12; i <= NF; i++) if ($i ~ /^RG:Z:/) $i = "RG:Z:" g
          print}' \
  | bcftools mpileup -f canFam4_plusY.fa -T "$OUT/$chr.sites" -q 30 -Q 20 -a AD -d 1000 -Ou - 2>/dev/null \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' | gzip > "$OUT/$chr.pileup.tsv.gz.tmp" \
  && mv "$OUT/$chr.pileup.tsv.gz.tmp" "$OUT/$chr.pileup.tsv.gz" || { echo "TRUTH-FAIL $chr pileup"; exit 1; }
fi
samples=ADMERA,UCLA,GEN50,OTHER
runs=""
for s in COSMO-SUPER COSMO3 COSMO3-DS0.1 COSMO3-DS0.25 COSMO3-DS0.5 COSMO3-DS1.0 admera-0001 UCLA-20001 UCLA-6164; do
  runs="$runs $s=$(ls work/$s/analysis/$G/*_imputed_dog10k.bcf)"
done
python analysis/cosmo_truth/score_chrom.py "$chr" "$OUT/$chr.pileup.tsv.gz" "$samples" "$OUT/$chr" $runs \
  && echo "TRUTH-DONE $chr" || { echo "TRUTH-FAIL $chr score"; exit 1; }
