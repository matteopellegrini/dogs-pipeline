#!/bin/bash
# COSMO-SUPER read truth in the regions the old panel never imputed (the 169 skipped chunks, failed_chunks.bed),
# at the v2 panel's sites, same rules as cluster/cosmo-truth.sh. One autosome per task.
#   sbatch -A pi_matteop -p campus -a 1-38 analysis/imputation_panel/truth_gap.sbatch
#SBATCH -J imp2_gap
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 8:00:00
#SBATCH -o logs/imp2_gap.%A.%a.log
set -uo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs; cd $D
export PATH=$D/envs/genomics/bin:$PATH
chr=chr${SLURM_ARRAY_TASK_ID:?}; OUT=$D/scratch/imp2_cosmo/gap; mkdir -p $OUT
awk -v c=$chr '$1==c' ref_panel_v2/failed_chunks_old.bed > $OUT/$chr.bed
[ -s $OUT/$chr.bed ] || { echo "GAP-NONE $chr"; exit 0; }
bcftools query -R $OUT/$chr.bed -f '%CHROM\t%POS\n' ref_panel_v2/AutoAndXPAR.Dog10K_DAP.bcf | uniq > $OUT/$chr.sites
samtools view -h /u/project/pellegrini_archive/data/dogs_bams/COSMO-SUPER.markdup.bam -L $OUT/$chr.bed \
  | awk 'BEGIN{OFS="\t"} /^@RG/{next}
         /^@/{print; next}
         !h{print "@RG\tID:ADMERA\tSM:ADMERA"; print "@RG\tID:UCLA\tSM:UCLA"; print "@RG\tID:GEN50\tSM:GEN50"; print "@RG\tID:OTHER\tSM:OTHER"; h=1}
         {g = $1 ~ /^ML15/ ? "ADMERA" : $1 ~ /^V350/ ? "UCLA" : $1 ~ /^LH00427/ ? "GEN50" : "OTHER"
          for (i = 12; i <= NF; i++) if ($i ~ /^RG:Z:/) $i = "RG:Z:" g
          print}' \
  | bcftools mpileup -f canFam4_plusY.fa -T $OUT/$chr.sites -q 30 -Q 20 -a AD -d 1000 -Ou - 2>/dev/null \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' | gzip > $OUT/$chr.pileup.tsv.gz
echo "GAP-DONE $chr $(zcat $OUT/$chr.pileup.tsv.gz | wc -l)"
