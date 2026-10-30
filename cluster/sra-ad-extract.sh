#!/bin/bash
# Direct allele counts (no imputation) at the breed-panel sites for the SRA top-up genomes, breed panel v2.x:
# one BAM per task -> sra_ad/<RUN>.tsv.gz (CHROM POS REF ALT AD). The imputed dosages (sra-ds-extract.sh) made
# those columns resemble every imputed customer dog; counts from the reads themselves do not.
#   qsub -t 1-N -tc 80 cluster/sra-ad-extract.sh logs/sra_ad_runs.txt      (one RUN accession per line)
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=2G,h_rt=0:30:00
#$ -pe shared 1
set -uo pipefail
D="${SGE_O_WORKDIR:?run via qsub from the repo root}"; cd "$D"
run=$(sed -n "${SGE_TASK_ID}p" "$1"); [ -n "$run" ] || exit 0
bam="work_sra/$run/analysis/markdup.bam"; [ -s "$bam" ] || { echo "SRAAD-SKIP $run (no bam)"; exit 0; }
mkdir -p sra_ad
[ -s logs/panel_pos_sorted.tsv ] || awk 'NR>1{print $1"\t"$2}' breed_panel/sites.tsv | sort -k1,1 -k2,2n > logs/panel_pos_sorted.tsv
envs/genomics/bin/bcftools mpileup -f canFam4_plusY.fa -T logs/panel_pos_sorted.tsv -q 30 -Q 20 -a AD -I -Ou "$bam" 2>/dev/null \
  | envs/genomics/bin/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' | gzip > "sra_ad/$run.tsv.gz.tmp" \
  && mv "sra_ad/$run.tsv.gz.tmp" "sra_ad/$run.tsv.gz" && echo "SRAAD-DONE $run $(zcat sra_ad/$run.tsv.gz | wc -l)"
