#!/bin/bash
# Panel-site dosage extract for the SRA top-up genomes (breed panel v2.1): one GLIMPSE BCF per task -> sra_ds/<RUN>.tsv.gz
# (CHROM POS REF ALT DS), the same extract validate_ucla.py reads for customer dogs.
#   ls work_sra/*/analysis/glimpse2/*_imputed_dog10k.bcf > logs/sra_bcfs.txt
#   qsub -t 1-N -tc 60 cluster/sra-ds-extract.sh logs/sra_bcfs.txt
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=2G,h_rt=0:30:00
#$ -pe shared 1
set -uo pipefail
D="${SGE_O_WORKDIR:?run via qsub from the repo root}"; cd "$D"
bcf=$(sed -n "${SGE_TASK_ID}p" "$1"); [ -n "$bcf" ] || exit 0
run=$(basename "$(dirname "$(dirname "$(dirname "$bcf")")")")
mkdir -p sra_ds
[ -s logs/panel_pos.tsv ] || awk 'NR>1{print $1"\t"$2}' breed_panel/sites.tsv > logs/panel_pos.tsv
envs/genomics/bin/bcftools query -R logs/panel_pos.tsv -f '%CHROM\t%POS\t%REF\t%ALT\t[%DS]\n' "$bcf" | gzip > "sra_ds/$run.tsv.gz.tmp" \
  && mv "sra_ds/$run.tsv.gz.tmp" "sra_ds/$run.tsv.gz" && echo "SRADS-DONE $run $(zcat sra_ds/$run.tsv.gz | wc -l)"
