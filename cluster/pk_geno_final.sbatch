#!/bin/bash
# Concatenate the 38 per-chromosome parts into one BCF, index it, write a PLINK
# 1.9 bed/bim/fam copy and the per-sample table. Runs after pk_geno_chr.
#SBATCH -J pkgenofin
#SBATCH -A pi_mcdb -p pi_shared --qos=pi_shared24
#SBATCH -c 4 --mem=48G -t 23:00:00
set -euo pipefail
module load bcftools
D=/u/project/pellegrini/gkislik/dogs
OUT=${OUT:-/u/project/pellegrini_archive/data/prosperkits/genotypes_dog10k}
NAME=prosperkits_dog10k_sites
cd "$OUT"
for i in $(seq 1 38); do
  [ -f parts/chr$i.done ] || { echo "missing chr$i" >&2; exit 1; }
  echo parts/chr$i.bcf
done > parts/concat.txt
bcftools concat --naive -f parts/concat.txt -Ob -o $NAME.bcf
bcftools index --threads 4 $NAME.bcf
bcftools index -n $NAME.bcf
"$HOME/plink" --bcf $NAME.bcf --dog --double-id --keep-allele-order \
  --vcf-idspace-to _ --threads 4 --memory 40000 --make-bed --out $NAME
OPENBLAS_NUM_THREADS=1 "$D/envs/genomics/bin/python" "$D/cluster/pk_geno_samples.py" \
  samples.txt "$D/results_prosper" "$D/sample_sheet.prosper.tsv" \
  "$D/sample_sheet.batch2.tsv" > ${NAME}_samples.tsv
md5sum $NAME.bcf $NAME.bed > md5.txt
echo FINAL-DONE
