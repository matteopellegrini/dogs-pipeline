#!/bin/bash
# Score old (production) vs v2-panel imputations of Cosmo against the COSMO-SUPER truth, one autosome per task.
#   sbatch -A pi_matteop -p campus -a 1-38 analysis/imputation_panel/score_cosmo.sbatch
#   then: python analysis/cosmo_truth/aggregate.py-style summary via analysis/imputation_panel/summarise_cosmo.py
#SBATCH -J imp2_score
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 6:00:00
#SBATCH -o logs/imp2_score.%A.%a.log
set -uo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs; cd $D
chr=chr${SLURM_ARRAY_TASK_ID:?}; OUT=$D/scratch/imp2_cosmo/score; mkdir -p $OUT
runs=""
for s in COSMO3-DS0.1 COSMO3-DS0.25 COSMO3-DS0.5 COSMO3-DS1.0 COSMO3 admera-0001; do
  runs="$runs old_$s=$(ls work/$s/analysis/glimpse2/*_imputed_dog10k.bcf)"
  runs="$runs new_$s=$(ls scratch/imp2_cosmo/$s/analysis/glimpse2/*_imputed_dog10k.bcf)"
done
envs/genomics/bin/python analysis/cosmo_truth/score_chrom.py "$chr" analysis/cosmo_truth/out/$chr.pileup.tsv.gz \
  ADMERA,UCLA,GEN50,OTHER "$OUT/$chr" $runs && echo "SCORE-DONE $chr"
