#!/bin/bash
# Score the v2-panel Cosmo imputations in the regions the old pipeline skipped (169 chunks), against the
# COSMO-SUPER truth from truth_gap.sbatch. The old runs have no genotypes there, so only the new runs are scored.
#   sbatch -A pi_matteop -p campus -a 1-38 analysis/imputation_panel/score_gap.sbatch
#SBATCH -J imp2_gapscore
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 6:00:00
#SBATCH -o logs/imp2_gapscore.%A.%a.log
set -uo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs; cd $D
chr=chr${SLURM_ARRAY_TASK_ID:?}; OUT=$D/scratch/imp2_cosmo/score_gap; mkdir -p $OUT
P=scratch/imp2_cosmo/gap/$chr.pileup.tsv.gz
[ -s $P ] || { echo "GAP-NONE $chr"; exit 0; }
runs=""
for s in COSMO3-DS0.1 COSMO3-DS0.25 COSMO3-DS0.5 COSMO3-DS1.0 COSMO3 admera-0001; do
  runs="$runs new_$s=$(ls scratch/imp2_cosmo/$s/analysis/glimpse2/*_imputed_dog10k.bcf)"
done
envs/genomics/bin/python analysis/cosmo_truth/score_chrom.py "$chr" $P ADMERA,UCLA,GEN50,OTHER "$OUT/$chr" $runs && echo "SCORE-DONE $chr"
