#!/bin/bash
# GLIMPSE2 needs AC/AN in the reference: add them to the test panels build_ref.sbatch stripped.
#   sbatch -A pi_matteop -p campus -a 1-7 analysis/imputation_panel/fill_acan.sbatch
#SBATCH -J panel_acan
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH -o logs/panel_acan.%A.%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs
B=$D/envs/genomics/bin
O=$D/scratch/panel_test/ref
f=$(echo both_chr18 both_chr27 both_chr38 dapx_chr18 dapx_chr27 dapx_chr38 d10k_chr27 | cut -d' ' -f${SLURM_ARRAY_TASK_ID:?})
$B/bcftools +fill-tags $O/$f.bcf -Ob -o $O/$f.acan.bcf -- -t AC,AN
$B/bcftools index -f $O/$f.acan.bcf
mv $O/$f.acan.bcf $O/$f.bcf; mv $O/$f.acan.bcf.csi $O/$f.bcf.csi
echo "ACAN-DONE $f"
