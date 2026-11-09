#!/bin/bash
# Breed composition for every DAP dog + training tables (build_dap.py), after extract.sbatch.
#   sbatch -A pi_matteop -p campus analysis/weight_model/dap/build.sbatch
#SBATCH -J dap_wb
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 6:00:00
#SBATCH -o logs/dap_wb.%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
D=/u/project/pellegrini/gkislik/dogs
$D/envs/genomics/bin/python3 $D/analysis/weight_model/dap/build_dap.py $D/dap_weight $D/breed_panel_v23 \
  /u/project/pellegrini/gkislik/dog_aging_project/data/pheno/dogs.tsv $D/logs/coat_b1101_list.tsv \
  $D/analysis_weights.tsv $D/dap_weight/cust_age.tsv $D/dap_weight/build
echo DAPWB-DONE
