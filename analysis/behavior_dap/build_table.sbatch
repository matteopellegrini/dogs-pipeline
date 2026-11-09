#!/bin/bash
#SBATCH -J dap_btab
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -o logs/dap_btab.%j.log
set -euo pipefail
D=/u/project/pellegrini/gkislik/dogs
export OPENBLAS_NUM_THREADS=2
$D/envs/genomics/bin/python3 $D/analysis/behavior_dap/build_table.py $D/dap_weight $D/breed_panel_v23 \
  $D/reference_json/prs_lmm_betas.json.gz /u/project/pellegrini/gkislik/dog_aging_project/data/pheno/dogs.tsv \
  $D/dap_behavior/table.npz
echo DAPBTAB-DONE
