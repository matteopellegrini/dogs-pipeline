#!/bin/bash
# Rare-breed check of the Dog10K+DAP panel: SRR40570403 (Finnish Lapphund, 21.6x, not in either panel) imputed from
# a 10M-spot low-pass slice at 0.25/0.5/1.0x with Dog10K alone (d10kx) and Dog10K+DAP (both); truth from its full gVCF.
#   sbatch -A pi_matteop -p campus -a 1-3 -c 1 --mem=8G -t 2:00:00 -J sra_score analysis/imputation_panel/sra_score.sbatch
#SBATCH -o logs/%x.%A.%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs; cd $D
T=$D/scratch/panel_test; R=SRR40570403
c=$(echo chr18 chr27 chr38 | cut -d' ' -f${SLURM_ARRAY_TASK_ID:?})
PY=$D/envs/genomics/bin/python3
mkdir -p $T/sra_truth $T/score_sra
$PY analysis/imputation_panel/sra_truth.py sra_full/$R.g.vcf.gz $c $T/ref/raf_$c.tsv.gz $T/sra_truth/${R}_$c.pileup.tsv.gz
runs=()
for d in 0.25x 0.5x 1.0x; do for p in d10kx both; do runs+=("${d}_$p=$T/imp/${c}_$R-${d}_$p.bcf"); done; done
TRUTH_PILEUP=$T/sra_truth/${R}_$c.pileup.tsv.gz TRUTH_SAMPLES=SRA TRUTH_DROP=NONE \
  $PY analysis/imputation_panel/score_panels.py $c $T/score_sra/$c.tsv "${runs[@]}"
