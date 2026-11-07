#!/bin/bash
#SBATCH -J gene_labels
#SBATCH -t 02:00:00
#SBATCH --mem=8G
#SBATCH -c 1
#SBATCH -o logs/%x.%j.log
# Rewrite CNV gene labels (cnv50_result / cnv_homdel / cnv_genes) in every result
# dir from the lifted canFam4 gene table (GENES; Hoffman reference_json still holds
# the ROS table until the branch is merged). No publish: republish is a separate step.
D=/u/project/pellegrini/gkislik/dogs
export OPENBLAS_NUM_THREADS=1
dirs=()
for d in "$D"/results_prosper/* "$D"/results/*; do
  [[ -s "$d/cnv50_result.json" || -s "$d/cnv_homdel.json" ]] && dirs+=("$d")
done
echo "result dirs: ${#dirs[@]}"
"$D/envs/genomics/bin/python3" "$D/cluster/backfill_gene_labels.py" "${GENES:-/u/scratch/g/gkislik/genefix/genes_canFam4.json}" "${dirs[@]}"
