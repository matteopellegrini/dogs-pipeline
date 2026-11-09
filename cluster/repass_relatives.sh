#!/bin/bash
# Stage-9b repass of relatives_result.json against the current relatives_ref,
# over a rerun list (name, barcode, pub dir, work dir), strided so ~100 tasks
# cover any list. cluster/relatives_9b.py is a verbatim extract of the stage-9b
# heredoc; this prelude mirrors the stage's bcftools dosage query. Keeps
# relatives_result.pre-<tag>.json once; republish with
#   FILES=relatives_result.json bash cluster/republish_files.sh <tag> <list>
#   sbatch -A pi_matteop -p campus -a 1-100 -c 1 --mem=12G -t 6:00:00 -J rel9b \
#          --export=ALL,STRIDE=100 cluster/repass_relatives.sh <list> <tag>
#SBATCH -o logs/%x.%A.%a.log
cd "$SLURM_SUBMIT_DIR"; D=$PWD; ENV=$D/envs/genomics
export LD_LIBRARY_PATH=$ENV/lib:${LD_LIBRARY_PATH:-} OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONIOENCODING=utf-8
LIST=${1:?list}; TAG=${2:?tag}; N=$(wc -l < "$LIST")
for ((i=SLURM_ARRAY_TASK_ID; i<=N; i+=${STRIDE:-100})); do
  row=$(sed -n "${i}p" "$LIST"); pk=$(echo "$row" | cut -f1); PUB=$(echo "$row" | cut -f3); WORK=$(echo "$row" | cut -f4)
  IMPUTED_BCF=$(ls "$WORK"/glimpse2/*_imputed_dog10k.bcf 2>/dev/null | head -1)
  [ -f "$IMPUTED_BCF" ] && [ -f "$PUB/qc_result.json" ] || { echo "REL9B-SKIP $pk (no bcf/qc)"; continue; }
  tmp=$(mktemp -d "$D/logs/rel9b.XXXXXX")
  awk 'NR>1{print $1"\t"$2}' "$D/breed_panel/sites.tsv" > "$tmp/pos"
  $ENV/bin/bcftools query -R "$tmp/pos" -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' "$IMPUTED_BCF" > "$tmp/ds" \
    || { echo "REL9B-FAIL $pk (bcftools)"; rm -rf "$tmp"; continue; }
  [ -f "$PUB/relatives_result.json" ] && [ ! -f "$PUB/relatives_result.pre-$TAG.json" ] \
    && cp "$PUB/relatives_result.json" "$PUB/relatives_result.pre-$TAG.json"
  REL_DS="$tmp/ds" REL_REF="$D/relatives_ref" PUB_DIR="$PUB" SAMPLE="$pk" \
    $ENV/bin/python3 cluster/relatives_9b.py && echo "REL9B-DONE $pk" || echo "REL9B-FAIL $pk (python)"
  rm -rf "$tmp"
done
