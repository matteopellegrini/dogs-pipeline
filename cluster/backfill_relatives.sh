#$ -cwd
#$ -j y
#$ -o logs/rel9b.$TASK_ID.log
#$ -l h_data=4G,h_rt=0:40:00
# Stage-9b repass against the v3 relatives reference (3,804 dogs, built
# 2026-09-08). cluster/relatives_9b.py is a VERBATIM extract of the stage-9b
# python heredoc in run_dog_pipeline.sh (quoted heredoc, so no shell
# interpolation to reproduce); the shell prelude below mirrors the stage's
# own bcftools dosage query. Only relatives_result.json is rewritten; the
# previous file is kept once as relatives_result.pre-v3.json.
cd $SGE_O_WORKDIR
D=$PWD
ENV=$D/envs/genomics
export LD_LIBRARY_PATH=$ENV/lib:${LD_LIBRARY_PATH:-}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # numpy on many-core nodes
export PYTHONIOENCODING=utf-8

row=$(sed -n "${SGE_TASK_ID}p" cluster/repass9b_list.txt)
[ -n "$row" ] || exit 0
pk=$(echo "$row" | cut -f1); PUB=$(echo "$row" | cut -f3); WORK=$(echo "$row" | cut -f4)
IMPUTED_BCF=$(ls "$WORK"/glimpse2/*_imputed_dog10k.bcf 2>/dev/null | head -1)
[ -f "$IMPUTED_BCF" ] && [ -f "$PUB/qc_result.json" ] || { echo "REL9B-SKIP $pk (no bcf/qc)"; exit 0; }

BREED_SITES=$D/breed_panel/sites.tsv
POS=$WORK/relatives_sites.pos
REL_DS=$WORK/relatives_ds.tsv
awk 'NR>1{print $1"\t"$2}' "$BREED_SITES" > "$POS"
$ENV/bin/bcftools query -R "$POS" -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' "$IMPUTED_BCF" > "$REL_DS" \
  || { echo "REL9B-FAIL $pk (bcftools)"; exit 1; }

[ -f "$PUB/relatives_result.json" ] && [ ! -f "$PUB/relatives_result.pre-v3.json" ] \
  && cp "$PUB/relatives_result.json" "$PUB/relatives_result.pre-v3.json"

REL_DS="$REL_DS" REL_REF="$D/relatives_ref" PUB_DIR="$PUB" SAMPLE="$pk" \
  $ENV/bin/python3 cluster/relatives_9b.py || { echo "REL9B-FAIL $pk (python)"; exit 1; }
rm -f "$REL_DS" "$POS"
echo "REL9B-DONE $pk"
