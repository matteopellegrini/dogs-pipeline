#$ -cwd
#$ -j y
#$ -o logs/rel_extract_b2.$TASK_ID.log
#$ -l h_data=4G,h_rt=1:00:00
cd $SGE_O_WORKDIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=$PWD/envs/genomics/lib:${LD_LIBRARY_PATH:-}
S=$PWD/envs/genomics/bin
s=$(sed -n "${SGE_TASK_ID}p" relatives/batch2_eligible.txt)
[ -n "$s" ] || exit 0
b=$(ls work_prosper/$s/analysis/glimpse2/*_imputed_dog10k.bcf 2>/dev/null | head -1)
[ -n "$b" ] || { echo "no bcf for $s"; exit 0; }
mkdir -p relatives/batch2
$S/bcftools query -R relatives/sites.pos -f "%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n" "$b" \
  | gzip > relatives/batch2/$s.ds.tsv.gz
echo B2-EXTRACT-DONE $s rows=$(zcat relatives/batch2/$s.ds.tsv.gz | wc -l)
