#$ -cwd
#$ -j y
#$ -o logs/lai_query.$TASK_ID.log
#$ -l h_data=4G,h_rt=2:00:00
cd $SGE_O_WORKDIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=$PWD/envs/genomics/lib
S=$PWD/envs/genomics/bin
mkdir -p lai/query
dogs=(COSMO3 COSMO3-DS1.0 COSMO3-DS0.5 COSMO3-DS0.25 COSMO3-DS0.1)
d=${dogs[$((SGE_TASK_ID-1))]}
b=$(ls work/$d/analysis/glimpse2/*_imputed_dog10k.bcf | head -1)
# phased imputed GT at the 131k panel sites, sample renamed to the dog
[ -s lai/query/sites.pos ] || awk "NR>1{print \$1\"\t\"\$2}" breed_panel/sites.tsv > lai/query/sites.pos
$S/bcftools view -R lai/query/sites.pos -Ou "$b" | $S/bcftools annotate -x INFO,^FORMAT/GT -Oz -o lai/query/$d.131k.vcf.gz
$S/bcftools index -f lai/query/$d.131k.vcf.gz
echo "LAI-QUERY-DONE $d sites=$($S/bcftools view -H lai/query/$d.131k.vcf.gz | wc -l)"
