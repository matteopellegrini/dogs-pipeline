#$ -cwd
#$ -j y
#$ -o logs/lai_phase_parker.log
#$ -l h_data=8G,h_rt=12:00:00
#$ -pe shared 4
cd $SGE_O_WORKDIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export PATH=/u/local/apps/java/jdk-17.0.12/bin:$PATH
export LD_LIBRARY_PATH=$PWD/envs/genomics/lib
S=$PWD/envs/genomics/bin
# wait for the Parker VCF builder (login node) to finish
until grep -q "^kept" lai/make_parker_vcf.log 2>/dev/null; do sleep 60; done
cd lai
$S/bcftools index -f -t dog10k_131k.vcf.gz 2>/dev/null || (zcat dog10k_131k.vcf.gz | $S/bgzip -c > dog10k_131k.bgz.vcf.gz && mv dog10k_131k.bgz.vcf.gz dog10k_131k.vcf.gz && $S/bcftools index -f -t dog10k_131k.vcf.gz)
zcat parker_131k.vcf.gz | $S/bgzip -c > parker_131k.bgz.vcf.gz && mv parker_131k.bgz.vcf.gz parker_131k.vcf.gz && $S/bcftools index -f -t parker_131k.vcf.gz
# Beagle 5.5: phase the unphased Parker genotypes using the phased Dog10K
# haplotypes as reference (same sites, same allele system). Constant-rate
# map if no canine map is available (map= omitted => 1 cM/Mb default).
java -Xmx24g -jar ../envs/lai/beagle.jar gt=parker_131k.vcf.gz ref=dog10k_131k.vcf.gz out=parker_131k.phased nthreads=4 impute=false
$S/bcftools index -f -t parker_131k.phased.vcf.gz
echo "LAI-PHASE-DONE $($S/bcftools query -l parker_131k.phased.vcf.gz | wc -l) samples"
