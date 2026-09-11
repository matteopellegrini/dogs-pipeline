#$ -cwd
#$ -j y
#$ -o logs/lai_flare.log
#$ -l h_data=8G,h_rt=12:00:00
#$ -pe shared 4
# Local ancestry pilot: FLARE on cosmo3 + its downsamples against the merged
# phased Parker + Dog10K reference (breeds with >= 5 dogs). Two runs per dog:
# FLARE's own EM (noprior) and the two-step design with the global lasso
# proportions as gt-ancestries prior (prior). Waits for the Beagle phasing
# job and the query extraction array.
cd $SGE_O_WORKDIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
module load java/jdk-17.0.12
export LD_LIBRARY_PATH=$PWD/envs/genomics/lib
S=$PWD/envs/genomics/bin
until grep -q "LAI-PHASE-DONE" logs/lai_phase_parker.log 2>/dev/null; do sleep 120; done
until [ "$(grep -l LAI-QUERY-DONE logs/lai_query.*.log 2>/dev/null | wc -l)" -ge 5 ]; do sleep 60; done
cd lai
$S/bcftools merge -Oz -o ref_merged.vcf.gz parker_131k.phased.vcf.gz dog10k_131k.vcf.gz
$S/bcftools view -S ref_keep_samples.txt --force-samples -Oz -o ref_panel.vcf.gz ref_merged.vcf.gz
$S/bcftools index -f -t ref_panel.vcf.gz
echo "reference: $($S/bcftools query -l ref_panel.vcf.gz | wc -l) samples, $($S/bcftools view -H ref_panel.vcf.gz | wc -l) sites"
mkdir -p flare_out priors
declare -A RES=([COSMO3]=cosmo3 [COSMO3-DS1.0]=cosmo3-ds1.0 [COSMO3-DS0.5]=cosmo3-ds0.5 [COSMO3-DS0.25]=cosmo3-ds0.25 [COSMO3-DS0.1]=cosmo3-ds0.1)
for d in COSMO3 COSMO3-DS1.0 COSMO3-DS0.5 COSMO3-DS0.25 COSMO3-DS0.1; do
  sid=$($S/bcftools query -l query/$d.131k.vcf.gz | head -1)
  ../envs/genomics/bin/python3 lai_priors.py ../results/${RES[$d]}/breed_result.json "$sid" priors/$d.tsv
  java -Xmx24g -jar ../envs/lai/flare.jar ref=ref_panel.vcf.gz ref-panel=ref_panel_map.tsv \
       gt=query/$d.131k.vcf.gz map=constant_1cM_per_Mb.map out=flare_out/$d.noprior probs=true nthreads=4 seed=1 \
    && echo "LAI-FLARE-DONE $d noprior" || echo "LAI-FLARE-FAIL $d noprior"
  java -Xmx24g -jar ../envs/lai/flare.jar ref=ref_panel.vcf.gz ref-panel=ref_panel_map.tsv \
       gt=query/$d.131k.vcf.gz map=constant_1cM_per_Mb.map out=flare_out/$d.prior probs=true \
       gt-ancestries=priors/$d.tsv nthreads=4 seed=1 \
    && echo "LAI-FLARE-DONE $d prior" || echo "LAI-FLARE-FAIL $d prior"
done
