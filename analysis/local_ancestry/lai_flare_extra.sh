#$ -cwd
#$ -j y
#$ -o logs/lai_flare_extra.log
#$ -l h_data=8G,h_rt=6:00:00
#$ -pe shared 4
# FLARE on additional queries (the two MGI replicates of Cosmo sequenced at
# UCLA: COSMO-ZYMO = ucla-20001, COSMO-GENOTEK = ucla-6164), reusing the
# reference built by lai_flare.sh. Same two runs per dog (noprior / prior).
cd $SGE_O_WORKDIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export PATH=/u/local/apps/java/jdk-17.0.12/bin:$PATH
export LD_LIBRARY_PATH=$PWD/envs/genomics/lib
S=$PWD/envs/genomics/bin
until [ -f lai/ref_panel.vcf.gz.tbi ]; do sleep 120; done
cd lai
mkdir -p flare_out priors
for d in $(cat query/extra_queries.txt); do
  sid=$($S/bcftools query -l query/$d.131k.vcf.gz | head -1)
  ../envs/genomics/bin/python3 lai_priors.py query/$d.breed_result.json "$sid" priors/$d.tsv
  java -Xmx24g -jar ../envs/lai/flare.jar ref=ref_panel.vcf.gz ref-panel=ref_panel_map.tsv \
       gt=query/$d.131k.vcf.gz map=constant_1cM_per_Mb.map out=flare_out/$d.noprior probs=true nthreads=4 seed=1 \
    && echo "LAI-FLARE-DONE $d noprior" || echo "LAI-FLARE-FAIL $d noprior"
  java -Xmx24g -jar ../envs/lai/flare.jar ref=ref_panel.vcf.gz ref-panel=ref_panel_map.tsv \
       gt=query/$d.131k.vcf.gz map=constant_1cM_per_Mb.map out=flare_out/$d.prior probs=true \
       gt-ancestries=priors/$d.tsv nthreads=4 seed=1 \
    && echo "LAI-FLARE-DONE $d prior" || echo "LAI-FLARE-FAIL $d prior"
done
