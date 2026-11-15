#!/bin/bash
# Imputation-panel v3 pilot, one chromosome per task: phase the added dogs against panel v2 (Beagle 5.5) and build test panels
# on panel v2's exact sites.
#   pv2  = panel v2 (Dog10K + Dog Aging Project), as used in production
#   pex  = pv2 + the existing high-depth SRA genomes (pilot_existing_dogs.tsv, from their sra_full gVCFs)
#   pall = pex + the newly streamed pilot dogs (pilot_v3_dogs.tsv role=panel, sra_panel/<BS>.panel_sites.bcf)
# Test dogs (role=test, SRR40570403) are in none of them. A BioSample listed under two breed groups enters once.
#   sbatch -A pi_matteop -p campus -a 18,27,38 -c 8 --mem=128G -t 12:00:00 analysis/imputation_panel/pilot_build.sbatch   (48G OOMs in SHAPEIT5 burn-in)
#SBATCH -J pilot_build
#SBATCH -o logs/pilot_build.%A.%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 THREADS=${SLURM_CPUS_PER_TASK:-8}
D=/u/project/pellegrini/gkislik/dogs; cd $D
B=$D/envs/genomics/bin; S=$D/envs/shapeit5/bin
T=$D/scratch/panel_test; P=$T/pilot; mkdir -p $P
c=chr${SLURM_ARRAY_TASK_ID:?}
REF=$D/ref_panel_v2/chr/$c.bcf

# Added dogs: existing gVCF genomes + every finished streamed panel dog.
cp analysis/imputation_panel/pilot_existing_dogs.tsv $P/dogs_$c.tsv
grep -v '^#' analysis/imputation_panel/pilot_v3_dogs.tsv | awk -F'\t' '$4=="panel" && !seen[$1]++{print $1}' | while read -r bs; do
  [ -s sra_panel/$bs.done ] && printf 'SRA-%s\tsites\tsra_panel/%s.panel_sites.bcf\tnew\n' "$bs" "$bs"
done >> $P/dogs_$c.tsv
cut -f1,2,3 $P/dogs_$c.tsv > $P/dogs_$c.in.tsv
echo "$c: $(grep -c gvcf $P/dogs_$c.tsv) existing + $(grep -c -w new $P/dogs_$c.tsv) new dogs"

$B/python3 analysis/imputation_panel/pilot_genotypes.py $c ref_panel_v2/sites/$c.tsv.gz $P/dogs_$c.in.tsv 2>$P/missing_$c.txt \
  | $B/bcftools +fill-tags - -Oz -o $P/add_$c.vcf.gz -- -t AC,AN   # SHAPEIT5 needs AC/AN in the input too
$B/bcftools index -f $P/add_$c.vcf.gz
awk '{s+=$3; n++} END{printf "mean missing per dog %.3f over %d dogs\n", s/n, n}' $P/missing_$c.txt

# Phase the added dogs against panel v2 with Beagle 5.5. SHAPEIT5 phase_common OOM-killed even at 128G here:
# ~45% of panel sites are rare and the added genotypes are ~14% missing, which it is not built for. Beagle phases
# every site and fills the missing genotypes. Its ref must be VCF/bref3, so panel v2 is written once as VCF (+AC/AN).
[ -s $P/refac_$c.vcf.gz.csi ] || { $B/bcftools +fill-tags $REF -Oz -o $P/refac_$c.vcf.gz -- -t AC,AN && $B/bcftools index -f $P/refac_$c.vcf.gz; }
java -Xmx110g -jar $D/envs/lai/beagle.jar gt=$P/add_$c.vcf.gz ref=$P/refac_$c.vcf.gz chrom=$c nthreads=$THREADS \
  impute=false out=$P/add_phased_$c > $P/phase_$c.log 2>&1
$B/bcftools view -Ob -o $P/add_phased_$c.bcf $P/add_phased_$c.vcf.gz && rm -f $P/add_phased_$c.vcf.gz
$B/bcftools index -f $P/add_phased_$c.bcf

cut -f1 analysis/imputation_panel/pilot_existing_dogs.tsv > $P/existing_samples.txt
[ -s $P/refac_$c.bcf.csi ] || { $B/bcftools view -Ob -o $P/refac_$c.bcf $P/refac_$c.vcf.gz && $B/bcftools index -f $P/refac_$c.bcf; }
ln -sf $P/refac_$c.bcf $T/ref/pv2_$c.bcf; ln -sf $P/refac_$c.bcf.csi $T/ref/pv2_$c.bcf.csi
$B/bcftools view -S $P/existing_samples.txt -Ob -o $P/add_ex_$c.bcf $P/add_phased_$c.bcf && $B/bcftools index -f $P/add_ex_$c.bcf
for pair in pex:$P/add_ex_$c.bcf pall:$P/add_phased_$c.bcf; do
  name=${pair%%:*}; add=${pair#*:}
  $B/bcftools merge -m none $REF $add -Ou | $B/bcftools +fill-tags -Ou -- -t AC,AN \
    | $B/bcftools view -Ob -o $T/ref/${name}_$c.bcf
  $B/bcftools index -f $T/ref/${name}_$c.bcf
  echo "$c $name: $($B/bcftools query -l $T/ref/${name}_$c.bcf | wc -l) dogs, $($B/bcftools index -n $T/ref/${name}_$c.bcf) sites (pv2 $($B/bcftools index -n $REF))"
done
echo "BUILD-DONE $c"
