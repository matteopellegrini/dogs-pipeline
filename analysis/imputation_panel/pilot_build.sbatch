#!/bin/bash
# Imputation-panel v3 pilot, one chromosome per task: phase the added dogs against panel v2 and build test panels
# on panel v2's exact sites.
#   pv2  = panel v2 (Dog10K + Dog Aging Project), as used in production
#   pex  = pv2 + the existing high-depth SRA genomes (pilot_existing_dogs.tsv, from their sra_full gVCFs)
#   pall = pex + the newly streamed pilot dogs (pilot_v3_dogs.tsv role=panel, sra_panel/<BS>.panel_sites.bcf)
# Test dogs (role=test, SRR40570403) are in none of them. A BioSample listed under two breed groups enters once.
#   sbatch -A pi_matteop -p campus -a 18,27,38 -c 8 --mem=48G -t 12:00:00 analysis/imputation_panel/pilot_build.sbatch
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
  | $B/bgzip -c > $P/add_$c.vcf.gz
$B/bcftools index -f $P/add_$c.vcf.gz
awk '{s+=$3; n++} END{printf "mean missing per dog %.3f over %d dogs\n", s/n, n}' $P/missing_$c.txt

# SHAPEIT5 wants AC/AN in the reference; panel v2's BCFs do not carry them
[ -s $P/refac_$c.bcf.csi ] || { $B/bcftools +fill-tags $REF -Ob -o $P/refac_$c.bcf -- -t AC,AN && $B/bcftools index -f $P/refac_$c.bcf; }
$S/SHAPEIT5_phase_common --input $P/add_$c.vcf.gz --reference $P/refac_$c.bcf --region $c --thread $THREADS \
  --output $P/add_phased_$c.bcf --log $P/phase_$c.log > /dev/null
$B/bcftools index -f $P/add_phased_$c.bcf

cut -f1 analysis/imputation_panel/pilot_existing_dogs.tsv > $P/existing_samples.txt
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
