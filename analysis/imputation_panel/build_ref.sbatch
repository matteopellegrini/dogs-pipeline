#!/bin/bash
# Test reference panels for one chromosome: Dog10K alone vs Dog10K + Dog Aging Project
# (and DAP alone), all restricted to the sites both panels carry, GT only.
#   sbatch -A pi_matteop -p campus -a 18,27,38 analysis/imputation_panel/build_ref.sbatch
#SBATCH -J panel_ref
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -o logs/panel_ref.%A.%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs
B=$D/envs/genomics/bin
O=$D/scratch/panel_test/ref
mkdir -p $O
c=chr${SLURM_ARRAY_TASK_ID:?}
[[ -s $O/both_$c.bcf.csi ]] && { echo "$c exists"; exit 0; }
$B/bcftools view -r $c -Ou $D/dog10k_panel/AutoAndXPAR.Dog10K.phased_plus_disease_rh.bcf 2>/dev/null \
  | $B/bcftools annotate -x INFO,^FORMAT/GT -Ob -o $O/d10k_$c.bcf && $B/bcftools index -f $O/d10k_$c.bcf
$B/bcftools view -r $c -Ou $D/dap_panel/dap.vcf.gz 2>/dev/null \
  | $B/bcftools annotate -x INFO,^FORMAT/GT -Ob -o $O/dap_$c.bcf && $B/bcftools index -f $O/dap_$c.bcf
# Records present in both (same REF/ALT) with no missing genotypes.
$B/bcftools merge -m none --force-samples $O/d10k_$c.bcf $O/dap_$c.bcf -Ou \
  | $B/bcftools view -i 'N_MISSING=0' -Ob -o $O/both_$c.tmp.bcf
$B/bcftools index -f $O/both_$c.tmp.bcf
$B/bcftools query -l $O/d10k_$c.bcf > $O/d10k_samples.txt
$B/bcftools query -l $O/both_$c.tmp.bcf | tail -n +$(( $(wc -l < $O/d10k_samples.txt) + 1 )) > $O/dap_samples_$c.txt
$B/bcftools view -S $O/d10k_samples.txt -Ob -o $O/d10kx_$c.bcf $O/both_$c.tmp.bcf && $B/bcftools index -f $O/d10kx_$c.bcf
$B/bcftools view -S $O/dap_samples_$c.txt -Ob -o $O/dapx_$c.bcf $O/both_$c.tmp.bcf && $B/bcftools index -f $O/dapx_$c.bcf
$B/bcftools query -f '%POS\t%REF\t%ALT\t%AC\t%AN\n' $O/d10kx_$c.bcf | gzip -1 > $O/raf_$c.tsv.gz
mv $O/both_$c.tmp.bcf $O/both_$c.bcf; mv $O/both_$c.tmp.bcf.csi $O/both_$c.bcf.csi
echo "$c: d10k $($B/bcftools index -n $O/d10k_$c.bcf) dap $($B/bcftools index -n $O/dap_$c.bcf) shared $($B/bcftools index -n $O/both_$c.bcf)"
echo "REF-DONE $c"
