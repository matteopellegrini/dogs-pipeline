#!/bin/bash
# Imputation reference v2 = Dog10K (1,929) + Dog Aging Project (7,828), one chromosome per task.
#   sbatch -A pi_matteop -p campus -a 1-39 analysis/imputation_panel/build_panel.sbatch   (39 = chrX PAR)
# then: bash analysis/imputation_panel/finish_panel.sh
#
# Per chromosome:
#  - Dog10K from the production panel, exact duplicates removed (the July concat -a doubled 1.2M records) and
#    every site with AC==0 or AC==AN removed: GLIMPSE2_phase aborts any chunk holding one ("AC/AN INFO fields
#    in VCF are inconsistent with GT field"), which is why 169 of 813 chunks never ran.
#  - DAP GT only, merged on identical CHROM/POS/REF/ALT; sites missing from either side are dropped. The 20
#    hand-added disease sites that DAP lacks are imputed separately from a Dog10K-only window panel.
#  - AC/AN recomputed; sites monomorphic in the merged panel dropped.
#SBATCH -J panel_v2
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 20:00:00
#SBATCH -o logs/panel_v2.%A.%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs
B=$D/envs/genomics/bin
O=$D/ref_panel_v2
mkdir -p $O/chr $O/tmp
i=${SLURM_ARRAY_TASK_ID:?}
c=chr$i; [[ $i == 39 ]] && c=chrX
[[ -s $O/chr/$c.bcf.csi ]] && { echo "$c exists"; exit 0; }
T=$O/tmp/$c; mkdir -p $T
$B/bcftools view -r $c -Ou $D/dog10k_panel/AutoAndXPAR.Dog10K.phased_plus_disease_rh.bcf 2>/dev/null \
  | $B/bcftools norm -d exact -Ou 2>/dev/null \
  | $B/bcftools annotate -x INFO,^FORMAT/GT -Ou \
  | $B/bcftools +fill-tags -Ou -- -t AC,AN \
  | $B/bcftools view -e 'INFO/AC==0 || INFO/AC==INFO/AN' -Ob -o $T/d10k.bcf
$B/bcftools index -f $T/d10k.bcf
n10=$($B/bcftools index -n $T/d10k.bcf)
if [[ $n10 == 0 ]]; then echo "$c: no Dog10K records"; exit 0; fi
$B/bcftools view -r $c -Ou $D/dap_panel/dap.vcf.gz 2>/dev/null \
  | $B/bcftools annotate -x INFO,^FORMAT/GT -Ob -o $T/dap.bcf
$B/bcftools index -f $T/dap.bcf
$B/bcftools merge -m none --force-samples $T/d10k.bcf $T/dap.bcf -Ou \
  | $B/bcftools view -i 'N_MISSING=0' -Ou \
  | $B/bcftools +fill-tags -Ou -- -t AC,AN \
  | $B/bcftools view -e 'INFO/AC==0 || INFO/AC==INFO/AN' -Ob -o $T/merged.bcf
$B/bcftools index -f $T/merged.bcf
nm=$($B/bcftools index -n $T/merged.bcf)
ns=$($B/bcftools query -l $T/merged.bcf | wc -l)
# Dog10K-only clean copy is kept for the added-site windows.
mv $T/d10k.bcf $O/chr/$c.dog10k_clean.bcf; mv $T/d10k.bcf.csi $O/chr/$c.dog10k_clean.bcf.csi
mv $T/merged.bcf $O/chr/$c.bcf; mv $T/merged.bcf.csi $O/chr/$c.bcf.csi
rm -rf $T
echo "$c: dog10k_clean $n10, merged $nm sites x $ns samples"
echo "PANEL-DONE $c"
