#!/bin/bash
# GLIMPSE2_phase aborts any chunk containing a reference site where every haplotype carries ALT
# (AC == AN; "AC/AN INFO fields in VCF are inconsistent with GT field"). Drop those from a test panel.
#   sbatch -A pi_matteop -p campus -a 1-4 analysis/imputation_panel/drop_fixed.sbatch
#SBATCH -J panel_fix
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH -o logs/panel_fix.%A.%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs
B=$D/envs/genomics/bin
O=$D/scratch/panel_test/ref
f=$(echo dapx_chr18 dapx_chr27 dapx_chr38 d10k_chr27 | cut -d' ' -f${SLURM_ARRAY_TASK_ID:?})
n=$($B/bcftools view -H -i 'INFO/AC==INFO/AN' $O/$f.bcf | wc -l)
$B/bcftools view -e 'INFO/AC==INFO/AN' -Ob -o $O/$f.fix.bcf $O/$f.bcf
$B/bcftools index -f $O/$f.fix.bcf
mv $O/$f.fix.bcf $O/$f.bcf; mv $O/$f.fix.bcf.csi $O/$f.bcf.csi
echo "FIX-DONE $f dropped $n all-ALT sites"
