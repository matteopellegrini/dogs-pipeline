#!/bin/bash
# Level 1 of the ProsperKits genotype matrix: one chromosome x one group of
# <=200 kits (opening all 1,775 kit BCFs at once needs >16 GB). Keeps phased GT
# + panel RAF and renames samples from "markdup" to pk-<kit>.
# Array index i: chromosome (i-1)/9+1, group (i-1)%9. Inputs from pk_geno_matrix.sh.
#SBATCH -J pkgeno
#SBATCH -A pi_mcdb -p pi_shared --qos=pi_shared24
#SBATCH -c 1 --mem=4G -t 6:00:00
#SBATCH --array=1-342
set -euo pipefail
module load bcftools
OUT=${OUT:-/u/project/pellegrini_archive/data/prosperkits/genotypes_dog10k}
i=$((SLURM_ARRAY_TASK_ID - 1))
c=chr$((i / 9 + 1)); g=$(printf %02d $((i % 9)))
p=$OUT/parts/$c.g$g
[ -f "$p.bcf.csi" ] && exit 0
bcftools merge --force-samples -r "$c" -l "$OUT/parts/files.g$g" -Ou \
  | bcftools annotate -x INFO/AF,INFO/INFO,FORMAT/DS,FORMAT/GP -Ob -o "$p.tmp.bcf"
bcftools reheader -s "$OUT/parts/samples.g$g" -o "$p.bcf" "$p.tmp.bcf"
rm -f "$p.tmp.bcf"
bcftools index -f "$p.bcf"
