#!/bin/bash
# Level 2: merge the 9 kit groups of one chromosome, add cohort AN/AC/AF.
#SBATCH -J pkgenochr
#SBATCH -A pi_mcdb -p pi_shared --qos=pi_shared24
#SBATCH -c 1 --mem=4G -t 6:00:00
#SBATCH --array=1-38
set -euo pipefail
module load bcftools
OUT=${OUT:-/u/project/pellegrini_archive/data/prosperkits/genotypes_dog10k}
c=chr$SLURM_ARRAY_TASK_ID
ls "$OUT"/parts/$c.g*.bcf > "$OUT/parts/$c.list"
[ "$(wc -l < "$OUT/parts/$c.list")" = 9 ] || { echo "$c: groups missing" >&2; exit 1; }
bcftools merge -l "$OUT/parts/$c.list" -Ou \
  | bcftools +fill-tags -Ob -o "$OUT/parts/$c.bcf" -- -t AN,AC,AF
echo "$c records $(bcftools view -H "$OUT/parts/$c.bcf" | wc -l) samples $(bcftools query -l "$OUT/parts/$c.bcf" | wc -l)"
touch "$OUT/parts/$c.done"
