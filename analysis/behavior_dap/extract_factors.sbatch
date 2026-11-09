#!/bin/bash
# Dog Aging Project -> per-chromosome Darwin's Ark behaviour-factor PRS (extract_factors.py).
#   sbatch -A pi_matteop -p campus -a 1-38 analysis/behavior_dap/extract_factors.sbatch
#SBATCH -J dap_bx
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 12:00:00
#SBATCH -o logs/dap_bx.%A.%a.log
set -euo pipefail
D=/u/project/pellegrini/gkislik/dogs
W=$D/dap_behavior
BIN=$D/envs/genomics/bin
DA=$D/reference_json/darwins_ark
T="human_sociability toy_motor_patterns biddability dog_sociability environmental_engagement proximity_seeking arousal_level agonistic_threshold"
c=chr${SLURM_ARRAY_TASK_ID:?}
mkdir -p $W/parts
[[ -s $W/parts/$c.npz ]] && { echo "$c exists"; exit 0; }
for t in $T; do zcat $DA/wts_$t.tsv.gz | awk -v c=$c '$1==c {print $1"\t"$2}'; done | sort -k2,2n -u > $W/parts/sites_$c.tsv
$BIN/bcftools query -r $c -T $W/parts/sites_$c.tsv -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' $D/dap_panel/dap.vcf.gz \
  | $BIN/python3 $D/analysis/behavior_dap/extract_factors.py $c $DA $W/parts/$c.tmp.npz $T
mv $W/parts/$c.tmp.npz $W/parts/$c.npz
echo "DAPBX-DONE $c"
