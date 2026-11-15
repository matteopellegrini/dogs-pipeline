#!/bin/bash
# Per-chromosome site lists of the Dog Aging Project VCF and the Dog10K imputation
# panel, for the panel-overlap comparison (summarise with overlap.py).
#   sbatch -A pi_matteop -p campus -a 1-39 analysis/imputation_panel/site_lists.sbatch
#SBATCH -J panel_sites
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 6:00:00
#SBATCH -o logs/panel_sites.%A.%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs
B=$D/envs/genomics/bin
O=$D/scratch/panel_overlap
mkdir -p $O
i=${SLURM_ARRAY_TASK_ID:?}
c=chr$i; [[ $i == 39 ]] && c=chrX
[[ -s $O/dap_$c.tsv.gz && -s $O/d10k_$c.tsv.gz ]] && { echo "$c exists"; exit 0; }
$B/bcftools query -r $c -f '%POS\t%REF\t%ALT\t%INFO/INFO\t%INFO/RAF\t%INFO/AF\n' $D/dap_panel/dap.vcf.gz 2>/dev/null | gzip -1 > $O/dap_$c.tmp.gz
$B/bcftools query -r $c -f '%POS\t%REF\t%ALT\t%AC\t%AN\n' $D/dog10k_panel/AutoAndXPAR.Dog10K.phased_plus_disease_rh.bcf 2>/dev/null | gzip -1 > $O/d10k_$c.tmp.gz
mv $O/dap_$c.tmp.gz $O/dap_$c.tsv.gz; mv $O/d10k_$c.tmp.gz $O/d10k_$c.tsv.gz
echo "SITES-DONE $c"
