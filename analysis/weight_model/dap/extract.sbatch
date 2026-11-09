#!/bin/bash
# Dog Aging Project -> per-chromosome breed-panel dosages + DA size PRS (extract_chr.py).
#   sbatch -A pi_matteop -p campus -a 1-38 analysis/weight_model/dap/extract.sbatch
#SBATCH -J dap_wx
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 8:00:00
#SBATCH -o logs/dap_wx.%A.%a.log
set -euo pipefail
D=/u/project/pellegrini/gkislik/dogs
W=$D/dap_weight
BIN=$D/envs/genomics/bin
c=chr${SLURM_ARRAY_TASK_ID:?}
mkdir -p $W/parts
[[ -s $W/parts/$c.npz ]] && { echo "$c exists"; exit 0; }
{ awk -v c=$c -F'\t' 'NR>1 && $1==c {print $1"\t"$2}' $D/breed_panel_v23/sites.tsv
  zcat $D/reference_json/darwins_ark_size_prs.tsv.gz | awk -v c=$c '$1==c {print $1"\t"$2}'
} | sort -k2,2n -u > $W/parts/sites_$c.tsv
$BIN/bcftools query -r $c -T $W/parts/sites_$c.tsv -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' $D/dap_panel/dap.vcf.gz \
  | $BIN/python3 $D/analysis/weight_model/dap/extract_chr.py $c $D/breed_panel_v23/sites.tsv \
      $D/reference_json/darwins_ark_size_prs.tsv.gz $W/parts/$c.tmp.npz
mv $W/parts/$c.tmp.npz $W/parts/$c.npz
echo "DAPWX-DONE $c"
