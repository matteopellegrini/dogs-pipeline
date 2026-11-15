#!/bin/bash
# After build_panel.sbatch: concatenate the per-chromosome Dog10K+DAP panel and build the Dog10K-only window panel
# for sites Dog10K carries but the merged panel lacks (the hand-added disease sites).
#   sbatch -A pi_matteop -p campus -c 4 --mem=16G -t 8:00:00 -o logs/panel_v2_finish.%j.log analysis/imputation_panel/finish_panel.sh
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs
B=$D/envs/genomics/bin
O=$D/ref_panel_v2
cd $O
chrs=$(for i in $(seq 1 38) X; do [[ -s chr/chr$i.bcf.csi ]] && echo chr/chr$i.bcf; done)
echo "chromosomes: $(echo $chrs | wc -w)"
$B/bcftools concat --threads 4 -Ob -o AutoAndXPAR.Dog10K_DAP.tmp.bcf $chrs
$B/bcftools index -f --threads 4 AutoAndXPAR.Dog10K_DAP.tmp.bcf
mv AutoAndXPAR.Dog10K_DAP.tmp.bcf AutoAndXPAR.Dog10K_DAP.bcf; mv AutoAndXPAR.Dog10K_DAP.tmp.bcf.csi AutoAndXPAR.Dog10K_DAP.bcf.csi
echo "panel records: $($B/bcftools index -n AutoAndXPAR.Dog10K_DAP.bcf)  samples: $($B/bcftools query -l AutoAndXPAR.Dog10K_DAP.bcf | wc -l)"
printf 'chrom\tpos\tref\talt\twin_start\twin_end\n' > added_sites.tsv
for f in $chrs; do
  c=$(basename $f .bcf)
  $B/bcftools isec -C -c none -w1 chr/$c.dog10k_clean.bcf chr/$c.bcf 2>/dev/null \
    | $B/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    | awk 'BEGIN{OFS="\t"} {lo=$2-1000000; if(lo<1) lo=1; print $1,$2,$3,$4,lo,$2+1000000}' >> added_sites.tsv
done
echo "added sites: $(($(wc -l < added_sites.tsv) - 1))"
# One window panel holding every window (overlapping windows merged), Dog10K-only clean haplotypes.
tail -n +2 added_sites.tsv | awk 'BEGIN{OFS="\t"} {print $1,$5-1,$6}' | sort -k1,1V -k2,2n \
  | awk 'BEGIN{OFS="\t"} $1==c && $2<=e {if($3>e) e=$3; next} {if(c!="") print c,s,e; c=$1; s=$2; e=$3} END{print c,s,e}' > added_windows.bed
: > added_windows.list
while read c s e; do
  $B/bcftools view -r $c:$((s+1))-$e -Ob -o win_${c}_$s.bcf chr/$c.dog10k_clean.bcf
  $B/bcftools index -f win_${c}_$s.bcf; echo win_${c}_$s.bcf >> added_windows.list
done < added_windows.bed
$B/bcftools concat -a -f added_windows.list -Ob -o added_windows.bcf
$B/bcftools index -f added_windows.bcf
rm -f win_chr*.bcf win_chr*.bcf.csi added_windows.list
echo "window panel records: $($B/bcftools index -n added_windows.bcf)"
: > failed_chunks.txt
echo "FINISH-DONE"
