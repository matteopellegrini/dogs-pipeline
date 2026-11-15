#!/bin/bash
# Impute Cosmo on one chromosome at one depth with one test panel (build_ref.sbatch output).
# Task list = jobs.tsv rows (chr depth panel); depths are COSMO3 (2.2x) subsamples or the Admera 5.3x run.
#   sbatch -A pi_matteop -p campus -a 1-N --constraint="$(cat cluster/avx2_constraint.txt)" analysis/imputation_panel/impute.sbatch
#SBATCH -J panel_imp
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 12:00:00
#SBATCH -o logs/panel_imp.%A.%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs
B=$D/envs/genomics/bin
G=$D/envs/glimpse/bin
T=$D/scratch/panel_test
read -r c depth panel < <(sed -n "${SLURM_ARRAY_TASK_ID:?}p" $T/jobs.tsv)
grep -q avx2 /proc/cpuinfo || { echo "no avx2 on $(hostname)"; exit 1; }
out=$T/imp/${c}_${depth}_${panel}.bcf
[[ -s $out.csi ]] && { echo "exists"; exit 0; }
W=$T/tmp/${c}_${depth}_${panel}; rm -rf $W; mkdir -p $W $T/imp
case $depth in
  SRR*) run=${depth%%-*}; d=${depth#*-}; S=$T/sra/$run/lowpass.bam
        full=$($B/samtools idxstats $S | awk '$1 ~ /^chr[0-9]+$/ {r+=$3; L+=$2} END{printf "%.4f", r*150/L}')
        frac=$(awk -v d=${d%x} -v f=$full 'BEGIN{x=d/f; if(x>=1) x=0.9999; printf "%.4f", x}'); frac=${frac#0}
        $B/samtools view -b -s 7$frac -o $W/in.bam $S $c ;;
  5.3x) $B/samtools view -b -o $W/in.bam $D/work/admera-0001/analysis/markdup.bam $c ;;
  2.2x) $B/samtools view -b -o $W/in.bam $D/work/COSMO3/analysis/markdup.bam $c ;;
  *)    frac=$(awk -v d=${depth%x} 'BEGIN{printf "%.4f", d/2.2}'); frac=${frac#0}
        $B/samtools view -b -s 7$frac -o $W/in.bam $D/work/COSMO3/analysis/markdup.bam $c ;;
esac
$B/samtools index $W/in.bam
ref=$T/ref/${panel}_$c.bcf
while IFS=$'\t' read -r id chrom ireg oreg rest; do
  echo "$G/GLIMPSE2_phase --bam-file $W/in.bam --reference $ref --fasta $D/canFam4.fa --input-region $ireg --output-region $oreg --threads 1 --output $W/chunk$id.bcf > $W/chunk$id.log 2>&1 && $B/bcftools index -f $W/chunk$id.bcf"
done < $D/COSMO/glimpse2_dog10k/chunks/$c.txt | xargs -P 8 -I{} bash -c '{}'
n=$(wc -l < $D/COSMO/glimpse2_dog10k/chunks/$c.txt)
ls $W/chunk*.bcf | sort -V > $W/list.txt
[[ $(wc -l < $W/list.txt) == $n ]] || { echo "only $(wc -l < $W/list.txt)/$n chunks"; tail -2 $W/chunk*.log; exit 1; }
$G/GLIMPSE2_ligate --input $W/list.txt --output $out > $W/ligate.log 2>&1
$B/bcftools index -f $out
rm -rf $W
echo "IMP-DONE $c $depth $panel"
