#!/bin/bash
# Imputation-panel v3 pilot: impute one test dog on one chromosome at one depth with one pilot panel
# (pilot_build.sbatch). Rows of $T/pilot/jobs.tsv (or $JOBS): chrom sample depth panel
#   sample COSMO: COSMO3 (2.2x) subsampled, or COSMO3 itself at 2.2x (truth = COSMO-SUPER, as impute.sbatch)
#   sample SRR*/SAMN*/SAMEA*: $T/sra/<id>/lowpass.bam subsampled to the depth (truth = its full-depth gVCF)
#   sbatch -A pi_matteop -p campus -a 1-N%60 -c 8 --mem=24G -t 12:00:00 --constraint="$(cat cluster/avx2_constraint.txt)" analysis/imputation_panel/pilot_impute.sbatch
#SBATCH -J pilot_imp
#SBATCH -o logs/pilot_imp.%A.%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs
B=$D/envs/genomics/bin; G=$D/envs/glimpse/bin
T=$D/scratch/panel_test
read -r c sample depth panel < <(sed -n "${SLURM_ARRAY_TASK_ID:?}p" ${JOBS:-$T/pilot/jobs.tsv})
grep -q avx2 /proc/cpuinfo || { echo "no avx2 on $(hostname)"; exit 1; }
out=$T/pilot/imp/${c}_${sample}_${depth}_${panel}.bcf
[[ -s $out.csi ]] && { echo "exists"; exit 0; }
W=$T/tmp/pilot_${c}_${sample}_${depth}_${panel}; rm -rf $W; mkdir -p $W $T/pilot/imp
# depth of a BAM over the autosomes = mapped reads x mean read length / autosome length
bamdepth() { local rl; rl=$($B/samtools view "$1" chr1:20000000-21000000 | head -5000 | awk '{s+=length($10)} END{print (NR ? s/NR : 100)}')
  $B/samtools idxstats "$1" | awk -v rl="$rl" '$1 ~ /^chr[0-9]+$/ {r+=$3; L+=$2} END{printf "%.4f", r*rl/L}'; }
if [[ $sample == COSMO ]]; then
  S=$D/work/COSMO3/analysis/markdup.bam; full=2.2
else
  S=$T/sra/$sample/lowpass.bam; full=$(bamdepth $S)
fi
frac=$(awk -v d=${depth%x} -v f=$full 'BEGIN{x=d/f; if(x>=1) x=0.9999; printf "%.4f", x}'); frac=${frac#0}
if [[ $frac == .9999 && $sample == COSMO ]]; then $B/samtools view -b -o $W/in.bam $S $c
# seed 11: the SRA lowpass BAMs were themselves drawn with seed 7, and the same seed keeps every read again
elif [[ $sample == COSMO ]]; then $B/samtools view -b -s 7$frac -o $W/in.bam $S $c
else $B/samtools view -b -s 11$frac -o $W/in.bam $S $c; fi
$B/samtools index $W/in.bam
ref=$T/ref/${panel}_$c.bcf
while IFS=$'\t' read -r id chrom ireg oreg rest; do
  echo "$G/GLIMPSE2_phase --bam-file $W/in.bam --reference $ref --fasta $D/canFam4.fa --input-region $ireg --output-region $oreg --threads 1 --output $W/chunk$id.bcf > $W/chunk$id.log 2>&1 && $B/bcftools index -f $W/chunk$id.bcf"
done < $D/COSMO/glimpse2_dog10k/chunks/$c.txt | xargs -P ${SLURM_CPUS_PER_TASK:-8} -I{} bash -c '{}'
n=$(wc -l < $D/COSMO/glimpse2_dog10k/chunks/$c.txt)
ls $W/chunk*.bcf | sort -V > $W/list.txt
[[ $(wc -l < $W/list.txt) == $n ]] || { echo "only $(wc -l < $W/list.txt)/$n chunks"; tail -2 $W/chunk*.log; exit 1; }
$G/GLIMPSE2_ligate --input $W/list.txt --output $out > $W/ligate.log 2>&1
$B/bcftools index -f $out
rm -rf $W
echo "IMP-DONE $c $sample $depth $panel (source ${full}x, frac $frac)"
