#!/bin/bash
# Low-pass reads for a rare-breed truth dog: the first 10M spots of an SRA run (a raw-FASTQ submission is in
# instrument order, so these are a random ~1.2x sample), aligned as run_dog_pipeline.sh stages 3+4.
# Truth comes from the full-depth gVCF already in sra_full/<RUN>.g.vcf.gz (sra-full-genotype.sh).
#   sbatch -A pi_matteop -p campus -a 1-2 analysis/imputation_panel/sra_lowpass.sbatch
#SBATCH -J sra_lowpass
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o logs/sra_lowpass.%A.%a.log
set -euo pipefail
D=/u/project/pellegrini/gkislik/dogs
export PATH=/u/home/g/gkislik/project-pellegrini/sratoolkit.3.1.1-ubuntu64/bin:$D/envs/genomics/bin:$PATH OPENBLAS_NUM_THREADS=1
run=$(echo SRR40570403 SRR40570422 | cut -d' ' -f${SLURM_ARRAY_TASK_ID:?})
O=$D/scratch/panel_test/sra/$run; mkdir -p $O; cd $O
[[ -s lowpass.bam.csi ]] && { echo exists; exit 0; }
fastq-dump -X 10000000 --split-spot -Z $run 2>fqd.log \
  | awk 'NR%4==1{k=$1; sub(/\/[12]$/,"",k); sub(/\.[12]$/,"",k)} {rec[NR%4]=$0}
         NR%4==0{r=rec[1]"\n"rec[2]"\n"rec[3]"\n"rec[0]; if(k==pk){print pr; print r; pk=""} else {pr=r; pk=k}}' \
  | bwa-mem2 mem -p -t 8 -R "@RG\tID:$run\tSM:$run\tPL:ILLUMINA\tLB:WGS" $D/canFam4_plusY.fa - 2>bwa.log \
  | samtools sort -n -m 1G -@ 8 -T tmp1 | samtools fixmate -m -@ 8 - - \
  | samtools sort -m 1G -@ 8 -T tmp2 | samtools markdup -@ 8 --write-index - lowpass.bam
samtools idxstats lowpass.bam | awk '$1 ~ /^chr[0-9]+$/ {n++; if($3>0) k++; t+=$3} END{print "autosomes with reads", k "/" n, "mapped", t}'
samtools idxstats lowpass.bam | awk '$1 ~ /^chr[0-9]+$/ {printf "%s:%.0f ", $1, $3/($2/1e6)} END{print ""}'
echo "LOWPASS-DONE $run"
