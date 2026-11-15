#!/bin/bash
# Genome-wide check of the v2 panel on Cosmo: stage 7 only, at the depths of the cosmo_truth runs, into
# scratch/imp2_cosmo/<run>. Score old vs new with cluster/cosmo-truth's score_chrom.py (score_cosmo.sbatch).
#   sbatch -A pi_matteop -p campus -a 1-6 -c 8 --mem=48G -t 10:00:00 --constraint="$(cat cluster/avx2_constraint.txt)" \
#          analysis/imputation_panel/validate_cosmo.sbatch
#SBATCH -J imp2_cosmo
#SBATCH -o logs/imp2_cosmo.%A.%a.log
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 DOGS_SITE=hoffman
D=/u/project/pellegrini/gkislik/dogs; cd $D
B=$D/envs/genomics/bin
run=$(echo COSMO3-DS0.1 COSMO3-DS0.25 COSMO3-DS0.5 COSMO3-DS1.0 COSMO3 admera-0001 | cut -d' ' -f${SLURM_ARRAY_TASK_ID:?})
grep -q avx2 /proc/cpuinfo || { echo "no avx2"; exit 1; }
W=$D/scratch/imp2_cosmo/$run; mkdir -p $W/analysis $W/pub
case $run in
  COSMO3-DS*) d=${run#COSMO3-DS}; frac=$(awk -v d=$d 'BEGIN{printf "%.4f", d/2.2}'); frac=${frac#0}
              [ -s $W/analysis/markdup.bam.csi ] || { $B/samtools view -b -s 7$frac -@ 8 -o $W/analysis/markdup.bam $D/work/COSMO3/analysis/markdup.bam && $B/samtools index -c $W/analysis/markdup.bam; } ;;
  *) ln -sf $(readlink -f $D/work/$run/analysis/markdup.bam) $W/analysis/markdup.bam
     for i in csi bai; do f=$(readlink -f $D/work/$run/analysis/markdup.bam).$i; [ -s $f ] && ln -sf $f $W/analysis/markdup.bam.$i; done ;;
esac
sheet=$W/sheet.tsv
printf 'sample_id\tfastq_dir\tage\toutput_name\twork_dir\tpub_dir\tfrom_stage\tsex\tnotes\n' > $sheet
printf '%s\t%s\t\t%s\t%s\t%s\t7\t\tv2 panel check\n' "$run" "$D/fastq/none" "$run" "$W/analysis" "$W/pub" >> $sheet
snap=$(mktemp "$D/.rdp_snapshot_v2c.XXXXXX.sh"); cp run_dog_pipeline.sh "$snap"
t0=$(date +%s)
TO_STAGE=7 PUBLISH_RESULTS=0 GLIMPSE_PARALLEL=8 bash "$snap" "$sheet" 2 7; rc=$?
rm -f "$snap"
echo "IMP2-COSMO-DONE $run rc=$rc $(( $(date +%s) - t0 ))s"
