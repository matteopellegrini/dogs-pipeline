#!/bin/bash
# F_ROH depth calibration: 19 distinct dogs sequenced at >=4x (incl. Cosmo 17.6x), each downsampled to
# 0.1/0.25/0.5/1/2x and imputed with the v2 panel (stage 7 only). Full-depth BCFs come from the cohort re-impute
# archive (dogs_imputed_v2/<name>.bcf). Keeps only the imputed BCF: scratch/froh_calib/<dog>/<depth>.bcf
#   sbatch -A pi_mcdb -p pi_shared --qos=pi_shared24 -a 1-95%60 -c 8 --mem=48G -t 20:00:00 \
#          --constraint="$(cat cluster/avx2_constraint.txt)" analysis/imputation_panel/froh_calib.sbatch
#SBATCH -J froh_calib
#SBATCH -o logs/froh_calib.%A.%a.log
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 DOGS_SITE=hoffman
D=/u/project/pellegrini/gkislik/dogs; cd $D
B=$D/envs/genomics/bin
LEVELS=(0.1 0.25 0.5 1 2)
i=$(( ${SLURM_ARRAY_TASK_ID:?} - 1 ))
row=$(sed -n "$(( i / 5 + 1 ))p" scratch/froh_calib_dogs.tsv); [ -n "$row" ] || exit 0
dog=$(echo "$row" | cut -f1); bam=$(echo "$row" | cut -f2); depth=$(echo "$row" | cut -f3)
lev=${LEVELS[$(( i % 5 ))]}
grep -q avx2 /proc/cpuinfo || { echo "no avx2"; exit 1; }
OUT=$D/scratch/froh_calib/$dog; mkdir -p $OUT
[ -s $OUT/$lev.bcf ] && { echo "CALIB-SKIP $dog $lev"; exit 0; }
W=$D/scratch/froh_calib/work/$dog-$lev; rm -rf $W; mkdir -p $W/analysis $W/pub
frac=$(awk -v d=$lev -v t=$depth 'BEGIN{printf "%.4f", d/t}'); frac=${frac#0}
$B/samtools view -b -s 11$frac -@ 8 -o $W/analysis/markdup.bam $bam && $B/samtools index -c $W/analysis/markdup.bam || { echo "CALIB-FAIL $dog $lev subsample"; exit 1; }
sheet=$W/sheet.tsv
printf 'sample_id\tfastq_dir\tage\toutput_name\twork_dir\tpub_dir\tfrom_stage\tsex\tnotes\n' > $sheet
printf '%s\t%s\t\t%s\t%s\t%s\t7\t\tF_ROH depth calibration\n' "$dog-$lev" "$D/fastq/none" "$dog-$lev" "$W/analysis" "$W/pub" >> $sheet
snap=$(mktemp "$D/.rdp_snapshot_frohc.XXXXXX.sh"); cp run_dog_pipeline.sh "$snap"
t0=$(date +%s)
TO_STAGE=7 PUBLISH_RESULTS=0 GLIMPSE_PARALLEL=8 bash "$snap" "$sheet" 2 7; rc=$?
rm -f "$snap"
bcf=$(compgen -G "$W/analysis/glimpse2/*_imputed_dog10k.bcf" | head -1)
if [ $rc -eq 0 ] && [ -n "$bcf" ]; then
  cp $bcf $OUT/$lev.bcf.tmp && cp $bcf.csi $OUT/$lev.bcf.csi && mv $OUT/$lev.bcf.tmp $OUT/$lev.bcf
  rm -rf $W
  echo "CALIB-DONE $dog $lev $(( $(date +%s) - t0 ))s"
else
  echo "CALIB-FAIL $dog $lev rc=$rc"; exit 1
fi
