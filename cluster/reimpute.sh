#!/bin/bash
# Re-impute one published kit against the imputation reference v2 (Dog10K + Dog Aging Project, 2026-11-11) and
# rerun every stage that reads the imputed genotypes (7-13, incl. 13c panel v3 calls), into a STAGING copy of the
# kit's results. Nothing published is touched: compare with analysis/imputation_panel/compare_staging.py, then
# promote + republish.
#
#   sbatch -A pi_matteop -p campus -a 1-N%40 -c 8 --mem=48G -t 10:00:00 \
#          --constraint="$(cat cluster/avx2_constraint.txt)" -J reimpute cluster/reimpute.sh <list>
#
# Rows: name, barcode, pub dir, work dir (the cohort list format, e.g. dap_behavior/kits.tsv deduplicated by pub dir).
# BAM: the archived markdup BAM (pellegrini_archive/data/dogs_bams/<name>.markdup.bam), else the work dir's.
# Outputs: results_imp2/<pub basename>/ (staging JSONs), the imputed BCF kept in
# pellegrini_archive/data/dogs_imputed_v2/<name>.bcf for later reruns.
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONIOENCODING=utf-8 DOGS_SITE=hoffman
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"; cd "$D"
LIST="${1:?usage: sbatch ... cluster/reimpute.sh <list>}"
A=/u/project/pellegrini_archive/data/dogs_bams
KEEP=/u/project/pellegrini_archive/data/dogs_imputed_v2
row=$(sed -n "${SLURM_ARRAY_TASK_ID:?}p" "$LIST"); [ -n "$row" ] || exit 0
name=$(echo "$row" | cut -f1); PUB=$(echo "$row" | cut -f3); WORK0=$(echo "$row" | cut -f4)
grep -q avx2 /proc/cpuinfo || { echo "REIMP-FAIL $name: no avx2 on $(hostname)"; exit 1; }
base=$(basename "$PUB")
STAGE=$D/results_imp2/$base
[ -s "$STAGE/.done" ] && { echo "REIMP-SKIP $name (done)"; exit 0; }
bam=""
for cand in "$A/$name.markdup.bam" "$A/$base.markdup.bam" "$A/${base#pk-}.markdup.bam" "$WORK0/markdup.bam"; do
  [ -s "$cand" ] && { bam=$cand; break; }
done
[ -n "$bam" ] || { echo "REIMP-SKIP $name (no BAM)"; exit 0; }
idx=""; for i in "$bam.csi" "$bam.bai" "${bam%.bam}.bai"; do [ -s "$i" ] && { idx=$i; break; }; done
[ -n "$idx" ] || { echo "REIMP-FAIL $name: BAM has no index"; exit 1; }

# Staging copy of the published results (report JSONs only, not the pre-<tag> backups).
mkdir -p "$STAGE"
find "$PUB" -maxdepth 1 -type f ! -name '*.pre-*' -exec cp -p {} "$STAGE/" \;
WORK=$D/work_imp2/$base/analysis
rm -rf "$WORK"; mkdir -p "$WORK"
ln -s "$bam" "$WORK/markdup.bam"
case $idx in *.csi) ln -s "$idx" "$WORK/markdup.bam.csi";; *) ln -s "$idx" "$WORK/markdup.bam.bai";; esac

sheet=$(mktemp "$D/logs/sheet_imp2.XXXXXX.tsv")
printf 'sample_id\tfastq_dir\tage\toutput_name\twork_dir\tpub_dir\tfrom_stage\tsex\tnotes\n' > "$sheet"
printf '%s\t%s\t\t%s\t%s\t%s\t7\t\tre-impute with Dog10K+DAP panel\n' "$name" "$D/fastq/none" "$name" "$WORK" "$STAGE" >> "$sheet"
# Snapshot in the repo root (the pipeline resolves its paths from its own location; edits mid-task poison bash).
snap=$(mktemp "$D/.rdp_snapshot_imp2.XXXXXX.sh"); cp run_dog_pipeline.sh "$snap"
t0=$(date +%s)
PANEL_V3=1 TO_STAGE=14 PUBLISH_RESULTS=0 GLIMPSE_PARALLEL=${SLURM_CPUS_PER_TASK:-8} bash "$snap" "$sheet" 2 7; rc=$?
rm -f "$sheet" "$snap"
bcf=$(compgen -G "$WORK/glimpse2/*_imputed_dog10k.bcf" | head -1)
if [ $rc -eq 0 ] && [ -n "$bcf" ]; then
  mkdir -p "$KEEP"
  cp "$bcf" "$KEEP/$name.bcf.tmp" && cp "$bcf.csi" "$KEEP/$name.bcf.csi" && mv "$KEEP/$name.bcf.tmp" "$KEEP/$name.bcf"
  echo "$name $(date -Is) $(( $(date +%s) - t0 ))s" > "$STAGE/.done"
  rm -rf "$D/work_imp2/$base"
  echo "REIMP-DONE $name $(( $(date +%s) - t0 ))s"
else
  echo "REIMP-FAIL $name rc=$rc"
  exit 1
fi
