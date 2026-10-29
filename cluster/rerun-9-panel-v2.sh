#!/bin/bash
# Re-run stage 9 (breed composition + 9b relatives + 9c chromosome painting +
# 9d Y line) on every published kit with the v2 breed panel (Parker + Dog10K +
# Hayward + Dog Aging Project). No alignment, no imputation: reads the kept
# imputed BCF in the kit's work dir. ~10 min per kit. PANEL/TAG pick another
# panel (v2.1, 2026-10-29: qsub -v PANEL=breed_panel_v21,TAG=v21 ...).
#
#   qsub -t 1-1872 -tc 50 cluster/rerun-9-panel-v2.sh cluster/repass9b_list.txt
#
# Rows: name, barcode, pub dir, work dir (same list as the LAI backfill).
# Nothing is uploaded here; republish from a login node afterwards.
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=8G,h_rt=2:00:00
#$ -pe shared 2
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONIOENCODING=utf-8
D="${SGE_O_WORKDIR:?run via qsub from the repo root}"
LIST="${1:?usage: qsub -t 1-N cluster/rerun-9-panel-v2.sh <list>}"
cd "$D"
PANEL="${PANEL:-breed_panel_v2}"; TAG="${TAG:-v2}"
row=$(sed -n "${SGE_TASK_ID}p" "$LIST"); [ -n "$row" ] || exit 0
pk=$(echo "$row" | cut -f1); PUB=$(echo "$row" | cut -f3); WORK=$(echo "$row" | cut -f4)
ls "$WORK"/glimpse2/*_imputed_dog10k.bcf >/dev/null 2>&1 || { echo "RERUN9-SKIP $pk (no bcf)"; exit 0; }
# keep the pre-<tag> breed result once, so the change is auditable per kit
[ -f "$PUB/breed_result.pre-$TAG.json" ] || cp "$PUB/breed_result.json" "$PUB/breed_result.pre-$TAG.json" 2>/dev/null || true
# one-row sample sheet in the pipeline's column order
sheet=$(mktemp "$D/logs/sheet9.XXXXXX.tsv")
printf 'sample_id\tfastq_dir\tage\toutput_name\twork_dir\tpub_dir\tfrom_stage\tsex\tnotes\n' > "$sheet"
printf '%s\t%s\t\t%s\t%s\t%s\t9\t\trerun stage 9 with breed panel %s\n' "$pk" "$D/fastq/none" "$pk" "$WORK" "$PUB" "$TAG" >> "$sheet"
# Snapshot the pipeline script: bash reads scripts incrementally, an edit mid-task would poison running tasks.
# The snapshot must live IN the repo root: the pipeline resolves PIPELINE_DIR (and site/hoffman.sh) from its own path.
snap=$(mktemp "$D/.rdp_snapshot9.XXXXXX.sh"); cp run_dog_pipeline.sh "$snap"
BREED_PANEL="$D/$PANEL" TO_STAGE=9 PUBLISH_RESULTS=0 bash "$snap" "$sheet" 2 9; rc=$?
rm -f "$sheet" "$snap"
echo "RERUN9-DONE $pk rc=$rc"
