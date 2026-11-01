#!/bin/bash
# Rerun stage 13 (coat colour + 13b single-gene traits) on published kits from a rerun list (name, barcode, pub dir,
# work dir — the same lists as rerun-8.sh / rerun-9-panel-v2.sh). Keeps coat_color.pre-<tag>.json and
# traits_result.pre-<tag>.json once per tag. Nothing is uploaded; republish with cluster/republish_files.sh.
#   qsub -t 1-N -tc 80 cluster/rerun-13-list.sh <list> <tag>
#   Slurm: sbatch -J rerun13 -a 1-N%20 -p campus -c 1 --mem=4G -t 1:00:00 cluster/rerun-13-list.sh <list> <tag>
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=4G,h_rt=1:00:00
#$ -pe shared 1
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
if [ -n "${SLURM_JOB_ID:-}" ]; then SGE_TASK_ID=$SLURM_ARRAY_TASK_ID; SGE_O_WORKDIR=$SLURM_SUBMIT_DIR; fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONIOENCODING=utf-8
D="${SGE_O_WORKDIR:?run via qsub from the repo root}"
LIST="${1:?usage: qsub -t 1-N cluster/rerun-13-list.sh <list> <tag>}"
TAG="${2:?tag required}"
cd "$D"
row=$(sed -n "${SGE_TASK_ID}p" "$LIST"); [ -n "$row" ] || exit 0
pk=$(echo "$row" | cut -f1); PUB=$(echo "$row" | cut -f3); WORK=$(echo "$row" | cut -f4)
bcf=$(compgen -G "$WORK/glimpse2/*_imputed_dog10k.bcf" | head -1)
[ -n "$bcf" ] || { echo "RERUN13-SKIP $pk (no bcf)"; exit 0; }
name=$(basename "$bcf" _imputed_dog10k.bcf)       # the pipeline finds the BCF by sample name, case-sensitively
for f in coat_color traits_result; do
  [ -f "$PUB/$f.pre-$TAG.json" ] || cp "$PUB/$f.json" "$PUB/$f.pre-$TAG.json" 2>/dev/null || true
done
sheet=$(mktemp "$D/logs/sheet13.XXXXXX.tsv")
printf 'sample_id\tfastq_dir\tage\toutput_name\twork_dir\tpub_dir\tfrom_stage\tsex\tnotes\n' > "$sheet"
printf '%s\t%s\t\t%s\t%s\t%s\t13\t\trerun stage 13 (%s)\n' "$name" "$D/fastq/none" "$name" "$WORK" "$PUB" "$TAG" >> "$sheet"
snap=$(mktemp "$D/.rdp_snapshot13.XXXXXX.sh"); cp run_dog_pipeline.sh "$snap"
[ -s "$snap" ] || { echo "RERUN13-DONE $pk rc=4 (empty snapshot)"; rm -f "$sheet" "$snap"; exit 0; }
TO_STAGE=13 PUBLISH_RESULTS=0 bash "$snap" "$sheet" 2 13; rc=$?
rm -f "$sheet" "$snap"
echo "RERUN13-DONE $pk rc=$rc"
