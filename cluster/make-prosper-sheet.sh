#!/usr/bin/env bash
# Build per-sample FASTQ symlink dirs + a pipeline sample sheet for the
# ProsperKits archive (~1,100 customer dogs at 0.1-1x).
#
#   bash cluster/make-prosper-sheet.sh
#
# Source layout: /u/project/pellegrini_archive/data/ProsperKits is FLAT —
# every dog's files sit in one directory, with several naming generations:
#   DT<kit>_R1.fastq.gz                 (no lane token — most files)
#   DT<kit>_S1_L001_R1.fastq.gz
#   DT<kit>-TR_R1.fastq.gz / _TR_...    (top-up resequencing of the same kit)
# Stage 1 requires *_R1_*.fastq.gz (trailing underscore) inside a per-sample
# dir, so each kit gets a dir of renamed symlinks; top-ups land in the same
# dir and stage 1's merge handles them.
#
# Ages come from prosperKitAgeInfo (customer-entered; 'NO AGE' rows skipped).
# Weights in that file are NOT used by the pipeline — they are the held-out
# validation set for the weight predictor.
set -euo pipefail

P=/u/project/pellegrini_archive/data/ProsperKits
D=/u/project/pellegrini/$USER/dogs
SHEET=$D/sample_sheet.prosper.tsv
FQ=$D/fastq_prosper

mkdir -p "$FQ"
printf 'sample_id\tfastq_dir\tage\toutput_name\twork_dir\tpub_dir\tfrom_stage\tsex\tnotes\n' > "$SHEET"

# kit -> age (first numeric age wins)
declare -A AGE
while IFS=$'\t' read -r kit age _weight _file; do
  [[ "$kit" =~ ^[0-9]+$ ]] || continue
  [[ "$age" =~ ^[0-9.]+$ ]] || continue
  [[ -z "${AGE[$kit]:-}" ]] && AGE[$kit]=$age
done < "$P/prosperKitAgeInfo"

made=0 skipped=0
while read -r kit r1file; do
  [[ -n "$kit" ]] || continue
  prefix="${r1file%%_R1*}"                 # e.g. DT31230712104028 or ...-TR
  base="${prefix%-TR}"; base="${base%_TR}" # collapse top-ups onto the base kit
  s="pk-${kit}"
  dir="$FQ/$s"
  [[ -d "$dir" ]] && continue              # manifest can list a kit twice
  # every R1 for this kit, any naming generation, plus its R2 twin
  mapfile -t r1s < <(ls "$P/${base}"*_R1*.fastq.gz "$P/${base}"*_R1.fastq.gz 2>/dev/null | sort -u)
  if (( ${#r1s[@]} == 0 )); then skipped=$((skipped+1)); continue; fi
  ok=1
  for r1 in "${r1s[@]}"; do
    r2="${r1/_R1/_R2}"
    [[ -e "$r2" ]] || { ok=0; break; }
  done
  (( ok )) || { skipped=$((skipped+1)); continue; }
  mkdir -p "$dir"
  i=0
  for r1 in "${r1s[@]}"; do
    i=$((i+1))
    ln -sf "$r1" "$dir/${s}-L${i}_S1_L$(printf '%03d' $i)_R1_001.fastq.gz"
    ln -sf "${r1/_R1/_R2}" "$dir/${s}-L${i}_S1_L$(printf '%03d' $i)_R2_001.fastq.gz"
  done
  printf '%s\t%s\t%s\t%s\t%s\t%s\t1\t\t\n' \
    "$s" "$dir" "${AGE[$kit]:-}" "$s" "$D/work_prosper/$s/analysis" "$D/results_prosper/$s" >> "$SHEET"
  made=$((made+1))
done < <(awk 'NR>0 {print $1, $2}' "$P/1k_input.txt")

echo "sheet: $SHEET — $made samples ($skipped skipped for missing files)"
echo "submit:  qsub -t 2-$((made+1)) -tc 40 cluster/submit-array.sh sample_sheet.prosper.tsv"
