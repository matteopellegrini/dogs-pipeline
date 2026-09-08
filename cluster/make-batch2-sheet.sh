#!/usr/bin/env bash
# Batch 2 sheet builder: the 2026-09-07 prosperkits drop (6.2TB raw).
# Differences from make-prosper-sheet.sh: no manifest (kits enumerated from
# the files themselves), new-drop naming variants (-BR top-ups, _S*_L00*
# multi-lane sets), and only kits with no existing work dir are included.
set -euo pipefail

P=/u/project/pellegrini_archive/data/prosperkits/raw
OLDP=/u/project/pellegrini_archive/data/ProsperKits
D=/u/project/pellegrini/$USER/dogs
SHEET=$D/sample_sheet.batch2.tsv
FQ=$D/fastq_prosper

mkdir -p "$FQ"
printf 'sample_id\tfastq_dir\tage\toutput_name\twork_dir\tpub_dir\tfrom_stage\tsex\tnotes\n' > "$SHEET"

declare -A AGE
while IFS=$'\t' read -r kit age _rest; do
  [[ "$kit" =~ ^[0-9]+$ ]] || continue
  [[ "$age" =~ ^[0-9.]+$ ]] || continue
  [[ -z "${AGE[$kit]:-}" ]] && AGE[$kit]=$age
done < "$OLDP/prosperKitAgeInfo"

made=0 skipped=0 have_age=0
# DW = dog WGS (analyze). DT = targeted BISULFITE sequencing (methylation
# product) — NOT WGS; running it through this pipeline produces garbage
# (GC ~25% vs 42%, fake aneuploidies). 27 DT samples slipped into the
# 2026-09-07 batch this way and were quarantined. DW only, forever.
for kit in $(ls "$P" | grep -oE '^DW[0-9]+' | sed -E 's/^DW//' | sort -u); do
  s="pk-${kit}"
  [[ -d "$D/work_prosper/$s" ]] && continue           # already processed (batch 1 or pilot)
  dir="$FQ/$s"
  [[ -d "$dir" ]] && continue
  mapfile -t r1s < <(ls "$P"/DW"${kit}"*_R1*.fastq.gz 2>/dev/null | sort -u)
  (( ${#r1s[@]} > 0 )) || { skipped=$((skipped+1)); continue; }
  ok=1
  for r1 in "${r1s[@]}"; do
    [[ -e "${r1/_R1/_R2}" ]] || { ok=0; break; }
  done
  (( ok )) || { skipped=$((skipped+1)); echo "SKIP $kit — missing R2 twin" >&2; continue; }
  mkdir -p "$dir"
  i=0
  for r1 in "${r1s[@]}"; do
    i=$((i+1))
    ln -sf "$r1" "$dir/${s}-L${i}_S1_L$(printf '%03d' $i)_R1_001.fastq.gz"
    ln -sf "${r1/_R1/_R2}" "$dir/${s}-L${i}_S1_L$(printf '%03d' $i)_R2_001.fastq.gz"
  done
  [[ -n "${AGE[$kit]:-}" ]] && have_age=$((have_age+1))
  printf '%s\t%s\t%s\t%s\t%s\t%s\t1\t\t\n' \
    "$s" "$dir" "${AGE[$kit]:-}" "$s" "$D/work_prosper/$s/analysis" "$D/results_prosper/$s" >> "$SHEET"
  made=$((made+1))
done

echo "sheet: $SHEET — $made samples ($skipped skipped; $have_age with ages)"
echo "submit:  cd $D && qsub -t 2-$((made+1)) -tc 40 cluster/submit-array.sh sample_sheet.batch2.tsv"
