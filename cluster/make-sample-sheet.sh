#!/usr/bin/env bash
# Build a sample sheet (and per-sample FASTQ directories) from a flat directory
# of DOGS-Gen-<N>-* files.
#
#   bash cluster/make-sample-sheet.sh <fastq-source-dir> [ages.tsv] > sample_sheet.gen.tsv
#
# Stage 1 globs "$FASTQ_DIR"/*_R1_*.fastq.gz, so every dog needs its own
# directory or one sample would swallow all 96. Rather than copy ~400GB, this
# creates $D/fastq/<sample>/ containing symlinks to the originals.
#
# Optional ages.tsv: two columns, <dog-number><TAB><age-years>. Dogs missing
# from it get age 0, which the pipeline treats as unknown for the microbiome age
# comparison. Real ages matter if these samples will retrain the age model.
#
# Deliberately avoids bash 4 features (associative arrays, ${v,,}, mapfile) and
# GNU-only find flags, so it behaves the same on macOS bash 3.2 and on Hoffman.
set -euo pipefail

SRC="${1:?usage: make-sample-sheet.sh <fastq-source-dir> [ages.tsv]}"
AGES="${2:-}"
[ -d "$SRC" ] || { echo "ERROR: not a directory: $SRC" >&2; exit 1; }
[ -n "$AGES" ] && [ ! -f "$AGES" ] && { echo "ERROR: ages file not found: $AGES" >&2; exit 1; }

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
DOGS_SITE="${DOGS_SITE:-hoffman}"
# shellcheck source=/dev/null
. "$PIPELINE_DIR/site/${DOGS_SITE}.sh"

age_of() {  # dog number -> age, or 0 when absent
  [ -n "$AGES" ] || { echo 0; return; }
  awk -F'\t' -v n="$1" '$1==n {print $2; found=1} END{if(!found) print 0}' "$AGES" | head -1
}

# Dog number is the integer after DOGS-Gen- and before the next hyphen.
NUMS=$(find "$SRC" -maxdepth 1 -name 'DOGS-Gen-*_R[12]_*.fastq.gz' 2>/dev/null \
        | sed 's|.*/||' | sed -n 's/^DOGS-Gen-\([0-9][0-9]*\)-.*/\1/p' | sort -n -u)
[ -n "$NUMS" ] || { echo "ERROR: no DOGS-Gen-<N>-*_R[12]_*.fastq.gz under $SRC" >&2; exit 1; }
echo "Found $(echo "$NUMS" | wc -l | tr -d ' ') dogs in $SRC" >&2

printf 'sample_id\tfastq_dir\tage\toutput_name\twork_dir\tpub_dir\tfrom_stage\tsex\tnotes\n'

kept=0; skipped=0
for n in $NUMS; do
  sample="DOGS-Gen-$n"
  lower=$(echo "$sample" | tr '[:upper:]' '[:lower:]')
  dir="$D/fastq/$sample"

  n_r1=$(find "$SRC" -maxdepth 1 -name "DOGS-Gen-$n-*_R1_*.fastq.gz" | wc -l | tr -d ' ')
  n_r2=$(find "$SRC" -maxdepth 1 -name "DOGS-Gen-$n-*_R2_*.fastq.gz" | wc -l | tr -d ' ')

  # An unpaired or empty dog would fail at Stage 1 an hour in; leave it out of
  # the sheet and report it rather than letting the array discover it.
  if [ "$n_r1" -eq 0 ] || [ "$n_r1" -ne "$n_r2" ]; then
    echo "  SKIP $sample: $n_r1 R1 vs $n_r2 R2 files" >&2
    skipped=$((skipped + 1))
    continue
  fi

  mkdir -p "$dir"
  find "$SRC" -maxdepth 1 -name "DOGS-Gen-$n-*_R[12]_*.fastq.gz" | while read -r f; do
    ln -sf "$f" "$dir/$(basename "$f")"
  done

  printf '%s\t%s\t%s\t%s\t%s\t%s\t1\t\t%s lane-pairs\n' \
    "$sample" "$dir" "$(age_of "$n")" "$sample" \
    "$D/work/$sample/analysis" "$D/results/$lower" "$n_r1"
  kept=$((kept + 1))
done

echo "Wrote sheet for $kept dogs (skipped $skipped)" >&2
echo "Submit with:  qsub -t 2-$((kept + 1)) -tc 40 cluster/submit-array.sh <sheet>" >&2
