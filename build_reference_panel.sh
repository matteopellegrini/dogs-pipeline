#!/usr/bin/env bash
# Build the CNV reference panel from the pooled reference dogs BAM.
#
# Run once (or whenever sample_dogs/ changes) to produce:
#   reference_panel/coverage_1mb.json   — per-1Mb normalised depth ratios
#
# Prerequisites:
#   - sample_dogs/male_panel_markdup.bam (pooled alignment of reference dogs)
#     To regenerate this BAM: merge all FASTQs in sample_dogs/, align with
#     bwa-mem2, sort, markdup (same pipeline as run_dog_pipeline.sh stages 1-4).
#   - canFam4.fa.fai (reference genome index)
#   - samtools on PATH
#
# Usage:
#   bash build_reference_panel.sh

set -euo pipefail

D="$(cd "$(dirname "$0")" && pwd)"
BAM="$D/sample_dogs/male_panel_markdup.bam"
FASTA_FAI="$D/canFam4.fa.fai"
OUT_DIR="$D/reference_panel"
OUT_JSON="$OUT_DIR/coverage_1mb.json"
WINDOWS_BED="$OUT_DIR/windows_1mb.bed"
COVERAGE_TSV="$OUT_DIR/coverage_1mb.tsv"

[[ -f "$BAM" ]]       || { echo "ERROR: BAM not found: $BAM"; exit 1; }
[[ -f "$FASTA_FAI" ]] || { echo "ERROR: FASTA index not found: $FASTA_FAI"; exit 1; }
mkdir -p "$OUT_DIR"

echo "[1/3] Generating 1Mb windows..."
awk 'BEGIN{OFS="\t"} $1~/^chr([0-9]+|X)$/ {
  for(s=0; s<$2; s+=1000000)
    print $1, s, (s+1000000<$2 ? s+1000000 : $2)
}' "$FASTA_FAI" > "$WINDOWS_BED"
echo "  $(wc -l < "$WINDOWS_BED") windows"

echo "[2/3] Running samtools bedcov..."
samtools bedcov "$WINDOWS_BED" "$BAM" > "$COVERAGE_TSV"
echo "  done"

echo "[3/3] Building coverage_1mb.json..."
python3 - << PYEOF
import json, statistics, collections

bed_path = "$COVERAGE_TSV"
out_path = "$OUT_JSON"

def chrom_key(c):
    c2 = c.replace('chr','')
    return (0, int(c2)) if c2.isdigit() else (1, c2)

data = collections.defaultdict(dict)
with open(bed_path) as f:
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) < 4: continue
        chrom, start, end, bases = cols[0], int(cols[1]), int(cols[2]), int(cols[3])
        size = end - start
        if size <= 0: continue
        data[chrom][start // 1_000_000] = round(bases / size, 4)

raw = {}
for chrom in sorted(data, key=chrom_key):
    pts = data[chrom]
    raw[chrom] = [pts.get(i, 0.0) for i in range(max(pts) + 1)]

auto_depths = [d for c, arr in raw.items() if c not in ('chrX', 'chrY') for d in arr if d > 0]
panel_median = statistics.median(auto_depths)
print(f"  Reference panel autosomal median depth: {panel_median:.2f}x")

result = {}
for chrom, arr in raw.items():
    if chrom == 'chrX':
        ratio_arr = [round(d / (panel_median * 0.5), 4) if d > 0 else 0.0 for d in arr]
    else:
        ratio_arr = [round(d / panel_median, 4) if d > 0 else 0.0 for d in arr]
    result[chrom] = {'ratio': ratio_arr}

result['_meta'] = {
    'source': 'reference_panel',
    'bam': '$BAM',
    'panel_median_depth': round(panel_median, 2),
    'n_dogs': 6,
    'dogs': ['DOGS-Gen-2', 'DOGS-Gen-3', 'DOGS-Gen-6', 'DOGS-Gen-9', 'DOGS-Gen-30', 'DOGS-Gen-47'],
    'note': ('Pooled alignment of 6 reference dogs. ratio = per-1Mb depth / autosomal median. '
             'Used in Stage 6 to compute ref_depth_pct for CNV artefact filtering.'),
}

with open(out_path, 'w') as f:
    json.dump(result, f)
print(f"  Written: {out_path} ({len(result)-1} chromosomes)")
PYEOF

echo "Done. Reference panel: $OUT_JSON"
