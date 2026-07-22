#!/usr/bin/env bash
# ============================================================
# run_dog_pipeline.sh  —  Full WGS → Dashboard pipeline
#
# Mode 1 — sample sheet (recommended):
#   bash run_dog_pipeline.sh sample_sheet.tsv [row] [from_stage]
#   Row defaults to 2 (first data row). TSV columns (tab-separated, header required):
#     sample_id  fastq_dir  age  output_name  work_dir  pub_dir  from_stage  sex  notes
#
# Mode 2 — positional args (legacy):
#   bash run_dog_pipeline.sh <DogName> [age] [from_stage] [fastq_dir]
#
# Sample sheet example (sample_sheet.tsv):
#   sample_id  fastq_dir                     age  output_name  work_dir                pub_dir                        from_stage
#   COSMO2     /path/to/COSMO2               3    cosmo2       /path/to/COSMO2/analysis /path/to/public/cosmo2        1
#   Kiki2      /path/to/Kiki                 7    kiki2        /path/to/Kiki2/analysis  /path/to/public/kiki2         1
# ============================================================
set -euo pipefail

D=/Users/matteopellegrini/Downloads/dogs   # base dir for default paths

# ── Parse arguments: sample sheet or legacy positional ───────
if [[ "${1:-}" == *.tsv ]]; then
    SHEET="$1"
    SHEET_ROW="${2:-2}"      # which data row to run (1-indexed including header, so 2 = first sample)
    FROM_STAGE_ARG="${3:-}"

    [[ -f "$SHEET" ]] || { echo "ERROR: sample sheet not found: $SHEET"; exit 1; }

    _row() { awk -F'\t' -v r="$SHEET_ROW" 'NR==r {print $'"$1"'}' "$SHEET"; }
    DOG_NAME=$(_row 1)
    FASTQ_DIR=$(_row 2)
    DOG_ACTUAL_AGE=$(_row 3)
    DOG_LOWER=$(_row 4)
    OUT=$(_row 5)
    PUB=$(_row 6)
    FROM_STAGE="${FROM_STAGE_ARG:-$(_row 7)}"
    FROM_STAGE="${FROM_STAGE:-1}"

    [[ -n "$DOG_NAME" ]]  || { echo "ERROR: empty sample_id in row $SHEET_ROW of $SHEET"; exit 1; }
    [[ -n "$FASTQ_DIR" ]] || { echo "ERROR: empty fastq_dir in row $SHEET_ROW of $SHEET"; exit 1; }
    [[ -n "$OUT" ]]       || OUT="$D/$DOG_NAME/analysis"
    [[ -n "$PUB" ]]       || PUB="$D/dogs-app/public/$DOG_LOWER"
else
    # Legacy positional mode
    DOG_NAME="${1:?Usage: $0 <sample_sheet.tsv> [row] [from_stage]  OR  $0 <DogName> [age] [from_stage] [fastq_dir]}"
    DOG_ACTUAL_AGE="${2:-}"
    FROM_STAGE="${3:-1}"
    FASTQ_SRC="${4:-$DOG_NAME}"
    DOG_LOWER=$(echo "$DOG_NAME" | tr '[:upper:]' '[:lower:]')
    FASTQ_DIR=$D/$FASTQ_SRC
    OUT=$D/$DOG_NAME/analysis
    PUB=$D/dogs-app/public/$DOG_LOWER
fi
REF=$D/canFam4_idx                   # BWA-MEM2 index prefix
FASTA=$D/canFam4.fa
VEP_CACHE=$D/vep_cache

# Shared reference data (same for every dog)
DOG10K_PANEL=$D/dog10k_panel/AutoAndXPAR.Dog10K.phased_plus_disease_rh.bcf
CHUNKS_DIR=$D/COSMO/glimpse2_dog10k/chunks   # reuse existing chunk definitions
OMIA_DB=$D/dogs-app/public/cosmo/omia_result.json  # OMIA variant database
COSMO_PUB=$D/dogs-app/public/cosmo             # reference JSONs to copy
SCOPE_P=$D/COSMO/analysis/cosmo_scope177Phat.txt     # (143933 SNPs × 177 breeds) full Parker panel allele freq matrix
SCOPE_CLUST=$D/COSMO/analysis/scope_clust.txt        # breed ordering for Phat columns
PARKER_BIM=$D/COSMO/analysis/cosmo_parker_full.bim
PARKER_FAM=$D/COSMO/analysis/cosmo_parker_full.fam
PARKER_BED=$D/COSMO/analysis/cosmo_parker_full.bed

SNPEFF_DB="ROS_Cfam_1.115"   # SnpEff database for canFam4 / ROS_Cfam_1.0

# Microbiome
METAPHLAN_BIN="${METAPHLAN_BIN:-/Users/matteopellegrini/Library/Python/3.9/bin/metaphlan}"
MICROBIOME_REF="$D/metagenome/merged_microbiome_age_weight_3.18_final.csv"

# Use env bin dirs directly — avoids micromamba lock contention when multiple
# tools run simultaneously in a pipe (bwa | samtools sort | fixmate | markdup).
ENV_GENOMICS="$HOME/micromamba/envs/genomics"
ENV_GLIMPSE="$HOME/micromamba/envs/glimpse_x86"
MM="env PATH=$ENV_GENOMICS/bin:$PATH LD_LIBRARY_PATH=$ENV_GENOMICS/lib"
MM_GLIMPSE="env PATH=$ENV_GLIMPSE/bin:$PATH LD_LIBRARY_PATH=$ENV_GLIMPSE/lib"
NPROC=8
GLIMPSE_PARALLEL=6

mkdir -p "$OUT" "$PUB"
LOG=$OUT/pipeline.log

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
die() { log "ERROR: $*"; exit 1; }

PIPELINE_START=$(date +%s)
PEAK_MEM_FILE=$(mktemp)
echo 0 > "$PEAK_MEM_FILE"

# Background process: poll RSS of this process tree every 10 s and record peak
(
  while kill -0 $$ 2>/dev/null; do
    # Sum RSS of all processes owned by current user (macOS-compatible)
    rss_kb=$(ps -u "$(id -u)" -o rss= 2>/dev/null | awk 'BEGIN{t=0} {t+=$1} END{print t}')
    cur=$(cat "$PEAK_MEM_FILE")
    if (( rss_kb > cur )); then echo "$rss_kb" > "$PEAK_MEM_FILE"; fi
    sleep 10
  done
) &
MEM_POLL_PID=$!

# Print runtime + peak memory summary on exit (normal or error)
_finish() {
  kill "$MEM_POLL_PID" 2>/dev/null || true
  local end=$(date +%s)
  local elapsed=$(( end - PIPELINE_START ))
  local h=$(( elapsed / 3600 ))
  local m=$(( (elapsed % 3600) / 60 ))
  local s=$(( elapsed % 60 ))
  local peak_kb=$(cat "$PEAK_MEM_FILE" 2>/dev/null || echo 0)
  local peak_mb=$(( peak_kb / 1024 ))
  local peak_gb
  peak_gb=$(awk "BEGIN{printf \"%.1f\", $peak_kb/1048576}")
  rm -f "$PEAK_MEM_FILE"
  log "========================================"
  log " Total runtime : ${h}h ${m}m ${s}s"
  log " Peak memory   : ${peak_mb} MB (${peak_gb} GB)"
  log "========================================"
}
trap _finish EXIT

log "========================================"
log " Pipeline start: $DOG_NAME"
log " FASTQ dir: $FASTQ_DIR"
log " Output:    $OUT"
log " Public:    $PUB"
log "========================================"

# Pre-derive path variables so they're available when skipping early stages
IMPUTED_BCF="$OUT/glimpse2/${DOG_LOWER}_imputed_dog10k.bcf"

if (( FROM_STAGE <= 9 )); then

# ── Stage 1: Merge FASTQ ─────────────────────────────────────
log "=== Stage 1: Merge FASTQ lanes ==="
R1_FILES=$(ls "$FASTQ_DIR"/*_R1_*.fastq.gz 2>/dev/null | sort -V) || die "No R1 FASTQ files in $FASTQ_DIR"
R2_FILES=$(ls "$FASTQ_DIR"/*_R2_*.fastq.gz 2>/dev/null | sort -V) || die "No R2 FASTQ files in $FASTQ_DIR"
log "R1 files: $(echo $R1_FILES | tr ' ' '\n' | wc -l)"
log "R2 files: $(echo $R2_FILES | tr ' ' '\n' | wc -l)"

cat $R1_FILES > "$OUT/merged_R1.fastq.gz"
cat $R2_FILES > "$OUT/merged_R2.fastq.gz"
log "Merged: R1=$(ls -lh $OUT/merged_R1.fastq.gz | awk '{print $5}'), R2=$(ls -lh $OUT/merged_R2.fastq.gz | awk '{print $5}')"

# ── Stage 2: Adapter trimming ─────────────────────────────────
log "=== Stage 2: Adapter trimming (fastp) ==="
$MM fastp \
  --in1  "$OUT/merged_R1.fastq.gz" \
  --in2  "$OUT/merged_R2.fastq.gz" \
  --out1 "$OUT/trimmed_R1.fastq.gz" \
  --out2 "$OUT/trimmed_R2.fastq.gz" \
  --detect_adapter_for_pe \
  --trim_poly_g --trim_poly_x \
  --length_required 36 \
  --qualified_quality_phred 20 \
  --thread $NPROC \
  --json "$OUT/fastp.json" \
  --html "$OUT/fastp.html" \
  2>"$OUT/fastp.log"
rm -f "$OUT/merged_R1.fastq.gz" "$OUT/merged_R2.fastq.gz"
log "Trimming done"

# ── Stages 3+4: Alignment → sort → fixmate → markdup (single pipe) ──
# Piping avoids writing the intermediate SAM (~30-50 GB) and namesorted/
# fixmate BAMs to disk — cuts disk I/O by ~3-4x and wall time by ~40%.
log "=== Stage 3: Alignment (bwa-mem2) ==="
log "=== Stage 4: Sort + markdup (piped) ==="
$MM bwa-mem2 mem \
  -t $NPROC \
  -R "@RG\tID:${DOG_NAME}\tSM:${DOG_NAME}\tPL:ILLUMINA\tLB:WGS" \
  "$REF" \
  "$OUT/trimmed_R1.fastq.gz" \
  "$OUT/trimmed_R2.fastq.gz" \
  2>"$OUT/bwa.log" \
| $MM samtools sort -n -@ $NPROC -T "$OUT/tmp_sort" \
| $MM samtools fixmate -m -@ $NPROC - - \
| $MM samtools sort -@ $NPROC -T "$OUT/tmp_sort2" \
| $MM samtools markdup -@ $NPROC --write-index - "$OUT/markdup.bam"
rm -f "$OUT/trimmed_R1.fastq.gz" "$OUT/trimmed_R2.fastq.gz"
$MM samtools flagstat "$OUT/markdup.bam" | tee -a "$LOG"
log "BAM ready: $OUT/markdup.bam"

# ── Stage 5: Coverage windows (1Mb for karyotype; adaptive for CNV) ─
log "=== Stage 5: Coverage windows ==="

# 5a — 1Mb windows for karyotype
awk 'BEGIN{OFS="\t"} $1~/^chr([0-9]+|X)$/ {
  for(s=0; s<$2; s+=1000000)
    print $1, s, (s+1000000<$2 ? s+1000000 : $2)
}' "$FASTA.fai" > "$OUT/windows_1mb.bed"
$MM samtools bedcov "$OUT/windows_1mb.bed" "$OUT/markdup.bam" > "$OUT/coverage_1mb.tsv"
log "1Mb coverage: $(wc -l < $OUT/coverage_1mb.tsv) windows"

# 5b — estimate mean depth from 1Mb data, then compute adaptive CNV window
# Artifact rejection is handled by the 6-dog reference panel (ref_depth_pct), so windows
# can be small enough to detect known deletions (~50-100kb at 2-6x depth).
# Formula: w = max(15000, 50000 / mean) → 23kb @ 2x, 15kb @ 3x+
# A deletion must occupy ~85% of the window to cross the <15% depth threshold.
CNV_WINDOW=$(awk '
  { bases+=$4; size+=($3-$2) }
  END {
    mean = bases/size
    w = int(50000 / mean + 0.5)
    if (w < 15000)  w = 15000
    if (w > 200000) w = 200000
    print w
  }
' "$OUT/coverage_1mb.tsv")
log "Adaptive CNV window: ${CNV_WINDOW} bp (50000/mean, min 15kb, max 200kb)"

awk -v w="$CNV_WINDOW" 'BEGIN{OFS="\t"} $1~/^chr([0-9]+|X)$/ {
  for(s=0; s<$2; s+=w)
    print $1, s, (s+w<$2 ? s+w : $2)
}' "$FASTA.fai" > "$OUT/windows_cnv.bed"
$MM samtools bedcov "$OUT/windows_cnv.bed" "$OUT/markdup.bam" > "$OUT/coverage_cnv.tsv"
log "CNV coverage: $(wc -l < $OUT/coverage_cnv.tsv) windows"

# ── Stage 6: Coverage + QC + CNV JSON ────────────────────────
log "=== Stage 6: Coverage / QC / CNV JSON ==="
python3 - << PYEOF
import json, collections, statistics, subprocess, re, os

tsv_1mb  = "$OUT/coverage_1mb.tsv"
tsv_cnv  = "$OUT/coverage_cnv.tsv"
pub      = "$PUB"
cnv_win  = int("$CNV_WINDOW")

def chrom_key(c):
    c2 = c.replace('chr','')
    return (0, int(c2)) if c2.isdigit() else (1, c2)

def load_tsv(path):
    rows = []
    with open(path) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 4: continue
            chrom, start, end, bases = cols[0], int(cols[1]), int(cols[2]), int(cols[3])
            size = end - start
            if size <= 0: continue
            rows.append((chrom, start, end, bases / size))
    return rows

# --- coverage_1mb.json (karyotype) ---
import statistics as _stats, json as _json
windows_1mb = load_tsv(tsv_1mb)
data = collections.defaultdict(dict)
for chrom, start, end, depth in windows_1mb:
    data[chrom][start // 1000000] = round(depth, 4)
raw = {}
for chrom in sorted(data, key=chrom_key):
    pts = data[chrom]
    raw[chrom] = [pts.get(i, 0) for i in range(max(pts) + 1)]

# Genome-wide median (autosomal) for ratio normalisation
auto_depths = [d for c, arr in raw.items() if c != 'chrX' for d in arr if d > 0]
kiki_median = _stats.median(auto_depths)

# Sex determination: male = one X (depth ~0.5× autosomal), female = two X (~1.0×)
x_depths = [d for d in raw.get('chrX', []) if d > 0]
x_auto_ratio = _stats.median(x_depths) / kiki_median if x_depths else 1.0
predicted_sex = 'male' if x_auto_ratio < 0.75 else 'female'
# For ratio normalisation of chrX: divide by 0.5 for males so hemizygous X → 1.0
chrx_norm = 0.5 if predicted_sex == 'male' else 1.0
print(f"Sex determination: chrX/auto ratio = {x_auto_ratio:.3f} → {predicted_sex} (chrX norm divisor: {chrx_norm})")

# Load COSMO reference (cosmo/panel reference tracks for tooltip)
cosmo_ref_path = "$COSMO_PUB/coverage_1mb.json"
try:
    with open(cosmo_ref_path) as _f:
        cosmo_ref = _json.load(_f)
except Exception:
    cosmo_ref = {}

result = {}
for chrom, arr in raw.items():
    ref   = cosmo_ref.get(chrom, {})
    cosmo_arr = ref.get('cosmo', []) if isinstance(ref, dict) else []
    panel_arr = ref.get('panel', []) if isinstance(ref, dict) else []
    n = len(arr)
    cosmo_arr = cosmo_arr[:n] + [0.0] * max(0, n - len(cosmo_arr))
    panel_arr = panel_arr[:n] + [0.0] * max(0, n - len(panel_arr))
    if chrom == 'chrX' and predicted_sex == 'male':
        # PAR1 windows are diploid in males (depth ≈ autosomal); non-PAR windows are hemizygous (depth ≈ 0.5×).
        # Use per-window normalization: above 0.75× autosomal → PAR1 → divide by kiki_median;
        # below 0.75× autosomal → non-PAR → divide by kiki_median × 0.5, so both map to ratio 1.0.
        par_thresh = kiki_median * 0.75
        ratio_arr = []
        for d in arr:
            if d == 0:
                ratio_arr.append(0.0)
            elif d >= par_thresh:
                ratio_arr.append(round(d / kiki_median, 4))        # PAR1: diploid norm
            else:
                ratio_arr.append(round(d / (kiki_median * 0.5), 4)) # non-PAR: hemizygous norm
    else:
        ratio_arr = [round(d / kiki_median, 4) if d > 0 else 0.0 for d in arr]
    result[chrom] = {
        'cosmo': [round(v, 4) for v in cosmo_arr],
        'panel': [round(v, 4) for v in panel_arr],
        'ratio': ratio_arr,
    }
result['_meta'] = {'predicted_sex': predicted_sex, 'chrx_auto_ratio': round(x_auto_ratio, 3)}
with open(f'{pub}/coverage_1mb.json', 'w') as f:
    json.dump(result, f)
print(f"coverage_1mb.json: {len(result)-1} chromosomes, median depth {kiki_median:.2f}×, sex={predicted_sex}")

# --- qc_result.json (derived from 1Mb windows) ---
depths = [w[3] for w in windows_1mb]
mean_d   = statistics.mean(depths)
median_d = statistics.median(depths)
std_d    = statistics.stdev(depths)
pct = {t: round(sum(1 for d in depths if d >= t)/len(depths)*100, 1) for t in [10,15,20,30]}
qc_status = 'PASS' if mean_d >= 20 else ('WARN' if mean_d >= 10 else 'FAIL')

chrom_data = collections.defaultdict(list)
for chrom, *_, d in windows_1mb: chrom_data[chrom].append(d)
chroms_out = []
for chrom in sorted(chrom_data, key=chrom_key):
    ds = chrom_data[chrom]
    chroms_out.append({'chrom': chrom, 'mean_depth': round(statistics.mean(ds),1),
        'median_depth': round(statistics.median(ds),1),
        'p10_depth': round(sorted(ds)[len(ds)//10],1),
        'n_bins': len(ds), 'low_bins': sum(1 for d in ds if d < 15)})
# --- read stats from fastp.json + samtools flagstat ---
read_stats = {}
fastp_path = "$OUT/fastp.json"
if os.path.exists(fastp_path):
    with open(fastp_path) as _f:
        fp = json.load(_f)
    bf = fp['summary']['before_filtering']
    af = fp['summary']['after_filtering']
    read_stats['total_reads_raw']              = bf['total_reads']
    read_stats['total_reads_after_qc']         = af['total_reads']
    read_stats['total_bases_raw_gb']           = round(bf['total_bases'] / 1e9, 2)
    read_stats['pct_q30_raw']                  = round(bf['q30_rate'] * 100, 1)
    read_stats['read_length_bp']               = bf.get('read1_mean_length', bf.get('read_mean_length'))
    read_stats['read_length_after_trimming_bp']= af.get('read1_mean_length', af.get('read_mean_length'))
    hist = fp.get('insert_size', {}).get('histogram', [])
    if hist:
        h = hist[1:]  # skip index 0 (undetermined)
        tot = sum(h)
        if tot > 0:
            read_stats['fragment_size_mean_bp'] = round(sum(i * c for i, c in enumerate(h, 1)) / tot)

flagstat = subprocess.run(
    ['samtools', 'flagstat', "$OUT/markdup.bam"],
    capture_output=True, text=True).stdout
for line in flagstat.splitlines():
    if 'primary mapped' in line and 'primary duplicate' not in line:
        m = re.match(r'(\d+)', line)
        if m: read_stats['reads_mapped'] = int(m.group(1))
    if 'primary duplicates' in line:
        m = re.match(r'(\d+)', line)
        if m and read_stats.get('total_reads_after_qc'):
            read_stats['duplication_rate_pct'] = round(
                int(m.group(1)) / read_stats['total_reads_after_qc'] * 100, 1)

qc = {'genome_mean_depth': round(mean_d,1), 'genome_median_depth': round(median_d,1),
    'genome_std_depth': round(std_d,1), 'uniformity_cv': round(std_d/mean_d,3),
    'pct_bins_gt10x': pct[10], 'pct_bins_gt15x': pct[15],
    'pct_bins_gt20x': pct[20], 'pct_bins_gt30x': pct[30],
    'n_low_bins': sum(1 for d in depths if d < 15), 'n_total_bins': len(depths),
    'chromosomes': chroms_out, 'qc_status': qc_status,
    'warning': None if pct[15] >= 95 else f"Only {pct[15]}% of 1Mb bins have ≥15x coverage.",
    'assessment': f"Mean genome coverage {mean_d:.1f}x across {len(depths)} 1Mb bins.",
    'method': 'samtools bedcov over 1Mb bins',
    **read_stats}
with open(f'{pub}/qc_result.json', 'w') as f:
    json.dump(qc, f, indent=2)
print(f"qc_result.json: {qc_status}, mean={mean_d:.1f}x, cnv_window={cnv_win}bp")

# --- cnv_homdel.json (adaptive-window data → structured for CnvTable) ---
windows_cnv = load_tsv(tsv_cnv)
mean_d2 = statistics.mean(d for _,_,_,d in windows_cnv)
hom_del_thresh = mean_d2 * 0.15
raw_dels = [{'chrom': c, 'start': s, 'end': e, 'depth': round(d,2),
             'norm': round(d/mean_d2,3), 'window_bp': cnv_win}
            for c,s,e,d in windows_cnv if d < hom_del_thresh]

# Load gene annotations — read from cosmo reference dir (always present),
# not from the sample pub dir which may not have it yet at Stage 6.
import json as _json, glob as _glob, os as _os
_gene_paths = [f'$COSMO_PUB/cnv_genes.json', f'{pub}/cnv_genes.json']
gene_map = {}
for _gp in _gene_paths:
    try:
        with open(_gp) as _f: gene_map = _json.load(_f); break
    except Exception: pass
try:
    with open(f'{pub}/coverage_1mb.json') as _f: cov_1mb = _json.load(_f)
except Exception: cov_1mb = {}

# Load reference panel coverage (pooled alignment of reference dogs)
ref_panel_path = "$D/reference_panel/coverage_1mb.json"
ref_panel = {}
try:
    with open(ref_panel_path) as _f:
        ref_panel = _json.load(_f)
    n_ref_dogs = ref_panel.get('_meta', {}).get('n_dogs', '?')
    print(f"Panel-of-normals: loaded {ref_panel_path} ({n_ref_dogs} reference dogs)")
except Exception as e:
    print(f"WARNING: could not load reference panel from {ref_panel_path}: {e}")
    print("ref_depth_pct will be None for all regions")

# All genes across all chroms
all_genes = [g for gs in gene_map.values() for g in gs]

# Merge adjacent windows (gap ≤ 1 window) into contiguous regions
raw_dels.sort(key=lambda w: (w['chrom'], w['start']))
merged = []
for w in raw_dels:
    if merged and merged[-1]['chrom'] == w['chrom'] and w['start'] <= merged[-1]['end'] + w['window_bp']:
        r = merged[-1]; r['end'] = max(r['end'], w['end']); r['norms'].append(w['norm'])
    else:
        merged.append({'chrom': w['chrom'], 'start': w['start'], 'end': w['end'], 'norms': [w['norm']]})

regions = []; disrupted_all = {}
for r in merged:
    size_bp = r['end'] - r['start']
    avg_norm = sum(r['norms']) / len(r['norms'])
    # Reference depth from pooled reference panel at overlapping 1Mb bins
    smb, emb = r['start']//1_000_000, r['end']//1_000_000
    ref_chrom_d = ref_panel.get(r['chrom'], {})
    ratio_arr = ref_chrom_d.get('ratio', []) if isinstance(ref_chrom_d, dict) else []
    rvals = [ratio_arr[i] for i in range(smb, min(emb+1, len(ratio_arr))) if ratio_arr[i] > 0]
    ref_depth_pct = round(sum(rvals)/len(rvals)*100) if rvals else None
    # Disrupted genes
    chrom_num = r['chrom'].replace('chr','')
    disrupted = []; disrupted_details = []
    for g in gene_map.get(chrom_num, []):
        if g['end'] < r['start'] or g['start'] > r['end']: continue
        ov = 'full' if g['start'] >= r['start'] and g['end'] <= r['end'] else 'partial'
        exon_ov = 'exonic' if any(es < r['end'] and ee > r['start']
                                  for es, ee in g.get('exons', [])) else 'intronic'
        detail = {'gene': g['name'], 'biotype': g['biotype'],
                  'chrom': r['chrom'], 'start': g['start'], 'end': g['end'],
                  'overlap': ov, 'exon_overlap': exon_ov}
        disrupted.append(g['name'])
        disrupted_details.append(detail)
        if g['name'] not in disrupted_all:
            disrupted_all[g['name']] = detail
    # Only keep ≥50kb regions or gene-disrupting ones
    if size_bp < 50_000 and not disrupted: continue
    size_str = (f"{size_bp/1e6:.2f}Mb" if size_bp >= 1_000_000
                else f"{size_bp//1000}kb" if size_bp >= 1000 else f"{size_bp}bp")
    # Classify: ref_depth_pct < 80 → mappability artefact in reference panel too
    is_artefact = ref_depth_pct is not None and ref_depth_pct < 80
    verdict = 'mappability_artefact' if is_artefact else 'putative_deletion'
    regions.append({'chrom': r['chrom'], 'start': r['start'], 'end': r['end'], 'size': size_str,
                    'sample_pct_mean': round(avg_norm*100), 'ref_depth_pct': ref_depth_pct,
                    'disrupted_genes': disrupted, 'disrupted_gene_details': disrupted_details,
                    'verdict': verdict})

real_regions = [r for r in regions if r['verdict'] != 'mappability_artefact']
artefact_regions = [r for r in regions if r['verdict'] == 'mappability_artefact']
# Only include disrupted genes from real (non-artefact) regions
real_gene_names = {g for r in real_regions for g in r['disrupted_genes']}
real_disrupted = {k: v for k, v in disrupted_all.items() if k in real_gene_names}

win_kb = round(cnv_win/1000, 1)
ref_meta = ref_panel.get('_meta', {})
n_ref_dogs = ref_meta.get('n_dogs', '?')
cnv_out = {
    'regions': real_regions, 'disrupted_genes': list(real_disrupted.values()),
    'artefact_regions': artefact_regions,
    'summary': {
        'total_regions': len(real_regions), 'unique_genes': len(real_disrupted),
        'method': (f'Adaptive CNV window ({win_kb}kb), depth normalised to genome-wide mean. '
                   f'Ratio<0.15 threshold for homozygous deletion. '
                   f'Ref depth from pooled reference panel ({n_ref_dogs} dogs); '
                   f'regions with ref_depth_pct<80% classified as mappability artefacts.'),
        'min_detectable_kb': round(win_kb*2), 'calling_resolution_kb': win_kb,
        'panel_note': (f'{len(real_regions)} putative deletion region(s) after artefact filtering '
                       f'({len(artefact_regions)} artefact region(s) excluded). '
                       f'Ref depth = normalised coverage in {n_ref_dogs}-dog reference panel.'),
        'artefact_note': (f'{len(artefact_regions)} region(s) with ref_depth_pct<80% are canFam4 '
                          f'mappability artefacts (low coverage in the reference panel too).'),
    }
}
with open(f'{pub}/cnv_homdel.json', 'w') as f:
    json.dump(cnv_out, f, indent=2)
print(f"cnv_homdel.json: {len(regions)} regions, {len(disrupted_all)} disrupted genes")
PYEOF

# ── Stage 7: Genotype estimation via GLIMPSE2 ───────────────
#
# The Dog10K panel ($DOG10K_PANEL) is a pre-phased reference panel —
# it is prepared ONCE and shared across all dogs. It is never modified here.
#
# GLIMPSE2_phase reads this dog's low-pass BAM and, for each genomic chunk,
# uses the reference panel's haplotype structure (LD) to estimate genotype
# posteriors (GP) at all 30.4M panel positions. The result is a per-dog BCF
# with genotype calls and posteriors — not a new phased panel.
#
# GLIMPSE2_ligate joins overlapping chunk BCFs into a single per-chromosome
# BCF by selecting the highest-confidence genotype call in the overlap regions.
# bcftools concat then merges all chromosomes into one genome-wide BCF.
log "=== Stage 7: Genotype estimation (GLIMPSE2 × Dog10K reference panel) ==="
GENO_DIR="$OUT/glimpse2/genotyped"
LIGATED_DIR="$OUT/glimpse2/ligated"
IMPUTED_BCF="$OUT/glimpse2/${DOG_LOWER}_imputed_dog10k.bcf"
mkdir -p "$GENO_DIR" "$LIGATED_DIR"

estimate_chunk() {
  local chr="$1" id="$2" ireg="$3" oreg="$4"
  local outfile="$GENO_DIR/${chr}_chunk${id}.bcf"
  [ -f "$outfile" ] && [ -f "${outfile}.csi" ] && { echo "SKIP ${chr}_chunk${id}"; return; }
  # GLIMPSE2_phase: estimates this dog's genotypes at panel positions using BAM reads + LD
  $MM_GLIMPSE GLIMPSE2_phase \
    --bam-file      "$OUT/markdup.bam" \
    --reference     "$DOG10K_PANEL" \
    --fasta         "$FASTA" \
    --input-region  "$ireg" \
    --output-region "$oreg" \
    --threads 2 \
    --output "$outfile" 2>&1 | grep -v "AC/AN INFO fields" | tail -1
  $MM_GLIMPSE bcftools index -f "$outfile"
  echo "DONE ${chr}_chunk${id}"
}
export -f estimate_chunk
export OUT DOG10K_PANEL FASTA GENO_DIR

log "Estimating genotypes across chunks (${GLIMPSE_PARALLEL} parallel jobs)..."
for chunkfile in "$CHUNKS_DIR"/*.txt; do
  chr=$(basename "$chunkfile" .txt)
  while IFS=$'\t' read -r id chrom ireg oreg rest; do
    estimate_chunk "$chr" "$id" "$ireg" "$oreg" &
    while [ "$(jobs -r | wc -l)" -ge "$GLIMPSE_PARALLEL" ]; do sleep 2; done
  done < "$chunkfile"
done
wait
log "Genotype estimation complete"

log "Ligating chunks per chromosome..."
for chr in $(ls "$CHUNKS_DIR"/*.txt | xargs -I{} basename {} .txt); do
  list=$(mktemp)
  ls "$GENO_DIR"/${chr}_chunk*.bcf 2>/dev/null | sort -V > "$list"
  [ -s "$list" ] || { rm "$list"; continue; }
  $MM_GLIMPSE GLIMPSE2_ligate --input "$list" --output "$LIGATED_DIR/${chr}.bcf" 2>&1 | tail -1
  $MM_GLIMPSE bcftools index -f "$LIGATED_DIR/${chr}.bcf"
  rm "$list"
  echo "Ligated $chr"
done

log "Merging chromosomes..."
$MM_GLIMPSE bcftools concat \
  $(ls "$LIGATED_DIR"/chr*.bcf | sort -V) \
  -O b -o "$IMPUTED_BCF"
$MM_GLIMPSE bcftools index -f "$IMPUTED_BCF"
TOTAL=$($MM_GLIMPSE bcftools stats "$IMPUTED_BCF" 2>/dev/null | grep "^SN.*number of records" | awk '{print $NF}')
log "Imputed BCF: $TOTAL variants → $IMPUTED_BCF"

# ── Stage 9: OMIA genotyping from Dog10K imputed panel ───────
log "=== Stage 8: OMIA genotyping (Dog10K imputed + BAM fallback) ==="
python3 - << PYEOF
import subprocess, pysam, json, re

BCF       = "$IMPUTED_BCF"
BAM       = "$OUT/markdup.bam"
OMIA      = "$OMIA_DB"
QC_JSON   = "$PUB/qc_result.json"
PUB       = "$PUB"
DOG       = "$DOG_LOWER"
GP_HIGH   = 0.90   # minimum max(GP) to report a GLIMPSE2 call
MIN_READS_BAM = 10  # minimum depth for direct BAM calling (≥10x threshold)

with open(QC_JSON) as f:
    mean_depth = json.load(f)['genome_mean_depth']
use_bam_fallback = mean_depth >= MIN_READS_BAM
print(f"Mean depth: {mean_depth:.1f}x → BAM fallback {'ENABLED' if use_bam_fallback else 'DISABLED (< 10x)'}")

with open(OMIA) as f:
    omia_ref = json.load(f)

def query_glimpse2(chrom, pos, ref, alt):
    """Return (gt, af, gp_list, in_panel) from imputed BCF."""
    result = subprocess.run(
        ['bcftools', 'query', '-r', f'{chrom}:{pos}-{pos}',
         '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RAF\t[%GT]\t[%GP]\n', BCF],
        capture_output=True, text=True)
    for line in result.stdout.strip().split('\n'):
        if not line: continue
        parts = line.split('\t')
        if len(parts) < 7: continue
        _, p, r, a, af_s, gt, gp_s = parts[:7]
        if int(p) != pos or r != ref or alt not in a.split(','): continue
        try:
            af_parts = af_s.split(',')
            alt_idx = a.split(',').index(alt)
            af = float(af_parts[min(alt_idx, len(af_parts)-1)])
        except Exception:
            af = None
        gp = None
        try:
            gp = [round(float(x), 4) for x in gp_s.split(',')]
        except Exception:
            pass
        return gt, af, gp, True   # in_panel = True
    return None, None, None, False  # not in Dog10K panel

def gt_from_bam(chrom, pos, ref, alt, min_bq=20, min_mq=20):
    counts = {}
    try:
        bam_fh = pysam.AlignmentFile(BAM, 'rb')
        for col in bam_fh.pileup(chrom, pos-1, pos, truncate=True,
                                  min_base_quality=min_bq, min_mapping_quality=min_mq,
                                  ignore_overlaps=True, ignore_orphans=True):
            if col.reference_pos != pos-1: continue
            for r in col.pileups:
                if not r.is_del and not r.is_refskip:
                    b = r.alignment.query_sequence[r.query_position].upper()
                    counts[b] = counts.get(b, 0) + 1
        bam_fh.close()
    except Exception:
        pass
    total = sum(counts.values())
    n_ref = counts.get(ref.upper(), 0)
    n_alt = counts.get(alt.upper(), 0)
    n_v = n_ref + n_alt
    if total < 5 or n_v < 5:
        return None
    f_alt = n_alt / n_v
    zyg = 'ref/ref' if f_alt < 0.1 else ('alt/alt' if f_alt > 0.9 else 'het')
    conf = 'high' if total >= 20 else ('medium' if total >= 10 else 'low')
    return {'zygosity': zyg, 'depth': total, 'ref_count': n_ref, 'alt_count': n_alt,
            'affected': zyg in ('alt/alt', 'het'), 'call_confidence': conf,
            'source': 'bam_direct'}

variants = []
n_panel = 0; n_bam = 0; n_indel = 0; n_not_callable = 0

for v in omia_ref.get('variants', []):
    chrom = v.get('chrom') or ''
    pos   = v.get('pos')
    ref   = v.get('ref') or ''
    alt   = v.get('alt') or ''
    new_v = {k: val for k, val in v.items() if k not in ('cosmo', 'nelk')}

    is_snv = pos and len(ref) == 1 and len(alt) == 1

    if not is_snv:
        # Indels: not resolvable by imputation or low-pass BAM
        new_v[DOG] = {'zygosity': 'indel_no_call', 'affected': False,
                      'call_confidence': 'none', 'source': 'indel'}
        n_indel += 1

    else:
        gt, af, gp, in_panel = query_glimpse2(chrom, int(pos), ref, alt)

        if in_panel:
            # Site is in Dog10K panel — only report if GP is high-confidence
            n_panel += 1
            max_gp = max(gp) if gp else 0.0
            if max_gp >= GP_HIGH:
                alleles = re.split(r'[|/]', gt)
                zyg = ('alt/alt' if set(alleles) == {'1'} else
                       'ref/ref' if set(alleles) == {'0'} else 'het')
                call = {'zygosity': zyg, 'affected': zyg in ('alt/alt', 'het'),
                        'call_confidence': 'high', 'glimpse2_gt': gt,
                        'glimpse2_gp': gp, 'source': 'dog10k_imputed'}
                if af is not None:
                    call['af_dog10k'] = round(af, 4)
            else:
                call = {'zygosity': 'low_gp_no_call', 'affected': False,
                        'call_confidence': 'low', 'glimpse2_gt': gt,
                        'glimpse2_gp': gp, 'source': 'dog10k_imputed',
                        'note': f'max GP {max_gp:.2f} below {GP_HIGH} threshold'}
            new_v[DOG] = call

        elif use_bam_fallback:
            # Not in panel, but depth ≥ 10x: direct BAM call
            bam_call = gt_from_bam(chrom, int(pos), ref, alt)
            if bam_call:
                new_v[DOG] = bam_call
                n_bam += 1
            else:
                new_v[DOG] = {'zygosity': 'no_call', 'affected': False,
                              'call_confidence': 'none',
                              'source': 'not_in_panel_insufficient_reads'}
                n_not_callable += 1
        else:
            # Not in panel, depth < 10x: cannot call reliably
            new_v[DOG] = {'zygosity': 'not_in_panel', 'affected': False,
                          'call_confidence': 'none', 'source': 'not_in_panel_low_coverage',
                          'note': f'Not in Dog10K panel; depth {mean_depth:.1f}x < 10x threshold'}
            n_not_callable += 1

    variants.append(new_v)

def is_snv_v(v): return len(v.get('ref') or '')==1 and len(v.get('alt') or '')==1
affected_snv = sum(1 for v in variants if (v.get(DOG) or {}).get('affected') and is_snv_v(v))
high_conf    = sum(1 for v in variants if (v.get(DOG) or {}).get('affected')
                   and (v.get(DOG) or {}).get('call_confidence') == 'high' and is_snv_v(v))

result = {
    'summary': {
        'total_screened': len(variants),
        'affected_snv': affected_snv,
        'affected_high_confidence': high_conf,
        'indel_unknown': n_indel,
        'unaffected': sum(1 for v in variants
                          if not (v.get(DOG) or {}).get('affected') and is_snv_v(v)),
        'in_dog10k_panel': n_panel,
        'called_from_bam': n_bam,
        'not_callable': n_not_callable,
        'mean_depth': mean_depth,
        'bam_fallback_used': use_bam_fallback,
    },
    'method': (
        f'Primary: GLIMPSE2 Dog10K imputed panel (30.4M SNPs); high-confidence calls require max GP ≥ {GP_HIGH}. '
        f'{"BAM direct call used for SNVs not in Dog10K panel (depth " + str(mean_depth) + "x ≥ 10x)." if use_bam_fallback else "BAM fallback disabled (depth < 10x)."}'
    ),
    'variants': variants,
}
with open(f'{PUB}/omia_result.json', 'w') as f:
    json.dump(result, f, indent=2)
print(f"omia_result.json: {len(variants)} variants | panel={n_panel} bam={n_bam} not_callable={n_not_callable} indels={n_indel}")
print(f"  affected SNVs: {affected_snv} ({high_conf} high confidence)")
PYEOF

# ── Stage 9: Breed prediction (GLIMPSE2 genotypes → supervised SCOPE) ──
#
# Step 1 — Infer genotypes at all ~143k Parker panel sites
#   Query the Dog10K imputed BCF (from Stage 7) at every Parker SNP position.
#   Extract GP-weighted posterior dosages: E[alt copies] = P(het)×1 + P(hom_alt)×2.
#   Handles allele-orientation mismatches between BIM (Parker) and BCF (Dog10K).
#
# Step 2 — Supervised SCOPE ancestry projection
#   Use the dosage vector as the sample's genotype input to SCOPE.
#   NNLS projection onto the 177-breed Parker 2017 reference Q matrix
#   estimates admixture proportions across all breeds.
log "=== Stage 9: Breed prediction (GLIMPSE2 genotypes at Parker sites → supervised SCOPE) ==="
python3 - << PYEOF
import subprocess, tempfile, os, numpy as np, json
from scipy.optimize import nnls

BCF     = "$IMPUTED_BCF"
BIM     = "$PARKER_BIM"
SCOPE_P = "$SCOPE_P"    # (103542 SNPs × 177 breeds) allele freq matrix
CLUST   = "$SCOPE_CLUST"
PUB     = "$PUB"
DOG     = "$DOG_NAME"

# Load Parker SNP positions and build lookup: (chr_with_prefix, pos) → row index
parker_snps = []
pos_index = {}
with open(BIM) as f:
    for line in f:
        p = line.strip().split()
        chrom_num, pos = p[0], int(p[3])
        parker_snps.append({'chrom_num': chrom_num, 'pos': pos, 'a1': p[4], 'a2': p[5]})
        pos_index[(f'chr{chrom_num}', pos)] = len(parker_snps) - 1
n_snps = len(parker_snps)
print(f"Parker panel: {n_snps} SNPs to query")

# Write a BED file for bcftools -R (avoids ARG_MAX on 100k+ positions)
bed_fh = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
for s in parker_snps:
    bed_fh.write(f"chr{s['chrom_num']}\t{s['pos']-1}\t{s['pos']}\n")
bed_fh.close()

# Extract GP posteriors for all Parker SNPs in a single bcftools call
result = subprocess.run(
    ['bcftools', 'query', '-R', bed_fh.name,
     '-f', '%CHROM\t%POS\t%REF\t%ALT\t[%GP]\n', BCF],
    capture_output=True, text=True)
os.unlink(bed_fh.name)

dosages = np.full(n_snps, np.nan)
genotyped = 0
for line in result.stdout.strip().split('\n'):
    if not line: continue
    parts = line.split('\t')
    if len(parts) < 5: continue
    c, pos_s, ref, alt, gp_s = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
    idx = pos_index.get((c, pos_s))
    if idx is None: continue
    s = parker_snps[idx]
    if ref not in (s['a1'], s['a2']) or alt not in (s['a1'], s['a2']): continue
    try:
        gp = [float(x) for x in gp_s.split(',')]
        if len(gp) < 3: continue
        # Dosage = E[copies of a1] to match Phat allele freq orientation
        # GP = [P(hom_ref), P(het), P(hom_alt)]
        if ref == s['a1']:   # REF=a1: E[a1] = 2*P(hom_ref) + P(het)
            dosages[idx] = 2.0 * gp[0] + gp[1]
        else:                # REF=a2: E[a1] = P(het) + 2*P(hom_alt)
            dosages[idx] = gp[1] + 2.0 * gp[2]
        genotyped += 1
    except Exception:
        continue

valid = ~np.isnan(dosages)
pct_covered = 100 * valid.sum() / n_snps
print(f"Imputed dosages: {valid.sum()}/{n_snps} Parker SNPs ({pct_covered:.1f}%)")

# Load SCOPE Phat: (n_snps=103542, K=177) — allele freq per SNP per breed component
P = np.loadtxt(SCOPE_P)   # (103542, 177)
print(f"Phat shape: {P.shape}")

# Breed ordering matches the column ordering in Phat (order of unique breeds in scope_clust.txt)
breed_labels = []
seen = set()
with open(CLUST) as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 3:
            b = parts[2]
            if b not in seen:
                seen.add(b)
                breed_labels.append(b)
print(f"Breed labels: {len(breed_labels)} (first 3: {breed_labels[:3]})")

# NNLS supervised projection: find q s.t. P_valid @ q ≈ dosages_valid, q ≥ 0
P_v = P[valid, :]          # restrict to genotyped sites
x_v = dosages[valid]

# Constrained NNLS with sum-to-1 row appended
A = np.vstack([P_v, np.ones((1, P_v.shape[1]))])
b = np.hstack([x_v, [1.0]])
q_raw, _ = nnls(A, b)
q_total = q_raw.sum()
q = q_raw / (q_total + 1e-12)  # normalize to sum to 1

# Build top breeds by proportion
breed_props = sorted(zip(breed_labels, q.tolist()), key=lambda x: -x[1])
top = [(b, s) for b, s in breed_props if s > 0.001][:20]
print("Top breeds:", [(b, round(s, 4)) for b, s in top[:5]])

# Parker code → human-readable breed name
PARKER_NAMES = {
    'ACKR':'American Cocker Spaniel','AFGH':'Afghan Hound','AIRT':'Airedale Terrier',
    'AKIT':'Akita','AMAL':'Alaskan Malamute','AMST':'American Staffordshire Terrier',
    'AHRT':'American Hairless Terrier','ANAT':'Anatolian Shepherd Dog',
    'AUSC':'Australian Cattle Dog','AUST':'Australian Shepherd',
    'AUST2':'Australian Terrier','BASS':'Basset Hound','BEAG':'Beagle',
    'BEAU':'Beauceron','BELS':'Belgian Malinois','BELG':'Belgian Sheepdog',
    'BELT':'Belgian Tervuren','BERN':'Bernese Mountain Dog',
    'BICH':'Bichon Frise','BLOO':'Bloodhound','BORD':'Border Collie',
    'BORZ':'Borzoi','BOST':'Boston Terrier','BOUX':'Bouvier des Flandres',
    'BOXR':'Boxer','BRIA':'Briard','BRIT':'Brittany',
    'BRUS':'Brussels Griffon','BULL':'Bulldog','BMAS':'Bullmastiff',
    'CAIR':'Cairn Terrier','CANE':'Cane Corso','CAVA':'Cavalier King Charles Spaniel',
    'CHIH':'Chihuahua','CHIN':'Chinese Crested','CHOW':'Chow Chow',
    'CLUM':'Clumber Spaniel','COCK':'Cocker Spaniel','COLL':'Collie',
    'COOK':'Cocker Spaniel','DACH':'Dachshund','DALM':'Dalmatian',
    'DAND':'Dandie Dinmont Terrier','DOBP':'Doberman Pinscher',
    'ECKR':'English Cocker Spaniel','EENGL':'English Foxhound',  # ECKR confirmed = English Cocker Spaniel
    'ESSP':'English Springer Spaniel','ESET':'English Setter',
    'ESKD':'American Eskimo Dog','FBUL':'French Bulldog',
    'FCR':'Flat-Coated Retriever','FINN':'Finnish Spitz',
    'FOXH':'Foxhound','FTRT':'Fox Terrier','GERM':'German Shepherd Dog',
    'GOLD':'Golden Retriever','GORD':'Gordon Setter',
    'GRDN':'Grand Danois','GREY':'Greyhound',
    'GSHP':'German Shorthaired Pointer','GSNAU':'Giant Schnauzer',
    'GWPG':'German Wirehaired Pointer','HAVA':'Havanese',
    'IBIZ':'Ibizan Hound','ICAL':'Icelandic Sheepdog',
    'IRIS':'Irish Setter','IRSW':'Irish Water Spaniel',
    'IRWT':'Irish Wolfhound','ISET':'Irish Setter',
    'ITAL':'Italian Greyhound','JACK':'Jack Russell Terrier',
    'JAP':'Japanese Chin','KEES':'Keeshond','KERRY':'Kerry Blue Terrier',
    'KOMN':'Komondor','KOMO':'Komondor','KUVZ':'Kuvasz',
    'LAB':'Labrador Retriever','LAKE':'Lakeland Terrier',
    'LHAP':'Lhasa Apso','MALT':'Maltese','MAST':'Mastiff',
    'MPIN':'Miniature Pinscher','MSCHN':'Miniature Schnauzer',
    'MPOO':'Miniature Poodle','NFLD':'Newfoundland',
    'NORW':'Norwegian Elkhound','NORB':'Norwich Terrier',
    'NOVA':'Nova Scotia Duck Tolling Retriever','OLDBS':'Old English Sheepdog',
    'OTTO':'Otterhound','PAPI':'Papillon','PEKE':'Pekingese',
    'PHAR':'Pharaoh Hound','PLSK':'Polish Lowland Sheepdog',
    'PNTG':'Pointer','POOD':'Poodle','PORT':'Portuguese Water Dog',
    'PRESA':'Dogo Canario','PUG':'Pug','PULI':'Puli',
    'ROTT':'Rottweiler','SALU':'Saluki','SAMO':'Samoyed',
    'SCHA':'Schapendoes','SCHN':'Schnauzer','SCOT':'Scottish Terrier',
    'SHAR':'Shar-Pei','SHED':'Shetland Sheepdog',
    'SHIB':'Shiba Inu','SHIH':'Shih Tzu','SILK':'Silky Terrier',
    'SLOU':'Sloughi','SMAL':'Small Munsterlander',
    'SOFT':'Soft Coated Wheaten Terrier','SPOO':'Standard Poodle',
    'SSKI':'Swedish Vallhund','STAF':'Staffordshire Bull Terrier',
    'SUSA':'Sussex Spaniel','TPOO':'Toy Poodle','TIBT':'Tibetan Mastiff',
    'TIBS':'Tibetan Spaniel','TIBT2':'Tibetan Terrier',
    'VISZL':'Vizsla','WEIM':'Weimaraner','WELCS':'Welsh Corgi',
    'WELSH':'Welsh Terrier','WEST':'West Highland White Terrier',
    'WHIP':'Whippet','WFOX':'Wire Fox Terrier',
    'WIRE':'Wirehaired Pointing Griffon','XOLO':'Xoloitzcuintli',
    'YORK':'Yorkshire Terrier',
    # Extended Parker panel codes not in original lookup
    'AESK':'American Eskimo Dog','AUCD':'Australian Cattle Dog',
    'AUSS':'Australian Shepherd','AZWK_Mali':'Azawakh (Mali)',
    'BEDT':'Bedlington Terrier','BERD':'Bergamasco Shepherd',
    'BLDH':'Bloodhound','BMAL':'Belgian Malinois',
    'BMD':'Bernese Mountain Dog','BOER':'Boerboel',
    'BORT':'Border Terrier','BOUV':'Bouvier des Flandres',
    'BOX':'Boxer','BPIC':'Berger Picard',
    'BRTR':'Brittany','BSJI':'Basenji',
    'BULD':'Bulldog','BULM':'Bullmastiff','BULT':'Bull Terrier',
    'CANE_Italy':'Cane Corso (Italy)','CARD':'Cardigan Welsh Corgi',
    'CCRT':'Curly-Coated Retriever','CIRN_Italy':'Cirneco dell\'Etna',
    'CKCS':'Cavalier King Charles Spaniel','COTO':'Coton de Tuléar',
    'CPAT_Italy':'Cane Pastore Abruzzese','CRES':'Chinese Crested',
    'DANE':'Great Dane','DDBX':'Dogue de Bordeaux',
    'DEER':'Scottish Deerhound','EURA':'Eurasier',
    'FIEL':'Field Spaniel','FINS':'Finnish Spitz',
    'GDJK':'Grand Basset Griffon Vendéen','GLEN':'Glen of Imaal Terrier',
    'GPYR':'Great Pyrenees','GREE':'Greenland Dog',
    'GSD':'German Shepherd Dog','GSNZ':'Giant Schnauzer',
    'GSMD':'Greater Swiss Mountain Dog','GWHP':'German Wirehaired Pointer',
    'HUSK':'Siberian Husky','ICES':'Icelandic Sheepdog',
    'INCA':'Peruvian Hairless Dog','IRIT':'Irish Terrier',
    'ITGY':'Italian Greyhound','IWOF':'Irish Wolfhound',
    'IWSP':'Irish Water Spaniel','KELP':'Australian Kelpie',
    'KERY':'Kerry Blue Terrier','LEON':'Leonberger',
    'LHSA':'Lhasa Apso','LMUN':'Large Munsterlander',
    'LVMD_Italy':'Levriero Meridionale (Italy)',
    'MAAB_Italy':'Maremma Abruzzese Sheepdog','MBLT':'Miniature Bull Terrier',
    'MNTY':'Montenegrin Mountain Hound','MSNZ':'Miniature Schnauzer',
    'MXOL':'Mexican Hairless Dog','NEAP':'Neapolitan Mastiff',
    'NELK':'Norrbottenspets','NEWF':'Newfoundland',
    'NORF':'Norfolk Terrier','NOWT':'Norwich Terrier',
    'NSDT':'Nova Scotia Duck Tolling Retriever','OES':'Old English Sheepdog',
    'OTTR':'Otterhound','PARS':'Parson Russell Terrier',
    'PBGV':'Petit Basset Griffon Vendéen','PEMB':'Pembroke Welsh Corgi',
    'POM':'Pomeranian','PTWD':'Portuguese Water Dog',
    'PUMI':'Pumi','RATT':'Rat Terrier',
    'REDB':'Redbone Coonhound','RHOD':'Rhodesian Ridgeback',
    'SALU_ArabPen':'Saluki (Arabian Peninsula)',
    'SALU_CentAsia':'Saluki (Central Asia)','SALU_Tribal':'Saluki (Tribal)',
    'SCWT':'Soft Coated Wheaten Terrier','SKIP':'Schipperke',
    'SLOU_NAfrica':'Sloughi (North Africa)','SPIN':'Spinone Italiano',
    'SSHP':'Smooth Collie','SSNZ':'Standard Schnauzer',
    'STBD':'Saint Bernard','SVAL':'Swedish Vallhund',
    'TIBM':'Tibetan Mastiff','TIBM_China':'Tibetan Mastiff (China)',
    'TURV':'Belgian Tervuren','TYFX':'Toy Fox Terrier',
    'VIZS':'Vizsla','VPIN_Italy':'Volpino Italiano',
    'WHPG':'Wirehaired Pointing Griffon','WHWT':'West Highland White Terrier',
    'WOLF-China':'Gray Wolf (China)','WOLF-Croatia':'Gray Wolf (Croatia)',
    'WOLF-India':'Gray Wolf (India)','WOLF-Israel':'Gray Wolf (Israel)',
    'WOLF-Italy':'Gray Wolf (Italy)','WOLF-Portugal':'Gray Wolf (Portugal)',
    'WOLF-Yellowstone':'Gray Wolf (Yellowstone)',
    'XIGO_China':'Xigou (China)',
}

breed_result = {
    # breed_composition: all 177 breeds sorted by proportion (dashboard shows top 6)
    'breed_composition': [{'breed': b,
                           'breed_name': PARKER_NAMES.get(b, b),
                           'proportion': round(s, 6)}
                          for b, s in sorted(zip(breed_labels, q.tolist()), key=lambda x: -x[1])],
    'snps_used': int(valid.sum()),
    'k': P.shape[1],
    'pct_parker_covered': round(pct_covered, 1),
    'reference_panel': 'Parker 2017 (Science) — 143,933 SNPs, 177 breeds',
    'method': ('Supervised SCOPE NNLS projection onto Parker 2017 allele frequency matrix. '
               f'P matrix: {P.shape[0]} SNPs × {P.shape[1]} breeds. '
               'Dosages from GLIMPSE2 Dog10K imputed BCF (posterior GP-weighted, E[a1]).')
}
with open(f'{PUB}/breed_result.json', 'w') as f:
    json.dump(breed_result, f, indent=2)
print("breed_result.json written")
PYEOF

fi # end stages 1–9

if (( FROM_STAGE <= 13 )); then

# ── Stage 10: Functional annotation (SnpEff) ────────────────
log "=== Stage 10: SnpEff annotation ==="
ANN_DIR="$OUT/snpeff"
mkdir -p "$ANN_DIR"

# snpEff wrapper (conda python script) needs python + java in PATH
SNPEFF_JAVA="$(find "$ENV_GENOMICS/lib/jvm/bin" -name java 2>/dev/null | head -1 || true)"
[ -z "$SNPEFF_JAVA" ] && SNPEFF_JAVA="$(command -v java 2>/dev/null || true)"
[ -z "$SNPEFF_JAVA" ] && die "java not found — required for SnpEff"
export PATH="$ENV_GENOMICS/bin:$(dirname "$SNPEFF_JAVA"):$PATH"
log "Java: $SNPEFF_JAVA"

# SnpEff works directly with BCF/VCF; no chr-prefix stripping needed.
# -canon: use only canonical transcripts (one effect per variant)
# -noStats: skip the HTML/CSV summary report (faster)
# -noLog: suppress usage reporting
# SnpEff cannot read BCF; pipe bcftools view to convert on the fly
# Use explicit paths to avoid $MM env-prefix issues with PATH inheritance
BCFTOOLS_BIN="$ENV_GENOMICS/bin/bcftools"
SNPEFF_BIN="$ENV_GENOMICS/bin/snpEff"
BGZIP_BIN="$ENV_GENOMICS/bin/bgzip"
"$BCFTOOLS_BIN" view "$IMPUTED_BCF" \
  | "$SNPEFF_BIN" ann \
      -canon \
      -noStats \
      -noLog \
      -v \
      "$SNPEFF_DB" \
      2>"$ANN_DIR/snpeff.log" \
  | "$BGZIP_BIN" -c > "$ANN_DIR/${DOG_LOWER}_annotated.vcf.gz"
"$BCFTOOLS_BIN" index -t "$ANN_DIR/${DOG_LOWER}_annotated.vcf.gz"
log "SnpEff done: $(wc -l < $ANN_DIR/snpeff.log) log lines"

# ── Rebuild cnv_genes.json from annotated VCF (genome-wide, this sample) ──
# Stage 6 uses a static reference cnv_genes.json that only covers chromosomes
# with CNVs in the reference dog. Re-annotate now using this sample's SnpEff
# output, which covers all chromosomes, then patch cnv_homdel.json.
log "  Rebuilding cnv_genes.json from SnpEff annotation…"
python3 - << PYEOF
import gzip, json, re
from collections import defaultdict

ANN_VCF  = "$ANN_DIR/${DOG_LOWER}_annotated.vcf.gz"
CNV_JSON = "$PUB/cnv_homdel.json"
PUB      = "$PUB"

# Effect types that indicate a variant falls in an exon
EXONIC_EFFECTS = {
    'exon_variant', 'missense_variant', 'synonymous_variant',
    'stop_gained', 'stop_lost', 'start_lost', 'frameshift_variant',
    'splice_acceptor_variant', 'splice_donor_variant',
    'protein_protein_contact', 'structural_interaction_variant',
    'inframe_insertion', 'inframe_deletion',
    'stop_retained_variant', 'start_retained_variant',
    'coding_sequence_variant', '5_prime_UTR_variant', '3_prime_UTR_variant',
}

# Parse ANN fields: build gene coordinate map and track per-gene exonic positions
gene_by_id = {}   # gene_id -> gene dict with exonic_positions set
with gzip.open(ANN_VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        cols = line.split('\t')
        if len(cols) < 8: continue
        chrom, pos = cols[0], int(cols[1])
        ann_match = re.search(r'ANN=([^;]+)', cols[7])
        if not ann_match: continue
        for ann in ann_match.group(1).split(','):
            parts = ann.split('|')
            if len(parts) < 8: continue
            effect    = parts[1]
            gene_name = parts[3]
            gene_id   = parts[4]
            biotype   = parts[7]
            if not gene_id or not gene_name: continue
            # Skip fusion/readthrough annotations (two gene names joined by '-')
            if '-' in gene_name:
                p = gene_name.split('-')
                if len(p) == 2 and all(len(x) > 4 for x in p):
                    continue
            is_exonic = any(e in effect for e in EXONIC_EFFECTS)
            if gene_id not in gene_by_id:
                gene_by_id[gene_id] = {
                    'gene_id': gene_id, 'name': gene_name, 'chrom': chrom,
                    'start': pos, 'end': pos, 'biotype': biotype,
                    'strand': '+', 'exons': [], 'cds': [],
                    '_exonic_pos': set()
                }
            g = gene_by_id[gene_id]
            g['start'] = min(g['start'], pos)
            g['end']   = max(g['end'],   pos)
            if is_exonic:
                g['_exonic_pos'].add(pos)

# Build chromosome-keyed map
gene_map = defaultdict(list)
for g in gene_by_id.values():
    chrom_num = g['chrom'].replace('chr', '')
    gene_map[chrom_num].append({
        'gene_id': g['gene_id'], 'name': g['name'],
        'start': g['start'], 'end': g['end'],
        'strand': g['strand'], 'biotype': g['biotype'],
        'exons': g['exons'], 'cds': g['cds'],
        '_exonic_pos': sorted(g['_exonic_pos'])
    })

print(f"cnv_genes.json: {len(gene_by_id)} genes across {len(gene_map)} chromosomes")
# Write without internal _exonic_pos (that's only for CNV re-annotation)
gene_map_out = {c: [{k: v for k, v in g.items() if k != '_exonic_pos'} for g in gs]
                for c, gs in gene_map.items()}
with open(f'{PUB}/cnv_genes.json', 'w') as f:
    json.dump(gene_map_out, f)

# Re-annotate CNV regions using gene map + exonic position evidence
with open(CNV_JSON) as f:
    cnv = json.load(f)

def find_genes(chrom, start, end):
    chrom_num = chrom.replace('chr', '')
    disrupted, details = [], []
    for g in gene_map.get(chrom_num, []):
        if g['end'] < start or g['start'] > end: continue
        ov = 'full' if g['start'] >= start and g['end'] <= end else 'partial'
        # Exon overlap: check if any exonic variant position falls in the CNV window
        exonic_in_region = [p for p in g.get('_exonic_pos', []) if start <= p <= end]
        exon_ov = 'exonic' if exonic_in_region else 'intronic'
        disrupted.append(g['name'])
        details.append({'gene': g['name'], 'biotype': g['biotype'],
                        'chrom': chrom, 'start': g['start'], 'end': g['end'],
                        'overlap': ov, 'exon_overlap': exon_ov})
    return disrupted, details

all_disrupted = {}
for r in cnv.get('regions', []):
    genes, dets = find_genes(r['chrom'], r['start'], r['end'])
    r['disrupted_genes'] = genes
    r['disrupted_gene_details'] = dets
    for d in dets:
        all_disrupted[d['gene']] = d

cnv['disrupted_genes'] = list(all_disrupted.values())
cnv['summary']['unique_genes'] = len(all_disrupted)
with open(CNV_JSON, 'w') as f:
    json.dump(cnv, f, indent=2)
print(f"cnv_homdel.json re-annotated: {len(all_disrupted)} disrupted genes")
for r in cnv.get('regions', []):
    print(f"  {r['chrom']}:{r['start']}-{r['end']} → "
          + ", ".join(f"{d['gene']} ({d['exon_overlap']})" for d in r['disrupted_gene_details']))
PYEOF

# Parse SnpEff ANN field → functional_variants.json
# ANN format (pipe-delimited per transcript, comma-separated per variant):
#   ALT | effect | impact | gene_name | gene_id | feature_type | feature_id |
#   biotype | rank | hgvsc | hgvsp | cdna_pos | cds_pos | aa_pos | distance | messages
python3 - << PYEOF
import gzip, re, subprocess, json

ANN_VCF = "$ANN_DIR/${DOG_LOWER}_annotated.vcf.gz"
BCF     = "$IMPUTED_BCF"
PUB     = "$PUB"

# Build AF lookup from imputed BCF (chr-prefixed keys)
print("Building AF lookup...")
af_lookup = {}
result = subprocess.run(
    ['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RAF\n', BCF],
    capture_output=True, text=True)
for line in result.stdout.strip().split('\n'):
    if not line: continue
    parts = line.split('\t')
    if len(parts) < 5: continue
    try:
        af_lookup[(parts[0], int(parts[1]), parts[3])] = float(parts[4].split(',')[0])
    except Exception:
        pass
print(f"  AF lookup: {len(af_lookup)} sites")

HIGH = []
HIGH_COUNTS = {}
MOD_BY_GENE = {}
seen = set()   # deduplicate: (chrom, pos, ref, alt, gene, impact)

with gzip.open(ANN_VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        cols = line.strip().split('\t')
        if len(cols) < 8: continue
        chrom, pos_s, _, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
        pos = int(pos_s)

        # Genotype
        gt_str = ''
        if len(cols) >= 10:
            fmt = cols[8].split(':')
            smp = cols[9].split(':')
            if 'GT' in fmt:
                gt_str = smp[fmt.index('GT')]
        alleles = re.split(r'[|/]', gt_str)
        if set(alleles) <= {'0', '.'}: continue
        zyg = ('hom_alt' if set(alleles) == {'1'} else
               'het'     if '0' in alleles and '1' in alleles else 'other')

        af = af_lookup.get((chrom, pos, alt))

        # Parse ANN field
        info = dict(x.split('=', 1) for x in cols[7].split(';') if '=' in x)
        ann_str = info.get('ANN', '')
        if not ann_str: continue

        best = {}   # gene → best (highest impact) annotation for this variant
        for ann in ann_str.split(','):
            fields = ann.split('|')
            if len(fields) < 4: continue
            ann_alt    = fields[0]
            effect     = fields[1]
            impact     = fields[2]
            gene       = fields[3]
            if not gene or impact not in ('HIGH', 'MODERATE'): continue
            if ann_alt != alt: continue   # skip if annotation is for a different ALT
            # keep highest-impact annotation per gene
            rank = {'HIGH': 0, 'MODERATE': 1}
            if gene not in best or rank[impact] < rank[best[gene]['impact']]:
                best[gene] = {'impact': impact, 'effect': effect}

        for gene, ann_data in best.items():
            impact = ann_data['impact']
            effect = ann_data['effect']
            key = (chrom, pos, ref, alt, gene, impact)
            if key in seen: continue
            seen.add(key)

            row = {'impact': impact, 'gene': gene,
                   'chr': chrom.replace('chr', ''), 'pos': str(pos),
                   'ref': ref, 'alt': alt, 'effect': effect, 'zygosity': zyg,
                   'af_dog10k': round(af, 6) if af is not None else None}

            if impact == 'HIGH':
                HIGH.append(row)
                base_effect = effect.split('&')[0]
                HIGH_COUNTS[base_effect] = HIGH_COUNTS.get(base_effect, 0) + 1
            else:
                MOD_BY_GENE.setdefault(gene, []).append(row)

def rare(thresh):
    return sum(1 for r in HIGH if r['zygosity'] == 'hom_alt'
               and r['af_dog10k'] is not None and r['af_dog10k'] < thresh)

h_hom     = [r for r in HIGH if r['zygosity'] == 'hom_alt']
h_het     = [r for r in HIGH if r['zygosity'] == 'het']
mod_total = sum(len(v) for v in MOD_BY_GENE.values())
mod_hom   = sum(sum(1 for r in v if r['zygosity'] == 'hom_alt') for v in MOD_BY_GENE.values())
mod_het   = sum(sum(1 for r in v if r['zygosity'] == 'het')     for v in MOD_BY_GENE.values())

high_sorted = sorted(HIGH, key=lambda r: (
    0 if r['zygosity'] == 'hom_alt' else 1,
    r['af_dog10k'] if r['af_dog10k'] is not None else 1.0))

mod_gene_list = []
for g, rows in sorted(MOD_BY_GENE.items(), key=lambda kv: -len(kv[1])):
    hom_rows = [r for r in rows if r['zygosity'] == 'hom_alt']
    het_rows = [r for r in rows if r['zygosity'] == 'het']
    effects = list(dict.fromkeys(r['effect'].split('&')[0] for r in rows))
    hom_afs = [r['af_dog10k'] for r in hom_rows if r['af_dog10k'] is not None]
    min_af = min(hom_afs) if hom_afs else None
    mod_gene_list.append({
        'gene': g,
        'n_moderate': len(rows),
        'hom_alt': len(hom_rows),
        'het': len(het_rows),
        'effects': effects,
        'min_af': round(min_af, 6) if min_af is not None else None,
    })

fv = {
    'summary': {
        'high_total':          len(HIGH),
        'high_hom_alt':        len(h_hom),
        'high_het':            len(h_het),
        'high_hom_rare_1pct':  rare(0.01),
        'high_hom_rare_5pct':  rare(0.05),
        'high_hom_rare_10pct': rare(0.10),
        'moderate_total':      mod_total,
        'moderate_hom_alt':    mod_hom,
        'moderate_het':        mod_het,
    },
    'high_effect_counts': HIGH_COUNTS,
    'high_variants':      high_sorted,
    'moderate_by_gene':   mod_gene_list,
    'source':   f'GLIMPSE2 Dog10K imputation + SnpEff ($SNPEFF_DB), canonical transcripts',
    'af_note':  'AF = allele frequency in Dog10K reference panel',
}
with open(f'{PUB}/functional_variants.json', 'w') as f:
    json.dump(fv, f)
print(f"functional_variants.json: {len(HIGH)} HIGH, {mod_total} MODERATE variants")
PYEOF

# ── Stage 11: PRS from imputed dosages ──────────────────────
log "=== Stage 11: PRS (imputed Parker dosages) ==="
python3 - << PYEOF
import subprocess, numpy as np, json, csv, io, tempfile, os

BCF      = "$IMPUTED_BCF"
BIM      = "$PARKER_BIM"
FAM      = "$PARKER_FAM"
BED      = "$PARKER_BED"
PUB      = "$PUB"
REF_PRS  = "$D/dogs-app/public/cosmo/prs_result.json"

# ── Parker breed code → AKC breed name (plural form used by kkakey/dog_traits_AKC) ──
PARKER_TO_AKC = {
    # Breed names match kkakey/dog_traits_AKC exactly (AKC uses "Retrievers (X)" etc.)
    'ACKR':'American Cocker Spaniels','AFGH':'Afghan Hounds','AIRT':'Airedale Terriers',
    'AKIT':'Akitas','AMAL':'Alaskan Malamutes','AMST':'American Staffordshire Terriers',
    'AHRT':'American Hairless Terriers','ANAT':'Anatolian Shepherd Dogs',
    'AESK':'American Eskimo Dogs','AUCD':'Australian Cattle Dogs',
    'AUST':'Australian Terriers','AUSS':'Australian Shepherds',
    'AZWK_Mali':'Azawakhs','BASS':'Basset Hounds','BEAG':'Beagles',
    'BEDT':'Bedlington Terriers','BELS':'Belgian Sheepdogs','BELM':'Belgian Malinois',
    'TURV':'Belgian Tervuren','BELT':'Belgian Tervuren','BERD':'Bergamasco Sheepdogs',
    'BPIC':'Berger Picards',
    'BMAL':'Bernese Mountain Dogs','BMD':'Bernese Mountain Dogs',
    'BICH':'Bichons Frises','BLDH':'Bloodhounds',
    'BOER':'Boerboels','BORD':'Border Collies',
    'BORT':'Border Terriers','BORZ':'Borzois',
    'BOST':'Boston Terriers','BOUV':'Bouviers des Flandres',
    'BOX':'Boxers','BRIA':'Briards',
    'BRIT':'Brittanys','BRTR':'Brittanys','BRUS':'Brussels Griffons',
    'BSJI':'Basenjis',
    'BULD':'Bulldogs','BULM':'Bullmastiffs','BULT':'Bull Terriers',
    'MBLT':'Miniature Bull Terriers','CAIR':'Cairn Terriers',
    'CANE':'Cane Corso','CANE_Italy':'Cane Corso',
    'CARD':'Cardigan Welsh Corgis',
    'CKCS':'Cavalier King Charles Spaniels','CCRT':'Retrievers (Curly-Coated)',
    'CHIH':'Chihuahuas',
    'CHIN':'Japanese Chin','CRES':'Chinese Crested','CHOW':'Chow Chows',
    'COOK':'Spaniels (Cocker)',
    'ESSP':'Spaniels (English Springer)',
    'COLL':'Collies','SSHP':'Shetland Sheepdogs',
    'COTO':'Coton de Tulear','DACH':'Dachshunds','DALM':'Dalmatians',
    'DANE':'Great Danes','DDBX':'Dogues de Bordeaux',
    'DEER':'Scottish Deerhounds','DOBP':'Doberman Pinschers',
    'ECKR':'English Cocker Spaniels','ESET':'Setters (English)',
    'FBUL':'French Bulldogs','FCR':'Retrievers (Flat-Coated)',
    'FIEL':'Spaniels (Field)','FINS':'Finnish Spitz',
    'WFOX':'Fox Terriers (Wire)',
    'GSD':'German Shepherd Dogs',
    'GSHP':'Pointers (German Shorthaired)','GSNZ':'Giant Schnauzers',
    'GOLD':'Retrievers (Golden)','GORD':'Setters (Gordon)',
    'GREY':'Greyhounds','GREE':'Greyhounds','GLEN':'Glen of Imaal Terriers',
    'GPYR':'Great Pyrenees','GSMD':'Greater Swiss Mountain Dogs',
    'GWHP':'Pointers (German Wirehaired)','HAVA':'Havanese',
    'HUSK':'Siberian Huskies','IBIZ':'Ibizan Hounds',
    'ICES':'Icelandic Sheepdogs','ISET':'Setters (Irish)',
    'IRIT':'Irish Terriers',
    'IWSP':'Spaniels (Irish Water)','IWOF':'Irish Wolfhounds','ITGY':'Italian Greyhounds',
    'JACK':'Russell Terriers',
    'KEES':'Keeshonden',
    'KERY':'Kerry Blue Terriers','KOMO':'Komondorok',
    'KUVZ':'Kuvaszok',
    'LAB':'Retrievers (Labrador)',
    'LHSA':'Lhasa Apsos',
    'LEON':'Leonbergers','MALT':'Maltese','MAST':'Mastiffs',
    'MPOO':'Poodles','MSNZ':'Miniature Schnauzers',
    'MPIN':'Miniature Pinschers',
    'NEAP':'Neapolitan Mastiffs','NEWF':'Newfoundlands',
    'NORF':'Norfolk Terriers',
    'NSDT':'Retrievers (Nova Scotia Duck Tolling)',
    'OES':'Old English Sheepdogs',
    'OTTR':'Otterhounds','PAPI':'Papillons',
    'PARS':'Parson Russell Terriers',
    'PBGV':'Petits Bassets Griffons Vendeens',
    'PEKE':'Pekingese','PEMB':'Pembroke Welsh Corgis',
    'PHAR':'Pharaoh Hounds',
    'SPOO':'Poodles','TPOO':'Poodles',
    'PTWD':'Portuguese Water Dogs',
    'POM':'Pomeranians','PUG':'Pugs','PULI':'Pulik','PUMI':'Pumik',
    'RATT':'Rat Terriers','RHOD':'Rhodesian Ridgebacks',
    'ROTT':'Rottweilers',
    'SALU':'Salukis','SALU_ArabPen':'Salukis','SALU_CentAsia':'Salukis','SALU_Tribal':'Salukis',
    'SAMO':'Samoyeds',
    'SCOT':'Scottish Terriers','SCWT':'Soft Coated Wheaten Terriers',
    'SHAR':'Chinese Shar-Pei',
    'SHIB':'Shiba Inu','SHIH':'Shih Tzu','SILK':'Silky Terriers',
    'SKIP':'Skye Terriers','SLOU_NAfrica':'Sloughis','SPIN':'Spinoni Italiani',
    'SSNZ':'Standard Schnauzers','STAF':'Staffordshire Bull Terriers',
    'TIBS':'Tibetan Spaniels','TIBT':'Tibetan Terriers',
    'TIBM':'Tibetan Mastiffs','TIBM_China':'Tibetan Mastiffs',
    'TYFX':'Fox Terriers (Smooth)',
    'YORK':'Yorkshire Terriers','VIZS':'Vizslas',
    'SVAL':'Swedish Vallhunds',
    'WEIM':'Weimaraners','WHWT':'West Highland White Terriers',
    'WELS':'Spaniels (Welsh Springer)',
    'WHIP':'Whippets','WHPG':'Wirehaired Pointing Griffons',
    'XOLO':'Xoloitzcuintli','MXOL':'Xoloitzcuintli',
}

TRAIT_COLS = [
    'Affectionate With Family','Good With Young Children','Good With Other Dogs',
    'Shedding Level','Coat Grooming Frequency','Drooling Level',
    'Openness To Strangers','Playfulness Level','Watchdog/Protective Nature',
    'Adaptability Level','Trainability Level','Energy Level',
    'Barking Level','Mental Stimulation Needs',
]

# ── Fetch AKC trait data ──────────────────────────────────────────────────
print("Fetching AKC trait data...")
r = subprocess.run(['curl', '-sL',
    'https://raw.githubusercontent.com/kkakey/dog_traits_AKC/main/data/breed_traits.csv'],
    capture_output=True, text=True)
akc_rows = list(csv.DictReader(io.StringIO(r.stdout)))
akc_by_breed = {row['Breed'].strip().replace('\xa0',' '): row for row in akc_rows}
TRAIT_COLS = [c for c in TRAIT_COLS if c in next(iter(akc_by_breed.values()))]
print(f"  {len(akc_by_breed)} AKC breeds, {len(TRAIT_COLS)} trait columns")

# ── Parker panel ──────────────────────────────────────────────────────────
breeds_fam = []
with open(FAM) as f:
    for line in f:
        breeds_fam.append(line.split()[0])
breeds_fam = np.array(breeds_fam)
n_ref = len(breeds_fam)

parker_snps = []
with open(BIM) as f:
    for line in f:
        p = line.strip().split()
        parker_snps.append({'chrom': p[0], 'pos': int(p[3]), 'a1': p[4], 'a2': p[5]})
n_snps = len(parker_snps)
print(f"Parker: {n_snps} SNPs, {n_ref} samples ({len(np.unique(breeds_fam))} breeds)")

def load_bed(bed_path, n_samples, n_snps):
    with open(bed_path, 'rb') as f:
        magic = f.read(3)
        assert magic == b'\x6c\x1b\x01', "Invalid BED file"
        bps = (n_samples + 3) // 4
        G = np.zeros((n_snps, n_samples), dtype=np.float32)
        for i in range(n_snps):
            rb = np.frombuffer(f.read(bps), dtype=np.uint8)
            for j in range(n_samples):
                g = (rb[j // 4] >> ((j % 4) * 2)) & 0x03
                G[i, j] = [0, np.nan, 1, 2][g]
    return G

print("Loading Parker BED...")
G_ref = load_bed(BED, n_ref, n_snps)

# ── Imputed dosages at Parker positions ───────────────────────────────────
print("Extracting imputed dosages at Parker SNP positions...")
_bed = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
for s in parker_snps:
    _bed.write(f"chr{s['chrom']}\t{s['pos']-1}\t{s['pos']}\n")
_bed.close()
result = subprocess.run(
    ['bcftools', 'query', '-R', _bed.name, '-f', '%CHROM\t%POS\t%REF\t%ALT\t[%GP]\n', BCF],
    capture_output=True, text=True)
os.unlink(_bed.name)

# Build position index for fast lookup
pos_index = {(s['chrom'], s['pos']): i for i, s in enumerate(parker_snps)}

# Orient dosages to Parker A2 convention (copies of A2 to match BED encoding)
sample_dosage = np.full(n_snps, np.nan)
genotyped = 0
for line in result.stdout.strip().split('\n'):
    if not line: continue
    parts = line.split('\t')
    if len(parts) < 5: continue
    chrom = parts[0].replace('chr', ''); pos = int(parts[1])
    ref, alt = parts[2], parts[3]
    idx = pos_index.get((chrom, pos))
    if idx is None: continue
    s = parker_snps[idx]
    if ref not in (s['a1'], s['a2']) or alt not in (s['a1'], s['a2']): continue
    try:
        gp = [float(x) for x in parts[4].split(',')]
        if len(gp) < 3: continue
        # Need E[copies of A2] to match Parker BED encoding
        if alt == s['a2']:
            sample_dosage[idx] = gp[1] + 2.0 * gp[2]
        else:  # alt == a1, ref == a2
            sample_dosage[idx] = 2.0 * gp[0] + gp[1]
        genotyped += 1
    except: pass

print(f"  {genotyped}/{n_snps} Parker SNPs covered (allele-oriented)")

# ── LD pruning from valid SNPs only + NaN fill ───────────────────────────
valid_all  = np.where(~np.isnan(sample_dosage))[0]
prune_idx  = valid_all[::5]           # every 5th valid SNP
valid_snps = np.ones(len(prune_idx), dtype=bool)  # all pruned are valid by construction

G_sub_raw = G_ref[prune_idx, :]
row_means  = np.nanmean(G_sub_raw, axis=1, keepdims=True)
G_sub      = np.where(np.isnan(G_sub_raw), row_means, G_sub_raw)
s_sub      = sample_dosage[prune_idx]
print(f"  LD-pruned valid SNPs: {len(prune_idx)}")

# ── Ridge-regularized GWAS PRS ────────────────────────────────────────────
def compute_prs_ridge(breed_scores_by_code, G_sub, s_sub, breeds, lambda_frac=0.1):
    y = np.array([breed_scores_by_code.get(b, np.nan) for b in breeds])
    valid = ~np.isnan(y)
    if valid.sum() < 30: return np.nan, np.nan
    y_v = y[valid]; G_v = G_sub[:, valid]
    G_c = G_v - G_v.mean(axis=1, keepdims=True)
    y_c = y_v - y_v.mean()
    var_j = np.sum(G_c**2, axis=1)
    beta  = np.dot(G_c, y_c) / (var_j + lambda_frac * var_j.mean())
    prs_raw = np.dot(beta, s_sub)
    prs_ref = np.dot(G_sub.T, beta)
    prs_ref_v = prs_ref[valid]
    prs_z = (prs_raw - prs_ref_v.mean()) / (prs_ref_v.std() + 1e-8)
    percentile = float(np.mean(prs_ref_v <= prs_raw) * 100)
    return prs_z, percentile

# ── Compute per trait ─────────────────────────────────────────────────────
print("Computing PRS per trait:")
traits_out = {}
for trait in TRAIT_COLS:
    breed_scores = {}
    for code in np.unique(breeds_fam):
        akc_name = PARKER_TO_AKC.get(code)
        if not akc_name or akc_name not in akc_by_breed: continue
        try:
            val = akc_by_breed[akc_name].get(trait, '').strip()
            if val: breed_scores[code] = float(val)
        except: pass
    if len(breed_scores) < 20: continue
    prs_z, pct = compute_prs_ridge(breed_scores, G_sub, s_sub, breeds_fam)
    if np.isnan(prs_z): continue
    all_scores = list(breed_scores.values())
    predicted  = float(np.clip(np.mean(all_scores) + prs_z * np.std(all_scores), 1, 5))
    traits_out[trait] = {
        'prs_z': round(float(prs_z), 3),
        'percentile': round(float(pct), 1),
        'predicted_score': round(predicted, 2),
        'n_ref_samples': int(len(prune_idx)),
        'description': f'GWAS prediction for {trait.lower()}',
    }
    print(f"  {trait}: z={prs_z:.3f}, pct={pct:.1f}, pred={predicted:.2f}")

# Carry heritability annotations from cosmo reference result
with open(REF_PRS) as f:
    ref_prs = json.load(f)
for t, v in traits_out.items():
    ref = ref_prs['traits'].get(t) or {}
    if ref.get('heritability'): v['heritability'] = ref['heritability']

# ── Physical traits PRS ───────────────────────────────────────────────────
print("Computing physical trait PRS...")

BREED_HEIGHT_CM = {
    'Afghan Hounds':68,'Airedale Terriers':58,'Akitas':65,'Alaskan Malamutes':61,
    'Australian Cattle Dogs':47,'Australian Shepherds':52,'Basenjis':42,
    'Basset Hounds':33,'Beagles':34,'Bearded Collies':53,'Belgian Malinois':60,
    'Belgian Sheepdogs':61,'Belgian Tervuren':61,'Bernese Mountain Dogs':64,
    'Bichons Frises':28,'Bloodhounds':62,'Border Collies':51,'Border Terriers':33,
    'Borzois':71,'Boston Terriers':38,'Bouviers des Flandres':64,'Boxers':58,
    'Briards':63,'Brittanys':50,'Bulldogs':38,'Bullmastiffs':66,'Cairn Terriers':28,
    'Cane Corso':66,'Cardigan Welsh Corgis':30,'Cavalier King Charles Spaniels':31,
    'Retrievers (Chesapeake Bay)':60,'Chihuahuas':18,'Chinese Crested':30,
    'Chinese Shar-Pei':48,'Chow Chows':48,'Spaniels (Cocker)':37,'Collies':61,
    'Retrievers (Curly-Coated)':65,'Dalmatians':56,'Doberman Pinschers':66,
    'Dogues de Bordeaux':63,'Spaniels (English Cocker)':40,'Setters (English)':63,
    'Spaniels (English Springer)':50,'Retrievers (Flat-Coated)':60,'French Bulldogs':30,
    'German Shepherd Dogs':62,'Pointers (German Shorthaired)':60,
    'Pointers (German Wirehaired)':63,'Retrievers (Golden)':58,'Setters (Gordon)':65,
    'Great Danes':79,'Great Pyrenees':71,'Greater Swiss Mountain Dogs':67,
    'Greyhounds':70,'Havanese':23,'Ibizan Hounds':60,'Setters (Irish)':67,
    'Irish Terriers':46,'Spaniels (Irish Water)':57,'Irish Wolfhounds':81,
    'Italian Greyhounds':35,'Keeshonden':46,'Kerry Blue Terriers':47,
    'Komondorok':70,'Kuvaszok':70,'Retrievers (Labrador)':57,'Leonbergers':72,
    'Lhasa Apsos':25,'Maltese':23,'Mastiffs':76,'Miniature Bull Terriers':33,
    'Miniature Pinschers':28,'Miniature Schnauzers':32,'Newfoundlands':69,
    'Norwegian Elkhounds':50,'Retrievers (Nova Scotia Duck Tolling)':50,
    'Old English Sheepdogs':56,'Papillons':23,'Parson Russell Terriers':33,
    'Pembroke Welsh Corgis':27,'Pomeranians':18,'Poodles':38,
    'Portuguese Water Dogs':52,'Pugs':30,'Rhodesian Ridgebacks':64,
    'Rottweilers':63,'Russell Terriers':28,'Salukis':66,'Samoyeds':55,
    'Scottish Deerhounds':76,'Shetland Sheepdogs':37,'Shiba Inu':38,'Shih Tzu':25,
    'Siberian Huskies':56,'Silky Terriers':23,'Soft Coated Wheaten Terriers':47,
    'Standard Schnauzers':47,'Tibetan Mastiffs':71,'Tibetan Spaniels':25,
    'Tibetan Terriers':38,'Vizslas':58,'Weimaraners':65,'Spaniels (Welsh Springer)':46,
    'West Highland White Terriers':27,'Whippets':51,'Fox Terriers (Wire)':38,
    'Yorkshire Terriers':18,'Entlebucher Mountain Dogs':50,
}
BREED_WEIGHT_KG = {
    'Afghan Hounds':25,'Airedale Terriers':24,'Akitas':40,'Alaskan Malamutes':38,
    'Australian Cattle Dogs':20,'Australian Shepherds':25,'Basenjis':10,
    'Basset Hounds':22,'Beagles':9,'Bearded Collies':22,'Belgian Malinois':28,
    'Belgian Sheepdogs':28,'Belgian Tervuren':28,'Bernese Mountain Dogs':40,
    'Bichons Frises':5.5,'Bloodhounds':45,'Border Collies':17,'Border Terriers':6,
    'Borzois':35,'Boston Terriers':8,'Bouviers des Flandres':38,'Boxers':30,
    'Briards':35,'Brittanys':16,'Bulldogs':23,'Bullmastiffs':55,'Cairn Terriers':6.5,
    'Cane Corso':50,'Cardigan Welsh Corgis':14,'Cavalier King Charles Spaniels':7,
    'Retrievers (Chesapeake Bay)':30,'Chihuahuas':2.5,'Chinese Crested':5,
    'Chinese Shar-Pei':22,'Chow Chows':28,'Spaniels (Cocker)':12,'Collies':28,
    'Retrievers (Curly-Coated)':32,'Dalmatians':25,'Doberman Pinschers':35,
    'Dogues de Bordeaux':55,'Spaniels (English Cocker)':13,'Setters (English)':28,
    'Spaniels (English Springer)':22,'Retrievers (Flat-Coated)':32,'French Bulldogs':11,
    'German Shepherd Dogs':32,'Pointers (German Shorthaired)':27,
    'Pointers (German Wirehaired)':32,'Retrievers (Golden)':30,'Setters (Gordon)':28,
    'Great Danes':65,'Great Pyrenees':45,'Greater Swiss Mountain Dogs':55,
    'Greyhounds':30,'Havanese':5,'Ibizan Hounds':22,'Setters (Irish)':30,
    'Irish Terriers':12,'Spaniels (Irish Water)':26,'Irish Wolfhounds':55,
    'Italian Greyhounds':4,'Keeshonden':18,'Kerry Blue Terriers':16,
    'Komondorok':50,'Kuvaszok':50,'Retrievers (Labrador)':30,'Leonbergers':55,
    'Lhasa Apsos':6,'Maltese':3,'Mastiffs':90,'Miniature Bull Terriers':11,
    'Miniature Pinschers':4,'Miniature Schnauzers':7,'Newfoundlands':60,
    'Norwegian Elkhounds':22,'Retrievers (Nova Scotia Duck Tolling)':21,
    'Old English Sheepdogs':34,'Papillons':4,'Parson Russell Terriers':6.5,
    'Pembroke Welsh Corgis':12,'Pomeranians':2.5,'Poodles':8,
    'Portuguese Water Dogs':21,'Pugs':8,'Rhodesian Ridgebacks':36,'Rottweilers':48,
    'Russell Terriers':5,'Salukis':23,'Samoyeds':24,'Scottish Deerhounds':45,
    'Shetland Sheepdogs':7,'Shiba Inu':9,'Shih Tzu':6,'Siberian Huskies':22,
    'Silky Terriers':4,'Soft Coated Wheaten Terriers':16,'Standard Schnauzers':17,
    'Tibetan Mastiffs':60,'Tibetan Spaniels':5,'Tibetan Terriers':11,'Vizslas':25,
    'Weimaraners':32,'Spaniels (Welsh Springer)':18,'West Highland White Terriers':8,
    'Whippets':12,'Wirehaired Pointing Griffons':27,'Fox Terriers (Wire)':7.5,
    'Yorkshire Terriers':3,'Entlebucher Mountain Dogs':25,
}
COAT_LEN_ORD = {'Short':1,'Medium':2,'Long':3}

def cont_scores(lookup):
    return {c: lookup[PARKER_TO_AKC[c]] for c in np.unique(breeds_fam)
            if PARKER_TO_AKC.get(c) in lookup}

def ord_scores_phys(col, omap):
    out = {}
    for c in np.unique(breeds_fam):
        akc = PARKER_TO_AKC.get(c)
        if not akc or akc not in akc_by_breed: continue
        v = akc_by_breed[akc].get(col,'').strip().replace('\xa0',' ')
        if v in omap: out[c] = float(omap[v])
    return out

def bin_scores_phys(col, target):
    out = {}
    for c in np.unique(breeds_fam):
        akc = PARKER_TO_AKC.get(c)
        if not akc or akc not in akc_by_breed: continue
        v = akc_by_breed[akc].get(col,'').strip().replace('\xa0',' ')
        if v: out[c] = 1.0 if v == target else 0.0
    return out

phys_traits = {}

# Height
h_sc = cont_scores(BREED_HEIGHT_CM)
z_h, pct_h = compute_prs_ridge(h_sc, G_sub, s_sub, breeds_fam)
if not np.isnan(z_h):
    vals = list(h_sc.values()); mu, sd = np.mean(vals), np.std(vals)
    pred_h = float(np.clip(mu + z_h * sd, 20, 110))
    phys_traits['height_cm'] = {
        'pred_cm': round(pred_h, 1),
        'prs_z': round(float(z_h), 3),
        'percentile': round(float(pct_h), 1),
        'n_ref_samples': int(len(prune_idx)),
        'description': 'Predicted adult height at withers.',
        'heritability': {'h2': 0.62, 'ci': '0.55–0.69', 'source': 'Hayward 2016 (Nat Gen)'},
    }
    print(f"  Height: {pred_h:.1f}cm (z={z_h:.3f}, pct={pct_h:.1f})")

# Weight
w_sc = cont_scores(BREED_WEIGHT_KG)
z_w, pct_w = compute_prs_ridge(w_sc, G_sub, s_sub, breeds_fam)
if not np.isnan(z_w):
    vals = list(w_sc.values()); mu, sd = np.mean(vals), np.std(vals)
    pred_w = float(np.clip(mu + z_w * sd, 1, 120))
    phys_traits['weight_kg'] = {
        'pred_kg': round(pred_w, 1),
        'pred_lbs': round(pred_w * 2.205, 1),
        'prs_z': round(float(z_w), 3),
        'percentile': round(float(pct_w), 1),
        'n_ref_samples': int(len(prune_idx)),
        'description': 'Predicted adult weight.',
        'heritability': {'h2': 0.60, 'ci': '0.52–0.68', 'source': 'Hayward 2016 (Nat Gen)'},
    }
    print(f"  Weight: {pred_w:.1f}kg / {pred_w*2.205:.1f}lbs (z={z_w:.3f}, pct={pct_w:.1f})")

# Coat type — per-category binary, pick highest z-score
coat_types = ['Double','Smooth','Wavy','Curly','Silky','Wiry','Rough']
ct_best, ct_best_z, ct_best_pct = None, -np.inf, 50.0
for ct in coat_types:
    sc = bin_scores_phys('Coat Type', ct)
    if len(sc) < 20: continue
    z, pct = compute_prs_ridge(sc, G_sub, s_sub, breeds_fam)
    if not np.isnan(z) and z > ct_best_z:
        ct_best, ct_best_z, ct_best_pct = ct, z, pct
if ct_best:
    phys_traits['coat_type'] = {
        'predicted': ct_best,
        'prs_z': round(float(ct_best_z), 3),
        'percentile': round(float(ct_best_pct), 1),
        'n_ref_samples': int(len(prune_idx)),
        'description': 'Predicted coat texture/type.',
        'heritability': {'h2': 0.32, 'ci': '0.18–0.46', 'source': 'Parker 2017 (Science)'},
    }
    print(f"  Coat type: {ct_best} (z={ct_best_z:.3f})")

# Coat length — ordinal Short=1, Medium=2, Long=3
cl_sc = ord_scores_phys('Coat Length', COAT_LEN_ORD)
z_cl, pct_cl = compute_prs_ridge(cl_sc, G_sub, s_sub, breeds_fam)
if not np.isnan(z_cl):
    vals = list(cl_sc.values()); mu, sd = np.mean(vals), np.std(vals)
    pred_ord = float(np.clip(mu + z_cl * sd, 1, 3))
    pred_cl = min(COAT_LEN_ORD, key=lambda k: abs(COAT_LEN_ORD[k] - pred_ord))
    phys_traits['coat_length'] = {
        'predicted': pred_cl,
        'pred_ordinal': round(pred_ord, 2),
        'prs_z': round(float(z_cl), 3),
        'percentile': round(float(pct_cl), 1),
        'n_ref_samples': int(len(prune_idx)),
        'description': 'Predicted coat length.',
        'heritability': {'h2': 0.63, 'ci': '0.52–0.74', 'source': 'Parker 2017 (Science)'},
    }
    print(f"  Coat length: {pred_cl} (ord={pred_ord:.2f}, z={z_cl:.3f})")

prs_result = {
    'traits': traits_out,
    'physical_traits': phys_traits,
    'method': (f'GWAS-based PRS using GLIMPSE2-imputed posterior dosages. '
               f'Parker et al. 2017 reference panel ({n_snps} SNPs, {len(np.unique(breeds_fam))} breeds). '
               f'{genotyped} Parker positions covered.'),
    'reference': 'AKC breed trait scores (kkakey/dog_traits_AKC)',
    'snps': n_snps,
    'snps_imputed': genotyped,
    'n_ref_breeds': int(len(np.unique(breeds_fam))),
}
with open(f'{PUB}/prs_result.json', 'w') as f:
    json.dump(prs_result, f, indent=2)
print(f"prs_result.json written ({len(traits_out)} traits)")
PYEOF

# ── Stage 13: Inbreeding (Dog10K ROH + F distribution) ──────
log "=== Stage 12: Inbreeding (Dog10K) ==="
python3 - << PYEOF
import subprocess, json, numpy as np, gzip, re

BCF     = "$IMPUTED_BCF"
AF_FILE = "$D/COSMO/glimpse2_dog10k/het_out/panel_af.tsv.gz"
HET_FILE= "$D/COSMO/glimpse2_dog10k/het_out/dog10k_het.het"
PUB     = "$PUB"
DOG     = "$DOG_LOWER"
AUTOSOMES = [str(i) for i in range(1,39)]

# ── ROH-based FROH (sliding window on imputed BCF) ───────────
print("Computing FROH from imputed BCF...")
WINDOW = 50; MAX_HET = 1; MIN_ROH_KB = 500; MAX_GAP_KB = 1000
AUTO_MB = 2200.0

roh_total_mb = 0.0
roh_segments = []

for chrom_num in AUTOSOMES:
    chrom = f'chr{chrom_num}'
    result = subprocess.run(
        ['bcftools', 'query', '-r', chrom, '-f', '[%GT]\t%POS\n', BCF],
        capture_output=True, text=True)
    sites = []
    for line in result.stdout.strip().split('\n'):
        if not line: continue
        parts = line.split('\t')
        if len(parts) < 2: continue
        gt, pos_s = parts[0], parts[1]
        alleles = re.split(r'[|/]', gt)
        if '.' in alleles: continue
        is_het = len(set(alleles)) > 1
        sites.append((int(pos_s), is_het))
    if len(sites) < WINDOW: continue

    i = 0
    while i <= len(sites) - WINDOW:
        window = sites[i:i+WINDOW]
        het_count = sum(1 for _,h in window if h)
        if het_count <= MAX_HET:
            j = i + WINDOW
            while j < len(sites):
                if not sites[j][1]:
                    j += 1
                elif sum(1 for _,h in sites[i:j+1] if h) <= MAX_HET:
                    j += 1
                else:
                    break
            roh_start = sites[i][0]; roh_end = sites[j-1][0]
            roh_kb = (roh_end - roh_start) / 1000
            if roh_kb >= MIN_ROH_KB:
                roh_total_mb += roh_kb / 1000
                roh_segments.append({'chrom': chrom, 'start': roh_start,
                    'end': roh_end, 'length_mb': round(roh_kb/1000,3)})
            i = j
        else:
            i += 1

froh = roh_total_mb / AUTO_MB
level = ('Very Low' if froh < 0.03125 else 'Low' if froh < 0.0625 else
         'Moderate' if froh < 0.125 else 'High' if froh < 0.25 else 'Very High')
print(f"FROH={froh:.4f} ({level}), {roh_total_mb:.1f} Mb in {len(roh_segments)} ROH segments")

inbreeding = {
    'f_roh': round(froh, 4), 'f_roh_pct': round(froh*100, 2),
    'roh_total_mb': round(roh_total_mb, 2), 'roh_n_segments': len(roh_segments),
    'level': level.replace(' ','_'),
    'autosomal_genome_mb': AUTO_MB,
    'roh_segments': roh_segments,
    'method': 'ROH from Dog10K GLIMPSE2-imputed panel; window=50 SNPs, max_het=1, min_roh=500kb'
}
with open(f'{PUB}/inbreeding_result.json', 'w') as f:
    json.dump(inbreeding, f, indent=2)
print("inbreeding_result.json written")

# ── Genotype F vs Dog10K distribution ────────────────────────
print("Computing genotype F vs Dog10K distribution...")
f_vals = []
with open(HET_FILE) as fh:
    next(fh)
    for line in fh:
        parts = line.split()
        if len(parts) > 4:
            f_vals.append(float(parts[4]))
f_arr = np.array(f_vals)

regions_arg = ','.join(f'chr{c}' for c in AUTOSOMES)
cosmo_proc = subprocess.Popen(
    ['bcftools', 'query', '-r', regions_arg, '-f', '%CHROM\t%POS\t[%GT]\n', BCF],
    stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)

panel_fh = gzip.open(AF_FILE, 'rt')

def parse_af(s):
    v = s.split(',')[0]
    return None if v == '.' else float(v)

def next_panel():
    line = panel_fh.readline()
    if not line: return None, None, None
    parts = line.rstrip().split('\t')
    return parts[0], int(parts[1]), parse_af(parts[2])

o_hom = 0; e_hom = 0.0; n_obs = 0
p_chrom, p_pos, p_af = next_panel()

for line in cosmo_proc.stdout:
    parts = line.rstrip().split('\t')
    c_chrom, c_pos, gt = parts[0], int(parts[1]), parts[2]
    c_num = int(c_chrom.replace('chr',''))
    while p_chrom is not None:
        if p_chrom == c_chrom and p_pos == c_pos: break
        p_num = int(p_chrom.replace('chr',''))
        if p_num < c_num or (p_chrom == c_chrom and p_pos < c_pos):
            p_chrom, p_pos, p_af = next_panel()
        else: break
    if p_chrom != c_chrom or p_pos != c_pos or p_af is None: continue
    af = p_af; q = 1.0 - af
    is_hom = gt in ('0|0','1|1','0/0','1/1')
    o_hom += int(is_hom); e_hom += af*af + q*q; n_obs += 1
    p_chrom, p_pos, p_af = next_panel()

cosmo_proc.wait(); panel_fh.close()
if n_obs == 0 or (n_obs - e_hom) == 0:
    print("WARNING: no sites matched for genotype F")
else:
    dog_F = (o_hom - e_hom) / (n_obs - e_hom)
    pct = float(np.mean(f_arr < dog_F) * 100)
    f_min, f_max = float(f_arr.min()), float(f_arr.max())
    hist_counts, hist_edges = np.histogram(f_arr, bins=40, range=(f_min, f_max+0.001))
    dog10k_dist = {
        'sample_froh': round(dog_F, 4),
        'sample_percentile': round(pct, 1),
        'n_samples': len(f_arr),
        'ref_froh_mean': round(float(f_arr.mean()), 4),
        'ref_froh_p25': round(float(np.percentile(f_arr, 25)), 4),
        'ref_froh_p50': round(float(np.percentile(f_arr, 50)), 4),
        'ref_froh_p75': round(float(np.percentile(f_arr, 75)), 4),
        'hist_counts': hist_counts.tolist(),
        'hist_edges': hist_edges.tolist(),
        'metric': 'genotype_F',
        'note': f'Genotype F from {n_obs/1e6:.1f}M autosomal SNPs vs Dog10K panel ({len(f_arr)} dogs).'
    }
    with open(f'{PUB}/inbreeding_froh_dog10k_result.json', 'w') as f:
        json.dump(dog10k_dist, f, indent=2)
    print(f"inbreeding_froh_dog10k_result.json: F={dog_F:.4f} ({pct:.1f}th pct)")
PYEOF

# ── Stage 13: Coat color (GLIMPSE2 imputed genotypes at causal loci) ─────
log "=== Stage 13: Coat color ==="
export IMPUTED_BCF MARKDUP_BAM="$OUT/markdup.bam" PUB DOG_LOWER
python3 - << 'PYEOF'
import subprocess, json, pysam, re, tempfile, os

BCF     = os.environ['IMPUTED_BCF']
BAM     = os.environ['MARKDUP_BAM']
PUB     = os.environ['PUB']
MIN_GP  = 0.80   # min max(GP) to trust a GLIMPSE2 call

# ── Causal variant table ──────────────────────────────────────────────────
# exp_ref/exp_alt: expected REF/ALT in canFam4; if swapped in BCF, n_alt is flipped.
# inheritance: 'recessive' (need 2 copies) | 'dominant' (1 copy sufficient)
KNOWN_VARIANTS = [
    # E locus: MC1R (chr5)
    # e1: p.Arg306* (c.916C>T) — major recessive-red allele (loss-of-function)
    #
    # NOTE: chr5:64186854 is NOT the MC1R p.Arg306* coding position in canFam4 (ROS_Cfam_1.0).
    # snpEff annotates it as an intron_variant in CPNE7 (c.1030-385G>A), located ~2–3 kb
    # downstream of the DPEP1–CPNE7 intergenic gap where MC1R sits in the canFam4 assembly.
    # The canFam4 annotation omits MC1R as a separate gene entry; the actual causal position
    # (c.916C>T p.Arg306*) falls in that gap (chr5:~64155000–64184000) and is NOT genotyped
    # by name in the Dog10K panel.
    #
    # This position (64186854) is in the Dog10K panel at T allele freq ≈28.6% and is in strong
    # LD with the e haplotype on most chromosomes — it acts as a PROXY SNP. When a dog carries
    # the T allele at this CPNE7 intronic site on both chromosomes without carrying the actual
    # causal MC1R mutation (e.g., a chocolate b/b dog on an e-haplotype background), the pipeline
    # incorrectly calls e/e. See call_locus() for the confidence-downgrade logic.
    #
    # TODO: identify the correct canFam4 coordinate for MC1R p.Arg306* from NCBI ROS_Cfam_1.0
    # annotation release (gene ID: LOC403803) and add it as a second E-locus entry. Until then,
    # all e/e calls from this position alone are reported at low confidence with a validation warning.
    # Actual MC1R p.Arg306* coding position in canFam4 (ROS_Cfam_1.0, Gene ID 489652).
    # MC1R is on the minus strand at chr5:63922271–63923224 (CDS). Position c.916 in the CDS
    # (codon 306, AGA→TGA = p.Arg306*) maps to genomic chr5:63922309. On the plus strand the
    # reference allele is T (= A on the coding/minus strand = Arg) and the e allele is A
    # (= T on the coding strand = TGA stop). This site is NOT in the Dog10K SNP panel; the
    # pipeline will fall back to BAM pileup. At low-pass depth (<5 quality reads) it reports
    # "not in panel / no BAM reads" which is expected for 2–6× samples.
    dict(locus='E', chrom='chr5', pos=63922309, exp_ref='T', exp_alt='A',
         allele='e', inheritance='recessive',
         effect='MC1R p.Arg306* — causal e allele (c.916A>T on coding strand, AGA→TGA stop; canFam4 chr5:63922309)'),
    # Proxy SNP: CPNE7 intron variant in LD with the e haplotype (in Dog10K panel at AF≈29%).
    # NOT the MC1R coding variant. Used as an imputation proxy because the actual position above
    # is absent from the Dog10K reference panel. e/e calls from this site alone are low-confidence.
    dict(locus='E', chrom='chr5', pos=64186854, exp_ref='C', exp_alt='T',
         allele='e', inheritance='recessive',
         effect='e haplotype proxy SNP (CPNE7 intron, chr5:64186854) — in LD with e allele in Dog10K panel but NOT the causal MC1R p.Arg306* coding variant'),
    # Em: melanistic mask — dominant; tagging SNP at chr5:64188070
    dict(locus='E', chrom='chr5', pos=64188070, exp_ref=None, exp_alt=None,
         allele='Em', inheritance='dominant', effect='MC1R — melanistic mask haplotype tag'),

    # K locus: CBD103 (chr16)
    # KB: p.Lys43Arg (c.128A>G) — dominant black
    dict(locus='K', chrom='chr16', pos=57074438, exp_ref='A', exp_alt='G',
         allele='KB', inheritance='dominant', effect='CBD103 p.Lys43Arg — dominant black'),
    dict(locus='K', chrom='chr16', pos=57036106, exp_ref=None, exp_alt=None,
         allele='KB', inheritance='dominant', effect='CBD103 — dominant black tagging SNP 2'),

    # A locus: ASIP (chr24)
    # ay (sable) and aw involve regulatory/structural variants — not detectable from SNP imputation.
    # These positions tag at/a coding alleles only.
    dict(locus='A', chrom='chr24', pos=23906214, exp_ref=None, exp_alt=None,
         allele='at_tag', inheritance='recessive', effect='ASIP — tan-points/recessive-black coding tag'),
    dict(locus='A', chrom='chr24', pos=23908000, exp_ref=None, exp_alt=None,
         allele='at_tag', inheritance='recessive', effect='ASIP — tan-points/recessive-black coding tag 2'),

    # B locus: TYRP1 (chr11)
    # b1: p.Arg345Cys (c.1033C>T); b2: p.Gln354* (c.1060C>T) — both recessive brown
    # b1 allele: p.Arg345Cys. BCF encodes on the + strand as ref=T, alt=A (AF≈8% in Dog10K panel).
    dict(locus='B', chrom='chr11', pos=33376317, exp_ref='T', exp_alt='A',
         allele='b', inheritance='recessive', effect='TYRP1 p.Arg345Cys (b1) — brown/liver'),
    # b2 allele: p.Gln354*. BCF encodes on the + strand as ref=C, alt=T (AF≈53% in Dog10K panel).
    dict(locus='B', chrom='chr11', pos=33440938, exp_ref='C', exp_alt='T',
         allele='b', inheritance='recessive', effect='TYRP1 p.Gln354* (b2) — brown/liver'),

    # D locus: MLPH (chr25)
    # d1: splice site c.123+1G>A — recessive dilute
    dict(locus='D', chrom='chr25', pos=48403161, exp_ref='G', exp_alt='A',
         allele='d', inheritance='recessive', effect='MLPH c.123+1G>A splice site — dilute (blue/isabella)'),
    dict(locus='D', chrom='chr25', pos=48431759, exp_ref=None, exp_alt=None,
         allele='d', inheritance='recessive', effect='MLPH — dilute tagging SNP 2'),

    # S locus: MITF (chr20) — piebald spotting
    dict(locus='S', chrom='chr20', pos=5711695, exp_ref=None, exp_alt=None,
         allele='sp', inheritance='recessive', effect='MITF — piebald white spotting'),
    # M (PMEL SINE insertion) and W (KIT structural) not callable from SNP imputation
]

ALLELES_REFERENCE = {
    'E': {'Em': 'Melanistic mask (dominant)', 'E': 'Wild type extension',
          'e':  'Recessive red/yellow — two copies needed'},
    'K': {'KB':  'Dominant black — one copy = solid black',
          'kbr': 'Brindle (incompletely dominant)',
          'ky':  'Non-black/agouti — A locus determines pattern'},
    'A': {'ay': 'Sable/fawn (dominant, regulatory variant)',
          'aw': 'Wild type agouti', 'at': 'Tan points / tricolor (recessive)',
          'a':  'Recessive black (recessive)'},
    'B': {'B': 'Black eumelanin (dominant)', 'b': 'Brown/liver eumelanin — two copies needed'},
    'D': {'D': 'Full pigment (dominant)', 'd': 'Dilute/blue — two copies needed'},
    'M': {'M': 'Merle (dominant, PMEL SINE insertion — not detectable from SNP data)', 'm': 'Non-merle (assumed)'},
    'S': {'S': 'Solid / minimal white', 'sp': 'Piebald spotting (recessive)', 'sw': 'Extreme white (recessive)'},
    'W': {'w': 'Non-white', 'W': 'Extreme white (dominant, KIT structural — not detectable from SNP data)'},
}

LOCUS_INFO = {
    'E': dict(gene='MC1R',  chrom='chr5',  name='Extension locus',
              role='Master pigment switch: eumelanin (black/brown) vs phaeomelanin (yellow/red)',
              phenotype_contribution='e/e → all coat pigment is yellow/red regardless of other loci'),
    'K': dict(gene='CBD103', chrom='chr16', name='Dominant black locus',
              role='KB locks melanocytes in eumelanin production, overriding the A locus',
              phenotype_contribution='KB/- → solid eumelanin; ky/ky → A locus controls patterning'),
    'A': dict(gene='ASIP',  chrom='chr24', name='Agouti locus',
              role='Controls eumelanin/phaeomelanin switching within the hair shaft',
              phenotype_contribution='Only expressed when ky/ky at K locus; determines sable/tan-points/solid pattern'),
    'B': dict(gene='TYRP1', chrom='chr11', name='Brown locus',
              role='Modifies eumelanin color: B → black, b/b → brown/liver/chocolate',
              phenotype_contribution='b/b converts all black pigment to brown; no effect on phaeomelanin'),
    'D': dict(gene='MLPH',  chrom='chr25', name='Dilution locus',
              role='Melanosome transport: d/d dilutes pigment (black → blue, brown → isabella)',
              phenotype_contribution='d/d lightens all eumelanin; phaeomelanin unaffected'),
    'M': dict(gene='PMEL',  chrom='chr10', name='Merle locus',
              role='SINE insertion causes mosaic pigment dilution producing merle pattern',
              phenotype_contribution='Not detectable from SNP imputation — requires PCR or long-read'),
    'S': dict(gene='MITF',  chrom='chr20', name='Spotting locus',
              role='Controls melanocyte migration extent → white spotting area',
              phenotype_contribution='sp/sp → piebald; limited resolution from single SNP'),
    'W': dict(gene='KIT',   chrom='chr13', name='White locus',
              role='Extreme white spotting; dominant W linked to deafness risk',
              phenotype_contribution='Not detectable from short-read SNP data'),
}

# ── Batch BCF query ───────────────────────────────────────────────────────
bed = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
for v in KNOWN_VARIANTS:
    bed.write(f"{v['chrom']}\t{v['pos']-1}\t{v['pos']}\n")
bed.close()

res = subprocess.run(
    ['bcftools', 'query', '-R', bed.name,
     '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RAF\t[%GT]\t[%GP]\n', BCF],
    capture_output=True, text=True)
os.unlink(bed.name)

bcf_hits = {}
for line in res.stdout.strip().split('\n'):
    if not line: continue
    p = line.split('\t')
    if len(p) < 6: continue
    chrom, pos, ref, alt, raf_s, gt_s = p[0], int(p[1]), p[2], p[3], p[4], p[5]
    gp_s = p[6] if len(p) > 6 else ''
    try:    gp = [float(x) for x in gp_s.split(',')]; max_gp = max(gp)
    except: gp = None; max_gp = 0.0
    try:    raf = float(raf_s.split(',')[0])
    except: raf = None
    n_alt = sum(1 for a in re.split(r'[|/]', gt_s) if a == '1')
    bcf_hits[(chrom, pos)] = dict(ref=ref, alt=alt, gt=gt_s, n_alt=n_alt,
                                   raf=raf, gp=gp, max_gp=max_gp)

print(f"BCF hits: {len(bcf_hits)}/{len(KNOWN_VARIANTS)} positions in Dog10K panel")

# ── BAM pileup fallback ───────────────────────────────────────────────────
def bam_pileup(chrom, pos):
    try:
        bam = pysam.AlignmentFile(BAM, 'rb')
        counts = {}
        for col in bam.pileup(chrom, pos-1, pos, truncate=True,
                               min_base_quality=20, min_mapping_quality=20,
                               ignore_overlaps=True, ignore_orphans=True):
            if col.reference_pos != pos - 1: continue
            for r in col.pileups:
                if not r.is_del and not r.is_refskip:
                    b = r.alignment.query_sequence[r.query_position].upper()
                    counts[b] = counts.get(b, 0) + 1
        bam.close()
        return counts if sum(counts.values()) >= 5 else None
    except Exception:
        return None

# ── Per-variant calling ───────────────────────────────────────────────────
variant_calls = []
for v in KNOWN_VARIANTS:
    key = (v['chrom'], v['pos'])
    hit = bcf_hits.get(key)
    if hit and hit['max_gp'] >= MIN_GP:
        ref, alt, n_alt = hit['ref'], hit['alt'], hit['n_alt']
        # Flip n_alt if BCF orientation is swapped vs expectation
        if v['exp_ref'] and v['exp_alt'] and ref == v['exp_alt'] and alt == v['exp_ref']:
            n_alt = 2 - n_alt
            ref, alt = v['exp_ref'], v['exp_alt']
        variant_calls.append({**v, 'found': True, 'source': 'Dog10K imputed',
            'n_alt': n_alt, 'ref': ref, 'alt': alt, 'gt': hit['gt'],
            'gp': hit['gp'], 'max_gp': hit['max_gp'], 'raf': hit['raf'],
            'conf': 'high' if hit['max_gp'] >= 0.90 else 'medium'})
    elif hit:
        variant_calls.append({**v, 'found': True, 'source': 'Dog10K imputed (low GP)',
            'n_alt': hit['n_alt'], 'ref': hit['ref'], 'alt': hit['alt'],
            'gt': hit['gt'], 'gp': hit['gp'], 'max_gp': hit['max_gp'],
            'raf': hit['raf'], 'conf': 'low'})
    else:
        counts = bam_pileup(v['chrom'], v['pos'])
        if counts:
            total = sum(counts.values())
            variant_calls.append({**v, 'found': True, 'source': f'BAM pileup ({total} reads)',
                'n_alt': None, 'bam_counts': counts, 'total_reads': total,
                'conf': 'medium' if total >= 15 else 'low'})
        else:
            variant_calls.append({**v, 'found': False, 'source': 'not in panel / no BAM reads',
                'n_alt': None, 'conf': 'none'})

# ── Per-locus diploid genotype calling ───────────────────────────────────
def n_copies(locus, allele, calls):
    """Max ALT copies seen across all variants for this locus+allele."""
    hits = [c for c in calls if c['locus'] == locus and c['allele'] == allele
            and c['found'] and c['n_alt'] is not None]
    return max((c['n_alt'] for c in hits), default=None)

def any_found(locus, calls):
    return any(c['locus'] == locus and c['found'] for c in calls)

def call_locus(locus, calls):
    """Returns (allele1, allele2, confidence, interpretation)."""

    if locus == 'E':
        n_e  = n_copies('E', 'e',  calls)
        n_em = n_copies('E', 'Em', calls)
        if not any_found('E', calls):
            return '?', '?', 'low', 'MC1R positions not found in Dog10K panel'
        if n_e == 2:
            # Determine whether the homozygous call is supported by the actual MC1R coding position
            # (chr5:63922309) or only by the proxy SNP (chr5:64186854, CPNE7 intron).
            ACTUAL_MC1R_POS = 63922309
            PROXY_POS       = 64186854
            e_hom_calls = [c for c in calls if c['locus'] == 'E' and c['allele'] == 'e'
                           and c['found'] and c.get('n_alt') == 2]
            # Check whether the actual coding position has any alt read evidence from BAM pileup.
            # bam_pileup returns counts when ≥5 quality reads are found; otherwise n_alt=None.
            actual_bam = next((c for c in calls if c['locus'] == 'E' and c['allele'] == 'e'
                               and c.get('pos') == ACTUAL_MC1R_POS and c['found']), None)
            actual_has_alt = (actual_bam is not None and
                              actual_bam.get('bam_counts', {}).get('A', 0) > 0)
            proxy_only = bool(e_hom_calls) and all(c['pos'] == PROXY_POS for c in e_hom_calls)
            if proxy_only and not actual_has_alt:
                bam_depth = (actual_bam.get('bam_counts', {}) if actual_bam else {})
                depth_str = f'{sum(bam_depth.values())} reads, all reference' if bam_depth else 'no reads'
                return ('e', 'e', 'low',
                    f'Homozygous e/e from proxy SNP only (chr5:64186854, CPNE7 intron, alt allele freq ≈29%). '
                    f'The actual MC1R p.Arg306* coding position (chr5:63922309) is not in the Dog10K panel; '
                    f'BAM pileup shows {depth_str} — insufficient coverage to confirm or exclude the causal allele. '
                    f'If this dog shows eumelanin pigmentation (black/brown coat), the proxy call is likely a false positive '
                    f'and the true genotype is E/e or E/E. If the dog is cream/yellow, e/e remains plausible.')
            if proxy_only and actual_has_alt:
                return ('e', 'e', 'medium',
                    'Homozygous e/e: proxy SNP homozygous (chr5:64186854) AND alt reads detected at actual '
                    'MC1R coding position (chr5:63922309) from BAM pileup — consistent with e/e but low '
                    'BAM depth limits confidence.')
            return 'e', 'e', 'high', 'Homozygous recessive red — all coat pigment is phaeomelanin (yellow/red/cream)'
        if n_e == 1 and n_em and n_em >= 1:
            return 'Em', 'e', 'medium', 'Melanistic mask carrier for recessive red (Em/e)'
        if n_e == 1:
            return 'E', 'e', 'medium', 'Carrier for recessive red (E/e) — normal extension expressed'
        if n_em and n_em >= 1:
            return 'Em', 'E', 'medium', 'Melanistic mask on wild-type extension background (Em/E)'
        return 'E', 'E', 'medium', 'No e or Em alleles detected — wild-type extension (E/E)'

    elif locus == 'K':
        n_kb = n_copies('K', 'KB', calls)
        if not any_found('K', calls):
            return '?', '?', 'low', 'CBD103 positions not found in Dog10K panel'
        if n_kb == 2:
            return 'KB', 'KB', 'high', 'Homozygous dominant black (KB/KB)'
        if n_kb == 1:
            return 'KB', 'ky', 'high', 'Dominant black carrier (KB/ky) — KB overrides A locus'
        return 'ky', 'ky', 'high', 'No KB allele (ky/ky) — A locus controls pattern'

    elif locus == 'A':
        at_calls = [c for c in calls if c['locus'] == 'A' and c['allele'] == 'at_tag'
                    and c['found'] and c['n_alt'] is not None]
        if not any_found('A', calls):
            return '?', '?', 'low', 'ASIP coding positions not found; ay/aw require structural variant analysis'
        n_at = max((c['n_alt'] for c in at_calls), default=0)
        n_het_sites = sum(1 for c in at_calls if c['n_alt'] == 1)
        any_hom = any(c['n_alt'] == 2 for c in at_calls)
        if any_hom or n_het_sites >= 2:
            return 'at', 'at', 'medium', 'Tan-points / tricolor (at/at) — recessive black or tricolor pattern'
        if n_at == 1:
            return 'ay/?', 'at', 'low', 'One putative at allele detected; other allele uncertain (ay/at or at/aw possible)'
        return 'ay/?', 'ay/?', 'low', ('No at/a coding variants detected. '
            'Sable (ay) or wild agouti (aw) likely but require structural variant analysis to confirm.')

    elif locus == 'B':
        b_by_pos = {c['pos']: c['n_alt'] for c in calls
                    if c['locus'] == 'B' and c['allele'] == 'b'
                    and c['found'] and c['n_alt'] is not None}
        if not b_by_pos:
            return '?', '?', 'low', 'TYRP1 brown alleles not found in Dog10K panel'
        any_hom = any(n == 2 for n in b_by_pos.values())
        n_het_sites = sum(1 for n in b_by_pos.values() if n == 1)
        if any_hom or n_het_sites >= 2:
            return 'b', 'b', 'high', 'Brown/liver eumelanin (b/b or compound b1/b2)'
        if n_het_sites == 1:
            # One b allele detected. The Dog10K panel only covers b1 (p.Arg345Cys) and
            # b2 (p.Gln354*). A dog that appears chocolate may carry b1/b3 or b1/b4
            # compound het where the second allele is absent from the panel. Cannot
            # distinguish B/b carrier from b/b compound het from panel data alone.
            return 'b', '?', 'low', ('One b allele detected (heterozygous). The Dog10K panel '
                'covers b1 (p.Arg345Cys) and b2 (p.Gln354*) only — a second b allele at '
                'an unqueried position (b3, b4, or other TYRP1 variant) cannot be excluded. '
                'Genotype is B/b (carrier, black eumelanin) OR b/b compound het (chocolate) '
                'if a second allele is present outside the panel.')
        return 'B', 'B', 'high', 'No brown alleles detected (B/B)'

    elif locus == 'D':
        d_by_pos = {c['pos']: c['n_alt'] for c in calls
                    if c['locus'] == 'D' and c['allele'] == 'd'
                    and c['found'] and c['n_alt'] is not None}
        if not d_by_pos:
            return '?', '?', 'low', 'MLPH dilute alleles not found in Dog10K panel'
        any_hom = any(n == 2 for n in d_by_pos.values())
        n_het_sites = sum(1 for n in d_by_pos.values() if n == 1)
        if any_hom or n_het_sites >= 2:
            return 'd', 'd', 'high', 'Dilute coat (d/d) — black→blue, brown→isabella'
        if n_het_sites == 1:
            return 'D', 'd', 'high', 'Carrier for dilute (D/d) — full pigment expressed'
        return 'D', 'D', 'high', 'No dilute alleles detected (D/D) — full pigment'

    elif locus == 'M':
        return 'm', 'm', 'low', 'Merle (PMEL SINE insertion) not detectable from SNP imputation — PCR required'

    elif locus == 'S':
        n_sp = n_copies('S', 'sp', calls)
        if not any_found('S', calls):
            return '?', '?', 'low', 'MITF spotting variant not found in Dog10K panel'
        if n_sp == 2:
            return 'sp', 'sp', 'medium', 'Piebald spotting (sp/sp) — white markings expected'
        if n_sp == 1:
            return 'S',  'sp', 'medium', 'Carrier for piebald (S/sp) — minimal or no white markings'
        return 'S', 'S', 'medium', 'No piebald allele at MITF queried position'

    elif locus == 'W':
        return 'w', 'w', 'low', 'KIT extreme white (structural variant) not detectable from SNP data'

    return '?', '?', 'low', 'Unknown locus'

loci_gt = {}
for locus in ['E', 'K', 'A', 'B', 'D', 'M', 'S', 'W']:
    a1, a2, conf, interp = call_locus(locus, variant_calls)
    loci_gt[locus] = dict(allele1=a1, allele2=a2, confidence=conf, interpretation=interp)

# ── Cross-locus E locus validation ────────────────────────────────────────
# When the E locus is called e/e from the proxy SNP only (low confidence),
# use K and B to compute what the coat would be if the proxy is a false positive.
# The proxy SNP (chr5:64186854, CPNE7 intron) tags the e haplotype but is not the
# causal MC1R coding variant; false positives are known (e.g. chocolate dogs).
# Knowing the K+B "eumelanic" prediction helps interpret ambiguity.
e_gt = loci_gt['E']
if (e_gt['allele1'] == 'e' and e_gt['allele2'] == 'e'
        and e_gt['confidence'] == 'low'):
    k_gt = loci_gt['K']
    b_gt = loci_gt['B']
    has_KB = 'KB' in (k_gt['allele1'], k_gt['allele2'])
    is_bb  = b_gt['allele1'] == 'b' and b_gt['allele2'] == 'b'
    if has_KB and is_bb:
        alt_color = 'chocolate'
    elif has_KB:
        alt_color = 'black'
    elif k_gt['allele1'] in ('ky', '?') and is_bb:
        alt_color = 'chocolate or sable (A locus)'
    else:
        alt_color = 'black or sable (A locus)'
    updated_interp = (
        e_gt['interpretation'] +
        f' If the proxy call is a false positive, the coat would be {alt_color} '
        f'based on K ({k_gt["allele1"]}/{k_gt["allele2"]}) and '
        f'B ({b_gt["allele1"]}/{b_gt["allele2"]}) loci. '
        f'Phenotype confirmation is required to distinguish e/e (cream/yellow) '
        f'from {alt_color}.'
    )
    loci_gt['E'] = {**e_gt, 'interpretation': updated_interp,
                    'proxy_false_positive_prediction': alt_color}

# ── Phenotype prediction: hierarchical epistasis protocol ─────────────────
# Implements the standard five-locus diagnostic hierarchy:
#   Step 2: B + D  → eumelanin base pigment (black / chocolate / blue / isabella)
#   Step 3 Tier 1: E locus  → e/e = cream/yellow, terminate
#   Step 3 Tier 2: K locus  → KB = solid eumelanin, terminate; kbr = brindle modifier
#   Step 3 Tier 3: A locus  → sable / agouti / tan-points / recessive black

def predict_phenotype(loci_gt):
    e1, e2 = loci_gt['E']['allele1'], loci_gt['E']['allele2']
    k1, k2 = loci_gt['K']['allele1'], loci_gt['K']['allele2']
    a1, a2 = loci_gt['A']['allele1'], loci_gt['A']['allele2']
    b1, b2 = loci_gt['B']['allele1'], loci_gt['B']['allele2']
    d1, d2 = loci_gt['D']['allele1'], loci_gt['D']['allele2']

    # Step 2: eumelanin base pigment (B + D)
    # b/? means one b allele confirmed, second allele unknown (could be B or another b).
    is_bb    = b1 == 'b' and b2 == 'b'
    is_b_unk = (b1 == 'b' and b2 == '?') or (b1 == '?' and b2 == 'b')
    is_dd    = d1 == 'd' and d2 == 'd'
    if is_bb and is_dd:
        eume_color = 'isabella/lilac'
        nose_color = 'isabella/lilac nose and pads'
    elif is_bb:
        eume_color = 'chocolate'
        nose_color = 'liver/brown nose and pads'
    elif is_b_unk and is_dd:
        eume_color = 'isabella/lilac or blue/grey'
        nose_color = 'isabella/lilac or blue/grey nose and pads (b allele status uncertain)'
    elif is_b_unk:
        eume_color = 'chocolate or black'
        nose_color = 'liver/brown or black nose and pads (b allele status uncertain)'
    elif is_dd:
        eume_color = 'blue/grey'
        nose_color = 'blue/grey nose and pads'
    else:
        eume_color = 'black'
        nose_color = 'black nose and pads'

    # Step 3 Tier 1: E locus — e/e = recessive red, overrides all other loci
    is_e_hom = e1 == 'e' and e2 == 'e'
    has_em   = 'Em' in (e1, e2)

    if is_e_hom:
        base_color = (f'Phaeomelanin — cream / yellow / red '
                      f'(e/e overrides K, A, B loci for coat; {nose_color} from B/D loci)')
        pattern = ('Solid phaeomelanin coat — no eumelanin in hair regardless of K or A locus. '
                   f'Skin pigment ({nose_color}) is determined by B and D loci independently.')
        dilution = 'Not applicable to coat (phaeomelanin unaffected by D locus)'
        return base_color, pattern, dilution

    # Step 3 Tier 2: K locus — KB = solid eumelanin
    has_KB  = 'KB' in (k1, k2)
    has_kbr = 'kbr' in (k1, k2)
    is_ky_hom = k1 == 'ky' and k2 == 'ky'
    mask_note = ' with melanistic mask (Em)' if has_em else ''

    if has_KB:
        base_color = f'Eumelanin — solid {eume_color}'
        pattern = (f'Solid {eume_color}{mask_note} — KB dominant black suppresses A locus entirely. '
                   f'{nose_color.capitalize()}.')
        dil_str = (f'Dilute (d/d) — {eume_color} coat' if is_dd
                   else f'Full pigment (D/D or D/d)')
        return base_color, pattern, dil_str

    # Step 3 Tier 3: A locus (reached only when ky/ky or kbr/ky)
    brindle = ' brindled' if has_kbr else ''

    has_Ay = 'ay' in (a1, a2)
    has_aw = 'aw' in (a1, a2)
    has_at = 'at' in (a1, a2)
    is_a   = a1 == 'a' and a2 == 'a'

    dil_str = (f'Dilute (d/d) — {eume_color} eumelanin' if is_dd
               else 'Full pigment (D/D or D/d)')

    if has_Ay:
        base_color = f'Eumelanin base — {eume_color}; phaeomelanin coat (sable/fawn)'
        pattern = (f'{eume_color.capitalize()}-based{brindle} sable/fawn{mask_note} — '
                   f'predominantly phaeomelanin (yellow/red/cream) coat with {eume_color}-tipped hairs. '
                   f'{nose_color.capitalize()}. '
                   f'Note: ay requires structural variant confirmation (not in Dog10K panel).')
    elif has_aw:
        base_color = f'Eumelanin base — {eume_color}; agouti banding (wolf sable)'
        pattern = (f'{eume_color.capitalize()}-based{brindle} wolf sable / agouti{mask_note} — '
                   f'individual hairs banded with alternating {eume_color} and phaeomelanin. '
                   f'{nose_color.capitalize()}.')
    elif has_at:
        base_color = f'Eumelanin — {eume_color} with tan points'
        pattern = (f'{eume_color.capitalize()}-based{brindle} tan points{mask_note} — '
                   f'{eume_color} body with phaeomelanin markings on muzzle, eyebrows, chest, '
                   f'inner ears, and lower legs. {nose_color.capitalize()}.')
    elif is_a:
        base_color = f'Eumelanin — solid {eume_color} (recessive black)'
        pattern = (f'Solid {eume_color}{mask_note} via recessive black (a/a) — '
                   f'A locus bypasses phaeomelanin expression entirely. {nose_color.capitalize()}.')
    else:
        base_color = f'Eumelanin base — {eume_color} (A locus undetermined)'
        pattern = (f'{eume_color.capitalize()}-based{brindle} coat{mask_note}; '
                   f'A locus pattern unknown (sable/agouti/tan-points require structural variant analysis). '
                   f'{nose_color.capitalize()}.')

    return base_color, pattern, dil_str

base_color, pattern, dilution = predict_phenotype(loci_gt)

# IRF4: check CNV data
try:
    with open(f'{PUB}/cnv_homdel.json') as _f:
        cnv = json.load(_f)
    irf4_dels = [g for r in cnv.get('regions', [])
                 for g in r.get('disrupted_genes', []) if 'IRF4' in g]
    if irf4_dels:
        irf4_note = ('IRF4 deletion detected — associated with progressive graying/silvering '
                     'of eumelanin pigment, particularly visible in dark-coated dogs')
    else:
        irf4_note = 'No IRF4 deletion detected in this sample.'
except Exception:
    irf4_note = 'IRF4 deletion status unknown (CNV data unavailable)'

overall_conf = ('medium'
    if all(loci_gt[l]['confidence'] in ('high', 'medium') for l in ['E', 'K', 'B', 'D'])
    else 'low')

e_gt = loci_gt['E']
validation_warning = None

# ── Build per-locus output ────────────────────────────────────────────────
loci_result = {}
for locus in ['E', 'K', 'A', 'B', 'D', 'M', 'S', 'W']:
    info = LOCUS_INFO[locus]
    g    = loci_gt[locus]
    obs  = []
    for vc in variant_calls:
        if vc['locus'] != locus or not vc['found']: continue
        ov = {'pos': vc['pos'], 'source': vc['source'], 'effect': vc['effect']}
        if vc.get('n_alt') is not None:
            ov.update(gt=vc.get('gt',''), n_alt=vc['n_alt'],
                      ref=vc.get('ref',''), alt=vc.get('alt',''))
        if vc.get('gp'):
            ov['gp'] = [round(x,3) for x in vc['gp']]
            ov['max_gp'] = round(vc['max_gp'], 3)
        if vc.get('raf') is not None:
            ov['af'] = round(1 - vc['raf'], 4)
        if vc.get('bam_counts'):
            ov['bam_counts'] = vc['bam_counts']
            ov['depth'] = vc['total_reads']
        obs.append(ov)

    locus_entry = {
        'gene': info['gene'], 'chrom': info['chrom'],
        'name': info['name'], 'role': info['role'],
        'phenotype_contribution': info['phenotype_contribution'],
        'alleles_reference': ALLELES_REFERENCE[locus],
        'predicted_alleles': [g['allele1'], g['allele2']],
        'confidence': g['confidence'],
        'interpretation': g['interpretation'],
        'observed_variants': obs,
    }
    if g.get('proxy_false_positive_prediction'):
        locus_entry['proxy_false_positive_prediction'] = g['proxy_false_positive_prediction']
    loci_result[locus] = locus_entry

coat = {
    'summary': {
        'predicted_base_color': base_color,
        'predicted_pattern': pattern,
        'predicted_dilution': dilution,
        'predicted_white': 'Not detectable from SNP data (S locus limited; W requires structural variant)',
        'predicted_merle': 'Not detectable from SNP imputation — requires PCR or long-read sequencing',
        'overall_confidence': overall_conf,
        **({'validation_warning': validation_warning} if validation_warning else {}),
        'caveat': ('E, K, B, D loci called from Dog10K GLIMPSE2 imputed BCF (causal SNPs). '
                   'A locus sable (ay/aw) requires structural variant analysis not available here. '
                   'Merle (M) and extreme white (W) require PCR or long-read. '
                   'Commercial tests (Embark, Wisdom Panel) cover additional alleles.'),
        'irf4_note': irf4_note,
    },
    'loci': loci_result,
    'method': (f'Coat color genotyping from GLIMPSE2 Dog10K imputed BCF (min GP={MIN_GP}). '
               'Causal SNPs queried at known canFam4 positions; BAM pileup fallback for sites not in panel. '
               'Compound heterozygosity handled for B (b1/b2) and D (d1/d2) loci.'),
}
with open(f'{PUB}/coat_color.json', 'w') as f:
    json.dump(coat, f, indent=2)
print('coat_color.json written')
for locus, g in loci_gt.items():
    print(f"  {locus}: {g['allele1']}/{g['allele2']} ({g['confidence']})")
PYEOF

fi # end stages 10–13

# ── Stage 15: Oral microbiome (MetaPhlAn4) ───────────────────
log "=== Stage 15: Oral microbiome (MetaPhlAn4) ==="
BAM_FOR_MICRO="$OUT/markdup.bam"
UNMAPPED_FQ="$OUT/${DOG_LOWER}_unmapped.fastq"
MICRO_OUT="$OUT/${DOG_LOWER}_metaphlan.txt"
MICRO_BT2="$OUT/${DOG_LOWER}_metaphlan.mapout.bz2"

    [[ -f "$BAM_FOR_MICRO" ]] || die "markdup.bam not found at $BAM_FOR_MICRO"
    [[ -f "$MICROBIOME_REF" ]] || die "Microbiome reference CSV not found at $MICROBIOME_REF"

    # Extract unmapped reads to avoid OOM on large BAMs
    log "  Extracting unmapped reads from BAM…"
    "$ENV_GENOMICS/bin/samtools" fastq -f 4 -@ 4 "$BAM_FOR_MICRO" > "$UNMAPPED_FQ"
    N_READS=$(wc -l < "$UNMAPPED_FQ")
    log "  Unmapped reads: $((N_READS/4)) ($(wc -c < "$UNMAPPED_FQ" | awk '{printf "%.1f", $1/1e6}') MB)"

    log "  Running MetaPhlAn4…"
    if [[ -f "$MICRO_BT2" ]]; then
        log "  Reusing existing mapout: $MICRO_BT2"
        "$METAPHLAN_BIN" "$MICRO_BT2" \
            --input_type mapout \
            --nproc 4 \
            -o "$MICRO_OUT"
    else
        "$METAPHLAN_BIN" "$UNMAPPED_FQ" \
            --input_type fastq \
            --mapout "$MICRO_BT2" \
            --nproc 4 \
            -o "$MICRO_OUT"
    fi

    log "  Computing microbiome JSONs…"
    python3 - << PYEOF
import json, re, math, datetime
import numpy as np
import pandas as pd
from scipy.stats import percentileofscore, entropy
from sklearn.linear_model import RidgeCV

PUB      = "$PUB"
OUT      = "$OUT"
DOG      = "$DOG_NAME"
MICRO_OUT = "$MICRO_OUT"
REF_CSV   = "$MICROBIOME_REF"
ACTUAL_AGE = "$DOG_ACTUAL_AGE"

# ── 1. Parse MetaPhlAn4 output ─────────────────────────────
RANK_MAP = {'k':'k','p':'p','c':'c','o':'o','f':'f','g':'g','s':'s','t':'t'}

def short_name(clade):
    parts = clade.split('|')
    last = parts[-1]
    return re.sub(r'^[a-z]__', '', last).replace('_', ' ')

taxa = {'kingdom':[],'phyla':[],'classes':[],'orders':[],'families':[],'genera':[],'species':[]}
rank_to_key = {'k':'kingdom','p':'phyla','c':'classes','o':'orders',
               'f':'families','g':'genera','s':'species'}
total_classified = 0.0

with open(MICRO_OUT) as fh:
    for line in fh:
        if line.startswith('#'): continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 2: continue
        clade = parts[0]
        # MetaPhlAn 4.2.4+ added NCBI_tax_id as col 1; abundance is col 2 if present, else col 1
        pct = float(parts[2]) if len(parts) >= 3 and parts[1].replace('|','').isdigit() else float(parts[1])
        ranks = [seg.split('__')[0] for seg in clade.split('|')]
        deepest = ranks[-1]
        if deepest not in rank_to_key: continue
        entry = {'clade': clade, 'rank': deepest,
                 'name': short_name(clade),
                 'relative_abundance': round(pct, 6),
                 'estimated_reads': None}
        taxa[rank_to_key[deepest]].append(entry)
        if deepest == 'k' and 'Bacteria' in clade:
            total_classified = pct

# Sort each rank by abundance descending
for k in taxa:
    taxa[k].sort(key=lambda x: -x['relative_abundance'])

n_species = len(taxa['species'])
n_reads_total = 0
with open("$OUT/pipeline.log") as lf:
    for line in lf:
        if 'reads_processed' in line.lower() or 'Total reads' in line:
            pass  # best-effort; leave as 0

# Scale all relative_abundance values to % of all reads
scale = total_classified / 100.0
for key in taxa:
    for entry in taxa[key]:
        entry['relative_abundance'] = round(entry['relative_abundance'] * scale, 6)

micro_result = {
    'sample': DOG.lower(),
    'run_date': datetime.date.today().isoformat(),
    'db_version': 'mpa_vJan25_CHOCOPhlAnSGB_202503',
    'total_classified_pct': round(total_classified, 4),
    **taxa,
}
with open(f'{PUB}/microbiome_result.json', 'w') as fh:
    json.dump(micro_result, fh, indent=2)
print(f"microbiome_result.json: {n_species} species, {total_classified:.2f}% classified")

# ── 2. Load reference dataset ──────────────────────────────
df = pd.read_csv(REF_CSV)
sp_cols   = [c for c in df.columns if '|s__' in c and '|t__' not in c]
prev      = (df[sp_cols] > 0).mean()
sp_filtered = prev[prev > 0.10].index.tolist()
print(f"Reference species features (>10% prevalence): {len(sp_filtered)}")

# ── 3. Kiki / dog species dict (% of classified bacteria) ──
kiki_species = {}
for sp in taxa['species']:
    kiki_species[sp['clade']] = round(sp['relative_abundance'] / scale, 6)

matched_features = [f for f in sp_filtered if f in kiki_species]
print(f"Matched features: {len(matched_features)}")

# ── 4. Age prediction (RidgeCV) ────────────────────────────
age_col = next((c for c in df.columns if c.lower() == 'age'), None)
if age_col is None:
    raise RuntimeError(f"No 'age' column found in {REF_CSV}. Columns: {df.columns.tolist()[:10]}")

X_ref = np.log10(df[sp_filtered].values + 1e-5)
y_ref = df[age_col].values

model = RidgeCV(alphas=[0.01,0.1,1,10,100], cv=5)
model.fit(X_ref, y_ref)

# Cross-val metrics
from sklearn.model_selection import cross_val_score
cv_r2  = cross_val_score(model, X_ref, y_ref, cv=5, scoring='r2').mean()
cv_mae = -cross_val_score(model, X_ref, y_ref, cv=5, scoring='neg_mean_absolute_error').mean()

# Predict for dog
kiki_vec = np.array([kiki_species.get(f, 0.0) for f in sp_filtered])
kiki_vec_log = np.log10(kiki_vec + 1e-5).reshape(1,-1)
pred_age = float(model.predict(kiki_vec_log)[0])

# Top species (by coefficient magnitude)
coef_pairs = sorted(zip(sp_filtered, model.coef_), key=lambda x: -abs(x[1]))[:10]
top_species = [{'name': re.sub(r'.*\|s__', '', f).replace('_', ' '), 'coefficient': round(c, 5)}
               for f, c in coef_pairs]

age_result = {
    'predicted_age_years': round(pred_age, 2),
    'cv_r2':               round(cv_r2,  3),
    'cv_mae_years':        round(cv_mae, 3),
    'n_training_samples':  len(df),
    'n_species_features':  len(sp_filtered),
    'n_features_matched': len(matched_features),
    'model':               'RidgeCV (log10-transformed, prevalence>10%)',
    'top_species':         top_species,
}
if ACTUAL_AGE:
    try:
        age_result['actual_age_years'] = float(ACTUAL_AGE)
    except ValueError:
        pass

with open(f'{PUB}/microbiome_age_result.json', 'w') as fh:
    json.dump(age_result, fh, indent=2)
print(f"microbiome_age_result.json: predicted={pred_age:.1f} yrs, matched={len(matched_features)} features")

# ── 5. Diversity & pathobiont health metrics ───────────────
# Diversity comparison uses genus-level aggregation because the reference was
# profiled with MetaPhlAn 3.18 (old NCBI taxonomy) while the sample uses
# MetaPhlAn4/Jan25 (GTDB taxonomy). Phylum/class/order names differ between
# versions, so full-clade species matching yields very few hits. Genus names
# are stable across versions and give representative percentiles.

# Build sample genus dict: sum species abundances per genus leaf name
sample_genus_abund = {}
for sp in taxa['species']:
    if '|g__' in sp['clade']:
        g_leaf = sp['clade'].split('|g__')[-1].split('|')[0]
        sample_genus_abund[g_leaf] = sample_genus_abund.get(g_leaf, 0.0) + sp['relative_abundance'] / scale

# Reference genus columns (leaf name → column); prevalence-filter at >10%
g_cols_all = [c for c in df.columns if '|g__' in c and '|s__' not in c and '|t__' not in c]
g_prev = (df[g_cols_all] > 0).mean()
g_filtered = g_prev[g_prev > 0.10].index.tolist()
g_leaf_map = {c: c.split('|g__')[-1] for c in g_filtered}   # col → leaf

matched_genera = [c for c in g_filtered if g_leaf_map[c] in sample_genus_abund]
print(f"Matched genera (>10% prevalence): {len(matched_genera)}")

kiki_abund = np.array([sample_genus_abund.get(g_leaf_map[c], 0.0) for c in matched_genera])
kiki_richness = int((kiki_abund > 0).sum())
nz = kiki_abund[kiki_abund > 0]
kiki_shannon = float(entropy(nz / nz.sum())) if len(nz) > 0 else 0.0

ref_richness_vec = (df[matched_genera] > 0).sum(axis=1).values
ref_shannon_vec  = np.array([
    entropy(row[row > 0] / row[row > 0].sum()) if (row > 0).any() else 0.0
    for row in df[matched_genera].values
])

richness_pct = round(percentileofscore(ref_richness_vec, kiki_richness, kind='rank'), 1)
shannon_pct  = round(percentileofscore(ref_shannon_vec,  kiki_shannon,  kind='rank'), 1)
r_p25, r_p50, r_p75 = [int(x) for x in np.percentile(ref_richness_vec, [25,50,75])]
s_p25, s_p50, s_p75 = [round(x,4) for x in np.percentile(ref_shannon_vec,  [25,50,75])]

# Pathobionts (key canine periodontal pathogens)
PATHOBIONTS = {
    's__Porphyromonas_gulae':         ('red',    'canine periodontal disease'),
    's__Tannerella_forsythia':         ('red',    'periodontal disease'),
    's__Porphyromonas_cangingivalis':  ('red',    'canine periodontitis'),
    's__Porphyromonas_canoris':        ('orange', 'canine oral disease'),
    's__Porphyromonas_gingivicanis':   ('red',    'canine periodontitis'),
    's__Treponema_denticola':          ('red',    'periodontal disease'),
    's__Fusobacterium_nucleatum':      ('orange', 'periodontal disease'),
    's__Prevotella_intermedia':        ('orange', 'periodontal disease'),
}

sp_pct_dict = {sp['clade']: sp['relative_abundance'] / scale
               for sp in taxa['species']}   # % of classified bacteria

pathobiont_hits = []
for clade_suffix, (color, assoc) in PATHOBIONTS.items():
    pct = 0.0
    for clade, v in sp_pct_dict.items():
        if clade_suffix in clade and '|t__' not in clade:
            pct += v
    if pct > 0:
        name_clean = clade_suffix.replace('s__','').replace('_',' ')
        pathobiont_hits.append({'name': name_clean, 'pct': round(pct,3),
                                'color': color, 'association': assoc})

pathobiont_hits.sort(key=lambda x: -x['pct'])
pathobiont_total = sum(h['pct'] for h in pathobiont_hits)
commensal_pct   = max(0.0, 100.0 - pathobiont_total)
dysbiosis_index = round(pathobiont_total / (commensal_pct + 1e-6), 3)

# Reference pathobiont distribution
path_cols = []
for clade_suffix in PATHOBIONTS:
    for c in df.columns:
        if clade_suffix in c and '|t__' not in c and '|s__' in c:
            path_cols.append(c)
            break
if path_cols:
    # Reference CSV is in fractions (0-1); convert to % to match pathobiont_total units
    ref_path_vec = df[path_cols].sum(axis=1).values * 100.0
else:
    ref_path_vec = np.zeros(len(df))
path_pct = round(percentileofscore(ref_path_vec, pathobiont_total, kind='rank'), 1)
r_pm, r_pmed, r_p75p, r_p90p = (round(float(x),2) for x in
    [ref_path_vec.mean(), np.median(ref_path_vec),
     np.percentile(ref_path_vec,75), np.percentile(ref_path_vec,90)])

health_result = {
    'sample_richness':         len(taxa['species']),
    'sample_shannon':          round(float(entropy(
        np.array([s['relative_abundance'] for s in taxa['species']]) /
        sum(s['relative_abundance'] for s in taxa['species'])
    )), 4) if taxa['species'] else 0.0,
    'sample_richness_matched': kiki_richness,
    'sample_shannon_matched':  round(kiki_shannon, 4),
    'richness_percentile':    richness_pct,
    'shannon_percentile':     shannon_pct,
    'ref_richness_p25':       r_p25,
    'ref_richness_p50':       r_p50,
    'ref_richness_p75':       r_p75,
    'ref_shannon_p25':        s_p25,
    'ref_shannon_p50':        s_p50,
    'ref_shannon_p75':        s_p75,
    'n_matched_genera':       len(matched_genera),
    'diversity_note':         'Genus-level comparison (ref: MetaPhlAn 3.18; sample: MetaPhlAn4 Jan25)',
    'pathobiont_burden_pct':  round(pathobiont_total, 3),
    'pathobiont_percentile':  path_pct,
    'commensal_pct':          round(commensal_pct, 3),
    'dysbiosis_index':        dysbiosis_index,
    'ref_pathobiont_mean':    r_pm,
    'ref_pathobiont_median':  r_pmed,
    'ref_pathobiont_p75':     r_p75p,
    'ref_pathobiont_p90':     r_p90p,
    'pathobiont_hits':        pathobiont_hits,
}
with open(f'{PUB}/microbiome_health_result.json', 'w') as fh:
    json.dump(health_result, fh, indent=2)
print(f"microbiome_health_result.json: genus_richness={kiki_richness}/{len(matched_genera)} ({richness_pct}th pct), "
      f"shannon={kiki_shannon:.3f} ({shannon_pct}th pct), pathobionts={pathobiont_total:.1f}%")
PYEOF

    log "  Microbiome stage complete."

# ── Stage 16: Copy reference JSONs ───────────────────────────
log "=== Stage 16: Copy reference JSONs ==="
for f in centromeres.json genes_1mb.json cnv_genes.json karyotype_zoom.json; do
    cp "$COSMO_PUB/$f" "$PUB/$f"
    log "  Copied $f"
done

# ── Stage 17: Commit & push to Vercel ───────────────────────
log "=== Stage 17: Commit and push ==="
cd "$D/dogs-app"
git add "public/$DOG_LOWER/"
git commit -m "Add $DOG_NAME genomic data (WGS pipeline)

Dog10K imputed BCF, VEP annotation, OMIA calls, breed prediction,
PRS, inbreeding, and coat color for $DOG_NAME.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
git push

log " Pipeline complete: $DOG_NAME"
log " Dashboard: dogs-app/public/$DOG_LOWER/"
echo "DONE" > "$OUT/pipeline.done"
