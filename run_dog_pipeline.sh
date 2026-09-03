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

# ── Site profile ─────────────────────────────────────────────
# One pipeline, several machines. A site profile supplies paths, tool locations
# and parallelism for wherever this is running. An explicit DOGS_SITE always
# wins; otherwise detect. $SGE_ROOT is set on Hoffman2 login and compute nodes.
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOGS_SITE="${DOGS_SITE:-$(
  if   [[ -n "${SGE_ROOT:-}" ]];      then echo hoffman
  elif [[ "$(uname -s)" == Darwin ]]; then echo mac
  else                                     echo generic
  fi)}"
SITE_FILE="$PIPELINE_DIR/site/${DOGS_SITE}.sh"
[[ -f "$SITE_FILE" ]] || { echo "ERROR: no site profile at $SITE_FILE (DOGS_SITE=$DOGS_SITE)"; exit 1; }
# shellcheck source=/dev/null
source "$SITE_FILE"

D="${D:-/Users/matteopellegrini/Downloads/dogs}"   # base dir for default paths

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
# Optional stop point. Re-running one stage matters now that intermediates are
# not all kept: a coverage re-analysis must start at stage 6 (the BAM stage 5
# needs is discarded by local-scratch runs) and must not continue into stage 7,
# which needs that same BAM. Default 99 = run everything.
TO_STAGE="${TO_STAGE:-99}"

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
REF_JSON=$D/reference_json                     # shared reference JSONs (genome annotations,
                                               # OMIA catalogue, PRS baseline); deliberately
                                               # outside public/ so no genomic data is
                                               # world-readable and the cluster needs no app checkout
OMIA_DB=$REF_JSON/omia_variants.json           # OMIA variant catalogue (479 variants)
# Merged Parker + Dog10K breed panel (see analysis/breed_accuracy/README.md).
# Built by make_merged_plink.py with harmonised labels, village dogs grouped by
# region, wolves pooled, min-n 6 -> 230 breeds over 131,353 SNPs.
# Coverage panels-of-normals built from 96 dogs by
# analysis/coverage_panel/build_panels.py. Thresholds are calibrated per scale
# against replicate agreement — see that directory's README.
COV_PANEL_1MB=$D/reference_panel/coverage_1mb_panel.json
COV_PANEL_CNV=$D/reference_panel/coverage_cnv_panel.json
COV_Z_1MB=5.0
COV_Z_CNV=6.0
BREED_PANEL=$D/breed_panel
BREED_SITES=$BREED_PANEL/sites.tsv     # chrom, pos, a1, a2 in Phat row order
BREED_PHAT=$BREED_PANEL/phat.npy       # (131353 x 230) float32 allele frequencies
BREED_LABELS=$BREED_PANEL/breeds.txt   # 230 canonical breed names, Phat column order
BREED_LASSO=0.3                        # non-negative L1 penalty; see lasso_sweep.py
SCOPE_P=$D/COSMO/analysis/cosmo_scope177Phat.txt     # (143933 SNPs × 177 breeds) full Parker panel allele freq matrix
SCOPE_CLUST=$D/COSMO/analysis/scope_clust.txt        # breed ordering for Phat columns
PARKER_BIM=$D/COSMO/analysis/cosmo_parker_full.bim
PARKER_FAM=$D/COSMO/analysis/cosmo_parker_full.fam
PARKER_BED=$D/COSMO/analysis/cosmo_parker_full.bed

SNPEFF_DB="ROS_Cfam_1.115"   # SnpEff database for canFam4 / ROS_Cfam_1.0

# Microbiome
METAPHLAN_BIN="${METAPHLAN_BIN:-$HOME/Library/Python/3.9/bin/metaphlan}"
MICROBIOME_REF="$D/reference_panel/microbiome_panel.json"

# Strip PATH entries with spaces (e.g. Claude plugin paths) that break $MM word-splitting
PATH=$(echo "$PATH" | tr ':' '\n' | grep -v ' ' | tr '\n' ':' | sed 's/:$//')
export PATH

# Use env bin dirs directly — avoids micromamba lock contention when multiple
# tools run simultaneously in a pipe (bwa | samtools sort | fixmate | markdup).
ENV_GENOMICS="${ENV_GENOMICS:-$HOME/micromamba/envs/genomics}"
ENV_GLIMPSE="${ENV_GLIMPSE:-$HOME/micromamba/envs/glimpse_x86}"
# Put the pinned genomics env first on PATH. Several inline Python blocks invoke
# `samtools`/`bcftools` by bare name; without this they resolve to whatever the
# host happens to provide — Homebrew on the Mac, nothing at all on Hoffman, where
# Stage 6 died with FileNotFoundError. Prepending makes every bare invocation use
# the same pinned build the rest of the pipeline uses via $MM.
export PATH="$ENV_GENOMICS/bin:$PATH"

MM="env PATH=$ENV_GENOMICS/bin:$PATH LD_LIBRARY_PATH=$ENV_GENOMICS/lib"
MM_GLIMPSE="env PATH=$ENV_GLIMPSE/bin:$PATH LD_LIBRARY_PATH=$ENV_GLIMPSE/lib"
NPROC="${NPROC:-8}"
GLIMPSE_PARALLEL="${GLIMPSE_PARALLEL:-8}"

# ── Node-local scratch ───────────────────────────────────────
# On a cluster, $OUT accumulates ~15-20GB of intermediates per sample (3GB BAM,
# merged/trimmed FASTQs, GLIMPSE chunks, SnpEff output). Ninety-six tasks doing
# that on shared project storage is slow and antisocial, so the working
# directory moves to node-local disk and only small artifacts are copied back.
#
# $PUB is NOT redirected — the result JSONs are the product and belong on
# shared storage.
#
# Consequence: intermediates vanish with the job, so --from-stage resume is not
# possible for local-scratch runs. Array jobs always start from stage 1 anyway.
FINAL_OUT=""
# Scratch is for FRESH runs only. A resume (--from-stage > 1) reads intermediates
# that a previous run left in the persistent work directory; redirecting $OUT to
# an empty scratch makes them invisible. That is not hypothetical — job 14285418
# resumed 96 dogs from stage 9, found no BCF, and Stage 9 wrote a uniform
# breed profile from zero coverage before Stage 10 died.
if (( FROM_STAGE > 1 )) && [[ "${USE_LOCAL_SCRATCH:-0}" == "1" ]]; then
  echo "NOTE: resuming from stage $FROM_STAGE — using $OUT directly, not scratch" >&2
  USE_LOCAL_SCRATCH=0
fi
if [[ "${USE_LOCAL_SCRATCH:-0}" == "1" ]]; then
  # Defaults to $TMPDIR (node-local, auto-cleaned by the scheduler). Set
  # LOCAL_SCRATCH_ROOT to a shared scratch filesystem instead if $TMPDIR is too
  # small — that trades local-disk speed for capacity, and we clean up ourselves.
  _scratch="${LOCAL_SCRATCH_ROOT:-${TMPDIR:-}}"
  # A site profile that sets LOCAL_SCRATCH_ROOT is stating where intermediates
  # belong, so create it rather than treating "does not exist yet" as "disabled".
  # Falling back silently is worse than failing: 96 samples x ~18GB of
  # intermediates then accumulate in $OUT forever instead of being cleaned up
  # after each sample, and nothing says so except one line of stderr.
  # $TMPDIR is different — it is a fallback guess, so a missing one just disables.
  if [[ -n "${LOCAL_SCRATCH_ROOT:-}" ]] && ! mkdir -p "$_scratch" 2>/dev/null; then
    echo "ERROR: LOCAL_SCRATCH_ROOT=$_scratch does not exist and cannot be created" >&2
    exit 1
  fi
  if [[ -n "$_scratch" && -d "$_scratch" ]]; then
    FINAL_OUT="$OUT"
    OUT="$_scratch/${DOG_LOWER}/analysis"
    mkdir -p "$FINAL_OUT"
  else
    echo "WARNING: USE_LOCAL_SCRATCH=1 but no usable scratch root" \
         "(LOCAL_SCRATCH_ROOT='${LOCAL_SCRATCH_ROOT:-}', TMPDIR='${TMPDIR:-}')" \
         "— intermediates stay in $OUT and will NOT be cleaned up" >&2
  fi
fi

mkdir -p "$OUT" "$PUB"
LOG=$OUT/pipeline.log

# A logging failure must never decide the pipeline's exit status — see the
# EXIT-trap bug above.
log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG" || true; }
die() { log "ERROR: $*"; exit 1; }

log "Site profile: $DOGS_SITE  (D=$D, NPROC=$NPROC, GLIMPSE_PARALLEL=$GLIMPSE_PARALLEL, publish=$PUBLISH_RESULTS)"
log "Run config: sample=$DOG_NAME from_stage=$FROM_STAGE out=$OUT"

# ── Preflight ────────────────────────────────────────────────
# Every prerequisite is checked up front. Previously a missing tool or reference
# surfaced as a crash deep into a 75-minute run; on a cluster that is 100 array
# tasks failing at 3am instead of one clear message before anything starts.
preflight() {
  local missing=0
  local t
  for t in "$ENV_GENOMICS/bin/bwa-mem2" "$ENV_GENOMICS/bin/samtools" \
           "$ENV_GENOMICS/bin/bcftools" "$ENV_GENOMICS/bin/fastp" \
           "$ENV_GLIMPSE/bin/GLIMPSE2_phase" "$ENV_GLIMPSE/bin/GLIMPSE2_ligate" \
           "$METAPHLAN_BIN"; do
    [[ -x "$t" ]] || { log "  MISSING TOOL: $t"; missing=1; }
  done
  # Every reference the pipeline reads. Derived from a sweep of $D/... paths in
  # this script — stages 8, 11 and 12 each died on a different missing file
  # because the transfer list was assembled by hand.
  for f in "$FASTA" "$DOG10K_PANEL" "$MICROBIOME_REF" "$PARKER_BED" "$PARKER_BIM" "$PARKER_FAM" \
           "$SCOPE_P" "$SCOPE_CLUST" "$OMIA_DB" "$REF_JSON/prs_reference.json" \
           "$REF_JSON/prs_lmm_betas.json.gz" \
           "$REF_JSON/darwins_ark_size_prs.tsv.gz" "$REF_JSON/darwins_ark_blend.json" \
           "$REF_JSON/darwins_ark/manifest.json.gz" "$REF_JSON/darwins_ark/wts_biddability.tsv.gz" \
           "$REF_JSON/mito_haplogroups.json.gz" "$REF_JSON/ortho_risk.json" \
           "$D/reference_panel/coverage_50kb_panel.npz" \
           "$D/relatives_ref/geno.npy" "$D/relatives_ref/meta.json.gz" "$D/relatives_ref/keys.tsv" \
           "$D/COSMO/glimpse2_dog10k/het_out/dog10k_het.het" \
           "$D/COSMO/glimpse2_dog10k/het_out/panel_af.tsv.gz"; do
    [[ -e "$f" ]] || { log "  MISSING REFERENCE: $f"; missing=1; }
  done
  [[ -d "$CHUNKS_DIR" ]] || { log "  MISSING REFERENCE: $CHUNKS_DIR"; missing=1; }
  # GLIMPSE2 is built with AVX2. On a node without it the binary prints its
  # banner and then dies with SIGILL *inside the compute loop* — so it cannot be
  # caught by running --help. Every chunk fails in ~0.2s, Stage 7 produces no
  # genotypes, and the run dies 45 minutes later when bcftools concat gets no
  # input. Job 14270502 task 19 lost DOGS-Gen-28 exactly this way on n6161.
  # A heterogeneous cluster will keep handing us such nodes; catch it in 1ms.
  # Only enforced when this run actually reaches Stage 7 — stage-6-only or
  # stage-9+ backfills never touch GLIMPSE2 and must not be refused a node
  # for a binary they will not execute.
  if (( FROM_STAGE <= 7 && TO_STAGE >= 7 )) \
     && [[ -r /proc/cpuinfo ]] && ! grep -qm1 '\bavx2\b' /proc/cpuinfo; then
    log "  CPU LACKS avx2 on $(hostname) — GLIMPSE2 would crash with SIGILL in Stage 7"
    missing=1
  fi
  # Inline Python calls these by bare name — make sure PATH resolves them into
  # the pinned env rather than to a host copy (or nothing).
  for b in samtools bcftools; do
    _r="$(command -v "$b" 2>/dev/null || true)"
    if [[ -z "$_r" ]]; then
      log "  '$b' not on PATH"; missing=1
    elif [[ "$_r" != "$ENV_GENOMICS/bin/$b" ]]; then
      log "  WARNING: bare '$b' resolves to $_r, not $ENV_GENOMICS/bin/$b"
    fi
  done
  # estimate_chunk indexes every chunk with `$MM_GLIMPSE bcftools index` and
  # DELETES the chunk on failure. The glimpse env has no bcftools of its own, so
  # if it is not inherited from PATH here, a whole Stage 7 silently produces
  # nothing after hours of work. Check the exact invocation the pipeline uses.
  $MM_GLIMPSE bcftools --version >/dev/null 2>&1 \
    || { log "  bcftools not reachable under \$MM_GLIMPSE — every chunk would be deleted after genotyping"; missing=1; }
  $MM_GLIMPSE GLIMPSE2_phase --help >/dev/null 2>&1 \
    || { log "  GLIMPSE2_phase not runnable under \$MM_GLIMPSE"; missing=1; }

  # Stage 10 needs a JVM for SnpEff.
  if [[ ! -x "$ENV_GENOMICS/bin/java" ]] && ! find "$ENV_GENOMICS/lib/jvm/bin" -name java >/dev/null 2>&1 \
     && ! command -v java >/dev/null 2>&1; then
    log "  java not found — Stage 10 (SnpEff) needs it"; missing=1
  fi

  # Stage 15 needs MetaPhlAn's ~8GB database. It is downloaded separately from
  # the tool, so the binary existing proves nothing; without it Stage 15 fails
  # around 40 minutes into a run.
  if [[ -n "${METAPHLAN_DB:-}" ]]; then
    [[ -d "$METAPHLAN_DB" ]] || { log "  METAPHLAN_DB is not a directory: $METAPHLAN_DB"; missing=1; }
    # A pinned index must actually be there, or MetaPhlAn quietly downloads ~40GB
    # of a different version instead of failing.
    if [[ -n "${METAPHLAN_INDEX:-}" && -d "$METAPHLAN_DB" ]]; then
      compgen -G "$METAPHLAN_DB/${METAPHLAN_INDEX}*.bt2l" >/dev/null \
        || { log "  METAPHLAN_INDEX $METAPHLAN_INDEX has no .bt2l files in $METAPHLAN_DB"; missing=1; }
    fi
  elif [[ -x "$METAPHLAN_BIN" ]]; then
    "$METAPHLAN_BIN" --version 2>&1 | grep -qi "installed databases" \
      || log "  WARNING: MetaPhlAn reports no installed database — Stage 15 will fail. Either run
        $METAPHLAN_BIN --install, or set METAPHLAN_DB to an existing database directory."
  fi

  # node only matters where this site actually publishes; the cluster stages
  # results instead and publishes later from a login node.
  if [[ "${PUBLISH_RESULTS:-0}" == "1" ]] && ! command -v node >/dev/null 2>&1; then
    log "  node not found — Stage 17 publishes with it (set PUBLISH_RESULTS=0 to stage instead)"; missing=1
  fi
  # Every inline Python block runs under DATA_PYTHON, so it needs the full set —
  # pysam included, which the stage-8/9 blocks import.
  if [[ -z "${DATA_PYTHON:-}" || ! -x "${DATA_PYTHON:-}" ]]; then
    log "  DATA_PYTHON not set or not executable: ${DATA_PYTHON:-<unset>}"; missing=1
  else
    for _m in numpy pysam pandas scipy sklearn; do
      "$DATA_PYTHON" -c "import $_m" 2>/dev/null || { log "  DATA_PYTHON ($DATA_PYTHON) lacks $_m"; missing=1; }
    done
    "$DATA_PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3,7) else 1)' \
      || { log "  DATA_PYTHON ($DATA_PYTHON) is $("$DATA_PYTHON" -V 2>&1) — 3.7+ required"; missing=1; }
  fi
  (( missing == 0 )) || die "Preflight failed for site '$DOGS_SITE' — see above. Nothing has been run."
  log "Preflight OK"
}
preflight

# Validate a site's setup without running anything:
#     DOGS_SITE=hoffman bash run_dog_pipeline.sh sample_sheet.tsv 2 1 --preflight-only
if [[ " $* " == *" --preflight-only "* ]]; then
  log "--preflight-only — exiting before any work"
  exit 0
fi

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
# Small, high-value artifacts worth keeping from a local-scratch run. The BAM,
# FASTQs, GLIMPSE chunks and SnpEff output are deliberately NOT kept — the JSONs
# in $PUB are the product, and realigning is cheaper than storing 3GB x 96.
#   coverage_*.tsv     feed the panel-of-normals work
#   *_metaphlan.*      let stage 15 be re-run without remapping
#   fastp.*            QC provenance
_keep_from_scratch() {
  [[ -n "$FINAL_OUT" && -d "$OUT" ]] || return 0
  mkdir -p "$FINAL_OUT"
  local f
  # markdup.bam is KEPT while the pipeline is still being iterated on —
  # deleting it forced three realignment campaigns (stage-8/13 reruns,
  # sites.bam, 50kb grid). BUT the per-user quota on the project filesystem
  # is nearly consumed by unrelated data (batch 14583198 went Eqw on 'Disk
  # quota exceeded'), so when BAM_ARCHIVE_DIR is set (site profile) the BAM
  # lands there — the archive filesystem has two orders of magnitude more
  # headroom — with a symlink in FINAL_OUT so stage 8/13 fallback still
  # resolves. Revisit (delete or CRAM) once the callers settle.
  for f in pipeline.log pipeline.done fastp.json fastp.html \
           coverage_1mb.tsv coverage_cnv.tsv coverage_50kb.tsv.gz \
           sites.bam sites.bam.bai sites.bam.csi sites.bed \
           "${DOG_LOWER}_metaphlan.txt" "${DOG_LOWER}_metaphlan.mapout.bz2"; do
    [[ -e "$OUT/$f" ]] && cp -p "$OUT/$f" "$FINAL_OUT/$f" 2>/dev/null || true
  done
  local bam
  for bam in markdup.bam markdup.bam.bai markdup.bam.csi; do
    [[ -e "$OUT/$bam" ]] || continue
    if [[ -n "${BAM_ARCHIVE_DIR:-}" ]] && mkdir -p "$BAM_ARCHIVE_DIR" 2>/dev/null; then
      cp -p "$OUT/$bam" "$BAM_ARCHIVE_DIR/${DOG_LOWER}.$bam" 2>/dev/null \
        && ln -sf "$BAM_ARCHIVE_DIR/${DOG_LOWER}.$bam" "$FINAL_OUT/$bam" 2>/dev/null || true
    else
      cp -p "$OUT/$bam" "$FINAL_OUT/$bam" 2>/dev/null || true
    fi
  done
  # The imputed BCF too. It was deliberately left out as an intermediate, but
  # that made stages 8-13 unrepeatable: re-scoring 96 dogs against a revised
  # breed panel needed a full re-run from stage 1 (~90 min each) instead of
  # minutes, and we have now wanted exactly that twice. ~200MB per dog, so ~19GB
  # for 96 against ~7TB free — cheap insurance for a pipeline whose downstream
  # models are still being iterated on.
  if [[ -d "$OUT/glimpse2" ]]; then
    mkdir -p "$FINAL_OUT/glimpse2"
    for f in "$OUT/glimpse2/"*_imputed_dog10k.bcf*; do
      [[ -e "$f" ]] && cp -p "$f" "$FINAL_OUT/glimpse2/" 2>/dev/null || true
    done
  fi
  # Remove the working directory ourselves, but ONLY on success. $TMPDIR would be
  # cleaned by the scheduler anyway; a shared scratch root would not, and 96 x
  # ~18GB left behind fills 2TB fast. On failure it is left in place so the run
  # can be debugged.
  if [[ -e "$OUT/pipeline.done" ]]; then
    echo "[$(date '+%H:%M:%S')] Kept artifacts in $FINAL_OUT; removing scratch $OUT" \
      | tee -a "$FINAL_OUT/pipeline.log" 2>/dev/null || true
    rm -rf "$OUT" 2>/dev/null || true
    # $LOG lived inside the directory just removed, and _finish still has its
    # summary to print. Without this the tee fails, and because we are inside
    # the EXIT trap under set -e that failure becomes the script's exit status:
    # every task in job 14306032 reported "pipeline exited 1" after completing
    # successfully. Repointing also means the runtime and peak-memory summary
    # lands in the log we keep, instead of being written to a deleted file.
    LOG="$FINAL_OUT/pipeline.log"
    # $OUT is <scratch>/<dog>/analysis — take the <dog> wrapper too, or scratch
    # accumulates one empty directory per sample (96 here, 1045 for the
    # microbiome reference rebuild). rmdir only ever removes an empty directory,
    # so this cannot touch a sibling run's data.
    rmdir "$(dirname "$OUT")" 2>/dev/null || true
  else
    echo "[$(date '+%H:%M:%S')] Kept artifacts in $FINAL_OUT; scratch RETAINED for debugging: $OUT" \
      | tee -a "$FINAL_OUT/pipeline.log" 2>/dev/null || true
  fi
}

_finish() {
  kill "$MEM_POLL_PID" 2>/dev/null || true
  _keep_from_scratch
  local end=$(date +%s)
  local elapsed=$(( end - PIPELINE_START ))
  local h=$(( elapsed / 3600 ))
  local m=$(( (elapsed % 3600) / 60 ))
  local s=$(( elapsed % 60 ))
  # An existing-but-empty peak file made $peak_kb empty, and the awk below then
  # parsed "/1048576" as a regex — failing the exit trap after a successful run.
  local peak_kb=$(cat "$PEAK_MEM_FILE" 2>/dev/null || echo 0)
  peak_kb=${peak_kb:-0}
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

if (( FROM_STAGE <= 1 && TO_STAGE >= 1 )); then

# ── Stage 1: Merge FASTQ ─────────────────────────────────────
log "=== Stage 1: Merge FASTQ lanes ==="
R1_FILES=$(ls "$FASTQ_DIR"/*_R1_*.fastq.gz 2>/dev/null | sort -V) || die "No R1 FASTQ files in $FASTQ_DIR"
R2_FILES=$(ls "$FASTQ_DIR"/*_R2_*.fastq.gz 2>/dev/null | sort -V) || die "No R2 FASTQ files in $FASTQ_DIR"
log "R1 files: $(echo $R1_FILES | tr ' ' '\n' | wc -l)"
log "R2 files: $(echo $R2_FILES | tr ' ' '\n' | wc -l)"

cat $R1_FILES > "$OUT/merged_R1.fastq.gz"
cat $R2_FILES > "$OUT/merged_R2.fastq.gz"
log "Merged: R1=$(ls -lh $OUT/merged_R1.fastq.gz | awk '{print $5}'), R2=$(ls -lh $OUT/merged_R2.fastq.gz | awk '{print $5}')"
fi # end stage 1

if (( FROM_STAGE <= 2 && TO_STAGE >= 2 )); then
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
fi # end stage 2

if (( FROM_STAGE <= 4 && TO_STAGE >= 4 )); then
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

# ── Permanent mini-BAM of known-variant regions ──────────────────────────
# Disk cleanup deletes markdup.bam eventually; when stage 8/13 was later rerun
# on 92 BAM-less dogs it silently degraded every read-based call to
# "insufficient reads". Keep a few-MB extract (catalogue sites, MC1R, TYRP1,
# PMEL merle junction, chrM) so read-based stages can always be rerun.
"$DATA_PYTHON" - << PYEOF > "$OUT/sites.bed"
import json
regs = [('chrM', 0, 17000), ('chr10', 643510, 645512),
        ('chr5', 64185000, 64189000), ('chr11', 33399000, 33402000)]
db = json.load(open("$OMIA_DB"))
for v in db['variants']:
    c, p = v.get('chrom'), v.get('pos')
    if not c or not p: continue
    e = v.get('pos_end') or p
    regs.append((c, max(0, int(p) - 500), int(e) + 500))
for c, s, e in regs: print(f'{c}\t{s}\t{e}')
PYEOF
$MM samtools view -b -M -L "$OUT/sites.bed" "$OUT/markdup.bam" > "$OUT/sites.bam"
$MM samtools index "$OUT/sites.bam"
log "sites.bam ready ($($MM samtools view -c "$OUT/sites.bam") reads at known-variant regions)"
fi # end stages 3+4

if (( FROM_STAGE <= 5 && TO_STAGE >= 5 )); then
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

# 5c — fixed 50kb grid for the unified CNV framework (planned replacement for
# the dual 1Mb/adaptive scheme): one window scale for every dog regardless of
# depth, so a cohort panel-of-normals can be built per window and segments can
# span 50kb to whole chromosomes. Emitted for every sample NOW so the batch in
# flight captures it while its BAM still exists (~2MB gz; kept from scratch).
awk 'BEGIN{OFS="\t"} $1~/^chr([0-9]+|X)$/ {
  for(s=0; s<$2; s+=50000)
    print $1, s, (s+50000<$2 ? s+50000 : $2)
}' "$FASTA.fai" > "$OUT/windows_50kb.bed"
$MM samtools bedcov "$OUT/windows_50kb.bed" "$OUT/markdup.bam" | gzip > "$OUT/coverage_50kb.tsv.gz"
log "50kb coverage: $(gunzip -c "$OUT/coverage_50kb.tsv.gz" | wc -l) windows"
fi # end stage 5

if (( FROM_STAGE <= 6 && TO_STAGE >= 6 )); then
# ── Stage 6: Coverage + QC + CNV JSON ────────────────────────
log "=== Stage 6: Coverage / QC / CNV JSON ==="

# CNV_WINDOW is computed in stage 5, so resuming at stage 6 left it unset and
# set -u killed the run. That mattered: stage 5 needs markdup.bam, which
# local-scratch runs discard, whereas coverage_cnv.tsv is kept — so stage 6 is
# exactly where a coverage re-analysis has to restart. Recover the window width
# from the data instead of depending on an earlier stage's variable.
if [[ -z "${CNV_WINDOW:-}" ]]; then
  [[ -s "$OUT/coverage_cnv.tsv" ]] \
    || die "stage 6 needs $OUT/coverage_cnv.tsv (run from stage 5 to regenerate it)"
  CNV_WINDOW=$(awk 'NR==1 {print $3-$2; exit}' "$OUT/coverage_cnv.tsv")
  log "Recovered CNV window from coverage_cnv.tsv: ${CNV_WINDOW} bp"
fi
"$DATA_PYTHON" - << PYEOF
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

COMMON_PCT = 5.0    # cohort share at or above which a region is "common"

def annotate_freq(events, idx, bin_bp, sex):
    """Tag each event with how often the reference cohort is flagged there.

    Recurrence analysis over 93 dogs found chr19:21.3-21.5Mb flagged as a gain
    in 42 of them and chr5:85.9Mb in 41, while NO bin was flagged in more than
    half. That pattern is the point: a reference artifact would hit nearly every
    dog, since all dogs map to the same reference, whereas hitting a fraction is
    what a real copy-number polymorphism looks like. The events are genuine and
    still not findings -- a region carried by 45% of dogs is common, and
    presenting it as a discovery about one dog is a category error.

    Takes the MAXIMUM frequency across the event's bins, not the peak bin's.
    Conservative on purpose: if any part of an event coincides with a region
    dogs commonly vary at, it is most likely that polymorphism, and the cost of
    wrongly calling something rare (alarming an owner) exceeds the cost of
    wrongly calling it common (a quieter listing).
    """
    for ev in events:
        key = 'all' if ev['chrom'] != 'chrX' else ('F' if sex == 'female' else 'M')
        fr = []
        for b in range(ev['start'], ev['end'], bin_bp):
            st = idx.get((ev['chrom'], b, key))
            if st and st.get('n'):
                c = st.get('ng', 0) if ev['direction'] == 'gain' else st.get('nl', 0)
                fr.append(c / st['n'])
        ev['cohort_freq'] = round(max(fr) * 100, 1) if fr else None
        ev['novelty'] = ('unknown' if ev['cohort_freq'] is None
                         else 'common' if ev['cohort_freq'] >= COMMON_PCT else 'rare')
    return events

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

# 96-dog panel-of-normals: per-window median and robust SD, so a window can be
# reported as "N SD from 96 dogs" rather than just "looks low".
cov_panel_path = "$COV_PANEL_1MB"
Z_CUT = float("$COV_Z_1MB")
try:
    with open(cov_panel_path) as _f:
        _pd = _json.load(_f)
    panel_idx = {(c, e['start'], k): e[k]
                 for c, es in _pd['panel'].items() for e in es
                 for k in ('all', 'F', 'M') if k in e}
    panel_n = _pd.get('meta', {}).get('n_samples')
    # Which panel produced this report. The panels are no longer tracked in git
    # (they are regenerated build artifacts, and the CNV one is 29MB), so the
    # identifying information lives here instead — enough to tell whether two
    # reports were scored against the same reference.
    import hashlib as _hl
    with open(cov_panel_path, 'rb') as _pf:
        panel_sha = _hl.sha256(_pf.read()).hexdigest()[:12]
except Exception as _e:
    panel_idx, panel_n = {}, None
    print(f"WARNING: no coverage panel at {cov_panel_path} ({_e}) — "
          "windows will carry no significance")

result = {}
for chrom, arr in raw.items():
    # Tooltip depth tracks. These used to be copied from a legacy reference
    # file that lost its 'cosmo'/'panel' keys in a rebuild, which silently
    # zero-padded both series in every report since. The sample's own depth
    # is already in hand, and the panel's expected depth is derived per
    # window from the panel stats in the loop below.
    cosmo_arr = arr
    panel_arr = []   # per window: panel median ratio × sample median depth, or None
    n = len(arr)
    if chrom == 'chrX' and predicted_sex == 'male':
        # PAR1 windows are diploid in males (depth ≈ autosomal); non-PAR windows are hemizygous (depth ≈ 0.5×).
        # Use per-window normalization: above 0.75× autosomal → PAR1 → divide by kiki_median;
        # below 0.75× autosomal → non-PAR → divide by kiki_median × 0.5, so both map to ratio 1.0.
        par_thresh = kiki_median * 0.75
        ratio_arr, scale_arr = [], []
        for d in arr:
            if d == 0:
                ratio_arr.append(0.0); scale_arr.append(1.0)
            elif d >= par_thresh:
                ratio_arr.append(round(d / kiki_median, 4))        # PAR1: diploid norm
                scale_arr.append(1.0)
            else:
                ratio_arr.append(round(d / (kiki_median * 0.5), 4)) # non-PAR: hemizygous norm
                scale_arr.append(2.0)
    else:
        ratio_arr = [round(d / kiki_median, 4) if d > 0 else 0.0 for d in arr]
        scale_arr = [1.0] * len(arr)
    # Robust z against the panel. Computed from the UNADJUSTED ratio
    # (depth / autosomal median), because that is the space the panel medians
    # live in — ratio_arr rescales male chrX so hemizygous reads as 1.0, which
    # would not line up with a male chrX panel median of ~0.5.
    pkey = 'all' if chrom != 'chrX' else ('F' if predicted_sex == 'female' else 'M')
    # Alongside this dog's z, carry how often the REFERENCE COHORT itself varies
    # at each window. A reader looking at a dip or a bump cannot tell whether it
    # matters; "34 of 93 healthy dogs are also low here" answers that directly,
    # and answers it for every window rather than only the flagged ones.
    #
    # Deliberately computed even where this dog's own window is unscoreable:
    # the frequency describes the cohort, not the sample, so "do other dogs vary
    # here?" has an answer either way.
    # The band is the more useful of the two for reading the chart. At 1Mb the
    # flag frequency is 0 almost everywhere (megabase averaging washes out focal
    # CNVs — only 8 of 2,498 windows are flagged in >=5% of the cohort), so on
    # its own it cannot tell a reader whether the bump they are looking at is
    # ordinary. The spread of the reference dogs at that window can: drawn
    # behind the bars, a bar inside the band is normal variation and a bar
    # outside it is not, with no numbers to interpret.
    #
    # Emitted in the SAME units as ratio_arr, hence scale_arr: on a male chrX
    # the ratio is rescaled so hemizygous reads as 1.0, while the panel median
    # sits at ~0.5, and a band in raw panel units would sit at half height and
    # look alarming on every male.
    BAND_SD = 2.0
    z_arr, pct_arr, lo_arr, hi_arr = [], [], [], []
    for i, d in enumerate(arr):
        st = panel_idx.get((chrom, i * 1000000, pkey))
        sc = scale_arr[i] if i < len(scale_arr) else 1.0
        # Expected depth of a typical panel dog here, expressed at this
        # sample's sequencing depth so Sample and Panel are comparable.
        panel_arr.append(round(st['median'] * kiki_median, 2) if st else None)
        # THE question a reader has: of the reference dogs, how many show this
        # much coverage here, or more? Counted from the panel's actual values at
        # this window, in the direction this dog deviates. A threshold-crossing
        # frequency cannot answer it — that asks how many dogs cross a fixed
        # line, not how many are as extreme as this one — and median/spread
        # cannot either, because the answer lives in the tail.
        vals = st.get('vals') if st else None
        if vals and d > 0:
            obs = d / kiki_median
            if obs >= st['median']:
                cnt = sum(1 for v in vals if v >= obs)
            else:
                cnt = sum(1 for v in vals if v <= obs)
            pct_arr.append(round(100 * cnt / len(vals), 1))
        else:
            pct_arr.append(None)
        if st and 0.004 <= st['mad_sd'] <= 0.30:
            lo_arr.append(round(max(0.0, st['median'] - BAND_SD * st['mad_sd']) * sc, 4))
            hi_arr.append(round((st['median'] + BAND_SD * st['mad_sd']) * sc, 4))
        else:
            lo_arr.append(None)
            hi_arr.append(None)
        if not st or d <= 0 or not (0.004 <= st['mad_sd'] <= 0.30):
            z_arr.append(None)
            continue
        z_arr.append(round(((d / kiki_median) - st['median']) / st['mad_sd'], 2))
    result[chrom] = {
        'cosmo': [round(v, 2) for v in cosmo_arr],
        'panel': panel_arr,
        'ratio': ratio_arr,
        'z': z_arr,
        'pct_dogs': pct_arr,
        'band_lo': lo_arr,
        'band_hi': hi_arr,
        'band_sd': BAND_SD,
    }

# Merge runs of adjacent flagged windows into events. "3 regions of unusual
# coverage" is a statement a dog owner can act on; "7 flagged windows" is not.
# Regions are the runs of windows the CHART draws in red or orange — i.e. the
# ones a reader can actually see — not runs defined by a z statistic.
#
# z was the wrong instrument here. It assumes the reference dogs are roughly
# normally distributed at a window, and at the scale of a real deletion they are
# not: at chr8's end the 93 panel dogs span 0.43-1.47, so a dog near the bottom
# scores z -1.1 and is "not flagged" while being in the lowest few percent of
# the population. That produced a chart shouting DELETION in red at windows the
# table said nothing about, and a tooltip quoting 2-4% beside a table row
# claiming "not seen".
#
# The honest statistic is the empirical one, and it is exactly a minor allele
# frequency: of the reference dogs, what share carry this much loss or gain
# here. No distributional assumption, and it means the same thing whether the
# deletion is 1Mb or 40.
DEL_T, DUP_T = 0.65, 1.35
events = []
for chrom in sorted(result):
    ratios = result[chrom].get('ratio') or []
    pcts   = result[chrom].get('pct_dogs') or []
    run = None
    for i, r in enumerate(list(ratios) + [None]):
        # ratio 0 means no coverage measured, not a deletion
        d = None
        if r is not None and r > 0:
            if r < DEL_T: d = 'loss'
            elif r > DUP_T: d = 'gain'
        if d and run is None:
            run = [i, i, d]
        elif d and run and run[2] == d:
            run[1] = i
        elif run is not None:
            lo, hi, direction = run
            seg = [pcts[j] for j in range(lo, hi + 1) if j < len(pcts) and pcts[j] is not None]
            events.append({'chrom': chrom, 'start': lo * 1000000,
                           'end': (hi + 1) * 1000000,
                           'windows': hi - lo + 1,
                           'direction': direction,
                           # range across the region, because that is what a
                           # reader sees when they hover along it
                           'pct_min': min(seg) if seg else None,
                           'pct_max': max(seg) if seg else None,
                           'pct_dogs': min(seg) if seg else None})
            run = [i, i, d] if d else None
events.sort(key=lambda e: (e['pct_dogs'] if e['pct_dogs'] is not None else 999, -e['windows']))

# Whole-chromosome events are aneuploidy, not copy-number variation, and must
# not be reported as "a gain on chrX". DOGS-Gen-111 is the case that forced
# this: 124 of her 125 chrX windows were elevated (X:autosome 1.461, i.e. three
# copies), which the event merger split into two large gains purely because one
# unscoreable window interrupted the run. Her breed-matched, same-batch control
# DOGS-Gen-110 measured 1.001, so the signal is hers and not the assay's.
#
# Detect per chromosome rather than by widening the merge gap: a chromosome
# qualifies when most of its SCOREABLE windows are flagged the same way. Using a
# fraction rather than a count makes it work on chr38 as well as chr1.
ANEUPLOIDY_FRAC = 0.80
aneuploidy = []
aneuploid_chroms = set()
for chrom in sorted(result):
    zs = result[chrom].get('z') or []
    scoreable = [i for i, z in enumerate(zs) if z is not None]
    if len(scoreable) < 10:
        continue                      # too little to judge a whole chromosome
    up   = [i for i in scoreable if zs[i] >=  Z_CUT]
    down = [i for i in scoreable if zs[i] <= -Z_CUT]
    hit, direction = (up, 'gain') if len(up) >= len(down) else (down, 'loss')
    frac = len(hit) / len(scoreable)
    if frac < ANEUPLOIDY_FRAC:
        continue
    # Copy number straight from the unadjusted ratio: 2 copies is the diploid
    # baseline, so a male's single X (~0.5) and a trisomy (~1.5) both fall out
    # of the same arithmetic without special-casing.
    rr = [raw[chrom][i] / kiki_median for i in scoreable if raw[chrom][i] > 0]
    med_ratio = _stats.median(rr) if rr else 0.0
    copies = round(med_ratio * 2)
    expected = 1 if (chrom == 'chrX' and predicted_sex == 'male') else 2
    if copies == expected:
        continue                      # flagged but not a copy-number change
    aneuploidy.append({
        'chrom': chrom, 'direction': direction,
        'call': {0: 'nullisomy', 1: 'monosomy',
                 3: 'trisomy', 4: 'tetrasomy'}.get(copies, f'{copies}-copy'),
        'copies_est': copies, 'expected_copies': expected,
        'median_ratio': round(med_ratio, 3),
        'fraction_flagged': round(frac, 3),
        'windows_flagged': len(hit), 'windows_scoreable': len(scoreable),
        'peak_z': max((zs[i] for i in hit), key=abs)})
    aneuploid_chroms.add(chrom)

# An aneuploid chromosome is fully described by its aneuploidy entry, so its
# windows must not ALSO appear in the event list — otherwise a trisomy reads as
# two large gains, which is both wrong in kind and needlessly alarming.
if aneuploid_chroms:
    events = [e for e in events if e['chrom'] not in aneuploid_chroms]
for a in aneuploidy:
    print(f"ANEUPLOIDY: {a['chrom']} {a['call']} — {a['copies_est']} copies "
          f"(expected {a['expected_copies']}), median ratio {a['median_ratio']}, "
          f"{a['windows_flagged']}/{a['windows_scoreable']} windows flagged")

annotate_freq(events, panel_idx, 1000000, predicted_sex)

# Per-dog confidence. A dog with 242 events should not produce 242 findings —
# it should say the coverage is not clean enough to call structural variation.
#
# Thresholds read off the 93-dog cohort, which separates cleanly: 93 dogs at
# 0-9 events, three at 55-270. A six-fold gap, so 20 is comfortably inside it.
# Depth is NOT the criterion — the lowest-coverage dog in the cohort (1.00x)
# produced 9 events, while a 2.40x dog produced 89 of which 85 were gains. That
# asymmetry is the second signal: real copy-number variation is not
# one-directional, so a gain-dominated profile is a mapping artifact even when
# the count alone would pass.
_gains = sum(1 for e in events if e['direction'] == 'gain')
_losses = len(events) - _gains
_reasons = []
# The old thresholds were calibrated against z-defined events and do not carry
# over. In particular the gain-dominated rule now misfires on clean samples:
# with regions defined by the ratio bands, gains above 1.35 are common at
# chromosome ends while losses below 0.65 are rare, so gain-dominance is the
# NORMAL pattern — Kiki2, platform-matched and clean, comes out 11 gains to 1
# loss. That rule is dropped rather than retuned; it was an artifact signature
# for a statistic no longer in use.
#
# What remains is a volume check. A clean dog sits around a dozen regions
# (Kiki2 12, Luna Genotek 13), so this only catches samples an order of
# magnitude worse. PROVISIONAL: wants recalibrating across the 96-dog cohort
# on this new basis before it is trusted as a QC gate.
if len(events) > 40:
    _reasons.append(f'{len(events)} flagged regions (a clean sample has roughly a dozen)')
cov_confidence = 'low' if _reasons else 'ok'
# Carry the verdict onto each aneuploidy entry. "Your dog has three copies of
# the X chromosome" is a strong claim to render, and whoever renders it should
# not have to remember to cross-check a sibling field before doing so.
for a in aneuploidy:
    a['confidence'] = cov_confidence
if _reasons:
    print(f"coverage confidence LOW: {'; '.join(_reasons)}")
print(f"coverage events at |z| >= {Z_CUT}: {len(events)} "
      f"({_gains} gain / {_losses} loss), confidence={cov_confidence}")

result['_meta'] = {'predicted_sex': predicted_sex,
                   'panel_sha': panel_sha,
                   'chrx_auto_ratio': round(x_auto_ratio, 3),
                   'panel_n': panel_n, 'z_threshold': Z_CUT,
                   'coverage_confidence': cov_confidence,
                   'confidence_reasons': _reasons,
                   'n_gain': _gains, 'n_loss': _losses,
                   'n_rare': sum(1 for e in events if e.get('novelty') == 'rare'),
                   'common_pct': COMMON_PCT,
                   'events': events,
                   'aneuploidy': aneuploidy}
with open(f'{pub}/coverage_1mb.json', 'w') as f:
    json.dump(result, f)
print(f"coverage_1mb.json: {len(result)-1} chromosomes, median depth {kiki_median:.2f}×, sex={predicted_sex}")

# --- qc_result.json (derived from 1Mb windows) ---
depths = [w[3] for w in windows_1mb]
mean_d   = statistics.mean(depths)
median_d = statistics.median(depths)
std_d    = statistics.stdev(depths)
pct = {t: round(sum(1 for d in depths if d >= t)/len(depths)*100, 1) for t in [10,15,20,30]}
# Thresholds for what this product actually does, not for base-resolution
# variant calling. 20x/15x are the bar when you are genotyping individual
# positions from the reads; nothing here does that. Genotypes come from
# GLIMPSE2 imputation against Dog10K, which is designed for low-pass input,
# and copy number is read at 1Mb, where a 7x sample puts ~70,000 reads in
# every window — a coefficient of variation under half a percent, against the
# 50% shift a real deletion produces. Calling that "FAIL" told the report, and
# the assistant reading it, that a perfectly good sample was inadequate.
# GLIMPSE2 is designed for exactly this input: published imputation accuracy
# holds well down to ~1x, so anything at or above that is a normal sample for
# this product and gets no flag at all. Below 1x accuracy starts to degrade
# (WARN); below 0.2x there is essentially no library and nothing downstream
# can be trusted (FAIL).
qc_status = 'PASS' if mean_d >= 1 else ('WARN' if mean_d >= 0.2 else 'FAIL')

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

# Reads per 1Mb window, and the Poisson precision that follows from it. This is
# the number that decides whether megabase copy number is trustworthy.
_rl = read_stats.get('read_length_bp') or 100
_reads_per_mb = (mean_d * 1_000_000) / _rl
_mb_cv_pct = round(100 / (_reads_per_mb ** 0.5), 2) if _reads_per_mb > 0 else None

qc = {'genome_mean_depth': round(mean_d,1), 'genome_median_depth': round(median_d,1),
    'genome_std_depth': round(std_d,1), 'uniformity_cv': round(std_d/mean_d,3),
    'pct_bins_gt10x': pct[10], 'pct_bins_gt15x': pct[15],
    'pct_bins_gt20x': pct[20], 'pct_bins_gt30x': pct[30],
    'n_low_bins': sum(1 for d in depths if d < 15), 'n_total_bins': len(depths),
    'chromosomes': chroms_out, 'qc_status': qc_status,
    'warning': None if mean_d >= 1 else
               f"Mean depth {mean_d:.1f}x is below 1x; genotype imputation accuracy degrades below this point."
               if mean_d >= 0.2 else
               f"Mean depth {mean_d:.1f}x — essentially no usable coverage; results cannot be trusted.",
    'assessment': (f"Mean genome coverage {mean_d:.1f}x across {len(depths)} 1Mb bins. "
                   f"At 1Mb resolution that is ~{_reads_per_mb:,.0f} reads per window "
                   f"(±{_mb_cv_pct}% counting precision), which is ample for copy-number "
                   f"analysis; a deletion shifts a window by ~50%."),
    'reads_per_mb_window': round(_reads_per_mb),
    'mb_precision_pct': _mb_cv_pct,
    'depth_note': ('Low-pass whole-genome sequencing. Genotypes are imputed against the '
                   'Dog10K panel rather than called from reads, so per-base depth '
                   'thresholds used for variant calling (20x, 30x) do not apply.'),
    'method': 'samtools bedcov over 1Mb bins',
    **read_stats}
with open(f'{pub}/qc_result.json', 'w') as f:
    json.dump(qc, f, indent=2)
print(f"qc_result.json: {qc_status}, mean={mean_d:.1f}x, cnv_window={cnv_win}bp")

# --- cnv_homdel.json (adaptive-window data → structured for CnvTable) ---
windows_cnv = load_tsv(tsv_cnv)
mean_d2 = statistics.mean(d for _,_,_,d in windows_cnv)
hom_del_thresh = mean_d2 * 0.15

# Load gene annotations — read from cosmo reference dir (always present),
# not from the sample pub dir which may not have it yet at Stage 6.
import json as _json, glob as _glob, os as _os
_gene_paths = [f'$REF_JSON/cnv_genes.json', f'{pub}/cnv_genes.json']
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

# ── CNV-scale panel of normals ──────────────────────────────────────────────
# The legacy reference above is 6 dogs pooled at 1Mb, used to judge regions as
# small as 50kb — 20x coarser than the calls it was vetting. This panel is 92
# dogs on a fixed 50kb grid with a per-bin median and robust SD, so a 50kb call
# is finally compared against many dogs at its own scale.
#
# Stage 5 sizes this dog's windows as 50000/mean_depth, so they never line up
# with the panel's grid. Re-bin by overlap, identically to build_panels.py —
# valid only downward, which is why a native window wider than the bin is
# refused rather than stretched.
cnv_panel_path = "$COV_PANEL_CNV"
Z_CUT_CNV = float("$COV_Z_CNV")
cnv_panel_idx, cnv_panel_n, cnv_bin = {}, None, None
try:
    with open(cnv_panel_path) as _f:
        _cp = _json.load(_f)
    cnv_bin = _cp['meta']['window_bp']
    cnv_panel_n = _cp['meta'].get('n_samples')
    cnv_panel_idx = {(c, e['start'], k): e[k]
                     for c, es in _cp['panel'].items() for e in es
                     for k in ('all', 'F', 'M') if k in e}
    print(f"CNV panel: {cnv_panel_n} dogs at {cnv_bin}bp, {len(cnv_panel_idx)} bins")
except Exception as _e:
    print(f"WARNING: no CNV panel at {cnv_panel_path} ({_e}) — "
          "ref_depth_pct falls back to the legacy 1Mb reference, no CNV z-scores")

def rebin_cnv(rows, bin_bp):
    acc, span = {}, {}
    for c, s, e, dpb in rows:
        for b in range(s // bin_bp, (e - 1) // bin_bp + 1):
            lo, hi = max(s, b * bin_bp), min(e, (b + 1) * bin_bp)
            if hi > lo:
                k = (c, b * bin_bp)
                acc[k] = acc.get(k, 0.0) + dpb * (hi - lo)
                span[k] = span.get(k, 0.0) + (hi - lo)
    # divide by length ACTUALLY covered, not bin_bp: a chromosome's final bin is
    # partial, and using the full width would make every chromosome end look
    # depleted — the bug that once flagged 39 windows in 96 of 96 dogs.
    return {k: v / span[k] for k, v in acc.items() if span[k] > 0}

cnv_z, cnv_events = {}, []
if cnv_panel_idx:
    _widest = max(e - s for _, s, e, _ in windows_cnv)
    if _widest > cnv_bin:
        print(f"WARNING: native CNV window {_widest}bp exceeds panel bin {cnv_bin}bp — "
              "skipping CNV z-scores (re-binning upward would invent resolution)")
    else:
        _binned = rebin_cnv(windows_cnv, cnv_bin)
        _auto = [v for (c, _), v in _binned.items() if c != 'chrX' and v > 0]
        _base = statistics.median(_auto) if _auto else 0.0
        if _base <= 0:
            print("WARNING: zero autosomal CNV coverage — skipping CNV z-scores")
        else:
            for (c, s), v in _binned.items():
                key = 'all' if c != 'chrX' else ('F' if predicted_sex == 'female' else 'M')
                st = cnv_panel_idx.get((c, s, key))
                if not st or v <= 0 or not (0.004 <= st['mad_sd'] <= 0.30):
                    continue
                cnv_z[(c, s)] = ((v / _base) - st['median']) / st['mad_sd']
            # Merge contiguous flagged bins into events. Unlike the homozygous
            # deletion path below, this sees GAINS too — the old detector could
            # only ever find losses, because its only test was norm < 0.15.
            for c, s in sorted(cnv_z):
                z = cnv_z[(c, s)]
                if abs(z) < Z_CUT_CNV:
                    continue
                if cnv_events and cnv_events[-1]['chrom'] == c and s == cnv_events[-1]['end']:
                    ev = cnv_events[-1]
                    ev['end'] = s + cnv_bin
                    ev['bins'] += 1
                    if abs(z) > abs(ev['peak_z']):
                        ev['peak_z'] = z
                else:
                    cnv_events.append({'chrom': c, 'start': s, 'end': s + cnv_bin,
                                       'bins': 1, 'peak_z': z})
            for ev in cnv_events:
                ev['peak_z'] = round(ev['peak_z'], 2)
                ev['direction'] = 'gain' if ev['peak_z'] > 0 else 'loss'
                ev['size_kb'] = round((ev['end'] - ev['start']) / 1000, 1)
            # A trisomic chromosome is elevated in every bin, so this detector
            # sees it as dozens of independent gains — 37 of them on a synthetic
            # trisomy X. The 1Mb path already suppresses those; without the same
            # filter here the aneuploidy gets reported twice, the second time in
            # entirely the wrong vocabulary.
            _sup = [e for e in cnv_events if e['chrom'] in aneuploid_chroms]
            if _sup:
                cnv_events = [e for e in cnv_events if e['chrom'] not in aneuploid_chroms]
                print(f"CNV: suppressed {len(_sup)} event(s) on aneuploid "
                      f"{', '.join(sorted(aneuploid_chroms))} — reported as aneuploidy instead")
            cnv_events.sort(key=lambda e: -abs(e['peak_z']))
            annotate_freq(cnv_events, cnv_panel_idx, cnv_bin, predicted_sex)
            _cg = sum(1 for e in cnv_events if e['direction'] == 'gain')
            _cr = sum(1 for e in cnv_events if e.get('novelty') == 'rare')
            print(f"CNV events at |z| >= {Z_CUT_CNV}: {len(cnv_events)} "
                  f"({_cg} gain / {len(cnv_events)-_cg} loss), {_cr} rare "
                  f"(<{COMMON_PCT}% of cohort), over {len(cnv_z)} scoreable bins")

# All genes across all chroms
all_genes = [g for gs in gene_map.values() for g in gs]

# Deletions are called on the PANEL's own 50kb grid, not on this dog's adaptive
# window. Stage 5 sizes that window as max(15000, 50000/depth), so a 7x dog is
# called at 15kb — but the reference it is judged against is binned at 50kb, and
# a 15kb call checked against a bin averaging 3.3x more sequence is checked
# against diluted evidence. The 15kb floor was never statistically motivated
# either: at 7x a 15kb window holds ~1,050 reads, putting the <15% threshold 28
# standard deviations below the mean, so resolution was never the limit.
#
# Calling at the panel's resolution means every call can be answered with the
# question that actually matters — how many reference dogs are this low here —
# instead of the one we were guessing at, which was whether low coverage in the
# panel meant bad mappability or a deletion other dogs also carry. Those two are
# indistinguishable in a median, and we were labelling both "mappability
# artefact".
raw_dels = []
if cnv_panel_idx and cnv_bin and windows_cnv:
    _widest_d = max(e - s for _, s, e, _ in windows_cnv)
    if _widest_d <= cnv_bin:
        _b = rebin_cnv(windows_cnv, cnv_bin)
        _auto_d = [v for (c, _), v in _b.items() if c != 'chrX' and v > 0]
        _base_d = statistics.median(_auto_d) if _auto_d else 0.0
        if _base_d > 0:
            for (c, st_bp), v in sorted(_b.items()):
                ratio_d = v / _base_d
                if ratio_d >= 0.15:
                    continue
                key = 'all' if c != 'chrX' else ('F' if predicted_sex == 'female' else 'M')
                pstat = cnv_panel_idx.get((c, st_bp, key))
                vals = pstat.get('vals') if pstat else None
                pct = (round(100 * sum(1 for x in vals if x <= ratio_d) / len(vals), 1)
                       if vals else None)
                raw_dels.append({'chrom': c, 'start': st_bp, 'end': st_bp + cnv_bin,
                                 'depth': round(v, 2), 'norm': round(ratio_d, 3),
                                 'window_bp': cnv_bin, 'pct_dogs': pct})
    else:
        print(f"WARNING: native CNV window {_widest_d}bp exceeds panel bin {cnv_bin}bp — "
              "no deletion calls")
else:
    print("WARNING: no CNV panel — no deletion calls")

# Merge adjacent windows (gap ≤ 1 window) into contiguous regions
raw_dels.sort(key=lambda w: (w['chrom'], w['start']))
merged = []
for w in raw_dels:
    if merged and merged[-1]['chrom'] == w['chrom'] and w['start'] <= merged[-1]['end'] + w['window_bp']:
        r = merged[-1]; r['end'] = max(r['end'], w['end'])
        r['norms'].append(w['norm']); r['pcts'].append(w['pct_dogs'])
    else:
        merged.append({'chrom': w['chrom'], 'start': w['start'], 'end': w['end'],
                       'norms': [w['norm']], 'pcts': [w['pct_dogs']]})

regions = []; disrupted_all = {}; _nongenic = 0
for r in merged:
    size_bp = r['end'] - r['start']
    avg_norm = sum(r['norms']) / len(r['norms'])
    # How many reference dogs are at least this low here. Same statistic as the
    # karyotype, and it replaces both ref_depth_pct and the artefact verdict.
    #
    # That verdict could not be supported: it called a region a mappability
    # artefact whenever the reference dogs were also low, which is equally what
    # a deletion many dogs carry looks like. On Luna Genotek seven of eight so
    # labelled had a tight panel spread (consistent with mappability) but
    # chr27:20.400Mb had a robust SD of 0.281 against a median of 0.356 — the
    # dogs disagree with each other there, which is a polymorphism, not a
    # mapping failure. Reporting the frequency states what we observe and
    # asserts no mechanism: a uniformly unmappable region comes out near 100%
    # and is self-evidently uninteresting.
    seg_p = [x for x in r['pcts'] if x is not None]
    pct_dogs = max(seg_p) if seg_p else None
    pct_min  = min(seg_p) if seg_p else None
    # Average coverage the reference dogs have here, as a percent of normal.
    # Distinct from the frequency and worth showing beside it: the frequency
    # says how many dogs are as low as this one, while this says what the
    # region looks like in a typical dog. A region where the panel sits at 30%
    # is hard to map for everyone, and a rare call there deserves more scepticism
    # than the same call where the panel sits at 100%.
    _k = 'all' if r['chrom'] != 'chrX' else ('F' if predicted_sex == 'female' else 'M')
    _meds = [cnv_panel_idx[(r['chrom'], b, _k)]['median']
             for b in range(r['start'], r['end'], cnv_bin)
             if (r['chrom'], b, _k) in cnv_panel_idx]
    ref_depth_pct = round(statistics.median(_meds) * 100) if _meds else None
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
    # Report only CNVs that disrupt a gene. The size bar let through losses in
    # gene deserts — Luna Genotek's sole reported CNV was 120kb at the tip of
    # chr8 hitting nothing at all. True, and not actionable: an owner cannot do
    # anything with it and a clinician would not either. What makes a copy
    # number change worth a customer's attention is the gene it removes.
    if not disrupted:
        _nongenic += 1
        continue
    size_str = (f"{size_bp/1e6:.2f}Mb" if size_bp >= 1_000_000
                else f"{size_bp//1000}kb" if size_bp >= 1000 else f"{size_bp}bp")
    regions.append({'chrom': r['chrom'], 'start': r['start'], 'end': r['end'], 'size': size_str,
                    'sample_pct_mean': round(avg_norm*100),
                    'pct_dogs': pct_dogs, 'pct_min': pct_min,
                    'ref_depth_pct': ref_depth_pct,
                    'disrupted_genes': disrupted, 'disrupted_gene_details': disrupted_details,
                    'n_named_genes': sum(1 for _g in disrupted
                                         if not re.match(r'^ENSCAFG\d+$', _g)),
                    })

# No artefact split any more: rarity is reported, not a mechanism we cannot
# observe. A region every reference dog is also missing comes out near 100% and
# speaks for itself.
real_regions = regions
artefact_regions = []
# Only include disrupted genes from real (non-artefact) regions
real_gene_names = {g for r in real_regions for g in r['disrupted_genes']}
real_disrupted = {k: v for k, v in disrupted_all.items() if k in real_gene_names}

win_kb = round(cnv_win/1000, 1)
ref_meta = ref_panel.get('_meta', {})
n_ref_dogs = ref_meta.get('n_dogs', '?')
# Say which reference actually produced ref_depth_pct, rather than naming one
# unconditionally — the fallback is silent otherwise.
_ref_src = (f'{cnv_panel_n}-dog panel at {round((cnv_bin or 0)/1000)}kb'
            if cnv_panel_idx else f'legacy {n_ref_dogs}-dog pooled reference at 1Mb')
cnv_out = {
    'regions': real_regions, 'disrupted_genes': list(real_disrupted.values()),
    'artefact_regions': artefact_regions,
    'events': cnv_events,
    'summary': {
        'total_regions': len(real_regions), 'unique_genes': len(real_disrupted),
        'n_nongenic_skipped': _nongenic,
        'n_events': len(cnv_events),
        'n_event_gain': sum(1 for e in cnv_events if e['direction'] == 'gain'),
        'n_event_loss': sum(1 for e in cnv_events if e['direction'] == 'loss'),
        'n_event_rare': sum(1 for e in cnv_events if e.get('novelty') == 'rare'),
        'common_pct': COMMON_PCT,
        'event_z_threshold': Z_CUT_CNV,
        'panel_n': cnv_panel_n, 'panel_bin_bp': cnv_bin,
        'method': (f'Deletions called on the reference panel grid ({round((cnv_bin or 0)/1000)}kb), '
                   f'depth normalised to genome-wide mean, ratio<0.15 for homozygous loss. '
                   f'Each region reports the share of the {cnv_panel_n}-dog panel that is at '
                   f'least as low there. Only regions disrupting a gene are listed.'),
        'event_note': (f'{len(cnv_events)} region(s) at |z| >= {Z_CUT_CNV} against the '
                       f'{cnv_panel_n}-dog panel, of which '
                       f'{sum(1 for e in cnv_events if e.get("novelty") == "rare")} are rare '
                       f'(flagged in <{COMMON_PCT}% of the reference cohort). Unlike the '
                       f'deletion list, this detects gains as well as losses. Regions common '
                       f'in the cohort are copy-number polymorphisms that many healthy dogs '
                       f'carry, not findings specific to this dog.' if cnv_panel_idx else
                       'No CNV panel available; gains were not screened for.'),
        'min_detectable_kb': round(win_kb*2), 'calling_resolution_kb': win_kb,
        'panel_note': (f'{len(real_regions)} gene-disrupting deletion region(s). Frequencies are '
                       f'the share of {cnv_panel_n} reference dogs at least as low at that point; '
                       f'a region most dogs also lack is common variation or unmappable sequence, '
                       f'not a finding about this dog.'),
    }
}
with open(f'{pub}/cnv_homdel.json', 'w') as f:
    json.dump(cnv_out, f, indent=2)
print(f"cnv_homdel.json: {len(regions)} regions, {len(disrupted_all)} disrupted genes")
PYEOF

# ── 6d: Mitochondrial lineage ───────────────────────────────
# Variants are called from a chrM-ONLY remap of the raw FASTQs, never from the
# whole-genome BAM's chrM slice: canFam4 carries 321 NUMTs covering the entire
# mito sequence, so ~98% of chrM reads in the whole-genome BAM are MAPQ<20 and
# a MAPQ-filtered caller is nearly blind there. Against a chrM-only index the
# same reads map at MAPQ 60. NUMT reads contaminate the remap instead, but mito
# depth (66-543x measured) dwarfs nuclear depth (1-7x), so they are a few
# percent of the pile — invisible to haploid consensus calling.
log "=== Stage 6d: Mitochondrial lineage (chrM remap) ==="
MITO_DIR=$OUT/mito
mkdir -p "$MITO_DIR"
$MM samtools faidx "$FASTA" chrM > "$MITO_DIR/chrM.fa"
$MM samtools faidx "$MITO_DIR/chrM.fa"
$MM bwa-mem2 index "$MITO_DIR/chrM.fa" >/dev/null 2>&1
MITO_R1=$(ls "$FASTQ_DIR"/*_R1_*.fastq.gz 2>/dev/null | sort -V)
MITO_R2=$(ls "$FASTQ_DIR"/*_R2_*.fastq.gz 2>/dev/null | sort -V)
[[ -n "$MITO_R1" && -n "$MITO_R2" ]] || die "Stage 6d: no _R1_/_R2_ FASTQs in $FASTQ_DIR"
# gunzip -c, not zcat: macOS zcat appends ".Z" to filenames and dies, which
# silently fed bwa an empty stream on Mac runs (0 reads, bogus haplogroup).
$MM bwa-mem2 mem -t "$NPROC" "$MITO_DIR/chrM.fa" <(gunzip -c $MITO_R1) <(gunzip -c $MITO_R2) 2>"$MITO_DIR/bwa.log" \
  | $MM samtools view -b -F 4 -q 20 - \
  | $MM samtools sort -o "$MITO_DIR/chrM.bam" -
$MM samtools index "$MITO_DIR/chrM.bam"
$MM bcftools mpileup -f "$MITO_DIR/chrM.fa" -d 4000 -q 20 -Q 20 -a AD,DP "$MITO_DIR/chrM.bam" 2>/dev/null \
  | $MM bcftools call -mv --ploidy 1 -Oz -o "$MITO_DIR/chrM.vcf.gz" 2>/dev/null
$MM bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\t[%AD]\n' "$MITO_DIR/chrM.vcf.gz" > "$MITO_DIR/snps.tsv"
MITO_DEPTH=$($MM samtools depth -a "$MITO_DIR/chrM.bam" | awk '{s+=$3; n++} END {printf "%.0f", (n?s/n:0)}')
log "  chrM mean depth ${MITO_DEPTH}x, $(wc -l < "$MITO_DIR/snps.tsv") raw variant rows"

MITO_DEPTH="$MITO_DEPTH" MITO_DIR="$MITO_DIR" PUB_DIR="$PUB" REFJ="$REF_JSON" SAMPLE="$DOG_NAME" \
"$DATA_PYTHON" - <<'PYEOF'
import gzip, json, os

# Assignment is SNP-set based (symmetric difference vs each labeled reference's
# substitutions relative to canFam4 chrM). Whole-sequence distance is broken
# here: a low-pass consensus inherits the reference everywhere reads don't
# disagree, so sequence distance is dominated by shared VNTR/reference noise.
refs = json.load(gzip.open(os.environ['REFJ'] + '/mito_haplogroups.json.gz', 'rt'))
mito_dir, pub = os.environ['MITO_DIR'], os.environ['PUB_DIR']
depth = int(os.environ['MITO_DEPTH'] or 0)

dog, het = set(), 0
for line in open(mito_dir + '/snps.tsv', encoding='utf-8'):
    f = line.rstrip('\n').split('\t')
    if len(f) < 5 or len(f[1]) != 1 or len(f[2]) != 1:
        continue
    pos, ref_a, alt, qual, ad = int(f[0]), f[1], f[2], float(f[3]), f[4]
    parts = [int(x) for x in ad.split(',') if x.isdigit()]
    af = parts[-1] / sum(parts) if parts and sum(parts) else 1.0
    if qual < 30:
        continue
    if 0.10 <= af < 0.90:
        het += 1
    if af >= 0.7:
        dog.add((pos - 1, alt))

if depth < 5 or (len(dog) == 0 and het == 0):
    # No usable mitochondrial reads (or an upstream failure produced an empty
    # call set): an empty variant set would spuriously 'match' the most
    # reference-like records, so refuse to assign rather than fabricate.
    out = {
        'haplogroup': None, 'nearest_haplotype': None, 'confidence': 'none',
        'suppressed': True, 'mean_chrm_depth': depth, 'n_snps_vs_reference': 0,
        'candidate_heteroplasmy_sites': 0,
        'group_story': '', 'group_frequencies': {},
        'summary': ('Maternal lineage could not be determined: no usable mitochondrial '
                    'sequence was recovered for this sample.'),
        'method': 'Suppressed — insufficient mitochondrial read data.',
    }
    with open(pub + '/mito_result.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2)
    print('mito_result.json: suppressed (depth {}x, {} snps)'.format(depth, len(dog)))
    raise SystemExit(0)

hits = sorted((len(dog ^ {tuple(s) for s in r['snps']}), r['haplotype'], r['acc'])
              for r in refs['refs'])
d0, hap0, acc0 = hits[0]
group = hap0[0]
margin = next(d for d, h, a in hits if h[0] != group) - d0
conf = 'high' if margin >= 20 else ('moderate' if margin >= 10 else 'low')

out = {
    'haplogroup': group,
    'nearest_haplotype': hap0,
    'nearest_acc': acc0,
    'snp_distance': d0,
    'margin_to_other_group': margin,
    'confidence': conf,
    'n_snps_vs_reference': len(dog),
    'mean_chrm_depth': depth,
    'candidate_heteroplasmy_sites': het,
    'heteroplasmy_note': ('Sites at intermediate allele fraction overlap the noise floor from '
                          'nuclear mitochondrial insertions (NUMTs) and are not individually '
                          'reported.'),
    'top_matches': [{'haplotype': h, 'snp_distance': d} for d, h, a in hits[:3]],
    'group_story': refs['groups'].get(group, ''),
    'group_frequencies': {'A': '~65-70% of dogs', 'B': '~20%', 'C': '~5-10%', 'D': 'rare'},
    'method': ('chrM-only remap of raw reads (whole-genome alignments are MAPQ-degraded by '
               'NUMTs), haploid consensus calling, SNP-set nearest neighbor against '
               '{} haplotype-labeled GenBank mitogenomes.'.format(refs['meta']['n_refs'])),
    'snps': sorted([p, b] for p, b in dog),
}
# ── Known mitochondrial disease variants ──
# The full mitogenome is called at 60-500x, so screening the (few) published
# canine mtDNA disease variants is essentially free. Heteroplasmy-aware: the
# report carries the fraction of reads supporting each variant.
try:
    dz = json.load(gzip.open(os.environ['REFJ'] + '/mito_disease.json.gz', 'rt'))
except Exception:
    dz = json.load(open(os.environ['REFJ'] + '/mito_disease.json'))
screen = []
rows = []
for line in open(os.environ['MITO_DIR'] + '/snps.tsv', encoding='utf-8'):
    f2 = line.rstrip('\n').split('\t')
    if len(f2) >= 5:
        rows.append(f2)
for dvar in dz['variants']:
    hit = None
    for f2 in rows:
        p2, r2, a2 = int(f2[0]), f2[1], f2[2]
        if dvar['type'] == 'del' and len(r2) > len(a2) and abs(p2 - dvar['chrm_pos']) <= dvar.get('window', 3):
            hit = f2; break
        if dvar['type'] == 'snv' and p2 == dvar['chrm_pos'] and r2 == dvar.get('ref') and a2 == dvar.get('alt'):
            hit = f2; break
    entry = {'name': dvar['name'], 'gene': dvar['gene'], 'breed': dvar['breed'],
             'description': dvar['description'], 'reference': dvar['reference'],
             'detected': bool(hit)}
    if hit:
        parts = [int(x) for x in hit[4].split(',') if x.isdigit()]
        entry['heteroplasmy_pct'] = round(100.0 * parts[-1] / sum(parts), 1) if parts and sum(parts) else None
        entry['note'] = 'Detected — discuss with your veterinarian; mitochondrial variants pass from mother to all offspring.'
    screen.append(entry)
out['disease_screen'] = screen
out['disease_screen_note'] = dz['note']

with open(pub + '/mito_result.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2)
print('mito_result.json: haplogroup {} ({}), dist {}, margin {}, depth {}x, screen {}'.format(
    group, hap0, d0, margin, depth, ['%s=%s' % (e['gene'], e['detected']) for e in screen]))
PYEOF

# ── Stage 6e: Unified 50kb CNV / aneuploidy (segmentation vs panel-of-normals)
# One framework across scales, replacing the dual 1Mb/adaptive analysis in the
# report. Requires the fixed 50kb grid (stage 5c) — absent on legacy runs, in
# which case a null result is written and the report falls back gracefully.
log "  Unified 50kb CNV segmentation…"
"$DATA_PYTHON" - << PYEOF
import gzip, json, os, sys
import numpy as np
sys.path.insert(0, "$D/analysis/cnv")
from segment_caller import call_sample

pub = "$PUB"
cov_path = "$OUT/coverage_50kb.tsv.gz"
panel_path = "$D/reference_panel/coverage_50kb_panel.npz"
out_path = pub + "/cnv50_result.json"
if not os.path.exists(cov_path):
    json.dump(None, open(out_path, "w"))
    print("cnv50_result.json: null (no 50kb grid for this run)")
else:
    panel = dict(np.load(panel_path, allow_pickle=True))
    counts = []
    with gzip.open(cov_path, "rt") as f:
        for l in f:
            counts.append(float(l.rstrip("\n").split("\t")[3]))
    v = np.array(counts, dtype=np.float64)
    chrom = panel["chrom"]
    auto = chrom != "chrX"
    med = np.median(v[auto][v[auto] > 0])
    v = v / med
    sex = "F" if float(np.median(v[~auto])) > 0.75 else "M"
    r = call_sample(v, sex, panel)

    # Gene overlap for reported segments (reference gene table; chrom keys
    # are bare numbers, X included).
    # Genome-wide gene spans (built once from the snpEff GTF — the per-sample
    # cnv_genes.json copies only cover chromosomes the legacy caller touched).
    gene_map = {}
    try:
        gene_map = json.load(open("$REF_JSON/genes_canFam4.json"))
    except Exception:
        pass
    for g in r["segments"]:
        c = g["chrom"].replace("chr", "")
        hits = [e["name"] for e in gene_map.get(c, [])
                if e.get("name") and e["start"] < g["end"] and e["end"] > g["start"]]
        g["genes"] = hits[:12]
        g["n_genes"] = len(hits)
    rare = [g for g in r["segments"] if not g.get("common_variant")]
    result = {
        "sex_inferred": sex,
        "confidence": r.get("confidence"),
        "sigma": r["sigma"],
        "aneuploidies": r["aneuploidies"],
        "segments": r["segments"],
        "n_rare": len(rare),
        "n_common": len(r["segments"]) - len(rare),
        "n_reference_dogs": 1257,
        "method": ("50kb-window coverage vs a 1,257-dog panel of normals; "
                   "recursive segmentation; chromosome-level shifts reported as "
                   "aneuploidy; events shared by >=5% of dogs marked as common "
                   "copy-number polymorphisms."),
    }
    if r.get("note"):
        result["note"] = r["note"]
    json.dump(result, open(out_path, "w"), indent=2)
    print(f"cnv50_result.json: conf={r.get('confidence')} aneuploidies={len(r['aneuploidies'])} rare={len(rare)} common={result['n_common']}")
PYEOF
fi # end stage 6

if (( FROM_STAGE <= 7 && TO_STAGE >= 7 )); then
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
  # GLIMPSE2_phase: estimates this dog's genotypes at panel positions using BAM reads + LD
  if ! $MM_GLIMPSE GLIMPSE2_phase \
    --bam-file      "$OUT/markdup.bam" \
    --reference     "$DOG10K_PANEL" \
    --fasta         "$FASTA" \
    --input-region  "$ireg" \
    --output-region "$oreg" \
    --threads 1 \
    --output "$outfile" 2>&1 | grep -v "AC/AN INFO fields" | tail -1; then
    echo "WARN ${chr}_chunk${id}: GLIMPSE2_phase failed — skipping"
    return 0
  fi
  if [ ! -f "$outfile" ]; then
    echo "WARN ${chr}_chunk${id}: no output produced — skipping"
    return 0
  fi
  $MM_GLIMPSE bcftools index -f "$outfile" || { echo "WARN ${chr}_chunk${id}: index failed — skipping"; rm -f "$outfile"; return 0; }
  echo "DONE ${chr}_chunk${id}"
}
export -f estimate_chunk
export OUT DOG10K_PANEL FASTA GENO_DIR

FAILED_CHUNKS="$(dirname "$CHUNKS_DIR")/failed_chunks.txt"
log "Estimating genotypes across chunks (${GLIMPSE_PARALLEL} parallel jobs)..."
GLIMPSE_PIDS_STR=""
glimpse_wait_slot() {
  while true; do
    local running=0 live_pids=""
    for pid in $GLIMPSE_PIDS_STR; do
      if kill -0 "$pid" 2>/dev/null; then
        running=$((running + 1))
        live_pids="${live_pids:+$live_pids }$pid"
      fi
    done
    GLIMPSE_PIDS_STR="$live_pids"
    if [ "$running" -lt "$GLIMPSE_PARALLEL" ]; then break; fi
    sleep 2
  done
}
for chunkfile in "$CHUNKS_DIR"/*.txt; do
  chr=$(basename "$chunkfile" .txt)
  while IFS=$'\t' read -r id chrom ireg oreg rest; do
    outfile="$GENO_DIR/${chr}_chunk${id}.bcf"
    if [ -f "$outfile" ] && [ -f "${outfile}.csi" ]; then
      echo "SKIP ${chr}_chunk${id}"; continue
    fi
    if grep -qxF "${chr}_chunk${id}" "$FAILED_CHUNKS" 2>/dev/null; then
      echo "SKIP(known-fail) ${chr}_chunk${id}"; continue
    fi
    glimpse_wait_slot
    estimate_chunk "$chr" "$id" "$ireg" "$oreg" &
    GLIMPSE_PIDS_STR="${GLIMPSE_PIDS_STR:+$GLIMPSE_PIDS_STR }$!"
  done < "$chunkfile"
done
kill "$MEM_POLL_PID" 2>/dev/null || true
wait "$MEM_POLL_PID" 2>/dev/null || true
wait
log "Genotype estimation complete"
( while kill -0 $$ 2>/dev/null; do
    rss_kb=$(ps -u "$(id -u)" -o rss= 2>/dev/null | awk 'BEGIN{t=0} {t+=$1} END{print t}')
    cur=$(cat "$PEAK_MEM_FILE")
    if (( rss_kb > cur )); then echo "$rss_kb" > "$PEAK_MEM_FILE"; fi
    sleep 10
  done ) &
MEM_POLL_PID=$!

log "Ligating chunks per chromosome..."
_ligated=0; _empty=""
for chr in $(ls "$CHUNKS_DIR"/*.txt | xargs -I{} basename {} .txt); do
  list=$(mktemp)
  # `|| true` matters: with `set -o pipefail`, an `ls` that matches nothing exits
  # non-zero and `set -e` kills the script before the emptiness guard below can
  # skip the chromosome. That happens whenever every chunk of a chromosome was
  # skipped or failed at runtime.
  ls "$GENO_DIR"/${chr}_chunk*.bcf 2>/dev/null | sort -V > "$list" || true
  [ -s "$list" ] || { rm -f "$list"; _empty="$_empty $chr"; continue; }
  $MM_GLIMPSE GLIMPSE2_ligate --input "$list" --output "$LIGATED_DIR/${chr}.bcf" \
    || die "GLIMPSE2_ligate failed for $chr"
  $MM_GLIMPSE bcftools index -f "$LIGATED_DIR/${chr}.bcf"
  rm -f "$list"
  _ligated=$((_ligated + 1))
done
log "Ligated $_ligated chromosomes"
if [ -n "$_empty" ]; then
  # Not fatal — the merge still succeeds — but the final BCF is missing these
  # chromosomes entirely, so say so rather than letting it pass unnoticed.
  log "WARNING: no genotypes for:$_empty — these are ABSENT from the imputed BCF"
fi

log "Merging chromosomes..."
$MM_GLIMPSE bcftools concat \
  $(ls "$LIGATED_DIR"/chr*.bcf | sort -V) \
  -O b -o "$IMPUTED_BCF"
$MM_GLIMPSE bcftools index -f "$IMPUTED_BCF"
TOTAL=$($MM_GLIMPSE bcftools stats "$IMPUTED_BCF" 2>/dev/null | grep "^SN.*number of records" | awk '{print $NF}')
log "Imputed BCF: $TOTAL variants → $IMPUTED_BCF"
fi # end stage 7

if (( FROM_STAGE <= 8 && TO_STAGE >= 8 )); then
# ── Stage 8: OMIA genotyping from Dog10K imputed panel ───────
log "=== Stage 8: OMIA genotyping (Dog10K imputed + BAM fallback) ==="
# markdup.bam may have been cleaned up; sites.bam is the permanent extract of
# the same reads at known-variant regions. NEVER run without one of them — a
# BAM-less rerun silently rewrites every read-based call as "insufficient
# reads" (happened to 92 published dogs on 2026-08-29).
CALL_BAM="$OUT/markdup.bam"
[[ -f "$CALL_BAM" ]] || CALL_BAM="$OUT/sites.bam"
[[ -f "$CALL_BAM" ]] || die "Stage 8 needs reads: neither markdup.bam nor sites.bam in $OUT — refusing to write a degraded report"
"$DATA_PYTHON" - << PYEOF
import subprocess, pysam, json, re

BCF       = "$IMPUTED_BCF"
BAM       = "$CALL_BAM"
OMIA      = "$OMIA_DB"
QC_JSON   = "$PUB/qc_result.json"
PUB       = "$PUB"
DOG       = "$DOG_LOWER"
GP_HIGH   = 0.90   # minimum max(GP) to report a GLIMPSE2 call
with open(QC_JSON) as f:
    mean_depth = json.load(f)['genome_mean_depth']
# The BAM fallback used to be gated on WHOLE-GENOME depth >= 10x, which turned
# it off for every dog in this low-pass product — the known-variants tab showed
# 430 of 479 variants as untestable on a 7x sample. The gate was redundant:
# gt_from_bam already refuses any SITE with fewer than 5 informative reads and
# grades each call high/medium/low by the reads actually present, which is the
# right granularity. At 7x mean, ~70% of sites clear 5 reads; the ones that do
# not are individually skipped rather than the whole assay being switched off.
print(f"Mean depth: {mean_depth:.1f}x — BAM fallback per-site (>=5 reads at the position)")

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

_fasta = pysam.FastaFile("$FASTA")
def gt_from_bam(chrom, pos, ref, alt, min_bq=20, min_mq=20):
    # Refuse to call unless the assembly base at this position IS the variant's
    # stated reference allele. The OMIA catalogue mixes strands and assembly
    # versions, and where its "alt" is actually the canFam4 reference base,
    # every read in every dog supports "alt" — counting bases then reports a
    # healthy animal as homozygous for haemophilia. Of the first 13 affected
    # calls this fallback produced, 12 were exactly that (0 ref reads, alt ==
    # assembly base); the 13th, MC1R e/e with the ref matching, was real.
    # The imputed path is immune because the BCF's ref/alt must match the
    # assembly; this check gives the direct path the same protection.
    try:
        base = _fasta.fetch(chrom, pos - 1, pos).upper()
    except Exception:
        return None
    if base != ref.upper():
        return {'zygosity': 'not_callable', 'affected': False,
                'call_confidence': 'none', 'source': 'assembly_mismatch',
                'note': (f'OMIA ref {ref} does not match canFam4 base {base} at this '
                         f'position (strand or assembly-version discrepancy); '
                         f'cannot be genotyped from reads')}
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
    # Floor lowered from 5 informative reads to 2 so a low-pass dog reports
    # most of the catalogue rather than listing it untested. 2-4 reads is
    # honest as a LOW-confidence call: with per-base error ~1% after the
    # BQ>=20 filter, two concordant reads mis-genotype at well under 1%,
    # while het vs hom remains genuinely uncertain — the zygosity_note says
    # so, and affected low-confidence calls tell the owner to confirm with a
    # targeted test before acting.
    if total < 2 or n_v < 2:
        return None
    f_alt = n_alt / n_v
    zyg = 'ref/ref' if f_alt < 0.1 else ('alt/alt' if f_alt > 0.9 else 'het')
    if n_v < 5:
        # At 2-4 reads a het can easily sample only one allele: an all-ref or
        # all-alt pile does not exclude het, so only the presence calls are
        # trustworthy; zygosity is indicative.
        conf = 'low'
    else:
        conf = 'high' if total >= 20 else ('medium' if total >= 10 else 'low')
    out = {'zygosity': zyg, 'depth': total, 'ref_count': n_ref, 'alt_count': n_alt,
           'affected': zyg in ('alt/alt', 'het'), 'call_confidence': conf,
           'source': 'bam_direct'}
    if n_v < 5 and out['affected']:
        out['note'] = ('Called from only {} reads — low confidence; confirm with a '
                       'targeted DNA test before acting on this.').format(n_v)
    return out

def gt_sv_from_bam(chrom, pos, svtype='', pos_end=None, bnd_window=15, min_mq=20):
    """Genotype a known structural variant (del/ins/dup/inv/untyped) at a
    known position from read alignments.

    Two orthogonal read signatures, per event class:
      - deletions shorter than a read appear as D operations in the CIGAR of
        every read crossing them;
      - any breakpoint (insertion, duplication junction, inversion, or a
        deletion longer than a read) truncates alignments, so reads soft-clip
        at the breakpoint coordinate.
    Reads that span the site cleanly with margin count as reference support.
    Catalogue coordinates are approximate to a few bases (target-site
    duplications, ambiguous naming), hence the +-bnd_window tolerance —
    same approach validated for MDR1 (4-bp del) and merle (PMEL SINE).
    """
    # Size-awareness matters: the catalogue mixes 1-4bp deletions (visible as
    # D operations inside reads) with deletions of tens of KILOBASES (never
    # visible as D — only as soft-clips at the two breakpoints). Counting any
    # overlapping D as the pathogenic allele mis-called a normal Poodle
    # homozygous for a 130kb dwarfism deletion off a benign 1bp indel.
    exp_len = (pos_end - pos) if (pos_end and pos_end > pos) else None
    large = exp_len is not None and exp_len > 60
    breakpoints = [pos, pos_end] if large else [pos]
    # Evidence is counted per FRAGMENT, not per alignment: overlapping mates of
    # one pair are a single molecule, and counting them twice let a lone
    # chimeric fragment (both mates 62S19M30S at the same spot) satisfy the
    # 2-read floor and call a healthy dog alt/alt for a RELN deletion.
    frag_del, frag_clip, frag_span = set(), set(), set()
    try:
        bam_fh = pysam.AlignmentFile(BAM, 'rb')
        for bp in breakpoints:
            lo, hi = bp - bnd_window, bp + bnd_window
            for r in bam_fh.fetch(chrom, max(0, bp - 200), bp + 200):
                if r.is_unmapped or r.mapping_quality < min_mq or r.is_secondary or r.is_supplementary:
                    continue
                start, end = r.reference_start + 1, r.reference_end
                cig = r.cigartuples or []
                aligned = sum(ln for op, ln in cig if op in (0, 7, 8))
                has_del = False
                if not large:
                    rp = start
                    for op, ln in cig:
                        if op in (0, 7, 8):
                            rp += ln
                        elif op == 2:
                            if rp - 1 <= hi and rp + ln >= lo:
                                # D length must be consistent with the expected
                                # deletion; unknown-size dels require >=2bp so a
                                # ubiquitous 1bp homopolymer slip cannot call.
                                if exp_len is not None:
                                    if exp_len // 2 <= ln <= exp_len * 2:
                                        has_del = True
                                elif ln >= 2:
                                    has_del = True
                            rp += ln
                        elif op == 3:
                            rp += ln
                clip_left  = cig and cig[0][0] in (4, 5) and cig[0][1] >= 8 and lo <= start <= hi
                clip_right = cig and cig[-1][0] in (4, 5) and cig[-1][1] >= 8 and lo <= end <= hi
                # A mostly-clipped alignment is mapping noise, not a breakpoint:
                # a <30bp anchor can land anywhere in the genome.
                if (clip_left or clip_right) and aligned < 30:
                    clip_left = clip_right = False
                name = r.query_name
                if has_del:
                    frag_del.add(name)
                elif clip_left or clip_right:
                    frag_clip.add(name)
                elif start <= lo - 5 and end >= hi + 5:
                    frag_span.add(name)
        bam_fh.close()
    except Exception:
        return None
    # A fragment with any alt-supporting alignment is alt; only fragments with
    # no alt evidence at all count as reference support.
    n_del  = len(frag_del)
    n_clip = len(frag_clip - frag_del)
    n_span = len(frag_span - frag_del - frag_clip)
    t = (svtype or '').lower()
    if 'del' in t:
        n_alt = n_del + n_clip          # short dels cigar; long dels clip
    elif t:
        n_alt = n_clip                  # ins/dup/inv: junction clips only
    else:
        n_alt = n_del + n_clip          # untyped: any structural evidence
    n_v = n_alt + n_span
    if n_v < 2:
        return None
    f_alt = n_alt / n_v
    zyg = 'ref/ref' if f_alt < 0.1 else ('alt/alt' if f_alt > 0.9 else 'het')
    conf = 'low' if n_v < 5 else ('high' if n_v >= 20 else ('medium' if n_v >= 10 else 'low'))
    if zyg != 'ref/ref' and n_alt < 2:
        return None                      # a single anomalous read is not a call
    out = {'zygosity': zyg, 'depth': n_v, 'ref_count': n_span, 'alt_count': n_alt,
           'affected': zyg in ('alt/alt', 'het'), 'call_confidence': conf,
           'source': 'bam_direct_sv', 'sv_evidence': {'del_reads': n_del, 'clip_reads': n_clip}}
    if out['affected'] and conf != 'high':
        out['note'] = ('Structural change seen in {} of {} reads — {} confidence; confirm with a '
                       'targeted DNA test before acting on this.').format(n_alt, n_v, conf)
    if not t:
        out['note'] = (out.get('note', '') + ' Variant type not specified in the catalogue; '
                       'this call reports structural evidence at the position.').strip()
    return out

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
        # Deletions ARE resolvable from reads: an aligner marks a deletion in
        # the CIGAR of every read spanning it, so counting deletion-carrying
        # vs clean-spanning reads genotypes it like a SNV pileup. Added for
        # MDR1 (ABCB1 4-bp del, chr14:13704489) — the most-requested variant
        # in the catalogue — but applies to any del with coordinates.
        # Insertions and complex events stay uncalled.
        del_call = None
        if pos and chrom:
            del_call = gt_sv_from_bam(chrom, int(pos), v.get('variant_type') or '', v.get('pos_end'))
        if del_call:
            new_v[DOG] = del_call
            n_bam += 1
        else:
            new_v[DOG] = {'zygosity': 'indel_no_call', 'affected': False,
                          'call_confidence': 'none', 'source': 'indel'}
            n_indel += 1

    else:
        gt, af, gp, in_panel = query_glimpse2(chrom, int(pos), ref, alt)

        if in_panel:
            # Site is in Dog10K panel — only report if GP is high-confidence
            n_panel += 1
            max_gp = max(gp) if gp else 0.0
            # Report the genotype at every panel site, graded by posterior
            # rather than hidden below a single cutoff: >=0.90 high,
            # >=0.70 medium, >=0.50 low. Below 0.50 the posterior barely
            # prefers one genotype over another and no call is honest.
            # (Previously anything under 0.90 was suppressed entirely, which
            # left most of the catalogue "untested" for sub-2x dogs.)
            if max_gp >= 0.50:
                alleles = re.split(r'[|/]', gt)
                zyg = ('alt/alt' if set(alleles) == {'1'} else
                       'ref/ref' if set(alleles) == {'0'} else 'het')
                conf = 'high' if max_gp >= GP_HIGH else ('medium' if max_gp >= 0.70 else 'low')
                call = {'zygosity': zyg, 'affected': zyg in ('alt/alt', 'het'),
                        'call_confidence': conf, 'glimpse2_gt': gt,
                        'glimpse2_gp': gp, 'source': 'dog10k_imputed'}
                if conf != 'high' and call['affected']:
                    call['note'] = ('Genotype posterior {:.2f} — below the high-confidence threshold; '
                                    'confirm with a targeted DNA test before acting on this.').format(max_gp)
                if af is not None:
                    call['af_dog10k'] = round(af, 4)
            else:
                call = {'zygosity': 'low_gp_no_call', 'affected': False,
                        'call_confidence': 'low', 'glimpse2_gt': gt,
                        'glimpse2_gp': gp, 'source': 'dog10k_imputed',
                        'note': f'max GP {max_gp:.2f} too uncertain to call'}
            new_v[DOG] = call

        else:
            # Not in the imputation panel: read the genotype directly from the
            # aligned reads at that position. gt_from_bam grades its own
            # confidence from the site depth and returns None below 5 reads.
            bam_call = gt_from_bam(chrom, int(pos), ref, alt)
            if bam_call and bam_call.get('source') == 'assembly_mismatch':
                new_v[DOG] = bam_call
                n_not_callable += 1
            elif bam_call:
                new_v[DOG] = bam_call
                n_bam += 1
            else:
                new_v[DOG] = {'zygosity': 'no_call', 'affected': False,
                              'call_confidence': 'none',
                              'source': 'not_in_panel_insufficient_reads',
                              'note': 'Fewer than 5 usable reads at this position'}
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
        'bam_fallback_used': True,
    },
    'method': (
        f'Primary: GLIMPSE2 Dog10K imputed panel (30.4M SNPs); high-confidence calls require max GP ≥ {GP_HIGH}. '
        f'SNVs not in the panel are read directly from the aligned reads at that position '
        f'(minimum 5 reads; confidence graded per site by depth — these calls may be lower '
        f'confidence than imputed ones).'
    ),
    'variants': variants,
}
with open(f'{PUB}/omia_result.json', 'w') as f:
    json.dump(result, f, indent=2)
print(f"omia_result.json: {len(variants)} variants | panel={n_panel} bam={n_bam} not_callable={n_not_callable} indels={n_indel}")
print(f"  affected SNVs: {affected_snv} ({high_conf} high confidence)")
PYEOF
fi # end stage 8

if (( FROM_STAGE <= 9 && TO_STAGE >= 9 )); then
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
"$DATA_PYTHON" - << PYEOF
import subprocess, tempfile, os, re, numpy as np, json
from scipy.optimize import nnls
from scipy.linalg import solve_triangular

BCF     = "$IMPUTED_BCF"
SITES   = "$BREED_SITES"
PHAT    = "$BREED_PHAT"
BREEDS  = "$BREED_LABELS"
LASSO   = float("$BREED_LASSO")
PUB     = "$PUB"
DOG     = "$DOG_NAME"

# Load Parker SNP positions and build lookup: (chr_with_prefix, pos) → row index
# Sites come from the panel itself, in Phat row order, carrying the a1/a2
# orientation its frequencies were computed against. Deriving them from a .bim
# instead would risk a silent row/allele mismatch against Phat.
parker_snps = []
pos_index = {}
with open(SITES) as f:
    next(f)
    for line in f:
        c, pos_s, a1, a2 = line.rstrip('\n').split('\t')
        pos = int(pos_s)
        parker_snps.append({'chrom': c, 'pos': pos, 'a1': a1, 'a2': a2})
        pos_index[(c, pos)] = len(parker_snps) - 1
n_snps = len(parker_snps)
print(f"Breed panel: {n_snps} SNPs to query")

# Write a BED file for bcftools -R (avoids ARG_MAX on 100k+ positions)
bed_fh = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
for s in parker_snps:
    bed_fh.write(f"{s['chrom']}\t{s['pos']-1}\t{s['pos']}\n")
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
# Refuse to produce a profile from too little data. With zero coverage NNLS
# happily returns a near-uniform vector, which looks like a real result and is
# not — that is exactly what got written for 96 dogs when a resume could not
# find its BCF. 20% is well below any healthy sample (real dogs run ~86%).
if pct_covered < 20:
    raise SystemExit(
        f"ERROR: only {valid.sum()}/{n_snps} panel SNPs genotyped "
        f"({pct_covered:.1f}%). Refusing to write a breed profile. "
        f"Check that {BCF} exists and is the imputed BCF for this dog.")
print(f"Imputed dosages: {valid.sum()}/{n_snps} Parker SNPs ({pct_covered:.1f}%)")

# Allele frequency matrix, one column per breed, rows aligned to SITES.
# Stored as .npy rather than text: 115MB and instant to load, against ~170MB
# and ~30s for np.loadtxt on the same numbers.
P = np.load(PHAT)
print(f"Phat shape: {P.shape}")

# Breed labels, one per Phat column. Already canonical and harmonised by
# analysis/breed_accuracy/harmonize.py — Parker and Dog10K code vocabularies
# collapsed onto one biology, village dogs grouped by region, wolves pooled.
# So no merging happens here any more; the panel arrives ready to use.
breed_labels = [l.strip() for l in open(BREEDS) if l.strip()]
breed_counts = {}
assert P.shape == (n_snps, len(breed_labels)), \
    f"panel mismatch: Phat {P.shape} vs {n_snps} sites x {len(breed_labels)} breeds"
print(f"Breed labels: {len(breed_labels)} (first 3: {breed_labels[:3]})")

# NNLS supervised projection: find q s.t. P_valid @ q ≈ dosages_valid, q ≥ 0
P_v = P[valid, :]          # restrict to genotyped sites
x_v = dosages[valid]

# Constrained NNLS with sum-to-1 row appended.
#
# That row is effectively inert — one row against ~131k data rows — and the
# simplex is really imposed by the normalisation two lines below. Do NOT
# "fix" this by weighting the row so it binds: P holds allele frequencies in
# [0,1] while the dosages are in [0,2], so the correct scale is sum(q) ~ 2 and
# forcing sum(q) = 1 is a misspecified constraint. Weighting it costs real
# accuracy on held-out dogs — 95.09% -> 92.55% at w=1e4.
#
# Fitting x/2 instead makes sum(q)=1 correct, and then the constraint weight
# makes no difference at all (95.09% at every weight tried).
#
# Regularisation, on the other hand, splits: RIDGE is neutral to harmful (it
# spreads mass across correlated breeds, which is exactly the leakage we want
# less of), but a non-negative LASSO helps materially on MIXED dogs — 27% lower
# proportion error against simulated crosses of known composition. Purebred
# top-1 is a metric that cannot see this, which is why the first sweep missed
# it. Adding the penalty is free: with q >= 0 it is linear, so
#     ||Aq-b||^2 + lam*1'q  =  q'Gq - 2(A'b - lam/2)'q + const
# i.e. the same Gram matrix with a shifted right-hand side. Beware a sharp
# cliff — lam=1 is near-optimal, lam=3 collapses the solution entirely.
# See analysis/breed_accuracy/{regularization_sweep,lasso_sweep}.py.
A = np.vstack([P_v, np.ones((1, P_v.shape[1]))])
b = np.hstack([x_v, [1.0]])
if LASSO:
    # Non-negative L1. With q >= 0 the penalty is linear, so it is a shift of
    # the normal-equation right-hand side rather than a different solver:
    #   ||Aq-b||^2 + lam*1'q = q'(A'A)q - 2(A'b - lam/2)'q + const
    G = (A.T @ A).astype(np.float64)
    rhs = (A.T @ b).astype(np.float64) - 0.5 * LASSO * float(np.mean(np.diag(P_v.T @ P_v)))
    G[np.diag_indices_from(G)] += 1e-6
    R = np.linalg.cholesky(G).T
    q_raw, _ = nnls(R, solve_triangular(R.T, rhs, lower=True))
else:
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
    # NELK is Norwegian ELKhound, not Norrbottenspets. Confirmed empirically:
    # a known Norwegian Elkhound sample projects onto NELK. The panel has no
    # Norrbottenspets population, and the code NORW ('Norwegian Elkhound' in
    # this table) is not in the panel at all.
    'NELK':'Norwegian Elkhound','NEWF':'Newfoundland',
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
    # SSHP is Shetland Sheepdog, not Smooth Collie — confirmed by Dog10K Table S1
    # (Nature/Genome Biology supplement), which lists every SSHP sample as
    # "Shetland Sheepdog". Same class of error as NELK: this name table carries
    # 236 codes for a 177-population panel and some panel codes were given the
    # wrong name out of that superset. Shetland Sheepdog is a top-30 AKC breed
    # whose apparent absence from the panel was the tell.
    'SSHP':'Shetland Sheepdog','SSNZ':'Standard Schnauzer',
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
    'WOLF':'Gray Wolf',   # the seven regional wolf populations, pooled
}

# ── AKC display grouping (PRESENTATION level) ────────────────────────
# Two separate problems, one mechanism: group by the name a customer should see
# and sum the proportions.
#
# 1. Same breed, two reference cohorts. AUSS/AUST are both "Australian
#    Shepherd", as are BELS/BMAL, BRIT/BRTR, CHIN/CRES and TIBM/TIBT. Showing a
#    customer two "Australian Shepherd" rows is simply wrong. These are caught
#    automatically by grouping on the display name.
#
# 2. AKC treats separately-bred varieties as ONE breed. Poodles are bred within
#    size variety so they form real genetic clusters, and Parker is right to
#    keep them apart for the MODEL — but a Labradoodle owner should see
#    "Poodle 59%", not Standard 27% / Miniature 20% / Toy 12%. Leave-one-out
#    supports this: Miniature Poodle reference dogs get called Toy Poodle, so
#    the panel cannot reliably separate the varieties anyway.
#
# Only merged where AKC calls it one breed. Schnauzers stay separate (Giant,
# Miniature and Standard are three distinct AKC breeds), as do Parson vs Russell
# Terrier and Bull Terrier vs Miniature Bull Terrier.
# Breed labels are already canonical names (STANDARD_POODLE, VILLAGE_EastAsia).
# This groups them into what a customer should actually read, and sums the
# proportions. Only where AKC treats the varieties as ONE breed — Schnauzers
# stay separate, and Shetland Sheepdog is NOT folded into Collie (SSHP was
# mislabelled "Smooth Collie" in the old Parker table; Dog10K Table S1 confirms
# it is a Sheltie).
AKC_DISPLAY = {
    'STANDARD_POODLE': 'Poodle', 'MINIATURE_POODLE': 'Poodle',
    'TOY_POODLE': 'Poodle',
    # Old German Shepherd is a coat variant, not a breed AKC or FCI recognises,
    # and the panel cannot separate it: its held-out dog is called German
    # Shepherd with 0% assigned to its own population.
    'OLD_GERMAN_SHEPHERD': 'German Shepherd', 'GERMAN_SHEPHERD': 'German Shepherd',
    # White Swiss Shepherd is deliberately NOT folded in with them. It is a
    # separate FCI breed and, unlike the Poodle varieties, the panel tells it
    # apart: held-out White Swiss dogs are called White Swiss 2/2 (74% own
    # share), held-out German Shepherds are called German Shepherd 2/2 (78%).
    # Merging would erase working discrimination. Note a pure GSD carries ~19%
    # White Swiss cross-loading, so only a materially higher share indicates
    # real White Swiss ancestry.
}

def _pretty(label):
    """canonical panel label -> the string a customer sees"""
    if label == 'GRAY_WOLF':
        return 'Gray Wolf'
    if label.startswith('VILLAGE_'):
        region = re.sub(r'(?<!^)(?=[A-Z])', ' ', label[len('VILLAGE_'):])
        return f'Village Dog ({region})'
    return label.replace('_', ' ').title()


def _display(code):
    return AKC_DISPLAY.get(code) or _pretty(code)

_grouped = {}
_best = {}
for _b, _s in zip(breed_labels, q.tolist()):
    d = _display(_b)
    g = _grouped.setdefault(d, {'breed': _b, 'breed_name': d, 'proportion': 0.0,
                                'components': []})
    g['proportion'] += _s
    g['components'].append({'code': _b, 'proportion': round(_s, 6)})
    if _s > _best.get(d, -1.0):
        _best[d] = _s
        g['breed'] = _b        # keep the dominant contributing code for provenance

# Drop components the fit assigned no mass. The lasso makes the solution
# genuinely sparse — typically ~14 non-zero of 230 — so without this the
# dashboard's expandable "Other ancestry" list renders 200+ rows of 0.000%.
# A component at exactly zero carries no information, in the report or in
# breed_composition_raw.
_MIN_REPORTED = 5e-5
_display_list = sorted((g for g in _grouped.values() if g['proportion'] >= _MIN_REPORTED),
                       key=lambda g: -g['proportion'])
for g in _display_list:
    g['proportion'] = round(g['proportion'], 6)
    g['components'] = sorted(g['components'], key=lambda c: -c['proportion'])
    if len(g['components']) == 1:
        del g['components']    # only carry provenance where something was merged

breed_result = {
    # Grouped for display: what the dashboard renders (it shows breed_name and
    # takes the top 6). breed_composition_raw keeps the per-population values.
    'breed_composition': _display_list,
    'breed_composition_raw': [{'breed': b,
                               'breed_name': _pretty(b),
                               'proportion': round(s, 6)}
                              for b, s in sorted(zip(breed_labels, q.tolist()), key=lambda x: -x[1])
                              if s >= _MIN_REPORTED],
    'snps_used': int(valid.sum()),
    'k': P.shape[1],
    'pct_parker_covered': round(pct_covered, 1),
    'reference_panel': (f'Parker 2017 (Cell Reports) — {P.shape[0]:,} SNPs, {P.shape[1]} populations '
                        'from 177 (regional Salukis pooled into one; the seven '
                        'geographic wolf populations pooled into one)'),
    'method': ('Supervised SCOPE NNLS projection onto Parker 2017 allele frequency matrix. '
               f'P matrix: {P.shape[0]} SNPs × {P.shape[1]} breeds. '
               'Dosages from GLIMPSE2 Dog10K imputed BCF (posterior GP-weighted, E[a1]).')
}
with open(f'{PUB}/breed_result.json', 'w') as f:
    json.dump(breed_result, f, indent=2)
print("breed_result.json written")
PYEOF

# ── 9b: Relatives (KING-robust kinship vs reference dogs) ───
# KING-robust, not a GRM: with global allele frequencies a GRM scores two
# unrelated Boxers at phi ~0.9 purely from breed homogeneity, while KING keeps
# same-breed background below ~0.05 — under the 2nd-degree threshold (0.0884).
# Validated on this cohort: same-dog kit pairs 0.499, cross-platform same-dog
# 0.489, known sib groups 0.20-0.26, Cosmo vs Luna -0.03 (unrelated).
log "=== Stage 9b: Relatives (KING kinship vs $((2031)) reference dogs) ==="
REL_DS="$OUT/relatives_ds.tsv"
awk 'NR>1{print $1"\t"$2}' "$BREED_SITES" > "$OUT/relatives_sites.pos"
$MM bcftools query -R "$OUT/relatives_sites.pos" -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' \
    "$IMPUTED_BCF" > "$REL_DS"
log "  dosages at $(wc -l < "$REL_DS") reference sites"

REL_DS="$REL_DS" REL_REF="$D/relatives_ref" PUB_DIR="$PUB" SAMPLE="$DOG_NAME" \
"$DATA_PYTHON" - <<'PYEOF'
import gzip, json, os
import numpy as np

# Depth gate, calibrated by downsampling cosmo3 (2026-08-28): at 0.3x KING
# still finds every true duplicate at phi 0.47 with zero false matches, but at
# 0.1x self-detection is lost entirely while FOUR spurious matches appear at up
# to phi 0.437 — imputation collapse fabricates near-identical relatives.
# Below the gate we publish an empty, honest result rather than noise.
MIN_DEPTH = 0.3
try:
    _qc = json.load(open(os.environ['PUB_DIR'] + '/qc_result.json'))
    _depth = float(_qc.get('genome_mean_depth') or 0)
except Exception:
    _depth = None
if _depth is not None and _depth < MIN_DEPTH:
    out = {
        'n_reference_dogs': 0, 'n_sites_used': 0, 'matches': [], 'n_matches': 0,
        'suppressed': True,
        'summary': ('Relative matching is not reported for this sample: sequencing depth '
                    '({:.1f}x) is below the {}x minimum at which kinship estimates are '
                    'reliable for imputed genotypes.').format(_depth, MIN_DEPTH),
        'method': 'Suppressed by depth gate (calibrated on downsampled truth data).',
    }
    with open(os.environ['PUB_DIR'] + '/relatives_result.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2)
    print('relatives_result.json: suppressed (depth {}x < {}x)'.format(_depth, MIN_DEPTH))
    raise SystemExit(0)

ref_dir = os.environ['REL_REF']
R = np.load(ref_dir + '/geno.npy')                       # sites x refs, int8, -1 = missing
meta = json.load(gzip.open(ref_dir + '/meta.json.gz', 'rt'))
keys = {}
with open(ref_dir + '/keys.tsv', encoding='utf-8') as f:
    for i, line in enumerate(f):
        keys[tuple(line.rstrip('\n').split('\t'))] = i

q = np.full(R.shape[0], -1, dtype=np.int8)
for line in open(os.environ['REL_DS'], encoding='utf-8'):
    p = line.rstrip('\n').split('\t')
    i = keys.get((p[0], p[1], p[2], p[3]))
    if i is None:
        continue
    ds = float(p[4])
    q[i] = 0 if ds < 0.5 else (1 if ds < 1.5 else 2)

use = q >= 0
Rs, qs = R[use], q[use]
valid = Rs >= 0                                           # per-ref non-missing mask
het_rows, homA_rows, homB_rows = qs == 1, qs == 0, qs == 2
Nhh   = (Rs[het_rows] == 1).sum(axis=0)
Nopp  = (Rs[homA_rows] == 2).sum(axis=0) + (Rs[homB_rows] == 0).sum(axis=0)
nhet_r = (Rs == 1).sum(axis=0)
nhet_q = valid[het_rows].sum(axis=0)                      # query hets over each ref's sites
denom = nhet_q + nhet_r
phi = np.where(denom > 0, (Nhh - 2.0 * Nopp) / np.maximum(denom, 1), -1.0)

t = meta['thresholds']
def category(v):
    # The 'distant' tier (0.045-0.088) is deliberately NOT reported: pilot runs
    # showed a recurring artifact trio of high-heterozygosity reference dogs
    # polluting that band even for 1x queries. Second degree and closer only.
    if v >= t['self']: return 'same_dog'
    if v >= t['first_degree']: return 'first_degree'
    if v >= t['second_degree']: return 'second_degree'
    return None

sample_name = os.environ.get('SAMPLE', '').lower()
matches = []
for j in np.argsort(phi)[::-1]:
    cat = category(float(phi[j]))
    if cat is None or len(matches) >= 10:
        break
    s = meta['samples'][j]
    # Reference dogs that ARE this sample are not a finding — skip the self
    # row by id, AND (user decision 2026-09-01) skip anything at self-level
    # kinship: a duplicate kit of the same dog listed as a "relative" is
    # confusing, and phi >= 0.45 is genetically the same dog regardless of id.
    if s['id'].lower() == sample_name or cat == 'same_dog':
        continue
    # Privacy: reference-dog identifiers never reach the report — only breed,
    # source and relationship tier. Panel dogs are public research data.
    entry = {
        'kinship': round(float(phi[j]), 3),
        'category': cat,
        'breed': s['label'],
        'source': 'Dog10K research panel' if s['source'] == 'dog10k' else 'ProsperK9 reference cohort',
    }
    # Reference dogs processed through the full pipeline carry their own
    # report-derived descriptors — still anonymous, but far more vivid than a
    # bare breed label ("a 24 kg black Labrador Retriever").
    if s.get('weight_kg'): entry['weight_kg'] = s['weight_kg']
    if s.get('coat'): entry['coat'] = s['coat']
    matches.append(entry)

out = {
    'n_reference_dogs': len(meta['samples']),
    'n_sites_used': int(use.sum()),
    'matches': matches,
    'n_matches': len(matches),
    'summary': ('No relatives detected among the reference dogs.' if not matches else
                '{} match(es) at third-degree kinship or closer.'.format(len(matches))),
    'thresholds': t,
    'method': meta['method'] + ' Same-breed background stays below the second-degree '
              'threshold, so listed matches reflect genuine kinship, not shared breed.',
}
with open(os.environ['PUB_DIR'] + '/relatives_result.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2)
top = matches[0] if matches else None
print('relatives_result.json: {} matches{}'.format(
    len(matches), '' if not top else '; top {} {} ({})'.format(top['kinship'], top['category'], top['breed'])))
PYEOF

fi # end stage 9

if (( FROM_STAGE <= 10 && TO_STAGE >= 10 )); then
# ── Stage 10: Functional annotation (SnpEff) ────────────────
log "=== Stage 10: SnpEff annotation ==="
ANN_DIR="$OUT/snpeff"
mkdir -p "$ANN_DIR"

# snpEff wrapper (conda python script) needs python + java in PATH
# openjdk lands in lib/jvm/bin on macOS conda but bin/ on Linux conda — check
# both inside the pinned env before falling back to whatever is on PATH.
SNPEFF_JAVA="$(find "$ENV_GENOMICS/lib/jvm/bin" -name java 2>/dev/null | head -1 || true)"
[ -z "$SNPEFF_JAVA" ] && [ -x "$ENV_GENOMICS/bin/java" ] && SNPEFF_JAVA="$ENV_GENOMICS/bin/java"
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
"$DATA_PYTHON" - << PYEOF
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
"$DATA_PYTHON" - << PYEOF
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
fi # end stage 10

if (( FROM_STAGE <= 11 && TO_STAGE >= 11 )); then
# ── Stage 11: PRS from imputed dosages ──────────────────────
log "=== Stage 11: PRS (imputed Parker dosages) ==="
"$DATA_PYTHON" - << PYEOF
import subprocess, numpy as np, json, csv, io, tempfile, os, gzip

BCF      = "$IMPUTED_BCF"
BIM      = "$PARKER_BIM"
FAM      = "$PARKER_FAM"
BED      = "$PARKER_BED"
PUB      = "$PUB"
REF_PRS  = "$REF_JSON/prs_reference.json"

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

# ── LMM betas (GEMMA) for the benchmarked traits ──────────────────────────
# Replaces on-the-fly marginal-ridge GWAS for the 14 AKC behaviour traits and
# height/weight. Benchmarked 2026-08-16, 5-fold breed-blocked CV on held-out
# breeds: LMM betas beat the ridge on 12-13/16 traits — the ridge's
# structure-confounded betas lost even on purebreds.
# DENSE (all ~131k SNPs, no p-cut): sparse CV-chosen cuts optimised
# breed-level ranking but were brittle for individual mixed-breed dogs
# (58-SNP weight score predicted a 13.6 kg dog at 44.6 kg; dense: 17.0 kg).
# Dense loses almost nothing at breed level (mean r 0.437 vs 0.444).
# Artifact: analysis/prs_lmm/ (GEMMA -lmm, GRM + genotype-source covariate,
# merged Parker+Dog10K panel; shared site table + per-trait beta vectors).
with gzip.open("$REF_JSON/prs_lmm_betas.json.gz", 'rt') as _f:
    LMM = json.load(_f)
LMM_SITES = LMM['sites']   # [chrom, pos, effect_allele, other_allele, panel_af]

# One BCF query for every artifact site, GP -> E[alt copies].
_bedf = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
for _sn in LMM_SITES:
    _bedf.write(f"{_sn[0]}\t{_sn[1]-1}\t{_sn[1]}\n")
_bedf.close()
_q = subprocess.run(['bcftools', 'query', '-R', _bedf.name,
                     '-f', '%CHROM\t%POS\t%REF\t%ALT\t[%GP]\n', BCF],
                    capture_output=True, text=True)
lmm_dose = {}
for _line in _q.stdout.splitlines():
    _f2 = _line.split('\t')
    if len(_f2) != 5: continue
    _gp = _f2[4].split(',')
    if len(_gp) != 3: continue
    try:
        _g = [float(x) for x in _gp]
    except ValueError:
        continue
    lmm_dose[(_f2[0], int(_f2[1]))] = (_f2[2], _f2[3], _g[1] + 2.0*_g[2])
print(f"LMM artifact: {len(LMM['traits'])} traits, {len(LMM_SITES)} sites, "
      f"{len(lmm_dose)} imputed in this dog")

# Precompute each site's dosage of its EFFECT allele once (shared site table,
# so this is trait-independent); None where the dog has no imputed genotype
# or the imputed alleles don't match the panel's.
_lmm_ea_dose = []
for _sn in LMM_SITES:
    _hit = lmm_dose.get((_sn[0], _sn[1]))
    if _hit and {_sn[2], _sn[3]} == {_hit[0], _hit[1]}:
        _lmm_ea_dose.append(_hit[2] if _sn[2] == _hit[1] else 2.0 - _hit[2])
    else:
        _lmm_ea_dose.append(None)

def compute_prs_lmm(tinfo):
    """Dense PRS from stored GEMMA betas; z and percentile against the stored
    reference distribution. Sites this dog cannot be imputed at contribute
    their panel-mean (beta * 2af), which keeps the dog's raw score on the
    same scale as ref_prs_sorted."""
    prs, matched = 0.0, 0
    for _sn, d, beta in zip(LMM_SITES, _lmm_ea_dose, tinfo['beta']):
        if beta == 0.0:
            continue
        if d is None:
            d = 2.0 * _sn[4]
        else:
            matched += 1
        prs += beta * d
    ref = tinfo['ref_prs_sorted']
    mu, sd = float(np.mean(ref)), float(np.std(ref)) + 1e-8
    z = (prs - mu) / sd
    pct = float(np.searchsorted(ref, prs) / len(ref) * 100)
    return z, pct, matched

# ── Compute per trait ─────────────────────────────────────────────────────
print("Computing PRS per trait (GEMMA LMM betas):")
traits_out = {}
for trait in TRAIT_COLS:
    tinfo = LMM['traits'].get(trait)
    if not tinfo: continue
    prs_z, pct, matched = compute_prs_lmm(tinfo)
    predicted = float(np.clip(tinfo['truth_mean'] + prs_z * tinfo['truth_sd'], 1, 5))
    traits_out[trait] = {
        'prs_z': round(float(prs_z), 3),
        'percentile': round(float(pct), 1),
        'predicted_score': round(predicted, 2),
        'n_ref_samples': int(LMM['meta']['n_dogs']),
        'n_snps': tinfo['n_snps'],
        'snps_matched': matched,
        'cv_r': tinfo['cv_r'],
        'description': f'LMM-PRS prediction for {trait.lower()}',
    }
    print(f"  {trait}: z={prs_z:.3f}, pct={pct:.1f}, pred={predicted:.2f} "
          f"({matched}/{tinfo['n_snps']} snps)")

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

# Height — dense GEMMA LMM betas (breed-blocked CV r 0.82 vs 0.71 for the ridge)
_th = LMM['traits'].get('height_cm')
z_h, pct_h, _mh = compute_prs_lmm(_th) if _th else (np.nan, np.nan, 0)
if not np.isnan(z_h):
    mu, sd = _th['truth_mean'], _th['truth_sd']
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

# Weight — dense GEMMA LMM betas (breed-blocked CV r 0.67 vs 0.58 for the
# ridge at breed level; dense chosen over sparse cuts for individual-dog
# robustness — see the artifact header above)
_tw = LMM['traits'].get('weight_kg')
z_w, pct_w, _mw = compute_prs_lmm(_tw) if _tw else (np.nan, np.nan, 0)
if not np.isnan(z_w):
    mu, sd = _tw['truth_mean'], _tw['truth_sd']
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

# ── Weight blend: dense breed-level pred + Darwin's Ark individual size PRS ──
# The two components err in OPPOSITE directions on chondrodysplastic breeds
# (dense +4.3 kg — FGF4 diluted across 131k SNPs; Darwin's Ark −4.3 kg), so a
# 2-coefficient blend cancels the bias. Validated 2026-08-18 on 94 cohort dogs
# with owner-reported weights: blend r=0.92 / MAE 5.1 kg vs 0.91/5.7 dense
# alone; dwarf-breed MAE 6.6 -> 3.8 kg. Coefficients fit on that cohort and
# stored with provenance in darwins_ark_blend.json; the raw Darwin's Ark PRS
# scale is platform-stable because unimputed sites fall back to beta*2af.
if 'weight_kg' in phys_traits:
    with open("$REF_JSON/darwins_ark_blend.json") as _f:
        _DAB = json.load(_f)
    _da_rows = []
    with gzip.open("$REF_JSON/darwins_ark_size_prs.tsv.gz", 'rt') as _f:
        for _line in _f:
            _da_rows.append(_line.split())
    _bedf2 = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
    for _r in _da_rows:
        _bedf2.write(f"{_r[0]}\t{int(_r[1])-1}\t{_r[1]}\n")
    _bedf2.close()
    _q2 = subprocess.run(['bcftools', 'query', '-R', _bedf2.name,
                          '-f', '%CHROM\t%POS\t%REF\t%ALT\t[%GP]\n', BCF],
                         capture_output=True, text=True)
    os.unlink(_bedf2.name)
    _da_dose = {}
    for _line in _q2.stdout.splitlines():
        _f2 = _line.split('\t')
        if len(_f2) != 5: continue
        _gp = _f2[4].split(',')
        if len(_gp) != 3: continue
        try:
            _g = [float(x) for x in _gp]
        except ValueError:
            continue
        _da_dose[(_f2[0], _f2[1])] = (_f2[2], _f2[3], _g[1] + 2.0*_g[2])
    _da_prs, _da_matched = 0.0, 0
    for _c, _p, _ea, _oa, _beta, _af in _da_rows:
        _beta, _af = float(_beta), float(_af)
        _hit = _da_dose.get((_c, _p))
        if _hit and {_ea, _oa} == {_hit[0], _hit[1]}:
            _d = _hit[2] if _ea == _hit[1] else 2.0 - _hit[2]
            _da_matched += 1
        else:
            _d = 2.0 * _af
        _da_prs += _beta * _d
    _dense_kg = phys_traits['weight_kg']['pred_kg']
    _blend_kg = float(np.clip(_DAB['intercept']
                              + _DAB['coef_dense_pred_kg'] * _dense_kg
                              + _DAB['coef_da_size_prs'] * _da_prs, 1, 120))
    phys_traits['weight_kg'].update({
        'pred_kg': round(_blend_kg, 1),
        'pred_lbs': round(_blend_kg * 2.205, 1),
        'pred_kg_dense': round(float(_dense_kg), 1),
        'da_size_prs': round(_da_prs, 2),
        'da_sites_matched': f'{_da_matched}/{len(_da_rows)}',
        'method_note': ('Blend of the breed-level dense LMM prediction with the '
                        "Darwin's Ark individual-level size PRS (Morrill 2022); "
                        'validated r=0.92, MAE 5.1 kg on 94 dogs with known weights.'),
    })
    print(f"  Weight blended: {_blend_kg:.1f}kg (dense {_dense_kg:.1f}kg, "
          f"DA prs {_da_prs:.1f}, {_da_matched}/{len(_da_rows)} sites)")

# ── Darwin's Ark individual-level traits (behaviour factors + physical) ────
# Trained on 2,155 dogs with INDIVIDUAL owner-reported phenotypes (Morrill
# 2022), so unlike the AKC breed-level scores these see within-breed
# variation. Weights are p<=0.1, platform-filtered to Dog10K-imputable sites
# (every site scores; no fallback dilution). Honest z uses the manifest's
# z_scale — in-sample references overstate spread ~25-40% (held-out CV;
# independently confirmed by 96 external dogs, sd ratio 0.795 vs 0.796).
# Behaviour factors additionally blend with the dog's breed-composition
# expectation (per-breed means from 14k survey-only purebreds) — held-out CV
# showed ancestry ~ PRS in signal with blend best on 8/9 traits.
import math as _math
with gzip.open("$REF_JSON/darwins_ark/manifest.json.gz", 'rt') as _f:
    _DAM = json.load(_f)
_da_all = {}      # trait -> list of (chr, pos, ea, oa, beta, af)
_da_union = {}    # (chr,pos) -> None (union of sites for ONE bcftools query)
for _t in _DAM['traits']:
    _rows = []
    with gzip.open(f"$REF_JSON/darwins_ark/wts_{_t}.tsv.gz", 'rt') as _f:
        for _line in _f:
            _p = _line.split()
            _rows.append((_p[0], _p[1], _p[2], _p[3], float(_p[4]), float(_p[5])))
            _da_union[(_p[0], _p[1])] = None
    _da_all[_t] = _rows
_bedf3 = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
for _c, _p in sorted(_da_union, key=lambda cp: (cp[0], int(cp[1]))):
    _bedf3.write(f"{_c}\t{int(_p)-1}\t{_p}\n")
_bedf3.close()
_q3 = subprocess.run(['bcftools', 'query', '-R', _bedf3.name,
                      '-f', '%CHROM\t%POS\t%REF\t%ALT\t[%GP]\n', BCF],
                     capture_output=True, text=True)
os.unlink(_bedf3.name)
_dam_dose = {}
for _line in _q3.stdout.splitlines():
    _f2 = _line.split('\t')
    if len(_f2) != 5: continue
    _gp = _f2[4].split(',')
    if len(_gp) != 3: continue
    try:
        _g = [float(x) for x in _gp]
    except ValueError:
        continue
    _dam_dose[(_f2[0], _f2[1])] = (_f2[2], _f2[3], _g[1] + 2.0*_g[2])
print(f"Darwin's Ark traits: {len(_DAM['traits'])} traits, "
      f"{len(_da_union)} union sites, {len(_dam_dose)} imputed")

# breed-composition expectation per trait (ancestry blend input)
_da_means = _DAM.get('breed_trait_means', {})
_da_comp = []
try:
    with open(f'{PUB}/breed_result.json') as _f:
        _brj = json.load(_f)
    _da_comp = [(_r['breed'], float(_r['proportion']))
                for _r in (_brj.get('breed_composition_raw') or _brj.get('breed_composition') or [])]
except Exception:
    pass
def _da_ancestry_expectation(trait):
    _tot, _s = 0.0, 0.0
    for _b, _pr in _da_comp:
        _v = _da_means.get(_b, {}).get(trait)
        if _v is not None:
            _tot += _pr; _s += _pr * _v
    return (_s / _tot) if _tot >= 0.5 else None

da_traits = {'factors': {}, 'physical': {}}
for _t, _info in _DAM['traits'].items():
    _prs, _matched = 0.0, 0
    for _c, _p, _ea, _oa, _beta, _af in _da_all[_t]:
        _hit = _dam_dose.get((_c, _p))
        if _hit and {_ea, _oa} == {_hit[0], _hit[1]}:
            _d = _hit[2] if _ea == _hit[1] else 2.0 - _hit[2]
            _matched += 1
        else:
            _d = 2.0 * _af
        _prs += _beta * _d
    _ref = _info['ref_prs_sorted']
    _rmu = float(np.mean(_ref)); _rsd = float(np.std(_ref)) + 1e-9
    _z = (_prs - _rmu) / (_rsd * _info['z_scale'])
    _pct = 50.0 * (1.0 + _math.erf(_z / _math.sqrt(2.0)))
    _entry = {'display': _info['display'],
              'prs_z': round(_z, 3),
              'percentile': round(_pct, 1),
              'h2_snp': _info['h2_snp'],
              'n_sites': _info['n_sites'],
              'sites_matched': _matched,
              'n_training_dogs': _info['n_pheno']}
    if 'cv_r_heldout' in _info:
        _entry['cv_r'] = _info['cv_r_heldout']
    if _info['kind'] == 'physical':
        _pred = _info['truth_mean'] + _z * _info['truth_sd']
        _anch = min(_info['anchors'], key=lambda a: abs(a[0] - _pred))
        _entry['predicted'] = _anch[1]
        da_traits['physical'][_t] = _entry
        print(f"  DA {_info['display']}: {_anch[1]} (z={_z:+.2f}, pct={_pct:.0f})")
    elif _info['kind'] == 'factor':
        _bl = _info.get('ancestry_blend')
        _anc = _da_ancestry_expectation(_t)
        if _bl and _anc is not None:
            _bv = _bl['w0'] + _bl['w_prs_z'] * _z + _bl['w_ancestry'] * _anc
            _zb = (_bv - _info['truth_mean']) / (_info['truth_sd'] + 1e-9)
            _entry['blended_percentile'] = round(50.0 * (1.0 + _math.erf(_zb / _math.sqrt(2.0))), 1)
            _entry['ancestry_expectation'] = round(_anc, 3)
            _entry['cv_r_blend'] = _bl['cv_r_blend']
        da_traits['factors'][_t] = _entry
        print(f"  DA {_info['display']}: z={_z:+.2f} pct={_pct:.0f}"
              + (f" blended_pct={_entry['blended_percentile']:.0f}" if 'blended_percentile' in _entry else ""))
    # 'size' kind feeds the weight blend above; skip separate reporting

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
        'heritability': {'h2': 0.32, 'ci': '0.18–0.46', 'source': 'Parker 2017 (Cell Reports)'},
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
        'heritability': {'h2': 0.63, 'ci': '0.52–0.74', 'source': 'Parker 2017 (Cell Reports)'},
    }
    print(f"  Coat length: {pred_cl} (ord={pred_ord:.2f}, z={z_cl:.3f})")

prs_result = {
    'traits': traits_out,
    'individual_traits': {
        **da_traits,
        'method': ("Individual-level PRS from Darwin's Ark (Morrill 2022, Science; 2,155 dogs "
                   "with owner-reported phenotypes; GEMMA LMM, p<=0.1 weights on "
                   "platform-imputable sites). Percentiles are held-out-calibrated (z_scale). "
                   "Behaviour factors blend PRS with the dog's breed-composition expectation."),
    },
    'physical_traits': phys_traits,
    'method': (f'GEMMA linear-mixed-model betas (GRM + genotype-source covariate), '
               f'trained on the merged Parker+Dog10K panel ({LMM["meta"]["n_dogs"]} dogs); '
               f'dense PRS = sum(beta x GLIMPSE2-imputed dosage) over all '
               f'{len(LMM_SITES)} panel SNPs (no p-value thresholding — sparse scores '
               f'are brittle for individual mixed-breed dogs). '
               f'Coat type and length remain marginal-ridge. '
               f'{len(lmm_dose)}/{len(LMM_SITES)} LMM sites imputed in this dog.'),
    'reference': 'AKC breed trait scores (kkakey/dog_traits_AKC) + breed-standard height/weight',
    'snps': n_snps,
    'snps_imputed': genotyped,
    'n_ref_breeds': int(len(np.unique(breeds_fam))),
}
with open(f'{PUB}/prs_result.json', 'w') as f:
    json.dump(prs_result, f, indent=2)
print(f"prs_result.json written ({len(traits_out)} traits)")
PYEOF
# ── 11c: Orthopedic risk (breed-ancestry-weighted epidemiology) ──────────
# Hip dysplasia and cranial cruciate ligament rupture are polygenic — no
# single variant to genotype — but both have solid breed-level epidemiology
# (OFA hip statistics; VetCompass/JAVMA CCL odds), and this dog's ancestry
# composition is known precisely. The report therefore gives an
# ancestry-weighted risk, clearly framed as breed statistics, not an
# individual DNA test.
log "  Orthopedic risk (ancestry-weighted)"
PUB_DIR="$PUB" REFJ="$REF_JSON" "$DATA_PYTHON" - <<'PYEOF'
import json, os, re

pub = os.environ['PUB_DIR']
art = json.load(open(os.environ['REFJ'] + '/ortho_risk.json'))
breed = json.load(open(pub + '/breed_result.json'))
bc = breed['breed_composition']
comp = (sorted(bc.items(), key=lambda x: -x[1]) if isinstance(bc, dict)
        else [(e.get('breed_name') or e.get('breed'), e['proportion']) for e in bc])

def norm(x):
    return re.sub(r'[^a-z0-9]+', '', str(x).replace('_', ' ').lower())

def weighted(table, default):
    total_w = 0.0; acc = 0.0; contrib = []
    for name, prop in comp:
        v = table.get(norm(name))
        if v is None:
            # panel size-variant labels: try dropping leading Standard/Miniature/Toy
            v = table.get(norm(re.sub(r'^(standard|miniature|toy)\s+', '', str(name), flags=re.I)))
        if v is None:
            v = default
        else:
            contrib.append({'breed': str(name).replace('_', ' ').title(), 'value': v, 'proportion': round(prop, 3)})
        acc += prop * v; total_w += prop
    return (acc / total_w if total_w else default), contrib

hip_tab = art['hip_dysplasia']['by_breed']
hip_mean = art['hip_dysplasia']['population_mean_pct']
hip, hip_contrib = weighted(hip_tab, hip_mean)
hip_ratio = hip / hip_mean
hip_cat = ('elevated' if hip_ratio >= 1.5 else 'somewhat elevated' if hip_ratio >= 1.15
           else 'reduced' if hip_ratio <= 0.7 else 'average')

ccl_tab = art['ccl']['by_breed']
ccl, ccl_contrib = weighted(ccl_tab, 1.0)
ccl_cat = ('elevated' if ccl >= 1.5 else 'somewhat elevated' if ccl >= 1.15
           else 'reduced' if ccl <= 0.7 else 'average')

out = {
  'hip_dysplasia': {
    'ancestry_weighted_pct': round(hip, 1),
    'population_average_pct': hip_mean,
    'relative_risk': round(hip_ratio, 2),
    'category': hip_cat,
    'top_contributors': sorted(hip_contrib, key=lambda c: -c['proportion'] * c['value'])[:4],
    'source': art['hip_dysplasia']['source'],
    'caveat': art['hip_dysplasia']['caveat'],
  },
  'ccl': {
    'relative_risk': round(ccl, 2),
    'category': ccl_cat,
    'top_contributors': sorted(ccl_contrib, key=lambda c: -c['proportion'] * abs(c['value'] - 1))[:4],
    'sources': art['ccl']['sources'],
  },
  'prevention_note': ('For both conditions, keeping a lean body weight is the single most '
                      'effective ownable factor; controlled exercise, and for at-risk dogs a '
                      'discussion with your vet about conditioning and neuter timing, also help.'),
  'method': art['method'],
}
with open(pub + '/ortho_result.json', 'w') as f:
    json.dump(out, f, indent=2)
print('ortho_result.json: hip {:.1f}% ({}), CCL RR {:.2f} ({})'.format(hip, hip_cat, ccl, ccl_cat))
PYEOF

fi # end stage 11

if (( FROM_STAGE <= 12 && TO_STAGE >= 12 )); then
# ── Stage 12: Inbreeding (Dog10K ROH + F distribution) ──────
log "=== Stage 12: Inbreeding (Dog10K) ==="
"$DATA_PYTHON" - << PYEOF
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
fi # end stage 12

if (( FROM_STAGE <= 13 && TO_STAGE >= 13 )); then
# ── Stage 13: Coat color (GLIMPSE2 imputed genotypes at causal loci) ─────
log "=== Stage 13: Coat color ==="
# Same BAM resolution as stage 8: fall back to the permanent sites.bam, refuse
# to run with neither (merle and B-locus phasing need reads).
CALL_BAM="$OUT/markdup.bam"
[[ -f "$CALL_BAM" ]] || CALL_BAM="$OUT/sites.bam"
[[ -f "$CALL_BAM" ]] || die "Stage 13 needs reads: neither markdup.bam nor sites.bam in $OUT — refusing to write a degraded report"
export IMPUTED_BCF MARKDUP_BAM="$CALL_BAM" PUB DOG_LOWER
"$DATA_PYTHON" - << 'PYEOF'
import subprocess, json, pysam, re, tempfile, os

BCF     = os.environ['IMPUTED_BCF']
BAM     = os.environ['MARKDUP_BAM']
PUB     = os.environ['PUB']
MIN_GP  = 0.80   # min max(GP) to trust a GLIMPSE2 call

# ── Causal variant table ──────────────────────────────────────────────────
# exp_ref/exp_alt: expected REF/ALT in canFam4; if swapped in BCF, n_alt is flipped.
# inheritance: 'recessive' (need 2 copies) | 'dominant' (1 copy sufficient)
KNOWN_VARIANTS = [
    # E locus: MC1R (chr5) — MC1R's CDS sits at chr5:64,186,690–64,187,643 on the MINUS
    # strand in canFam4 (verified by translating the reverse complement: the frame opens
    # with MSGQGPQRRLLGSLN and codon 306 is CGA = Arg exactly at 64,186,728–). Earlier
    # revisions placed the gene at chr5:63.92Mb (another assembly's coordinate) and called
    # 64186854 a "CPNE7 intron proxy" — both wrong.
    #
    # e1 — the causal recessive-red stop, MC1R c.916 p.Arg306* (codon CGA→TGA). Plus
    # strand: ref G, alt A. NOT in the Dog10K panel, so only direct reads can show it;
    # min_reads=1 because even a single read at the stop codon is informative context.
    dict(locus='E', chrom='chr5', pos=64186728, exp_ref='G', exp_alt='A',
         allele='e', inheritance='recessive', min_reads=1,
         effect='MC1R p.Arg306* — causal recessive-red (e) stop, c.916C>T coding = chr5:64186728 G>A plus strand'),
    # e_tag — MC1R c.790 M264V (codon 264). The canFam4 REFERENCE (a German Shepherd, a
    # masked breed) carries V264 = the Em melanistic-mask allele; the ALT (plus-strand T,
    # coding A) is the ordinary non-mask M264 that the e stop arose on. So alt copies
    # bound the possible e copies: alt/alt is NECESSARY for e/e but far from sufficient
    # (alt AF ≈71% in Dog10K), and ref copies are mask-type functional MC1R.
    dict(locus='E', chrom='chr5', pos=64186854, exp_ref='C', exp_alt='T',
         allele='e_tag', inheritance='recessive',
         effect='MC1R c.790 M264V (chr5:64186854) — REF C = V264 Em mask haplotype; ALT T = non-mask M264 background that e arises on'),

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
    # Second brown allele: the bc-region haplotype at the exon-1 end of TYRP1
    # (minus strand, so codon 41 sits at the HIGH coordinate end ~33.40Mb).
    # Chosen empirically from the Dog10K panel: carried by the b1-heterozygous
    # obligate-brown dogs (Vizslas, Weimaraners) and by zero Boxers.
    #
    # The previous entry here (chr11:33440938, labeled "b2 p.Gln354*") was
    # WRONG: its ALT allele runs at 53% panel-wide, is carried by Boxers and
    # absent from obligate browns — a common linked variant, not a brown
    # allele. Combined with phase-blind compound-het logic it called a third
    # of the reference cohort chocolate (b/b), including at least one dog
    # whose owner confirms a black nose.
    dict(locus='B', chrom='chr11', pos=33400944, exp_ref='A', exp_alt='G',
         allele='b', inheritance='recessive', effect='TYRP1 bc-region brown haplotype tag'),

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
          'e':  'Recessive red/yellow — two copies needed',
          'e?': 'Possibly recessive red — e-compatible background but causal stop not covered by reads'},
    'K': {'KB':  'Dominant black — one copy = solid black',
          'kbr': 'Brindle (incompletely dominant)',
          'ky':  'Non-black/agouti — A locus determines pattern'},
    'A': {'ay': 'Sable/fawn (dominant, regulatory variant)',
          'aw': 'Wild type agouti', 'at': 'Tan points / tricolor (recessive)',
          'a':  'Recessive black (recessive)'},
    'B': {'B': 'Black eumelanin (dominant)', 'b': 'Brown/liver eumelanin — two copies needed'},
    'D': {'D': 'Full pigment (dominant)', 'd': 'Dilute/blue — two copies needed'},
    'M': {'M': 'Merle (dominant, PMEL SINE insertion — detected from reads at the insertion site)', 'm': 'Non-merle'},
    'S': {'S': 'Solid / minimal white', 'sp': 'Piebald spotting (recessive)', 'sw': 'Extreme white (recessive)'},
    'W': {'w': 'Non-white', 'W': 'Extreme white (dominant, KIT structural — not detectable from SNP data)'},
}

LOCUS_INFO = {
    'E': dict(gene='MC1R',  chrom='chr5',  name='Extension locus',
              role='Master pigment switch: eumelanin (black/brown) vs phaeomelanin (yellow/red)',
              phenotype_contribution='e/e → all coat pigment yellow/red; called from reads at the causal MC1R stop plus the M264V haplotype; low-pass depth can leave it unresolved'),
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
              phenotype_contribution='Detected from sequencing reads at the PMEL insertion site; the merle class (cryptic to harlequin) needs a specialized length test'),
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
def bam_pileup(chrom, pos, min_reads=5):
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
        return counts if sum(counts.values()) >= min_reads else None
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
        counts = bam_pileup(v['chrom'], v['pos'], v.get('min_reads', 5))
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

# ── Merle (PMEL SINE insertion) from read alignments ─────────────────────
# The merle allele is a ~250-280bp SINE inserted at the intron10/exon11
# boundary of PMEL. Coordinate derived by mapping the RefSeq (ROS_Cfam
# NC_051814.1:299657, NM_001103216.1 exon structure) junction flank onto
# canFam4: unique 800bp MAPQ60 hit -> junction at chr10:644,511 (15bp TSD).
# Short reads cannot span the insertion, so a merle chromosome shows up as
# reads SOFT-CLIPPED at the junction; a non-merle chromosome as reads
# spanning it cleanly. The poly-A length that separates cryptic from full
# merle is NOT measurable from short reads — a positive call is therefore
# "a merle-type insertion is present", class unknown.
MERLE_CHROM, MERLE_POS, MERLE_TSD = 'chr10', 644511, 15

def detect_merle_sine():
    try:
        bam_fh = pysam.AlignmentFile(BAM, 'rb')
    except Exception:
        return None
    n_clip = 0; n_span = 0
    lo, hi = MERLE_POS - MERLE_TSD, MERLE_POS + MERLE_TSD
    try:
        for r in bam_fh.fetch(MERLE_CHROM, MERLE_POS - 200, MERLE_POS + 200):
            if r.is_unmapped or r.mapping_quality < 20 or r.is_secondary or r.is_supplementary:
                continue
            start, end = r.reference_start + 1, r.reference_end  # 1-based inclusive
            cig = r.cigartuples or []
            clip_left  = cig and cig[0][0] in (4, 5)
            clip_right = cig and cig[-1][0] in (4, 5)
            # a clip whose breakpoint falls inside the TSD window is insertion evidence
            if (clip_right and lo <= end <= hi and cig[-1][1] >= 8) or                (clip_left and lo <= start <= hi and cig[0][1] >= 8):
                n_clip += 1
            elif start <= lo - 5 and end >= hi + 5:
                n_span += 1
    except Exception:
        return None
    finally:
        bam_fh.close()
    return {'n_clip': n_clip, 'n_span': n_span}

MERLE_EVIDENCE = detect_merle_sine()

def call_locus(locus, calls):
    """Returns (allele1, allele2, confidence, interpretation)."""

    if locus == 'E':
        # Two informative sites, both inside MC1R (see the variant table):
        #   e1  chr5:64186728 G>A  — the causal p.Arg306* stop; BAM reads only.
        #   tag chr5:64186854 C>T  — c.790 M264V; REF = Em mask-type functional
        #       haplotype, ALT = the non-mask background e arises on, so the
        #       imputed alt count is an upper bound on possible e copies.
        e1  = next((c for c in calls if c['locus'] == 'E' and c.get('pos') == 64186728
                    and c.get('bam_counts')), None)
        tag = next((c for c in calls if c['locus'] == 'E' and c.get('pos') == 64186854
                    and c['found'] and c.get('n_alt') is not None), None)
        e1_ref = e1['bam_counts'].get('G', 0) if e1 else 0
        e1_alt = e1['bam_counts'].get('A', 0) if e1 else 0
        if not tag and not e1:
            return '?', '?', 'low', 'MC1R not covered by the imputation panel or by reads at this depth'
        # Direct reads at the stop codon outrank the tag.
        if e1_alt >= 2 and e1_ref == 0:
            return ('e', 'e', 'medium' if e1_alt >= 3 else 'low',
                f'Recessive red (e/e): {e1_alt} reads carry the MC1R p.Arg306* stop and none the functional '
                f'allele — all coat pigment is phaeomelanin (cream/yellow/red)')
        if e1_ref >= 1 and e1_alt >= 1:
            return ('E', 'e', 'low',
                f'Carrier for recessive red (E/e): reads at MC1R p.306 show both the functional allele '
                f'({e1_ref}) and the p.Arg306* stop ({e1_alt}) — eumelanin coat is expressed')
        if e1_ref >= 2:
            return ('E', 'E' if e1_ref >= 3 else '?',
                'medium' if e1_ref >= 3 else 'low',
                f'No recessive-red stop in {e1_ref} reads at MC1R p.306 — a functional copy is present, '
                f'eumelanin coat is expressed')
        # Stop codon not covered decisively — bound e by the M264V tag.
        n_tag = tag['n_alt'] if tag else None
        if n_tag == 0:
            return ('Em', 'Em', 'medium',
                'Both MC1R copies are the V264 mask-type functional haplotype (Em/Em) — eumelanin coat '
                'expressed; a melanistic mask may show where the coat is phaeomelanin-patterned')
        if n_tag == 1:
            return ('E', 'Em', 'medium',
                'One mask-type (V264) and one standard MC1R copy — at most one recessive-red allele '
                'possible, so a functional copy is present and eumelanin coat is expressed')
        if n_tag == 2:
            hint = (f' One read at the stop codon does carry e — consistent with E/e or e/e.'
                    if e1_alt == 1 else '')
            return ('e?', 'e?', 'low',
                'Both MC1R copies are the non-mask haplotype that recessive red arises on, but the causal '
                'p.Arg306* stop itself is not covered by enough reads to decide. The dog may be e/e '
                '(cream/yellow/red coat) or E/E / E/e (eumelanin coat expressed).' + hint +
                ' The dog\'s actual coat color resolves this.')
        return '?', '?', 'low', 'E locus not resolvable at this sequencing depth'

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
        b_calls = [c for c in calls
                   if c['locus'] == 'B' and c['allele'] == 'b'
                   and c['found'] and c['n_alt'] is not None]
        if not b_calls:
            return '?', '?', 'low', 'TYRP1 brown alleles not found in Dog10K panel'
        if any(c['n_alt'] == 2 for c in b_calls):
            return 'b', 'b', 'high', 'Brown/liver eumelanin (b/b) — black pigment becomes brown; nose and pads liver/brown'
        het = [c for c in b_calls if c['n_alt'] == 1]
        if len(het) >= 2:
            # Two different brown variants, each heterozygous. Phase decides:
            # on OPPOSITE haplotypes (trans) both gene copies are broken and
            # the dog is brown; on the SAME haplotype (cis) one copy is intact
            # and the dog is a black-nosed carrier. GLIMPSE genotypes are
            # phased, so read the haplotype side of each ALT.
            sides = set()
            phased = True
            for c in het:
                gt = str(c.get('gt', ''))
                if '|' not in gt:
                    phased = False
                    break
                sides.add(gt.split('|').index('1'))
            if phased and len(sides) > 1:
                return 'b', 'b', 'medium', ('Brown/liver eumelanin — two different brown variants on opposite '
                    'chromosome copies (compound heterozygous). Phasing is statistical at low pass, so treat '
                    'with moderate confidence.')
            if phased:
                return 'B', 'b', 'medium', ('Carrier (B/b): two brown variants detected but on the SAME '
                    'chromosome copy, so one intact copy remains — black pigment, black nose. '
                    'Puppies may inherit the brown haplotype.')
            return 'b', '?', 'low', ('Two brown variants detected but phase is unavailable; genotype is '
                'B/b (carrier, black nose) or b/b (brown) — a DNA test with parental phasing would resolve it.')
        if len(het) == 1:
            return 'B', 'b', 'medium', ('Carrier (B/b): one brown allele detected. Pigment stays black '
                '(black nose and pads); puppies may inherit brown if the other parent also carries it. '
                'Rare brown alleles outside the panel cannot be fully excluded.')
        return 'B', 'B', 'high', 'No brown alleles detected (B/B) — black eumelanin, black nose and pads'

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
        ev = MERLE_EVIDENCE
        if not ev:
            return 'm', 'm', 'low', 'Merle (PMEL SINE insertion) not assessable — no read data at the PMEL junction'
        nc, ns = ev['n_clip'], ev['n_span']
        if nc >= 2:
            if ns == 0 and nc >= 3:
                return 'M*', 'M*', 'low', ('Merle-type SINE insertion detected on both read sets ({} clipped, 0 spanning reads) — '
                    'likely two merle-family alleles. The class (cryptic to full merle) requires a specialized '
                    'length test.').format(nc)
            return 'M*', 'm', 'medium' if nc >= 3 else 'low', ('Merle-type SINE insertion detected ({} clipped vs {} clean reads at the '
                'PMEL junction). The dog carries a merle-family allele; whether it shows as merle depends on the '
                'insertion length (cryptic merle looks solid), which requires a specialized test.').format(nc, ns)
        if nc == 1:
            return 'm', '?', 'low', ('One read hints at a merle-type insertion ({} clean reads) — inconclusive at this '
                'sequencing depth; a targeted merle test would resolve it.').format(ns)
        if ns >= 5:
            return 'm', 'm', 'medium', 'No merle insertion seen across {} reads spanning the PMEL junction (m/m).'.format(ns)
        if ns >= 2:
            return 'm', 'm', 'low', ('No merle insertion seen, but only {} reads span the PMEL junction — a merle allele '
                'could be missed at this depth.').format(ns)
        return 'm', '?', 'low', 'Too few reads at the PMEL junction to assess merle at this sequencing depth.'

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
# When the E locus is unresolved ('e?' — both copies on the non-mask M264V
# background but the causal p.Arg306* stop not covered by reads), use K and B
# to name what the coat would be if a functional copy is present, so the report
# states both scenarios concretely instead of asserting cream/red.
e_gt = loci_gt['E']
if e_gt['allele1'] == 'e?' and e_gt['allele2'] == 'e?':
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
        f' If a functional MC1R copy is present, the coat would be {alt_color} '
        f'based on K ({k_gt["allele1"]}/{k_gt["allele2"]}) and '
        f'B ({b_gt["allele1"]}/{b_gt["allele2"]}) loci.'
    )
    loci_gt['E'] = {**e_gt, 'interpretation': updated_interp,
                    'eumelanic_alternative': alt_color}

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

    # Unresolved E ('e?'): both MC1R copies on the e-compatible background but
    # the causal stop not covered by reads. State both scenarios rather than
    # asserting either one — the owner's eyes resolve it instantly.
    if e1 == 'e?':
        alt_solid = f'solid {eume_color}' if 'KB' in (k1, k2) else f'{eume_color} (patterned by the A locus)'
        base_color = (f'Not resolved at this sequencing depth — either cream/yellow/red (if e/e) '
                      f'or {alt_solid} (if a functional MC1R copy is present). '
                      f'The dog\'s actual coat color tells you which; {nose_color} from B/D loci either way')
        pattern = (f'If e/e: solid phaeomelanin (cream/yellow/red). '
                   f'Otherwise: {alt_solid}. Skin pigment ({nose_color}) comes from B and D loci independently.')
        dilution = 'See B/D loci — applies to the eumelanin scenario'
        return base_color, pattern, dilution

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
    if g.get('eumelanic_alternative'):
        locus_entry['eumelanic_alternative'] = g['eumelanic_alternative']
    loci_result[locus] = locus_entry

coat = {
    'summary': {
        'predicted_base_color': base_color,
        'predicted_pattern': pattern,
        'predicted_dilution': dilution,
        'predicted_white': 'Not detectable from SNP data (S locus limited; W requires structural variant)',
        'predicted_merle': loci_gt['M']['interpretation'],
        'overall_confidence': overall_conf,
        **({'validation_warning': validation_warning} if validation_warning else {}),
        'caveat': ('E, K, B, D loci called from Dog10K GLIMPSE2 imputed BCF (causal SNPs). '
                   'A locus sable (ay/aw) requires structural variant analysis not available here. '
                   'Merle (M) is screened from reads at the PMEL insertion site (the exact merle class needs a length test); extreme white (W) requires PCR or long-read. '
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

fi # end stage 13

# Stages 15-17 had no FROM_STAGE guard at all, so they ran on every resume
# regardless of what was asked. TO_STAGE could not stop them either.
if (( FROM_STAGE <= 15 && TO_STAGE >= 15 )); then
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

    METAPHLAN_LOG="$OUT/${DOG_LOWER}_metaphlan_stderr.log"
    # A site may point at an existing database (a shared cluster copy, or one
    # downloaded separately) instead of MetaPhlAn's default per-install location.
    # The database is ~34GB extracted, so this avoids fetching it per environment.
    MPA_DB_ARG=()
    if [[ -n "${METAPHLAN_DB:-}" ]]; then
        [[ -d "$METAPHLAN_DB" ]] || die "METAPHLAN_DB set but not a directory: $METAPHLAN_DB"
        # MetaPhlAn 4.2 renamed --bowtie2db to --db_dir; the old name is rejected.
        MPA_DB_ARG=(--db_dir "$METAPHLAN_DB")
        # Pin the index too, or MetaPhlAn ignores db_dir and fetches "latest".
        [[ -n "${METAPHLAN_INDEX:-}" ]] && MPA_DB_ARG+=(--index "$METAPHLAN_INDEX")
        log "  Using MetaPhlAn database: $METAPHLAN_DB"
    fi
    # Ensure MetaPhlAn4's Python bin is in PATH so its helper scripts (read_fastx.py) can be found
    export PATH="$(dirname "$METAPHLAN_BIN"):$PATH"
    log "  Running MetaPhlAn4…"
    if [[ -f "$MICRO_BT2" ]]; then
        if bzip2 -t "$MICRO_BT2" 2>/dev/null; then
            log "  Reusing existing mapout: $MICRO_BT2"
            "$METAPHLAN_BIN" "$MICRO_BT2" ${MPA_DB_ARG[@]+"${MPA_DB_ARG[@]}"} \
                --input_type mapout \
                --nproc 4 \
                -o "$MICRO_OUT" 2>"$METAPHLAN_LOG" \
            || { log "  MetaPhlAn4 error:"; cat "$METAPHLAN_LOG" | head -20 | while read -r l; do log "    $l"; done; die "MetaPhlAn4 failed"; }
        else
            log "  Mapout truncated — removing and re-running from FASTQ"
            rm -f "$MICRO_BT2"
            "$METAPHLAN_BIN" "$UNMAPPED_FQ" ${MPA_DB_ARG[@]+"${MPA_DB_ARG[@]}"} \
                --input_type fastq \
                --mapout "$MICRO_BT2" \
                --nproc 4 \
                -o "$MICRO_OUT" 2>"$METAPHLAN_LOG" \
            || { log "  MetaPhlAn4 error:"; cat "$METAPHLAN_LOG" | head -20 | while read -r l; do log "    $l"; done; die "MetaPhlAn4 failed"; }
        fi
    else
        "$METAPHLAN_BIN" "$UNMAPPED_FQ" ${MPA_DB_ARG[@]+"${MPA_DB_ARG[@]}"} \
            --input_type fastq \
            --mapout "$MICRO_BT2" \
            --nproc 4 \
            -o "$MICRO_OUT" 2>"$METAPHLAN_LOG" \
        || { log "  MetaPhlAn4 error:"; cat "$METAPHLAN_LOG" | head -20 | while read -r l; do log "    $l"; done; die "MetaPhlAn4 failed"; }
    fi

    log "  Computing microbiome JSONs…"
    # Pin to an interpreter that has the data-science stack. Bare `python3` here
    # resolves to the genomics env (prepended to PATH by Stage 10), which lacks
    # pandas — so a fresh run reaches this point without it while a resume-from-15
    # does not. Probe explicit candidates for the full stack.
    _site_python="${DATA_PYTHON:-}"
    DATA_PYTHON=""
    for cand in "$_site_python" /usr/bin/python3 "$(command -v python3 || true)"; do
        if [ -n "$cand" ] && [ -x "$cand" ] && "$cand" -c 'import pandas,numpy,scipy,sklearn' 2>/dev/null; then
            DATA_PYTHON="$cand"; break
        fi
    done
    [ -n "$DATA_PYTHON" ] || die "No python3 with pandas/numpy/scipy/sklearn found for microbiome analysis"
    log "  Using python: $DATA_PYTHON"
    "$DATA_PYTHON" - << PYEOF 2>&1 | while IFS= read -r l; do log "  [py] $l"; done; [ "${PIPESTATUS[0]}" -eq 0 ] || die "Microbiome Python block failed"
import json, re, math, datetime
import numpy as np
from scipy.stats import percentileofscore, entropy
from sklearn.linear_model import RidgeCV, ElasticNetCV

PUB      = "$PUB"
OUT      = "$OUT"
DOG      = "$DOG_NAME"
MICRO_OUT = "$MICRO_OUT"
REF_PANEL = "$MICROBIOME_REF"
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

# ── 2. Reference panel: the 96 cohort dogs, processed by THIS pipeline ─────
# Replaces the 1,045-sample MetaPhlAn 3.18 CSV. That reference was larger but
# taxonomically incompatible: the version gap meant the age model matched only
# ~16 of its 62 features for a typical sample, and diversity percentiles were
# meaningless. Here every feature matches by construction — same MetaPhlAn,
# same database, same read handling.
with open(REF_PANEL) as fh:
    _panel = json.load(fh)
panel_dogs = _panel['dogs']
sp_filtered = sorted({c for d in panel_dogs for c in d['species']})
print(f"Panel: {len(panel_dogs)} dogs, {len(sp_filtered)} species features "
      f"(db {_panel['meta'].get('db_version','?')})")

# ── 3. This dog's species dict (% of classified bacteria) ──
kiki_species = {}
for sp in taxa['species']:
    kiki_species[sp['clade']] = round(sp['relative_abundance'] / scale, 6)
matched_features = [f for f in sp_filtered if f in kiki_species]
print(f"Matched features: {len(matched_features)}")

# ── 4. Age prediction (RidgeCV on the cohort) ──────────────
# Panel v2 contains this cohort's own dogs: never train on the sample being
# scored (self-inclusion leaks the label and flatters the prediction).
aged = [d for d in panel_dogs
        if d.get('age') is not None
        and len([v for v in d['species'].values() if v]) >= 5
        and str(d.get('sample', '')).lower() != DOG_LOWER.lower()]
# Two conditions under which the prediction must not be made.
#
# Too few aged dogs: the model is fiction.
#
# Platform mismatch: validated directly on this cohort. Four MGI (100bp)
# samples scored against the Illumina (151bp) panel all predicted ~+5 years
# over baseline, with near-identical per-feature contributions for a
# 2.5-year-old and a 9-year-old — the model reads the batch shift as age, and
# no per-sample statistic we tried (range, centroid or nearest-neighbour
# distance) can detect the shift. Read length is the fingerprint we have.
# The UI shows no age card when the file carries null.
MIN_AGED = 20
_panel_rl = set(_panel['meta'].get('read_lengths_bp') or [])
_sample_rl = None
try:
    with open(f'{PUB}/qc_result.json') as _qf:
        _sample_rl = json.load(_qf).get('read_length_bp')
except Exception:
    pass
_platform_ok = not _panel_rl or _sample_rl is None or int(_sample_rl) in _panel_rl
skip_reason = None
_sample_n_sp = len([v for v in kiki_species.values() if v])
if len(aged) < MIN_AGED:
    skip_reason = f"only {len(aged)} panel dogs have ages (<{MIN_AGED})"
elif _sample_n_sp < 5:
    skip_reason = f"only {_sample_n_sp} species detected in this sample (<5) — profile too sparse to score"
elif not _platform_ok:
    skip_reason = (f"read length {_sample_rl}bp vs panel {sorted(_panel_rl)} — "
                   f"cross-platform age prediction reads batch shift as age")
if skip_reason:
    print(f"WARNING: age prediction skipped — {skip_reason}")
    with open(f'{PUB}/microbiome_age_result.json', 'w') as fh:
        fh.write('null')
else:
    X_ref = np.log10(np.array([[d['species'].get(f, 0.0) for f in sp_filtered]
                               for d in aged]) + 1e-5)
    y_ref = np.array([d['age'] for d in aged])

    # ElasticNet for the pooled 600+-dog panel (validated 2026-09-01: CV MAE
    # 2.12y r 0.67 vs Ridge 2.30y r 0.61); RidgeCV kept for small panels
    # where ElasticNet's sparsity is unstable.
    if len(aged) >= 100:
        try:
            model = ElasticNetCV(l1_ratio=[.1, .5, .9], n_alphas=20, max_iter=5000, cv=5)
        except TypeError:  # sklearn >=1.9 renamed n_alphas -> alphas=<int>
            model = ElasticNetCV(l1_ratio=[.1, .5, .9], alphas=20, max_iter=5000, cv=5)
        model_name = 'ElasticNetCV (log10, prevalence>10%, pooled panel)'
    else:
        model = RidgeCV(alphas=[0.01,0.1,1,10,100], cv=5)
        model_name = 'RidgeCV (log10, prevalence>10%, cohort panel)'
    model.fit(X_ref, y_ref)

    from sklearn.model_selection import cross_val_score
    cv_r2  = cross_val_score(model, X_ref, y_ref, cv=5, scoring='r2').mean()
    cv_mae = -cross_val_score(model, X_ref, y_ref, cv=5, scoring='neg_mean_absolute_error').mean()

    kiki_vec = np.array([kiki_species.get(f, 0.0) for f in sp_filtered])
    kiki_vec_log = np.log10(kiki_vec + 1e-5).reshape(1,-1)
    pred_age = float(model.predict(kiki_vec_log)[0])

    coef_pairs = sorted(zip(sp_filtered, model.coef_), key=lambda x: -abs(x[1]))[:10]
    top_species = [{'name': re.sub(r'.*\|s__', '', f).replace('_', ' '), 'coefficient': round(c, 5)}
                   for f, c in coef_pairs]

    age_result = {
        'predicted_age_years': round(pred_age, 2),
        'cv_r2':               round(cv_r2,  3),
        'cv_mae_years':        round(cv_mae, 3),
        'n_training_samples':  len(aged),
        'n_species_features':  len(sp_filtered),
        'n_features_matched': len(matched_features),
        'model':               model_name,
        'reference':           f"{len(aged)} reference dogs with known ages, same pipeline and database",
        'top_species':         top_species,
    }
    if ACTUAL_AGE:
        try:
            age_result['actual_age_years'] = float(ACTUAL_AGE)
        except ValueError:
            pass

    with open(f'{PUB}/microbiome_age_result.json', 'w') as fh:
        json.dump(age_result, fh, indent=2)
    print(f"microbiome_age_result.json: predicted={pred_age:.1f} yrs "
          f"(cv_r2={cv_r2:.2f}, mae={cv_mae:.2f} on {len(aged)} dogs)")

# ── 5. Pathobionts and diversity vs the cohort ─────────────
PATHOBIONTS = {
    's__Porphyromonas_gulae':          ('red',    'canine periodontal disease'),
    's__Tannerella_forsythia':         ('red',    'periodontal disease'),
    's__Porphyromonas_cangingivalis':  ('red',    'canine periodontitis'),
    's__Porphyromonas_canoris':        ('orange', 'canine oral disease'),
    's__Porphyromonas_gingivicanis':   ('red',    'canine periodontitis'),
    's__Treponema_denticola':          ('red',    'periodontal disease'),
    's__Fusobacterium_nucleatum':      ('orange', 'periodontal disease'),
    's__Prevotella_intermedia':        ('orange', 'periodontal disease'),
}

sp_pct_dict = {sp['clade']: sp['relative_abundance'] / scale
               for sp in taxa['species']}

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

ref_path_vec = np.array([d['pathobiont_pct'] for d in panel_dogs])
path_pct = round(percentileofscore(ref_path_vec, pathobiont_total, kind='rank'), 1)
r_pm, r_pmed, r_p75p, r_p90p = (round(float(x),2) for x in
    [ref_path_vec.mean(), np.median(ref_path_vec),
     np.percentile(ref_path_vec,75), np.percentile(ref_path_vec,90)])

# Diversity vs the cohort. Directly comparable now — same database — so the
# genus-level cross-version matching machinery is gone. The dashboard no longer
# displays these, but they are kept in the JSON for the record.
sample_richness = len(taxa['species'])
_ab = np.array([s['relative_abundance'] for s in taxa['species']])
sample_shannon = round(float(entropy(_ab / _ab.sum())), 4) if len(_ab) else 0.0
ref_rich = np.array([d['richness'] for d in panel_dogs])
ref_shan = np.array([d['shannon'] for d in panel_dogs])

health_result = {
    'sample_richness':        sample_richness,
    'sample_shannon':         sample_shannon,
    'richness_percentile':    round(percentileofscore(ref_rich, sample_richness, kind='rank'), 1),
    'shannon_percentile':     round(percentileofscore(ref_shan, sample_shannon, kind='rank'), 1),
    'ref_richness_p50':       int(np.median(ref_rich)),
    'ref_shannon_p50':        round(float(np.median(ref_shan)), 4),
    'reference':              f"{len(panel_dogs)} cohort dogs, same pipeline and database",
    'pathobiont_burden_pct':  round(pathobiont_total, 3),
    'pathobiont_percentile':  path_pct,
    'commensal_pct':          round(commensal_pct, 3),
    'ref_pathobiont_mean':    r_pm,
    'ref_pathobiont_median':  r_pmed,
    'ref_pathobiont_p75':     r_p75p,
    'ref_pathobiont_p90':     r_p90p,
    'pathobiont_hits':        pathobiont_hits,
}
with open(f'{PUB}/microbiome_health_result.json', 'w') as fh:
    json.dump(health_result, fh, indent=2)
print(f"microbiome_health_result.json: pathobionts={pathobiont_total:.1f}% "
      f"({path_pct}th pct of cohort), richness={sample_richness}")
PYEOF

    log "  Microbiome stage complete."

fi # end stage 15
if (( FROM_STAGE <= 16 && TO_STAGE >= 16 )); then
# ── Stage 16: Copy reference JSONs ───────────────────────────
log "=== Stage 16: Copy reference JSONs ==="
# NB: cnv_genes.json is NOT copied here — Stage 10 rebuilds it genome-wide for
# this sample, and copying cosmo's static version over it replaced each dog's
# own CNV gene calls with cosmo's.
for f in centromeres.json genes_1mb.json karyotype_zoom.json; do
    cp "$REF_JSON/$f" "$PUB/$f"
    log "  Copied $f"
done

fi # end stage 16
if (( FROM_STAGE <= 17 && TO_STAGE >= 17 )); then
# ── Stage 17: Publish results ───────────────────────────────
log "=== Stage 17: Publish results ==="

# Every sample publishes the same way: results go to private Blob storage keyed
# by the kit barcode, and the kit is flipped to complete. No git, no site
# rebuild, and the genomic data is never world-readable. The barcode is the
# sample's output_name upper-cased, so the sample sheet needs no extra column.
if (( PUBLISH_RESULTS )); then
    log "  Publishing $DOG_NAME to Blob storage"
    ( cd "$D/dogs-app" && node scripts/publish-results.mjs "$DOG_NAME" "$PUB" ) \
        || die "Publishing results for $DOG_NAME failed"
else
    # Compute nodes have no outbound internet. Results are staged; publish later
    # from a login node with cluster/publish-pending.sh.
    log "  PUBLISH_RESULTS=0 — results staged at $PUB, not published"
    echo "$DOG_NAME" > "$PUB/.pending-publish"
fi

log " Pipeline complete: $DOG_NAME"
log " Dashboard: kit $(echo "$DOG_NAME" | tr '[:lower:]' '[:upper:]')"
fi # end stage 17
echo "DONE" > "$OUT/pipeline.done"
