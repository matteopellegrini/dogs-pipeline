#!/usr/bin/env bash
# Resume pipeline for a dog from Stage 8 onwards (GLIMPSE2 already done).
# Usage: bash resume_from_stage8.sh <DogName>
set -euo pipefail

DOG_NAME="${1:?Usage: $0 <DogName>}"
DOG_LOWER=$(echo "$DOG_NAME" | tr '[:upper:]' '[:lower:]')

D=/Users/matteopellegrini/Downloads/dogs
OUT=$D/$DOG_NAME/analysis
PUB=$D/dogs-app/public/$DOG_LOWER
FASTA=$D/canFam4.fa
DOG10K_PANEL=$D/dog10k_panel/AutoAndXPAR.Dog10K.phased_plus_disease_rh.bcf
OMIA_DB=$D/dogs-app/public/cosmo/omia_result.json
COSMO_PUB=$D/dogs-app/public/cosmo
SCOPE_P=$D/COSMO/analysis/cosmo_scope_fullPhat.txt
SCOPE_CLUST=$D/COSMO/analysis/scope_clust.txt
PARKER_BIM=$D/COSMO/analysis/cosmo_parker_scope.bim
PARKER_FAM=$D/COSMO/analysis/cosmo_parker_scope.fam
SNPEFF_DB="ROS_Cfam_1.115"
MM="micromamba run -n genomics"
MM_GLIMPSE="micromamba run -n glimpse_x86"
NPROC=8
IMPUTED_BCF="$OUT/glimpse2/${DOG_LOWER}_imputed_dog10k.bcf"
ANN_DIR="$OUT/snpeff"

mkdir -p "$PUB" "$ANN_DIR"
LOG=$OUT/pipeline_resume.log

# Define log() so it's available in this shell and any sourced blocks
log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
export -f log

# Export every variable that the stage scripts reference via $VAR expansion
export DOG_NAME DOG_LOWER D OUT PUB FASTA DOG10K_PANEL OMIA_DB COSMO_PUB \
       SCOPE_P SCOPE_CLUST PARKER_BIM PARKER_FAM SNPEFF_DB \
       MM MM_GLIMPSE NPROC IMPUTED_BCF ANN_DIR LOG

START_LINE="${2:-291}"   # default Stage 8 (291); pass 468 for Stage 9, 583 for Stage 10

log "========================================"
log "Resuming $DOG_NAME pipeline from line $START_LINE"
log "========================================"

TMP=$(mktemp /tmp/resume_stages_XXXXXX)  # no .sh suffix — macOS mktemp needs Xs at end
sed -n "${START_LINE},1099p" "$D/run_dog_pipeline.sh" > "$TMP"
bash "$TMP"
rm -f "$TMP"
