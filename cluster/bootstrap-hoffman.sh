#!/usr/bin/env bash
# One-time Hoffman2 setup: tool environments + reference data.
#
#   bash cluster/bootstrap-hoffman.sh            # everything
#   bash cluster/bootstrap-hoffman.sh --tools    # environments only
#   bash cluster/bootstrap-hoffman.sh --check    # verify, install nothing
#
# Run from a LOGIN node — this needs outbound internet. No root required:
# micromamba unpacks into your own $SCRATCH, which is how installs work on a
# cluster you don't administer.
#
# Sizes: tools ~6GB, reference data ~37GB. The Dog10K panel (23GB) is the bulk
# and must be copied from your Mac; everything else downloads.
set -euo pipefail

# Everything lives together on persistent lab project storage. /scratch is purged
# and $HOME is too small for the ~43GB of reference data.
D="${D:-/u/project/pellegrini/dogs}"
ENVS="$D/envs"
MAMBA="$D/bin/micromamba"
export MAMBA_ROOT_PREFIX="$D/micromamba"

MODE="${1:-all}"

step() { printf '\n\033[1m== %s\033[0m\n' "$*"; }
ok()   { printf '   ok   %s\n' "$*"; }
bad()  { printf '   MISSING %s\n' "$*"; }

# ── verify ───────────────────────────────────────────────────
check() {
  local missing=0 t f
  step "Tools"
  for t in "$ENVS/genomics/bin/bwa-mem2" "$ENVS/genomics/bin/samtools" \
           "$ENVS/genomics/bin/bcftools" "$ENVS/genomics/bin/fastp" \
           "$ENVS/genomics/bin/snpEff"   "$ENVS/genomics/bin/metaphlan" \
           "$ENVS/glimpse/bin/GLIMPSE2_phase" "$ENVS/glimpse/bin/GLIMPSE2_ligate"; do
    [[ -x "$t" ]] && ok "$t" || { bad "$t"; missing=1; }
  done
  step "Python data stack"
  if "$ENVS/genomics/bin/python3" -c 'import pandas,numpy,scipy,sklearn' 2>/dev/null; then
    ok "pandas/numpy/scipy/sklearn"
  else bad "pandas/numpy/scipy/sklearn in $ENVS/genomics"; missing=1; fi
  step "Reference data under $D"
  for f in canFam4.fa canFam4_idx.bwt.2bit.64 \
           dog10k_panel/AutoAndXPAR.Dog10K.phased_plus_disease_rh.bcf \
           COSMO/glimpse2_dog10k/chunks COSMO/analysis/cosmo_parker_full.bed \
           COSMO/analysis/cosmo_scope177Phat.txt \
           metagenome/merged_microbiome_age_weight_3.18_final.csv \
           reference_json/centromeres.json; do
    [[ -e "$D/$f" ]] && ok "$f" || { bad "$f"; missing=1; }
  done
  return $missing
}

if [[ "$MODE" == "--check" ]]; then
  check && { echo; echo "Bootstrap complete — ready to submit jobs."; exit 0; } \
        || { echo; echo "Incomplete: see MISSING above."; exit 1; }
fi

# ── tools ────────────────────────────────────────────────────
step "micromamba"
if [[ -x "$MAMBA" ]]; then
  ok "already at $MAMBA"
else
  mkdir -p "$D/bin"
  ( cd "$D" && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba )
  ok "installed"
fi

step "genomics env"
# Versions pinned to what the Mac install runs today, so the cluster reproduces
# it rather than whatever the solver picks months from now.
"$MAMBA" create -y -p "$ENVS/genomics" -c conda-forge -c bioconda \
  bwa-mem2=2.3 samtools=1.23.1 bcftools=1.23.1 fastp=1.3.6 \
  snpeff=5.4.0c openjdk=25 metaphlan=4.2.4 \
  python=3.11 pandas numpy scipy scikit-learn

step "glimpse env"
"$MAMBA" create -y -p "$ENVS/glimpse" -c conda-forge -c bioconda glimpse-bio=2.0.1

[[ "$MODE" == "--tools" ]] && { echo; echo "Tools installed. Re-run without --tools for reference data."; exit 0; }

# ── reference data ───────────────────────────────────────────
step "MetaPhlAn database (~8GB, downloads)"
if [[ -d "$ENVS/genomics/lib/python3.11/site-packages/metaphlan/metaphlan_databases" ]]; then
  ok "present"
else
  "$ENVS/genomics/bin/metaphlan" --install --nproc 4 || \
    echo "   WARNING: MetaPhlAn DB install failed — retry from a login node with internet"
fi

step "SnpEff dog database"
"$ENVS/genomics/bin/snpEff" download -v ROS_Cfam_1.115 || \
  echo "   WARNING: snpEff download failed — copy \$ENV/share/snpeff*/data from the Mac instead"

step "Reference data that must be copied from the Mac"
mkdir -p "$D"
cat <<EOF

   These cannot be downloaded — copy them from your Mac (~26GB total).
   From the Mac, with <user> your Hoffman2 username:

     rsync -avP --info=progress2 \\
       ~/Downloads/dogs/canFam4.fa \\
       ~/Downloads/dogs/canFam4.fa.fai \\
       ~/Downloads/dogs/canFam4_idx* \\
       ~/Downloads/dogs/dog10k_panel \\
       ~/Downloads/dogs/reference_json \\
       ~/Downloads/dogs/metagenome \\
       <user>@hoffman2.idre.ucla.edu:$D/

     rsync -avR --info=progress2 \\
       ~/Downloads/dogs/./COSMO/glimpse2_dog10k/chunks \\
       ~/Downloads/dogs/./COSMO/analysis/cosmo_parker_full.* \\
       ~/Downloads/dogs/./COSMO/analysis/cosmo_scope177Phat.txt \\
       ~/Downloads/dogs/./COSMO/analysis/scope_clust.txt \\
       <user>@hoffman2.idre.ucla.edu:$D/

   The Dog10K panel is 23GB of that and is the slow part.

EOF

step "Verifying"
check || { echo; echo "Still incomplete — finish the rsync above, then: bash cluster/bootstrap-hoffman.sh --check"; exit 1; }
echo
echo "Bootstrap complete. Test one sample before any array job:"
echo "  qsub -t 2-2 cluster/submit-array.sh sample_sheet.tsv"
