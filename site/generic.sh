# Site profile: fallback for anywhere that isn't the Mac or Hoffman2.
#
# Deliberately assumes nothing beyond $HOME. If the auto-detection picked this
# and it's wrong, set DOGS_SITE explicitly:
#     DOGS_SITE=hoffman bash run_dog_pipeline.sh sample_sheet.tsv 2 1
#
# Preflight will name whatever is missing before any work starts.

D="${D:-$HOME/dogs}"

ENV_GENOMICS="${ENV_GENOMICS:-$HOME/micromamba/envs/genomics}"
ENV_GLIMPSE="${ENV_GLIMPSE:-$HOME/micromamba/envs/glimpse}"
METAPHLAN_BIN="${METAPHLAN_BIN:-$ENV_GENOMICS/bin/metaphlan}"
DATA_PYTHON="${DATA_PYTHON:-$ENV_GENOMICS/bin/python3}"

NPROC="${NPROC:-4}"
GLIMPSE_PARALLEL="${GLIMPSE_PARALLEL:-2}"

# Conservative: stage results rather than publishing from an unknown machine.
PUBLISH_RESULTS="${PUBLISH_RESULTS:-0}"
