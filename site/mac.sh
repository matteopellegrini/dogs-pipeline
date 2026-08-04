# Site profile: Matteo's Mac (interactive, one sample at a time).
#
# Sourced by run_dog_pipeline.sh before anything else. Every value here is a
# default the main script honours via ${VAR:-...}, so this profile reproduces
# exactly the behaviour the pipeline had before site profiles existed.

D="${D:-$HOME/Downloads/dogs}"

# Tool environments (micromamba, created by hand on this machine)
ENV_GENOMICS="${ENV_GENOMICS:-$HOME/micromamba/envs/genomics}"
ENV_GLIMPSE="${ENV_GLIMPSE:-$HOME/micromamba/envs/glimpse_x86}"
METAPHLAN_BIN="${METAPHLAN_BIN:-$HOME/Library/Python/3.9/bin/metaphlan}"

# Interpreter holding pandas/numpy/scipy/sklearn. The genomics env's python
# does NOT have them, and Stage 10 prepends that env to PATH — so this must be
# an absolute path, not a bare `python3`.
DATA_PYTHON="${DATA_PYTHON:-/usr/bin/python3}"

# One sample at a time, so parallelism lives *within* the sample.
NPROC="${NPROC:-8}"
GLIMPSE_PARALLEL="${GLIMPSE_PARALLEL:-8}"

# This machine has outbound internet, so Stage 17 publishes directly.
PUBLISH_RESULTS="${PUBLISH_RESULTS:-1}"
