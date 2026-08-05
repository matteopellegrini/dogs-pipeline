# Site profile: UCLA Hoffman2 (UGE/SGE cluster, many samples in parallel).
#
# Set up once with cluster/bootstrap-hoffman.sh, which creates the environments
# and stages the ~37GB of reference data these paths point at.

# Lab project storage: persistent (unlike /scratch, which is purged) and large.
# $HOME cannot be used — its 60GB quota is mostly consumed, leaving far less than
# the ~43GB of reference data this pipeline needs.
D="${D:-/u/project/pellegrini/$USER/dogs}"

ENV_GENOMICS="${ENV_GENOMICS:-$D/envs/genomics}"
ENV_GLIMPSE="${ENV_GLIMPSE:-$D/envs/glimpse}"
METAPHLAN_BIN="${METAPHLAN_BIN:-$D/envs/genomics/bin/metaphlan}"
DATA_PYTHON="${DATA_PYTHON:-$D/envs/genomics/bin/python3}"

# Threads come from the scheduler: `qsub -pe shared N` exports NSLOTS.
NPROC="${NPROC:-${NSLOTS:-8}}"

# GLIMPSE parallelism is 1 on purpose. Each array task is one sample, so
# parallelism belongs ACROSS samples — 100 tasks x 8 GLIMPSE jobs would
# oversubscribe the node badly and slow everything down.
GLIMPSE_PARALLEL="${GLIMPSE_PARALLEL:-1}"

# Compute nodes generally have no outbound internet, so the pipeline stops after
# staging results locally. Publish afterwards from a login node:
#     bash cluster/publish-pending.sh
PUBLISH_RESULTS="${PUBLISH_RESULTS:-0}"
