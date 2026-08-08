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

# MetaPhlAn's database is ~34GB extracted and is downloaded separately from the
# tool, so it lives beside the other reference data rather than inside the conda
# env. Passed to MetaPhlAn as --bowtie2db, so it is never re-downloaded.
METAPHLAN_DB="${METAPHLAN_DB:-$D/metaphlan_db}"

# Pin the index version. Without this MetaPhlAn queries the server for the
# "latest" index, decides whatever is in db_dir is "not present or partially
# present", and silently starts downloading a different ~40GB database — which
# would also make results incomparable with the Mac, which runs vJan25.
METAPHLAN_INDEX="${METAPHLAN_INDEX:-mpa_vJan25_CHOCOPhlAnSGB_202503}"
DATA_PYTHON="${DATA_PYTHON:-$D/envs/genomics/bin/python3}"

# Threads come from the scheduler: `qsub -pe shared N` exports NSLOTS.
NPROC="${NPROC:-${NSLOTS:-8}}"

# Match GLIMPSE parallelism to the slots the scheduler already gave us.
#
# Stage 4 (bwa-mem2) uses all NPROC slots regardless, so the job holds them for
# its whole life. Running GLIMPSE on one core left 7 of 8 slots idle for ~3.8h
# per sample — allocated to us and doing nothing. Using them costs no extra
# resource and cuts per-sample wall-clock roughly 4x (3h49m -> ~30min observed
# on the Mac at 8-way).
#
# Aggregate throughput across 100 samples is about the same either way; this
# just avoids wasting slots we are already holding. Set GLIMPSE_PARALLEL=1
# explicitly if you ever request a single-slot job.
GLIMPSE_PARALLEL="${GLIMPSE_PARALLEL:-${NSLOTS:-8}}"

# Work on node-local disk, not shared project storage. Each sample generates
# ~15-20GB of intermediates; 96 concurrent tasks writing that to /u/project would
# be slow for us and disruptive for everyone else on the filesystem. Only the
# result JSONs, coverage tracks, MetaPhlAn profile and logs are copied back.
USE_LOCAL_SCRATCH="${USE_LOCAL_SCRATCH:-1}"

# Compute nodes generally have no outbound internet, so the pipeline stops after
# staging results locally. Publish afterwards from a login node:
#     bash cluster/publish-pending.sh
PUBLISH_RESULTS="${PUBLISH_RESULTS:-0}"
