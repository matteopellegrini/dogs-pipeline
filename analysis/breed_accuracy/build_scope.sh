#!/usr/bin/env bash
# Build SCOPE (https://github.com/sriramlab/SCOPE) for supervised breed inference.
#
#   bash analysis/breed_accuracy/build_scope.sh [install-dir]
#
# On Linux (Hoffman) the upstream instructions work as published and this script
# is just a wrapper. On an Apple Silicon Mac three things need fixing, none of
# them obvious from the error output:
#
#   1. cmake >= 4 rejects the project's old cmake_minimum_required.
#   2. include/mailman.h uses SSE2 intrinsics (__m128d, _mm_loadu_pd) directly,
#      so it cannot compile for arm64 at all. Building x86_64 under Rosetta is
#      the path of least resistance and matches how GLIMPSE2 already runs here.
#   3. The vendored Eigen 3.3.3 has a known incompatibility with modern clang in
#      Transpositions.h: it calls trt.derived() on a Transpose&, which has no
#      such member. Upstream Eigen fixed this by passing trt directly.
set -euo pipefail

DEST="${1:-$HOME/opt}"
mkdir -p "$DEST"
cd "$DEST"

if [[ ! -d SCOPE ]]; then
  git clone --depth 1 https://github.com/sriramlab/SCOPE.git
fi
cd SCOPE

# Fix 3 — vendored Eigen vs modern clang. Idempotent.
EIGEN_FILE=include/Eigen/src/Core/Transpositions.h
if grep -q 'matrix.derived(), trt.derived()' "$EIGEN_FILE"; then
  perl -pi -e 's/matrix\.derived\(\), trt\.derived\(\)/matrix.derived(), trt/' "$EIGEN_FILE"
  echo "patched vendored Eigen ($EIGEN_FILE)"
fi

CMAKE_ARGS=(-DCMAKE_POLICY_VERSION_MINIMUM=3.5)          # fix 1
if [[ "$(uname -s)" == Darwin && "$(uname -m)" == arm64 ]]; then
  CMAKE_ARGS+=(-DCMAKE_OSX_ARCHITECTURES=x86_64)         # fix 2
  echo "Apple Silicon detected — building x86_64 (runs under Rosetta)"
fi

rm -rf build && mkdir build && cd build
cmake "${CMAKE_ARGS[@]}" ..
make -j4

echo
echo "built: $PWD/scope"
"$PWD/scope" 2>&1 | head -3 || true
echo
echo "Supervised usage:"
echo "  scope -g <plink-prefix> -freq <plink .frq.strat> -k <n-breeds> -o <out-prefix> -nt 4"
echo "SNP order and allele coding MUST match between the frequency file and the"
echo "genotypes — analysis/breed_accuracy/make_merged_plink.py writes both together"
echo "for exactly that reason."
