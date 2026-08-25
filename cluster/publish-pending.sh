#!/usr/bin/env bash
# Publish every result set the cluster staged but could not upload.
#
#   bash cluster/publish-pending.sh [--dry-run]
#
# Run this from a Hoffman2 LOGIN node (compute nodes have no outbound internet,
# which is why the pipeline stages results and sets PUBLISH_RESULTS=0 there).
#
# Each staged sample leaves a .pending-publish marker containing its barcode.
# Publishing removes the marker, so re-running only picks up what is still
# outstanding and is safe to run repeatedly.
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOGS_SITE="${DOGS_SITE:-hoffman}"
source "$PIPELINE_DIR/site/${DOGS_SITE}.sh"

DRY=""
[[ "${1:-}" == "--dry-run" ]] && DRY="--dry-run"

APP="$D/dogs-app"
[[ -d "$APP" ]] || { echo "ERROR: no app checkout at $APP — clone dogs-app there first (it is private)"; exit 1; }

# publish-results.mjs needs global fetch (node >= 18); Hoffman's system node is
# v16 and there is no nodejs module, so bootstrap installs one at $D/envs/node.
[[ -x "$D/envs/node/bin/node" ]] && export PATH="$D/envs/node/bin:$PATH"
node -e 'if (typeof fetch !== "function") process.exit(1)' \
  || { echo "ERROR: node $(node --version 2>/dev/null) lacks fetch — need >= 18 (mamba create -p \$D/envs/node -c conda-forge nodejs=20)"; exit 1; }

# Results may be staged anywhere under $D (the cluster keeps them outside the app
# checkout, unlike the Mac), so search rather than assuming one layout.
mapfile -t markers < <(find "$D" -maxdepth 4 -name .pending-publish -type f 2>/dev/null)
if (( ${#markers[@]} == 0 )); then
  echo "Nothing pending — no .pending-publish markers found."
  exit 0
fi

echo "Found ${#markers[@]} sample(s) awaiting publish."
failed=0
for m in "${markers[@]}"; do
  dir="$(dirname "$m")"
  barcode="$(tr -d '[:space:]' < "$m")"
  [[ -n "$barcode" ]] || { echo "  SKIP $dir — empty marker"; continue; }
  printf '  %-16s ' "$barcode"
  if ( cd "$APP" && node scripts/publish-results.mjs "$barcode" "$dir" $DRY >/dev/null 2>&1 ); then
    echo "published"
    [[ -z "$DRY" ]] && rm -f "$m"
  else
    echo "FAILED — rerun for detail: (cd $APP && node scripts/publish-results.mjs $barcode $dir)"
    failed=$((failed + 1))
  fi
done

(( failed == 0 )) || { echo "$failed sample(s) failed; markers left in place."; exit 1; }
echo "All pending samples published."
