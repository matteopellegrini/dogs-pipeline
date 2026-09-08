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

# publish-results.mjs uses node:https (system node v16 suffices; fetch/undici
# cannot allocate its WASM parser under the 4GB login-node vmem cap). Prefer a
# bootstrapped $D/envs/node if present, else the system node.
[[ -x "$D/envs/node/bin/node" ]] && export PATH="$D/envs/node/bin:$PATH"
command -v node >/dev/null || { echo "ERROR: node not found"; exit 1; }

# Results may be staged anywhere under $D (the cluster keeps them outside the app
# checkout, unlike the Mac), so search rather than assuming one layout.
# QC holds first: samples the pipeline refused to queue (mean depth below
# QC_PUBLISH_MIN). They need an operator decision — resequence, or publish
# deliberately with: node scripts/publish-results.mjs <BARCODE> <dir>
mapfile -t holds < <(find "$D" -maxdepth 4 -name .qc-hold -type f 2>/dev/null)
if (( ${#holds[@]} > 0 )); then
  echo "=========================================================="
  echo "  QC HOLDS — ${#holds[@]} sample(s) NOT queued for publish:"
  for h in "${holds[@]}"; do printf '    %s  (%s)\n' "$(cat "$h")" "$(dirname "$h")"; done
  echo "  Review each: publish manually with publish-results.mjs, or resequence."
  echo "=========================================================="

  # Email the operator a triage report for each NEW hold, once. Only holds
  # with a .qc-hold.json (written by the pipeline's triage step) qualify;
  # a .qc-hold-notified sibling marks it as already reported. Bulk/backfill
  # runs must set SKIP_QC_EMAILS=1 so old batches never flood the inbox —
  # this notification path is for the ongoing new-customer flow.
  if [[ "${SKIP_QC_EMAILS:-0}" != "1" ]]; then
    QC_EMAIL_TO="${QC_EMAIL_TO:-matteope@gmail.com}"
    API_KEY="$(grep -m1 '^PIPELINE_API_KEY=' "$APP/.env.local" 2>/dev/null | cut -d= -f2- | tr -d '"' || true)"
    if [[ -z "$API_KEY" ]]; then
      echo "  (QC emails skipped: PIPELINE_API_KEY not found in $APP/.env.local)"
    else
      for h in "${holds[@]}"; do
        hd="$(dirname "$h")"
        hj="$hd/.qc-hold.json"
        smp="$(cut -f1 "$h" | head -1)"
        [[ -f "$hj" && ! -f "$hd/.qc-hold-notified" ]] || continue
        payload="$(python3 - "$hj" "$QC_EMAIL_TO" <<'PYEOF'
import json, sys
j = json.load(open(sys.argv[1]))
mf = j.get('dog_mapped_fraction')
body = f"""Sample {j['sample']} failed sequencing QC and was NOT published.

Verdict: {j['verdict']}
Mean depth: {j['mean_depth_x']}x (threshold {j['threshold_x']}x)
Duplication: {j['duplication_pct']:.1f}%
Read length: {j['read_length_raw_bp']:.0f} bp raw -> {j['read_length_trimmed_bp']:.0f} bp after trimming
Reads after QC: {j['total_reads']:,}""" + (f"\nDog-mapped fraction: {100*mf:.0f}%" if mf is not None else "") + f"""

Recommendation: {j['recommendation']}

The report is held on the cluster. To publish anyway:
  cd $D/dogs-app && node scripts/publish-results.mjs {j['barcode']} <results-dir>
"""
print(json.dumps({'to': sys.argv[2],
                  'subject': f"QC hold: {j['sample']} — {j['verdict']} ({j['mean_depth_x']}x)",
                  'text': body}))
PYEOF
)" || { echo "  (QC email compose failed for $hd)"; continue; }
        if [[ -n "$DRY" ]]; then
          echo "  [dry-run] would email QC hold for $smp"
        elif curl -sf -X POST "https://my.prosperk9.com/api/admin/send-email" \
               -H "x-api-key: $API_KEY" -H 'Content-Type: application/json' \
               -d "$payload" >/dev/null; then
          touch "$hd/.qc-hold-notified"
          echo "  QC email sent for $smp"
        else
          echo "  QC email FAILED for $smp (will retry next run)"
        fi
      done
    fi
  fi
fi

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
