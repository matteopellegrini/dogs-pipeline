#!/bin/bash
# Silent republish of the kits in cluster/coatconf_kits.txt after the stage-13
# confidence-cap rerun (2026-09-13). Only rows whose task log says RERUN13-DONE.
D=/u/project/pellegrini/gkislik/dogs; APP=$D/dogs-app
[[ -x "$D/envs/node/bin/node" ]] && export PATH="$D/envs/node/bin:$PATH"
done_kits=$(grep -h "RERUN13-DONE" logs/coat13-*.1473339[23].*.log 2>/dev/null | awk "{print \$2}" | sort -u)
ok=0; fail=0; skip=0
while read -r pk; do
  [[ -n "$pk" ]] || continue
  grep -qx "$pk" <<<"$done_kits" || { skip=$((skip+1)); echo "SKIP $pk (no RERUN13-DONE)"; continue; }
  barcode=${pk#pk-}; pub=$D/results_prosper/$pk
  if ( cd "$APP" && node scripts/publish-results.mjs "$barcode" "$pub" >/dev/null 2>&1 ); then ok=$((ok+1)); else fail=$((fail+1)); echo "FAIL $pk"; fi
done < "$D/cluster/coatconf_kits.txt"
echo "coat republish done: ok=$ok fail=$fail skipped=$skip"
