#!/bin/bash
# Silent bulk publish of batch-2 DW dogs as unowned holding kits.
# List-driven (cluster/publish_batch2.txt); removes .pending-publish markers
# on success so publish-pending.sh stays clean. No emails: silent by default.
D=/u/project/pellegrini/gkislik/dogs
APP=$D/dogs-app
[[ -x "$D/envs/node/bin/node" ]] && export PATH="$D/envs/node/bin:$PATH"
ok=0; fail=0
while IFS=$'\t' read -r pk barcode; do
  [[ -n "$pk" && -n "$barcode" ]] || continue
  if ( cd "$APP" && node scripts/publish-results.mjs "$barcode" "$D/results_prosper/$pk" >/dev/null 2>&1 ); then
    ok=$((ok+1)); rm -f "$D/results_prosper/$pk/.pending-publish"
  else
    fail=$((fail+1)); echo "FAIL $pk $barcode"
  fi
  (( (ok+fail) % 50 == 0 )) && echo "progress: $((ok+fail)) done (ok=$ok fail=$fail)"
done < "$D/cluster/publish_batch2.txt"
echo "batch2 publish done: ok=$ok fail=$fail"
