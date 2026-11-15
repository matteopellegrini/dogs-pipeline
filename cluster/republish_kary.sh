#!/bin/bash
# Silent republish after the 1,283-dog karyotype backfill. No emails: publish
# is silent by default; --notify is never passed here.
D=/u/project/pellegrini/gkislik/dogs
APP=$D/dogs-app
[[ -x "$D/envs/node/bin/node" ]] && export PATH="$D/envs/node/bin:$PATH"
ok=0; fail=0
while IFS=$'\t' read -r pk barcode; do
  [[ -n "$pk" && -n "$barcode" ]] || continue
  if ( cd "$APP" && node scripts/publish-results.mjs "$barcode" "$D/results_prosper/$pk" >/dev/null 2>&1 ); then
    ok=$((ok+1))
  else
    fail=$((fail+1)); echo "FAIL $pk $barcode"
  fi
done < "$D/cluster/publish_prosper.txt"
while IFS=$'\t' read -r low barcode; do
  [[ -n "$low" && -n "$barcode" ]] || continue
  if ( cd "$APP" && node scripts/publish-results.mjs "$barcode" "$D/results/$low" >/dev/null 2>&1 ); then
    ok=$((ok+1))
  else
    fail=$((fail+1)); echo "FAIL $low $barcode"
  fi
done < "$D/kit_map_96.tsv"
if ( cd "$APP" && node scripts/publish-results.mjs COSMO3 "$D/results/cosmo3" >/dev/null 2>&1 ); then
  ok=$((ok+1))
else
  fail=$((fail+1)); echo "FAIL cosmo3 COSMO3"
fi
echo "republish done: ok=$ok fail=$fail"
