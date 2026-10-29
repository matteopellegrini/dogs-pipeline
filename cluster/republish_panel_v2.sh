#!/bin/bash
# Silent republish after the breed panel v2 rerun: every row of the two rerun lists whose
# breed_result.json is newer than its pre-v2 copy (i.e. stage 9 succeeded). No emails
# (publish-results.mjs is silent unless --notify). Run from a login node.
#   bash cluster/republish_panel_v2.sh cluster/repass9b_list.txt cluster/rerun9_results_list.txt
# TAG=v21 for the v2.1 rerun (compares against breed_result.pre-v21.json).
D=/u/project/pellegrini/gkislik/dogs
APP=$D/dogs-app
[[ -x "$D/envs/node/bin/node" ]] && export PATH="$D/envs/node/bin:$PATH"
TAG="${TAG:-v2}"
ok=0; fail=0; skip=0
for LIST in "$@"; do
  while IFS=$'\t' read -r pk barcode pub work; do
    [[ -n "$pk" && -n "$barcode" ]] || continue
    if [[ ! -f "$pub/breed_result.json" || ! -f "$pub/breed_result.pre-$TAG.json" || ! "$pub/breed_result.json" -nt "$pub/breed_result.pre-$TAG.json" ]]; then skip=$((skip+1)); echo "SKIP $pk (no fresh $TAG result)"; continue; fi
    if ( cd "$APP" && node scripts/publish-results.mjs "$barcode" "$pub" >/dev/null 2>&1 ); then ok=$((ok+1)); else fail=$((fail+1)); echo "FAIL $pk $barcode"; fi
    (( (ok+fail) % 100 == 0 )) && echo "progress: $((ok+fail)) (ok=$ok fail=$fail)"
  done < "$LIST"
done
echo "panel $TAG republish done: ok=$ok fail=$fail skipped=$skip"
