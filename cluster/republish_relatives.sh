#!/bin/bash
# Silent republish after the stage-9b v3 repass: every row of
# cluster/repass9b_list.txt whose relatives_result.json was rewritten
# (REL9B-DONE in its task log). Silent by default — no --notify.
D=/u/project/pellegrini/gkislik/dogs
APP=$D/dogs-app
[[ -x "$D/envs/node/bin/node" ]] && export PATH="$D/envs/node/bin:$PATH"
ok=0; fail=0; skip=0
while IFS=$'\t' read -r pk barcode pub work; do
  [[ -n "$pk" && -n "$barcode" ]] || continue
  if [[ ! -f "$pub/relatives_result.pre-v3.json" ]]; then skip=$((skip+1)); continue; fi
  if ( cd "$APP" && node scripts/publish-results.mjs "$barcode" "$pub" >/dev/null 2>&1 ); then
    ok=$((ok+1))
  else
    fail=$((fail+1)); echo "FAIL $pk $barcode"
  fi
  (( (ok+fail) % 100 == 0 )) && echo "progress: $((ok+fail)) (ok=$ok fail=$fail)"
done < "$D/cluster/repass9b_list.txt"
echo "relatives republish done: ok=$ok fail=$fail skipped(not repassed)=$skip"
