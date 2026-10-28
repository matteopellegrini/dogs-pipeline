#!/bin/bash
# Silent republish after a stage-8 catalogue rerun (cluster/rerun-8.sh): every row whose omia_result.json is newer
# than its omia_result.pre-<tag>.json backup. No emails (publish-results.mjs is silent unless --notify).
# Run from a login node:  bash cluster/republish_omia.sh <tag> <list> [<list> ...]
D=/u/project/pellegrini/gkislik/dogs
APP=$D/dogs-app
[[ -x "$D/envs/node/bin/node" ]] && export PATH="$D/envs/node/bin:$PATH"
TAG="${1:?tag}"; shift
ok=0; fail=0; skip=0
for LIST in "$@"; do
  while IFS=$'\t' read -r pk barcode pub work; do
    [[ -n "$pk" && -n "$barcode" ]] || continue
    if [[ ! -f "$pub/omia_result.pre-$TAG.json" || ! "$pub/omia_result.json" -nt "$pub/omia_result.pre-$TAG.json" ]]; then skip=$((skip+1)); continue; fi
    if ( cd "$APP" && node scripts/publish-results.mjs "$barcode" "$pub" >/dev/null 2>&1 ); then ok=$((ok+1)); else fail=$((fail+1)); echo "FAIL $pk $barcode"; fi
    (( (ok+fail) % 100 == 0 )) && echo "progress: $((ok+fail)) (ok=$ok fail=$fail)"
  done < "$LIST"
done
echo "omia republish ($TAG) done: ok=$ok fail=$fail skipped=$skip"
