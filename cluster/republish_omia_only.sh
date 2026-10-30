#!/bin/bash
# Silent republish of omia_result.json ALONE (health conditions) after a stage-8 rerun, leaving every other published
# file as it is: each kit's omia_result.json is copied into an empty staging dir and published from there
# (publish-results.mjs uploads per file and never deletes). Used 2026-10-29 when the pub dirs also held an
# unpublished stage-9 rerun. Only rows whose omia_result.json is newer than omia_result.pre-<tag>.json.
#   bash cluster/republish_omia_only.sh <tag> <list> [<list> ...]
D=/u/project/pellegrini/gkislik/dogs
APP=$D/dogs-app
[[ -x "$D/envs/node/bin/node" ]] && export PATH="$D/envs/node/bin:$PATH"
TAG="${1:?tag}"; shift
stage=$(mktemp -d "$D/logs/omia_only.XXXXXX")
ok=0; fail=0; skip=0
for LIST in "$@"; do
  while IFS=$'\t' read -r pk barcode pub work; do
    [[ -n "$pk" && -n "$barcode" ]] || continue
    if [[ ! -f "$pub/omia_result.pre-$TAG.json" || ! "$pub/omia_result.json" -nt "$pub/omia_result.pre-$TAG.json" ]]; then skip=$((skip+1)); continue; fi
    rm -f "$stage"/*.json; cp "$pub/omia_result.json" "$stage/"
    if ( cd "$APP" && node scripts/publish-results.mjs "$barcode" "$stage" >/dev/null 2>&1 ); then ok=$((ok+1)); else fail=$((fail+1)); echo "FAIL $pk $barcode"; fi
    (( (ok+fail) % 100 == 0 )) && echo "progress: $((ok+fail)) (ok=$ok fail=$fail)"
  done < "$LIST"
done
rm -rf "$stage"
echo "omia-only republish ($TAG) done: ok=$ok fail=$fail skipped=$skip"
