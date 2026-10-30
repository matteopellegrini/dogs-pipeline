#!/bin/bash
# Silent republish of NAMED result files only (others already published stay untouched): each kit's files are copied
# into an empty staging dir and published from there (publish-results.mjs uploads per file, never deletes). Only kits
# whose first named file is newer than its <file>.pre-<tag>.json backup.
#   FILES="coat_color.json traits_result.json" bash cluster/republish_files.sh <tag> <list> [<list> ...]
D=/u/project/pellegrini/gkislik/dogs
APP=$D/dogs-app
[[ -x "$D/envs/node/bin/node" ]] && export PATH="$D/envs/node/bin:$PATH"
TAG="${1:?tag}"; shift
read -r -a NAMES <<< "${FILES:?FILES=\"a.json b.json\" required}"
first=${NAMES[0]%.json}
stage=$(mktemp -d "$D/logs/files_only.XXXXXX")
ok=0; fail=0; skip=0
for LIST in "$@"; do
  while IFS=$'\t' read -r pk barcode pub work; do
    [[ -n "$pk" && -n "$barcode" ]] || continue
    if [[ ! -f "$pub/$first.pre-$TAG.json" || ! "$pub/$first.json" -nt "$pub/$first.pre-$TAG.json" ]]; then skip=$((skip+1)); continue; fi
    rm -f "$stage"/*.json
    for f in "${NAMES[@]}"; do [[ -f "$pub/$f" ]] && cp "$pub/$f" "$stage/"; done
    if ( cd "$APP" && node scripts/publish-results.mjs "$barcode" "$stage" >/dev/null 2>&1 ); then ok=$((ok+1)); else fail=$((fail+1)); echo "FAIL $pk $barcode"; fi
    (( (ok+fail) % 100 == 0 )) && echo "progress: $((ok+fail)) (ok=$ok fail=$fail)"
  done < "$LIST"
done
rm -rf "$stage"
echo "files-only republish ($TAG: ${NAMES[*]}) done: ok=$ok fail=$fail skipped=$skip"
