#!/bin/bash
# Silent republish after the stage-11 weight-model rerun (2026-09-13): every row
# in cluster/w11_{prosper,batch2}_rows.txt whose prs_result.json now carries the
# rebuilt model note. Usage: bash cluster/republish_w11.sh > logs/republish_w11.log
D=/u/project/pellegrini/gkislik/dogs; APP=$D/dogs-app
[[ -x "$D/envs/node/bin/node" ]] && export PATH="$D/envs/node/bin:$PATH"
ok=0; fail=0; skip=0
for pair in sample_sheet.prosper.tsv:cluster/w11_prosper_rows.txt sample_sheet.batch2.tsv:cluster/w11_batch2_rows.txt; do
  sheet=${pair%%:*}; rows=${pair##*:}
  while read -r r; do
    [[ -n "$r" ]] || continue
    line=$(awk -v i="$r" "NR==i{print; exit}" "$D/$sheet")
    name=$(cut -f4 <<<"$line"); pub=$(cut -f6 <<<"$line"); barcode=${name#pk-}
    if ! grep -q "577 ProsperK9 customer-reported" "$pub/prs_result.json" 2>/dev/null; then skip=$((skip+1)); echo "SKIP $name (old model still)"; continue; fi
    if ( cd "$APP" && node scripts/publish-results.mjs "$barcode" "$pub" >/dev/null 2>&1 ); then ok=$((ok+1)); else fail=$((fail+1)); echo "FAIL $name"; fi
    (( (ok+fail) % 100 == 0 )) && echo "progress: $((ok+fail)) (ok=$ok fail=$fail)"
  done < "$D/$rows"
done
echo "w11 republish done: ok=$ok fail=$fail skipped=$skip"
