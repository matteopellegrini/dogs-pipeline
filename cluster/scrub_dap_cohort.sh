#!/bin/bash
# Remove the Dog Aging Project name from published report text (2026-11-11). Originals of every touched file go
# to logs/nodap_originals.<time>.tgz; breed_result.json is also copied to breed_result.pre-nodap.json, the marker
# republish_files.sh uses. Then scrub breed/prs/relatives JSONs in place (string values only).
#   DRY=1 bash cluster/scrub_dap_cohort.sh <list>     then
#   FILES="breed_result.json prs_result.json relatives_result.json" bash cluster/republish_files.sh nodap <list>
D=/u/project/pellegrini/gkislik/dogs; cd $D
PY=$D/envs/genomics/bin/python3
LIST=${1:?list}
hits=$(mktemp)
cut -f3 "$LIST" | sort -u | while read -r pub; do
  [ -d "$pub" ] || continue
  grep -l -i "dog aging project" "$pub"/breed_result.json "$pub"/prs_result.json "$pub"/relatives_result.json 2>/dev/null
done > "$hits"
echo "kits with the name: $(xargs -n1 dirname < "$hits" | sort -u | wc -l)  files: $(wc -l < "$hits")"
[ -n "${DRY:-}" ] && { rm -f "$hits"; exit 0; }
tar czf "logs/nodap_originals.$(date +%Y%m%d%H%M).tgz" -T "$hits" 2>/dev/null
for pub in $(xargs -n1 dirname < "$hits" | sort -u); do
  [ -f "$pub/breed_result.pre-nodap.json" ] || cp -p "$pub/breed_result.json" "$pub/breed_result.pre-nodap.json"
done
xargs $PY -I analysis/imputation_panel/scrub_dap.py < "$hits" | awk '{n+=$2} END {print "replacements:", n, "in", NR, "files"}'
for pub in $(xargs -n1 dirname < "$hits" | sort -u); do touch "$pub/breed_result.json"; done
left=$(xargs grep -l -i "dog aging project" < "$hits" | wc -l)
echo "files still naming it: $left"
rm -f "$hits"
