#!/bin/bash
# Promote panel v2 re-imputed reports from results_imp2 staging into the published result dirs.
# Only the files stages 7-14 rewrite are copied (a survey of 1,875 staged kits found exactly these differ);
# everything else in staging is the copy reimpute.sh took at job start, and must not overwrite anything
# published since (karyotype, microbiome). Each replaced file is kept once as <file>.pre-<tag>.json, and
# breed_result.json is always backed up and touched for a promoted kit, so
#   FILES="$PROMOTE_FILES" bash cluster/republish_files.sh <tag> <list>
# republishes every promoted kit (republish_files.sh keys on the first named file).
# Kits with a .qc-hold.json in staging, or without .done, are skipped.
#   DRY=1 bash cluster/promote_imp2.sh <tag> <list>     rows: name, barcode, pub dir, work dir
D=/u/project/pellegrini/gkislik/dogs; cd $D
TAG=${1:?tag}; LIST=${2:?list}
PROMOTE_FILES="breed_result.json local_ancestry.json omia_result.json coat_color.json traits_result.json prs_result.json relatives_result.json inbreeding_result.json inbreeding_froh_dog10k_result.json mhc_result.json functional_variants.json ortho_result.json cnv_genes.json cnv_homdel.json"
kits=0; files=0; hold=0; skip=0
while IFS=$'\t' read -r name barcode pub work; do
  [[ -n "$name" && -d "$pub" ]] || continue
  st=$D/results_imp2/$(basename "$pub")
  [[ -s "$st/.done" ]] || { skip=$((skip+1)); echo "SKIP $name (not re-imputed)"; continue; }
  [[ -e "$st/.qc-hold.json" ]] && { hold=$((hold+1)); echo "HOLD $name (.qc-hold.json in staging)"; continue; }
  n=0
  for f in $PROMOTE_FILES; do
    [[ -f "$st/$f" ]] || continue
    cmp -s "$st/$f" "$pub/$f" && continue
    n=$((n+1)); [[ -n "${DRY:-}" ]] && continue
    b="$pub/${f%.json}.pre-$TAG.json"
    [[ -f "$pub/$f" && ! -f "$b" ]] && cp -p "$pub/$f" "$b"
    cp "$st/$f" "$pub/$f"
  done
  (( n > 0 )) || continue
  kits=$((kits+1)); files=$((files+n))
  if [[ -z "${DRY:-}" ]]; then
    [[ -f "$pub/breed_result.pre-$TAG.json" ]] || cp -p "$pub/breed_result.json" "$pub/breed_result.pre-$TAG.json"
    touch "$pub/breed_result.json"
  fi
done < "$LIST"
echo "promote $TAG${DRY:+ (dry run)}: kits $kits, files $files, held $hold, not re-imputed $skip"
echo "republish with: FILES=\"$PROMOTE_FILES\" bash cluster/republish_files.sh $TAG $LIST"
