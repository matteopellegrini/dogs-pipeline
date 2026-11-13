#!/bin/bash
# Rerun stage 12b (mhc_result.json) in results_imp2 staging after the MHC reference refit
# (analysis/mhc/refit_from_staging.py), from each kit's archived panel v2 BCF, its staged breed_result.json and
# its depth-corrected inbreeding_result.json. ~1 s per dog; run on the login node with nohup.
#   bash cluster/mhc_restage_imp2.sh <list>      rows: name, barcode, pub dir, work dir
D=/u/project/pellegrini/gkislik/dogs; cd $D; export LD_LIBRARY_PATH=$D/envs/genomics/lib
ok=0; fail=0; skip=0
while IFS=$'\t' read -r name barcode pub work; do
  [[ -n "$name" ]] || continue
  bcf=/u/project/pellegrini_archive/data/dogs_imputed_v2/$name.bcf
  st=$D/results_imp2/$(basename "$pub")
  [[ -s "$bcf" && -s "$st/.done" ]] || { skip=$((skip+1)); echo "SKIP $name (not re-imputed)"; continue; }
  if PYTHONIOENCODING=utf-8 envs/genomics/bin/python3 analysis/mhc/mhc_stage.py "$bcf" reference_json/mhc_reference.json "$st/mhc_result.json" --breed "$st/breed_result.json" --inbreeding "$st/inbreeding_result.json" >/dev/null 2>&1; then
    ok=$((ok+1)); else fail=$((fail+1)); echo "FAIL $name"; fi
  (( (ok+fail) % 200 == 0 )) && echo "progress $((ok+fail)) ok=$ok fail=$fail"
done < "${1:?list}"
echo "MHC-RESTAGE-DONE ok=$ok fail=$fail skip=$skip"
