#!/bin/bash
# Stage-12-only pass (depth-corrected inbreeding, analysis/inbreeding/inbreeding_stage.py) over re-imputed kits,
# from the archived v2 BCFs into their staging dirs (results_imp2/<pub basename>, or $OUT_ROOT/<basename>).
#   sbatch -A pi_mcdb -p pi_shared --qos=pi_shared24 -a 1-N%100 -c 2 --mem=8G -t 2:00:00 cluster/inbreeding_pass.sh <list>
# Rows: the reimpute list format (name, barcode, pub dir, work dir). Kits not yet re-imputed are skipped.
#SBATCH -J inbreeding_pass
#SBATCH -o logs/%x.%A.%a.log
export OPENBLAS_NUM_THREADS=1 PYTHONIOENCODING=utf-8
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"; cd "$D"
row=$(sed -n "${SLURM_ARRAY_TASK_ID:?}p" "${1:?list}"); [ -n "$row" ] || exit 0
name=$(echo "$row" | cut -f1); base=$(basename "$(echo "$row" | cut -f3)")
bcf=/u/project/pellegrini_archive/data/dogs_imputed_v2/$name.bcf
STAGE=$D/results_imp2/$base
[ -s "$STAGE/.done" ] && [ -s "$bcf" ] || { echo "IB-SKIP $name (not re-imputed)"; exit 0; }
OUT=${OUT_ROOT:+$OUT_ROOT/$base}; OUT=${OUT:-$STAGE}; mkdir -p "$OUT"
[ "$OUT" = "$STAGE" ] || cp -p "$STAGE/qc_result.json" "$OUT/"
H=$D/COSMO/glimpse2_dog10k/het_out
envs/genomics/bin/python3 analysis/inbreeding/inbreeding_stage.py "$bcf" "$OUT" --bcftools envs/genomics/bin/bcftools \
  --af-file $H/panel_af.tsv.gz --het-file $H/dog10k_het.het --calibration reference_json/inbreeding_depth_calibration.json \
  && echo "IB-DONE $name" || { echo "IB-FAIL $name"; exit 1; }
