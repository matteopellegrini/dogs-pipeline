#!/bin/bash
# Panel v3, step 3a: structural-variant read evidence (analysis/omia/sv_evidence.py) for every kit.
#   sbatch -A pi_matteop -p campus -a 1-99%60 -c 1 --mem=4G -t 2:00:00 -J pv3sv cluster/panel_v3_sv.sh work_panel_v3/manifest.tsv 20
# Output: work_panel_v3/sv/<pk>.sv.tsv (kit, site_id, bp_alt, bp_ref, disc, depth_in, depth_flank, len); kits whose
# output exists are skipped, so a resubmit only redoes failures.
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"
cd "$D"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PATH="$D/envs/genomics/bin:$PATH"
MANIFEST="${1:?manifest}"; PER="${2:-20}"
O="$D/work_panel_v3"; mkdir -p "$O/sv"
start=$(( (SLURM_ARRAY_TASK_ID - 1) * PER + 1 )); end=$(( start + PER - 1 ))
ok=0; fail=0; skip=0
while IFS=$'\t' read -r pk barcode pub work; do
  [[ -n "$pk" ]] || continue
  out="$O/sv/$pk.sv.tsv"
  if [[ -s "$out" ]]; then skip=$((skip+1)); continue; fi
  bam="$work/markdup.bam"
  if [[ ! -s "$bam" ]]; then echo "PV3SV-NOBAM $pk"; fail=$((fail+1)); continue; fi
  if python3 analysis/omia/sv_evidence.py count "$O/sv_sites.tsv" "$pk" "$bam" > "$out.tmp" 2> "$out.err"; then
    mv "$out.tmp" "$out"; rm -f "$out.err"; ok=$((ok+1))
  else
    echo "PV3SV-FAIL $pk $(tail -1 "$out.err")"; rm -f "$out.tmp"; fail=$((fail+1))
  fi
done < <(sed -n "${start},${end}p" "$MANIFEST")
echo "PV3SV-DONE task=$SLURM_ARRAY_TASK_ID rows=$start-$end ok=$ok fail=$fail skip=$skip"
