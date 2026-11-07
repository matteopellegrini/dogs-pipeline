#!/bin/bash
# Panel v3: depth evidence at the depth-only CNV sites (known_variants/panel_v3/depth_only_cnv.tsv) for every kit,
# appended to the cohort SV counts. Same per-kit loop as panel_v3_sv.sh, restricted to the depth-only rows of
# sv_sites.tsv, so the 84 breakpoint sites are not recounted.
#   sbatch -A pi_matteop -p campus -a 1-99%60 -c 1 --mem=2G -t 1:00:00 -J pv3cnv cluster/panel_v3_cnv_depth.sh \
#       work_panel_v3/manifest.tsv 20 work_panel_v3/cnv_depth
# DIR holds sv_sites.depth_only.tsv and sv_evidence.py (the version with depth-only support);
# output DIR/sv/<pk>.sv.tsv, same columns as panel_v3_sv.sh. Kits whose output exists are skipped.
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"
cd "$D"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PATH="$D/envs/genomics/bin:$PATH"
MANIFEST="${1:?manifest}"; PER="${2:-20}"; O="${3:?work dir}"
mkdir -p "$O/sv"
start=$(( (SLURM_ARRAY_TASK_ID - 1) * PER + 1 )); end=$(( start + PER - 1 ))
ok=0; fail=0; skip=0
while IFS=$'\t' read -r pk barcode pub work; do
  [[ -n "$pk" ]] || continue
  out="$O/sv/$pk.sv.tsv"
  if [[ -s "$out" ]]; then skip=$((skip+1)); continue; fi
  bam="$work/markdup.bam"
  if [[ ! -s "$bam" ]]; then echo "PV3CNV-NOBAM $pk"; fail=$((fail+1)); continue; fi
  if python3 "$O/sv_evidence.py" count "$O/sv_sites.depth_only.tsv" "$pk" "$bam" > "$out.tmp" 2> "$out.err"; then
    mv "$out.tmp" "$out"; rm -f "$out.err"; ok=$((ok+1))
  else
    echo "PV3CNV-FAIL $pk $(tail -1 "$out.err")"; rm -f "$out.tmp"; fail=$((fail+1))
  fi
done < <(sed -n "${start},${end}p" "$MANIFEST")
echo "PV3CNV-DONE task=$SLURM_ARRAY_TASK_ID rows=$start-$end ok=$ok fail=$fail skip=$skip"
