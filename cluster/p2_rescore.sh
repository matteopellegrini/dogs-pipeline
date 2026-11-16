#!/bin/bash
# Rescore one published kit's panel v3 calls (stage 13c inputs + score, no apply) so its score files carry p_two,
# the chance of two copies. The 2026-11-13 re-imputation scored every kit this way and then deleted its work dir, so
# the live calls can only be matched by scoring again with the same frozen model, BAM and imputed BCF.
#
#   CHUNK=20 sbatch -A pi_matteop -p campus -a 1-<ceil(N/20)> -c 1 --mem=4G -t 2:00:00 -J p2score cluster/p2_rescore.sh <list>
#
# Rows: name, barcode, pub dir, work dir (logs/nd1115_list.tsv); task t takes rows (t-1)*CHUNK+1 .. t*CHUNK (~1 min
# a kit; one task per kit hits campus's 500-job submit cap). Output: work_panel_v3/p2score/<pub basename>/score.
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
export OPENBLAS_NUM_THREADS=1
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"; cd "$D"
export PATH="$D/envs/genomics/bin:$PATH" DATA_PYTHON="$D/envs/genomics/bin/python3" FASTA="$D/canFam4_plusY.fa"
C=${CHUNK:-1}; fail=0
for r in $(seq $(( (${SLURM_ARRAY_TASK_ID:?} - 1) * C + 1 )) $(( SLURM_ARRAY_TASK_ID * C ))); do
row=$(sed -n "${r}p" "${1:?list}"); [ -n "$row" ] || break
name=$(echo "$row" | cut -f1); PUB=$(echo "$row" | cut -f3); WORK0=$(echo "$row" | cut -f4); base=$(basename "$PUB")
A=/u/project/pellegrini_archive/data/dogs_bams
bam=""; for c in "$A/$name.markdup.bam" "$A/$base.markdup.bam" "$A/${base#pk-}.markdup.bam" "$WORK0/markdup.bam"; do [ -s "$c" ] && { bam=$c; break; }; done
[ -n "$bam" ] || { echo "P2-SKIP $name no BAM"; continue; }
O=$D/work_panel_v3/p2score/$base; [ -s "$O/score/bp_p_carrier.npz" ] && { echo "P2-SKIP $name done"; continue; }
rm -rf "$O"; mkdir -p "$O"
bash analysis/omia/panel_v3_kit.sh "$name" "$bam" "/u/project/pellegrini_archive/data/dogs_imputed_v2/$name.bcf" "$PUB" \
     reference_panel/panel_v3 "$O/kit" \
  && "$DATA_PYTHON" analysis/omia/panel_v3_breed_prior.py score reference_panel/panel_v3 "$O/kit" "$O/score" \
  && echo "P2-DONE $name" || { echo "P2-FAIL $name"; fail=1; }
done
exit $fail
