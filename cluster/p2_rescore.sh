#!/bin/bash
# Rescore one published kit's panel v3 calls (stage 13c inputs + score, no apply) so its score files carry p_two,
# the chance of two copies. The 2026-11-13 re-imputation scored every kit this way and then deleted its work dir, so
# the live calls can only be matched by scoring again with the same frozen model, BAM and imputed BCF.
#
#   sbatch -A pi_mcdb -p pi_shared --qos=pi_shared24 -a 1-N%100 -c 1 --mem=4G -t 1:00:00 -J p2score cluster/p2_rescore.sh <list>
#
# Rows: name, barcode, pub dir, work dir (logs/nd1115_list.tsv). Output: work_panel_v3/p2score/<pub basename>/score.
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
export OPENBLAS_NUM_THREADS=1
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"; cd "$D"
export PATH="$D/envs/genomics/bin:$PATH" DATA_PYTHON="$D/envs/genomics/bin/python3" FASTA="$D/canFam4_plusY.fa"
row=$(sed -n "${SLURM_ARRAY_TASK_ID:?}p" "${1:?list}"); [ -n "$row" ] || exit 0
name=$(echo "$row" | cut -f1); PUB=$(echo "$row" | cut -f3); WORK0=$(echo "$row" | cut -f4); base=$(basename "$PUB")
A=/u/project/pellegrini_archive/data/dogs_bams
bam=""; for c in "$A/$name.markdup.bam" "$A/$base.markdup.bam" "$A/${base#pk-}.markdup.bam" "$WORK0/markdup.bam"; do [ -s "$c" ] && { bam=$c; break; }; done
[ -n "$bam" ] || { echo "P2-SKIP $name no BAM"; exit 0; }
O=$D/work_panel_v3/p2score/$base; rm -rf "$O"; mkdir -p "$O"
bash analysis/omia/panel_v3_kit.sh "$name" "$bam" "/u/project/pellegrini_archive/data/dogs_imputed_v2/$name.bcf" "$PUB" \
     reference_panel/panel_v3 "$O/kit" \
  && "$DATA_PYTHON" analysis/omia/panel_v3_breed_prior.py score reference_panel/panel_v3 "$O/kit" "$O/score" \
  && echo "P2-DONE $name" || { echo "P2-FAIL $name"; exit 1; }
