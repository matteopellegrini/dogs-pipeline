#!/bin/bash
# Rebuild functional_variants.json on published kits from a rerun list (name, barcode, pub dir, work dir — the same
# lists as rerun-13-list.sh) after the Dog10K panel annotation changed (reference_panel/dog10k_ann/dog10k_sites.ann.bcf).
# Only the two stage-10 steps that depend on it: copy ANN/LOF/NMD onto the kit's imputed BCF, then
# analysis/snpeff/functional_variants.py. No imputation, and cnv_genes.json / cnv_homdel.json are left alone.
# Keeps functional_variants.pre-<tag>.json once per tag. Nothing is uploaded; republish with
#   FILES="functional_variants.json" bash cluster/republish_files.sh <tag> <list>
#   sbatch -A pi_matteop -J rerun10fv -a 1-N%100 -p campus -c 1 --mem=6G -t 1:00:00 cluster/rerun-10-fv-list.sh <list> <tag>
# CHUNK=k makes each array task do k consecutive rows (the user's job submit limit is a few hundred): -a 1-ceil(N/k)
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"
LIST="${1:?usage: sbatch -a 1-N cluster/rerun-10-fv-list.sh <list> <tag>}"
TAG="${2:?tag required}"
cd "$D"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONIOENCODING=utf-8
E="$D/envs/genomics"; export PATH="$E/bin:$PATH"
SNPEFF_DB="${SNPEFF_DB:-canFam4_NCBI_AR106}"
ANN_REF="${SNPEFF_ANN_REF:-$D/reference_panel/dog10k_ann/dog10k_sites.ann.bcf}"
CHUNK="${CHUNK:-1}"
rc=0
for ((r = (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1; r <= SLURM_ARRAY_TASK_ID * CHUNK; r++)); do
row=$(sed -n "${r}p" "$LIST"); [ -n "$row" ] || break
pk=$(echo "$row" | cut -f1); PUB=$(echo "$row" | cut -f3); WORK=$(echo "$row" | cut -f4)
bcf=$(compgen -G "$WORK/glimpse2/*_imputed_dog10k.bcf" | head -1)
[ -n "$bcf" ] || { echo "RERUN10FV-SKIP $pk (no bcf)"; continue; }
[ -f "$PUB/functional_variants.pre-$TAG.json" ] || cp "$PUB/functional_variants.json" "$PUB/functional_variants.pre-$TAG.json" 2>/dev/null || true
tmp=$(mktemp -d "${TMPDIR:-/tmp}/fv.XXXXXX")
if ! bcftools annotate -a "$ANN_REF" -c INFO/ANN,INFO/LOF,INFO/NMD -Oz -o "$tmp/annotated.vcf.gz" "$bcf"; then
  echo "RERUN10FV-FAIL $pk (annotate)"; rc=1
elif ! "$E/bin/python3" analysis/snpeff/functional_variants.py "$tmp/annotated.vcf.gz" "$bcf" "$PUB" "$SNPEFF_DB"; then
  cp "$PUB/functional_variants.pre-$TAG.json" "$PUB/functional_variants.json" 2>/dev/null; echo "RERUN10FV-FAIL $pk (parse)"; rc=1
else
  echo "RERUN10FV-DONE $pk"
fi
rm -rf "$tmp"
done
exit $rc
