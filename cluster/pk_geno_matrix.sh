#!/bin/bash
# Build one multi-sample genotype matrix (Dog10K+DAP imputation-panel sites,
# 29.09M SNVs/indels on chr1-38) for every ProsperKits dog with an imputed
# genome in dogs_imputed_v2. Prepares inputs and submits the three stages:
# pk_geno_merge (chromosome x kit group) -> pk_geno_chr -> pk_geno_final.
# Run from $D on hoffman-slurm.
set -euo pipefail
IMP=/u/project/pellegrini_archive/data/dogs_imputed_v2
OUT=${OUT:-/u/project/pellegrini_archive/data/prosperkits/genotypes_dog10k}
mkdir -p "$OUT/parts" logs
ls "$IMP" | grep -E '^pk-[0-9]+\.bcf$' | sort | sed "s|^|$IMP/|" > "$OUT/parts/files.txt"
sed 's|.*/||; s|\.bcf$||' "$OUT/parts/files.txt" > "$OUT/samples.txt"
split -d -l 200 "$OUT/parts/files.txt" "$OUT/parts/files.g"
split -d -l 200 "$OUT/samples.txt" "$OUT/parts/samples.g"
[ "$(ls "$OUT"/parts/files.g* | wc -l)" = 9 ] || { echo "expected 9 groups" >&2; exit 1; }
echo "kits: $(wc -l < "$OUT/samples.txt")"
[ "${DRY:-0}" = 1 ] && exit 0
a=$(sbatch --parsable --export=ALL,OUT="$OUT" -o "logs/pkgeno_%A_%a.out" cluster/pk_geno_merge.sbatch)
b=$(sbatch --parsable --export=ALL,OUT="$OUT" -d afterok:"${a%%;*}" -o "logs/pkgenochr_%A_%a.out" cluster/pk_geno_chr.sbatch)
f=$(sbatch --parsable --export=ALL,OUT="$OUT" -d afterok:"${b%%;*}" -o "logs/pkgenofin_%j.out" cluster/pk_geno_final.sbatch)
echo "groups $a, chromosomes $b, final $f"
