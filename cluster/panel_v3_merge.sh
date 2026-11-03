#!/bin/bash
# Panel v3, step 2a: merge the per-kit outputs of cluster/panel_v3_genotype.sh into cohort files and
# run the carrier scan (analysis/omia/panel_v3_carriers.py).
#   sbatch -A pi_matteop -p campus -c 2 --mem=16G -t 2:00:00 -J pv3merge cluster/panel_v3_merge.sh
#SBATCH -o logs/%x.%j.log
set -euo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"
cd "$D"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PATH="$D/envs/genomics/bin:$PATH"
O="$D/work_panel_v3"
ulimit -n 8192 2>/dev/null || true
ls "$O"/kits/*.reads.bcf > "$O/reads.list"; ls "$O"/kits/*.imp.bcf > "$O/imp.list"
bcftools merge -m none --file-list "$O/reads.list" -Ob -o "$O/cohort_reads.bcf" --threads 2
bcftools index -f "$O/cohort_reads.bcf"
bcftools merge -m none --file-list "$O/imp.list" -Ob -o "$O/cohort_imputed.bcf" --threads 2
bcftools index -f "$O/cohort_imputed.bcf"
{ printf 'CHROM\tPOS\tREF\tALT\t'; bcftools query -l "$O/cohort_reads.bcf" | paste -sd '\t'; } > "$O/cohort_reads.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$O/cohort_reads.bcf" >> "$O/cohort_reads.tsv"
gzip -f "$O/cohort_reads.tsv"
python3 analysis/omia/panel_v3_carriers.py "$O/cohort_reads.tsv.gz" "$O/manifest.tsv" \
  known_variants/panel_v3/panel_v3_sites.tsv "$O/v3_carriers"
echo "PV3-MERGE-DONE $(bcftools query -l "$O/cohort_reads.bcf" | wc -l) kits"
