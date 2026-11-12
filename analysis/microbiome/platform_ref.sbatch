#!/bin/bash
# New-process microbiome reference: one swab per array task, profiled the way the old-lab cohort was
# (unmapped reads capped at CAP with the pipeline's fixed seed, MetaPhlAn4, HUMAnN 3.9 EC CPM, scored
# against reference_panel/oral_function_panel.json). The cap matters: richness and, mildly, pathobiont
# burden depend on how many bacterial reads a sample has, and the new-process swabs carry 12-52M against
# ~1-1.5M for a typical cohort kit. fit_platform_offset.py turns the outputs into the offset file.
#   sbatch -A pi_mcdb -p pi_shared --qos=pi_shared24 -a 1-N -c 8 --mem=64G -t 12:00:00 \
#          analysis/microbiome/platform_ref.sbatch analysis/microbiome/platform_ref_samples.tsv
# Rows: sample (BAM name in the archive), dog id, kit/platform label.
#SBATCH -J micro_platform_ref
#SBATCH -o logs/%x.%A.%a.log
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONIOENCODING=utf-8
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"; cd "$D"
source site/hoffman.sh
CAP=${CAP:-2000000}; NPROC=${SLURM_CPUS_PER_TASK:-8}
row=$(grep -v '^#' "${1:?list}" | sed -n "${SLURM_ARRAY_TASK_ID:?}p"); [ -n "$row" ] || exit 0
s=$(echo "$row" | cut -f1); lc=$(echo "$s" | tr A-Z a-z)
O=$D/scratch/micro_platform/$lc; mkdir -p "$O"; cd "$O"
fq=$(ls "$D"/work/{$s,$lc}/analysis/*_unmapped.fastq 2>/dev/null | head -1)
if [ -z "$fq" ]; then
  "$ENV_GENOMICS/bin/samtools" fastq -f 4 -@ "$NPROC" "$BAM_ARCHIVE_DIR/$s.markdup.bam" > all_unmapped.fastq 2>/dev/null
  fq=$O/all_unmapped.fastq
fi
n=$(( $(wc -l < "$fq") / 4 ))
keep=$(awk -v n="$n" -v cap="$CAP" 'BEGIN{p=cap/n; if(p>1)p=1; printf "%.6f", p}')
awk -v p="$keep" 'BEGIN{srand(11)} NR%4==1{k=(rand()<p)} k' "$fq" > capped.fastq
rm -f all_unmapped.fastq
echo "$s unmapped=$n keep=$keep capped=$(( $(wc -l < capped.fastq) / 4 ))"
export PATH="$(dirname "$METAPHLAN_BIN"):$PATH"
"$METAPHLAN_BIN" capped.fastq --db_dir "$METAPHLAN_DB" --index "$METAPHLAN_INDEX" --input_type fastq \
    --mapout "${lc}.mapout.bz2" --nproc "$NPROC" -o "${lc}_metaphlan.txt" 2> metaphlan.log || { echo "PREF-FAIL $s metaphlan"; exit 1; }
reads=$(grep -m1 -o '^#[0-9]* reads processed' "${lc}_metaphlan.txt" | tr -dc '0-9')
H=$D/envs/humann; HDB=$D/humann_db
sed 's/mpa_vJan25_CHOCOPhlAnSGB_202503/mpa_vJun23_CHOCOPhlAnSGB_202403/' "${lc}_metaphlan.txt" > profile_for_humann.tsv
if ( export PATH="$H/bin:$ENV_GENOMICS/bin:$PATH"; \
     humann --input capped.fastq --output humann --threads "$NPROC" --nucleotide-database "$HDB/chocophlan" \
       --protein-database "$HDB/uniref" --output-basename "$lc" --remove-temp-output --taxonomic-profile profile_for_humann.tsv > humann.log 2>&1 \
     && humann_regroup_table -i "humann/${lc}_genefamilies.tsv" -g uniref90_level4ec -o "humann/${lc}_ec.tsv" >> humann.log 2>&1 \
     && humann_renorm_table -i "humann/${lc}_ec.tsv" -u cpm -o "humann/${lc}_ec_cpm.tsv" >> humann.log 2>&1 ); then
  python3 "$D/analysis/oral_function/score_oral_function.py" --ec "humann/${lc}_ec_cpm.tsv" --profile "${lc}_metaphlan.txt" \
      --reads "$reads" --panel "$D/reference_panel/oral_function_panel.json" --out oral_function_result.json
  gzip -f "humann/${lc}_genefamilies.tsv"; rm -rf humann/*_humann_temp_*
else
  echo "PREF-WARN $s humann failed"; tail -5 humann.log
fi
rm -f capped.fastq
echo "PREF-DONE $s reads=$reads"
