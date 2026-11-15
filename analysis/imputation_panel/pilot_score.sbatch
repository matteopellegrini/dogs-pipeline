#!/bin/bash
# Imputation-panel v3 pilot: score every pilot imputation of one chromosome (pilot_impute.sbatch) against truth.
# Cosmo: COSMO-SUPER read truth (score_panels.py defaults). SRA test dogs: their full-depth gVCF (sra_truth.py).
# Writes $T/pilot/score/<sample>_<chrom>.tsv; summarise with
#   python3 analysis/imputation_panel/summarise.py $T/pilot/score/<sample>_chr*.tsv
#   sbatch -A pi_matteop -p campus -a 18,27,38 -c 1 --mem=8G -t 6:00:00 analysis/imputation_panel/pilot_score.sbatch
#SBATCH -J pilot_score
#SBATCH -o logs/pilot_score.%A.%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
D=/u/project/pellegrini/gkislik/dogs; cd $D
T=$D/scratch/panel_test; P=$T/pilot; PY=$D/envs/genomics/bin/python3
c=chr${SLURM_ARRAY_TASK_ID:?}
mkdir -p $P/score $P/truth
for sample in $(awk -v c=$c '$1==c{print $2}' $P/jobs.tsv | sort -u); do
  runs=()
  while read -r cc s depth panel; do
    f=$P/imp/${c}_${s}_${depth}_${panel}.bcf
    [[ $cc == $c && $s == $sample && -s $f.csi ]] && runs+=("${depth}_${panel}=$f")
  done < $P/jobs.tsv
  [[ ${#runs[@]} -gt 0 ]] || continue
  if [[ $sample == COSMO ]]; then
    $PY analysis/imputation_panel/score_panels.py $c $P/score/${sample}_$c.tsv "${runs[@]}"
  else
    tr=$P/truth/${sample}_$c.pileup.tsv.gz
    [[ -s $tr ]] || $PY analysis/imputation_panel/sra_truth.py sra_full/$sample.g.vcf.gz $c $T/ref/raf_$c.tsv.gz $tr
    TRUTH_PILEUP=$tr TRUTH_SAMPLES=SRA TRUTH_DROP=NONE \
      $PY analysis/imputation_panel/score_panels.py $c $P/score/${sample}_$c.tsv "${runs[@]}"
  fi
done
echo "SCORE-DONE $c"
