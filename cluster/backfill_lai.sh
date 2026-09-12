#$ -cwd
#$ -j y
#$ -o logs/lai9c.$TASK_ID.log
#$ -l h_data=4G,h_rt=0:30:00
# Stage-9c backfill: chromosome painting + two-step breed_composition for
# every published report. Rows from cluster/repass9b_list.txt (name, barcode,
# pub dir, work dir). Re-queries the dosages exactly as stage 9b does, runs
# analysis/local_ancestry/lai_hmm.py, writes local_ancestry.json and
# rewrites breed_result.json (original kept once as breed_result.pre-lai.json).
cd $SGE_O_WORKDIR
D=$PWD
ENV=$D/envs/genomics
export LD_LIBRARY_PATH=$ENV/lib:${LD_LIBRARY_PATH:-}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONIOENCODING=utf-8

row=$(sed -n "${SGE_TASK_ID}p" cluster/repass9b_list.txt)
[ -n "$row" ] || exit 0
pk=$(echo "$row" | cut -f1); PUB=$(echo "$row" | cut -f3); WORK=$(echo "$row" | cut -f4)
IMPUTED_BCF=$(ls "$WORK"/glimpse2/*_imputed_dog10k.bcf 2>/dev/null | head -1)
[ -f "$IMPUTED_BCF" ] && [ -f "$PUB/breed_result.json" ] || { echo "LAI9C-SKIP $pk (no bcf/breed_result)"; exit 0; }

POS=$WORK/relatives_sites.pos; DS=$WORK/relatives_ds.tsv
awk 'NR>1{print $1"\t"$2}' "$D/breed_panel/sites.tsv" > "$POS"
$ENV/bin/bcftools query -R "$POS" -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' "$IMPUTED_BCF" > "$DS" \
  || { echo "LAI9C-FAIL $pk (bcftools)"; exit 1; }
[ -f "$PUB/breed_result.pre-lai.json" ] || cp "$PUB/breed_result.json" "$PUB/breed_result.pre-lai.json"
# always start from the lasso composition so a re-run is idempotent
cp "$PUB/breed_result.pre-lai.json" "$PUB/breed_result.json"
BREED_PANEL=$D/breed_panel $ENV/bin/python3 analysis/local_ancestry/lai_hmm.py \
    "$DS" "$PUB/breed_result.json" "$WORK/local_ancestry" --ds --panel "$D/breed_panel" --write-breed-result \
  && cp "$WORK/local_ancestry.lai.json" "$PUB/local_ancestry.json" \
  || { cp "$PUB/breed_result.pre-lai.json" "$PUB/breed_result.json"; echo "LAI9C-FAIL $pk (hmm)"; exit 1; }
rm -f "$DS" "$POS"
echo "LAI9C-DONE $pk"
