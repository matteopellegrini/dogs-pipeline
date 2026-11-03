#!/bin/bash
# Panel v3: allele counts at the mitochondrial catalogue sites, plus each kit's inferred sex (for X ploidy).
# chrM catalogue sites sit in NUMT-shadowed stretches where every read in markdup.bam has MAPQ 0 (e.g. chrM:1900),
# so the MAPQ filter of panel_v3_genotype.sh drops them all; mtDNA outnumbers the nuclear copies ~100:1, so here
# reads are counted with -q 0. Output: work_panel_v3/chrM_counts.tsv (kit, chrom, pos, ref, alt, ref_reads,
# alt_reads) and work_panel_v3/kit_sex.tsv (kit, M/F from cnv50_result.json).
#   sbatch -A pi_matteop -p campus -c 1 --mem=4G -t 4:00:00 -J pv3chrM cluster/panel_v3_chrM.sh work_panel_v3/manifest.tsv
#SBATCH -o logs/%x.%j.log
set -uo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"
cd "$D"
export OPENBLAS_NUM_THREADS=1
PATH="$D/envs/genomics/bin:$PATH"
O="$D/work_panel_v3"
zcat "$O/v3_small.alleles.tsv.gz" | awk '$1=="chrM"' > "$O/chrM.alleles.tsv"
: > "$O/chrM_counts.tsv"; : > "$O/kit_sex.tsv"
while IFS=$'\t' read -r pk barcode pub work; do
  s=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1])).get('sex_inferred',''))" "$pub/cnv50_result.json" 2>/dev/null)
  printf '%s\t%s\n' "$pk" "${s:-NA}" >> "$O/kit_sex.tsv"
  bam="$work/markdup.bam"; [[ -s "$bam" ]] || continue
  while IFS=$'\t' read -r c p al; do
    ref=${al%,*}; alt=${al#*,}
    bcftools mpileup -f "$D/canFam4_plusY.fa" -r "$c:$p-$p" -q 0 -Q 20 -B -m 1 -d 100000 -a FORMAT/AD "$bam" 2>/dev/null \
      | bcftools query -f '%REF\t%ALT\t[%AD]\n' 2>/dev/null \
      | awk -v k="$pk" -v c="$c" -v p="$p" -v r="$ref" -v a="$alt" -F'\t' '
          $1==r { n=split($2,al,","); split($3,ad,","); x=0; for(i=1;i<=n;i++) if(al[i]==a) x=ad[i+1]; print k"\t"c"\t"p"\t"r"\t"a"\t"ad[1]"\t"x; found=1 }
          END { if(!found) print k"\t"c"\t"p"\t"r"\t"a"\t0\t0" }' >> "$O/chrM_counts.tsv"
  done < "$O/chrM.alleles.tsv"
done < "${1:?manifest}"
echo "PV3-CHRM-DONE $(wc -l < "$O/kit_sex.tsv") kits"
