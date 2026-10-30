#!/bin/bash
# Breed panel v2.3: full-depth genotypes for the SRA reference genomes that the 50M-spot download
# (sra-download3.sh) left partial — coordinate-sorted submissions kept only the first chromosomes.
# One run per task: prefetch the whole run -> fasterq-dump -> bwa-mem2 -> markdup (as stages 3+4)
# -> keep only the genotypes, then delete the .sra, FASTQ and BAM:
#   sra_full/<RUN>.g.vcf.gz      genome-wide gVCF (chr1-38, X, Y, M), for later panels/traits/relatives
#   sra_full/<RUN>.panel.tsv.gz  breed-panel sites: CHROM POS REF ALT GT GQ AD DP (REF/ALT from canFam4;
#                                hom-ref rows carry no GQ, so judge them by DP)
#   sra_full/<RUN>.done          label + stats; written last, so its presence means the task finished
#   qsub -t 1-N -tc 20 cluster/sra-full-genotype.sh logs/sra_full_runs.txt     (RUN<TAB>LABEL per line)
#$ -cwd
#$ -j y
#$ -r y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=4G,h_rt=24:00:00
#$ -pe shared 8
set -uo pipefail
D="${SGE_O_WORKDIR:?run via qsub from the repo root}"; cd "$D"
IFS=$'\t' read -r run label < <(sed -n "${SGE_TASK_ID}p" "$1"); [ -n "${run:-}" ] || exit 0
export PATH=/u/home/g/gkislik/project-pellegrini/sratoolkit.3.1.1-ubuntu64/bin:$D/envs/genomics/bin:$PATH
NP=${NSLOTS:-8}; REF=canFam4_plusY.fa; OUT=work_sra_full/$run; ALLELES=logs/panel_alleles.tsv.gz
mkdir -p sra_full "$OUT"
[ -s "sra_full/$run.done" ] && { echo "SRAFULL-SKIP $run (done)"; exit 0; }
[ -s "$ALLELES.tbi" ] || { echo "SRAFULL-FAIL $run (no $ALLELES)"; exit 1; }
fail() { echo "SRAFULL-FAIL $run $label: $1"; exit 1; }
t0=$(date +%s); echo "[$(date +%T)] $run $label on $(hostname) with $NP slots"

# 1. Whole run (no spot cap), unpacked next to the .sra; reads unpaired after filtering are dropped.
if [ ! -s "$OUT/markdup.bam" ]; then
  prefetch --max-size 200G -O "$OUT" "$run" >/dev/null || fail prefetch
  fasterq-dump --split-3 -e "$NP" -t "$OUT/tmp" -O "$OUT" "$OUT/$run/$run.sra" || fail fasterq-dump
  rm -rf "$OUT/$run" "$OUT/tmp" "$OUT/$run.fastq"
  [ -s "$OUT/${run}_1.fastq" ] && [ -s "$OUT/${run}_2.fastq" ] || fail "no read pairs"
  echo "[$(date +%T)] reads: $(du -sh "$OUT" | cut -f1)"

  # 2. Align -> sort -> fixmate -> markdup, as run_dog_pipeline.sh stages 3+4 (smaller sort buffers: 4G/slot).
  bwa-mem2 mem -t "$NP" -R "@RG\tID:SRA-$run\tSM:SRA-$run\tPL:ILLUMINA\tLB:WGS" "$REF" \
    "$OUT/${run}_1.fastq" "$OUT/${run}_2.fastq" 2>"$OUT/bwa.log" \
  | samtools sort -n -m 400M -@ "$NP" -T "$OUT/tmp_sort" \
  | samtools fixmate -m -@ "$NP" - - \
  | samtools sort -m 400M -@ "$NP" -T "$OUT/tmp_sort2" \
  | samtools markdup -@ "$NP" --write-index - "$OUT/markdup.bam.tmp" \
  && mv "$OUT/markdup.bam.tmp" "$OUT/markdup.bam" && mv "$OUT/markdup.bam.tmp.csi" "$OUT/markdup.bam.csi" \
  || { tail -5 "$OUT/bwa.log"; fail align; }
  rm -f "$OUT/${run}_1.fastq" "$OUT/${run}_2.fastq"
  echo "[$(date +%T)] BAM: $(du -sh "$OUT/markdup.bam" | cut -f1)"
fi

# 3. Genome-wide gVCF, one chromosome per process.
contigs=$(awk '$1 ~ /^chr([0-9]+|X|Y|M)$/ {print $1}' "$REF.fai")
echo "$contigs" | xargs -P "$NP" -I{} sh -c \
  "bcftools mpileup -f $REF -r {} -q 30 -Q 20 -a AD,DP -g 5 -Ou $OUT/markdup.bam 2>/dev/null \
   | bcftools call -m -g 5 -Oz -o $OUT/g.{}.vcf.gz 2>/dev/null" || fail gvcf
bcftools concat -Oz -o "sra_full/$run.g.vcf.gz.tmp" $(for c in $contigs; do echo "$OUT/g.$c.vcf.gz"; done) 2>/dev/null \
  && mv "sra_full/$run.g.vcf.gz.tmp" "sra_full/$run.g.vcf.gz" && tabix -f -p vcf "sra_full/$run.g.vcf.gz" || fail concat

# 4. Breed-panel sites, constrained to the panel alleles.
bcftools mpileup -f "$REF" -T "$ALLELES" -q 30 -Q 20 -a AD,DP -Ou "$OUT/markdup.bam" 2>/dev/null \
  | bcftools call -m -C alleles -T "$ALLELES" -f GQ -Ou 2>/dev/null \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%GQ\t%AD\t%DP]\n' | gzip > "sra_full/$run.panel.tsv.gz.tmp" \
  && mv "sra_full/$run.panel.tsv.gz.tmp" "sra_full/$run.panel.tsv.gz" || fail panel

# 5. Check, record, delete the BAM.
n=$(zcat "sra_full/$run.panel.tsv.gz" | wc -l)
stats=$(zcat "sra_full/$run.panel.tsv.gz" | awk -F'\t' '{s+=$8; if($8>=8) c++; ch[$1]=1} END{n=0; for(k in ch) n++; printf "mean_dp=%.1f sites_dp8=%d chroms=%d", s/NR, c, n}')
(( n > 100000 )) || fail "only $n panel sites"
echo -e "$label\t$run\tpanel_sites=$n\t$stats\tgvcf=$(du -h "sra_full/$run.g.vcf.gz" | cut -f1)\tmin=$(( ($(date +%s)-t0)/60 ))" > "sra_full/$run.done"
rm -rf "$OUT"
echo "SRAFULL-DONE $run $label panel_sites=$n $stats min=$(( ($(date +%s)-t0)/60 ))"
