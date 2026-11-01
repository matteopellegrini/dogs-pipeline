#!/bin/bash
# Breed panel v2.3: full-depth genotypes for the SRA reference genomes that the 50M-spot download
# (sra-download3.sh) left partial — coordinate-sorted submissions kept only the first chromosomes.
# Since 2026-10-30 every SRA reference genome is processed this way (Matteo: full files, no caps).
# One run per task: prefetch the whole run -> stream into bwa-mem2 -> markdup (as stages 3+4)
# -> keep only the genotypes, then delete the .sra and BAM:
#   sra_full/<RUN>.g.vcf.gz      genome-wide gVCF (chr1-38, X, Y, M), for later panels/traits/relatives
#   sra_full/<RUN>.panel.tsv.gz  breed-panel sites: CHROM POS REF ALT GT GQ AD DP (REF/ALT from canFam4;
#                                hom-ref rows carry no GQ, so judge them by DP)
#   sra_full/<RUN>.done          label + stats; written last, so its presence means the task finished
#   qsub -t 1-N -tc 20 cluster/sra-full-genotype.sh logs/sra_full_runs.txt     (RUN<TAB>LABEL per line)
#   highp copy alongside it: qsub -N sra_full_hp -t 1-N -tc 4 -l h_data=4G,h_rt=48:00:00,highp ...
#   Slurm (login.hoffman2.oarc.ucla.edu): sbatch -J sra_full -a 3-N%20 -p campus -c 16 --mem=60G -t 1-00:00:00 \
#     cluster/sra-full-genotype.sh logs/sra_full_runs_smallfirst.txt     (MCDB nodes: --qos=highp -A pi_mcdb)
#$ -cwd
#$ -j y
#$ -r y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=4G,h_rt=24:00:00
#$ -pe shared 8
#SBATCH -o logs/%x.%A.%a.log
#SBATCH --requeue
set -uo pipefail
if [ -n "${SLURM_JOB_ID:-}" ]; then   # Slurm: map its variables onto the SGE names used below
  JOB_ID=${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}; SGE_TASK_ID=$SLURM_ARRAY_TASK_ID
  SGE_O_WORKDIR=$SLURM_SUBMIT_DIR; NSLOTS=$SLURM_CPUS_PER_TASK
fi
D="${SGE_O_WORKDIR:?run via qsub from the repo root}"; cd "$D"
IFS=$'\t' read -r run label < <(sed -n "${SGE_TASK_ID}p" "$1"); [ -n "${run:-}" ] || exit 0
export PATH=/u/home/g/gkislik/project-pellegrini/sratoolkit.3.1.1-ubuntu64/bin:$D/envs/genomics/bin:$PATH
NP=${NSLOTS:-8}; REF=canFam4_plusY.fa; OUT=work_sra_full/$run; ALLELES=logs/panel_alleles.tsv.gz
mkdir -p sra_full "$OUT"
[ -s "sra_full/$run.done" ] && { echo "SRAFULL-SKIP $run (done)"; exit 0; }
[ -s "$ALLELES.tbi" ] || { echo "SRAFULL-FAIL $run (no $ALLELES)"; exit 1; }
# Two arrays (general pool + highp) walk the same list: a lock directory claims each run. The holder touches it
# every 10 min; a lock untouched for 30 min was left by a killed task and is taken over.
LOCK=work_sra_full/$run.lock
if ! mkdir "$LOCK" 2>/dev/null; then
  if [ -n "$(find "$LOCK" -maxdepth 0 -mmin -30)" ]; then echo "SRAFULL-BUSY $run (held by $(ls "$LOCK"))"; exit 0; fi
  echo "taking over stale lock from $(ls "$LOCK")"; rm -rf "$LOCK"; mkdir "$LOCK" || { echo "SRAFULL-BUSY $run"; exit 0; }
fi
touch "$LOCK/$JOB_ID.$SGE_TASK_ID"
( while sleep 600; do touch "$LOCK"; done ) & HB=$!
trap 'kill $HB 2>/dev/null; rm -rf "$LOCK"' EXIT
fail() { echo "SRAFULL-FAIL $run $label: $1"; exit 1; }
t0=$(date +%s); echo "[$(date +%T)] $run $label on $(hostname) with $NP slots"

# 1+2. Whole run (no spot cap), streamed: fasterq-dump interleaves the mates, awk drops spots whose mate was
# filtered out, bwa-mem2 -p aligns pairs -> sort -> fixmate -> markdup as run_dog_pipeline.sh stages 3+4.
# No FASTQ on disk (a 59 GB .sra unpacks to ~300 GB); peak disk = .sra + BAM. Sort buffers sized for 4G/slot.
if [ ! -s "$OUT/markdup.bam" ]; then
  rm -rf "$OUT"/tmp_sort* "$OUT"/markdup.bam.tmp*   # a killed attempt's leftovers
  prefetch --max-size 200G -O "$OUT" "$run" >/dev/null || fail prefetch
  echo "[$(date +%T)] .sra: $(du -sh "$OUT/$run" | cut -f1)"
  fasterq-dump --split-spot --stdout -e 2 -t "$OUT/tmp" --seq-defline '@$ac.$si/$ri' --qual-defline '+' \
      "$OUT/$run/$run.sra" 2>"$OUT/fqd.log" \
  | awk 'NR%4==1{k=$1; sub(/\/[12]$/,"",k)} {rec[NR%4]=$0}
         NR%4==0{r=rec[1]"\n"rec[2]"\n"rec[3]"\n"rec[0]; if(k==pk){print pr; print r; n++; pk=""} else {pr=r; pk=k}}
         END{print n+0 > "'"$OUT"'/pairs.txt"}' \
  | bwa-mem2 mem -p -t "$NP" -R "@RG\tID:SRA-$run\tSM:SRA-$run\tPL:ILLUMINA\tLB:WGS" "$REF" - 2>"$OUT/bwa.log" \
  | samtools sort -n -m 400M -@ "$NP" -T "$OUT/tmp_sort" \
  | samtools fixmate -m -@ "$NP" - - \
  | samtools sort -m 400M -@ "$NP" -T "$OUT/tmp_sort2" \
  | samtools markdup -@ "$NP" --write-index - "$OUT/markdup.bam.tmp" \
  && [ "$(cat "$OUT/pairs.txt" 2>/dev/null || echo 0)" -gt 1000000 ] \
  && mv "$OUT/markdup.bam.tmp" "$OUT/markdup.bam" && mv "$OUT/markdup.bam.tmp.csi" "$OUT/markdup.bam.csi" \
  || { tail -3 "$OUT/fqd.log" "$OUT/bwa.log"; fail "align (pairs=$(cat "$OUT/pairs.txt" 2>/dev/null))"; }
  rm -rf "$OUT/$run" "$OUT/tmp"
  echo "[$(date +%T)] pairs: $(cat "$OUT/pairs.txt")  BAM: $(du -sh "$OUT/markdup.bam" | cut -f1)"
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
