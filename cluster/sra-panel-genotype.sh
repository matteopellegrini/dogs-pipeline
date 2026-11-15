#!/bin/bash
# Imputation-panel v3 pilot: stream a public dog genome from SRA and keep only its genotypes.
# One BioSample per task (all its WGS runs): prefetch one run at a time -> fasterq-dump -> random read sample
# -> bwa-mem2 -> markdup (as run_dog_pipeline.sh stages 3+4) -> genotypes -> delete the .sra and BAM.
#   role=panel: reads sampled down to ~CAP_X (default 10x) by a random pair filter, so coordinate-sorted
#               submissions are capped evenly (taking the first N spots would keep only the first chromosomes).
#               Keeps sra_panel/<BS>.panel_sites.bcf  GT:DP:GQ at every panel v2 site (ref_panel_v2/sites/<chr>.tsv.gz,
#               alleles forced to the panel's), sra_panel/<BS>.panel.tsv.gz (breed-panel sites, as sra-full-genotype.sh)
#               and sra_slices/<BS>.bam (candidate-gene reads).
#   role=test:  full depth, no cap; keeps sra_full/<BS>.g.vcf.gz (truth for scoring) and a ~1x random subsample
#               scratch/panel_test/sra/<BS>/lowpass.bam for imputation tests. Test dogs never enter a panel.
# Rows (analysis/imputation_panel/pilot_v3_dogs.tsv): biosample group breed role depth_x bases runs(comma)
#   sbatch -J sra_panel -a 1-N%40 -A pi_matteop -p campus -c 8 --mem=64G -t 1-00:00:00 cluster/sra-panel-genotype.sh <list>
#   (48G: ~23% of tasks died at the limit in the pipe bwa-mem2 | sort -n | fixmate | sort | markdup, 2026-11-15)
# A second array (e.g. pi_shared) can walk the same list: a lock directory claims each BioSample.
# FAST=1 (panel role): when the first run's reads are in random genomic order (first 2,000 spots land on >= 20
#   chromosomes), read only the first ~CAP_X worth of spots of each run straight from SRA (fastq-dump -X, no prefetch)
#   instead of downloading and reading everything; position-sorted submissions fall back to the full read.
#SBATCH -o logs/%x.%A.%a.log
#SBATCH --requeue
set -uo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"; cd "$D"
row=$(grep -v '^#' "${1:?list}" | sed -n "${SLURM_ARRAY_TASK_ID:?}p"); [ -n "$row" ] || exit 0
IFS=$'\t' read -r bs group breed role depth bases runs <<< "$row"
export PATH=/u/home/g/gkislik/project-pellegrini/sratoolkit.3.1.1-ubuntu64/bin:$D/envs/genomics/bin:$PATH
NP=${SLURM_CPUS_PER_TASK:-8}; REF=canFam4_plusY.fa; OUT=work_sra_panel/$bs; SITES=ref_panel_v2/sites
CAP_X=${CAP_X:-10}; ALLELES=logs/panel_alleles.tsv.gz
mkdir -p sra_panel sra_full sra_slices "$OUT"
DONE=sra_panel/$bs.done
[ -s "$DONE" ] && { echo "SRAPANEL-SKIP $bs (done)"; exit 0; }
[ -s "$SITES/chr38.tsv.gz.tbi" ] || { echo "SRAPANEL-FAIL $bs (no $SITES)"; exit 1; }
LOCK=work_sra_panel/$bs.lock
if ! mkdir "$LOCK" 2>/dev/null; then
  if [ -n "$(find "$LOCK" -maxdepth 0 -mmin -30)" ]; then echo "SRAPANEL-BUSY $bs"; exit 0; fi
  rm -rf "$LOCK"; mkdir "$LOCK" || { echo "SRAPANEL-BUSY $bs"; exit 0; }
fi
touch "$LOCK/${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}.$SLURM_ARRAY_TASK_ID"
( while sleep 600; do touch "$LOCK"; done ) & HB=$!
trap 'kill $HB 2>/dev/null; rm -rf "$LOCK"' EXIT
# A failed task drops its downloads and partial BAMs (a failed task kept ~30 GB each; 163 of them filled 4.7 TB),
# keeping only a finished markdup.bam, which a rerun resumes from.
fail() { echo "SRAPANEL-FAIL $bs $group: $1"; rm -rf "$OUT/sra" "$OUT/tmp" "$OUT"/tmp_sort* "$OUT"/markdup.bam.tmp*; exit 1; }
t0=$(date +%s); echo "[$(date +%T)] $bs $group ($breed) role=$role ${depth}x runs=$runs on $(hostname)"

# Fraction of read pairs to keep: CAP_X over the raw depth (bases / 2.48 Gb); test dogs keep everything.
keep=1
[ "$role" = panel ] && keep=$(awk -v b="$bases" -v c="$CAP_X" 'BEGIN{k=c*2.48e9/b; if(k>1)k=1; printf "%.5f", k}')

# Paired runs only when the sample has any (a few PAIRED submissions store each mate as its own spot: align those
# single-end, as sra-full-genotype.sh does).
pe=0; use=()
for r in ${runs//,/ }; do
  p=$(fastq-dump -X 1000 --split-spot -Z "$r" 2>/dev/null | awk 'END{print (NR/4 > 1500) ? 1 : 0}')
  [ "$p" = 1 ] && pe=1
  use+=("$r:$p")
done
sel=(); for x in "${use[@]}"; do [ "$pe" = 0 ] || [ "${x#*:}" = 1 ] && sel+=("${x%%:*}"); done
echo "[$(date +%T)] keep=$keep $([ "$pe" = 1 ] && echo paired || echo single-end) runs used: ${sel[*]}"

# FAST: spots needed for CAP_X (bases per spot from the first 2,000 spots), if those spots are spread over the genome
nspots=0
if [ "${FAST:-0}" = 1 ] && [ "$role" = panel ] && [ "$keep" != 1.00000 ]; then
  fastq-dump -X 2000 --split-spot -Z "${sel[0]}" 2>/dev/null > "$OUT/head.fq"
  bps=$(awk -v pe="$pe" 'NR%4==2{b+=length($0)} END{print (NR ? b/(NR/4/(pe?2:1)) : 0)}' "$OUT/head.fq")
  nchr=$(bwa-mem2 mem $([ "$pe" = 1 ] && echo -p) -t "$NP" "$REF" "$OUT/head.fq" 2>/dev/null \
         | samtools view -F 0x904 -q 20 - | cut -f3 | grep -E '^chr[0-9]+$' | sort -u | wc -l)
  rm -f "$OUT/head.fq"
  if [ "${nchr:-0}" -ge 20 ] && awk -v b="$bps" 'BEGIN{exit !(b>50)}'; then
    nspots=$(awk -v c="$CAP_X" -v b="$bps" 'BEGIN{printf "%d", c*2.48e9/b}')
    keep=1
  fi
  echo "[$(date +%T)] FAST: first spots on $nchr autosomes, $bps bases/spot -> $([ "$nspots" -gt 0 ] && echo "read first $nspots spots" || echo "sorted or unknown: full read")"
fi

if [ ! -s "$OUT/markdup.bam" ]; then
  rm -rf "$OUT"/tmp_sort* "$OUT"/markdup.bam.tmp* "$OUT"/sra
  left=$nspots
  for r in "${sel[@]}"; do
    if [ "$nspots" -gt 0 ]; then
      [ "$left" -gt 0 ] || break
      fastq-dump --split-spot -Z -X "$left" --defline-seq '@$ac.$si/$ri' --defline-qual '+' "$r" 2>>"$OUT/fqd.log" \
        | awk -v f="$OUT/n.$r" '{print} END{print NR/4 > f}' || { echo "fastq-dump $r failed" >&2; exit 4; }
      left=$(( left - $(awk -v pe="$pe" '{printf "%d", $1/(pe?2:1)}' "$OUT/n.$r") ))
      continue
    fi
    prefetch --max-size 400G -O "$OUT/sra" "$r" >/dev/null 2>>"$OUT/prefetch.log" || { echo "prefetch $r failed" >&2; exit 3; }
    fasterq-dump --split-spot --stdout -e 2 -t "$OUT/tmp" --seq-defline '@$ac.$si/$ri' --qual-defline '+' \
        "$OUT/sra/$r/$r.sra" 2>>"$OUT/fqd.log" || { echo "fasterq-dump $r failed" >&2; exit 4; }
    rm -rf "$OUT/sra/$r" "$OUT/tmp"
  done \
  | awk -v pe="$pe" -v keep="$keep" 'BEGIN{srand(11)}
         NR%4==1{k=$1; sub(/\/[12]$/,"",k)} {rec[NR%4]=$0}
         NR%4==0{r=rec[1]"\n"rec[2]"\n"rec[3]"\n"rec[0]
                 if(!pe){if(rand()<keep){print r; n++}}
                 else if(k==pk){if(rand()<keep){print pr; print r; n++}; pk=""} else {pr=r; pk=k}}
         END{print n+0 > "'"$OUT"'/pairs.txt"}' \
  | bwa-mem2 mem $([ "$pe" = 1 ] && echo -p) -t "$NP" -R "@RG\tID:SRA-$bs\tSM:SRA-$bs\tPL:ILLUMINA\tLB:WGS" "$REF" - 2>"$OUT/bwa.log" \
  | samtools sort -n -m 400M -@ "$NP" -T "$OUT/tmp_sort" \
  | samtools fixmate -m -@ "$NP" - - \
  | samtools sort -m 400M -@ "$NP" -T "$OUT/tmp_sort2" \
  | samtools markdup -@ "$NP" --write-index - "$OUT/markdup.bam.tmp" \
  && [ "$(cat "$OUT/pairs.txt" 2>/dev/null || echo 0)" -gt 1000000 ] \
  && mv "$OUT/markdup.bam.tmp" "$OUT/markdup.bam" && mv "$OUT/markdup.bam.tmp.csi" "$OUT/markdup.bam.csi" \
  || { tail -n 3 "$OUT/prefetch.log" "$OUT/fqd.log" "$OUT/bwa.log" 2>/dev/null; fail "align (pairs=$(cat "$OUT/pairs.txt" 2>/dev/null))"; }
  rm -rf "$OUT/sra" "$OUT/tmp"
  echo "[$(date +%T)] pairs: $(cat "$OUT/pairs.txt")  BAM: $(du -sh "$OUT/markdup.bam" | cut -f1)"
fi
mdp=$(samtools idxstats "$OUT/markdup.bam" | awk -v rl="$(samtools view "$OUT/markdup.bam" chr1:10000000-10100000 | head -2000 | awk '{s+=length($10)} END{print (NR?s/NR:100)}')" \
      '$1 ~ /^chr[0-9]+$/ {r+=$3; L+=$2} END{printf "%.2f", r*rl/L}')
echo "[$(date +%T)] mapped depth ~${mdp}x"
contigs=$(awk '$1 ~ /^chr([0-9]+|X)$/ {print $1}' "$REF.fai")

if [ "$role" = panel ]; then
  # Panel v2 sites, alleles forced to the panel's, one chromosome per process.
  echo "$contigs" | xargs -P "$NP" -I{} sh -c \
    "[ -s $SITES/{}.tsv.gz ] || exit 0; bcftools mpileup -f $REF -r {} -T $SITES/{}.tsv.gz -q 30 -Q 20 -a AD,DP -Ou $OUT/markdup.bam 2>/dev/null \
     | bcftools call -m -C alleles -T $SITES/{}.tsv.gz -f GQ -Ou 2>/dev/null \
     | bcftools annotate -x INFO,FORMAT/PL,FORMAT/AD -Ob -o $OUT/p.{}.bcf" || fail "panel sites"
  bcftools concat -Ob -o "sra_panel/$bs.panel_sites.bcf.tmp" $(for c in $contigs; do [ -s "$OUT/p.$c.bcf" ] && echo "$OUT/p.$c.bcf"; done) 2>/dev/null \
    && bcftools index -f "sra_panel/$bs.panel_sites.bcf.tmp" \
    && mv "sra_panel/$bs.panel_sites.bcf.tmp" "sra_panel/$bs.panel_sites.bcf" && mv "sra_panel/$bs.panel_sites.bcf.tmp.csi" "sra_panel/$bs.panel_sites.bcf.csi" \
    || fail concat
  bcftools mpileup -f "$REF" -T "$ALLELES" -q 30 -Q 20 -a AD,DP -Ou "$OUT/markdup.bam" 2>/dev/null \
    | bcftools call -m -C alleles -T "$ALLELES" -f GQ -Ou 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%GQ\t%AD\t%DP]\n' | gzip > "sra_panel/$bs.panel.tsv.gz" || fail "breed panel"
  out="sites=$(bcftools index -n "sra_panel/$bs.panel_sites.bcf") bcf=$(du -h "sra_panel/$bs.panel_sites.bcf" | cut -f1)"
else
  # Truth: genome-wide gVCF at full depth (as sra-full-genotype.sh), plus a ~1x random subsample for imputation.
  echo "$contigs" | xargs -P "$NP" -I{} sh -c \
    "bcftools mpileup -f $REF -r {} -q 30 -Q 20 -a AD,DP -g 5 -Ou $OUT/markdup.bam 2>/dev/null \
     | bcftools call -m -g 5 -Oz -o $OUT/g.{}.vcf.gz 2>/dev/null" || fail gvcf
  bcftools concat -Oz -o "sra_full/$bs.g.vcf.gz.tmp" $(for c in $contigs; do echo "$OUT/g.$c.vcf.gz"; done) 2>/dev/null \
    && mv "sra_full/$bs.g.vcf.gz.tmp" "sra_full/$bs.g.vcf.gz" && tabix -f -p vcf "sra_full/$bs.g.vcf.gz" || fail concat
  L=scratch/panel_test/sra/$bs; mkdir -p "$L"
  frac=$(awk -v d="$mdp" 'BEGIN{x=1/d; if(x>=1) x=0.9999; printf "%.4f", x}'); frac=${frac#0}
  samtools view -@ "$NP" -b -s 7$frac -o "$L/lowpass.bam" "$OUT/markdup.bam" && samtools index -c "$L/lowpass.bam" || fail lowpass
  echo "$mdp" > "$L/full_depth.txt"
  out="gvcf=$(du -h "sra_full/$bs.g.vcf.gz" | cut -f1) lowpass=$(du -h "$L/lowpass.bam" | cut -f1)"
fi

SLICES=logs/candidate_slices.bed
if [ -s "$SLICES" ]; then
  samtools view -@ "$NP" -b -L "$SLICES" -o "sra_slices/$bs.bam.tmp" "$OUT/markdup.bam" \
    && samtools index "sra_slices/$bs.bam.tmp" "sra_slices/$bs.bam.bai" && mv "sra_slices/$bs.bam.tmp" "sra_slices/$bs.bam" || fail slice
fi
echo -e "$bs\t$group\t$breed\t$role\traw=${depth}x\tmapped=${mdp}x\t$out\tmin=$(( ($(date +%s)-t0)/60 ))" > "$DONE"
rm -rf "$OUT"
echo "SRAPANEL-DONE $(cat "$DONE")"
