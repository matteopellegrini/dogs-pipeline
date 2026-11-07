#!/bin/bash
# Build the canFam4 snpEff database (canFam4_NCBI_AR106) the pipeline annotates with. The ROS_Cfam_1.115 database
# snpEff ships is Ensembl ROS_Cfam_1.0, a different assembly from the canFam4 (UU_Cfam_GSD_1.0) positions the pipeline
# calls, so gene names and effects from it were computed against gene models at the wrong place.
# Source: NCBI RefSeq annotation release 106 of GCF_011100685.1 (= canFam4), chromosomes renamed RefSeq -> UCSC with
# the assembly report so they match canFam4.fa; gene records (transcript_id "") dropped, snpEff builds genes from
# transcripts. protein.fa / cds.fa are keyed by transcript ID (protein_id -> transcript_id from the GTF CDS lines) so
# snpEff's build checks compare its own translation of canFam4.fa against the NCBI proteins.
#   sbatch -A pi_matteop -p campus -c 2 --mem=24G -t 4:00:00 -J snpeffdb cluster/build_snpeff_canfam4.sh
# Output: envs/genomics/share/snpeff-*/data/canFam4_NCBI_AR106/snpEffectPredictor.bin (+ the genome entry in
# snpEff.config); build log snpeff_build/build.log.
#SBATCH -o logs/%x.%j.log
set -euo pipefail
D="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$D"
E="$D/envs/genomics"
JAVA="$(find "$E/lib/jvm/bin" -name java 2>/dev/null | head -1 || true)"; [[ -n "$JAVA" ]] || JAVA="$E/bin/java"
export PATH="$E/bin:$(dirname "$JAVA"):$PATH"
DB=canFam4_NCBI_AR106
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/100/685/GCF_011100685.1_UU_Cfam_GSD_1.0
P=GCF_011100685.1_UU_Cfam_GSD_1.0
B="$D/snpeff_build"; mkdir -p "$B"; cd "$B"
for f in genomic.gtf.gz assembly_report.txt protein.faa.gz; do [[ -s "${P}_$f" ]] || curl -sSf -O "$NCBI/${P}_$f"; done

SE="$(ls -d "$E"/share/snpeff-*/ | head -1)"
OUT="$SE/data/$DB"; mkdir -p "$OUT"
# RefSeq accession -> UCSC name (assembly report columns 7 and 10)
grep -v '^#' "${P}_assembly_report.txt" | awk -F'\t' '$7!="na"{sub(/\r$/,"",$10); print $7"\t"$10}' > refseq2ucsc.tsv
zcat "${P}_genomic.gtf.gz" | awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{m[$1]=$2; next}
  /^#/ || !($1 in m) || $9 ~ /transcript_id ""/ {next} {$1=m[$1]; print}' refseq2ucsc.tsv - > "$OUT/genes.gtf"
# protein.fa keyed by transcript ID
zcat "${P}_genomic.gtf.gz" | awk -F'\t' '$3=="CDS"' \
  | sed -n 's/.*transcript_id "\([^"]*\)".*protein_id "\([^"]*\)".*/\2\t\1/p' | sort -u > prot2tx.tsv
zcat "${P}_protein.faa.gz" | awk 'NR==FNR{t[$1]=t[$1] ? t[$1]" "$2 : $2; next}
  /^>/{id=substr($1,2); keep=(id in t); if(keep){n=split(t[id],a," "); hdr=a[1]} ; if(keep) print ">"hdr; next} keep' \
  prot2tx.tsv - > "$OUT/protein.fa"
ln -sf "$D/canFam4.fa" "$OUT/sequences.fa"
grep -q "^$DB.genome" "$SE/snpEff.config" \
  || printf '\n# canFam4 (UU_Cfam_GSD_1.0), NCBI RefSeq annotation release 106 - cluster/build_snpeff_canfam4.sh\n%s.genome : Dog_canFam4\n%s.reference : %s\n' \
       "$DB" "$DB" "$NCBI" >> "$SE/snpEff.config"
echo "genes.gtf: $(wc -l < "$OUT/genes.gtf") lines, $(cut -f1 "$OUT/genes.gtf" | sort -u | wc -l) contigs; protein.fa: $(grep -c '>' "$OUT/protein.fa") transcripts"
snpEff -Xmx20g build -gtf22 -v -noCheckCds "$DB" > build.log 2>&1
rm -f "$OUT/sequences.fa"   # the predictor stores the sequence it needs; canFam4.fa stays the one copy
grep -iE 'protein check|errors|ok:|error rate|warnings' build.log | tail -20
echo "SNPEFF-DB-DONE $DB"
