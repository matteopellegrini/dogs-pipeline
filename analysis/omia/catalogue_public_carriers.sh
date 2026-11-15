#!/bin/bash
# How often is each Panel v3 catalogue variant (reference_panel/panel_v3/bp_sites.tsv, small variants) seen in the
# non-customer dogs we hold? Three sources, each written as long tables to work_panel_v3/public/:
#   panel.tsv.gz  imputation reference panel v2 (Dog10K + US cohort), only the catalogue sites it contains: key sample GT
#   slices.tsv.gz public genomes with a read slice (sra_slices/*.bam, ~10x), forced-allele calls: key dog GT DP
#   gvcf.tsv.gz   full-depth public gVCFs (sra_full/*.g.vcf.gz) whose dog has no slice: key dog GT DP (ALT "." = hom-ref)
# key = CHROM:POS:REF:ALT after bcftools norm; sites.tsv maps each key to its catalogue site_id.
# Summarised locally by analysis/omia/catalogue_public_carriers.py. Rerun as the streamed dogs add slices.
#   sbatch -A pi_matteop -p campus -c 8 --mem=16G -t 4:00:00 -J catpub analysis/omia/catalogue_public_carriers.sh
#SBATCH -o logs/%x.%j.log
set -uo pipefail
D="${SLURM_SUBMIT_DIR:?run via sbatch from the repo root}"; cd "$D"
export OPENBLAS_NUM_THREADS=1 PATH=$D/envs/genomics/bin:$PATH
P=reference_panel/panel_v3; O=work_panel_v3/public; FASTA=canFam4_plusY.fa; NP=${SLURM_CPUS_PER_TASK:-8}
mkdir -p $O

# Catalogue small variants -> normalised keys (bp_sites keeps original coordinates; calls are left-normalised)
{ printf '##fileformat=VCFv4.2\n'; awk '{printf "##contig=<ID=%s,length=%s>\n",$1,$2}' $FASTA.fai
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  awk -F'\t' -v OFS='\t' 'NR>1 && $4~/^[ACGT]+$/ && $5~/^[ACGT]+$/{print $2,$3,$1,$4,$5,".",".","."}' $P/bp_sites.tsv | sort -k1,1V -k2,2n
} | bcftools norm -f $FASTA -c w 2>$O/norm.log | bcftools query -f '%ID\t%CHROM:%POS:%REF:%ALT\t%CHROM\t%POS\n' > $O/sites.tsv
cut -f3,4 $O/sites.tsv | sort -u -k1,1V -k2,2n > $O/regions.tsv
echo "sites: $(wc -l < $O/sites.tsv)"

# 1. Reference panel: only non-ref genotypes are kept, plus the number of called dogs per site
bcftools view -R $O/regions.tsv ref_panel_v2/AutoAndXPAR.Dog10K_DAP.bcf -Ou | bcftools norm -f $FASTA -m -any -Ou \
  | bcftools query -f '%CHROM:%POS:%REF:%ALT[\t%SAMPLE=%GT]\n' \
  | awk -F'\t' -v OFS='\t' '{n=0; for(i=2;i<=NF;i++){split($i,a,"="); g=a[2]; if(g ~ /\./) continue; n++; if(g ~ /1/) print $1,a[1],g}
                             print $1,"__called__",n}' | gzip > $O/panel.tsv.gz
echo "panel done"

# 2. Read slices: forced-allele calls at every catalogue site (as cluster/panel_v3_genotype.sh)
zcat $P/v3_small.alleles.tsv.gz > $O/alleles.tsv; bgzip -f $O/alleles.tsv; tabix -f -s1 -b2 -e2 $O/alleles.tsv.gz
ls sra_slices/*.bam > $O/slices.list
bcftools mpileup --threads $NP -f $FASTA -R $P/v3_small.regions.tsv -q 20 -Q 20 -B -m 1 -d 1000 -a FORMAT/AD,FORMAT/DP \
    -b $O/slices.list -Ou 2>$O/mpileup.log \
  | bcftools call -m -A -C alleles -T $O/alleles.tsv.gz -Ou 2>/dev/null | bcftools norm -f $FASTA -m -any -Ou 2>/dev/null \
  | bcftools query -f '%CHROM:%POS:%REF:%ALT[\t%SAMPLE=%GT=%DP]\n' \
  | awk -F'\t' -v OFS='\t' '{for(i=2;i<=NF;i++){split($i,a,"="); print $1,a[1],a[2],a[3]}}' | gzip > $O/slices.tsv.gz
echo "slices done: $(wc -l < $O/slices.list) BAMs"

# 3. Full-depth gVCFs for dogs without a slice (bcftools call -g gVCFs: hom-ref sites have ALT ".")
: > $O/gvcf.tsv
for g in sra_full/*.g.vcf.gz; do
  id=$(basename $g .g.vcf.gz); [ -s sra_slices/$id.bam ] && continue
  bcftools view -R $O/regions.tsv $g 2>/dev/null | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t$id[\t%GT\t%DP]\n" >> $O/gvcf.tsv
done
gzip -f $O/gvcf.tsv
echo "CATPUB-DONE"
