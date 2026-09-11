# Phased Dog10K haplotypes at the 131k Parker sites -> VCF (from the
# relatives extraction, which kept the 0|1 phase from the GLIMPSE2 panel).
import gzip, sys
samples=[l.strip() for l in open("../relatives/panel_samples.txt") if l.strip()]
out=gzip.open("dog10k_131k.vcf.gz","wt")
out.write("##fileformat=VCFv4.2\n##source=dog10k_panel/AutoAndXPAR.Dog10K.phased.bcf at breed_panel/sites.tsv\n")
out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype\">\n")
chroms=[]
n=0
for l in gzip.open("../relatives/panel_gt.tsv.gz","rt"):
    p=l.rstrip("\n").split("\t")
    c=p[0]
    if not chroms or chroms[-1]!=c: chroms.append(c)
    n+=1
for c in chroms: out.write(f"##contig=<ID={c}>\n")
out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(samples)+"\n")
for l in gzip.open("../relatives/panel_gt.tsv.gz","rt"):
    p=l.rstrip("\n").split("\t")
    gts=[g if g else "./." for g in p[4:4+len(samples)]]
    out.write(f"{p[0]}\t{p[1]}\t{p[0]}_{p[1]}\t{p[2]}\t{p[3]}\t.\tPASS\t.\tGT\t"+"\t".join(gts)+"\n")
out.close()
print("rows:", n, "samples:", len(samples))
