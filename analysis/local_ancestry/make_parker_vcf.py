# Parker array genotypes (PLINK bed, 1,356 dogs, 143,933 sites) -> unphased
# VCF at the 131k Parker panel sites, with alleles oriented to the Dog10K
# VCF REF/ALT so Beagle/FLARE see one allele system. Sites whose alleles
# do not match Dog10K (after allowing a swap) are dropped and counted.
import gzip, sys, time, os, numpy as np
while not os.path.exists("make_dog10k_vcf.log") or "rows:" not in open("make_dog10k_vcf.log").read():
    time.sleep(30)
prefix="../COSMO/analysis/cosmo_parker_full"
fam=[l.split() for l in open(prefix+".fam")]; ids=[f[1] for f in fam]; fids=[f[0] for f in fam]
bim=[l.split() for l in open(prefix+".bim")]
n_ind=len(fam); n_snp=len(bim); bps=(n_ind+3)//4
data=np.fromfile(prefix+".bed",dtype=np.uint8,offset=3).reshape(n_snp,bps)
codes=np.empty((n_snp,bps*4),dtype=np.uint8)
for k in range(4): codes[:,k::4]=(data>>(2*k))&3
# PLINK bed 2-bit: 0=hom A1, 1=missing, 2=het, 3=hom A2  -> count of A1 alleles
a1count=np.array([2,-1,1,0],dtype=np.int8)[codes[:,:n_ind]]
bykey={}
for i,b in enumerate(bim): bykey[(("chr"+b[0]) if not b[0].startswith("chr") else b[0], b[3])]=i
out=gzip.open("parker_131k.vcf.gz","wt")
out.write("##fileformat=VCFv4.2\n##source=COSMO/analysis/cosmo_parker_full (Parker array) at breed_panel/sites.tsv, alleles oriented to Dog10K\n")
out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
contigs=[]; rows=[]
for l in gzip.open("dog10k_131k.vcf.gz","rt"):
    if l.startswith("##contig"): contigs.append(l); continue
    if l.startswith("#"): continue
    p=l.split("\t",5); rows.append((p[0],p[1],p[3],p[4]))
for c in contigs: out.write(c)
out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(ids)+"\n")
kept=drop=miss=0
gtmap={2:"0/0",1:"0/1",0:"1/1",-1:"./."}   # keyed by ALT-allele count -> written below
for chrom,pos,ref,alt in rows:
    i=bykey.get((chrom,pos))
    if i is None: miss+=1; continue
    A1,A2=bim[i][4],bim[i][5]
    if A1==alt and A2==ref:   altcount=a1count[i]
    elif A1==ref and A2==alt: altcount=np.where(a1count[i]<0,-1,2-a1count[i])
    else: drop+=1; continue
    out.write(f"{chrom}\t{pos}\t{chrom}_{pos}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t"+"\t".join({2:"1/1",1:"0/1",0:"0/0",-1:"./."}[int(x)] for x in altcount)+"\n")
    kept+=1
out.close()
open("parker_samples.tsv","w").write("".join(f"{s}\t{f}\n" for s,f in zip(ids,fids)))
print(f"kept {kept} sites, dropped(allele mismatch) {drop}, not in array {miss}, samples {n_ind}")
