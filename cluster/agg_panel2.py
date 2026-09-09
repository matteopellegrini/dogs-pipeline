import glob
from collections import defaultdict, Counter
# Per-contig bins: key on (species, contig, bin) so multi-plasmid Borrelia
# breadth is counted genome-wide, and conserved cross-map loci are masked
# per contig.
sb_dogs=defaultdict(set)
dogdata=defaultdict(lambda: defaultdict(dict))  # dog -> sp -> (contig,bin) -> reads
dogs=set()
for f in glob.glob("/u/project/pellegrini/gkislik/dogs/pathogen_screen/panel_out/*.tsv"):
    dog=f.split("/")[-1][:-4]; dogs.add(dog)
    for l in open(f):
        p=l.rstrip("\n").split("\t")
        if len(p)!=5: continue
        d,sp,ctg,b,n=p[0],p[1],p[2],int(p[3]),int(p[4])
        dogdata[dog][sp][(ctg,b)]=dogdata[dog][sp].get((ctg,b),0)+n
        sb_dogs[(sp,ctg,b)].add(dog)
N=len(dogs)
conserved={k for k,v in sb_dogs.items() if len(v) > 0.05*N}
print("dogs:", N, "conserved(masked) loci:", len(conserved))
for sp,c in sorted(Counter(sp for sp,ctg,b in conserved).items()): print("  masked", sp, c, "bins")
tops=defaultdict(list); pos=[]
for dog in dogs:
    for sp,bins in dogdata[dog].items():
        info={k:v for k,v in bins.items() if (sp,k[0],k[1]) not in conserved}
        nb=len(info); reads=sum(info.values())
        tops[sp].append((nb,reads,dog))
        if nb>=15 and reads>=40: pos.append((dog,sp,nb,reads))
print("=== POSITIVES (>=15 informative bins & >=40 reads):", len(pos), "===")
for dog,sp,nb,reads in sorted(pos,key=lambda x:-x[2]): print(f"  {dog} {sp} bins={nb} reads={reads}")
print("=== max informative breadth per species ===")
for sp in sorted(tops):
    nb,reads,dog=max(tops[sp]); print(f"  {sp:26s} max_bins={nb:3d} reads={reads:4d} ({dog})")
