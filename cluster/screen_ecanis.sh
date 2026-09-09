#$ -cwd
#$ -j y
#$ -o logs/ecanis_screen.$TASK_ID.log
#$ -l h_data=4G,h_rt=1:00:00
cd $SGE_O_WORKDIR
export LD_LIBRARY_PATH=$PWD/envs/genomics/lib:${LD_LIBRARY_PATH:-}
S=$PWD/envs/genomics/bin
dog=$(sed -n "${SGE_TASK_ID}p" pathogen_screen/bam_list.txt | cut -f1)
B=$(sed -n "${SGE_TASK_ID}p" pathogen_screen/bam_list.txt | cut -f2)
[ -n "$dog" ] && [ -f "$B" ] || exit 0
mkdir -p pathogen_screen/out
# unmapped reads -> E. canis; per-hit identity and 10kb-bin position.
# Verdict logic lives in the aggregator: true positives need genome BREADTH
# (many distinct bins at high identity), not just read count — conserved
# rRNA loci (bins 28/109) cross-map in every dog.
$S/samtools fastq -f 4 "$B" 2>/dev/null \
 | $S/minimap2 -ax sr pathogen_screen/ecanis.fna - 2>/dev/null \
 | $S/samtools view -q 30 -F 0x904 - \
 | awk '{nm=0; for(i=12;i<=NF;i++) if($i ~ /^NM:i:/){split($i,a,":"); nm=a[3]}
        len=length($10); if(len>0 && (len-nm)/len>=0.97) print int($4/10000)}' \
 | sort -n | uniq -c | awk -v d="$dog" '{print d"\t"$2"\t"$1}' > pathogen_screen/out/$dog.tsv
echo ECANIS-DONE $dog bins=$(wc -l < pathogen_screen/out/$dog.tsv)
