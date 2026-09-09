#$ -cwd
#$ -j y
#$ -o logs/pathogen_panel.$TASK_ID.log
#$ -l h_data=4G,h_rt=1:00:00
cd $SGE_O_WORKDIR
export LD_LIBRARY_PATH=$PWD/envs/genomics/lib:${LD_LIBRARY_PATH:-}
S=$PWD/envs/genomics/bin
dog=$(sed -n "${SGE_TASK_ID}p" pathogen_screen/bam_list.txt | cut -f1)
B=$(sed -n "${SGE_TASK_ID}p" pathogen_screen/bam_list.txt | cut -f2)
[ -n "$dog" ] && [ -f "$B" ] || exit 0
mkdir -p pathogen_screen/panel_out
# Competitive mapping vs the 7-genome tick/zoonotic panel; per-hit species
# + 10kb bin at >=97% identity. Verdicts need genome BREADTH (aggregator
# masks bins present across the cohort = conserved cross-map loci).
$S/samtools fastq -f 4 "$B" 2>/dev/null \
 | $S/minimap2 -ax sr pathogen_screen/panel.fna - 2>/dev/null \
 | $S/samtools view -q 30 -F 0x904 - \
 | awk '{nm=0; for(i=12;i<=NF;i++) if($i ~ /^NM:i:/){split($i,a,":"); nm=a[3]}
        len=length($10); if(len>0 && (len-nm)/len>=0.97){split($3,r,"|"); print r[1]"\t"r[2]"\t"int($4/10000)}}' \
 | sort | uniq -c | awk -v d="$dog" '{print d"\t"$2"\t"$3"\t"$4"\t"$1}' > pathogen_screen/panel_out/$dog.tsv
echo PANEL-DONE $dog rows=$(wc -l < pathogen_screen/panel_out/$dog.tsv)
