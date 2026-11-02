"""Curate the 65 catalogue entries left without alleles (one-off, 2026-11-01; chains/OMIA download lived in that session's scratchpad). Small variants get explicit left-normalised alleles;
SVs get symbolic records (SVTYPE, END, SVLEN) with breakpoints lifted from the source assembly."""
import csv,json,re,subprocess
from pyliftover import LiftOver
S='/private/tmp/claude-501/-Users-matteopellegrini-Downloads-dogs/0fb26461-2afd-4142-936c-a58295d35959/scratchpad'
FA='/Users/matteopellegrini/Downloads/dogs/canFam4.fa'
L={'cf3':LiftOver(S+'/chains/canFam3ToCanFam4.over.chain.gz'),'cf6':LiftOver(S+'/chains/canFam6ToCanFam4.over.chain.gz')}
def seq(c,s,e): return ''.join(subprocess.run(['samtools','faidx',FA,f'{c}:{s}-{e}'],capture_output=True,text=True).stdout.split('\n')[1:]).upper()
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def lift(a,c,p):
    if a=='cf4': return c,p,'+'
    r=L[a].convert_coordinate(c,p-1); assert r, (a,c,p); return r[0][0],r[0][1]+1,r[0][2]
V=json.load(open('/Users/matteopellegrini/Downloads/dogs/reference_json/omia_variants.json'))['variants']
SV=[ # sid, asm, chrom, start, end(or None for INS), svtype, svlen, source, note
('cat014','cf3','chr28',23380074,23483377,'DEL',None,'OMIA 1259',''),
('cat015','cf3','chr3',40812274,45768123,'INV',None,'OMIA 1086','~5 Mb inversion'),
('cat020','cf3','chr18',44791414,44890185,'DUP',None,'OMIA 1129','blue-eye 98.8 kb duplication'),
('cat025','cf3','chr5',32551064,None,'INS',227,'OMIA 734','SINE (LT796559.1:50-276) in exon 2'),
('cat026','cf3','chr26',8200944,None,'INS',205,'OMIA 1192','~205 bp insertion at c.2098-3'),
('cat027','cf3','chr7',41325010,41340731,'DEL',None,'OMIA 635','15.7 kb deletion replaced by AAGGCC'),
('cat032','cf4','chr29',33340028,33745116,'DEL',None,'catalogue','whole-gene CNGB3 deletion (~404 kb); OMIA 631 gives no coordinates'),
('cat035','cf3','chr20',40524267,40524380,'CPX',None,'OMIA 1276','7-bp del + 114-bp dup'),
('cat038','cf3','chrX',27631972,33069482,'INV',None,'OMIA 750','~5.4 Mb inversion'),
('cat042','cf3','chr4',80439639,86910352,'INV',None,'OMIA 752','6.47 Mb inversion, breakpoint in FAM134B intron 3'),
('cat044','cf4','chr18',48828546,48962003,'DUP',None,'catalogue','ridge 133 kb duplication (FGF3/FGF4/FGF19/ORAOV1)'),
('cat050','cf3','chr1',37564500,None,'INS',62,'OMIA 693',''),
('cat058','cf4','chr17',36516912,None,'INS',6401,'catalogue + OMIA 819','full-length LINE-1; OMIA gives no coordinates'),
('cat060','cf3','chr9',20407670,None,'INS',None,'OMIA 956','poly-A (40-70) + 11-bp duplication; length varies'),
('cat068','cf3','chr2',19371988,None,'INS',249,'OMIA 924','SINE 236 bp + 13-bp duplication'),
('cat069','cf4','chr8',14468376,14473136,'DEL',None,'OMIA 1613','catalogue position was 163 kb off'),
('cat077','cf3','chr27',20447905,None,'INS',233,'OMIA 700','SINE'),
('cat078','cf3','chr8',4236320,None,'INS',1980,'OMIA 706','LINE-1'),
('cat088','cf3','chrX',27622834,29823788,'INV',None,'OMIA 1234',''),
('cat089','cf3','chrX',27721607,None,'INS',4800,'OMIA 708','~4.8 kb insertion'),
('cat090','cf3','chrX',27851768,28247504,'DUP',None,'OMIA 1492',''),
('cat091','cf3','chr10',61822372,None,'INS',230,'OMIA 925','SINE'),
('cat093','cf4','chr12',33934716,None,'INS',3209,'catalogue + OMIA 853/855','FGF4 retrogene CFA12 (MF040221)'),
('cat094','cf3','chrX',109521130,None,'INS',1500,'OMIA 705','~1.5 kb insertion'),
('cat095','cf4','chrX',124076026,None,'INS',None,'catalogue','Havanese F8 insertion; length not in OMIA'),
('cat096','cf3','chr16',44477343,None,'INS',90,'OMIA 702',''),
('cat102','cf3','chr19',37903870,None,'INS',218,'OMIA 730','LINE fragment LN864704:123-340'),
('cat193','cf3','chr14',60628774,60758561,'DEL',None,'OMIA 642',''),
('cat195','cf3','chr4',73867275,73871357,'DEL',None,'OMIA 675',''),
('cat198','cf3','chr21',42583698,42587925,'DEL',None,'OMIA 638',''),
('cat227','cf3','chr37',30852988,30902901,'DEL',None,'OMIA 690',''),
('cat301','cf3','chr37',25698028,25705826,'DEL',None,'OMIA 632','CEA 7.8 kb deletion'),
('cat328','cf3','chr25',40417857,40417862,'CPX',None,'OMIA 1250','6 bp replaced by 35 bp of unknown sequence'),
('cat379','cf3','chrX',27615280,27982912,'DEL',None,'OMIA 1249',''),
('cat396','cf3','chr3',91523238,91523524,'DEL',None,'OMIA 1237',''),
('cat406','cf3','chr1',67734331,67736575,'DEL',None,'OMIA 1399',''),
('cat408','cf3','chr15',52905336,52927296,'DEL',None,'OMIA 943',''),
('cat418','cf3','chrX',120749179,120763577,'DEL',None,'OMIA 790',''),
('cat430','cf3','chr4',66946539,66963863,'DEL',None,'OMIA 1334','whole-gene SEPP1 deletion'),
('cat453','cf4','chr15',2716453,2717727,'DEL',None,'catalogue + OMIA 641','~1.27 kb deletion'),
('cat456','cf3','chr12',22603767,None,'INS',226,'OMIA 419','SINE'),
('cat469','cf4','chr18',20723884,None,'INS',5000,'catalogue + OMIA 694','FGF4 retrogene CFA18 (~5 kb)'),
('cat472','cf4','chr7',64457504,None,'INS',6500,'catalogue + OMIA 701','~6.5 kb satellite insertion, intron 35'),
]
# small variants: (sid, chrom, pos, ref, alt, source, note) with canFam4 coordinates already derived and checked
def ins_after(c,p,s): return (c,p,seq(c,p,p),seq(c,p,p)+s)
SMALL=[]
def add(sid,rec,src,note): SMALL.append((sid,)+rec+(src,note))
c,p,_=lift('cf3','chr16',6344748); add('cat030',ins_after(c,p,'T'),'OMIA 609','c.2647_2648insA (NM_001003124.2, codon 883)')
add('cat039',('chr18',52123541,'AG','A'),'OMIA 1481','c.3120delC (XM_038563927.1); OMIA g. is the VCF anchor')
add('cat041',('chr18',26126094,seq('chr18',26126094,26126106),seq('chr18',26126094,26126094)),'PMID 35385474','12-bp in-frame deletion removing VHFL (XM_038562964.1 c.1033_1044); catalogue position was 132 kb off')
c,p,s=lift('cf3','chr32',4528617); c2,p2,_=lift('cf3','chr32',4528632); a,b=sorted([p,p2])
add('cat046',(c,a-1,seq(c,a-1,b),seq(c,a-1,a-1)),'OMIA 952','lh3 c.556_571del (16 bp)')
c,p,_=lift('cf3','chr9',20134857); add('cat047',ins_after(c,p,'T'*60+seq(c,p-15,p)),'OMIA 1361','poly-T(60) + 16-bp duplication; poly-A length may vary')
add('cat071',('chr28',22371378,seq('chr28',22371378,22371400),seq('chr28',22371378,22371378)),'catalogue + OMIA 1323','22-bp deletion in exon 11 (checked inside exon 11)')
c,p,_=lift('cf3','chr15',18332036); add('cat073',ins_after(c,p,'A'*29+'GGAAGCAACAGGATG'),'OMIA 699','44-bp insertion c.142_143ins (codon 48 checked); poly-A tract length varies')
c,p,_=lift('cf3','chr11',19841331); add('cat070',ins_after(c,p,rc('TTTTTTTTTTTTTTTTTTTTTAGCCCTTGAAATTTT')),'OMIA 1138','36-bp insertion (gene on minus strand, sequence reverse-complemented)')
add('cat074',('chr1',48141426,'CACTT','C'),'Donner 2023 S1 + XM_038526433.1','c.182+1_182+4delGTAA (same as the T4 curated row)')
add('cat075',('chr5',43854679,seq('chr5',43854679,43854683),seq('chr5',43854679,43854679)),'PMID 32033218','c.1021_1024delTCAG')
add('cat164',('chr9',20763542,'GC','G'),'OMIA 1502','Dalmatian c.1107del; catalogue position was 540 kb off')
c,p,_=lift('cf3','chr5',82198114); add('cat173',ins_after(c,p,'G'),'OMIA 568','c.971_972insC (JHC)')
c,p,_=lift('cf3','chr4',73868000-725+0) if False else (None,None,None)
add('cat250',('chrX',124075494,seq('chrX',124075494,124075495),seq('chrX',124075494,124075494)),'OMIA 1655','Border Collie c.3206delA')
c,p,_=lift('cf3','chr24',23393552); add('cat361',(c,p,'C','T'),'OMIA 30','a allele R96C (NM_001007263.1 codon 96 checked); catalogue position was 14.5 kb off')
c,p,_=lift('cf3','chr4',52737279); add('cat417',(c,p-1,seq(c,p-1,p),seq(c,p-1,p-1)),'OMIA 563','c.744del (codon 248 checked)')
add('cat423',('chr17',19497253,seq('chr17',19497253,19497267),seq('chr17',19497253,19497253)),'catalogue + PMID 27160313','POMC 14-bp deletion (Raffan 2016); OMIA 809 g. spans 15 positions, catalogue 14-bp span kept')
c,p,_=lift('cf3','chr3',91747728); add('cat450',(c,p-1,seq(c,p-1,p+2),seq(c,p-1,p-1)),'OMIA 528','c.2407_2409del')
add('cat451',('chrX',33219537,seq('chrX',33219537,33219542),seq('chrX',33219537,33219537)),'catalogue + OMIA 480','XLPRA1 deletion of GAGAA (gene strand) checked')
add('cat452',('chrX',33219486,seq('chrX',33219486,33219488),seq('chrX',33219486,33219486)),'catalogue + OMIA 481','XLPRA2 deletion of GA (gene strand) checked')
add('cat466',('chr17',29430064,seq('chr17',29430064,29430072),'T'),'PMID 20596727','8-bp deletion CATACACA in the 3\' UTR')
add('cat471',('chr31',24668878,'ACTCCACAAA','GACT'),'PMID 28158183','c.383_392ACTCCACAAA>GACT (5\' non-coding exon)')
add('cat467',('chr14',20789292,seq('chr14',20789292,20789308),seq('chr14',20789292,20789292)),'catalogue (hgvs: intron 10 donor 16-bp deletion)','removes the intron-10 GT donor (exon 10 ends at 20789309, minus strand)')
OPEN=[]
_x=[('cat467','PDK4 16-bp deletion: catalogue position is at the 3\' end of the gene, not the intron-1 donor the paper describes (Meurs 2012, paywalled); needs the deleted sequence')]
w=csv.writer(open(S+'/t4/curated65.tsv','w'),delimiter='\t')
w.writerow(['site_id','gene','chrom','pos','ref','alt','kind','svtype','end','svlen','source','note','cat_pos','delta'])
for sid,ch,pos,ref,alt,src,note in SMALL:
    v=V[int(sid[3:])]; w.writerow([sid,v['gene'],ch,pos,ref,alt,'small','','','',src,note,v['pos'],pos-int(v['pos'])])
for sid,a,c,s,e,t,ln,src,note in SV:
    v=V[int(sid[3:])]; c4,s4,st=lift(a,c,s)
    if e: c4b,e4,_=lift(a,c,e); s4,e4=sorted([s4,e4])
    else: e4=s4
    ln=ln if ln else (e4-s4+1 if t in('DEL','DUP','INV','CPX') else '')
    w.writerow([sid,v['gene'],c4,s4,seq(c4,s4,s4),f'<{t}>','sv',t,e4,ln,src,note,v['pos'],s4-int(v['pos'])])
for sid,note in OPEN:
    v=V[int(sid[3:])]; w.writerow([sid,v['gene'],v['chrom'],v['pos'],'.','.','open','','','','',note,v['pos'],0])
print(len(SMALL),'small',len(SV),'sv',len(OPEN),'open')
