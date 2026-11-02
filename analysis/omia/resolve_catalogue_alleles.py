"""One-off (2026-11-01): fill alleles for catalogue entries that carry ref/alt None, from OMIA g. notation
(gene match, lifted to canFam4, reference base checked, within 15 bp and same variant type). Inputs were the OMIA
download and chains in that session's scratchpad; the committed output is known_variants/panel_v3/catalogue_alleles_resolved.tsv."""
import csv,json,re,subprocess,collections
from pyliftover import LiftOver
S='/private/tmp/claude-501/-Users-matteopellegrini-Downloads-dogs/0fb26461-2afd-4142-936c-a58295d35959/scratchpad'
FA='/Users/matteopellegrini/Downloads/dogs/canFam4.fa'
LO={'CanFam3.1':LiftOver(S+'/chains/canFam3ToCanFam4.over.chain.gz'),'Dog10K_Boxer_Tasha':LiftOver(S+'/chains/canFam6ToCanFam4.over.chain.gz')}
def seq(c,s,e): return ''.join(subprocess.run(['samtools','faidx',FA,f'{c}:{s}-{e}'],capture_output=True,text=True).stdout.split('\n')[1:]).upper()
comp=str.maketrans('ACGT','TGCA'); rc=lambda s:s.translate(comp)[::-1]
def nasm(a):
    return 'UU' if 'Cfam_GSD' in a else 'CanFam3.1' if 'CanFam3' in a else 'Dog10K_Boxer_Tasha' if ('Tasha' in a or 'Dog10K' in a) else None
def conv(asm,chrom,g):
    """simple single-event g. -> (chrom,pos,ref,alt) on canFam4, or None"""
    gg=g.split(':')[-1].replace(',','').replace(' ','')
    m=re.fullmatch(r'[gG]\.(\d+)(?:_(\d+))?(del[ACGT]*|dup[ACGT]*|ins[ACGT]+|delins[ACGT]+|[ACGT]>[ACGT])',gg)
    if not m: return None
    s=int(m.group(1)); e=int(m.group(2) or s); op=m.group(3)
    if asm=='UU': S_,E_,st=s,e,'+'
    else:
        a=LO[asm].convert_coordinate(chrom,s-1); b=LO[asm].convert_coordinate(chrom,e-1)
        if not a or not b or a[0][0]!=chrom: return None
        S_,E_,st=a[0][1]+1,b[0][1]+1,a[0][2]
        if st=='-': S_,E_=E_,S_
        if abs(E_-S_)!=e-s: return None
    if '>' in op:
        r,x=op[0],op[2]
        if st=='-': r,x=rc(r),rc(x)
        return (chrom,S_,r,x) if seq(chrom,S_,S_)==r else ('MISMATCH',)
    if op.startswith('delins'):
        ins=op[6:]; ins=ins if st=='+' else rc(ins)
        return (chrom,S_-1,seq(chrom,S_-1,E_),seq(chrom,S_-1,S_-1)+ins)
    if op.startswith('del'):
        d=op[3:]; ref=seq(chrom,S_-1,E_)
        if d and (d if st=='+' else rc(d))!=ref[1:]: return ('MISMATCH',)
        return (chrom,S_-1,ref,ref[0])
    if op.startswith('dup'): return (chrom,E_,seq(chrom,E_,E_),seq(chrom,E_,E_)+seq(chrom,S_,E_))
    if op.startswith('ins'):
        ins=op[3:]; ins=ins if st=='+' else rc(ins)
        if st=='-': S_=E_
        return (chrom,S_,seq(chrom,S_,S_),seq(chrom,S_,S_)+ins)
om=collections.defaultdict(list)
for r in csv.DictReader(open(S+'/omia_dl/omia_dog_variants.tsv'),delimiter='\t'):
    for g in re.split(r'[ ,;/]+',r['Gene']): om[g.upper()].append(r)
V=json.load(open('/Users/matteopellegrini/Downloads/dogs/reference_json/omia_variants.json'))['variants']
res=[]; stat=collections.Counter()
for i,v in enumerate(V):
    if str(v['ref'])!='None': continue
    pos=int(v['pos']); pend=int(v.get('pos_end') or pos); best=None; why='no_omia_gene'
    for r in om.get(v['gene'].upper(),[]):
        asm=nasm(r['Reference Sequence'])
        if not asm or not r['g. or m.']: why='omia_no_g' if why=='no_omia_gene' else why; continue
        out=conv(asm,v['chrom'],r['g. or m.'])
        if not out: why='unparsed' if why in('no_omia_gene','omia_no_g') else why; continue
        if out[0]=='MISMATCH': why='ref_mismatch'; continue
        if abs(out[1]-pos)<=15 or (pos-15<=out[1]<=pend+15):
            best=(out,r['OMIA Variant ID'],r['g. or m.']); break
        else: why='pos_far' if why!='ref_mismatch' else why
    stat['resolved' if best else why]+=1
    res.append([f'cat{i:03d}',v['gene'],v['chrom'],pos,v.get('variant_type',''),*(list(best[0])+[best[1],best[2]] if best else ['','','','','',why])])
csv.writer(open(S+'/t4/none_resolved.tsv','w'),delimiter='\t').writerows(res)
print(stat)
