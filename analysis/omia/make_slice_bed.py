import gzip,re,json,csv,collections
span={}
for l in gzip.open('annot/canFam4.ncbiRefSeq.gtf.gz','rt'):
    f=l.split('\t')
    if f[2]!='transcript' or '_' in f[0]: continue
    g=re.search(r'gene_name "([^"]+)"',f[8]).group(1).upper()
    s,e=int(f[3]),int(f[4])
    if g in span and span[g][0]==f[0]: span[g]=(f[0],min(s,span[g][1]),max(e,span[g][2]))
    elif g not in span: span[g]=(f[0],s,e)
genes=set(); pts=[]
cat=json.load(open('/Users/matteopellegrini/Downloads/dogs/reference_json/omia_variants.json'))['variants']
for v in cat:
    genes.add((v.get('gene') or '').upper()); pts.append((v['chrom'],v['pos']))
for r in csv.DictReader(open('omia_dl/omia_dog_variants.tsv'),delimiter='\t'): genes.add(r['Gene'].upper())
for l in open('omia_dl/dog_genes.tsv'): genes.add(l.split('\t')[0].upper())
tr=json.load(open('/Users/matteopellegrini/Downloads/dogs/reference_json/trait_catalog.json'))['traits']
for t in tr:
    genes.add((t.get('gene') or '').upper())
    for L in t.get('loci',[]):
        if L.get('chrom') and L.get('pos'): pts.append((L['chrom'],int(L['pos'])))
for r in csv.DictReader(open('omia_missing_triage.tsv'),delimiter='\t'):
    if r['canfam4_pos']: pts.append((r['canfam4_chrom'],int(r['canfam4_pos'])))
alias={'PTPLA':'HACD1'}
iv=[]; found=0; missing=[]
for g in sorted(genes):
    if not g: continue
    s=span.get(g) or span.get(alias.get(g,''))
    if s and s[2]-s[1]<3_000_000: iv.append((s[0],max(0,s[1]-10000),s[2]+10000)); found+=1
    else: missing.append(g)
for c,p in pts: iv.append((c,max(0,p-50000),p+50000))
# merge
iv.sort(); m=[]
for c,s,e in iv:
    if m and m[-1][0]==c and s<=m[-1][2]: m[-1][2]=max(m[-1][2],e)
    else: m.append([c,s,e])
open('candidate_slices.bed','w').write(''.join(f'{c}\t{s}\t{e}\n' for c,s,e in m))
print('genes',len(genes),'found',found,'not in GTF',len(missing),missing[:25])
print('intervals',len(m),'total Mb',round(sum(e-s for c,s,e in m)/1e6,1))
