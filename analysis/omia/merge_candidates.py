import csv, re, collections, json
rd=lambda f: list(csv.DictReader(open(f),delimiter='\t'))
G=lambda s: re.split(r'[ (/,;]',(s or '').strip().upper())[0]
rows=[]
def add(**k):
    base=dict(gene='',condition='',breeds='',variant='',canfam4='',vclass='',evidence='causal',sources=set(),freq='',omia_id='',pmid='',year='',notes='')
    base.update(k); base['sources']=set(base['sources']); rows.append(base); return base
# A. OMIA gaps
for r in rd('omia_missing_triage.tsv'):
    pos=f"{r['canfam4_chrom']}:{r['canfam4_pos']}" if r['canfam4_pos'] else ''
    ev='trait' if r['deleterious']=='no' else 'causal'
    add(gene=r['gene'],condition=r['phenotype'],breeds=r['breeds'],variant=' '.join(x for x in [r['c'],r['p']] if x)[:120] or r['g'][:80],
        canfam4=pos,vclass=r['class'],evidence=ev,sources={'OMIA'},omia_id=r['omia_variant_id'],pmid=r['pmid'],year=r['year'],
        notes=r['coord_status'] if r['coord_status'] not in ('native','lifted') else '')
byg=collections.defaultdict(list)
for x in rows: byg[G(x['gene'])].append(x)
def attach_or_add(gene,src,**k):
    hits=byg.get(G(gene),[])
    if hits:
        for h in hits: h['sources'].add(src)
        return hits[0]
    x=add(gene=gene,sources={src},**k); byg[G(gene)].append(x); return x
# B. commercial
for r in rd('commercial_missing.tsv'):
    srcs=[s.strip() for s in re.split(r'[;,]',r['sources']) if s.strip()]
    g=r['gene'] if r['gene'].lower()!='not stated' else ''
    if g:
        x=attach_or_add(g,'commercial',condition=r['condition'],breeds=r['breed(s)'],variant=r['variant_detail_if_any'][:120])
    else:
        x=add(gene='(not stated)',condition=r['condition'],breeds=r['breed(s)'],sources={'commercial'},notes='no gene/variant published')
    x['sources'].update(srcs)
    if not x['omia_id'] and ('DLA' in g.upper() or 'marker' in r['condition'].lower() or 'risk' in r['condition'].lower() or 'suscept' in r['condition'].lower()): x['evidence']='risk'
# C. screening not in catalogue
for r in rd('screening_frequencies.tsv'):
    if r['in_our_catalogue']!='no': continue
    x=attach_or_add(r['gene'] or '(not stated)','Donner screening',condition=r['disease'],variant=r['variant'][:120],canfam4=r['canfam4'])
    x['freq']=r['freq_overall']; x['sources'].add('Donner screening')
# D. literature
for r in rd('lit_2025_2026.tsv'):
    if r['status'] not in ('NEW','known_variant_not_in_OMIA_or_catalogue'): continue
    ev='risk' if re.search(r'risk|assoc|GWAS|severity|suscept|metaboli',r['evidence']+r['disease'],re.I) else 'causal'
    if re.search(r'weak|unresolved',r['evidence'],re.I): ev='weak'
    x=add(gene=r['gene'],condition=r['disease'],breeds=r['breed'],variant=r['variant'][:120],sources={'PubMed 2025-26' if r['status']=='NEW' else 'PubMed (pre-2025)'},
          pmid=r['pmid'],year=r['year'],evidence=ev,notes=('assembly '+r['assembly']) if r['assembly'] else 'coordinates not given')
# E. HGTD
for r in rd('hgtd_missing.tsv'):
    if r['status']=='known_gap':
        for h in byg.get(G(r['gene']),[]): h['sources'].add('HGTD')
        continue
    ev='risk' if 'risk' in r['condition'].lower() else ('trait' if 'tan' in r['condition'].lower() else 'causal')
    if r['status']=='possible_additional_allele':
        attach_or_add(r['gene'],'HGTD',condition=r['condition'],breeds=r['breeds'])['notes']+=' HGTD: '+r['breeds'][:40]
    else:
        add(gene=r['gene'] or '(not stated)',condition=r['condition'],breeds=r['breeds'],variant=r['variant_detail'][:120],sources={'HGTD'},evidence=ev,
            notes='' if r['gene'] else 'no gene/variant published')
# F. major coat/trait loci absent from the catalogue (VGL/Orivet trait tests)
for g,c in [('PMEL','Merle (SINE insertion, allele lengths)'),('MITF','Piebald / white spotting (SINE + length polymorphism)'),('RSPO2','Furnishings (3\' UTR insertion)')]:
    attach_or_add(g,'UCDavisVGL',condition=c,evidence='trait',vclass='structural')['evidence']='trait'
# hand corrections after review (2026-11-01)
for x in rows:
    g=x['gene'].upper()
    if g=='ABCB1': x['evidence']='causal'                      # MDR1 alleles: OMIA marks 'not deleterious', but drug sensitivity is reported
    if g=='RETN': x['evidence']='risk'                         # copper toxicosis modifier
    if g.startswith('HS RISK') or 'histiocytic' in x['condition'].lower(): x['evidence']='risk'
    if g=='SIRT6': x['evidence']='causal'                      # dominant hyperparathyroidism, Keeshond
# tiers
def tier(x):
    if x['evidence']=='trait': return 'T6 trait/coat'
    if x['gene'] in ('(not stated)','') or 'no gene' in x['notes']: return 'T7 no published variant'
    if x['evidence'] in ('risk','weak'): return 'T5 risk/association (panel only)'
    if not x['canfam4'] and x['sources']<= {'OMIA'} : return 'T4 needs coordinates'
    pri = bool(x['freq']) or len(x['sources']-{'OMIA'})>=1
    if x['vclass']=='structural': return 'T3 structural' + (' (priority)' if pri else '')
    if not x['canfam4']: return 'T4 needs coordinates' + (' (priority)' if pri else '')
    return 'T1 add now (priority)' if pri else 'T2 add now'
for x in rows: x['tier']=tier(x); x['n_sources']=len(x['sources']); x['sources']='; '.join(sorted(x['sources']))
rows.sort(key=lambda x:(x['tier'],-x['n_sources'],x['gene']))
cols=['tier','gene','condition','breeds','variant','canfam4','vclass','evidence','sources','n_sources','freq','omia_id','pmid','year','notes']
with open('catalogue_expansion_candidates.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
c=collections.Counter(x['tier'] for x in rows)
for k in sorted(c): print(f'{c[k]:4d}  {k}')
print('total',len(rows))
