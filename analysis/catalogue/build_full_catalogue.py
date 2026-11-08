#!/usr/bin/env python3
"""Full catalogue of what a report screens: every condition/trait variant, its variant class (SNV, indel,
structural/CNV, mitochondrial) and how it is called, plus the 25 single-gene traits, the coat-colour loci and the
genome-wide scans. Joins reference_json/omia_variants_v3.json (the 601 report entries), the panel v3 site tables
(engine + cohort call counts), dogs-app/lib/conditionCatalog.json (condition grouping, as the app shows it) and
reference_json/trait_catalog.json.
  python3 analysis/catalogue/build_full_catalogue.py OUT.xlsx
"""
import csv, json, sys, collections, pathlib
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

R = pathlib.Path(__file__).resolve().parents[2]
OUT = sys.argv[1]
V3 = json.load(open(R / 'reference_json/omia_variants_v3.json'))
CC = json.load(open(R / 'dogs-app/lib/conditionCatalog.json'))
BP = {r['site_id']: r for r in csv.DictReader(open(R / 'reference_panel/panel_v3/bp_sites.tsv'), delimiter='\t')}
BS = {r['site_id']: r for r in csv.DictReader(open(R / 'reference_panel/panel_v3/bp_sv_sites.tsv'), delimiter='\t')}
SV = {r['site_id']: r for r in csv.DictReader(open(R / 'reference_panel/panel_v3/sv_sites.tsv'), delimiter='\t')}
TR = json.load(open(R / 'reference_json/trait_catalog.json'))
HELD, DROPPED = set(V3['held']), set(CC.get('dropped') or [])
RETRO = {'cat093', 'cat469'}
MODE = {'AR': 'Recessive', 'AD': 'Dominant', 'SD': 'One copy raises risk', 'XR': 'X-linked recessive',
        'XD': 'X-linked dominant', 'ML': 'Risk marker in males', 'complex': 'Risk variant', 'mixed': 'Varies by variant'}
CATN = {'disease': 'Health', 'coat': 'Coat colour', 'trait': 'Trait'}

M_IMP = 'Imputed (GLIMPSE2 on the Dog10K panel) x reads'
M_RDS = 'Direct reads x breed-ancestry prior'
M_HAP = 'Haplotype imputation (Dog10K indel haplotypes, Beagle-phased) x reads'
M_MT = 'Mitochondrial allele fraction (reads, MAPQ 0 allowed)'
M_JD = 'Junction reads (split / clipped / discordant pairs) + read depth'
M_J = 'Junction reads (split / clipped / discordant pairs)'
M_JU = 'Junction reads, copy not confirmable (reports clear or uncertain only)'
M_D = 'Read depth only (vs distant reference windows; no published breakpoints)'
M_RT = 'Retrogene insertion: mate pairs landing in the parent gene'
M_NONE = 'Not genotyped (no reliable call at this site)'

def vkey(e):
    return f"{e.get('chrom')}:{e.get('pos')}:{str(e.get('gene') or '').strip().upper()}"

def cond_of(e):
    hit = CC['variants'].get(f"{vkey(e)}:{e.get('hgvs_c')}") or CC['variants'].get(vkey(e))
    if not hit:
        hit = next((v for k, v in CC['variants'].items() if k.startswith(f"{e.get('chrom')}:{e.get('pos')}:")), None)
    if not hit:
        return None, None, None
    c = CC['conditions'].get(hit['c'], {})
    return c.get('name'), hit['cat'], hit['m']

def sv_confirmable(s):
    n = int(s.get('svlen') or 0)
    if s['site_id'] in RETRO: return True
    if s['type'] == 'INS': return 2 <= n <= 60
    return int(s['end']) > int(s['start']) or 2 <= n <= 60

rows = []
for e in V3['variants']:
    sid = e['site_id']
    name, cat, mode = cond_of(e)
    name = name or (e.get('phene_name') or e.get('trait') or f"Inherited variant in {e.get('gene')}").strip()
    cat = cat or e.get('category') or 'disease'
    ref, alt = e.get('ref'), e.get('alt')
    allele = f"{ref}>{alt}" if ref and alt and ref not in ('.',) else (e.get('hgvs_c') or e.get('variant_type') or '')
    cohort = {}
    if sid in BP:
        b = BP[sid]
        size = abs(len(b['ref']) - len(b['alt']))
        if b['chrom'] == 'chrM':
            vclass = 'Mitochondrial SNV' if size == 0 else 'Mitochondrial indel'
        elif b['kind'] == 'snv':
            vclass = 'SNV'
        else:
            vclass = f"Indel ({'+' if len(b['alt']) > len(b['ref']) else '-'}{size} bp)" + (' in homopolymer' if b['kind'] == 'homopolymer' else '')
        src = b['source']
        method = {'GLIMPSE (Dog10K)': M_IMP, 'reads + breed prior': M_RDS, 'chrM allele fraction': M_MT}.get(src, M_HAP if 'Haplotype' in src else src)
        allele = f"{b['ref']}>{b['alt']}" if len(b['ref']) + len(b['alt']) < 30 else f"{len(b['ref'])}bp>{len(b['alt'])}bp"
        cohort = b; rclass = b['report_class']; extra = b['documented_breeds']
    elif sid in BS:
        b, s = BS[sid], SV[sid]
        n = int(s['svlen'] or 0) or (int(s['end']) - int(s['start']))
        t = s['type']
        tname = {'DEL': 'deletion', 'DUP': 'duplication', 'INS': 'insertion', 'INV': 'inversion', 'CPX': 'complex rearrangement'}[t]
        if t in ('DEL', 'DUP') and n >= 50:
            vclass = f"CNV: {tname} ({n:,} bp)"
        elif t in ('DEL', 'DUP', 'INS') and 0 < n < 50:
            vclass = f"Indel ({'-' if t == 'DEL' else '+'}{n} bp, typed as structural)"
        elif t == 'INS':
            vclass = f"Structural: insertion ({n:,} bp)" if n else 'Structural: insertion (size not catalogued)'
        else:
            vclass = f"Structural: {tname} ({n:,} bp)"
        if s.get('flanks'):
            method = M_D
        elif sid in RETRO:
            method = M_RT
        elif not sv_confirmable(s):
            method = M_JU
        elif b['depth_used'] != '0':
            method = M_JD
        else:
            method = M_J
        allele = f"{t} {s['chrom']}:{int(s['start']):,}-{int(s['end']):,}"
        cohort = b; rclass = b['report_class']; extra = b['documented_breeds']
    else:
        h = (e.get('hgvs_c') or '')
        vclass = 'Indel' if ('del' in h or 'dup' in h or 'ins' in h or 'dup' in str(e.get('variant_type'))) else 'SNV'
        method = M_IMP if sid == 'cat353' else M_NONE
        rclass = 'report'; extra = ''
    status = ('Dropped (not shown)' if vkey(e) in DROPPED else 'Reported (earlier call kept, held from panel v3)' if sid in HELD
              else 'Panel only (not shown)' if rclass == 'panel_only' else 'Listed, not genotyped' if method == M_NONE else 'Reported')
    if sid in HELD:
        method = 'Held from panel v3: earlier per-kit call kept'
    carriers = sum(int(cohort.get(k) or 0) for k in ('carrier', 'two_copies', 'hemizygous')) if cohort else ''
    rows.append(dict(
        condition=name, category=CATN.get(cat, cat), gene=str(e.get('gene') or '').strip(),
        location=f"{e.get('chrom')}:{e.get('pos')}" if e.get('pos') else '', allele=allele,
        hgvs=' '.join(x for x in (e.get('hgvs_c'), e.get('hgvs_p')) if x), vclass=vclass,
        vgroup=('SNV' if vclass.startswith('SNV') else 'Mitochondrial' if vclass.startswith('Mito') else 'Indel'
                if vclass.startswith('Indel') else 'CNV' if vclass.startswith('CNV') else 'Structural (other)'),
        method=method, inheritance=MODE.get(mode, ''), status=status,
        breeds=(extra or e.get('variant_breed') or '').replace('_', ' ').title().replace(';', ', ') if (extra or e.get('variant_breed')) else '',
        c_copy=carriers, c_two=cohort.get('two_copies', '') if cohort else '', c_unc=cohort.get('uncertain', '') if cohort else '',
        c_clear=cohort.get('clear', '') if cohort else '', site_id=sid,
        ref=e.get('reference') or e.get('references') or ''))

# A condition sits in the category most of its variants carry: the app's public list counts Brown colour (bc) and
# Colour dilution (d1) as health conditions because one variant in each lacks its coat tag (variantCategories.json).
maj = collections.defaultdict(collections.Counter)
for r in rows: maj[r['condition']][r['category']] += 1
for r in rows: r['category'] = maj[r['condition']].most_common(1)[0][0]

# ── traits ──
T_IMP = 'Imputed (Dog10K panel, GP>=0.8), else direct reads'
T_RDS = 'Direct reads (pileup, MQ/BQ>=20)'
PANEL = {'chr32:35494497', 'chr32:35475211', 'chr27:44113063', 'chr1:24541931', 'chr2:74142969', 'chr15:30145498',
         'chr18:13147324', 'chr20:56252402', 'chr21:10955561', 'chr32:34768246', 'chr1:54742499', 'chr16:20112737',
         'chr16:20112974', 'chr15:41513523', 'chr3:42234279', 'chr4:67960775', 'chr4:67960816', 'chr10:49679610',
         'chr3:32515468', 'chr4:40070215', 'chr10:8703415'}   # checked against phased_plus_disease_rh.bcf 2026-11-07
GROUP = {'coat': 'Coat type & pigment', 'body': 'Body features', 'size': 'Body size genes', 'performance': 'Performance'}
trows = []
for t in TR['traits']:
    for l in t['loci']:
        ty = l['type']
        if ty == 'snp':
            loc = f"{l['chrom']}:{l['pos']}"; vc = 'SNV'; al = f"{l['ref']}>{l['alt']}"
            m = T_IMP if loc in PANEL else T_RDS
        elif ty in ('del', 'ins'):
            loc = f"{l['chrom']}:{l['pos']}"; n = l['len']; al = f"{'-' if ty == 'del' else '+'}{n} bp"
            vc = f"Indel ({al})" if n < 50 else f"Structural: insertion ({n} bp)"
            m = ('Direct reads: deletion in the read alignment at the site' if ty == 'del' else
                 'Direct reads: insertion in the alignment or reads soft-clipped at the breakpoint, vs reads spanning it')
        else:
            loc = f"{l['chrom']}:{l['start']}-{l['end']}"; n = l['end'] - l['start']
            vc = f"CNV: {'deletion' if ty == 'cnv_del' else 'duplication'} ({n:,} bp)"; al = ''
            m = 'Read depth over the event vs flanks (copy ratio 1.0 / 1.5 / 0.5 / 0)'
        trows.append(dict(trait=t['name'], group=GROUP.get(t['group'], t['group']), gene=t['gene'], locus=l.get('id'),
                          location=loc, allele=al, vclass=vc, method=m))

# ── coat colour (stage 13) ──
coat = [
    ('E (extension)', 'MC1R', 'e: p.Arg306* stop', 'chr5:64186728 G>A', 'SNV', 'Direct reads (not on the panel)'),
    ('E (extension)', 'MC1R', 'Em mask / e background: M264V haplotype', 'chr5:64186854 C>T', 'SNV', 'Imputed (Dog10K panel)'),
    ('K (dominant black)', 'CBD103', 'KB: 3-bp deletion', 'chr16:55468988-55468990', 'Indel (-3 bp)', 'Direct reads spanning the deletion'),
    ('K (dominant black)', 'CBD103', 'KB haplotype tags (fallback only)', 'chr16:55464672, 55464673', 'SNV', 'Imputed; evidence only, never a genotype on their own'),
    ('A (agouti)', 'ASIP', 'at / a coding tags', 'chr24:23906214, 23908000', 'SNV', 'Imputed (23906214) / direct reads (23908000)'),
    ('A (agouti)', 'ASIP', 'ay (sable), aw (wild type)', 'regulatory / structural', 'Structural', 'Not determined from DNA'),
    ('B (brown)', 'TYRP1', 'bs', 'chr11:33385200 C>T', 'SNV', 'Direct reads'),
    ('B (brown)', 'TYRP1', 'bd 3-bp deletion', 'chr11:33385242-33385244', 'Indel (-3 bp)', 'Direct reads'),
    ('B (brown)', 'TYRP1', 'bc p.Cys41Ser', 'chr11:33376317 T>A', 'SNV', 'Imputed, or direct reads at >=5x; brown copies = bs + max(bd, bc)'),
    ('D (dilute)', 'MLPH', 'd1 splice c.123+1G>A + tag SNP', 'chr25:48403161, 48431759', 'SNV', 'Imputed (Dog10K panel)'),
    ('M (merle)', 'PMEL', 'SINE insertion', 'chr10:644511 junction', 'Structural: insertion (~250-280 bp)', 'Direct reads: clipped vs spanning at the junction (merle class/length not determined)'),
    ('H (harlequin)', 'PSMB7', 'p.Val49Gly', 'chr9:58614853 T>G', 'SNV', 'Direct reads (too rare for imputation); shown only with merle'),
    ('S (spotting)', 'MITF', 'piebald tag', 'chr20:5711695', 'SNV', 'Tag SNP only; causal MITF insertion not typed (limited)'),
    ('W (extreme white)', 'KIT', 'structural', '', 'Structural', 'Not determined from DNA'),
]

# ── genome-wide scans ──
scans = [
    ('Predicted gene-disrupting changes', '2,400+ disease-linked genes', 'SNV + indel',
     'Dog10K site annotation (snpEff canFam4_NCBI_AR106) looked up for every called variant (stage 10); listed only when two copies, allele frequency <1%, gene recessive/X-linked/unstated and never homozygous in Dog10K. Predictions, not counted as conditions.'),
    ('Large DNA changes (CNV scan)', 'Whole genome, 50-kb windows', 'CNV',
     'Read depth vs a 93-dog panel of normals; gains/losses and homozygous gene deletions labelled with canFam4 genes.'),
    ('Sex chromosomes / karyotype', 'X, Y and autosomes', 'Aneuploidy',
     'Chromosome-level depth (e.g. trisomy X); never rendered as CNV gains.'),
    ('AMY2B amylase copy number', 'AMY2B 15-kb unit', 'CNV',
     'Depth over the unit vs GC-matched chr6 two-copy windows (stage 13b2).'),
    ('Immune-gene (DLA class II) diversity', 'DLA-DRB1/DQA1/DQB1 region', 'SNV haplotypes',
     'Imputed SNP heterozygosity across the MHC region (stage 12b).'),
    ('Maternal line', 'Mitochondrial genome', 'mtDNA SNVs', 'chrM-only remap, haplogroup call.'),
    ('Mitochondrial disease screen', 'Sensory ataxic neuropathy (MT-TY m.5304del = canFam4 chrM:11072 TA>T)', 'mtDNA indel',
     'Share of chrM reads carrying the deletion (heteroplasmy %), from mito_result.json disease_screen; shown as a card in Genetic Health, separate from the conditions list (where panel v3 holds the same site). CYTB m.14474G>A excluded: a lineage polymorphism in ~7-11% of healthy dogs.'),
    ('Paternal line', 'Y chromosome', 'Y SNVs', 'canFam4_plusY reference, Y haplogroup (new kits).'),
]

# ── workbook ──
wb = Workbook()
H = Font(bold=True, color='FFFFFF'); HF = PatternFill('solid', fgColor='3540CA'); B = Font(bold=True)
def sheet(ws, header, data, widths):
    ws.append(header)
    for c in ws[1]: c.font, c.fill = H, HF
    for d in data: ws.append(d)
    for i, w in enumerate(widths, 1): ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = 'A2'; ws.auto_filter.ref = ws.dimensions
    for r in ws.iter_rows(min_row=2):
        for c in r: c.alignment = Alignment(vertical='top', wrap_text=False)

shown = [r for r in rows if r['status'].startswith('Reported')]
conds = collections.defaultdict(set)
for r in shown: conds[r['category']].add(r['condition'])
ws = wb.active; ws.title = 'Summary'
ws.append(['ProsperK9 report catalogue', '', f"built from omia_variants_v3.json + trait catalogue {TR['version']}"]); ws['A1'].font = Font(bold=True, size=14)
ws.append([])
ws.append(['Conditions and variants shown in reports']); ws['A3'].font = B
ws.append(['Category', 'Conditions', 'Variants'])
for cat in ('Health', 'Coat colour', 'Trait'):
    ws.append([cat, len(conds[cat]), sum(r['category'] == cat for r in shown)])
ws.append(['Total', sum(len(v) for v in conds.values()), len(shown)])
ws.append(['Catalogue entries not shown', '', len(rows) - len(shown)])
ws.append([])
ws.append(['Reported variants by type and calling method']); ws[f'A{ws.max_row}'].font = B
meths = [m for m, _ in collections.Counter(r['method'] for r in shown).most_common()]
groups = ['SNV', 'Indel', 'CNV', 'Structural (other)', 'Mitochondrial']
ws.append(['Calling method'] + groups + ['Total'])
for m in meths:
    c = collections.Counter(r['vgroup'] for r in shown if r['method'] == m)
    ws.append([m] + [c.get(g, 0) for g in groups] + [sum(c.values())])
c = collections.Counter(r['vgroup'] for r in shown)
ws.append(['Total'] + [c.get(g, 0) for g in groups] + [len(shown)])
ws.append([])
ws.append(['Plus', '25 single-gene traits (sheet Traits)', '9 coat-colour loci (sheet Coat colour)', 'genome-wide scans (sheet Genome-wide)'])
ws.column_dimensions['A'].width = 62
for col in 'BCDEFG': ws.column_dimensions[col].width = 16

sheet(wb.create_sheet('Conditions (variants)'),
      ['Condition', 'Category', 'Gene', 'Location (canFam4)', 'Allele', 'HGVS', 'Variant type', 'Type group', 'Calling method',
       'Inheritance', 'Status', 'Breeds where described', 'Cohort: dogs with a copy', 'Cohort: two copies',
       'Cohort: uncertain', 'Cohort: clear', 'Site id', 'Reference'],
      [[r[k] for k in ('condition', 'category', 'gene', 'location', 'allele', 'hgvs', 'vclass', 'vgroup', 'method', 'inheritance',
                       'status', 'breeds', 'c_copy', 'c_two', 'c_unc', 'c_clear', 'site_id', 'ref')]
       for r in sorted(rows, key=lambda r: (not r['status'].startswith('Reported'), r['category'] != 'Health', r['condition'].lower()))],
      [44, 11, 12, 20, 18, 26, 30, 14, 50, 18, 20, 30, 12, 10, 10, 10, 16, 40])
sheet(wb.create_sheet('Traits (25)'), ['Trait', 'Group', 'Gene', 'Locus', 'Location (canFam4)', 'Allele', 'Variant type', 'Calling method'],
      [[r[k] for k in ('trait', 'group', 'gene', 'locus', 'location', 'allele', 'vclass', 'method')] for r in trows],
      [30, 20, 18, 34, 24, 12, 28, 60])
sheet(wb.create_sheet('Coat colour'), ['Locus', 'Gene', 'Allele / marker', 'Location (canFam4)', 'Variant type', 'Calling method'],
      coat, [20, 10, 40, 26, 30, 70])
sheet(wb.create_sheet('Genome-wide'), ['Result', 'Scope', 'Variant type', 'Method'], scans, [36, 30, 14, 110])
sheet(wb.create_sheet('Methods'), ['Calling method', 'What it means'], [
    (M_IMP, 'GLIMPSE2 imputes the genotype from the Dog10K reference haplotypes using the low-pass reads around the site; the genotype probability is combined with any reads covering the site. Reliable at 1-2x.'),
    (M_RDS, 'The site is not on the imputation panel. Reads covering it are counted (allowing for sequencing error) and combined with a prior from how common the variant is in the dog\'s breeds. Needs reads at the site: at ~1x many dogs stay "uncertain".'),
    (M_HAP, 'Disease indels not on the panel but carried on recognisable Dog10K haplotypes: the haplotype is imputed from nearby SNPs (GLIMPSE2) and multiplied by the reads likelihood. MDR1, GH1, PDK4, STRN, FGF5.'),
    (M_MT, 'Allele fraction among mitochondrial reads (heteroplasmy-aware).'),
    (M_JD, 'Reads crossing the breakpoints (split, soft-clipped, discordant pairs) plus depth inside vs outside the event (DEL/DUP >= 2 kb), with shared depth biases removed.'),
    (M_J, 'Reads crossing the breakpoints only (insertions, inversions, short events).'),
    (M_JU, 'Insertions of unknown sequence or events catalogued without a size: junction reads cannot confirm a copy, so the report only says clear or uncertain (p capped at 50%).'),
    (M_D, 'CNVs with no published breakpoints: depth over a curated interval vs distant two-copy reference windows.'),
    (M_RT, 'FGF4 retrogene insertions (CDDY, CDPA): read pairs at the insertion site whose mates map to the parent FGF4 gene.'),
    ('Call thresholds', 'clear if P(copy) <= 0.2%; carrier if >= 95%; two copies if P(2) >= 90%; otherwise uncertain. Frozen cohort fit, so new kits are scored the same way.'),
], [62, 140])
wb.save(OUT)
print('rows', len(rows), 'reported', len(shown), {k: len(v) for k, v in conds.items()}, 'traits loci', len(trows))
print(collections.Counter(r['status'] for r in rows))
for m in meths: print(f"{sum(r['method']==m for r in shown):4d}", m, dict(collections.Counter(r['vgroup'] for r in shown if r['method'] == m)))
print(dict(collections.Counter(r['vgroup'] for r in shown)))
