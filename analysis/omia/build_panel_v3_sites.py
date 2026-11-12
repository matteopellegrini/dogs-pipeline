"""Freeze the panel v3 site list: the current catalogue plus the curated expansion candidates (2026-11-01).

Inputs (known_variants/):
  ../reference_json/omia_variants.json          current catalogue (474 variants)
  catalogue_expansion_candidates.tsv            merged five-source table with tiers
  panel_v3/t1_t2_norm.vcf                       T1/T2 small variants, canFam4, left-normalised (ids GENE|omia_id)
  panel_v3/t4_priority_curated.tsv              T4 priority rows curated from the papers
  panel_v3/t5_curated.tsv                       T5 (panel-only risk) rows with coordinates
  panel_v3/t3_t6_sv_coords.tsv                  T3 structural + T6 trait rows (coords or gene span)
  panel_v3/depth_only_cnv.tsv                   CNVs without published breakpoints, genotyped from depth alone
                                                (canFam4 interval + distant reference windows; analysis/cnv catalogue)
Outputs (known_variants/panel_v3/):
  panel_v3_sites.tsv        every site, with report class and genotyping route
  panel_v3_small.vcf.gz     small variants with resolved alleles (imputation / direct-call targets)
Run from the pipeline root: python3 analysis/omia/build_panel_v3_sites.py
"""
import csv, json, re, subprocess, collections, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]
KV = ROOT / 'known_variants'; P = KV / 'panel_v3'; FA = ROOT / 'canFam4.fa'
COLS = ['site_id', 'origin', 'report_class', 'gene', 'condition', 'chrom', 'pos', 'end', 'ref', 'alt', 'kind',
        'genotyping', 'status', 'refs', 'note']
rows = []

def kind_of(ref, alt):
    if not re.fullmatch(r'[ACGTN]+', ref or '') or not re.fullmatch(r'[ACGTN,]+', alt or ''): return 'sv'
    if max(len(ref), len(alt)) > 50: return 'sv'
    return 'snv' if len(ref) == len(alt) == 1 else 'indel'

def add(**k):
    k.setdefault('end', k.get('pos')); rows.append({c: k.get(c, '') for c in COLS})

# 1. current catalogue; 217 entries carry no alleles, 152 of them resolved from OMIA g. notation
#    (analysis/omia/resolve_catalogue_alleles.py -> panel_v3/catalogue_alleles_resolved.tsv)
RES = {r['site_id']: r for r in csv.DictReader(open(P / 'catalogue_alleles_resolved.tsv'), delimiter='\t')}
# the other 65 curated by hand from OMIA/papers (analysis/omia/curate_catalogue_alleles.py): 22 small, 43 SV
CUR = {r['site_id']: r for r in csv.DictReader(open(P / 'catalogue_alleles_curated.tsv'), delimiter='\t')}
CUR_OMIA = {m for r in CUR.values() for m in re.findall(r'OMIA (\d+)', r['source'])}
for i, v in enumerate(json.load(open(ROOT / 'reference_json/omia_variants.json'))['variants']):
    sid = f'cat{i:03d}'; ref, alt = str(v['ref']).upper(), str(v['alt']).upper(); pos = v['pos']; status = 'live'; chrom = v['chrom']
    if ref == 'NONE':
        r = RES.get(sid, {})
        cu = CUR.get(sid)
        if r.get('ref'): pos, ref, alt, status = r['pos4'], r['ref'], r['alt'], f"live; alleles from OMIA variant {r['omia_variant_id']}"
        elif cu:
            chrom = cu['chrom']; pos, ref, alt = cu['pos'], cu['ref'], cu['alt']
            status = f"live; curated 2026-11-01 from {cu['source']}" + (f"; {cu['note']}" if cu['note'] else '')
            if cu['kind'] == 'sv':
                add(site_id=sid, origin='catalogue', report_class={'coat': 'trait', 'trait': 'trait'}.get(v.get('category'), 'report'),
                    gene=v['gene'], condition=v.get('trait', ''), chrom=chrom, pos=pos, end=cu['end'], ref='.', alt=alt, kind='sv',
                    genotyping='read_based', status=status + f"; SVLEN={cu['svlen']}", refs=v.get('reference', ''))
                continue
        else: ref, alt, status = '.', v.get('variant_type') or 'SV', f"live; needs alleles ({r.get('omia_g_or_reason', '')})"
    kd = kind_of(ref, alt)
    cat = v.get('category', 'disease')
    add(site_id=sid, origin='catalogue', report_class={'coat': 'trait', 'trait': 'trait'}.get(cat, 'report'),
        gene=v['gene'], condition=v.get('trait', ''), chrom=chrom, pos=pos, ref=ref, alt=alt, kind=kd,
        genotyping=('needs_alleles' if ref == '.' else 'read_based') if kd == 'sv' else 'imputable_or_direct',
        status=status, refs=v.get('reference', ''))

cand = list(csv.DictReader(open(KV / 'catalogue_expansion_candidates.tsv'), delimiter='\t'))
by_omia = {r['omia_id']: r for r in cand if r['omia_id']}
by_gene = {r['gene']: r for r in cand}
EXCLUDE = {'93': 'BRAF V595E is a somatic tumour mutation, not inherited'}

# 2. T1/T2 small variants
for l in open(P / 't1_t2_norm.vcf'):
    if l.startswith('#'): continue
    c, pos, vid, ref, alt = l.split('\t')[:5]
    gene, oid = vid.split('|'); r = by_omia.get(oid) or by_gene[gene]
    if oid in EXCLUDE: continue
    add(site_id=f'new_{gene}_{oid or pos}', origin=r['tier'], report_class='report', gene=gene, condition=r['condition'],
        chrom=c, pos=pos, ref=ref, alt=alt, kind=kind_of(ref, alt), genotyping='imputable_or_direct', status='verified',
        refs=r['pmid'])

# 3. T4 priority and T5 curated
for fn, origin, rc in [('t4_priority_curated.tsv', 'T4 needs coordinates (priority)', 'report'),
                       ('t5_curated.tsv', 'T5 risk/association (panel only)', 'panel_only')]:
    for r in csv.DictReader(open(P / fn), delimiter='\t'):
        kd = r['kind']; end = r['pos']
        m = re.search(r'END=(\d+)', r['note'])
        if m: end = m.group(1)
        add(site_id=f"new_{r['gene']}_{r['pos']}", origin=origin, report_class=rc, gene=r['gene'], condition=r['condition'],
            chrom=r['chrom'], pos=r['pos'], end=end, ref=r['ref'], alt=r['alt'], kind=kd,
            genotyping='read_based' if kd == 'sv' or 'read-based' in r['note'] else
                       'pending_alleles' if r['alt'] == 'N' else 'imputable_or_direct',
            status=r['status'], refs=r['pmid'], note=r['note'])

# 4. T3 structural and T6 traits
# CNVs whose breakpoints were never published get a curated canFam4 interval instead of the gene span; they are
# genotyped from depth against distant reference windows (sv_evidence.py reads the windows from the same table).
DEPTH = {r['omia_id']: r for r in csv.DictReader(open(P / 'depth_only_cnv.tsv'), delimiter='\t')}
for r in csv.DictReader(open(P / 't3_t6_sv_coords.tsv'), delimiter='\t'):
    if r['omia_id'] in CUR_OMIA: continue  # already a catalogue entry, curated above
    t3 = r['tier'].startswith('T3')
    rc = 'trait' if (not t3 or 'merle' in r['condition'].lower()) else 'report'
    d = DEPTH.get(r['omia_id'])
    add(site_id=f"new_{r['gene']}_{r['omia_id'] or r['start']}", origin=r['tier'], report_class=rc, gene=r['gene'],
        condition=r['condition'], chrom=r['chrom'], pos=d['start'] if d else r['start'], end=d['end'] if d else r['end'],
        ref='.', alt=r['type'] or 'SV', kind='sv', genotyping='read_based',
        status=f"depth_only; {d['basis']}" if d else r['status'], refs=r['omia_id'] and f"OMIA variant {r['omia_id']}",
        note='trait loci already genotyped in stages 13/13b are not duplicated there' if not t3 else '')

# 5. T5 rows without coordinates and T4 non-priority: listed so nothing is lost, not genotyped yet
have = {(x['gene'].split('_')[0].upper(), x['origin'][:2]) for x in rows}
for r in cand:
    t = r['tier'][:2]
    g = re.split(r'[ (/;,+]', r['gene'].strip())[0].upper()
    if t == 'T5' and (g, 'T5') not in have:
        add(site_id=f"todo_{g}", origin=r['tier'], report_class='panel_only', gene=r['gene'], condition=r['condition'],
            chrom='', pos='', kind='', genotyping='needs_coordinates', status='later', refs=r['pmid'], note=r['variant'][:150])
    if r['tier'] == 'T4 needs coordinates' and r['omia_id'] in DEPTH:
        d = DEPTH[r['omia_id']]
        add(site_id=f"new_{d['gene']}_{r['omia_id']}", origin=r['tier'], report_class='report', gene=d['gene'],
            condition=r['condition'], chrom=d['chrom'], pos=d['start'], end=d['end'], ref='.',
            alt={'DEL': 'deletion, gross (>20)', 'DUP': 'duplication'}[d['type']], kind='sv', genotyping='read_based',
            status=f"depth_only; {d['basis']}", refs=f"OMIA variant {r['omia_id']}")
    elif r['tier'] == 'T4 needs coordinates':
        add(site_id=f"todo_{g}_{len(rows)}", origin=r['tier'], report_class='report', gene=r['gene'], condition=r['condition'],
            chrom='', pos='', kind='', genotyping='needs_coordinates', status='later', refs=r['pmid'], note=r['variant'][:150])

# 6. literature search of 2026-11-11 (~/Projects/dog literature, results/canfam4_variants.csv): abstract-mined
#     variants not in OMIA, placed on canFam4 with REF checked; report_class per row (Mendelian -> report,
#     associations -> panel_only). Rows with no alleles, an ambiguous ALT or 'not usable' confidence were left out.
#     Last, so the todo_* ids above (numbered by row) do not shift.
for r in csv.DictReader(open(P / 'literature_2026-11-11.tsv'), delimiter='\t'):
    add(site_id=f"new_{r['gene']}_{r['pos']}", origin='Literature 2026-11-11', report_class=r['report_class'],
        gene=r['gene'], condition=r['condition'], chrom=r['chrom'], pos=r['pos'], ref=r['ref'], alt=r['alt'],
        kind=kind_of(r['ref'], r['alt']), genotyping='imputable_or_direct', status=r['status'],
        refs=f"PMID {r['pmid']}", note=r['note'])

# dedupe exact small-variant duplicates (keep the catalogue copy)
seen = set(); out = []
for x in rows:
    k = (x['chrom'], str(x['pos']), x['ref'], x['alt']) if x['kind'] in ('snv', 'indel') else None
    if k and k in seen: continue
    if k: seen.add(k)
    out.append(x)
# fixes made on the frozen file 2026-11-03 (sv_evidence.py reads SVLEN from status)
FIXES = {
    'cat227': dict(note='END fixed 2026-11-03: CanFam3 g.30852988_30902901del runs to the chr37 end (30,785,915 in '
                        'canFam4); the lifted END 60261342 was past the chromosome'),
    'new_SPAG1_1922': dict(status='coords_native; SVLEN=29',
                           note='2026-11-03: OMIA gives the allele, a 29-bp GC-rich duplication (c.1190_1218dup, '
                                'insGGCGGCGGCAAGCGGCCGGAGAGGGGCGC); typed by size, not as an unknown gross insertion'),
}
for x in out:
    x.update(FIXES.get(x['site_id'], {}))
with open(P / 'panel_v3_sites.tsv', 'w') as f:
    w = csv.DictWriter(f, COLS, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(out)

# small-variant VCF
fai = [l.split('\t')[:2] for l in open(str(FA) + '.fai')]
order = {c: i for i, (c, _) in enumerate(fai)}
small = sorted([x for x in out if x['kind'] in ('snv', 'indel') and x['alt'] not in ('N', '')],
               key=lambda x: (order.get(x['chrom'], 999), int(x['pos'])))
vcf = P / 'panel_v3_small.vcf'
with open(vcf, 'w') as f:
    f.write('##fileformat=VCFv4.2\n##reference=canFam4 (UU_Cfam_GSD_1.0)\n')
    f.writelines(f'##contig=<ID={c},length={n}>\n' for c, n in fai)
    f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">\n##INFO=<ID=RC,Number=1,Type=String,Description="report_class">\n')
    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for x in small:
        f.write(f"{x['chrom']}\t{x['pos']}\t{x['site_id']}\t{x['ref']}\t{x['alt']}\t.\t.\tGENE={re.sub(r'[^A-Za-z0-9._-]', '_', x['gene'])};RC={x['report_class']}\n")
subprocess.run(f'bcftools norm -f {FA} -c x {vcf} -Oz -o {vcf}.gz && bcftools index -t {vcf}.gz && rm {vcf}', shell=True, check=True)

cnt = collections.Counter((x['origin'][:2] if x['origin'] != 'catalogue' else 'catalogue', x['report_class'], x['genotyping']) for x in out)
for k in sorted(cnt): print(*k, cnt[k], sep='\t')
print('total sites', len(out), '| small VCF records', len(small))
