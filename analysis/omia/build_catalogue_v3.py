#!/usr/bin/env python3
"""Panel v3 catalogue: the live catalogue (reference_json/omia_variants.json) plus the new panel v3 sites that the
cohort model can call, as reference_json/omia_variants_v3.json.

Every entry gains a site_id: catNNN is the live catalogue's index (unchanged order, so old reports still line up
entry by entry), new_* are the panel v3 additions. A new site enters the catalogue only when
  * it has calls (it is in bp_sites.tsv or bp_sv_sites.tsv, i.e. it has coordinates and alleles),
  * it is reportable (report_class report or trait; T5 risk associations stay panel-only),
  * it is not held (held sites have no trustworthy call),
  * it is not a live catalogue site under another name (same chrom, pos and gene).
A new site in a gene whose live entries are all coat or trait variants takes that category (FGF5 long hair, TYR
Himalayan), so a new allele does not turn an existing coat or trait condition into a disease.
Live entries in DROPPED leave the screen (their catNNN ids stay, so the remaining ids keep the live index).
Sites in MANUAL_HOLD are held after review even though the model calls them. Names, modes, categories and
plain-language descriptions for the new sites come from known_variants/panel_v3/new_site_text.tsv (by gene; modes
checked against the source papers, since OMIA lists none for most of these recent conditions).
Metadata (breeds, OMIA ids, PubMed ids, HGVS) comes from known_variants/catalogue_expansion_candidates.tsv, joined by
OMIA variant id (the site_id suffix), else by gene and canFam4 position, else by gene.

  python3 analysis/omia/build_catalogue_v3.py [WORK_DIR]
"""
import csv, json, pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
P = ROOT / 'known_variants/panel_v3'
W = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'work_panel_v3'
DOG_KEYS = ('cosmo', 'nelk')
# Held after review: called by the model, but the call cannot be trusted yet.
MANUAL_HOLD = {
    # chrM 1-bp deletion: high alt fractions could be heteroplasmy or nuclear (NUMT) reads; one kit came out at risk
    # (90% of 48 reads) and another 50% uncertain (2026-11-05).
    'new_MT-TY_11072',
}

# Removed from the screen: not inherited, so a cheek-swab genotype says little (the app hides them on older reports).
DROPPED = {
    'cat275',  # JAK2 V617F polycythemia: arises in the bone marrow during life (2026-11-05)
}

def tsv(path):
    return list(csv.DictReader(open(path), delimiter='\t'))

def main():
    live = json.load(open(ROOT / 'reference_json/omia_variants.json'))
    sites = {r['site_id']: r for r in tsv(P / 'panel_v3_sites.tsv')}
    called, held = set(), set()
    for f in ('bp_sites.tsv', 'bp_sv_sites.tsv'):
        for r in tsv(W / f):
            called.add(r['site_id'])
            if r['held'] != '0' and int(r['held']) >= int(r.get('kits_with_reads') or r.get('kits_with_evidence') or 1):
                held.add(r['site_id'])
    held |= MANUAL_HOLD & called
    text = {r['gene']: r for r in tsv(P / 'new_site_text.tsv')}
    cand = tsv(ROOT / 'known_variants/catalogue_expansion_candidates.tsv')
    by_omia = {c['omia_id']: c for c in cand if c['omia_id']}
    by_pos = {(c['gene'], c['canfam4']): c for c in cand if c['canfam4']}
    by_gene = {}
    for c in cand:
        by_gene.setdefault(c['gene'], c)

    cats = json.load(open(ROOT / 'dogs-app/lib/variantCategories.json'))
    live_sites, gene_cats = set(), {}
    for v in live['variants']:
        live_sites.add((v['chrom'], int(v['pos']), v['gene'].strip().upper()))
        gene_cats.setdefault(v['gene'].strip().upper(), set()).add(cats.get(f"{v['chrom']}:{v['pos']}", v.get('category', 'disease')))
    out = []
    for i, v in enumerate(live['variants']):
        e = {k: x for k, x in v.items() if k not in DOG_KEYS}
        e['site_id'] = f'cat{i:03d}'
        if e['site_id'] not in DROPPED:
            out.append(e)
    n_new = 0
    for sid, s in sites.items():
        if sid.startswith('cat') or sid not in called or sid in held or s['report_class'] not in ('report', 'trait'):
            continue
        if (s['chrom'], int(s['pos']), s['gene'].strip().upper()) in live_sites:
            continue
        gc = gene_cats.get(s['gene'].strip().upper(), set())
        category = next(iter(gc)) if len(gc) == 1 and gc != {'disease'} else 'disease' if s['report_class'] == 'report' else 'trait'
        suffix = sid.rsplit('_', 1)[-1]
        c = (by_omia.get(suffix) or by_pos.get((s['gene'], f"{s['chrom']}:{s['pos']}"))
             or by_gene.get(s['gene']) or {})
        sv = s['kind'] == 'sv'
        e = {'site_id': sid, 'gene': s['gene'], 'chrom': s['chrom'], 'pos': int(s['pos']),
             'ref': s['ref'] if not sv else None, 'alt': s['alt'] if not sv else None,
             'variant_type': s['alt'] if sv else f"{s['ref']}>{s['alt']}",
             'trait': s['condition'], 'phene_name': s['condition'],
             'variant_breed': re.sub(r'\s*\(Dog\)', '', c.get('breeds', '')),
             'category': category,
             'source': f"panel_v3 {s['origin']}"}
        t = text.get(s['gene'], {})
        if t.get('name'):
            e['trait'] = e['phene_name'] = t['name']
        for f in ('mode', 'category', 'description'):
            if t.get(f):
                e[f] = t[f]
        if sv:
            e['pos_end'] = int(s['end'])
        hg = c.get('variant', '')
        m = re.search(r'(c\.[^\s]+)', hg)
        if m:
            e['hgvs_c'] = m.group(1)
        m = re.search(r'(p\.[^\s]+)', hg)
        if m:
            e['hgvs_p'] = m.group(1)
        if c.get('omia_id'):
            e['omia_id'] = c['omia_id']
        refs = [f"https://pubmed.ncbi.nlm.nih.gov/{p}/" for p in re.findall(r'\d{6,9}', c.get('pmid', '') + ' ' + s['refs'])]
        if refs:
            e['reference'] = refs[0]
        out.append(e)
        n_new += 1
    path = ROOT / 'reference_json/omia_variants_v3.json'
    json.dump({'method': 'Panel v3 catalogue: live catalogue (site_id catNNN = live index) + panel v3 additions (new_*)',
               'held': sorted(held), 'dropped': sorted(DROPPED), 'variants': out}, open(path, 'w'), indent=1)
    print(f'{len(live["variants"]) - len(DROPPED)} live + {n_new} new = {len(out)} entries -> {path}; held: {sorted(held)}')

if __name__ == '__main__':
    main()
