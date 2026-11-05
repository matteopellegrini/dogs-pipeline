#!/usr/bin/env python3
"""Write panel v3 calls into a dog's omia_result.json (the file the Genetic Health tab reads).

The report keeps its shape: one entry per catalogue variant, the dog's genotype under its own sample key. Entries
follow reference_json/omia_variants_v3.json (live catalogue in its original order, then the panel v3 additions).
For every site the cohort model calls (work_panel_v3/bp_p_carrier.npz, bp_sv.npz), the dog's genotype is replaced:

  call         zygosity    meaning
  clear        ref/ref     P(>=1 copy) <= CALL_LO (0.2%)
  carrier      het         P(>=1 copy) >= 95%
  two_copies   alt/alt     P(2 copies) >= 90%
  hemizygous   alt/alt     a male's single X copy carries it
  uncertain    uncertain   in between: too few reads to decide; p_carrier says how likely a copy is

Each panel v3 genotype has source 'panel_v3', p_carrier, and basis: 'imputed' (GLIMPSE posterior at a Dog10K site),
'reads' (the dog's own reads at the site, combined with the breed prior), 'breed_prior' (no reads at the site, so
the call rests on the dog's breed mix), 'allele_fraction' (mitochondrial) or 'structural' (breakpoint reads, read
pairs and depth). Held sites and sites the model does not call keep the previous call.

  python3 panel_v3_apply.py WORK_DIR CATALOGUE SAMPLE OLD_OMIA|- OUT_OMIA
"""
import csv, json, sys
import numpy as np

CALL_LO = 0.002  # must match panel_v3_breed_prior.py
ZYG = {'clear': 'ref/ref', 'carrier': 'het', 'two_copies': 'alt/alt', 'hemizygous': 'alt/alt', 'uncertain': 'uncertain'}
DOG_KEYS = ('cosmo', 'nelk')

def is_gt(x):
    return isinstance(x, dict) and 'zygosity' in x

def genotype(call, p, basis, extra):
    g = {'zygosity': ZYG[call], 'affected': None if call == 'uncertain' else call != 'clear',
         'call_confidence': 'low' if call == 'uncertain' else 'high' if call == 'clear' or p >= 0.99 else 'medium',
         'source': 'panel_v3', 'p_carrier': round(float(p), 4), 'basis': basis}
    if call == 'hemizygous':
        g['hemizygous'] = True
    g.update(extra)
    return g

def main(work, catalogue, sample, old_path, out_path):
    cat = json.load(open(catalogue))
    small = np.load(f'{work}/bp_p_carrier.npz')
    sv = np.load(f'{work}/bp_sv.npz')
    col = {s: i for i, s in enumerate(small['samples'])}
    if sample not in col:
        sys.exit(f'{sample} is not in the panel v3 cohort run')
    j = col[sample]
    calls = {}
    sites = list(csv.DictReader(open(f'{work}/bp_sites.tsv'), delimiter='\t'))
    source = {r['site_id']: r['source'] for r in sites}
    # Dog10K panel frequency at the imputed sites (glimpse_sites_af.tsv: CHROM POS REF ALT AC AN from the panel BCF):
    # the report reads a variant carried by >=5% of reference dogs as common, not as a finding.
    af = {}
    for r in csv.reader(open(f'{work}/glimpse_sites_af.tsv'), delimiter='\t'):
        if int(r[5]) > 0:
            af[tuple(r[:4])] = int(r[4]) / int(r[5])
    site_af = {r['site_id']: af.get((r['chrom'], r['pos'], r['ref'], r['alt'])) for r in sites}
    for i, sid in enumerate(small['sid']):
        c, p, n = str(small['call'][i, j]), float(small['p'][i, j]), int(small['reads'][i, j])
        src = source.get(str(sid), '')
        basis = ('imputed' if src.startswith('GLIMPSE') else 'allele_fraction' if src.startswith('chrM')
                 else 'reads' if n > 0 else 'breed_prior')
        extra = {'site_reads': n}
        if basis == 'imputed' and site_af.get(str(sid)) is not None:
            extra['af_dog10k'] = round(site_af[str(sid)], 4)
        calls[str(sid)] = (c, p, basis, extra)
    svj = {s: i for i, s in enumerate(sv['samples'])}[sample]
    ev = {}
    for r in csv.reader(open(f'{work}/sv_counts.tsv'), delimiter='\t'):
        if r[0] == sample:
            ev[r[1]] = {'breakpoint_alt_reads': int(r[2]), 'breakpoint_ref_reads': int(r[3]), 'pair_reads': int(r[4])}
    for i, sid in enumerate(sv['sid']):
        calls[str(sid)] = (str(sv['call'][i, svj]), float(sv['p'][i, svj]), 'structural', ev.get(str(sid), {}))

    old = json.load(open(old_path)) if old_path != '-' else None
    old_v = old['variants'] if old else []
    key = next((k for v in old_v for k, x in v.items() if is_gt(x) and k not in DOG_KEYS), sample.lower())
    variants, n = [], {'panel_v3': 0, 'kept': 0, 'none': 0}
    for i, e in enumerate(cat['variants']):
        sid = e['site_id']
        base = dict(e)
        if sid.startswith('cat') and i < len(old_v) and old_v[i].get('gene') == e['gene'] and str(old_v[i].get('pos')) == str(e['pos']):
            base = dict(old_v[i]); base['site_id'] = sid
        c = calls.get(sid)
        if c and c[0] != 'held':
            base[key] = genotype(*c)
            n['panel_v3'] += 1
        elif is_gt(base.get(key)):
            n['kept'] += 1
        else:
            n['none'] += 1
        variants.append(base)
    res = dict(old or {})
    res['variants'] = variants
    res['panel_v3'] = {'call_lo': CALL_LO, 'catalogue_entries': len(variants), **n}
    res['method'] = ((old or {}).get('method', '') + ' Panel v3 (2026-11): each catalogue site is called from the dog\'s '
                     'own reads combined with a breed-aware prior (or the GLIMPSE posterior at Dog10K sites); structural '
                     'variants from breakpoint reads, read pairs and depth. A site is reported clear only when the chance '
                     f'of a copy is below {CALL_LO:.1%}, a copy only above 95%, and uncertain in between.').strip()
    json.dump(res, open(out_path, 'w'), indent=1)
    print(f'{sample}: {len(variants)} entries, panel_v3 {n["panel_v3"]}, kept {n["kept"]}, no call {n["none"]} -> {out_path}')

if __name__ == '__main__':
    main(*sys.argv[1:6])
