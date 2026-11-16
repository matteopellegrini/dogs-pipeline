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

Each genotype also carries p_two, the chance of two copies (for a male's X, of the one copy carrying it): what a
recessive condition turns on, so the app flags an undecided recessive site only when p_two is meaningful.

Each panel v3 genotype has source 'panel_v3', p_carrier, and basis: 'imputed' (GLIMPSE posterior at a Dog10K site),
'haplotype' (a disease indel imputed from its Dog10K haplotype, combined with the dog's reads at it),
'reads' (the dog's own reads at the site, combined with the breed prior), 'breed_prior' (no reads at the site, so
the call rests on the dog's breed mix), 'allele_fraction' (mitochondrial) or 'structural' (breakpoint reads, read
pairs and depth). Held sites and sites the model does not call keep the previous call.

Sites in NOT_DETERMINED are shown as not determined (zygosity 'not_callable') whatever the model or an older stage
says: their short-read calls are known to miss carriers, so a 'clear' there would be false reassurance. --hold
rewrites only those entries in existing reports (no cohort model needed).

The CBD103 KB deletion (cat218) is the exception: it sits at the end of a GGG run, where the pileup counts reads that
stop inside the run as reference (Cosmo at 18x: 16 alt + 2 ref against 18 of 18 spanning reads deleted), so its fitted
model never calls clear and rarely two copies. Its genotype comes from the coat stage's K locus call instead, which
counts only reads spanning the deletion (k_locus_genotype; pass the dog's coat_color.json as COAT_JSON).

  python3 panel_v3_apply.py WORK_DIR CATALOGUE SAMPLE OLD_OMIA|- OUT_OMIA [COAT_JSON]
  python3 panel_v3_apply.py --hold OMIA_JSON...
  python3 panel_v3_apply.py --p2 WORK_DIR LIST   # LIST: sample<TAB>omia_result.json; adds p_two to existing calls
"""
import csv, json, sys
import numpy as np

CALL_LO = 0.002  # must match panel_v3_breed_prior.py
ZYG = {'clear': 'ref/ref', 'carrier': 'het', 'two_copies': 'alt/alt', 'hemizygous': 'alt/alt', 'uncertain': 'uncertain'}
DOG_KEYS = ('cosmo', 'nelk')
# Calls that miss carriers: the cohort and 400+ public genomes show none, while Donner et al. 2023 (1,054,293 dogs)
# found carriers at frequencies that predict 11-127 among our customers (2026-11-15). Held until a caller that reads
# the variant's own sequence is validated on carrier-breed genomes.
NOT_DETERMINED = {
    'cat073': 'RPGRIP1 cord1: 44-bp insertion carrying an A29 run; aligners soft-clip or drop the insertion reads',
    'cat147': 'PNPLA1: complex deletion-insertion in a TAC repeat; the pileup does not see the inserted sequence',
    'cat010': 'ABCA4: 1-bp duplication in a C run; short-read pileups miss it',
    'cat446': 'KLKB1: no carriers called where ~0.6% of mixed-breed dogs carry it; catalogue allele under review',
}

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

KB_SITE = 'cat218'

def held_genotype(sid):
    return {'zygosity': 'not_callable', 'affected': None, 'call_confidence': 'low', 'source': 'panel_v3',
            'basis': 'held', 'held_reason': NOT_DETERMINED[sid]}

def gt_key(v):
    keys = [k for k, x in v.items() if is_gt(x) and k not in DOG_KEYS]
    return keys[-1] if keys else None

def hold(paths):
    """Rewrite the NOT_DETERMINED entries of existing reports in place; prints the calls it replaced."""
    for path in paths:
        res = json.load(open(path))
        changed = []
        for v in res['variants']:
            sid, k = v.get('site_id'), gt_key(v)
            if sid in NOT_DETERMINED and k and v[k].get('basis') != 'held':
                changed.append(f"{sid}:{v[k].get('zygosity')}")
                v[k] = held_genotype(sid)
        if changed:
            json.dump(res, open(path, 'w'), indent=1)
        print(path, ' '.join(changed) or 'unchanged')

def add_p2(work, lst):
    """Add p_two to the panel v3 calls of existing reports, in place, from a cohort run that saved it (bp_*.npz 'p2').
    Only where the report's p_carrier is the run's, so a call the run did not make (K locus, sequence caller, held)
    is left alone; no call changes."""
    p2 = {}
    for f in ('bp_p_carrier.npz', 'bp_sv.npz'):
        z = dict(np.load(f'{work}/{f}'))
        col = {str(x): j for j, x in enumerate(z['samples'])}
        p2[f] = (z, col, {str(x): i for i, x in enumerate(z['sid'])})
    for line in open(lst):
        sample, path = line.rstrip('\n').split('\t')[:2]
        res = json.load(open(path))
        n = {'added': 0, 'skipped': 0}
        for v in res['variants']:
            sid, k = v.get('site_id'), gt_key(v)
            g = v.get(k) if k else None
            if not g or g.get('source') != 'panel_v3' or g.get('p_carrier') is None or 'p_two' in g:
                continue
            for z, col, row in p2.values():
                if sample in col and sid in row:
                    i, j = row[sid], col[sample]
                    if abs(round(float(z['p'][i, j]), 4) - g['p_carrier']) < 1e-4:
                        g['p_two'] = round(float(z['p2'][i, j]), 4); n['added'] += 1
                    else:
                        n['skipped'] += 1
                    break
        json.dump(res, open(path, 'w'), indent=1)
        print(sample, path, n['added'], n['skipped'])

def k_locus_genotype(coat, model_g):
    """The KB deletion genotype from the coat stage's K locus read call (stage 13), or None to keep model_g.
    KB/KB -> two copies, KB/ky -> one copy, high-confidence ky/ky (>=6 spanning reads, none deleted) -> clear;
    KB/? and lower-confidence ky/ky stay uncertain with the model's p_carrier. Tag-only calls (no spanning read)
    keep the model call."""
    k = ((coat or {}).get('loci') or {}).get('K') or {}
    ev = k.get('read_evidence') or {}
    nd, nr = int(ev.get('n_del') or 0), int(ev.get('n_ref') or 0)
    if nd + nr == 0:
        return None
    a = sorted(k.get('predicted_alleles') or [])
    conf = k.get('confidence') or 'low'
    extra = {'site_reads': nd + nr, 'k_locus_reads': {'deleted': nd, 'not_deleted': nr}}
    if a == ['KB', 'KB']:
        call, p = 'two_copies', 1.0
    elif a == ['KB', 'ky']:
        call, p = 'carrier', 1.0
    elif a == ['ky', 'ky'] and conf == 'high':
        call, p = 'clear', 0.0
    else:
        call, p = 'uncertain', float((model_g or {}).get('p_carrier', 0.5))
        conf = 'low'
    g = genotype(call, p, 'reads', extra)
    g['call_confidence'] = conf
    return g

def main(work, catalogue, sample, old_path, out_path, coat_path=None):
    cat = json.load(open(catalogue))
    # dict(): an NpzFile decompresses an array again on every item access, i.e. per site
    small = dict(np.load(f'{work}/bp_p_carrier.npz'))
    sv = dict(np.load(f'{work}/bp_sv.npz'))
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
    p2s = small.get('p2')
    for i, sid in enumerate(small['sid']):
        c, p, n = str(small['call'][i, j]), float(small['p'][i, j]), int(small['reads'][i, j])
        src = source.get(str(sid), '')
        basis = ('imputed' if src.startswith('GLIMPSE') else 'haplotype' if src.startswith('Haplotype')
                 else 'allele_fraction' if src.startswith('chrM') else 'reads' if n > 0 else 'breed_prior')
        extra = {'site_reads': n}
        if p2s is not None:
            extra['p_two'] = round(float(p2s[i, j]), 4)
        if basis in ('imputed', 'haplotype') and site_af.get(str(sid)) is not None:
            extra['af_dog10k'] = round(site_af[str(sid)], 4)
        calls[str(sid)] = (c, p, basis, extra)
    svj = {s: i for i, s in enumerate(sv['samples'])}[sample]
    ev = {}
    for r in csv.reader(open(f'{work}/sv_counts.tsv'), delimiter='\t'):
        if r[0] == sample:
            ev[r[1]] = {'breakpoint_alt_reads': int(r[2]), 'breakpoint_ref_reads': int(r[3]), 'pair_reads': int(r[4])}
    for i, sid in enumerate(sv['sid']):
        e = dict(ev.get(str(sid), {}))
        if 'p2' in sv:
            e['p_two'] = round(float(sv['p2'][i, svj]), 4)
        calls[str(sid)] = (str(sv['call'][i, svj]), float(sv['p'][i, svj]), 'structural', e)

    old = json.load(open(old_path)) if old_path != '-' else None
    try:
        coat = json.load(open(coat_path)) if coat_path else None
    except (OSError, ValueError):
        coat = None
    old_v = old['variants'] if old else []
    key = next((k for v in old_v for k, x in v.items() if is_gt(x) and k not in DOG_KEYS), sample.lower())
    variants, n = [], {'panel_v3': 0, 'kept': 0, 'none': 0}
    for e in cat['variants']:
        sid = e['site_id']
        base = dict(e)
        # catNNN is the live catalogue index, which is also the entry's position in a report written from it.
        k = int(sid[3:]) if sid.startswith('cat') else -1
        if 0 <= k < len(old_v) and old_v[k].get('gene') == e['gene'] and str(old_v[k].get('pos')) == str(e['pos']):
            base = dict(old_v[k]); base['site_id'] = sid
        c = calls.get(sid)
        if sid in NOT_DETERMINED:
            base[key] = held_genotype(sid)
            n['none'] += 1
        elif c and c[0] != 'held':
            base[key] = genotype(*c)
            if sid == KB_SITE and coat is not None:
                base[key] = k_locus_genotype(coat, base[key]) or base[key]
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
    if sys.argv[1] == '--hold':
        hold(sys.argv[2:])
    elif sys.argv[1] == '--p2':
        add_p2(*sys.argv[2:4])
    else:
        main(*sys.argv[1:7])
