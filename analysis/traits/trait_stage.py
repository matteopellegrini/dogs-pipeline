#!/usr/bin/env python3
"""
Stage 13b: single-gene physical traits -> traits_result.json (Embark-style
"Traits" section: coat type, body features, body size genes, performance).

  trait_stage.py <trait_catalog.json> <out.json> --bcf <imputed.bcf> --bam <markdup.bam>
                 [--sex male|female]

Genotyping engines (chosen per locus by its `type`):
  snp / small indel   Dog10K-imputed GT when the site is in the panel and
                      max(GP) >= 0.8, else a BAM pileup (MQ>=20, BQ>=20);
                      indels are read from CIGAR at the position.
  ins (structural)    reads carrying an insertion of ~len at the breakpoint OR
                      soft-clipped at it, vs reads spanning it cleanly.
  cnv_dup / cnv_del   mean depth over the event vs the flanks (BAM), copy-number
                      ratio -> 1.0 normal, ~1.5 het dup, ~0.5 het del, ~0 hom del.
Each trait's interpretation is a rule in TRAIT_RULES keyed by trait id; the
catalogue carries coordinates, alleles and copy. Calls carry a confidence and
the evidence (source, reads, GP) so the report can say how sure we are.
"""
import json, sys, math
import numpy as np
import pysam

MIN_GP = 0.8
MIN_BQ, MIN_MQ = 20, 20


# ───────────────────────── genotyping engines ─────────────────────────
def imputed_gt(bcf, chrom, pos, ref, alt):
    """-> (n_alt, max_gp, gt_str) or None if site absent / different alleles"""
    if bcf is None: return None
    try:
        for rec in bcf.fetch(chrom, pos - 1, pos):
            if rec.pos != pos: continue
            if not (rec.ref == ref and alt in rec.alts):
                # allow the swapped orientation
                if rec.ref == alt and ref in rec.alts:
                    s = rec.samples[0]; gt = s.get('GT'); gp = s.get('GP')
                    if gt is None or None in gt: return None
                    n_alt = 2 - sum(1 for a in gt if a and a > 0)
                    return n_alt, (max(gp) if gp else None), '|'.join(str(a) for a in gt) + ' (swapped)'
                continue
            s = rec.samples[0]; gt = s.get('GT'); gp = s.get('GP')
            if gt is None or None in gt: return None
            n_alt = sum(1 for a in gt if a and a > 0)
            return n_alt, (max(gp) if gp else None), '|'.join(str(a) for a in gt)
    except ValueError:
        return None
    return None


def pileup_bases(bam, chrom, pos):
    counts = {}
    for col in bam.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=MIN_BQ,
                          min_mapping_quality=MIN_MQ, ignore_overlaps=True):
        if col.reference_pos != pos - 1: continue
        for r in col.pileups:
            if r.is_del or r.is_refskip or r.query_position is None: continue
            b = r.alignment.query_sequence[r.query_position].upper()
            counts[b] = counts.get(b, 0) + 1
    return counts


def cigar_indel_counts(bam, chrom, pos, kind, length, tol=None):
    """Reads with a deletion ('del') or insertion ('ins') of ~length near pos, OR soft-clipped
    inside the event window (an insertion the aligner would not gap), vs reads that span the
    whole event window [pos-5, pos+length+5] cleanly. -> (n_event, n_span)"""
    if tol is None: tol = max(3, length + 2) if kind == 'ins' else max(3, min(length, 20))
    w0, w1 = pos - 1 - tol, pos - 1 + (length if kind == 'del' else length) + tol   # 0-based window
    n_event = n_span = 0
    for r in bam.fetch(chrom, max(0, pos - 300), pos + length + 300):
        if r.mapping_quality < MIN_MQ or r.is_secondary or r.is_supplementary or r.is_duplicate: continue
        if r.reference_end <= w0 or r.reference_start >= w1: continue
        rp = r.reference_start; hit = False; first = True
        for op, ln in r.cigartuples or []:
            if op in (0, 7, 8):        # M/=/X
                rp += ln
            elif op == 2:              # D
                if kind == 'del' and abs(rp - (pos - 1)) <= tol and abs(ln - length) <= max(2, length // 4): hit = True
                rp += ln
            elif op == 1:              # I
                if kind == 'ins' and abs(rp - (pos - 1)) <= tol and abs(ln - length) <= max(2, length // 4): hit = True
            elif op == 3:              # N
                rp += ln
            elif op == 4:              # S: soft clip at the breakpoint = evidence only for insertions the aligner
                bp = r.reference_start if first else rp   # will not gap (>= 10 bp); never for deletions, where clips
                if kind == 'ins' and length >= 10 and ln >= 5 and w0 - 2 <= bp <= w1 + 2: hit = True   # come from divergent reference sequence
            first = False
        if hit: n_event += 1
        elif r.reference_start <= w0 - 3 and r.reference_end >= w1 + 3: n_span += 1
    return n_event, n_span


def region_depth(bam, chrom, start, end):
    n = bam.count(chrom, start, end, read_callback=lambda r: r.mapping_quality >= MIN_MQ and not r.is_duplicate and not r.is_secondary)
    return n / max(1, end - start)


def cnv_ratio(bam, chrom, start, end, flank=None):
    flank = flank or max(20000, end - start)
    inside = region_depth(bam, chrom, start, end)
    left = region_depth(bam, chrom, max(0, start - flank), start)
    right = region_depth(bam, chrom, end, end + flank)
    base = (left + right) / 2
    return (inside / base if base > 0 else None), inside, base


# ───────────────────────── per-locus call ─────────────────────────
def call_locus(loc, bcf, bam):
    """-> dict(n_alt, source, confidence, evidence) ; n_alt None = no call"""
    t = loc.get('type', 'snp'); chrom = loc['chrom']; pos = int(loc.get('pos') or loc.get('start'))
    if t in ('snp', 'del', 'ins') and len(loc.get('ref', '')) == 1 and len(loc.get('alt', '')) == 1 and t == 'snp':
        r = imputed_gt(bcf, chrom, pos, loc['ref'], loc['alt'])
        if r and r[1] is not None and r[1] >= MIN_GP:
            return dict(n_alt=r[0], source='Dog10K imputed', confidence='high' if r[1] >= 0.95 else 'medium', evidence={'gt': r[2], 'max_gp': round(r[1], 3)})
        if bam is not None:
            c = pileup_bases(bam, chrom, pos); nr, na = c.get(loc['ref'], 0), c.get(loc['alt'], 0); tot = nr + na
            if tot >= 3 or (tot >= 1 and not (r and r[1] is not None and r[1] >= 0.5)):
                frac = na / tot
                n_alt = 2 if frac >= 0.85 else (0 if frac <= 0.15 else 1)
                conf = 'high' if tot >= 12 else ('medium' if tot >= 5 else 'low')
                return dict(n_alt=n_alt, source=f'reads ({tot})', confidence=conf, evidence={'ref_reads': nr, 'alt_reads': na, 'other': {k: v for k, v in c.items() if k not in (loc['ref'], loc['alt'])}})
        if r and r[1] is not None and r[1] >= 0.5:   # imputed call below the confident threshold
            return dict(n_alt=r[0], source='Dog10K imputed (low GP)', confidence='low', evidence={'gt': r[2], 'max_gp': round(r[1], 3)})
        return dict(n_alt=None, source='no data', confidence='none', evidence={})
    if t in ('del', 'ins'):
        if bam is None: return dict(n_alt=None, source='no BAM', confidence='none', evidence={})
        ne, ns = cigar_indel_counts(bam, chrom, pos, t, int(loc['len']))
        tot = ne + ns
        if tot == 0: return dict(n_alt=None, source='no reads', confidence='none', evidence={})
        frac = ne / tot
        n_alt = 2 if frac >= 0.8 else (0 if frac <= 0.1 else 1)
        conf = 'high' if tot >= 12 else ('medium' if tot >= 5 else 'low')
        return dict(n_alt=n_alt, source=f'reads ({tot})', confidence=conf, evidence={'event_reads': ne, 'spanning_reads': ns})
    if t in ('cnv_dup', 'cnv_del'):
        if bam is None: return dict(n_alt=None, source='no BAM', confidence='none', evidence={})
        ratio, inside, base = cnv_ratio(bam, chrom, int(loc['start']), int(loc['end']))
        if ratio is None or base * (int(loc['end']) - int(loc['start'])) < 30:
            return dict(n_alt=None, source='insufficient coverage', confidence='none', evidence={'ratio': ratio})
        if t == 'cnv_dup':
            n_alt = 0 if ratio < 1.25 else (1 if ratio < 1.75 else 2)
        else:
            n_alt = 0 if ratio > 0.75 else (1 if ratio > 0.25 else 2)
        conf = 'high' if base * (int(loc['end']) - int(loc['start'])) >= 300 else 'medium'
        return dict(n_alt=n_alt, source='read depth', confidence=conf, evidence={'depth_ratio': round(ratio, 2), 'depth_inside': round(inside, 2), 'depth_flanks': round(base, 2)})
    return dict(n_alt=None, source='unsupported', confidence='none', evidence={})


# ───────────────────────── trait interpretation rules ─────────────────────────
def gstr(a, b): return f'{a}{b}'

def rule_recessive(calls, allele_alt, allele_ref, pheno_alt2, pheno_carrier, pheno_ref, label_locus=None):
    """single-locus recessive (two alt copies -> phenotype)"""
    c = calls[0]
    if c['n_alt'] is None: return None
    g = [allele_alt] * c['n_alt'] + [allele_ref] * (2 - c['n_alt'])
    txt = pheno_alt2 if c['n_alt'] == 2 else (pheno_carrier if c['n_alt'] == 1 else pheno_ref)
    return ''.join(sorted(g)), txt, c['confidence']

def rule_dominant(calls, allele_alt, allele_ref, pheno_any, pheno_ref):
    c = calls[0]
    if c['n_alt'] is None: return None
    g = [allele_alt] * c['n_alt'] + [allele_ref] * (2 - c['n_alt'])
    return ''.join(sorted(g)), (pheno_any if c['n_alt'] >= 1 else pheno_ref), c['confidence']

def rule_additive3(calls, allele_alt, allele_ref, p2, p1, p0):
    c = calls[0]
    if c['n_alt'] is None: return None
    g = [allele_alt] * c['n_alt'] + [allele_ref] * (2 - c['n_alt'])
    return ''.join(sorted(g)), [p0, p1, p2][c['n_alt']], c['confidence']


def interpret(trait, calls, ctx):
    tid = trait['id']; loci = trait['loci']
    known = [c for c in calls if c['n_alt'] is not None]
    if not known: return None
    conf = min((c['confidence'] for c in known), key=lambda x: ['none', 'low', 'medium', 'high'].index(x))
    if tid == 'coat_length':
        # any two long-coat alleles across the FGF5 variants -> long
        n_lh = sum(c['n_alt'] for c in known)
        g = f'{"Lh" * min(2, n_lh)}{"Sh" * max(0, 2 - min(2, n_lh))}' if n_lh <= 2 else 'LhLh'
        return g, ('Likely long coat' if n_lh >= 2 else ('Likely short coat (carries one long-coat copy)' if n_lh == 1 else 'Likely short coat')), conf
    if tid == 'furnishings':
        return rule_dominant(calls, 'F', 'I', 'Likely furnished (mustache, beard and/or eyebrows)', 'Likely unfurnished')
    if tid == 'curl':
        n = sum(c['n_alt'] for c in known); furnished = ctx.get('furnished')
        g = 'C' * (2 - min(2, n)) + 'T' * min(2, n)
        if n >= 1: txt = 'Likely wavy or curly coat' if n == 1 else 'Likely curly coat'
        else: txt = 'Likely wavy coat (furnishings add curl)' if furnished else 'Likely straight coat'
        return g, txt, conf
    if tid == 'shedding':
        c = calls[0]; n = c['n_alt']; furnished = ctx.get('furnished')
        g = 'C' * (2 - n) + 'T' * n
        if furnished: txt = 'Likely light shedding (furnished coats shed little regardless)'
        else: txt = 'Likely light shedding' if n == 2 else 'Likely heavy or seasonal shedding'
        return g, txt, conf
    if tid == 'hairless_xolo':
        return rule_dominant(calls, 'Dup', 'N', 'Likely hairless (one copy)', 'Very unlikely to be hairless')
    if tid == 'hairless_terrier':
        return rule_recessive(calls, 'D', 'N', 'Likely hairless', 'Normal coat, carries the hairless variant', 'Very unlikely to be hairless')
    if tid == 'albinism':
        return rule_recessive(calls, 'D', 'N', 'Likely albino (oculocutaneous albinism)', 'Not albino, carries the variant', 'Likely not albino')
    if tid == 'saddle_tan':
        # I = the 16-bp duplication (Embark's letter); II favours black-and-tan points, N/I or N/N saddle tan
        return rule_additive3(calls, 'I', 'N', 'Black-and-tan points more likely than a saddle (only shows in tan-pointed dogs)',
                              'Saddle tan pattern possible (only shows in tan-pointed dogs)', 'Saddle tan pattern possible (only shows in tan-pointed dogs)')
    if tid == 'red_intensity':
        # additive over the intensity loci: count intensifying (red) alleles, which may be REF or ALT per locus
        n = 0
        for l, c in zip(loci, calls):
            if c['n_alt'] is None: continue
            n += c['n_alt'] if l.get('red_allele') == l['alt'] else 2 - c['n_alt']
        m = 2 * len(known)
        frac = n / m
        txt = 'Intense red pigment' if frac >= 0.67 else ('Intermediate red pigment (tan to red)' if frac >= 0.34 else 'Dilute red pigment (cream to yellow)')
        return f'{n}/{m} intensifying alleles', txt, conf
    if tid == 'muzzle':
        return rule_additive3(calls, 'A', 'C', 'Likely short muzzle', 'Likely medium or long muzzle', 'Likely medium or long muzzle')
    if tid == 'bobtail':
        return rule_dominant(calls, 'G', 'C', 'Likely natural bobtail (one copy)', 'Likely normal-length tail')
    if tid == 'hind_dewclaws':
        n = min(2, sum(c['n_alt'] for c in known))
        return 'T' * n + 'C' * (2 - n), ('About a 50% chance of hind dewclaws' if n >= 1 else 'Unlikely to have hind dewclaws'), conf
    if tid == 'muscling':
        c = calls[0]; n = c['n_alt']; sex = ctx.get('sex')
        g = 'T' * n + 'C' * (2 - n)
        if n == 2 or (n == 1 and sex == 'male'): txt = 'Likely heavy back muscling (in bulky large breeds)'
        elif n == 1: txt = 'Carries one copy of the muscling variant (heavy muscling needs two in females)'
        else: txt = 'Likely normal muscling'
        return g, txt, conf
    if tid == 'eye_color':
        return rule_dominant(calls, 'Dup', 'N', 'More likely to have blue eyes', 'Less likely to have blue eyes')
    if tid.startswith('size_'):
        return rule_additive3(calls, trait['alt_label'], trait['ref_label'], trait['pheno_alt2'], trait['pheno_het'], trait['pheno_ref2'])
    if tid == 'altitude':
        return rule_dominant(calls, 'A', 'G', 'Enhanced tolerance of low oxygen (high altitude)', 'Normal altitude tolerance')
    if tid == 'appetite':
        return rule_dominant(calls, 'D', 'N', 'Higher food motivation (POMC deletion) — prone to overeating and weight gain', 'Normal food motivation')
    return None


def main():
    args = sys.argv[1:]; cat_path, out_path = args[:2]
    def opt(n): return args[args.index(n) + 1] if n in args else None
    cat = json.load(open(cat_path))
    bcf = pysam.VariantFile(opt('--bcf')) if opt('--bcf') else None
    bam = pysam.AlignmentFile(opt('--bam'), 'rb') if opt('--bam') else None
    ctx = {'sex': opt('--sex')}
    traits_out = []
    # furnishings first: it modifies curl and shedding
    order = sorted(cat['traits'], key=lambda t: 0 if t['id'] == 'furnishings' else 1)
    for tr in order:
        calls = [dict(locus=l['id'], **call_locus(l, bcf, bam)) for l in tr['loci']]
        res = interpret(tr, calls, ctx)
        entry = {'id': tr['id'], 'name': tr['name'], 'group': tr['group'], 'gene': tr['gene'],
                 'description': tr.get('description', ''), 'loci': calls}
        if res:
            g, txt, conf = res
            entry.update({'genotype': g, 'result': txt, 'confidence': conf, 'status': 'called'})
            if tr['id'] == 'furnishings': ctx['furnished'] = g.count('F') >= 1
        else:
            entry.update({'genotype': None, 'result': 'Not resolvable at this sequencing depth', 'confidence': 'none', 'status': 'no_call'})
        traits_out.append(entry)
    out = {'traits': traits_out, 'method': cat.get('method', ''), 'catalog_version': cat.get('version', '')}
    json.dump(out, open(out_path, 'w'), indent=2)
    print('traits: ' + '; '.join(f"{t['name']}: {t['genotype']} -> {t['result']} [{t['confidence']}]" for t in traits_out))


if __name__ == '__main__':
    main()
