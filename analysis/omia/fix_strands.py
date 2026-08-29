#!/usr/bin/env python3
"""One-time strand repair for reference_json/omia_variants.json.

A block of catalogue SNVs (VWF including vWD Type 1, DNM1/EIC, FVIII, and
others) carried ref/alt on the OPPOSITE strand from canFam4: stage 8's
assembly-mismatch guard correctly refused to genotype them (the guard exists
because where alt==assembly base, every healthy dog reads as affected).

The fix is principled, not a blind flip. An entry is complemented only when
ALL of:
  1. its stated ref does NOT match the canFam4 base,
  2. the COMPLEMENT of its ref DOES match the canFam4 base,
  3. the SNV is not strand-ambiguous (A/T and C/G pairs are skipped —
     complementing those cannot be verified),
Entries that fail these stay as they are, still guarded.

Verification anchors (Dog10K panel carries the complemented allele pair at
the exact position): DNM1 chr9:55370803 C>A (AC=12/3858), vWD1
chr27:7140281 C>T (AC=35/3858).

Usage: python3 analysis/omia/fix_strands.py   (rewrites the JSON in place,
prints every change; needs canFam4.fa + .fai in the repo root)
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FASTA = os.path.join(ROOT, 'canFam4.fa')
DB = os.path.join(ROOT, 'reference_json', 'omia_variants.json')
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def fetch_base(chrom, pos, fai, fh):
    length, offset, linebases, linewidth = fai[chrom]
    if pos < 1 or pos > length: return None
    fh.seek(offset + (pos - 1) // linebases * linewidth + (pos - 1) % linebases)
    return fh.read(1).upper()

def main():
    fai = {}
    for line in open(FASTA + '.fai'):
        p = line.split('\t')
        fai[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
    fh = open(FASTA)
    db = json.load(open(DB))
    vs = db.get('variants', db if isinstance(db, list) else [])
    fixed, ambiguous, still_bad = 0, 0, 0
    for v in vs:
        ref, alt = v.get('ref'), v.get('alt')
        chrom, pos = v.get('chrom'), v.get('pos')
        if not (ref and alt and chrom and pos): continue
        if len(ref) != 1 or len(alt) != 1: continue
        if chrom not in fai: continue
        base = fetch_base(chrom, int(pos), fai, fh)
        if base is None or ref.upper() == base: continue
        if COMP.get(ref.upper()) == base:
            if COMP.get(ref.upper()) == alt.upper():
                ambiguous += 1     # A/T or C/G — unverifiable, leave guarded
                continue
            v['ref'], v['alt'] = COMP[ref.upper()], COMP[alt.upper()]
            v['strand_fixed'] = True
            fixed += 1
            print('FLIPPED %-8s %s:%s %s>%s -> %s>%s  (%s)' % (
                v.get('gene'), chrom, pos, ref, alt, v['ref'], v['alt'],
                str(v.get('trait'))[:40]))
        else:
            still_bad += 1
            print('STILL MISMATCHED %-8s %s:%s ref=%s canFam4=%s (%s)' % (
                v.get('gene'), chrom, pos, ref, base, str(v.get('trait'))[:40]))
    json.dump(db, open(DB, 'w'), indent=1)
    print('\nflipped: %d · ambiguous(skipped): %d · still mismatched: %d' % (fixed, ambiguous, still_bad))

if __name__ == '__main__':
    sys.exit(main())
