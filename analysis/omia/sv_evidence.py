#!/usr/bin/env python3
"""Panel v3, step 3a: read evidence at the catalogue's structural variants, per kit (counts, not calls).

For each SV with known breakpoints (known_variants/panel_v3/panel_v3_sites.tsv, kind == 'sv'):
  bp_alt / bp_ref   fragments with breakpoint evidence (soft-clip at the breakpoint, or a D/I operation of the
                    expected size) vs fragments spanning every breakpoint cleanly. Same fragment-level rules as
                    the live caller gt_sv_from_bam in run_dog_pipeline.sh (clips need a >= 30 bp anchor; a
                    fragment with any alt evidence is alt).
  disc              read pairs whose mates straddle the event with the orientation it predicts (DEL: FR pair
                    spanning start..end; DUP: RF pair inside; INV: same-strand pair across one breakpoint).
  depth_in / depth_flank
                    MAPQ>=20 read starts inside the event and in equal-sized flanking windows (DEL/DUP >= 2 kb),
                    so a deletion shows as depth_in/depth_flank ~0.5 (one copy) or ~0, a duplication as ~1.5.
The likelihood model is in panel_v3_breed_prior.py (sv mode); this file only counts.

  python3 sv_evidence.py sites  > sv_sites.tsv                 # the SV list with types and breakpoints
  python3 sv_evidence.py count  sv_sites.tsv KIT BAM > KIT.sv.tsv
"""
import csv, os, re, sys

MIN_MQ, WINDOW, MIN_CLIP, MIN_ANCHOR = 20, 15, 8, 30
DEPTH_MIN_LEN, FLANK_MAX = 2000, 50000

def sv_type(alt, status):
    a = (alt or '').lower()
    if 'del' in a and 'ins' in a: return 'DEL'          # delins, gross: the deleted segment dominates
    for k, t in (('<del>', 'DEL'), ('deletion', 'DEL'), ('<ins>', 'INS'), ('insertion', 'INS'), ('<dup>', 'DUP'),
                 ('duplication', 'DUP'), ('<inv>', 'INV'), ('inversion', 'INV'), ('<cpx>', 'CPX')):
        if k in a: return t
    if re.fullmatch(r'[acgtn]+', a): return 'SEQ'       # explicit long alleles (written out in the VCF sense)
    return None                                         # haplotype, repeat variation: not breakpoint-typable

def sites():
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    w = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    w.writerow(['site_id', 'gene', 'chrom', 'start', 'end', 'type', 'svlen', 'ref', 'alt'])
    for s in csv.DictReader(open(os.path.join(root, 'known_variants/panel_v3/panel_v3_sites.tsv')), delimiter='\t'):
        if s['kind'] != 'sv' or s['status'] == 'gene_span_only' or not s['pos'].isdigit():
            continue
        t = sv_type(s['alt'], s['status'])
        if t is None:
            continue
        start, end = int(s['pos']), int(s['end']) if s['end'].isdigit() else int(s['pos'])
        m = re.search(r'SVLEN=(-?\d+)', s['status'])
        svlen = abs(int(m.group(1))) if m else (abs(len(s['alt']) - len(s['ref'])) if t == 'SEQ' else end - start)
        if t == 'SEQ':
            t = 'DEL' if len(s['ref']) > len(s['alt']) else 'INS'
            end = start + len(s['ref']) - 1
        w.writerow([s['site_id'], s['gene'], s['chrom'], start, end, t, svlen, s['ref'][:60], s['alt'][:60]])

def breakpoint_counts(bam, chrom, bps, t, svlen):
    """Fragments with alt evidence vs fragments spanning all breakpoints cleanly."""
    small = svlen <= 60 and t in ('DEL', 'INS')
    alt, span_at = set(), {}
    for bp in bps:
        lo, hi = bp - WINDOW, bp + WINDOW
        span_at[bp] = set()
        for r in bam.fetch(chrom, max(0, bp - 300), bp + 300):
            if r.is_unmapped or r.mapping_quality < MIN_MQ or r.is_secondary or r.is_supplementary:
                continue
            cig = r.cigartuples or []
            start, end = r.reference_start + 1, r.reference_end
            aligned = sum(ln for op, ln in cig if op in (0, 7, 8))
            hit = False
            rp = start
            for op, ln in cig:  # D/I of the expected size near the breakpoint (events that fit inside a read)
                if op in (0, 7, 8, 3):
                    rp += ln
                elif op == 2:
                    if small and t == 'DEL' and rp - 1 <= hi and rp + ln >= lo and svlen // 2 <= ln <= svlen * 2:
                        hit = True
                    rp += ln
                elif op == 1 and small and t == 'INS' and lo <= rp <= hi and max(1, svlen // 2) <= ln <= svlen * 2:
                    hit = True
            clip_l = cig and cig[0][0] in (4, 5) and cig[0][1] >= MIN_CLIP and lo <= start <= hi
            clip_r = cig and cig[-1][0] in (4, 5) and cig[-1][1] >= MIN_CLIP and lo <= end <= hi
            if (clip_l or clip_r) and aligned >= MIN_ANCHOR:
                hit = True
            if hit:
                alt.add(r.query_name)
            elif start <= lo - 5 and end >= hi + 5:
                span_at[bp].add(r.query_name)
    # a reference fragment must span a breakpoint without alt evidence anywhere; for two-breakpoint events either
    # breakpoint counts (each spans independently on a reference haplotype)
    ref = set().union(*span_at.values()) - alt
    return len(alt), len(ref)

def discordant(bam, chrom, start, end, t, insert_max=1000):
    if t not in ('DEL', 'DUP', 'INV') or end - start < 300:
        return 0
    names = set()
    for r in bam.fetch(chrom, max(0, start - insert_max), start + insert_max):
        if (r.is_unmapped or r.mate_is_unmapped or r.mapping_quality < MIN_MQ or r.is_secondary
                or r.is_supplementary or r.next_reference_name != chrom):
            continue
        mpos = r.next_reference_start
        if t == 'DEL':   # forward read left of start, reverse mate right of end
            ok = (not r.is_reverse and r.mate_is_reverse and r.reference_start < start and end <= mpos <= end + insert_max)
        elif t == 'DUP':  # tandem junction: reverse read near start, forward mate near end
            ok = (r.is_reverse and not r.mate_is_reverse and r.reference_start >= start - 50
                  and end - insert_max <= mpos <= end + 50)
        else:            # inversion: both mates on the same strand, one at each breakpoint
            ok = (r.is_reverse == r.mate_is_reverse and abs(mpos - end) <= insert_max)
        if ok:
            names.add(r.query_name)
    return len(names)

def depth(bam, chrom, start, end, t):
    if t not in ('DEL', 'DUP') or end - start < DEPTH_MIN_LEN:
        return '', ''
    def starts(a, b):
        return sum(1 for r in bam.fetch(chrom, max(0, a), b) if not r.is_unmapped and r.mapping_quality >= MIN_MQ
                   and not r.is_secondary and not r.is_supplementary and a <= r.reference_start < b)
    flank = min(end - start, FLANK_MAX)
    return starts(start, end), starts(start - flank, start) + starts(end, end + flank)

def count(sites_tsv, kit, bam_path):
    import pysam
    bam = pysam.AlignmentFile(bam_path, 'rb')
    w = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    for s in csv.DictReader(open(sites_tsv), delimiter='\t'):
        c, start, end, t, svlen = s['chrom'], int(s['start']), int(s['end']), s['type'], int(s['svlen'])
        bps = [start] if t == 'INS' or end - start <= 60 else [start, end]
        try:
            a, r = breakpoint_counts(bam, c, bps, t, svlen)
            d = discordant(bam, c, start, end, t)
            di, df = depth(bam, c, start, end, t)
        except (ValueError, KeyError):   # contig absent from this BAM
            continue
        w.writerow([kit, s['site_id'], a, r, d, di, df, end - start])

if __name__ == '__main__':
    if sys.argv[1] == 'sites':
        sites()
    else:
        count(*sys.argv[2:5])
