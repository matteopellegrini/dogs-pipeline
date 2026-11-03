#!/usr/bin/env python3
"""Panel v3, step 3a: read evidence at the catalogue's structural variants, per kit (counts, not calls).

For each SV with known breakpoints (known_variants/panel_v3/panel_v3_sites.tsv, kind == 'sv'):
  bp_alt / bp_ref   fragments with breakpoint evidence (a D/I operation of the expected size, or a soft-clip at
                    the breakpoint whose clipped bases read as the alt haplotype) vs fragments spanning every
                    breakpoint cleanly. A clip counts only if its bases match the reference at the event's other
                    breakpoint better than they match the read's own continuation (either strand for inversions);
                    for insertions of unknown sequence, only if they are foreign to the locus (+-300 bp, both
                    strands) and not low-complexity. Bare clips, as the live caller gt_sv_from_bam counts them,
                    turned up in 2-3 fragments per dog at several loci in every breed. Clips need a >= 30 bp anchor;
                    a fragment with any alt evidence is alt.
  disc              read pairs whose mates straddle the event with the orientation it predicts (DEL: FR pair
                    spanning start..end; DUP: RF pair inside; INV: same-strand pair across one breakpoint); for
                    insertions copied from a known gene (SOURCE), pairs pointing at the insertion with the mate in
                    the source gene. Pairs straddling a >= 1 kb insertion with a normal insert add to bp_ref.
  depth_in / depth_flank
                    MAPQ>=20 read starts inside the event and in equal-sized flanking windows (DEL/DUP >= 2 kb),
                    so a deletion shows as depth_in/depth_flank ~0.5 (one copy) or ~0, a duplication as ~1.5.
The likelihood model is in panel_v3_breed_prior.py (sv mode); this file only counts.

  python3 sv_evidence.py sites  > sv_sites.tsv                 # the SV list with types and breakpoints
  python3 sv_evidence.py count  sv_sites.tsv KIT BAM [FASTA] > KIT.sv.tsv   # FASTA defaults to canFam4_plusY.fa
"""
import csv, os, re, sys

MIN_MQ, WINDOW, MIN_CLIP, MIN_ANCHOR = 20, 15, 8, 30
DEPTH_MIN_LEN, FLANK_MAX = 2000, 50000
# insertions copied from a known gene (retrogenes): read pairs with one read beside the insertion and the mate in the
# source gene are alt evidence. Parental FGF4 in canFam4: Cosmo's 18 retrogene mates land at chr18:48,869,203-48,873,308.
SOURCE = {'cat093': ('chr18', 48868700, 48873800), 'cat469': ('chr18', 48868700, 48873800)}
PAIR_MAX = 1000          # largest insert of a normal pair; also how far from the insertion a pair's anchor read may sit
SPAN_MIN_INS = 1000      # pairs straddling an insertion this long can only come from the reference allele

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
    w.writerow(['site_id', 'gene', 'chrom', 'start', 'end', 'type', 'svlen', 'ref', 'alt', 'source'])
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
        if t in ('INS', 'DUP') and svlen <= 1:
            svlen = 0   # "gross insertion" / duplication of unrecorded size: breakpoint clips only, no size match
        src = '%s:%d-%d' % SOURCE[s['site_id']] if s['site_id'] in SOURCE else ''
        w.writerow([s['site_id'], s['gene'], s['chrom'], start, end, t, svlen, s['ref'][:60], s['alt'][:60], src])

def mismatches(a, b):
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))

def best_hit(seg, ref):
    """Fewest mismatches of seg anywhere in ref (ungapped)."""
    if len(ref) < len(seg):
        return len(seg)
    return min(mismatches(seg, ref[i:i + len(seg)]) for i in range(len(ref) - len(seg) + 1))

RC = str.maketrans('ACGTN', 'TGCAN')
def rc(x):
    return x.translate(RC)[::-1]

def low_complexity(seg):
    """Homopolymer or short tandem repeat (period <= 6) over >= 70% of the tail: repeat-length variation, not a
    junction (SPAG1: GGC-repeat expansions in a Dachshund passed every other test)."""
    return any(sum(seg[i] == seg[i + k] for i in range(len(seg) - k)) >= 0.7 * (len(seg) - k)
               for k in range(1, 7) if len(seg) > k + 4)

def junction_ok(fa, chrom, seg, own_next, partner, t, local_win, local_lo):
    """A clipped tail is breakpoint evidence only if it reads as the alt haplotype.
    seg: clipped bases next to the junction; own_next: reference continuing past the read's own alignment
    (what a clip from sequencing error or misalignment would match); partner: (lo, hi) window around the
    event's other breakpoint, or None for insertions of unknown sequence; local_win: reference +-300 bp around
    this breakpoint, starting at local_lo. A tail that also matches the locus outside the partner window is a
    tandem repeat shifting, not a junction (IDUA: GC-rich repeat at both ends of the deletion)."""
    tol = max(1, len(seg) // 8)
    local = mismatches(seg, own_next)
    if low_complexity(seg):
        return False
    if partner is None:  # insertion: the clip must be foreign to the locus
        return local > len(seg) // 3 and min(best_hit(seg, local_win), best_hit(rc(seg), local_win)) > tol
    win = fa.fetch(chrom, max(0, partner[0]), partner[1]).upper()
    hit = best_hit(seg, win)
    if t == 'INV':
        hit = min(hit, best_hit(rc(seg), win))
    a, b = max(0, partner[0] - local_lo), max(0, partner[1] - local_lo)
    elsewhere = local_win[:a] + 'N' * max(0, min(b, len(local_win)) - a) + local_win[b:]
    return hit <= tol and hit + 2 <= local and min(best_hit(seg, elsewhere), best_hit(rc(seg), elsewhere)) > tol

def breakpoint_counts(bam, fa, chrom, bps, t, svlen, start, end):
    """Fragments with junction evidence vs fragments spanning all breakpoints cleanly."""
    small = 2 <= svlen <= 60 and t in ('DEL', 'INS')   # svlen 0 = size unknown
    # the catalogue size, +-10%: the live caller's 0.5-2x window let a 24-bp GGC-repeat expansion in a Dachshund
    # pass for SPAG1's 29-bp duplication
    size_tol = max(2, svlen // 10)
    alt, span_at = set(), {}
    for bp in bps:
        lo, hi = bp - WINDOW, bp + WINDOW
        span_at[bp] = set()
        # the other side of the event: for a two-breakpoint event the other breakpoint; for a short deletion
        # written at one position, its end; insertions of unknown sequence have none
        if t == 'INS':
            partner = None
        elif len(bps) == 2:
            o = bps[1] if bp == bps[0] else bps[0]
            partner = (o - 80, o + 80)
        else:
            partner = (end - 20, end + 20) if t == 'DEL' else None
        local_lo = max(0, bp - 300)
        local_win = fa.fetch(chrom, local_lo, bp + 300).upper()
        for r in bam.fetch(chrom, max(0, bp - 300), bp + 300):
            if r.is_unmapped or r.mapping_quality < MIN_MQ or r.is_secondary or r.is_supplementary:
                continue
            cig = r.cigartuples or []
            rs, re_ = r.reference_start + 1, r.reference_end
            aligned = sum(ln for op, ln in cig if op in (0, 7, 8))
            hit = False
            rp = rs
            for op, ln in cig:  # D/I of the expected size near the breakpoint (events that fit inside a read)
                if op in (0, 7, 8, 3):
                    rp += ln
                elif op == 2:
                    if small and t == 'DEL' and rp - 1 <= hi and rp + ln >= lo and abs(ln - svlen) <= size_tol:
                        hit = True
                    rp += ln
                elif op == 1 and small and t == 'INS' and lo <= rp <= hi and abs(ln - svlen) <= size_tol:
                    hit = True
            q = r.query_sequence
            if not hit and q and aligned >= MIN_ANCHOR:
                if cig and cig[0][0] == 4 and cig[0][1] >= MIN_CLIP and lo <= rs <= hi:
                    seg = q[:cig[0][1]][-25:].upper()
                    own = fa.fetch(chrom, max(0, rs - 1 - len(seg)), rs - 1).upper()
                    hit = junction_ok(fa, chrom, seg, own, partner, t, local_win, local_lo)
                if not hit and cig and cig[-1][0] == 4 and cig[-1][1] >= MIN_CLIP and lo <= re_ <= hi:
                    seg = q[len(q) - cig[-1][1]:][:25].upper()
                    own = fa.fetch(chrom, re_, re_ + len(seg)).upper()
                    hit = junction_ok(fa, chrom, seg, own, partner, t, local_win, local_lo)
            if hit:
                alt.add(r.query_name)
            elif rs <= lo - 5 and re_ >= hi + 5:
                span_at[bp].add(r.query_name)
    # a reference fragment must span a breakpoint without alt evidence anywhere; for two-breakpoint events either
    # breakpoint counts (each spans independently on a reference haplotype)
    ref = set().union(*span_at.values()) - alt
    return alt, ref

def insertion_pairs(bam, chrom, bp, svlen, source):
    """Read-pair evidence at an insertion. alt: pairs anchored within PAIR_MAX of the insertion, pointing at it,
    whose mate lies in the source gene. ref: normal-insert pairs straddling the insertion point, for insertions of
    >= SPAN_MIN_INS (a fragment cannot reach across one on the alt allele)."""
    alt, ref = set(), set()
    for r in bam.fetch(chrom, max(0, bp - PAIR_MAX), bp + PAIR_MAX):
        if (r.is_unmapped or r.mate_is_unmapped or r.mapping_quality < MIN_MQ or r.is_secondary
                or r.is_supplementary):
            continue
        left = r.reference_end is not None and r.reference_end <= bp + 10
        toward = (left and not r.is_reverse) or (not left and r.is_reverse and r.reference_start >= bp - 10)
        if source and toward and r.next_reference_name == source[0] and source[1] <= r.next_reference_start <= source[2]:
            alt.add(r.query_name)
        elif (svlen >= SPAN_MIN_INS and r.next_reference_name == chrom and not r.is_reverse and r.mate_is_reverse
              and r.reference_end is not None and r.reference_end < bp - 10 and r.next_reference_start > bp + 10
              and 0 < r.template_length <= PAIR_MAX):
            ref.add(r.query_name)
    return alt, ref - alt

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

def count(sites_tsv, kit, bam_path, fasta=None):
    import pysam
    bam = pysam.AlignmentFile(bam_path, 'rb')
    fa = pysam.FastaFile(fasta or os.path.join(os.path.dirname(__file__), '..', '..', 'canFam4_plusY.fa'))
    w = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    for s in csv.DictReader(open(sites_tsv), delimiter='\t'):
        c, start, end, t, svlen = s['chrom'], int(s['start']), int(s['end']), s['type'], int(s['svlen'])
        bps = [start] if t == 'INS' or end - start <= 60 else [start, end]
        t_bp = 'INS' if t == 'DUP' and 2 <= svlen <= 60 else t   # a short tandem dup reads as an insertion
        try:
            alt, ref = breakpoint_counts(bam, fa, c, bps, t_bp, svlen, start, end)
            if t == 'INS':   # pair evidence goes in 'disc'; bp_alt keeps the junction-only fragments (no double count)
                src = s.get('source', '')
                src = (src.split(':')[0], *map(int, src.split(':')[1].split('-'))) if src else None
                palt, pref = insertion_pairs(bam, c, start, svlen, src)
                alt -= palt
                ref = (ref | pref) - alt - palt
                a, r, d = len(alt), len(ref), len(palt)
            else:
                a, r, d = len(alt), len(ref), discordant(bam, c, start, end, t)
            di, df = depth(bam, c, start, end, t)
        except (ValueError, KeyError):   # contig absent from this BAM
            continue
        w.writerow([kit, s['site_id'], a, r, d, di, df, end - start])

if __name__ == '__main__':
    if sys.argv[1] == 'sites':
        sites()
    else:
        count(*sys.argv[2:6])
