#!/usr/bin/env python3
"""
AMY2B (pancreatic amylase) copy number from read depth -> amy2b_result.json.

  amy2b.py <markdup.bam> <out.json> [--fasta canFam4_plusY.fa]

The amylase gene sits on a ~15 kb unit (canFam4 chr6:47,375,000-47,390,500)
that dogs carry in 2 to ~30 copies per diploid genome (wolves ~2; Arendt 2014,
Freedman 2014). canFam4 holds one uniquely mappable copy, so every extra copy
piles its reads onto it and diploid copies = 2 x (depth over the unit / depth
of a two-copy region). The unit is AT-rich (GC 33%), so the two-copy depth is
taken from GC-matched control windows elsewhere on chr6 rather than the
flanks: GC dropout differs by chemistry, and a 1.9 Mb segmental duplication
around the gene (Arendt 2014) would inflate local flanks.

Uncertainty: read starts over the unit are Poisson with the overdispersion
seen across the control windows of the same kit; the control mean has its own
standard error. The interval is 95% on the log scale.
"""
import json, math, sys, argparse, zlib
import numpy as np
import pysam

CHROM = 'chr6'
UNIT = (47_375_500, 47_390_000)            # inner unit, 0-based half-open; edges trimmed
CTRL_SPAN = [(30_000_000, 46_500_000), (48_500_000, 65_000_000)]
WIN = UNIT[1] - UNIT[0]
GC_TOL = 0.02
MIN_MQ = 20
REF_COPIES = 2                             # control windows are two-copy (diploid)


def gc(fa, a, b):
    s = fa.fetch(CHROM, a, b).upper()
    n = len(s) - s.count('N')
    return (s.count('G') + s.count('C')) / n if n > 0.99 * len(s) else None


def control_windows(fa):
    target = gc(fa, *UNIT)
    out = []
    for a0, b0 in CTRL_SPAN:
        for a in range(a0, b0 - WIN, WIN):
            g = gc(fa, a, a + WIN)
            if g is not None and abs(g - target) <= GC_TOL:
                out.append((a, a + WIN))
    return target, out


def starts(bam, a, b, keep=1.0):
    cut = int(keep * 2 ** 32)
    return sum(1 for r in bam.fetch(CHROM, a, b)
               if (keep >= 1 or zlib.crc32(r.query_name.encode()) < cut)
               and r.mapping_quality >= MIN_MQ and not (r.is_secondary or r.is_supplementary
                                                       or r.is_duplicate or r.is_qcfail or r.is_unmapped)
               and a <= r.reference_start < b)


def estimate(bam_path, fasta, keep=1.0):
    """keep < 1 thins read pairs deterministically (validation of low-depth behaviour)."""
    fa = pysam.FastaFile(fasta)
    bam = pysam.AlignmentFile(bam_path, 'rb')
    unit_gc, wins = control_windows(fa)
    c = np.array([starts(bam, a, b, keep) for a, b in wins], float)
    # Drop control windows that are themselves CNVs or dropouts (outside 0.5-1.5x the median).
    med = np.median(c)
    ok = (c > 0.5 * med) & (c < 1.5 * med) if med > 0 else np.ones_like(c, bool)
    c = c[ok]
    n_unit = starts(bam, *UNIT, keep)
    half = (UNIT[0] + UNIT[1]) // 2
    halves = (starts(bam, UNIT[0], half, keep), starts(bam, half, UNIT[1], keep))
    mu = float(c.mean()) if len(c) else 0.0
    res = {
        'gene': 'AMY2B', 'locus': f'{CHROM}:{UNIT[0] + 1}-{UNIT[1]}', 'reference': 'canFam4',
        'unit_gc': round(unit_gc, 3), 'control_windows': int(len(c)),
        'unit_reads': n_unit, 'control_mean_reads': round(mu, 1),
        'unit_halves_reads': list(halves),
    }
    if mu < 5 or n_unit == 0:
        res.update(status='insufficient_depth', copies=None, ci_low=None, ci_high=None)
        return res
    phi = max(1.0, float(c.var(ddof=1)) / mu) if len(c) > 2 else 1.0   # overdispersion vs Poisson
    ratio = n_unit / mu
    copies = REF_COPIES * ratio
    se = math.sqrt(phi / n_unit + float(c.var(ddof=1)) / len(c) / mu ** 2)
    lo, hi = copies * math.exp(-1.96 * se), copies * math.exp(1.96 * se)
    res.update(status='called', copies=round(copies, 2), ci_low=round(lo, 2), ci_high=round(hi, 2),
               overdispersion=round(phi, 2), mean_depth_reads_per_window=round(mu, 1))
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('bam'); ap.add_argument('out')
    ap.add_argument('--fasta', default='canFam4_plusY.fa')
    ap.add_argument('--subsample', type=float, default=1.0)
    a = ap.parse_args()
    r = estimate(a.bam, a.fasta, a.subsample)
    json.dump(r, open(a.out, 'w'), indent=1)
    print(json.dumps(r), file=sys.stderr)
