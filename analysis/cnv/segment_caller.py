#!/usr/bin/env python3
"""Unified CNV/aneuploidy caller: segmentation of 50kb panel-normalized coverage.

One framework across scales, replacing the dual 1Mb-karyotype / adaptive-CNV
split: per-window log2 ratios against the 1,257-dog panel-of-normals are
segmented per chromosome by recursive binary splitting (dependency-free CBS
analogue), and segments classify by size — whole-chromosome segments are
ANEUPLOIDY findings, sub-chromosomal ones are CNVs.

Noise model: per-sample residual SD estimated robustly from that sample's own
autosomal log-ratios (low-pass depth varies 0.1-6x across the cohort, so a
fixed threshold would over-call shallow dogs). A split is accepted when the
two sides differ by > SPLIT_Z * se; a segment is REPORTED when its mean
log-ratio is significant (|mean|/se > CALL_Z) AND biologically meaningful
(|mean lr| >= 0.32, i.e. >=25% coverage change).

Usage (library): call_sample(v, sex, panel) -> dict
Usage (CLI, dev): segment_caller.py <matrix.npz> <panel.npz> [sample ...]
"""
import sys, json
import numpy as np

MIN_WINDOWS   = 2      # smallest reportable event: 100kb
SPLIT_Z       = 6.0    # split acceptance
CALL_Z        = 8.0    # segment report threshold (mean-level significance)
MIN_ABS_LR    = 0.32   # >=25% coverage shift — below this it's noise/waves
ANEU_FRACTION = 0.85   # segment covering >=85% of a chromosome = aneuploidy


def _split_scores(x):
    """For each interior split point, |mean(left)-mean(right)| in units of se."""
    n = len(x)
    cs = np.cumsum(x)
    i = np.arange(1, n)
    ml = cs[:-1] / i
    mr = (cs[-1] - cs[:-1]) / (n - i)
    se = np.sqrt(1.0 / i + 1.0 / (n - i))
    return np.abs(ml - mr) / se


def _segment(x, sigma, lo, out):
    """Recursive binary segmentation of x (log-ratios); indices offset by lo."""
    n = len(x)
    if n < 2 * MIN_WINDOWS:
        out.append((lo, lo + n)); return
    sc = _split_scores(x) / max(sigma, 1e-6)
    # forbid splits that create fragments shorter than MIN_WINDOWS
    sc[: MIN_WINDOWS - 1] = 0
    if MIN_WINDOWS > 1:
        sc[-(MIN_WINDOWS - 1):] = 0
    k = int(np.argmax(sc))
    if sc[k] < SPLIT_Z:
        out.append((lo, lo + n)); return
    _segment(x[: k + 1], sigma, lo, out)
    _segment(x[k + 1:], sigma, lo + k + 1, out)


def call_sample(v, sex, panel):
    """v: normalized coverage vector (sample/autosomal-median), one per window.

    Returns {'segments': [...], 'aneuploidies': [...], 'sigma': float}."""
    chrom = panel['chrom']
    bad = panel['bad']
    med = np.where(chrom == 'chrX',
                   panel['medF'] if sex == 'F' else panel['medM'],
                   panel['med']).astype(np.float64)
    good = (~bad) & (med > 0.25) & np.isfinite(v)
    lr = np.full(len(v), np.nan)
    lr[good] = np.log2(np.maximum(v[good], 1e-3) / med[good])
    auto_good = good & (chrom != 'chrX')
    resid = lr[auto_good]
    sigma = float(1.4826 * np.median(np.abs(resid - np.median(resid))))

    # Sample quality gate: wavy/degraded libraries (residual sigma >= 0.4)
    # produce genome-wide pseudo-events in both directions — report nothing
    # rather than fiction, mirroring the old coverage-confidence gate.
    confidence = 'high' if sigma < 0.2 else ('medium' if sigma < 0.4 else 'low')
    if confidence == 'low':
        return {'segments': [], 'aneuploidies': [], 'sigma': round(sigma, 4),
                'confidence': 'low',
                'note': 'Coverage too uneven for reliable copy-number analysis on this sample.'}

    segments, aneuploidies = [], []
    for c in list(dict.fromkeys(chrom)):     # keep genomic order
        cidx = np.where((chrom == c) & good)[0]
        if len(cidx) < MIN_WINDOWS:
            continue
        x = lr[cidx]
        n_chrom_windows = int((chrom == c).sum())

        # Aneuploidy is a CHROMOSOME-level question, decided before
        # segmentation: a uniformly shifted chromosome (e.g. trisomy X at
        # +0.55) can fragment under segmentation and never yield one >=85%
        # segment. Median shift + majority agreement is the robust test.
        cmed = float(np.median(x))
        frac_shifted = float(np.mean(np.sign(x) == np.sign(cmed))) if cmed != 0 else 0.0
        zc = abs(cmed) / (sigma / np.sqrt(len(x)))
        if abs(cmed) >= 0.18 and zc >= CALL_Z and frac_shifted >= 0.75 and len(x) >= 0.5 * n_chrom_windows:
            aneuploidies.append({
                'chrom': c, 'ratio': round(2 ** cmed, 3),
                'copy_estimate': round(2 * 2 ** cmed, 2) if c != 'chrX' else round((2 if sex == 'F' else 1) * 2 ** cmed, 2),
                'mean_log2': round(cmed, 3), 'z': round(zc, 1),
                'n_windows': len(x)})
            x = x - cmed   # segment the residual for additional sub-events

        # Per-chromosome noise floor: chrX (and occasionally others) carries
        # more panel-relative variance than the autosomal average.
        sigma_c = max(sigma, float(1.4826 * np.median(np.abs(x - np.median(x)))))
        segs = []
        _segment(x, sigma_c, 0, segs)
        for a, b in segs:
            m = float(np.mean(x[a:b]))
            nw = b - a
            z = abs(m) / (sigma_c / np.sqrt(nw))
            if nw / n_chrom_windows >= ANEU_FRACTION:
                continue   # chromosome-level already handled above
            if nw < MIN_WINDOWS or abs(m) < MIN_ABS_LR or z < CALL_Z:
                continue
            w0, w1 = int(cidx[a]), int(cidx[b - 1])
            # Population frequency from the cohort recurrence track: events
            # shared by many healthy dogs are normal copy-number polymorphisms
            # (olfactory clusters etc.), presented separately from rare events.
            freq = None
            if 'seg_freq' in panel:
                freq = float(np.median(panel['seg_freq'][w0:w1 + 1]))
            # windows are laid out per-chromosome in fixed 50kb steps, so the
            # genomic start is the window's rank within its chromosome * 50kb
            first = int(np.where(chrom == c)[0][0])
            segments.append({
                'chrom': c,
                'start': (w0 - first) * 50000,
                'end': (w1 - first) * 50000 + 50000,
                'window_start': w0, 'window_end': w1,
                'n_windows': nw,
                'mean_log2': round(m, 3),
                'fold': round(2 ** m, 2),
                'type': 'loss' if m < 0 else 'gain',
                'copy_estimate': round(2 * 2 ** m, 2),
                'z': round(z, 1),
                **({'panel_frequency': round(freq, 3), 'common_variant': freq >= 0.05}
                   if freq is not None else {})})
    # A profile throwing dozens of large 'events' is instability, not biology.
    if len(segments) > 25:
        return {'segments': [], 'aneuploidies': aneuploidies, 'sigma': round(sigma, 4),
                'confidence': 'low',
                'note': f'{len(segments)} candidate segments — unstable profile, CNV calls suppressed.'}
    return {'segments': segments, 'aneuploidies': aneuploidies, 'sigma': round(sigma, 4),
            'confidence': confidence}


def main():
    mat = np.load(sys.argv[1], allow_pickle=True)
    panel = dict(np.load(sys.argv[2], allow_pickle=True))
    M, names, sexes = mat['M'], list(mat['names']), mat['sex']
    targets = sys.argv[3:] or names
    for s in targets:
        i = names.index(s)
        r = call_sample(M[:, i].astype(np.float64), str(sexes[i]), panel)
        print(json.dumps({'sample': s, **r}))


if __name__ == '__main__':
    main()
