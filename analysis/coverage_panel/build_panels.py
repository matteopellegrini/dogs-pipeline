#!/usr/bin/env python3
"""
Build coverage panels-of-normals at both scales from many dogs.

    python3 build_panels.py OUT_DIR WORKDIR_GLOB [--cnv-bin 50000] [--min-n 20]

    e.g. python3 build_panels.py reference_panel 'work/DOGS-Gen-*/analysis'

Produces
    coverage_1mb_panel.json    2,373 fixed 1Mb windows
    coverage_cnv_panel.json    fixed CNV-scale bins (default 50kb)

Why two scales, and why they need different handling
----------------------------------------------------
The 1Mb windows come off a fixed grid derived from the .fai, so every dog has
the same 2,373 windows and a panel is a straight per-window summary.

The CNV windows do not. Stage 5 sizes them as 50000/mean_depth clamped to
15kb-200kb, so each dog is on its OWN grid — observed across 96 dogs: min
15,000 (the floor), median 17,953, max 51,102. There is no shared window to
take a mean over, which is why a CNV panel was never built.

The fix is to re-bin each dog onto one fixed grid. samtools bedcov gives summed
base depth per window, so a source window overlapping two target bins
contributes in proportion to the overlap. That assumes depth is uniform within
a source window — exactly the assumption the original binning already makes, so
it adds no new error. It is only valid DOWNWARD: the target bin must be at
least as wide as the widest source window, or we would be inventing resolution
that was never measured. 50kb clears every observed window bar one at 51.1kb.

Statistics: median and MAD*1.4826, not mean and SD
--------------------------------------------------
The cohort is unscreened, so some dogs genuinely carry large CNVs. One dog with
a 3x amplification would inflate that window's SD enough to hide the same event
in every other dog — which defeats the point. The median ignores a minority of
outliers and MAD*1.4826 estimates the SD of the underlying normal part.

chrX is computed per inferred sex, since pooling ~1x females with ~0.5x males
gives a bimodal distribution whose median describes no actual dog.
"""
import glob
import json
import os
import statistics
import sys
from collections import defaultdict

# Frequency counting must use exactly the same scoreability rules and thresholds
# the scorer applies, or the stored frequencies describe a different test than
# the one that produces events. Import rather than restate them.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_sample import MIN_SD, MAX_SD, REPORT_Z, REPORT_Z_BY_WINDOW

MAD_TO_SD = 1.4826
SEX_THRESHOLD = 0.75          # chrX:autosome ratio; ~1.0 female, ~0.5 male
COMMON_PCT = 5.0              # >= this share of the cohort = common polymorphism


def read_bedcov(path):
    """-> list of (chrom, start, end, summed_depth)"""
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 4:
                rows.append((f[0], int(f[1]), int(f[2]), float(f[3])))
    return rows


def rebin(rows, bin_bp):
    """Per-bp depth on a fixed grid, weighting each source window by overlap.

    Divides by the length ACTUALLY covered, not by bin_bp. A chromosome's final
    bin is partial, and dividing it by the full bin width would make it look
    proportionally depleted — which is exactly the bug that made every
    chromosome end flag in 96/96 dogs.
    """
    acc = defaultdict(float)
    span = defaultdict(float)
    for chrom, s, e, total in rows:
        width = e - s
        if width <= 0:
            continue
        depth_per_bp = total / width
        b0, b1 = s // bin_bp, (e - 1) // bin_bp
        for b in range(b0, b1 + 1):
            lo = max(s, b * bin_bp)
            hi = min(e, (b + 1) * bin_bp)
            if hi > lo:
                acc[(chrom, b * bin_bp)] += depth_per_bp * (hi - lo)
                span[(chrom, b * bin_bp)] += (hi - lo)
    return {k: v / span[k] for k, v in acc.items() if span[k] > 0}


def to_ratio(depths):
    """summed depth per window -> per-bp depth normalised by autosomal median"""
    per_bp = {}
    for (chrom, start), total in depths.items():
        per_bp[(chrom, start)] = total
    auto = [v for (c, _), v in per_bp.items() if c != 'chrX' and v > 0]
    if not auto:
        return None, None
    base = statistics.median(auto)
    ratios = {k: v / base for k, v in per_bp.items()}
    x = [v for (c, _), v in ratios.items() if c == 'chrX' and v > 0]
    sex = 'F' if (statistics.median(x) if x else 1.0) >= SEX_THRESHOLD else 'M'
    return ratios, sex


def summarise(samples, sexes, min_n, label, store_values=False):
    windows = set()
    for r in samples.values():
        windows.update(r)
    panel = defaultdict(list)
    dropped = defaultdict(int)
    for chrom, start in sorted(windows, key=lambda w: (w[0], w[1])):
        if chrom == 'chrX':
            groups = {s: [samples[n][(chrom, start)] for n in samples
                          if sexes[n] == s and (chrom, start) in samples[n]]
                      for s in ('F', 'M')}
        else:
            groups = {'all': [samples[n][(chrom, start)] for n in samples
                              if (chrom, start) in samples[n]]}
        entry = {'start': start}
        for key, vals in groups.items():
            if len(vals) < min_n:
                dropped[(chrom, key)] += 1
                continue
            med = statistics.median(vals)
            mad = statistics.median([abs(v - med) for v in vals]) if len(vals) > 1 else 0.0
            entry[key] = {'median': round(med, 4),
                          'mad_sd': round(mad * MAD_TO_SD, 4),
                          'n': len(vals)}
            if store_values:
                # The sorted per-dog values at this window, so a scorer can ask
                # the question a reader actually has: "how many reference dogs
                # are as far from normal here as mine?" Median and spread cannot
                # answer that — they describe the centre, and the answer lives in
                # the tail. Stored only for the 1Mb panel; at 50kb this would be
                # 47k windows x 92 dogs and the file would be unusable.
                entry[key]['vals'] = [round(v, 3) for v in sorted(vals)]
        if len(entry) > 1:
            panel[chrom].append(entry)
    total = sum(len(v) for v in panel.values())
    print(f"  {label}: {total} windows across {len(panel)} chromosomes")
    for (chrom, key), n in sorted(dropped.items()):
        grp = '' if key == 'all' else f' [{key}]'
        print(f"     DROPPED {chrom}{grp}: {n} windows with < {min_n} samples")
    spreads = sorted(w['all']['mad_sd'] for c in panel for w in panel[c] if 'all' in w)
    if spreads:
        print(f"     autosomal robust SD: median {spreads[len(spreads)//2]:.3f}  "
              f"90th {spreads[int(len(spreads)*0.9)]:.3f}  max {spreads[-1]:.3f}")
    return dict(panel)


def add_frequencies(panel, samples, sexes, z_cut):
    """Count, per bin, how many cohort dogs are flagged there, by direction.

    Why this belongs in the panel
    -----------------------------
    Recurrence over the 93 clean dogs found chr19:21.3-21.5Mb flagged as a gain
    in 42 of them, chr5:85.9Mb in 41, and a chr6:45.6-47.0Mb cluster in 19. No
    bin was flagged in more than half, which is the important part: a reference
    artifact would hit nearly EVERY dog, because every dog is mapped against the
    same reference. Hitting a fraction is what a genuine copy-number
    polymorphism looks like.

    So these events are real, and still not findings. A variant carried by 45%
    of the population is common, and reporting it as a discovery about one dog
    is a category error. Storing the frequency lets the report rank by rarity,
    exactly as allele frequency is used for small variants.

    Frequencies are computed against the panel's own cohort, so a bin's count
    includes every dog that contributed to its median. That is the intended
    denominator: we are asking "how often is this region flagged among normal
    dogs", not estimating a population rate from an independent sample.
    """
    idx = {}
    for chrom, entries in panel.items():
        for e in entries:
            for k in ('all', 'F', 'M'):
                if k in e:
                    idx[(chrom, e['start'], k)] = e[k]
    for name, ratios in samples.items():
        for (chrom, start), v in ratios.items():
            key = 'all' if chrom != 'chrX' else sexes[name]
            st = idx.get((chrom, start, key))
            if st is None or v <= 0:
                continue
            sd = st['mad_sd']
            if not (MIN_SD <= sd <= MAX_SD):
                continue
            z = (v - st['median']) / sd
            if z >= z_cut:
                st['ng'] = st.get('ng', 0) + 1
            elif z <= -z_cut:
                st['nl'] = st.get('nl', 0) + 1
    flagged = [(e[k].get('ng', 0) + e[k].get('nl', 0)) / e[k]['n']
               for c in panel for e in panel[c]
               for k in ('all', 'F', 'M') if k in e and e[k].get('n')]
    common = sum(1 for f in flagged if f >= COMMON_PCT / 100)
    print(f"     frequencies at |z| >= {z_cut}: {common} of {len(flagged)} bins "
          f"flagged in >= {COMMON_PCT:.0f}% of the cohort")
    return panel


def main():
    out_dir = sys.argv[1]
    pattern = sys.argv[2]
    args = sys.argv[3:]
    cnv_bin = int(args[args.index('--cnv-bin') + 1]) if '--cnv-bin' in args else 50000
    min_n = int(args[args.index('--min-n') + 1]) if '--min-n' in args else 20
    # A panel of normals should not contain dogs whose coverage is not normal.
    # Median/MAD limits their influence but does not remove it, and a sample
    # that cannot be scored reliably has no business defining what "reliable"
    # means. Identify candidates by event count and gain/loss balance, not by
    # depth alone — DOGS-Gen-29 is the lowest-coverage dog in the cohort at
    # 1.00x and produced only 10 events, while DOGS-Gen-54 at 2.40x produced 86,
    # of which 85 were gains. That asymmetry is a mapping artifact signature;
    # real copy-number variation is not one-directional.
    exclude = set()
    if '--exclude' in args:
        exclude = {x.strip() for x in args[args.index('--exclude') + 1].split(',') if x.strip()}

    # Comma-separated globs: the panel now spans the study cohort AND the
    # ProsperKits cohort (same Illumina low-pass platform; UCLA MGI dogs are
    # deliberately not included — cross-platform panels mis-score everything).
    dirs = sorted(d for pat in pattern.split(',') for d in glob.glob(pat))
    print(f"{len(dirs)} sample directories"
          + (f", excluding {len(exclude)}: {', '.join(sorted(exclude))}" if exclude else "")
          + "\n")

    mb, mb_sex = {}, {}
    cnv, cnv_sex = {}, {}
    widest = 0
    skipped = []
    for d in dirs:
        name = os.path.basename(os.path.dirname(d.rstrip('/')))
        if name in exclude:
            skipped.append((name, 'excluded by --exclude'))
            continue
        f1 = os.path.join(d, 'coverage_1mb.tsv')
        f2 = os.path.join(d, 'coverage_cnv.tsv')
        if os.path.exists(f1):
            rows = read_bedcov(f1)
            # per-bp depth, not summed: the last window of each chromosome is
            # truncated, and summed depth would make it look depleted in the
            # panel while the scorer (which divides by width) sees it as normal.
            acc = {(c, s): t / (e - s) for c, s, e, t in rows if e > s}
            r, sx = to_ratio(acc)
            if r:
                mb[name], mb_sex[name] = r, sx
            else:
                skipped.append((name, '1mb: zero autosomal coverage'))
        if os.path.exists(f2):
            rows = read_bedcov(f2)
            w = max((e - s) for _, s, e, _ in rows) if rows else 0
            widest = max(widest, w)
            if w > cnv_bin:
                # Re-binning upward would invent resolution that was never
                # measured, so refuse rather than quietly produce a panel.
                skipped.append((name, f'cnv: native window {w}bp > bin {cnv_bin}bp'))
                continue
            r, sx = to_ratio(rebin(rows, cnv_bin))
            if r:
                cnv[name], cnv_sex[name] = r, sx

    print(f"1Mb samples: {len(mb)}   CNV samples: {len(cnv)}   "
          f"widest native CNV window seen: {widest} bp\n")
    os.makedirs(out_dir, exist_ok=True)

    for name, samples, sexes, binbp, fn in (
            ('1Mb', mb, mb_sex, 1000000, 'coverage_1mb_panel.json'),
            (f'CNV@{cnv_bin//1000}kb', cnv, cnv_sex, cnv_bin, 'coverage_cnv_panel.json')):
        if len(samples) < min_n:
            print(f"  {name}: only {len(samples)} samples, need >= {min_n} — skipped")
            continue
        panel = summarise(samples, sexes, min_n, name, store_values=True)
        z_cut = REPORT_Z_BY_WINDOW.get(binbp, REPORT_Z)
        panel = add_frequencies(panel, samples, sexes, z_cut)
        doc = {'meta': {'n_samples': len(samples),
                        'freq_z_cut': z_cut,
                        'common_pct': COMMON_PCT,
                        'n_female': sum(1 for s in sexes.values() if s == 'F'),
                        'n_male': sum(1 for s in sexes.values() if s == 'M'),
                        'window_bp': binbp,
                        'dispersion': 'MAD * 1.4826 (robust SD estimate)',
                        'normalisation': 'per-sample autosomal median window depth',
                        'samples': sorted(samples)},
               'panel': panel}
        p = os.path.join(out_dir, fn)
        with open(p, 'w') as fh:
            json.dump(doc, fh, separators=(',', ':'))
        print(f"     wrote {p} ({os.path.getsize(p)/1e6:.1f} MB)\n")

    for n, why in skipped:
        print(f"  SKIPPED {n}: {why}")


if __name__ == '__main__':
    main()
