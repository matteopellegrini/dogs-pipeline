#!/usr/bin/env python3
"""
Build a coverage panel-of-normals from many samples' coverage_1mb.tsv files.

    python3 analysis/coverage_panel/build_panel.py OUT.json DIR [DIR ...]

Each DIR is a sample's analysis directory (containing coverage_1mb.tsv), or a
glob of them. The output JSON gives, for every 1Mb window, the expected
normalised coverage and how much it varies across the cohort — so a new dog's
window can be called unusual against a real distribution instead of a single
reference ratio.

Why median/MAD and not mean/SD
------------------------------
The cohort is unscreened: some dogs genuinely carry large CNVs, and a single
sample with a 3x amplification would inflate that window's SD enough to hide
the same event in everyone else. The median is unmoved by a minority of
outliers, and MAD*1.4826 estimates the SD of the underlying normal part. That
matters because the panel's whole job is to make outliers detectable.

Sex and chrX
------------
chrX coverage is ~1x in females and ~0.5x in males, so pooling them produces a
bimodal distribution whose median describes nobody. Sex is inferred from each
sample's own chrX:autosome ratio and chrX statistics are computed per sex.
Autosomes are pooled across all samples.
"""
import glob
import json
import os
import statistics
import sys

MAD_TO_SD = 1.4826       # scale MAD so it estimates the SD of a normal
MIN_SAMPLES = 5          # below this a per-window distribution is meaningless
SEX_CALL_THRESHOLD = 0.75  # chrX:autosome ratio; ~1.0 female, ~0.5 male


def read_coverage(path):
    """coverage_1mb.tsv -> {(chrom, start): mean depth over the window}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 4:
                continue
            chrom, start, end, total = f[0], int(f[1]), int(f[2]), float(f[3])
            width = end - start
            if width > 0:
                out[(chrom, start)] = total / width
    return out


def autosomal_median(depths):
    vals = [d for (c, _), d in depths.items() if c != 'chrX' and d > 0]
    return statistics.median(vals) if vals else 0.0


def mad(values, med):
    if len(values) < 2:
        return 0.0
    return statistics.median([abs(v - med) for v in values])


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__.strip())
    out_path, patterns = sys.argv[1], sys.argv[2:]

    dirs = []
    for p in patterns:
        dirs.extend(sorted(glob.glob(p)) if any(ch in p for ch in '*?[') else [p])

    # ── load and normalise each sample ───────────────────────────
    samples = {}       # name -> {(chrom,start): ratio}
    sex = {}           # name -> 'M' | 'F'
    skipped = []
    for d in dirs:
        tsv = d if d.endswith('.tsv') else os.path.join(d, 'coverage_1mb.tsv')
        if not os.path.exists(tsv):
            skipped.append((d, 'no coverage_1mb.tsv'))
            continue
        name = os.path.basename(os.path.dirname(os.path.abspath(tsv)))
        if name == 'analysis':
            name = os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(tsv))))

        depths = read_coverage(tsv)
        base = autosomal_median(depths)
        if base <= 0:
            skipped.append((d, 'zero autosomal coverage'))
            continue

        ratios = {k: v / base for k, v in depths.items()}
        x = [v for (c, _), v in ratios.items() if c == 'chrX' and v > 0]
        x_ratio = statistics.median(x) if x else 1.0
        sex[name] = 'F' if x_ratio >= SEX_CALL_THRESHOLD else 'M'
        samples[name] = ratios

    if len(samples) < MIN_SAMPLES:
        sys.exit(f"ERROR: only {len(samples)} usable samples; need >= {MIN_SAMPLES}")

    # ── per-window statistics ────────────────────────────────────
    windows = set()
    for r in samples.values():
        windows.update(r)

    panel = {}
    n_windows = 0
    for chrom, start in sorted(windows, key=lambda w: (w[0], w[1])):
        if chrom == 'chrX':
            groups = {s: [samples[n][(chrom, start)]
                          for n in samples if sex[n] == s and (chrom, start) in samples[n]]
                      for s in ('F', 'M')}
        else:
            groups = {'all': [samples[n][(chrom, start)]
                              for n in samples if (chrom, start) in samples[n]]}

        entry = {'start': start}
        for key, vals in groups.items():
            if len(vals) < MIN_SAMPLES:
                continue
            med = statistics.median(vals)
            entry[key] = {
                'median': round(med, 4),
                'mad_sd': round(mad(vals, med) * MAD_TO_SD, 4),
                'n': len(vals),
            }
        if len(entry) > 1:
            panel.setdefault(chrom, []).append(entry)
            n_windows += 1

    result = {
        'meta': {
            'n_samples': len(samples),
            'n_female': sum(1 for s in sex.values() if s == 'F'),
            'n_male': sum(1 for s in sex.values() if s == 'M'),
            'n_windows': n_windows,
            'window_bp': 1000000,
            'dispersion': 'MAD * 1.4826 (robust SD estimate)',
            'normalisation': 'per-sample autosomal median window depth',
            'samples': sorted(samples),
        },
        'panel': panel,
    }
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, 'w') as fh:
        json.dump(result, fh, separators=(',', ':'))

    print(f"samples   : {len(samples)} "
          f"({result['meta']['n_female']}F / {result['meta']['n_male']}M)")
    print(f"windows   : {n_windows} across {len(panel)} chromosomes")
    print(f"written   : {out_path}")
    for d, why in skipped:
        print(f"  SKIPPED {d}: {why}")

    # A quick sense of how tight the panel is. Windows with a large spread are
    # the ones where no per-sample call will ever be confident.
    spreads = [w['all']['mad_sd'] for c in panel for w in panel[c] if 'all' in w]
    if spreads:
        spreads.sort()
        print(f"autosomal robust SD: median {spreads[len(spreads)//2]:.3f}, "
              f"90th pct {spreads[int(len(spreads)*0.9)]:.3f}, max {spreads[-1]:.3f}")


if __name__ == '__main__':
    main()
