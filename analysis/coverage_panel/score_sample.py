#!/usr/bin/env python3
"""
Score coverage windows against a panel-of-normals.

    # one dog
    python3 score_sample.py PANEL.json SAMPLE/coverage_1mb.tsv

    # calibrate a threshold across the whole cohort
    python3 score_sample.py PANEL.json --cohort DIR/coverage_1mb.tsv [...]

Reports each window as a robust z: (ratio - panel median) / panel robust SD.

Choosing a threshold
--------------------
2373 windows x 96 dogs is ~228k tests. At |z| > 3 a normal distribution alone
would throw ~600 false flags, which would bury every real event. But coverage
ratios are not normal — they have heavy tails from mappability, GC and repeat
content — so a theoretical correction is the wrong instrument too.

--cohort therefore reports the OBSERVED tail. Read it as: "at |z| > 5, the
median dog has N flagged windows". Pick the threshold where the median dog has
approximately zero, because most dogs should not carry large CNVs; whatever
remains above it is either real or a property of that window in that dog.

Windows the panel cannot judge
------------------------------
A window whose panel SD is ~0 makes every deviation look infinite, and one with
a huge SD can never be called at all. Both are excluded and counted, because
silently scoring them produces confident nonsense.
"""
import json
import statistics
import sys

# MIN_SD guards against dividing by a spread that is really zero. It was 0.02,
# chosen when the only panel available came from 8 replicate-heavy samples. On
# the real 96-dog panel that excluded 1,702 of 2,248 windows — three quarters of
# the genome — because a well-measured 1Mb window genuinely has SD ~0.015. With
# n=96 the sampling error on that estimate is ~0.001, so it is a measurement,
# not noise. The floor now only catches degenerate windows.
#
# For scale: observed SDs match Poisson counting noise at both resolutions
# (1Mb ~13k reads -> 0.009 expected vs 0.015 seen; 50kb ~650 reads -> 0.039 vs
# 0.044), so the panel is measuring the genome rather than batch effects.
MIN_SD = 0.004
MAX_SD = 0.30        # above this the window is uninformative (mappability etc.)
REPORT_Z = 5.0       # default per-window reporting threshold


def read_coverage(path, bin_bp=None):
    """Per-window mean depth, keyed by (chrom, start).

    With bin_bp set, the sample is first re-binned onto that fixed grid. This is
    required for the CNV panel: Stage 5 sizes CNV windows per dog as
    50000/mean_depth, so a sample's native windows never line up with a shared
    grid and every one of them would score as uncallable. Re-binning splits each
    source window's summed depth across target bins by overlap — the same
    uniform-depth assumption the original binning makes.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 4:
                continue
            chrom, start, end, total = f[0], int(f[1]), int(f[2]), float(f[3])
            if end > start:
                rows.append((chrom, start, end, total))
    if bin_bp is None:
        return {(c, s): t / (e - s) for c, s, e, t in rows}

    widest = max((e - s) for c, s, e, t in rows) if rows else 0
    if widest > bin_bp:
        raise SystemExit(
            f"ERROR: {path} has native windows up to {widest}bp, wider than the "
            f"panel's {bin_bp}bp grid. Re-binning upward would invent resolution "
            f"that was never measured.")
    acc, span = {}, {}
    for chrom, s, e, total in rows:
        per_bp = total / (e - s)
        for b in range(s // bin_bp, (e - 1) // bin_bp + 1):
            lo, hi = max(s, b * bin_bp), min(e, (b + 1) * bin_bp)
            if hi > lo:
                k = (chrom, b * bin_bp)
                acc[k] = acc.get(k, 0.0) + per_bp * (hi - lo)
                span[k] = span.get(k, 0.0) + (hi - lo)
    # divide by the length actually covered, matching build_panels.rebin
    return {k: v / span[k] for k, v in acc.items() if span[k] > 0}


def normalise(depths):
    vals = [d for (c, _), d in depths.items() if c != 'chrX' and d > 0]
    base = statistics.median(vals) if vals else 0.0
    if base <= 0:
        return None, None
    ratios = {k: v / base for k, v in depths.items()}
    x = [v for (c, _), v in ratios.items() if c == 'chrX' and v > 0]
    sex = 'F' if (statistics.median(x) if x else 1.0) >= 0.75 else 'M'
    return ratios, sex


def load_panel(path):
    with open(path) as fh:
        doc = json.load(fh)
    index = {}
    for chrom, entries in doc['panel'].items():
        for e in entries:
            for key in ('all', 'F', 'M'):
                if key in e:
                    index[(chrom, e['start'], key)] = e[key]
    return doc.get('meta', {}), index


def score(ratios, sex, index):
    """-> (list of (chrom, start, ratio, z), n_uncallable)"""
    out, uncallable = [], 0
    for (chrom, start), ratio in ratios.items():
        key = sex if chrom == 'chrX' else 'all'
        stats = index.get((chrom, start, key))
        if stats is None:
            uncallable += 1
            continue
        sd = stats['mad_sd']
        if sd < MIN_SD or sd > MAX_SD:
            uncallable += 1
            continue
        out.append((chrom, start, ratio, (ratio - stats['median']) / sd))
    return out, uncallable


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__.strip())
    panel_path = sys.argv[1]
    args = sys.argv[2:]
    cohort = args[0] == '--cohort'
    if cohort:
        args = args[1:]

    meta, index = load_panel(panel_path)
    bin_bp = meta.get('window_bp')
    if bin_bp == 1000000:
        bin_bp = None          # 1Mb windows already share a fixed grid
    print(f"panel: {meta.get('n_samples','?')} samples, "
          f"{meta.get('n_windows','?')} windows "
          f"({meta.get('n_female','?')}F/{meta.get('n_male','?')}M)")

    if not cohort:
        ratios, sex = normalise(read_coverage(args[0], bin_bp))
        if ratios is None:
            sys.exit("ERROR: zero autosomal coverage")
        scored, uncallable = score(ratios, sex, index)
        flagged = sorted((w for w in scored if abs(w[3]) >= REPORT_Z),
                         key=lambda w: -abs(w[3]))
        print(f"sex: {sex}   scored: {len(scored)}   uncallable: {uncallable}")
        print(f"flagged at |z| >= {REPORT_Z}: {len(flagged)}")
        for chrom, start, ratio, z in flagged[:40]:
            kind = 'GAIN' if z > 0 else 'LOSS'
            print(f"  {kind} {chrom}:{start//1000000}Mb  ratio={ratio:.3f}  z={z:+.1f}")
        return

    # ── cohort calibration ───────────────────────────────────────
    thresholds = [3, 4, 5, 6, 8, 10]
    per_dog = []
    recurrence = {}
    unscoreable = []
    for path in args:
        # Refusing to re-bin upward is correct for a single sample, but in a
        # cohort one such sample must not abort the other 94 — skip and report.
        try:
            cov = read_coverage(path, bin_bp)
        except SystemExit as e:
            unscoreable.append((path, str(e).split(' has ')[-1].strip()))
            continue
        ratios, sex = normalise(cov)
        if ratios is None:
            continue
        scored, uncallable = score(ratios, sex, index)
        counts = [sum(1 for w in scored if abs(w[3]) >= t) for t in thresholds]
        for c, st, _, z in scored:
            if abs(z) >= 5:
                recurrence[(c, st)] = recurrence.get((c, st), 0) + 1
        per_dog.append((path, sex, counts, uncallable))

    for p_, why in unscoreable:
        print(f"  SKIPPED {p_}: {why}")
    print(f"\ncalibration across {len(per_dog)} dogs — flagged windows per dog")
    print("      " + "".join(f"|z|>={t:<7}" for t in thresholds))
    for label, agg in (('median', statistics.median), ('max', max)):
        row = [agg([d[2][i] for d in per_dog]) for i in range(len(thresholds))]
        print(f"{label:<6}" + "".join(f"{v:<11.0f}" for v in row))

    print("\nPick the lowest threshold whose MEDIAN dog is ~0: most dogs should")
    print("carry no large CNVs, so a nonzero median is the false-positive rate.")

    # Windows flagged in most dogs are not per-dog events — they are windows the
    # panel cannot score. Recurrence separates them from real CNVs cleanly: a
    # true copy-number variant appears in a minority of dogs, an unscoreable
    # window appears in nearly all of them.
    print("\n--- recurrent windows (flagged in many dogs at |z| >= 5) ---")
    rec = sorted(((n, w) for w, n in recurrence.items() if n >= 0.5 * len(per_dog)),
                 reverse=True)
    print(f"  {len(rec)} windows flagged in >=50% of dogs; they account for "
          f"{sum(n for n, _ in rec)} of {sum(d[2][2] for d in per_dog)} total flags")
    for n, (c, st) in rec[:15]:
        pm = index.get((c, st, 'all'))
        extra = f"panel median {pm['median']:.3f} sd {pm['mad_sd']:.4f}" if pm else ''
        print(f"     {c}:{st//10**6}Mb  flagged in {n}/{len(per_dog)} dogs   {extra}")

    noisy = sorted(per_dog, key=lambda d: -d[2][2])[:5]
    print(f"\nmost-flagged dogs at |z| >= {thresholds[2]} (check these for QC problems):")
    for path, sex, counts, uncallable in noisy:
        print(f"  {counts[2]:>5} windows  {sex}  uncallable={uncallable}  {path}")


if __name__ == '__main__':
    main()
