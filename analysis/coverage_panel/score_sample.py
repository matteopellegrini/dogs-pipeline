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

MIN_SD = 0.02        # below this the panel is too tight to divide by
MAX_SD = 0.30        # above this the window is uninformative (mappability etc.)
REPORT_Z = 5.0       # default per-window reporting threshold


def read_coverage(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 4:
                continue
            chrom, start, end, total = f[0], int(f[1]), int(f[2]), float(f[3])
            if end > start:
                out[(chrom, start)] = total / (end - start)
    return out


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
    print(f"panel: {meta.get('n_samples','?')} samples, "
          f"{meta.get('n_windows','?')} windows "
          f"({meta.get('n_female','?')}F/{meta.get('n_male','?')}M)")

    if not cohort:
        ratios, sex = normalise(read_coverage(args[0]))
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
    for path in args:
        ratios, sex = normalise(read_coverage(path))
        if ratios is None:
            continue
        scored, uncallable = score(ratios, sex, index)
        counts = [sum(1 for w in scored if abs(w[3]) >= t) for t in thresholds]
        per_dog.append((path, sex, counts, uncallable))

    print(f"\ncalibration across {len(per_dog)} dogs — flagged windows per dog")
    print("      " + "".join(f"|z|>={t:<7}" for t in thresholds))
    for label, agg in (('median', statistics.median), ('max', max)):
        row = [agg([d[2][i] for d in per_dog]) for i in range(len(thresholds))]
        print(f"{label:<6}" + "".join(f"{v:<11.0f}" for v in row))

    print("\nPick the lowest threshold whose MEDIAN dog is ~0: most dogs should")
    print("carry no large CNVs, so a nonzero median is the false-positive rate.")

    noisy = sorted(per_dog, key=lambda d: -d[2][2])[:5]
    print(f"\nmost-flagged dogs at |z| >= {thresholds[2]} (check these for QC problems):")
    for path, sex, counts, uncallable in noisy:
        print(f"  {counts[2]:>5} windows  {sex}  uncallable={uncallable}  {path}")


if __name__ == '__main__':
    main()
