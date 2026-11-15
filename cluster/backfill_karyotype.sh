#!/bin/bash
# Backfill the 1Mb karyotype (coverage_1mb.json) for every dog against the
# pooled panel-of-normals (1,283 Illumina dogs, built 2026-09-03 by
# analysis/coverage_panel/build_panels.py job 14650653), replacing scores
# computed against the original 93-dog panel.
#
# The python below is a VERBATIM extract of the stage-6a block of
# run_dog_pipeline.sh (lines 569-923 at commit time) — same heredoc
# interpolation, same logic — so backfilled files are identical to what a
# full stage-6 rerun would produce, without touching mito (6d needs a FASTQ
# remap) or the 50kb CNV caller (6e, already backfilled separately).
# Only coverage_1mb.json is rewritten; qc_result.json and everything after
# it in stage 6 are left alone.
set -u
D=/u/project/pellegrini/gkislik/dogs
COV_PANEL_1MB=$D/reference_panel/coverage_1mb_panel.json
COV_Z_1MB=5.0
CNV_WINDOW=0   # referenced by the extracted block but unused before the 1Mb write
DATA_PYTHON=$(command -v python3)
# Non-tty ssh sessions default to an ASCII locale, and the extracted block
# prints a literal arrow character — without this the block dies mid-run.
export PYTHONIOENCODING=utf-8

run_one() {
  local OUT=$1 PUB=$2
  [[ -s "$OUT/coverage_1mb.tsv" ]] || { echo "SKIP no coverage_1mb.tsv: $OUT"; return 1; }
  [[ -d "$PUB" ]] || { echo "SKIP no results dir: $PUB"; return 1; }
"$DATA_PYTHON" - << PYEOF
import json, collections, statistics, subprocess, re, os

tsv_1mb  = "$OUT/coverage_1mb.tsv"
tsv_cnv  = "$OUT/coverage_cnv.tsv"
pub      = "$PUB"
cnv_win  = int("$CNV_WINDOW")

def chrom_key(c):
    c2 = c.replace('chr','')
    return (0, int(c2)) if c2.isdigit() else (1, c2)

def load_tsv(path):
    rows = []
    with open(path) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 4: continue
            chrom, start, end, bases = cols[0], int(cols[1]), int(cols[2]), int(cols[3])
            size = end - start
            if size <= 0: continue
            rows.append((chrom, start, end, bases / size))
    return rows

COMMON_PCT = 5.0    # cohort share at or above which a region is "common"

def annotate_freq(events, idx, bin_bp, sex):
    """Tag each event with how often the reference cohort is flagged there.

    Recurrence analysis over 93 dogs found chr19:21.3-21.5Mb flagged as a gain
    in 42 of them and chr5:85.9Mb in 41, while NO bin was flagged in more than
    half. That pattern is the point: a reference artifact would hit nearly every
    dog, since all dogs map to the same reference, whereas hitting a fraction is
    what a real copy-number polymorphism looks like. The events are genuine and
    still not findings -- a region carried by 45% of dogs is common, and
    presenting it as a discovery about one dog is a category error.

    Takes the MAXIMUM frequency across the event's bins, not the peak bin's.
    Conservative on purpose: if any part of an event coincides with a region
    dogs commonly vary at, it is most likely that polymorphism, and the cost of
    wrongly calling something rare (alarming an owner) exceeds the cost of
    wrongly calling it common (a quieter listing).
    """
    for ev in events:
        key = 'all' if ev['chrom'] != 'chrX' else ('F' if sex == 'female' else 'M')
        fr = []
        for b in range(ev['start'], ev['end'], bin_bp):
            st = idx.get((ev['chrom'], b, key))
            if st and st.get('n'):
                c = st.get('ng', 0) if ev['direction'] == 'gain' else st.get('nl', 0)
                fr.append(c / st['n'])
        ev['cohort_freq'] = round(max(fr) * 100, 1) if fr else None
        ev['novelty'] = ('unknown' if ev['cohort_freq'] is None
                         else 'common' if ev['cohort_freq'] >= COMMON_PCT else 'rare')
    return events

# --- coverage_1mb.json (karyotype) ---
import statistics as _stats, json as _json
windows_1mb = load_tsv(tsv_1mb)
data = collections.defaultdict(dict)
for chrom, start, end, depth in windows_1mb:
    data[chrom][start // 1000000] = round(depth, 4)
raw = {}
for chrom in sorted(data, key=chrom_key):
    pts = data[chrom]
    raw[chrom] = [pts.get(i, 0) for i in range(max(pts) + 1)]

# Genome-wide median (autosomal) for ratio normalisation
auto_depths = [d for c, arr in raw.items() if c != 'chrX' for d in arr if d > 0]
kiki_median = _stats.median(auto_depths)

# Sex determination: male = one X (depth ~0.5× autosomal), female = two X (~1.0×)
x_depths = [d for d in raw.get('chrX', []) if d > 0]
x_auto_ratio = _stats.median(x_depths) / kiki_median if x_depths else 1.0
predicted_sex = 'male' if x_auto_ratio < 0.75 else 'female'
# For ratio normalisation of chrX: divide by 0.5 for males so hemizygous X → 1.0
chrx_norm = 0.5 if predicted_sex == 'male' else 1.0
print(f"Sex determination: chrX/auto ratio = {x_auto_ratio:.3f} → {predicted_sex} (chrX norm divisor: {chrx_norm})")

# 96-dog panel-of-normals: per-window median and robust SD, so a window can be
# reported as "N SD from 96 dogs" rather than just "looks low".
cov_panel_path = "$COV_PANEL_1MB"
Z_CUT = float("$COV_Z_1MB")
try:
    with open(cov_panel_path) as _f:
        _pd = _json.load(_f)
    panel_idx = {(c, e['start'], k): e[k]
                 for c, es in _pd['panel'].items() for e in es
                 for k in ('all', 'F', 'M') if k in e}
    panel_n = _pd.get('meta', {}).get('n_samples')
    # Which panel produced this report. The panels are no longer tracked in git
    # (they are regenerated build artifacts, and the CNV one is 29MB), so the
    # identifying information lives here instead — enough to tell whether two
    # reports were scored against the same reference.
    import hashlib as _hl
    with open(cov_panel_path, 'rb') as _pf:
        panel_sha = _hl.sha256(_pf.read()).hexdigest()[:12]
except Exception as _e:
    panel_idx, panel_n = {}, None
    print(f"WARNING: no coverage panel at {cov_panel_path} ({_e}) — "
          "windows will carry no significance")

result = {}
for chrom, arr in raw.items():
    # Tooltip depth tracks. These used to be copied from a legacy reference
    # file that lost its 'cosmo'/'panel' keys in a rebuild, which silently
    # zero-padded both series in every report since. The sample's own depth
    # is already in hand, and the panel's expected depth is derived per
    # window from the panel stats in the loop below.
    cosmo_arr = arr
    panel_arr = []   # per window: panel median ratio × sample median depth, or None
    n = len(arr)
    if chrom == 'chrX' and predicted_sex == 'male':
        # PAR1 windows are diploid in males (depth ≈ autosomal); non-PAR windows are hemizygous (depth ≈ 0.5×).
        # Use per-window normalization: above 0.75× autosomal → PAR1 → divide by kiki_median;
        # below 0.75× autosomal → non-PAR → divide by kiki_median × 0.5, so both map to ratio 1.0.
        par_thresh = kiki_median * 0.75
        ratio_arr, scale_arr = [], []
        for d in arr:
            if d == 0:
                ratio_arr.append(0.0); scale_arr.append(1.0)
            elif d >= par_thresh:
                ratio_arr.append(round(d / kiki_median, 4))        # PAR1: diploid norm
                scale_arr.append(1.0)
            else:
                ratio_arr.append(round(d / (kiki_median * 0.5), 4)) # non-PAR: hemizygous norm
                scale_arr.append(2.0)
    else:
        ratio_arr = [round(d / kiki_median, 4) if d > 0 else 0.0 for d in arr]
        scale_arr = [1.0] * len(arr)
    # Robust z against the panel. Computed from the UNADJUSTED ratio
    # (depth / autosomal median), because that is the space the panel medians
    # live in — ratio_arr rescales male chrX so hemizygous reads as 1.0, which
    # would not line up with a male chrX panel median of ~0.5.
    pkey = 'all' if chrom != 'chrX' else ('F' if predicted_sex == 'female' else 'M')
    # Alongside this dog's z, carry how often the REFERENCE COHORT itself varies
    # at each window. A reader looking at a dip or a bump cannot tell whether it
    # matters; "34 of 93 healthy dogs are also low here" answers that directly,
    # and answers it for every window rather than only the flagged ones.
    #
    # Deliberately computed even where this dog's own window is unscoreable:
    # the frequency describes the cohort, not the sample, so "do other dogs vary
    # here?" has an answer either way.
    # The band is the more useful of the two for reading the chart. At 1Mb the
    # flag frequency is 0 almost everywhere (megabase averaging washes out focal
    # CNVs — only 8 of 2,498 windows are flagged in >=5% of the cohort), so on
    # its own it cannot tell a reader whether the bump they are looking at is
    # ordinary. The spread of the reference dogs at that window can: drawn
    # behind the bars, a bar inside the band is normal variation and a bar
    # outside it is not, with no numbers to interpret.
    #
    # Emitted in the SAME units as ratio_arr, hence scale_arr: on a male chrX
    # the ratio is rescaled so hemizygous reads as 1.0, while the panel median
    # sits at ~0.5, and a band in raw panel units would sit at half height and
    # look alarming on every male.
    BAND_SD = 2.0
    z_arr, pct_arr, lo_arr, hi_arr = [], [], [], []
    for i, d in enumerate(arr):
        st = panel_idx.get((chrom, i * 1000000, pkey))
        sc = scale_arr[i] if i < len(scale_arr) else 1.0
        # Expected depth of a typical panel dog here, expressed at this
        # sample's sequencing depth so Sample and Panel are comparable.
        panel_arr.append(round(st['median'] * kiki_median, 2) if st else None)
        # THE question a reader has: of the reference dogs, how many show this
        # much coverage here, or more? Counted from the panel's actual values at
        # this window, in the direction this dog deviates. A threshold-crossing
        # frequency cannot answer it — that asks how many dogs cross a fixed
        # line, not how many are as extreme as this one — and median/spread
        # cannot either, because the answer lives in the tail.
        vals = st.get('vals') if st else None
        if vals and d > 0:
            obs = d / kiki_median
            if obs >= st['median']:
                cnt = sum(1 for v in vals if v >= obs)
            else:
                cnt = sum(1 for v in vals if v <= obs)
            pct_arr.append(round(100 * cnt / len(vals), 1))
        else:
            pct_arr.append(None)
        if st and 0.004 <= st['mad_sd'] <= 0.30:
            lo_arr.append(round(max(0.0, st['median'] - BAND_SD * st['mad_sd']) * sc, 4))
            hi_arr.append(round((st['median'] + BAND_SD * st['mad_sd']) * sc, 4))
        else:
            lo_arr.append(None)
            hi_arr.append(None)
        if not st or d <= 0 or not (0.004 <= st['mad_sd'] <= 0.30):
            z_arr.append(None)
            continue
        z_arr.append(round(((d / kiki_median) - st['median']) / st['mad_sd'], 2))
    result[chrom] = {
        'cosmo': [round(v, 2) for v in cosmo_arr],
        'panel': panel_arr,
        'ratio': ratio_arr,
        'z': z_arr,
        'pct_dogs': pct_arr,
        'band_lo': lo_arr,
        'band_hi': hi_arr,
        'band_sd': BAND_SD,
    }

# Merge runs of adjacent flagged windows into events. "3 regions of unusual
# coverage" is a statement a dog owner can act on; "7 flagged windows" is not.
# Regions are the runs of windows the CHART draws in red or orange — i.e. the
# ones a reader can actually see — not runs defined by a z statistic.
#
# z was the wrong instrument here. It assumes the reference dogs are roughly
# normally distributed at a window, and at the scale of a real deletion they are
# not: at chr8's end the 93 panel dogs span 0.43-1.47, so a dog near the bottom
# scores z -1.1 and is "not flagged" while being in the lowest few percent of
# the population. That produced a chart shouting DELETION in red at windows the
# table said nothing about, and a tooltip quoting 2-4% beside a table row
# claiming "not seen".
#
# The honest statistic is the empirical one, and it is exactly a minor allele
# frequency: of the reference dogs, what share carry this much loss or gain
# here. No distributional assumption, and it means the same thing whether the
# deletion is 1Mb or 40.
DEL_T, DUP_T = 0.65, 1.35
events = []
for chrom in sorted(result):
    ratios = result[chrom].get('ratio') or []
    pcts   = result[chrom].get('pct_dogs') or []
    run = None
    for i, r in enumerate(list(ratios) + [None]):
        # ratio 0 means no coverage measured, not a deletion
        d = None
        if r is not None and r > 0:
            if r < DEL_T: d = 'loss'
            elif r > DUP_T: d = 'gain'
        if d and run is None:
            run = [i, i, d]
        elif d and run and run[2] == d:
            run[1] = i
        elif run is not None:
            lo, hi, direction = run
            seg = [pcts[j] for j in range(lo, hi + 1) if j < len(pcts) and pcts[j] is not None]
            events.append({'chrom': chrom, 'start': lo * 1000000,
                           'end': (hi + 1) * 1000000,
                           'windows': hi - lo + 1,
                           'direction': direction,
                           # range across the region, because that is what a
                           # reader sees when they hover along it
                           'pct_min': min(seg) if seg else None,
                           'pct_max': max(seg) if seg else None,
                           'pct_dogs': min(seg) if seg else None})
            run = [i, i, d] if d else None
events.sort(key=lambda e: (e['pct_dogs'] if e['pct_dogs'] is not None else 999, -e['windows']))

# Whole-chromosome events are aneuploidy, not copy-number variation, and must
# not be reported as "a gain on chrX". DOGS-Gen-111 is the case that forced
# this: 124 of her 125 chrX windows were elevated (X:autosome 1.461, i.e. three
# copies), which the event merger split into two large gains purely because one
# unscoreable window interrupted the run. Her breed-matched, same-batch control
# DOGS-Gen-110 measured 1.001, so the signal is hers and not the assay's.
#
# Detect per chromosome rather than by widening the merge gap: a chromosome
# qualifies when most of its SCOREABLE windows are flagged the same way. Using a
# fraction rather than a count makes it work on chr38 as well as chr1.
ANEUPLOIDY_FRAC = 0.80
aneuploidy = []
aneuploid_chroms = set()
for chrom in sorted(result):
    zs = result[chrom].get('z') or []
    scoreable = [i for i, z in enumerate(zs) if z is not None]
    if len(scoreable) < 10:
        continue                      # too little to judge a whole chromosome
    up   = [i for i in scoreable if zs[i] >=  Z_CUT]
    down = [i for i in scoreable if zs[i] <= -Z_CUT]
    hit, direction = (up, 'gain') if len(up) >= len(down) else (down, 'loss')
    frac = len(hit) / len(scoreable)
    if frac < ANEUPLOIDY_FRAC:
        continue
    # Copy number straight from the unadjusted ratio: 2 copies is the diploid
    # baseline, so a male's single X (~0.5) and a trisomy (~1.5) both fall out
    # of the same arithmetic without special-casing.
    rr = [raw[chrom][i] / kiki_median for i in scoreable if raw[chrom][i] > 0]
    med_ratio = _stats.median(rr) if rr else 0.0
    copies = round(med_ratio * 2)
    expected = 1 if (chrom == 'chrX' and predicted_sex == 'male') else 2
    if copies == expected:
        continue                      # flagged but not a copy-number change
    aneuploidy.append({
        'chrom': chrom, 'direction': direction,
        'call': {0: 'nullisomy', 1: 'monosomy',
                 3: 'trisomy', 4: 'tetrasomy'}.get(copies, f'{copies}-copy'),
        'copies_est': copies, 'expected_copies': expected,
        'median_ratio': round(med_ratio, 3),
        'fraction_flagged': round(frac, 3),
        'windows_flagged': len(hit), 'windows_scoreable': len(scoreable),
        'peak_z': max((zs[i] for i in hit), key=abs)})
    aneuploid_chroms.add(chrom)

# An aneuploid chromosome is fully described by its aneuploidy entry, so its
# windows must not ALSO appear in the event list — otherwise a trisomy reads as
# two large gains, which is both wrong in kind and needlessly alarming.
if aneuploid_chroms:
    events = [e for e in events if e['chrom'] not in aneuploid_chroms]
for a in aneuploidy:
    print(f"ANEUPLOIDY: {a['chrom']} {a['call']} — {a['copies_est']} copies "
          f"(expected {a['expected_copies']}), median ratio {a['median_ratio']}, "
          f"{a['windows_flagged']}/{a['windows_scoreable']} windows flagged")

annotate_freq(events, panel_idx, 1000000, predicted_sex)

# Per-dog confidence. A dog with 242 events should not produce 242 findings —
# it should say the coverage is not clean enough to call structural variation.
#
# Thresholds read off the 93-dog cohort, which separates cleanly: 93 dogs at
# 0-9 events, three at 55-270. A six-fold gap, so 20 is comfortably inside it.
# Depth is NOT the criterion — the lowest-coverage dog in the cohort (1.00x)
# produced 9 events, while a 2.40x dog produced 89 of which 85 were gains. That
# asymmetry is the second signal: real copy-number variation is not
# one-directional, so a gain-dominated profile is a mapping artifact even when
# the count alone would pass.
_gains = sum(1 for e in events if e['direction'] == 'gain')
_losses = len(events) - _gains
_reasons = []
# The old thresholds were calibrated against z-defined events and do not carry
# over. In particular the gain-dominated rule now misfires on clean samples:
# with regions defined by the ratio bands, gains above 1.35 are common at
# chromosome ends while losses below 0.65 are rare, so gain-dominance is the
# NORMAL pattern — Kiki2, platform-matched and clean, comes out 11 gains to 1
# loss. That rule is dropped rather than retuned; it was an artifact signature
# for a statistic no longer in use.
#
# What remains is a volume check. A clean dog sits around a dozen regions
# (Kiki2 12, Luna Genotek 13), so this only catches samples an order of
# magnitude worse. PROVISIONAL: wants recalibrating across the 96-dog cohort
# on this new basis before it is trusted as a QC gate.
if len(events) > 40:
    _reasons.append(f'{len(events)} flagged regions (a clean sample has roughly a dozen)')
cov_confidence = 'low' if _reasons else 'ok'
# Carry the verdict onto each aneuploidy entry. "Your dog has three copies of
# the X chromosome" is a strong claim to render, and whoever renders it should
# not have to remember to cross-check a sibling field before doing so.
for a in aneuploidy:
    a['confidence'] = cov_confidence
if _reasons:
    print(f"coverage confidence LOW: {'; '.join(_reasons)}")
print(f"coverage events at |z| >= {Z_CUT}: {len(events)} "
      f"({_gains} gain / {_losses} loss), confidence={cov_confidence}")

result['_meta'] = {'predicted_sex': predicted_sex,
                   'panel_sha': panel_sha,
                   'chrx_auto_ratio': round(x_auto_ratio, 3),
                   'panel_n': panel_n, 'z_threshold': Z_CUT,
                   'coverage_confidence': cov_confidence,
                   'confidence_reasons': _reasons,
                   'n_gain': _gains, 'n_loss': _losses,
                   'n_rare': sum(1 for e in events if e.get('novelty') == 'rare'),
                   'common_pct': COMMON_PCT,
                   'events': events,
                   'aneuploidy': aneuploidy}
with open(f'{pub}/coverage_1mb.json', 'w') as f:
    json.dump(result, f)
print(f"coverage_1mb.json: {len(result)-1} chromosomes, median depth {kiki_median:.2f}×, sex={predicted_sex}")
PYEOF
}

if [[ $# -ge 2 ]]; then
  # single-dog mode: backfill_karyotype.sh <out_dir> <pub_dir>
  run_one "$1" "$2"
  exit $?
fi

ok=0; fail=0
for w in "$D"/work_prosper/pk-*/analysis; do
  s=$(basename "$(dirname "$w")")
  if run_one "$w" "$D/results_prosper/$s" > /dev/null 2>"$D/logs/backfill_kary.err.tmp"; then
    ok=$((ok+1))
  else
    fail=$((fail+1)); echo "FAIL $s"; cat "$D/logs/backfill_kary.err.tmp"
  fi
done
for w in "$D"/work/*/analysis; do
  s=$(basename "$(dirname "$w")")
  low=$(echo "$s" | tr '[:upper:]' '[:lower:]')
  [[ -d "$D/results/$low" ]] || continue
  if run_one "$w" "$D/results/$low" > /dev/null 2>"$D/logs/backfill_kary.err.tmp"; then
    ok=$((ok+1))
  else
    fail=$((fail+1)); echo "FAIL $s"; cat "$D/logs/backfill_kary.err.tmp"
  fi
done
echo "backfill done: ok=$ok fail=$fail"
