#!/usr/bin/env python3
"""Recompute karyotype event cohort_freq / novelty in published coverage_1mb.json.

Stage 6 used to count the cohort with the panel's |z| >= 5 flags while the
events themselves are ratio-band runs (< 0.65 / > 1.35). Where the panel is
widely spread z never reaches 5, so a region a quarter of dogs share came out
"cohort_freq 0.0, rare": chr8 74-77 Mb (the IGH locus, a breed-structured
copy-number polymorphism; canFam4's German Shepherd carries the long form) was
headlined as "a large single-copy deletion" on 176 of 1,878 published kits.
This applies stage 6's corrected annotate_freq (band=(DEL_T, DUP_T)) to the
events already in each file; ratios, pct_dogs and everything else are untouched.

The original is kept once as coverage_1mb.pre-<tag>.json beside it, which is
what cluster/republish_files.sh keys on. Idempotent: always reads the backup.

Usage: backfill_kary_freq.py <tag> <coverage_1mb_panel.json> <result_dir> [...]
"""
import json
import os
import shutil
import sys

DEL_T, DUP_T, COMMON_PCT = 0.65, 1.35, 5.0
tag, panel_path, dirs = sys.argv[1], sys.argv[2], sys.argv[3:]
pd = json.load(open(panel_path))
idx = {(c, e['start'], k): e[k] for c, es in pd['panel'].items() for e in es
       for k in ('all', 'F', 'M') if k in e}
del pd


def cohort_freq(ev, sex):
    key = 'all' if ev['chrom'] != 'chrX' else ('F' if sex == 'female' else 'M')
    fr = []
    for b in range(ev['start'], ev['end'], 1000000):
        st = idx.get((ev['chrom'], b, key))
        if not st or not st.get('n'):
            continue
        fr.append((st.get('ng', 0) if ev['direction'] == 'gain' else st.get('nl', 0)) / st['n'])
        vals = [v for v in st.get('vals') or [] if v > 0]
        if vals:
            sc = 2.0 if key == 'M' and st['median'] < 0.75 else 1.0
            hit = (sum(1 for v in vals if v * sc > DUP_T) if ev['direction'] == 'gain'
                   else sum(1 for v in vals if v * sc < DEL_T))
            fr.append(hit / len(vals))
    return round(max(fr) * 100, 1) if fr else None


changed = same = missing = 0
for d in dirs:
    cur = os.path.join(d, 'coverage_1mb.json')
    bak = os.path.join(d, f'coverage_1mb.pre-{tag}.json')
    if not os.path.exists(cur):
        missing += 1
        continue
    if not os.path.exists(bak):
        shutil.copy2(cur, bak)
    cov = json.load(open(bak))
    meta = cov.get('_meta') or {}
    evs = meta.get('events') or []
    live = (json.load(open(cur)).get('_meta') or {}).get('events') or []
    before = [(e.get('cohort_freq'), e.get('novelty')) for e in live]
    for e in evs:
        f = cohort_freq(e, meta.get('predicted_sex'))
        e['cohort_freq'] = f
        e['novelty'] = 'unknown' if f is None else 'common' if f >= COMMON_PCT else 'rare'
    if [(e['cohort_freq'], e['novelty']) for e in evs] == before:
        same += 1
        continue          # leave the file (and its mtime) alone: nothing to republish
    meta['n_rare'] = sum(1 for e in evs if e.get('novelty') == 'rare')
    meta['freq_rule'] = 'ratio-band'
    with open(cur, 'w') as fh:
        json.dump(cov, fh)
    changed += 1
print(f'{tag}: changed {changed}, unchanged {same}, missing {missing}')
