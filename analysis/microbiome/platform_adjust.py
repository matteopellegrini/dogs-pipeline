#!/usr/bin/env python3
"""Place a new-process (MGI / Complete Genomics) swab on the old-lab cohort's scale.

The cohort references were sequenced on Illumina; new-process swabs differ in composition, and there are
too few new-process dogs for a reference of their own. So: shift = average new-process offset from the
reference swabs in reference_json/microbiome_platform_offset.json (dog means, depth-matched, leaving out
the dog being scored), spread = the cohort. The adjusted value gets a cohort percentile, reported only as
a band (low <= 25th < typical < 75th <= high) because a single swab varies a lot from day to day.

Rewrites in place, keeping the originals:
  microbiome_health_result.json  pathobiont_percentile (+ _unadjusted), pathobiont_band, richness and
                                 Shannon percentiles nulled (not comparable across processes)
  oral_function_result.json      per-module percentile/level/tone (+ percentile_unadjusted), note
Old-process samples (read length in the cohort panel's list) are left untouched.

usage: platform_adjust.py PUB --sample NAME [--force] [--offset F] [--panel F] [--oral-panel F]
"""
import argparse, bisect, json, math, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ap = argparse.ArgumentParser()
ap.add_argument('pub'); ap.add_argument('--sample', required=True)
ap.add_argument('--force', action='store_true', help='treat as new-process whatever the read length')
ap.add_argument('--offset', default=os.path.join(ROOT, 'reference_json', 'microbiome_platform_offset.json'))
ap.add_argument('--panel', default=os.path.join(ROOT, 'reference_panel', 'microbiome_panel.json'))
ap.add_argument('--oral-panel', default=os.path.join(ROOT, 'reference_panel', 'oral_function_panel.json'))
a = ap.parse_args()

panel = json.load(open(a.panel))
qc = os.path.join(a.pub, 'qc_result.json')
rl = json.load(open(qc)).get('read_length_bp') if os.path.exists(qc) else None
if not a.force and (rl is None or int(rl) in set(panel['meta'].get('read_lengths_bp') or [])):
    print(f'old process (read length {rl}bp) - no adjustment'); sys.exit(0)

ref = json.load(open(a.offset))
LO, HI = ref['bands']['low_max_percentile'], ref['bands']['high_min_percentile']
# With few reference dogs the shift itself is uncertain, so a high or low call needs the unadjusted comparison to
# agree; otherwise the result reads typical (2026-11-12: a single swab moved nitrate 14th -> 87th on a 2-dog shift).
MIN_REF_DOGS = 10
me = a.sample.lower()
own_dogs = {r['dog'] for r in ref['samples'] if me == r['sample'] or me in r['aliases']}
pool = [r for r in ref['samples'] if r['dog'] not in own_dogs]
dogs = sorted({r['dog'] for r in pool})
NOTE = ('This sample was sequenced with our newer lab process, which reads mouth bacteria a little differently '
        'from the process behind most of our reference dogs. Comparisons are adjusted for that and shown as '
        'low, typical or high rather than exact percentiles.')

def pct(ref_sorted, v):
    lo, hi = bisect.bisect_left(ref_sorted, v), bisect.bisect_right(ref_sorted, v)
    return 100.0 * (lo + 0.5 * (hi - lo)) / len(ref_sorted)

def band(p): return 'low' if p <= LO else 'high' if p >= HI else 'typical'

def agreed(p_adj, p_raw):
    b = band(p_adj)
    if len(dogs) >= MIN_REF_DOGS or p_raw is None: return b
    return b if band(p_raw) == b else 'typical'

def shift(values_by_dog, cohort_t):
    """Mean of dog means (transformed) minus the cohort median (transformed); None without reference dogs."""
    means = [sum(v) / len(v) for v in values_by_dog.values() if v]
    if not means: return None
    c = sorted(cohort_t)
    return sum(means) / len(means) - c[len(c) // 2]

def by_dog(get):
    out = {}
    for r in pool:
        v = get(r)
        if v is not None: out.setdefault(r['dog'], []).append(v)
    return out

# Gum-disease bacteria share
tp = lambda x: math.log(x + 0.5)
hp = os.path.join(a.pub, 'microbiome_health_result.json')
if os.path.exists(hp):
    h = json.load(open(hp))
    raw_pct = h.get('pathobiont_percentile_unadjusted', h.get('pathobiont_percentile'))
    cohort = [d['pathobiont_pct'] for d in panel['dogs']]
    s = shift(by_dog(lambda r: tp(r['pathobiont_pct'])), [tp(x) for x in cohort])
    h['pathobiont_percentile_unadjusted'] = raw_pct
    for k in ('richness_percentile', 'shannon_percentile'):
        h[k + '_unadjusted'] = h.get(k + '_unadjusted', h.get(k)); h[k] = None
    if s is None:
        h['pathobiont_percentile'] = None; h['pathobiont_band'] = None
    else:
        adj = math.exp(tp(h['pathobiont_burden_pct']) - s) - 0.5
        p = round(pct(sorted(cohort), adj), 1)
        h.update(pathobiont_percentile=p, pathobiont_band=agreed(p, raw_pct), pathobiont_adjusted_pct=round(max(adj, 0.0), 3))
    h.update(platform_adjusted=True, platform_note=NOTE,
             platform_reference=f'{len(dogs)} new-process reference dogs ({len(pool)} swabs), shift = '
                                f'{"x%.2f" % math.exp(s) if s is not None else "n/a"} on gum-disease bacteria share')
    json.dump(h, open(hp, 'w'), indent=2)
    print(f"pathobiont {h['pathobiont_burden_pct']:.1f}%: pct {raw_pct} -> {h['pathobiont_percentile']} "
          f"({h['pathobiont_band']}); reference dogs {dogs}")

# Bacterial activity (HUMAnN modules)
op = os.path.join(a.pub, 'oral_function_result.json')
if os.path.exists(op):
    o = json.load(open(op)); oral = json.load(open(a.oral_panel))
    if o.get('status') == 'ok':
        tc = lambda x: math.log(x + 1.0)
        pm = {m['key']: m for m in oral['modules']}
        for m in o['modules']:
            if not m.get('measurable') and m.get('percentile_unadjusted') is None: continue
            m['measurable'] = True; m.pop('message', None)
            refm = pm[m['key']]
            s = shift(by_dog(lambda r: tc(r['oral_cpm'][m['key']]) if m['key'] in r.get('oral_cpm', {}) else None),
                      [tc(x) for x in refm['ref_cpm_sorted']])
            m['percentile_unadjusted'] = m.get('percentile_unadjusted', m.get('percentile'))
            if s is None:  # no reference swab measured this module: show it, but not as a comparison
                m.update(measurable=False, message='Not compared: too few reference samples from our newer lab process yet.')
                for k in ('percentile', 'level', 'tone'): m.pop(k, None)
                continue
            adj = math.exp(tc(m['cpm']) - s) - 1.0
            p = int(round(pct(refm['ref_cpm_sorted'], round(adj, 3))))
            b = agreed(p, m['percentile_unadjusted'])
            hi, lo = b == 'high', b == 'low'
            if m['direction'] == 'neutral' or not (hi or lo): level, tone = 'typical', 'neutral'
            else:
                level = 'high' if hi else 'low'
                tone = 'concern' if (m['direction'] == 'bad' and hi) or (m['direction'] == 'good' and lo) else 'good'
            m.update(percentile=p, level=level, tone=tone, cpm_adjusted=round(max(adj, 0.0), 2))
        concern = [m for m in o['modules'] if m.get('tone') == 'concern']
        if any(m['key'] in ('sulfur', 'ging') for m in concern):
            o['note'] = ("Worth raising at your next vet visit: the bacteria behind gum disease are more active than in most dogs. "
                         "Ask about a dental check. At home, daily tooth brushing helps most, and dental chews or a dental diet with the "
                         "VOHC seal (vohc.org) are shown to reduce plaque and tartar.")
            o['note_kind'] = 'vet'
        elif concern:
            o['note'] = ("Fewer of the bacteria linked to healthier mouths than in most dogs. Daily tooth brushing and VOHC-approved "
                         "dental chews (vohc.org) help keep the balance healthy.")
            o['note_kind'] = 'watch'
        else:
            o['note'] = "The mouth bacteria look typical. Keep up what you're doing; daily brushing and VOHC-approved chews help keep it that way."
            o['note_kind'] = 'typical'
        o.update(platform_adjusted=True, platform_note=NOTE)
        json.dump(o, open(op, 'w'), indent=1)
        print('oral function: ' + ' '.join(f"{m['key']}={m.get('percentile_unadjusted')}->{m.get('percentile', 'NA')}({m.get('level', '-')})"
                                            for m in o['modules']))
