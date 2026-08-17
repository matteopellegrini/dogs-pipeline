#!/usr/bin/env python3
"""
Validate PRS/trait predictions across a results cohort against each dog's
own breed call — the check the training CV cannot do, because these dogs
are real individuals, not panel rows.

    python3 cluster/check-prs-cohort.py <results-root> [run_dog_pipeline.sh]

For every <results-root>/*/prs_result.json with a sibling breed_result.json:
  - purebred dogs (top-1 breed fraction >= 0.80) with a breed-standard
    height/weight: predicted vs breed-standard value
  - the same dogs, behaviour traits: predicted score vs the breed's AKC truth,
    correlated per trait across the cohort

Mixed dogs are listed with their predictions only — there is no per-dog truth
for them, which is exactly why the purebred subset carries the validation.

Breed tables are parsed out of the pipeline script (same extraction as
analysis/prs_lmm/build_inputs.py), so this cannot drift from what Stage 11
ships. AKC behaviour truth comes from the cached breed_traits.csv next to the
pipeline if present, else is fetched.
"""
import csv
import json
import os
import re
import subprocess
import sys

PUREBRED_MIN = 0.80
AKC_URL = 'https://raw.githubusercontent.com/kkakey/dog_traits_AKC/main/data/breed_traits.csv'

TRAIT_COLS = [
    'Affectionate With Family', 'Good With Young Children', 'Good With Other Dogs',
    'Shedding Level', 'Coat Grooming Frequency', 'Drooling Level',
    'Openness To Strangers', 'Playfulness Level', 'Watchdog/Protective Nature',
    'Adaptability Level', 'Trainability Level', 'Energy Level',
    'Barking Level', 'Mental Stimulation Needs',
]


def extract_dict(src, name):
    m = re.search(rf'^{name} = \{{', src, re.M)
    if not m:
        raise SystemExit(f'{name} not found in pipeline script')
    depth, i = 0, m.end() - 1
    while i < len(src):
        if src[i] == '{':
            depth += 1
        elif src[i] == '}':
            depth -= 1
            if depth == 0:
                break
        i += 1
    ns = {}
    exec(src[m.start():i + 1], {}, ns)
    return ns[name]


def tokens(name):
    return frozenset(w.rstrip('s') for w in re.sub(r'[^a-z ]', ' ', name.lower()).split() if w)


def main():
    root = sys.argv[1]
    pipeline_sh = sys.argv[2] if len(sys.argv) > 2 else \
        os.path.join(os.path.dirname(__file__), '..', 'run_dog_pipeline.sh')
    src = open(pipeline_sh, encoding='utf-8').read()
    height_cm = extract_dict(src, 'BREED_HEIGHT_CM')
    weight_kg = extract_dict(src, 'BREED_WEIGHT_KG')
    parker_to_akc = extract_dict(src, 'PARKER_TO_AKC')
    hw_by_tokens = {tokens(n): n for n in height_cm}

    def breed_longname(breed):
        # older results carry short Parker codes (SPOO); newer ones the merged
        # panel's long names (STANDARD_POODLE) — normalise to an AKC-ish name
        return parker_to_akc.get(breed, breed).replace('_', ' ')

    cache = os.path.join(os.path.dirname(pipeline_sh), 'breed_traits.csv')
    if not os.path.exists(cache):
        r = subprocess.run(['curl', '-sL', AKC_URL], capture_output=True, text=True, check=True)
        open(cache, 'w', encoding='utf-8').write(r.stdout)
    akc = {row['Breed'].strip().replace('\xa0', ' '): row
           for row in csv.DictReader(open(cache, encoding='utf-8'))}
    akc_by_tokens = {tokens(n): n for n in akc}

    rows = []
    for d in sorted(os.listdir(root)):
        pd = os.path.join(root, d)
        try:
            prs = json.load(open(os.path.join(pd, 'prs_result.json'), encoding='utf-8'))
            br = json.load(open(os.path.join(pd, 'breed_result.json'), encoding='utf-8'))
        except Exception:
            continue
        comp = br.get('breed_composition') or []
        top = comp[0] if comp else {'breed': '?', 'proportion': 0.0}
        rows.append((d, top['breed'], float(top['proportion']), prs))
    if not rows:
        raise SystemExit(f'no prs_result.json + breed_result.json pairs under {root}')

    print(f'{len(rows)} dogs under {root}\n')

    # ── height / weight, purebreds ────────────────────────────────────────
    print(f"{'dog':<14} {'top breed':<28} {'frac':>5}  "
          f"{'pred_kg':>7} {'std_kg':>7}  {'pred_cm':>7} {'std_cm':>7}")
    hw_pairs_w, hw_pairs_h = [], []
    n_mixed = 0
    for d, breed, frac, prs in rows:
        phys = prs.get('physical_traits', {})
        pw = phys.get('weight_kg', {}).get('pred_kg')
        ph = phys.get('height_cm', {}).get('pred_cm')
        if frac < PUREBRED_MIN:
            n_mixed += 1
            continue
        name = hw_by_tokens.get(tokens(breed_longname(breed)))
        sw = weight_kg.get(name) if name else None
        sh = height_cm.get(name) if name else None
        if pw is not None and sw is not None:
            hw_pairs_w.append((pw, sw))
        if ph is not None and sh is not None:
            hw_pairs_h.append((ph, sh))
        print(f"{d:<14} {breed[:27]:<28} {frac:>5.2f}  "
              f"{pw if pw is not None else '—':>7} {sw if sw is not None else '—':>7}  "
              f"{ph if ph is not None else '—':>7} {sh if sh is not None else '—':>7}")
    print(f"({n_mixed} dogs below top-1 fraction {PUREBRED_MIN} not shown; "
          f"no per-dog truth for mixes)\n")

    def corr(pairs):
        import statistics as st
        if len(pairs) < 3:
            return float('nan'), float('nan')
        xs, ys = [p[0] for p in pairs], [p[1] for p in pairs]
        mx, my = st.mean(xs), st.mean(ys)
        num = sum((x - mx) * (y - my) for x, y in pairs)
        den = (sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys)) ** 0.5
        mae = st.mean(abs(x - y) for x, y in pairs)
        return (num / den if den else float('nan')), mae

    for label, pairs, unit in (('weight', hw_pairs_w, 'kg'), ('height', hw_pairs_h, 'cm')):
        r, mae = corr(pairs)
        print(f'purebred {label:<6}: n={len(pairs):<3} r={r:+.3f}  MAE={mae:.1f} {unit}')

    # ── behaviour traits, purebreds ───────────────────────────────────────
    print(f"\n{'trait':<28} {'n':>3} {'r':>7} {'MAE':>6}   (predicted_score vs top breed's AKC value)")
    for trait in TRAIT_COLS:
        pairs = []
        for d, breed, frac, prs in rows:
            if frac < PUREBRED_MIN:
                continue
            v = prs.get('traits', {}).get(trait, {}).get('predicted_score')
            name = akc_by_tokens.get(tokens(breed_longname(breed)))
            row = akc.get(name) if name else None
            tv = (row.get(trait) or '').strip() if row else ''
            try:
                tvf = float(tv)
            except ValueError:
                continue
            if v is not None:
                pairs.append((float(v), tvf))
        r, mae = corr(pairs)
        print(f'{trait:<28} {len(pairs):>3} {r:>7.3f} {mae:>6.2f}')


if __name__ == '__main__':
    main()
