#!/usr/bin/env python3
"""Sex-chromosome karyotype from coverage: X dosage, the X/Y shared tip (PAR) and Y.

    sex_karyotype.py PUB_DIR [--y y_result.json] [--dry-run]

Reads PUB_DIR/coverage_1mb.json, decides the dog's sex chromosomes and rewrites
its _meta in place:
  predicted_sex    'male' if a Y is present, else 'female' (the dog's sex: weight,
                   lifespan, pronouns, the paternal line)
  x_dosage_sex     what the X copy number alone says (1 X = 'male', 2+ = 'female').
                   Stage 6 normalised chrX and chose its reference panel with this,
                   and X-linked conditions follow it: an XXY dog carries two X.
  sex_chromosomes  karyotype string + the three measurements behind it
  aneuploidy       gains a chrX entry for an XXY or X0 dog (XXX is already called
                   by stage 6, which compares a female's X with two copies)
Prints the karyotype; exit status 0 whether or not anything changed.

Why: stage 6 called sex from X dosage alone (X/auto < 0.75 = male). An XXY dog has
female X dosage, so DOGS-Gen-15, a male chocolate Lab, was reported as a female
with a normal karyotype (2026-11-14). canFam4 chrX starts with the pseudoautosomal
region (PAR, 0 to ~6.6 Mb), shared with the Y: Y-PAR reads land on the X-PAR, so
PAR depth counts X + Y copies while the rest of the X counts X copies alone.

  karyotype    non-PAR X   PAR    Y/auto
  XY           0.5         1.0    ~0.43
  XX           1.0         1.0    0
  XXY          1.0         1.5    ~0.43
  XXX          1.5         1.5    0
  X0           0.5         0.5    0

A female sample contaminated with male DNA, or an XX/XY chimera, never lifts PAR
above 1.0, so PAR separates XXY from contamination. On 1,871 kits (2026-11-14):
PAR 1.016 +/- 0.046 at >=2x; male Y/auto min 0.27, female max 0.03.
"""
import argparse, json, os, statistics, sys

PAR_WINDOWS = range(0, 6)     # 1 Mb windows 0-6 Mb; male depth halves from 7 Mb
NONPAR_FROM = 7
Y_PRESENT = 0.15              # Y/autosome; male min 0.27 / female max 0.03 in the cohort
Y_ABSENT_MIN_DEPTH = 0.5      # below this a real male's Y median can read ~0: never infer X0
PAR_ONLY_EXCESS = 0.35        # PAR - nonPAR needed to call a Y without Y reads (XY ~0.5)
PAR_ONLY_MIN = 1.35           # and PAR itself (sd 0.046 at >=2x; XXY 1.54)


def measure(cov):
    auto = [d for c, v in cov.items() if c != '_meta' and c[3:].isdigit() for d in v['cosmo'] if d > 0]
    if not auto:
        return None, None
    med = statistics.median(auto)
    x = cov['chrX']['cosmo']
    par = [x[i] / med for i in PAR_WINDOWS if i < len(x) and x[i] > 0]
    nonpar = [d / med for d in x[NONPAR_FROM:] if d > 0]
    return (round(statistics.mean(par), 3) if par else None,
            round(statistics.median(nonpar), 3) if nonpar else None)


def classify(par, nonpar, x_copies, confidence, y=None, depth=None):
    """-> (karyotype, predicted_sex, basis). x_copies comes from stage 6 (its X/auto
    call plus any chrX aneuploidy it made, e.g. XXX); this adds only the Y."""
    y_known = (y is not None and y.get('y_to_autosomal') is not None
               and y.get('status') != 'no_chrY_in_reference')
    ya = y.get('y_to_autosomal') if y_known else None
    if y_known and ya >= Y_PRESENT:
        has_y, basis = True, f'Y reads (Y/autosome {ya:.2f})'
    elif y_known and (depth or 0) >= Y_ABSENT_MIN_DEPTH:
        has_y, basis = False, f'no Y reads (Y/autosome {ya:.2f})'
    elif (x_copies >= 2 and par is not None and confidence != 'low'
          and par >= PAR_ONLY_MIN and par - nonpar >= PAR_ONLY_EXCESS):
        has_y, basis = True, f'shared X/Y region depth {par:.2f} vs X {nonpar:.2f}'
    else:
        has_y, basis = x_copies == 1, 'X dosage'
    kary = 'X' * x_copies + ('Y' if has_y else '')
    return ('X0' if kary == 'X' else kary), ('male' if has_y else 'female'), basis


def apply(pub, y_path=None, dry=False):
    path = os.path.join(pub, 'coverage_1mb.json')
    cov = json.load(open(path))
    meta = cov.setdefault('_meta', {})
    x_sex = meta.get('x_dosage_sex') or meta.get('predicted_sex')
    if x_sex not in ('male', 'female') or 'chrX' not in cov:
        return None, False
    par, nonpar = measure(cov)
    if nonpar is None:
        print('sex karyotype: no coverage, left as is')
        return None, False
    y = json.load(open(y_path)) if y_path and os.path.exists(y_path) else None
    depth = None
    try:
        depth = json.load(open(os.path.join(pub, 'qc_result.json'))).get('genome_mean_depth')
    except Exception:
        pass
    aneu = [a for a in meta.get('aneuploidy', []) if a.get('source') != 'sex_karyotype']
    stage6_x = [a for a in aneu if a.get('chrom') == 'chrX']
    x_copies = stage6_x[0]['copies_est'] if stage6_x else (1 if x_sex == 'male' else 2)
    kary, sex, basis = classify(par, nonpar, x_copies, meta.get('coverage_confidence'), y, depth)
    if kary not in ('XX', 'XY') and not stage6_x:
        expected = 1 if sex == 'male' else 2
        if x_copies != expected:
            aneu.append({
                'chrom': 'chrX', 'direction': 'gain' if x_copies > expected else 'loss',
                'call': {'XXY': 'XXY', 'X0': 'monosomy'}.get(kary, kary),
                'karyotype': kary,
                'label': {'XXY': 'XXY (a male with an extra X chromosome)',
                          'X0': 'X0 (a single X chromosome and no Y)'}.get(kary, kary),
                'copies_est': x_copies, 'expected_copies': expected,
                'median_ratio': nonpar,
                'confidence': meta.get('coverage_confidence', 'ok'),
                'source': 'sex_karyotype'})
    new = {'predicted_sex': sex, 'x_dosage_sex': x_sex, 'aneuploidy': aneu,
           'sex_chromosomes': {'karyotype': kary, 'basis': basis,
                               'par_auto_ratio': par,
                               'nonpar_x_ratio': nonpar,
                               'y_to_autosomal': y.get('y_to_autosomal') if y else None}}
    changed = meta.get('aneuploidy', []) != aneu
    sex_changed = meta.get('predicted_sex') != sex
    if not dry:
        meta.update(new)
        tmp = path + '.tmp'
        with open(tmp, 'w') as f:
            json.dump(cov, f)
        os.replace(tmp, path)
    print(f"sex karyotype: {kary} -> {sex} ({basis}; PAR {par}, X {nonpar})"
          + (' [CHANGED predicted_sex]' if sex_changed else ' [aneuploidy updated]' if changed else ''))
    return new, sex_changed


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('pub')
    ap.add_argument('--y', help='y_result.json (Y/autosome ratio); omit when the reference has no chrY')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()
    _, sex_changed = apply(a.pub, a.y, a.dry_run)
    sys.exit(0)
