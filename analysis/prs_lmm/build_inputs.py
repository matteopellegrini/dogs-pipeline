#!/usr/bin/env python3
"""
Build GEMMA inputs for the PRS-LMM experiment.

    python3 build_inputs.py WORKDIR BFILE_PREFIX PIPELINE_SH

Produces, in WORKDIR:
    pheno.tsv    one column per (trait x fold) — fold columns have the held-out
                 breeds' dogs set to NA, so GEMMA's own missing-phenotype
                 handling does the masking and every run reads the same file
    covar.txt    intercept + genotype-source indicator (Parker array = 0,
                 Dog10K WGS = 1). The two cohorts were genotyped by different
                 technologies, and without this the LMM can absorb
                 platform-correlated allele-frequency error into trait betas.
    columns.tsv  trait, fold, 1-based column index for GEMMA -n
    folds.tsv    breed -> fold assignment

Design notes
------------
Phenotypes are breed-level (AKC trait scores and breed-standard heights and
weights) assigned to every dog of the breed. That makes dogs pseudo-replicates:
the honest cross-validation unit is the BREED, so folds split breeds, never
dogs. A dog-level split would leak the training breeds' scores into the test
set and reward memorising breed membership — the exact failure this experiment
exists to measure.

The Parker-code -> AKC-name map, the height table and the weight table are
parsed out of run_dog_pipeline.sh rather than copied, so this harness cannot
drift from what Stage 11 ships.
"""
import csv
import hashlib
import io
import os
import re
import subprocess
import sys

N_FOLDS = 5
AKC_URL = 'https://raw.githubusercontent.com/kkakey/dog_traits_AKC/main/data/breed_traits.csv'

TRAIT_COLS = [
    'Affectionate With Family', 'Good With Young Children', 'Good With Other Dogs',
    'Shedding Level', 'Coat Grooming Frequency', 'Drooling Level',
    'Openness To Strangers', 'Playfulness Level', 'Watchdog/Protective Nature',
    'Adaptability Level', 'Trainability Level', 'Energy Level',
    'Barking Level', 'Mental Stimulation Needs',
]


def extract_dict(src, name):
    """Pull an inline `NAME = { ... }` python dict out of the pipeline script."""
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


def main():
    workdir, bfile, pipeline_sh = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(workdir, exist_ok=True)
    src = open(pipeline_sh).read()
    parker_to_akc = extract_dict(src, 'PARKER_TO_AKC')
    height_cm = extract_dict(src, 'BREED_HEIGHT_CM')
    weight_kg = extract_dict(src, 'BREED_WEIGHT_KG')

    # final6.fam carries the MERGED panel's harmonized long names
    # (ENGLISH_SPRINGER_SPANIEL), not the short Parker codes PARKER_TO_AKC is
    # keyed by, so join on normalised token sets instead: AKC's
    # "Spaniels (English Springer)" and the panel's ENGLISH_SPRINGER_SPANIEL
    # both reduce to {english, springer, spaniel} once plurals are stripped.
    def tokens(name):
        return frozenset(w.rstrip('s') for w in re.sub(r'[^a-z ]', ' ', name.lower()).split() if w)

    cache = os.path.join(workdir, 'breed_traits.csv')
    if not os.path.exists(cache):
        r = subprocess.run(['curl', '-sL', AKC_URL], capture_output=True, text=True, check=True)
        open(cache, 'w').write(r.stdout)
    akc = {row['Breed'].strip().replace('\xa0', ' '): row
           for row in csv.DictReader(open(cache))}

    fam = [l.split() for l in open(bfile + '.fam')]
    breeds = [f[0] for f in fam]
    iids = [f[1] for f in fam]
    # Genotype source by ID convention: Parker IDs carry an underscore
    # (ESSP_04198), Dog10K renames are CODE + 6 digits (ESSP000001). Verified a
    # clean 1306/1589 split with zero ambiguous on final6.
    source = [0 if '_' in i else 1 for i in iids]

    akc_by_tokens = {tokens(n): n for n in akc}
    hw_by_tokens  = {tokens(n): n for n in height_cm}

    def akc_name_for(fid):
        t = tokens(fid.replace('_', ' '))
        return akc_by_tokens.get(t)

    def hw_name_for(fid):
        t = tokens(fid.replace('_', ' '))
        return hw_by_tokens.get(t)

    # per-breed truth
    def breed_value(code, trait):
        if trait in ('height_cm', 'weight_kg'):
            name = hw_name_for(code)
        else:
            name = akc_name_for(code)
        if not name:
            return None
        if trait == 'height_cm':
            return height_cm.get(name)
        if trait == 'weight_kg':
            return weight_kg.get(name)
        row = akc.get(name)
        if not row:
            return None
        v = (row.get(trait) or '').strip()
        try:
            return float(v)
        except ValueError:
            return None

    all_traits = TRAIT_COLS + ['height_cm', 'weight_kg']

    # deterministic breed folds, only over breeds that have any phenotype
    scored = sorted({b for b in set(breeds)
                     if any(breed_value(b, t) is not None for t in all_traits)})
    fold_of = {b: int(hashlib.sha1(b.encode()).hexdigest(), 16) % N_FOLDS + 1
               for b in scored}
    with open(os.path.join(workdir, 'folds.tsv'), 'w') as fh:
        for b in scored:
            fh.write(f'{b}\t{fold_of[b]}\n')

    cols, header = [], []
    for t in all_traits:
        vals = [breed_value(b, t) for b in breeds]
        n_ok = sum(v is not None for v in vals)
        if n_ok < 500:
            print(f'  skip {t}: only {n_ok} dogs with truth')
            continue
        for fold in ['full'] + [str(k) for k in range(1, N_FOLDS + 1)]:
            col = []
            for b, v in zip(breeds, vals):
                masked = v is None or (fold != 'full' and fold_of.get(b) == int(fold))
                col.append('NA' if masked else f'{v}')
            cols.append(col)
            header.append((t, fold))

    with open(os.path.join(workdir, 'pheno.tsv'), 'w') as fh:
        for row in zip(*cols):
            fh.write('\t'.join(row) + '\n')
    with open(os.path.join(workdir, 'columns.tsv'), 'w') as fh:
        for i, (t, fold) in enumerate(header, 1):
            fh.write(f'{t}\t{fold}\t{i}\n')
    with open(os.path.join(workdir, 'covar.txt'), 'w') as fh:
        for s in source:
            fh.write(f'1\t{s}\n')

    n_traits = len({t for t, _ in header})
    print(f'{len(fam)} dogs, {len(scored)} scored breeds, '
          f'{n_traits} traits x {N_FOLDS + 1} folds = {len(header)} phenotype columns')


if __name__ == '__main__':
    main()
