#!/usr/bin/env python3
"""
Build the oral-microbiome reference panel from the cohort's own pipeline output.

    python3 build_microbiome_panel.py OUT.json 'results/dogs-gen-*/microbiome_result.json' sample_sheet.gen.tsv

Why this replaces the 1,045-sample CSV
--------------------------------------
The previous reference (merged_microbiome_age_weight_3.18_final.csv) was
profiled with MetaPhlAn 3.18 under the old NCBI taxonomy, while samples run
MetaPhlAn4 with the Jan25 GTDB database. The taxonomies disagree enough that
the age model matched only 16 of its 62 species features for a typical sample,
and the alpha-diversity percentiles were meaningless (a dog at "percentile
100.0" against a reference median of 18 genera). A reference built from the
96 cohort dogs processed by THIS pipeline — same MetaPhlAn version, same
database, same read handling — is smaller but actually comparable: every
feature matches by construction.

Ages come from the sample sheet (column 3), where the cohort's owner-reported
ages already live.

The output is deliberately small (per-dog species profiles, not read data), so
unlike the coverage panels it is committed to git.
"""
import glob
import json
import math
import sys
from datetime import date

# Species below this cohort prevalence are dropped from the panel: they cannot
# be model features, and keeping them only bloats the file.
MIN_PREVALENCE = 0.10

# Must stay in lockstep with the PATHOBIONTS table in Stage 16 of
# run_dog_pipeline.sh — the percentile is only meaningful if the panel's totals
# were computed from the same list as the sample's.
PATHOBIONTS = [
    's__Porphyromonas_gulae', 's__Tannerella_forsythia',
    's__Porphyromonas_cangingivalis', 's__Porphyromonas_canoris',
    's__Porphyromonas_gingivicanis', 's__Treponema_denticola',
    's__Fusobacterium_nucleatum', 's__Prevotella_intermedia',
]


def load_ages(sheet_path):
    ages = {}
    with open(sheet_path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 4:
                continue
            # key by output_name (col 4, lower-cased pipeline name)
            try:
                ages[f[3].lower()] = float(f[2])
            except (ValueError, IndexError):
                continue
    return ages


def main():
    out_path, pattern, sheet = sys.argv[1], sys.argv[2], sys.argv[3]
    ages = load_ages(sheet)
    print(f"ages for {len(ages)} samples from {sheet}")

    dogs = []
    db_versions = set()
    read_lengths = set()
    for p in sorted(glob.glob(pattern)):
        d = json.load(open(p))
        # The platform fingerprint. An age model trained on this panel must not
        # be applied to samples from another platform: the four MGI (100bp)
        # samples validated against the Illumina (151bp) panel all came out
        # ~+5 years over baseline with near-identical feature contributions
        # for a 2.5-year-old and a 9-year-old — batch shift read as age. No
        # cheap per-sample statistic detects it (range and distance checks
        # both fail); read length at least catches cross-platform use.
        try:
            import os as _os
            q = json.load(open(_os.path.join(_os.path.dirname(p), 'qc_result.json')))
            if q.get('read_length_bp'):
                read_lengths.add(int(q['read_length_bp']))
        except Exception:
            pass
        sample = d.get('sample', p)
        db_versions.add(d.get('db_version', '?'))
        # species abundances as % of classified bacteria — the same units
        # Stage 16 uses for the sample side. microbiome_result.json stores
        # % of ALL reads (scaled by total_classified), so unscale.
        scale = (d.get('total_classified_pct') or 100.0) / 100.0
        species = {}
        for s in d.get('species', []):
            v = s['relative_abundance'] / scale if scale > 0 else 0.0
            if v > 0:
                species[s['clade']] = round(v, 6)
        patho = 0.0
        for clade, v in species.items():
            if any(px in clade for px in PATHOBIONTS) and '|t__' not in clade:
                patho += v
        vals = list(species.values())
        tot = sum(vals)
        shannon = -sum((v / tot) * math.log(v / tot) for v in vals) if tot > 0 else 0.0
        dogs.append({'sample': sample, 'age': ages.get(str(sample).lower()),
                     'richness': len(species), 'shannon': round(shannon, 4),
                     'pathobiont_pct': round(patho, 4), 'species': species})

    if len(db_versions) > 1:
        raise SystemExit(f"ERROR: mixed MetaPhlAn databases in cohort: {db_versions} — "
                         "a panel across versions rebuilds the exact problem this replaces")

    # prevalence filter across the cohort
    from collections import Counter
    prev = Counter()
    for d in dogs:
        for c in d['species']:
            prev[c] += 1
    keep = {c for c, n in prev.items() if n / len(dogs) > MIN_PREVALENCE}
    for d in dogs:
        d['species'] = {c: v for c, v in d['species'].items() if c in keep}

    n_aged = sum(1 for d in dogs if d['age'] is not None)
    doc = {'meta': {'n_samples': len(dogs), 'n_with_age': n_aged,
                    'db_version': next(iter(db_versions)),
                    'built': str(date.today()),
                    'min_prevalence': MIN_PREVALENCE,
                    'n_species_kept': len(keep),
                    'ages_source': sheet.split('/')[-1],
                    'read_lengths_bp': sorted(read_lengths),
                    'pathobionts': PATHOBIONTS},
           'dogs': dogs}
    with open(out_path, 'w') as fh:
        json.dump(doc, fh, separators=(',', ':'))
    import os
    print(f"panel: {len(dogs)} dogs ({n_aged} with age), {len(keep)} species kept, "
          f"{os.path.getsize(out_path)/1e6:.2f} MB -> {out_path}")


if __name__ == '__main__':
    main()
