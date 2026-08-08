#!/usr/bin/env python3
"""
Print a canonical one-line-per-metric summary of a results directory.

    python3 cluster/summarize-results.py <results-dir>

Run it on two sites and diff the output to confirm they produce the same
science, not merely that both exited cleanly. Values are formatted to fixed
precision so trivial float noise does not show up as a difference.
"""
import json, os, sys

d = sys.argv[1] if len(sys.argv) > 1 else '.'

def load(name):
    try:
        with open(os.path.join(d, name)) as fh:
            return json.load(fh)
    except Exception:
        return None

def emit(label, value):
    print(f"{label:<34} {value}")

def num(x, places=4):
    return f"{x:.{places}f}" if isinstance(x, (int, float)) else str(x)

# ── breed ────────────────────────────────────────────────────
b = load('breed_result.json')
if b and b.get('breed_composition'):
    for row in b['breed_composition'][:5]:
        emit(f"breed.{row['breed']}", num(row['proportion']))
    emit("breed.n_components", len(b['breed_composition']))

# ── known variants ───────────────────────────────────────────
o = load('omia_result.json')
if o and 'summary' in o:
    for k in ('total_screened', 'affected_snv', 'affected_high_confidence',
              'in_dog10k_panel', 'called_from_bam', 'not_callable'):
        if k in o['summary']:
            emit(f"omia.{k}", o['summary'][k])
    emit("omia.mean_depth", num(o['summary'].get('mean_depth'), 2))

# ── inbreeding ───────────────────────────────────────────────
i = load('inbreeding_result.json')
if i:
    for k in ('f_roh', 'roh_total_mb', 'roh_n_segments'):
        if k in i:
            emit(f"inbreeding.{k}", num(i[k]))
i2 = load('inbreeding_froh_dog10k_result.json')
if i2 and 'sample_froh' in i2:
    emit("inbreeding.sample_froh_dog10k", num(i2['sample_froh']))

# ── qc ───────────────────────────────────────────────────────
q = load('qc_result.json')
if q:
    for k in ('genome_mean_depth', 'total_reads_raw', 'total_reads_after_qc'):
        if k in q:
            emit(f"qc.{k}", num(q[k], 2))

# ── microbiome ───────────────────────────────────────────────
m = load('microbiome_result.json')
if m:
    sp = m.get('species')
    if isinstance(sp, list):
        emit("microbiome.n_species", len(sp))
a = load('microbiome_age_result.json')
if a and 'predicted_age_years' in a:
    emit("microbiome.predicted_age_years", num(a['predicted_age_years'], 2))
h = load('microbiome_health_result.json')
if h:
    for k in ('richness_percentile', 'shannon_percentile'):
        if k in h:
            emit(f"microbiome.{k}", num(h[k], 1))

# ── coverage / cnv, as counts only ───────────────────────────
for name, label in (('coverage_1mb.json', 'coverage.windows'),
                    ('cnv_homdel.json', 'cnv.homdel_entries'),
                    ('functional_variants.json', 'functional.entries')):
    v = load(name)
    if isinstance(v, dict):
        emit(label, len(v))
    elif isinstance(v, list):
        emit(label, len(v))

emit("files.json_count", len([f for f in os.listdir(d) if f.endswith('.json')]))
