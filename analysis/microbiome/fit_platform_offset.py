#!/usr/bin/env python3
"""Collect the new-process microbiome reference (platform_ref.sbatch outputs) into
reference_json/microbiome_platform_offset.json.

The old-lab cohort (reference_panel/microbiome_panel.json, oral_function_panel.json) was sequenced on
Illumina; swabs from the new lab (MGI/DNBSEQ, Complete Genomics) differ systematically in composition.
Too few new-process dogs exist to build a reference of their own, so a new-process sample is scored by
removing the average new-process shift (learned from these reference swabs, depth-matched to the cohort)
and placing the result in the cohort's spread. This script only records the per-swab metrics; the shift
is computed at scoring time by platform_adjust.py, leaving out the dog being scored.

usage: fit_platform_offset.py SAMPLES_TSV SCRATCH_DIR > reference_json/microbiome_platform_offset.json
"""
import json, math, sys, datetime, os

PATHOBIONTS = ['s__Porphyromonas_gulae', 's__Tannerella_forsythia', 's__Porphyromonas_cangingivalis',
               's__Porphyromonas_canoris', 's__Porphyromonas_gingivicanis', 's__Treponema_denticola',
               's__Fusobacterium_nucleatum', 's__Prevotella_intermedia']

def profile(path):
    sp, bact, reads = {}, 100.0, None
    for line in open(path):
        if line.startswith('#'):
            if 'reads processed' in line: reads = int(''.join(c for c in line.split('reads')[0] if c.isdigit()))
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 2: continue
        pct = float(p[2]) if len(p) >= 3 and p[1].replace('|', '').isdigit() else float(p[1])
        last = p[0].split('|')[-1]
        if p[0] == 'k__Bacteria': bact = pct
        if last.startswith('s__'): sp[p[0]] = pct
    return sp, bact, reads

samples = []
for line in open(sys.argv[1]):
    if line.startswith('#') or not line.strip(): continue
    s, dog, process, aliases = line.rstrip('\n').split('\t')
    d = os.path.join(sys.argv[2], s.lower())
    mp = os.path.join(d, f'{s.lower()}_metaphlan.txt')
    if not os.path.exists(mp):
        print(f'missing {mp}', file=sys.stderr); continue
    sp, bact, reads = profile(mp)
    path = sum(v for c, v in sp.items() if any(c.endswith('|' + x) or ('|' + x + '|') in c for x in PATHOBIONTS) and '|t__' not in c)
    ab = [v for v in sp.values() if v > 0]; tot = sum(ab)
    shannon = -sum(v / tot * math.log(v / tot) for v in ab) if tot else 0.0
    rec = {'sample': s.lower(), 'dog': dog, 'process': process, 'aliases': aliases.split(','),
           'microbial_reads': reads, 'total_classified_pct': round(bact, 3),
           'pathobiont_pct': round(path, 4), 'richness': len(ab), 'shannon': round(shannon, 4)}
    of = os.path.join(d, 'oral_function_result.json')
    if os.path.exists(of):
        o = json.load(open(of))
        rec['oral_cpm'] = {m['key']: m['cpm'] for m in o.get('modules', []) if m.get('measurable')}
    samples.append(rec)

json.dump({'version': 1, 'built': datetime.date.today().isoformat(),
           'description': 'New-process (MGI/Complete Genomics) reference swabs, bacterial reads capped at 2M to match '
                          'old-lab cohort depth. Shift = mean over reference dogs (excluding the dog scored) of the '
                          'dog-mean transformed value, minus the cohort median.',
           'transforms': {'pathobiont': 'log(pct + 0.5), pct of all reads as the pipeline and cohort panel compute it', 'oral_cpm': 'log(cpm + 1)'},
           'bands': {'low_max_percentile': 25, 'high_min_percentile': 75},
           'samples': samples}, sys.stdout, indent=1)
print(f'{len(samples)} reference swabs from {len({r["dog"] for r in samples})} dogs', file=sys.stderr)
