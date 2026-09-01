#!/usr/bin/env python3
"""Microbiome panel v2: the 96 study dogs + 1,167 ProsperKits dogs, pooled.

Age labels: 643 dogs (96 cohort + 547 prosper with customer-reported ages at
collection). Validation (2026-09-01, clade-keyed features, 5-fold CV):
  ElasticNet(log10, prevalence>=10%): MAE 2.12y, r 0.67  (old 96-dog Ridge:
  held-out MAE 2.54y, r 0.55 on 555 prosper dogs).
Read length is NOT a batch effect within Illumina (train-151 -> test-101:
MAE 2.06, bias +0.07); the earlier "cross-platform failure" was a feature-key
mismatch in the test harness. The MGI/DNBSEQ gate remains (chemistry shift).

Usage: python3 analysis/microbiome/build_panel_v2.py
Reads  reference_panel/microbiome_panel.json (v1, kept as .v1 backup)
       analysis/microbiome_training_table_clade.json.gz
Writes reference_panel/microbiome_panel.json
"""
import json, gzip, math, os, shutil, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PANEL = os.path.join(ROOT, 'reference_panel', 'microbiome_panel.json')
TABLE = os.path.join(ROOT, 'analysis', 'microbiome_training_table_clade.json.gz')

v1 = json.load(open(PANEL))
if not os.path.exists(PANEL + '.v1'):
    shutil.copy(PANEL, PANEL + '.v1')
pathobionts = v1['meta']['pathobionts']

def stats(abund):
    vals = [v for v in abund.values() if v and v > 0]
    total = sum(vals)
    rich = len(vals)
    shannon = 0.0
    if total > 0:
        for v in vals:
            p = v / total
            shannon -= p * math.log(p)
    patho = sum(v for k, v in abund.items()
                if v and any(k.endswith('|' + p) or k.endswith(p) for p in pathobionts))
    return rich, round(shannon, 4), round(patho, 3)

dogs = []
for d in v1['dogs']:
    d = dict(d)
    d.setdefault('read_length_bp', 151)
    d['source'] = 'cohort96'
    dogs.append(d)

recs = json.load(gzip.open(TABLE, 'rt'))
for r in recs:
    if r['n_species'] < 1:
        continue
    rich, shannon, patho = stats(r['abund'])
    dogs.append({
        'sample': r['sample'], 'age': r['age'], 'species': r['abund'],
        'richness': rich, 'shannon': shannon, 'pathobiont_pct': patho,
        'read_length_bp': r.get('read_length_bp'), 'source': 'prosperkits',
    })

meta = dict(v1['meta'])
meta.update({
    'n_samples': len(dogs),
    'n_with_age': sum(1 for d in dogs if d.get('age') is not None),
    'built': datetime.date.today().isoformat(),
    'read_lengths_bp': sorted({int(d['read_length_bp']) for d in dogs if d.get('read_length_bp')}),
    'ages_source': 'sample_sheet.gen.tsv (96) + sample_sheet.prosper.tsv / prosperKitAgeInfo (547)',
    'version': 2,
    'note': ('Pooled panel: read length is not a batch effect within Illumina '
             '(validated 151->101); the read-length gate now guards against '
             'other chemistries (MGI/DNBSEQ) whose lengths are absent here.'),
})
json.dump({'meta': meta, 'dogs': dogs}, open(PANEL, 'w'))
print(f"panel v2: {len(dogs)} dogs, {meta['n_with_age']} aged, read lengths {meta['read_lengths_bp']}")
print(f"size: {os.path.getsize(PANEL)/1e6:.1f} MB")
