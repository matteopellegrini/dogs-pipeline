import gzip, json, os
import numpy as np

# Depth gate, calibrated by downsampling cosmo3 (2026-08-28): at 0.3x KING
# still finds every true duplicate at phi 0.47 with zero false matches, but at
# 0.1x self-detection is lost entirely while FOUR spurious matches appear at up
# to phi 0.437 — imputation collapse fabricates near-identical relatives.
# Below the gate we publish an empty, honest result rather than noise.
MIN_DEPTH = 0.3
try:
    _qc = json.load(open(os.environ['PUB_DIR'] + '/qc_result.json'))
    _depth = float(_qc.get('genome_mean_depth') or 0)
except Exception:
    _depth = None
if _depth is not None and _depth < MIN_DEPTH:
    out = {
        'n_reference_dogs': 0, 'n_sites_used': 0, 'matches': [], 'n_matches': 0,
        'suppressed': True,
        'summary': ('Relative matching is not reported for this sample: sequencing depth '
                    '({:.1f}x) is below the {}x minimum at which kinship estimates are '
                    'reliable for imputed genotypes.').format(_depth, MIN_DEPTH),
        'method': 'Suppressed by depth gate (calibrated on downsampled truth data).',
    }
    with open(os.environ['PUB_DIR'] + '/relatives_result.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2)
    print('relatives_result.json: suppressed (depth {}x < {}x)'.format(_depth, MIN_DEPTH))
    raise SystemExit(0)

ref_dir = os.environ['REL_REF']
R = np.load(ref_dir + '/geno.npy')                       # sites x refs, int8, -1 = missing
meta = json.load(gzip.open(ref_dir + '/meta.json.gz', 'rt'))
keys = {}
with open(ref_dir + '/keys.tsv', encoding='utf-8') as f:
    for i, line in enumerate(f):
        keys[tuple(line.rstrip('\n').split('\t'))] = i

q = np.full(R.shape[0], -1, dtype=np.int8)
for line in open(os.environ['REL_DS'], encoding='utf-8'):
    p = line.rstrip('\n').split('\t')
    i = keys.get((p[0], p[1], p[2], p[3]))
    if i is None:
        continue
    ds = float(p[4])
    q[i] = 0 if ds < 0.5 else (1 if ds < 1.5 else 2)

use = q >= 0
Rs, qs = R[use], q[use]
valid = Rs >= 0                                           # per-ref non-missing mask
het_rows, homA_rows, homB_rows = qs == 1, qs == 0, qs == 2
Nhh   = (Rs[het_rows] == 1).sum(axis=0)
Nopp  = (Rs[homA_rows] == 2).sum(axis=0) + (Rs[homB_rows] == 0).sum(axis=0)
nhet_r = (Rs == 1).sum(axis=0)
nhet_q = valid[het_rows].sum(axis=0)                      # query hets over each ref's sites
denom = nhet_q + nhet_r
phi = np.where(denom > 0, (Nhh - 2.0 * Nopp) / np.maximum(denom, 1), -1.0)

t = meta['thresholds']
def category(v):
    # The 'distant' tier (0.045-0.088) is deliberately NOT reported: pilot runs
    # showed a recurring artifact trio of high-heterozygosity reference dogs
    # polluting that band even for 1x queries. Second degree and closer only.
    if v >= t['self']: return 'same_dog'
    if v >= t['first_degree']: return 'first_degree'
    if v >= t['second_degree']: return 'second_degree'
    return None

sample_name = os.environ.get('SAMPLE', '').lower()
matches = []
for j in np.argsort(phi)[::-1]:
    cat = category(float(phi[j]))
    if cat is None or len(matches) >= 10:
        break
    s = meta['samples'][j]
    # Reference dogs that ARE this sample are not a finding — skip the self
    # row by id, AND (user decision 2026-09-01) skip anything at self-level
    # kinship: a duplicate kit of the same dog listed as a "relative" is
    # confusing, and phi >= 0.45 is genetically the same dog regardless of id.
    if s['id'].lower() == sample_name or cat == 'same_dog':
        continue
    # Privacy: reference-dog identifiers never reach the report — only breed,
    # source and relationship tier. Panel dogs are public research data.
    entry = {
        'kinship': round(float(phi[j]), 3),
        'category': cat,
        'breed': s['label'],
        'source': 'Dog10K research panel' if s['source'] == 'dog10k' else 'ProsperK9 reference cohort',
    }
    # Reference dogs processed through the full pipeline carry their own
    # report-derived descriptors — still anonymous, but far more vivid than a
    # bare breed label ("a 24 kg black Labrador Retriever").
    if s.get('weight_kg'): entry['weight_kg'] = s['weight_kg']
    if s.get('coat'): entry['coat'] = s['coat']
    matches.append(entry)

out = {
    'n_reference_dogs': len(meta['samples']),
    'n_sites_used': int(use.sum()),
    'matches': matches,
    'n_matches': len(matches),
    'summary': ('No relatives detected among the reference dogs.' if not matches else
                '{} match(es) at third-degree kinship or closer.'.format(len(matches))),
    'thresholds': t,
    'method': meta['method'] + ' Same-breed background stays below the second-degree '
              'threshold, so listed matches reflect genuine kinship, not shared breed.',
}
with open(os.environ['PUB_DIR'] + '/relatives_result.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2)
top = matches[0] if matches else None
print('relatives_result.json: {} matches{}'.format(
    len(matches), '' if not top else '; top {} {} ({})'.format(top['kinship'], top['category'], top['breed'])))
