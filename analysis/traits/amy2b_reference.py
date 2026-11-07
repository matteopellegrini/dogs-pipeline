#!/usr/bin/env python3
"""
Cohort reference + validation for AMY2B copy number (analysis/traits/amy2b.py).

  amy2b_reference.py <manifest.tsv> <out_dir_of_kit_jsons> <reference.json> [--report report.tsv]

- One estimate per dog (kits sharing a pub dir counted once), status 'called'.
- Cohort quantiles (0..100) and median -> the app's lib/amy2bReference.json.
- Breed means for dogs >=75% one breed, compared with the ddPCR breed means of
  Arendt et al. 2014 (Anim Genet 45:716, Table 1 + text; diploid copies).
"""
import csv, json, os, sys, argparse
import numpy as np

# Arendt 2014 ddPCR means (diploid). Display order first: well-known breeds spanning the range.
PUBLISHED = [
    ('German Shepherd', 15.71), ('Golden Retriever', 12.91), ('Labrador Retriever', 11.86), ('Poodle', 9.45),
    ('Samoyed', 6.88), ('English Springer Spaniel', 17.28), ('Nova Scotia Duck Tolling Retriever', 14.25),
    ('Rottweiler', 13.66), ('Australian Terrier', 13.68), ('Border Collie', 12.74), ('Norwegian Elkhound', 11.96),
    ('Beagle', 11.88), ('Boxer', 11.55), ('Bearded Collie', 10.05), ('Greenland Sledge Dog', 4.3),
]
SHORT = {'Labrador Retriever': 'Labrador', 'Nova Scotia Duck Tolling Retriever': 'Toller'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('manifest'); ap.add_argument('outdir'); ap.add_argument('reference')
    ap.add_argument('--report')
    a = ap.parse_args()
    seen, dogs = set(), []
    for pk, bc, pub, work in csv.reader(open(a.manifest), delimiter='\t'):
        f = os.path.join(a.outdir, pk + '.json')
        if pub in seen or not os.path.exists(f) or pk.upper().startswith('COSMO3-DS'):
            continue
        seen.add(pub)
        r = json.load(open(f))
        breed, prop = None, 0
        try:
            comp = json.load(open(os.path.join(pub, 'breed_result.json'))).get('breed_composition', [])
            if comp:
                breed, prop = comp[0].get('breed_name', comp[0].get('breed')), comp[0]['proportion']
        except Exception:
            pass
        dogs.append((pk, r, breed, prop))
    called = [d for d in dogs if d[1].get('status') == 'called']
    x = np.array([d[1]['copies'] for d in called])
    q = [round(float(v), 2) for v in np.percentile(x, range(101))]
    width = np.array([d[1]['ci_high'] / d[1]['ci_low'] for d in called])
    print(f'dogs {len(dogs)}  called {len(called)}  median {np.median(x):.2f}  p5 {q[5]}  p95 {q[95]}  '
          f'min {x.min():.2f}  max {x.max():.2f}  precise(<=1.65) {np.mean(width <= 1.65):.1%}', file=sys.stderr)

    by = {}
    for pk, r, b, p in called:
        if b and p >= 0.75:
            by.setdefault(b, []).append(r['copies'])
    ours = {b: (float(np.mean(v)), float(np.std(v)), len(v)) for b, v in by.items()}
    pairs = []
    print('\nbreed                      ours (sd, n)        Arendt 2014', file=sys.stderr)
    for b, m in PUBLISHED:
        o = ours.get(b)
        if o and o[2] >= 2:
            pairs.append((o[0], m))
        print(f'{b:26s} ' + (f'{o[0]:5.1f} ({o[1]:.1f}, {o[2]:3d})' if o else ' ' * 17) + f'   {m}', file=sys.stderr)
    if len(pairs) >= 3:
        o, m = np.array(pairs).T
        print(f'\nbreed-mean r = {np.corrcoef(o, m)[0, 1]:.2f} over {len(pairs)} breeds; mean ours-published = {np.mean(o - m):+.2f}', file=sys.stderr)
    print('\nlargest breed groups in our cohort:', file=sys.stderr)
    for b, (mu, sd, n) in sorted(ours.items(), key=lambda kv: -kv[1][2])[:30]:
        print(f'  {b:30s} {mu:5.1f} sd {sd:4.1f} n {n}', file=sys.stderr)

    ref = {
        'cohort': {'n': len(called), 'median': int(round(float(np.median(x)))), 'quantiles': q},
        'published': [{'breed': SHORT.get(b, b), 'mean': round(m, 1)} for b, m in PUBLISHED],
        'source': 'Cohort: Prosper K9 genomes (depth estimate, analysis/traits/amy2b.py). Breed means: Arendt et al. 2014 Anim Genet 45:716 (ddPCR).',
    }
    json.dump(ref, open(a.reference, 'w'), indent=1)
    if a.report:
        with open(a.report, 'w') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['pk', 'status', 'copies', 'ci_low', 'ci_high', 'unit_reads', 'top_breed', 'top_prop'])
            for pk, r, b, p in dogs:
                w.writerow([pk, r.get('status'), r.get('copies'), r.get('ci_low'), r.get('ci_high'), r.get('unit_reads'), b, round(p, 3)])


if __name__ == '__main__':
    main()
