#!/usr/bin/env python3
# Two-step ancestry (user design, 2026-09-10): the global lasso proportions
# become FLARE's per-sample ancestry prior (gt-ancestries), the local step
# refines them, and FLARE's .global.anc.gz (genome-wide mean of the local
# posteriors) is the final breed prediction.
#
#   python3 lai_priors.py <breed_result.json> <query sample id> <out.tsv>
#
# Format (FLARE global-ancestry file): header "SAMPLE <anc names...>", then
# one line per sample with proportions. Ancestry names = the reference
# panel names in ref_panel_map.tsv. Breeds the lasso did not use get a
# small floor so the local step can still discover them.
import json, sys
br, sample, out = sys.argv[1:4]
FLOOR = 0.002
ref = sorted({l.split('\t')[1].strip() for l in open('ref_panel_map.tsv') if l.strip()})
b = json.load(open(br))
prior = {}
for e in b.get('breed_composition', []):
    comps = e.get('components') or [{'code': e.get('breed'), 'proportion': e.get('proportion', 0)}]
    for c in comps:
        code = c.get('code') or e.get('breed')
        prior[code] = prior.get(code, 0.0) + float(c.get('proportion', 0))
used = {k: v for k, v in prior.items() if k in ref}
missing = {k: v for k, v in prior.items() if k not in ref and v >= 0.01}
vec = {a: max(used.get(a, 0.0), FLOOR) for a in ref}
tot = sum(vec.values())
vec = {a: v / tot for a, v in vec.items()}
with open(out, 'w') as f:
    f.write('SAMPLE ' + ' '.join(ref) + '\n')
    f.write(sample + ' ' + ' '.join(f'{vec[a]:.6f}' for a in ref) + '\n')
top = sorted(vec.items(), key=lambda kv: -kv[1])[:5]
print(f'{sample}: prior mass on lasso breeds {sum(used.values()):.2f}; lasso breeds absent from reference (>=1%): {missing}; top: {[(a, round(v,3)) for a, v in top]}')
