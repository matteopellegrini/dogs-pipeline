#!/usr/bin/env python3
# FLARE reference panel map: sample<TAB>breed for every Parker + Dog10K dog
# whose harmonised canonical breed has >= 8 individuals across both sources.
# Codes -> canonical names come from analysis/breed_accuracy/harmonize.py
# (labels.json). Run from $D/lai with envs/genomics python.
import json, re
from collections import Counter
L = json.load(open('labels.json'))
rows = []
for s in open('../relatives/panel_samples.txt'):
    s = s.strip()
    if not s:
        continue
    m = re.match(r'^([A-Za-z\-]+?)\d+$', s)
    code = m.group(1) if m else s
    rows.append((s, 'dog10k', code, L.get(code)))
for l in open('../COSMO/analysis/cosmo_parker_full.fam'):
    fid, iid = l.split()[:2]
    rows.append((iid, 'parker', fid, L.get(fid)))
cnt = Counter(b for _, _, _, b in rows if b)
keep = {b for b, n in cnt.items() if n >= 8}
with open('ref_map_all.tsv', 'w') as f:
    for s, src, code, b in rows:
        f.write(f"{s}\t{src}\t{code}\t{b or ''}\n")
with open('ref_panel_map.tsv', 'w') as f:
    for s, src, code, b in rows:
        if b in keep:
            f.write(f"{s}\t{b}\n")
with open('ref_keep_samples.txt', 'w') as f:
    for s, src, code, b in rows:
        if b in keep:
            f.write(s + '\n')
unl = [r for r in rows if not r[3]]
print('samples:', len(rows), '| unlabeled:', len(unl), Counter(r[2] for r in unl).most_common(5))
print('breeds with >=8 dogs:', len(keep), '| dogs in those breeds:', sum(1 for r in rows if r[3] in keep))
print('smallest kept:', sorted(((n, b) for b, n in cnt.items() if n >= 8))[:3], '| largest:', cnt.most_common(3))
