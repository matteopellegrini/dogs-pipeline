#!/usr/bin/env python3
"""Side-by-side breed composition of two breed_result.json files (grouped by display name).
    python3 compare_results.py OLD.json NEW.json [label]"""
import json, sys
def comp(p):
    d = json.load(open(p)); out = {}
    for b in d.get('breed_composition') or []:
        name = b.get('breed_name') or b.get('breed'); out[name] = out.get(name, 0) + b.get('proportion', 0)
    return out
c0, c1 = comp(sys.argv[1]), comp(sys.argv[2])
label = sys.argv[3] if len(sys.argv) > 3 else ''
keys = sorted(set(c0) | set(c1), key=lambda k: -max(c0.get(k, 0), c1.get(k, 0)))[:12]
print(f'== {label}: {"breed":34s} {"old":>7s} {"new":>7s}')
for k in keys: print(f'   {k:34s} {100*c0.get(k,0):6.1f}% {100*c1.get(k,0):6.1f}%')
print(f'   components old/new: {len(c0)}/{len(c1)}')
