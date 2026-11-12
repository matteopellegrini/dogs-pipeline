#!/usr/bin/env python3
"""Remove the Dog Aging Project name from customer-visible report text (Matteo 2026-11-11: it must not be
mentioned on reports or the website). Only the full name: DAP is also a gene symbol. Rewrites string values in report JSONs in place; keys are left alone.
usage: scrub_dap.py [--dry] FILE [...]   prints FILE<TAB>n_replacements for each changed file."""
import json, re, sys

RULES = [
    (r'Dog Aging Project HLES owner survey', 'US pet-dog owner survey'),
    (r'Dog Aging Project owner-reported weights \(HLES dd_weight_lbs\)', 'owner-reported weights of US pet dogs'),
    (r'adult dogs in the Dog Aging Project with owner-reported weights', 'adult US pet dogs with owner-reported weights'),
    (r'Dog Aging Project adults', 'US pet-dog adults'),
    (r'Dog Aging Project reference dogs', 'US pet-dog research cohort'),
    (r'Dog Aging Project dogs', 'US pet dogs'),
    (r'the Dog Aging Project', 'a US pet-dog research cohort'),
    (r'Dog Aging Project', 'US pet-dog research cohort'),
]
RX = [(re.compile(a, re.I), b) for a, b in RULES]

def scrub(o, n):
    if isinstance(o, str):
        for rx, b in RX:
            o, k = rx.subn(b, o); n[0] += k
        return o
    if isinstance(o, list):
        return [scrub(x, n) for x in o]
    if isinstance(o, dict):
        return {k: scrub(v, n) for k, v in o.items()}
    return o

if __name__ == '__main__':
    dry = '--dry' in sys.argv
    for f in (a for a in sys.argv[1:] if a != '--dry'):
        with open(f) as fh:
            raw = fh.read()
        if not re.search(r'(?i)dog aging project', raw):
            continue
        n = [0]
        d = scrub(json.loads(raw), n)
        if n[0] and not dry:
            compact = '\n' not in raw.strip()
            with open(f + '.tmp', 'w') as fh:
                json.dump(d, fh, ensure_ascii='\\u' in raw, **({'separators': (',', ':')} if compact else {'indent': 2}))
            import os; os.replace(f + '.tmp', f)
        print(f'{f}\t{n[0]}')
