#!/usr/bin/env python3
"""
Compare two runs of the pipeline on the same samples, JSON-aware.

    python3 analysis/compare_runs.py MAC_DIR HOFFMAN_DIR

    e.g. python3 analysis/compare_runs.py \
             dogs-app/public/ucla-6198  /tmp/results-ucla-check/ucla-6198

Byte-diffing JSON reports is the wrong test: run_date differs by design, float
formatting can differ across numpy/sklearn builds, and key order is arbitrary.
This walks both structures and reports only differences that matter, with
numeric tolerance — the question is whether the two machines computed the same
biology, not whether they serialised it identically.
"""
import json
import math
import sys

REL_TOL = 1e-4          # floats within 0.01% are the same measurement
IGNORE_KEYS = {'run_date', 'built'}   # differ by design


def walk(a, b, path, diffs):
    if type(a) is not type(b) and not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
        diffs.append(f"{path}: type {type(a).__name__} vs {type(b).__name__}")
        return
    if isinstance(a, dict):
        for k in sorted(set(a) | set(b)):
            if k in IGNORE_KEYS:
                continue
            if k not in a:
                diffs.append(f"{path}.{k}: only in second")
            elif k not in b:
                diffs.append(f"{path}.{k}: only in first")
            else:
                walk(a[k], b[k], f"{path}.{k}", diffs)
    elif isinstance(a, list):
        if len(a) != len(b):
            diffs.append(f"{path}: length {len(a)} vs {len(b)}")
            return
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f"{path}[{i}]", diffs)
    elif isinstance(a, float) or isinstance(b, float):
        if not math.isclose(float(a), float(b), rel_tol=REL_TOL, abs_tol=1e-6):
            diffs.append(f"{path}: {a} vs {b}")
    elif a != b:
        diffs.append(f"{path}: {a!r} vs {b!r}")


def main():
    import os
    d1, d2 = sys.argv[1], sys.argv[2]
    files = sorted(f for f in os.listdir(d1) if f.endswith('.json'))
    total = 0
    for f in files:
        p2 = os.path.join(d2, f)
        if not os.path.exists(p2):
            print(f"  {f}: MISSING in {d2}")
            total += 1
            continue
        try:
            a = json.load(open(os.path.join(d1, f)))
            b = json.load(open(p2))
        except Exception as e:
            print(f"  {f}: unreadable ({e})")
            total += 1
            continue
        diffs = []
        walk(a, b, f, diffs)
        if diffs:
            total += len(diffs)
            print(f"  {f}: {len(diffs)} difference(s)")
            for d in diffs[:8]:
                print(f"     {d}")
            if len(diffs) > 8:
                print(f"     ... and {len(diffs) - 8} more")
        else:
            print(f"  {f}: identical (within tolerance)")
    print(f"\n{'REPRODUCED' if total == 0 else f'{total} differences'} across {len(files)} files")


if __name__ == '__main__':
    main()
