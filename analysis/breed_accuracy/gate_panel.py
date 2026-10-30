#!/usr/bin/env python3
"""
Customer-kit gate for a candidate breed panel. For each customer kit (the pipeline's relatives_ds.tsv: CHROM POS REF
ALT DS) and each panel, run Stage 9 as the report does — the lasso picks the candidate breeds and prior, then the
local-ancestry HMM (analysis/local_ancestry/lai_hmm.py, stage 9c) produces the percentages the customer sees — and
count, per population, the kits where it reaches >= 5%. A candidate fails when a population multiplies against
the live panel. Lasso-only scoring passed v2.1 (Dutch Shepherd 2 -> 4 of 300 kits) while the reports, i.e. the HMM,
took Dutch Shepherd from 9 to 39 of 1,877 kits and Sloughi from 0 to 13.

    python3 gate_panel.py REL_DS_DIR LIVE_PANEL CANDIDATE [CANDIDATE ...] [--jobs 8]
"""
import glob, gzip, os, sys
from collections import Counter
from multiprocessing import Pool
import numpy as np
from scipy.optimize import nnls
from scipy.linalg import solve_triangular
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'local_ancestry'))
import lai_hmm

args = sys.argv[1:]
jobs = int(args[args.index('--jobs') + 1]) if '--jobs' in args else 8
args = [a for i, a in enumerate(args) if a != '--jobs' and (i == 0 or args[i - 1] != '--jobs')]
ds_dir, panels = args[0], args[1:]
paths = sorted(glob.glob(f'{ds_dir}/*.tsv.gz'))
PANEL = {}

def lasso(P, x):
    valid = ~np.isnan(x); Pv = P[valid]; xv = x[valid].astype(np.float32)
    A = np.vstack([Pv, np.ones((1, Pv.shape[1]), dtype=np.float32)]); b = np.hstack([xv, [1.0]]).astype(np.float32)
    G = (A.T @ A).astype(np.float64); rhs = (A.T @ b).astype(np.float64) - 0.5 * 0.3 * float(np.mean(np.diag(Pv.T @ Pv)))
    G[np.diag_indices_from(G)] += 1e-6; R = np.linalg.cholesky(G).T
    q, _ = nnls(R, solve_triangular(R.T, rhs, lower=True)); return q / q.sum()

def one(job):
    pd, path = job
    if pd not in PANEL:
        P, breeds, sites = lai_hmm.load_panel(pd); chroms, pos = lai_hmm.site_coords(pd); PANEL[pd] = (P, breeds, sites, chroms, pos)
    P, breeds, sites, chroms, pos = PANEL[pd]
    x = np.full(len(sites), np.nan)
    for l in gzip.open(path, 'rt'):
        p = l.rstrip('\n').split('\t')
        if len(p) < 5: continue
        hit = sites.get((p[0], p[1]))
        if not hit: continue
        try: d = float(p[4])
        except ValueError: continue
        i, a1, a2 = hit
        if p[2] == a1 and p[3] == a2: x[i] = 2 - d
        elif p[2] == a2 and p[3] == a1: x[i] = d
    q = lasso(P, x)
    breed = {'breed_composition_raw': [{'breed': b, 'breed_name': b, 'proportion': float(v)} for b, v in zip(breeds, q) if v >= 0.005],
             'breed_composition': []}
    keep, pretty, group = lai_hmm.lasso_prior(breed, breeds, 0.02)
    cands = list(keep); idx = [breeds.index(b) for b in cands]; pi = np.array([keep[b] for b in cands])
    dos = lai_hmm.load_query_ds(path, sites); F = P[:, idx].astype(np.float64)
    best = None
    for g in (1, 2, 3, 5, 8, 12):
        ll, gamma = lai_hmm.run_hmm(F, dos, chroms, pos, pi, g, 0.05, 10.0)
        if best is None or ll > best[0]: best = (ll, gamma)
    gamma = best[1]; marg = (gamma.sum(axis=2) + gamma.sum(axis=1)) / 2.0
    return pd, os.path.basename(path), dict(zip(cands, marg.mean(axis=0).tolist()))

if __name__ == '__main__':
    with Pool(jobs) as pool:
        res = pool.map(one, [(pd, p) for pd in panels for p in paths], chunksize=4)
    count = {pd: Counter() for pd in panels}
    for pd, kit, props in res:
        for b, v in props.items():
            if v >= 0.05: count[pd][b] += 1
    live = panels[0]
    print(f'{len(paths)} customer kits, reported (HMM) composition; populations at >=5%: {os.path.basename(live)} -> candidate')
    for pd in panels[1:]:
        rows = sorted(set(count[live]) | set(count[pd]), key=lambda b: -(count[pd][b] - count[live][b]))
        grown = [b for b in rows if count[pd][b] >= max(3, 2 * count[live][b] + 2)]
        print(f'== {os.path.basename(pd)}: {"FAIL" if grown else "PASS"} ({len(grown)} populations multiplied: {", ".join(grown)})')
        for b in rows[:12]: print(f'   {b:34s} {count[live][b]:4d} -> {count[pd][b]:4d}{"   <<" if b in grown else ""}')
