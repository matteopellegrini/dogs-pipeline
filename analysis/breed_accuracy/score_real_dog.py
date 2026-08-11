#!/usr/bin/env python3
"""
Score a real dog's imputed BCF against a merged panel, the way Stage 9 would.

    python3 score_real_dog.py PANEL_PREFIX IMPUTED.bcf [IMPUTED.bcf ...]

Reproduces Stage 9's projection exactly: dosage is the GP-weighted expected
count of a1, restricted to the panel's SNPs, projected by NNLS with a sum-to-1
row and then normalised.

These dogs are not in the panel, so no holdout is needed — the test is honest by
construction. Replicate runs of the same dog (kiki/kiki2, ferdl/ferdl2) also
show how stable the call is across independent sequencing and imputation.
"""
import re
import subprocess
import sys
import tempfile
import os

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import nnls


def load_panel(pre):
    breeds = [l.strip() for l in open(pre + '.breeds')]
    bim = [l.split() for l in open(pre + '.bim')]
    n_snp = len(bim)
    bidx = {b: j for j, b in enumerate(breeds)}
    P = np.zeros((n_snp, len(breeds)), dtype=np.float32)
    with open(pre + '.frq.strat') as fh:
        fh.readline()
        i, seen = -1, 0
        for line in fh:
            f = line.split()
            j = bidx.get(f[2])
            if j is None:
                continue
            if seen % len(breeds) == 0:
                i += 1
            P[i, j] = float(f[5])
            seen += 1
    return breeds, bim, P


def dosages_from_bcf(bcf, bim, bed_path):
    """GP-weighted E[count of a1] at each panel SNP; NaN where not called"""
    idx = {}
    for i, b in enumerate(bim):
        c = b[0] if b[0].startswith('chr') else 'chr' + b[0]
        idx[(c, b[3])] = i
    out = np.full(len(bim), np.nan, dtype=np.float32)
    p = subprocess.run(['bcftools', 'query', '-R', bed_path,
                        '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%GP]\n', bcf],
                       capture_output=True, text=True)
    for line in p.stdout.split('\n'):
        if not line:
            continue
        f = line.split('\t')
        if len(f) < 5:
            continue
        i = idx.get((f[0], f[1]))
        if i is None:
            continue
        a1, a2 = bim[i][4], bim[i][5]
        ref, alt = f[2], f[3]
        if ref not in (a1, a2) or alt not in (a1, a2):
            continue
        try:
            gp = [float(x) for x in f[4].split(',')]
        except ValueError:
            continue
        if len(gp) < 3:
            continue
        # GP = [P(hom ref), P(het), P(hom alt)] -> expected copies of a1
        out[i] = (2.0 * gp[0] + gp[1]) if ref == a1 else (gp[1] + 2.0 * gp[2])
    return out


def main():
    pre = sys.argv[1]
    args = sys.argv[2:]
    top_n = int(args[args.index('--top') + 1]) if '--top' in args else 8
    LASSO = float(args[args.index('--lasso') + 1]) if '--lasso' in args else 0.3
    min_q = float(args[args.index('--min') + 1]) if '--min' in args else 0.005
    bcfs = [a for a in args if a.endswith('.bcf')]
    breeds, bim, P = load_panel(pre)
    K = len(breeds)
    exp = (2.0 * P.mean(axis=1)).astype(np.float32)
    print(f"panel: {P.shape[0]} SNPs x {K} breeds\n")

    fd, bed = tempfile.mkstemp(suffix='.bed')
    with os.fdopen(fd, 'w') as fh:
        rows = sorted(((b[0] if b[0].startswith('chr') else 'chr' + b[0]), int(b[3]))
                      for b in bim)
        for c, pos in rows:
            fh.write(f"{c}\t{pos-1}\t{pos}\n")

    for bcf in bcfs:
        x = dosages_from_bcf(bcf, bim, bed)
        cov = np.isfinite(x)
        # RESTRICT to covered sites rather than mean-imputing the rest, which is
        # what Stage 9 does. Imputing with the panel average pulls the fit toward
        # the mean and inflates spurious components — on PK9-0002 (13.7%
        # uncovered) it moved German Shepherd 50.8 -> 45.2% and doubled a
        # spurious village-dog component. Verified to reproduce the pipeline to
        # two decimals.
        Pv = P[cov, :]
        A = np.vstack([Pv, np.ones((1, Pv.shape[1]))])
        b = np.hstack([x[cov], [1.0]])
        Gv = (A.T @ A).astype(np.float64)
        rhs = (A.T @ b).astype(np.float64)
        if LASSO:
            rhs = rhs - 0.5 * LASSO * float(np.mean(np.diag(Pv.T @ Pv)))
        Gv[np.diag_indices_from(Gv)] += 1e-6
        Rv = np.linalg.cholesky(Gv).T
        q, _ = nnls(Rv, solve_triangular(Rv.T, rhs, lower=True))
        q = q / (q.sum() + 1e-12)
        order = np.argsort(-q)
        name = os.path.basename(bcf).replace('_imputed_dog10k.bcf', '')
        print(f"=== {name}   ({100*cov.mean():.1f}% of panel SNPs covered) ===")
        for i in order[:top_n]:
            if q[i] < min_q:
                break
            print(f"    {100*q[i]:5.1f}%  {breeds[i].replace('_',' ').title()}")
        print()
    os.unlink(bed)


if __name__ == '__main__':
    main()
