#!/usr/bin/env python3
"""
Package the full-data GEMMA betas into the reference artifact Stage 11 reads.

    python3 package_betas.py SWEEP_WORKDIR BFILE_PREFIX OUT.json

For each trait:
  - choose the p-value cut by CROSS-VALIDATION (the cut whose held-out-breed
    correlation, pooled over the 5 folds, is best) — never on the full data,
    which would be selection on the test set
  - take the full-data run's betas at that cut, with explicit effect alleles
  - compute the reference PRS for all panel dogs, stored sorted, so Stage 11
    turns a dog's score into a z and a percentile without any genotype panel
    at run time
  - store the truth mean/sd across breeds, so predicted scores stay on the
    same scale the current pipeline reports (pred = mu + z * sd)

The artifact replaces on-the-fly GWAS in Stage 11: deterministic, faster, and
the betas are structure-corrected (GEMMA LMM, GRM + genotype-source covariate)
where the old marginal ridge encoded breed membership.
"""
import json
import os
import sys
import numpy as np

P_CUTS = [1e-2, 1e-4, 1e-6]


def read_bed(prefix):
    fam = [l.split() for l in open(prefix + '.fam')]
    n = len(fam)
    bim = [l.split() for l in open(prefix + '.bim')]
    data = open(prefix + '.bed', 'rb').read()
    assert data[:3] == b'\x6c\x1b\x01'
    bpf = (n + 3) // 4
    raw = np.frombuffer(data, dtype=np.uint8, offset=3).reshape(len(bim), bpf)
    codes = np.zeros((len(bim), n), dtype=np.int8)
    for shift in range(4):
        cols = np.arange(bpf) * 4 + shift
        keep = cols < n
        codes[:, cols[keep]] = (raw >> (2 * shift))[:, keep] & 0b11
    dose = np.where(codes == 0, 2, np.where(codes == 2, 1,
                    np.where(codes == 3, 0, -1))).astype(np.float32)
    for j in np.where((dose < 0).any(axis=1))[0]:
        row = dose[j]
        mu = row[row >= 0].mean() if (row >= 0).any() else 0.0
        row[row < 0] = mu
    return fam, bim, dose


def load_assoc(fn, snp_idx, m):
    beta = np.zeros(m, dtype=np.float32)
    pval = np.ones(m, dtype=np.float32)
    a1 = {}
    with open(fn) as fh:
        next(fh)
        for line in fh:
            f = line.split('\t')
            i = snp_idx.get(f[1])
            if i is not None:
                beta[i] = float(f[7])
                pval[i] = float(f[11])
                # (effect allele, other allele, panel frequency of effect allele)
                # af lets Stage 11 substitute the panel-mean contribution
                # (beta * 2af) at sites a dog cannot be imputed at, keeping the
                # dog's score on the same scale as ref_prs_sorted.
                a1[i] = (f[4], f[5], float(f[6]))
    return beta, pval, a1


def main():
    work, bfile, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    fam, bim, dose = read_bed(bfile)
    breeds = np.array([f[0] for f in fam])
    snp_idx = {b[1]: i for i, b in enumerate(bim)}
    m = len(bim)

    cols = [l.rstrip('\n').split('\t') for l in open(f'{work}/columns.tsv')]
    pheno = np.loadtxt(f'{work}/pheno.tsv', dtype=str)
    traits = sorted({t for t, f, c in cols})

    artifact = {'meta': {'model': 'GEMMA -lmm 1, GRM + genotype-source covariate',
                         'panel': os.path.basename(bfile),
                         'n_dogs': len(fam),
                         'p_cut_selection': '5-fold breed-blocked CV',
                         'effect_allele': 'a1 (beta is per copy of a1)'},
                'traits': {}}

    for trait in traits:
        full_col = next(int(c) for t, f, c in cols if t == trait and f == 'full')
        yfull = pheno[:, full_col - 1]
        have = yfull != 'NA'
        yv = np.full(len(fam), np.nan)
        yv[have] = yfull[have].astype(float)

        # CV: pooled held-out predictions per cut
        pooled = {cut: ([], []) for cut in P_CUTS}
        for fold in '12345':
            col = next((int(c) for t, f, c in cols if t == trait and f == fold), None)
            if col is None:
                continue
            tr = pheno[:, col - 1] != 'NA'
            te = have & ~tr
            if te.sum() == 0:
                continue
            slug = ''.join(ch for ch in (trait + '_' + fold).replace(' ', '_').replace('/', '_')
                           if ch.isalnum() or ch == '_')
            fn = f'{work}/output/lmm_{slug}_col{col}.assoc.txt'
            if not os.path.exists(fn):
                continue
            beta, pval, _ = load_assoc(fn, snp_idx, m)
            bt = breeds[te]
            for cut in P_CUTS:
                sel = pval <= cut
                if sel.sum() < 5:
                    continue
                pl = beta[sel] @ dose[np.ix_(sel, te)]
                for b in sorted(set(bt)):
                    mask = bt == b
                    pooled[cut][0].append(float(pl[mask].mean()))
                    pooled[cut][1].append(float(yv[te][mask][0]))

        cv = {}
        for cut, (p, t) in pooled.items():
            if len(t) >= 10:
                cv[cut] = float(np.corrcoef(p, t)[0, 1])
        if not cv:
            print(f'  skip {trait}: no CV signal')
            continue
        best_cut = max(cv, key=cv.get)

        slug = ''.join(ch for ch in (trait + '_full').replace(' ', '_').replace('/', '_')
                       if ch.isalnum() or ch == '_')
        fn = f'{work}/output/lmm_{slug}_col{full_col}.assoc.txt'
        beta, pval, a1 = load_assoc(fn, snp_idx, m)
        sel = np.where(pval <= best_cut)[0]

        ref_prs = beta[sel] @ dose[np.ix_(sel, np.arange(len(fam)))]
        # truth stats across BREEDS (matching current Stage 11's scale mapping)
        bvals = {}
        for b in set(breeds[have]):
            bvals[b] = float(yv[breeds == b][0])
        vals = list(bvals.values())

        snps = []
        for i in sel:
            chrom, snpid, _, pos = bim[i][0], bim[i][1], bim[i][2], bim[i][3]
            ea, oa, af = a1.get(int(i), (bim[i][4], bim[i][5], 0.0))
            snps.append([f'chr{chrom}' if not str(chrom).startswith('chr') else chrom,
                         int(pos), ea, oa, round(float(beta[i]), 6), round(af, 4)])

        artifact['traits'][trait] = {
            'p_cut': best_cut,
            'cv_r': round(cv[best_cut], 4),
            'n_snps': len(snps),
            'truth_mean': round(float(np.mean(vals)), 4),
            'truth_sd': round(float(np.std(vals)), 4),
            'n_breeds': len(vals),
            'snps': snps,
            'ref_prs_sorted': [round(float(x), 4) for x in np.sort(ref_prs)],
        }
        print(f'  {trait:<28} cut={best_cut:<8} cv_r={cv[best_cut]:+.3f}  snps={len(snps)}')

    with open(out_path, 'w') as fh:
        json.dump(artifact, fh, separators=(',', ':'))
    print(f'\nwrote {out_path} ({os.path.getsize(out_path)/1e6:.1f} MB, '
          f'{len(artifact["traits"])} traits)')


if __name__ == '__main__':
    main()
