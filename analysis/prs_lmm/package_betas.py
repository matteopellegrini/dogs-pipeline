#!/usr/bin/env python3
"""
Package the full-data GEMMA betas into the reference artifact Stage 11 reads.

    python3 package_betas.py SWEEP_WORKDIR BFILE_PREFIX OUT.json.gz

DENSE format (v2): every panel SNP contributes, no p-value thresholding.

Why dense. The first version chose a per-trait p-cut by breed-blocked CV,
which optimises breed-level RANKING — and sparse scores (58 SNPs for weight)
are brittle for an individual mixed-breed dog, where segregation luck at a
handful of loci dominates: a 13.6 kg dog was predicted at 44.6 kg. Re-running
the CV with dense cuts showed all-SNP scores lose almost nothing at breed
level (mean r 0.437 vs 0.444 at p<=0.1, and 14/16 traits are as good or
better dense than at the old sparse cuts), while the same dog's prediction
moves to 17.0 kg. A dense score behaves like a breed-mix average — exactly
the right prior for an individual dog.

Artifact layout: one shared `sites` table ([chrom, pos, effect_allele,
other_allele, panel_af]) plus a per-trait `beta` vector aligned to it
(0.0 where a trait's GEMMA run dropped the SNP). Stored gzipped — dense
betas for 16 traits x 131k SNPs are ~20 MB raw.

Per trait we still store: cv_r at the dense setting (5-fold breed-blocked,
honesty metadata only), ref_prs_sorted over all panel dogs (z / percentile
at run time with no genotype panel), and truth mean/sd across breeds
(pred = mu + z * sd, same scale the pipeline reports).
"""
import gzip
import json
import os
import sys
import numpy as np


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
    """Dense beta vector aligned to the bim; a1 = (effect allele, other, af)."""
    beta = np.zeros(m, dtype=np.float32)
    seen = np.zeros(m, dtype=bool)
    a1 = {}
    with open(fn) as fh:
        next(fh)
        for line in fh:
            f = line.split('\t')
            i = snp_idx.get(f[1])
            if i is not None:
                beta[i] = float(f[7])
                seen[i] = True
                a1[i] = (f[4], f[5], float(f[6]))
    return beta, seen, a1


def main():
    work, bfile, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    fam, bim, dose = read_bed(bfile)
    breeds = np.array([f[0] for f in fam])
    snp_idx = {b[1]: i for i, b in enumerate(bim)}
    m = len(bim)

    cols = [l.rstrip('\n').split('\t') for l in open(f'{work}/columns.tsv')]
    pheno = np.loadtxt(f'{work}/pheno.tsv', dtype=str)
    traits = sorted({t for t, f, c in cols})

    # site table: alleles/af from the first full-data assoc that reports the
    # SNP (GEMMA's per-run MAF filter can drop a site for a trait whose
    # phenotyped-dog subset differs — that trait just gets beta 0.0 there)
    site_ea = [None] * m
    trait_data = {}

    for trait in traits:
        full_col = next(int(c) for t, f, c in cols if t == trait and f == 'full')
        slug = ''.join(ch for ch in (trait + '_full').replace(' ', '_').replace('/', '_')
                       if ch.isalnum() or ch == '_')
        fn = f'{work}/output/lmm_{slug}_col{full_col}.assoc.txt'
        if not os.path.exists(fn):
            print(f'  skip {trait}: no full-data assoc')
            continue
        beta, seen, a1 = load_assoc(fn, snp_idx, m)
        for i, ea in a1.items():
            if site_ea[i] is None:
                site_ea[i] = ea
        trait_data[trait] = (full_col, beta, seen)

    keep = [i for i in range(m) if site_ea[i] is not None]
    print(f'{len(keep)} sites in union across {len(trait_data)} traits')

    sites = []
    for i in keep:
        chrom, pos = bim[i][0], int(bim[i][3])
        ea, oa, af = site_ea[i]
        sites.append([f'chr{chrom}' if not str(chrom).startswith('chr') else chrom,
                      pos, ea, oa, round(af, 4)])

    artifact = {'meta': {'model': 'GEMMA -lmm 1, GRM + genotype-source covariate',
                         'panel': os.path.basename(bfile),
                         'n_dogs': len(fam),
                         'format': 'dense-v2',
                         'p_cut': 1.0,
                         'effect_allele': 'a1 (beta is per copy of a1)'},
                'sites': sites,
                'traits': {}}

    for trait, (full_col, beta, seen) in trait_data.items():
        yfull = pheno[:, full_col - 1]
        have = yfull != 'NA'
        yv = np.full(len(fam), np.nan)
        yv[have] = yfull[have].astype(float)

        # honesty metric: 5-fold breed-blocked CV at the dense setting
        preds, truth = [], []
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
            bf, _, _ = load_assoc(fn, snp_idx, m)
            pl = bf @ dose[:, te]
            bt = breeds[te]
            for b in sorted(set(bt)):
                mask = bt == b
                preds.append(float(pl[mask].mean()))
                truth.append(float(yv[te][mask][0]))
        cv_r = float(np.corrcoef(preds, truth)[0, 1]) if len(truth) >= 10 else float('nan')

        # round exactly as stored, then build the reference distribution from
        # the rounded betas so a dog's score and ref_prs_sorted share a scale
        bvec = [float(f'{beta[i]:.4g}') for i in keep]
        bround = np.array(bvec, dtype=np.float32)
        ref_prs = bround @ dose[keep, :]

        bvals = {}
        for b in set(breeds[have]):
            bvals[b] = float(yv[breeds == b][0])
        vals = list(bvals.values())

        artifact['traits'][trait] = {
            'cv_r': round(cv_r, 4),
            'n_snps': int(seen[keep].sum()),
            'truth_mean': round(float(np.mean(vals)), 4),
            'truth_sd': round(float(np.std(vals)), 4),
            'n_breeds': len(vals),
            'beta': bvec,
            'ref_prs_sorted': [round(float(x), 4) for x in np.sort(ref_prs)],
        }
        print(f'  {trait:<28} cv_r={cv_r:+.3f}  snps={int(seen[keep].sum())}')

    with gzip.open(out_path, 'wt') as fh:
        json.dump(artifact, fh, separators=(',', ':'))
    print(f'\nwrote {out_path} ({os.path.getsize(out_path)/1e6:.1f} MB gz, '
          f'{len(artifact["traits"])} traits)')


if __name__ == '__main__':
    main()
