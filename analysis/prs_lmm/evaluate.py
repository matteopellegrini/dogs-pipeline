#!/usr/bin/env python3
"""
Score LMM-based PRS against the current marginal-ridge PRS, breed-blocked.

    python3 evaluate.py WORKDIR BFILE_PREFIX

For every trait and every fold: train both models on the fold's training
breeds, predict the held-out breeds' dogs, aggregate to breed means, and
correlate with the breed's true score. The unit of validation is the BREED —
dogs are pseudo-replicates of their breed's phenotype, so dog-level metrics
would flatter whichever model best memorises breed membership.

Old model  = the pipeline's compute_prs_ridge: per-SNP marginal covariance
             with ridge shrinkage, no structure control.
LMM model  = GEMMA -lmm betas from the same training dogs (structure absorbed
             by the GRM, genotype-source covariate), PRS = sum(beta_j * g_j)
             over SNPs passing a p-value cut, with the cut chosen among a few
             candidates on TRAINING breeds only.
"""
import os
import sys
import numpy as np

P_CUTS = [1e-2, 1e-4, 1e-6]


def read_bed(prefix):
    fam = [l.split() for l in open(prefix + '.fam')]
    n = len(fam)
    snps = [l.split()[1] for l in open(prefix + '.bim')]
    m = len(snps)
    data = open(prefix + '.bed', 'rb').read()
    assert data[:3] == b'\x6c\x1b\x01', 'not a SNP-major bed'
    bpf = (n + 3) // 4
    raw = np.frombuffer(data, dtype=np.uint8, offset=3).reshape(m, bpf)
    codes = np.zeros((m, n), dtype=np.int8)
    for shift in range(4):
        cols = np.arange(bpf) * 4 + shift
        keep = cols < n
        codes[:, cols[keep]] = (raw >> (2 * shift))[:, keep] & 0b11
    # PLINK: 00=hom A1(2 copies), 10=het(1), 11=hom A2(0), 01=missing
    dose = np.where(codes == 0, 2, np.where(codes == 2, 1, np.where(codes == 3, 0, -1))).astype(np.float32)
    # mean-impute missing per SNP
    for j in np.where((dose < 0).any(axis=1))[0]:
        row = dose[j]
        mu = row[row >= 0].mean() if (row >= 0).any() else 0.0
        row[row < 0] = mu
    return fam, snps, dose


def ridge_betas(G, y, lambda_frac=0.1):
    """The pipeline's compute_prs_ridge, verbatim in the parts that matter."""
    Gc = G - G.mean(axis=1, keepdims=True)
    yc = y - y.mean()
    var_j = np.sum(Gc**2, axis=1)
    return np.dot(Gc, yc) / (var_j + lambda_frac * var_j.mean())


def main():
    work, bfile = sys.argv[1], sys.argv[2]
    # Optional third arg: breed list to restrict EVALUATION to, so different
    # panels can be compared on identical held-out breeds. Training is
    # untouched — only which test breeds are scored.
    only = None
    if len(sys.argv) > 3:
        only = {l.split()[0] for l in open(sys.argv[3])}
    fam, snps, dose = read_bed(bfile)   # dose: m x n
    breeds = np.array([f[0] for f in fam])
    snp_idx = {s: i for i, s in enumerate(snps)}

    cols = [l.rstrip('\n').split('\t') for l in open(f'{work}/columns.tsv')]
    pheno = np.loadtxt(f'{work}/pheno.tsv', dtype=str)
    fold_of = dict(l.split() for l in open(f'{work}/folds.tsv'))

    traits = sorted({t for t, f, c in cols})
    print(f"{'trait':<28} {'r_ridge':>8} {'r_lmm':>8} {'breeds':>7}")
    summary = []
    for trait in traits:
        full_col = next(int(c) for t, f, c in cols if t == trait and f == 'full')
        yfull = pheno[:, full_col - 1]
        have = yfull != 'NA'
        yv = np.full(len(fam), np.nan)
        yv[have] = yfull[have].astype(float)

        preds_r, preds_l, truth, breed_names = [], [], [], []
        for fold in '12345':
            col = next((int(c) for t, f, c in cols if t == trait and f == fold), None)
            if col is None:
                continue
            ytr = pheno[:, col - 1]
            tr = ytr != 'NA'
            te = have & ~tr        # dogs with truth but masked in this fold
            if te.sum() == 0 or tr.sum() < 200:
                continue
            # old model
            br = ridge_betas(dose[:, tr], yv[tr])
            pr = br @ dose[:, te]
            # LMM model
            slug = (trait + '_' + fold).replace(' ', '_').replace('/', '_')
            slug = ''.join(ch for ch in slug if ch.isalnum() or ch == '_')
            fn = f'{work}/output/lmm_{slug}_col{col}.assoc.txt'
            if not os.path.exists(fn):
                continue
            bl = np.zeros(len(snps), dtype=np.float32)
            pv = np.ones(len(snps), dtype=np.float32)
            with open(fn) as fh:
                next(fh)
                for line in fh:
                    f2 = line.split('\t')
                    i = snp_idx.get(f2[1])
                    if i is not None:
                        bl[i] = float(f2[7]); pv[i] = float(f2[11])
            # choose the p-cut on TRAINING breeds only
            best_cut, best_r = P_CUTS[0], -2
            for cut in P_CUTS:
                sel = pv <= cut
                if sel.sum() < 5:
                    continue
                pl_tr = bl[sel] @ dose[np.ix_(sel, tr)]
                bt = breeds[tr]
                bm = {b: pl_tr[bt == b].mean() for b in set(bt)}
                tv = {b: yv[tr][bt == b][0] for b in set(bt)}
                r = np.corrcoef([bm[b] for b in bm], [tv[b] for b in bm])[0, 1]
                if r > best_r:
                    best_r, best_cut = r, cut
            sel = pv <= best_cut
            pl = bl[sel] @ dose[np.ix_(sel, te)]

            bt = breeds[te]
            for b in sorted(set(bt)):
                if only is not None and b not in only:
                    continue
                mask = bt == b
                preds_r.append(pr[mask].mean())
                preds_l.append(pl[mask].mean())
                truth.append(yv[te][mask][0])
                breed_names.append(b)

        if len(truth) < 10:
            continue
        r_r = np.corrcoef(preds_r, truth)[0, 1]
        r_l = np.corrcoef(preds_l, truth)[0, 1]
        summary.append((trait, r_r, r_l, len(truth)))
        print(f"{trait:<28} {r_r:>8.3f} {r_l:>8.3f} {len(truth):>7}")

    if summary:
        mr = np.mean([s[1] for s in summary]); ml = np.mean([s[2] for s in summary])
        wins = sum(1 for s in summary if s[2] > s[1])
        print(f"\nmean breed-level r:  ridge {mr:.3f}   lmm {ml:.3f}   "
              f"(lmm better on {wins}/{len(summary)} traits)")


if __name__ == '__main__':
    main()
