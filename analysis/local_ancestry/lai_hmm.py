#!/usr/bin/env python3
"""
Local ancestry from breed ALLELE FREQUENCIES: a diploid HMM over the production
breed panel (breed_panel/phat.npy, the same per-breed frequencies the global
lasso uses), with the lasso proportions as the ancestry prior.

Why frequencies, not haplotype copying (FLARE): the panel has 5-15 dogs per
breed across 231 breeds. Haplotype-copying models need hundreds of haplotypes
per ancestry and, on this panel, collapse onto one breed with clean data or
onto a diverse "sink" breed with noisy data (simulation, 2026-09-11).
Frequencies are stable at small n (the lasso is 95% top-1 on them), and the
model works on genotype DOSAGES, so low-pass phase-switch errors are not an
error source at all.

Model
  candidates  K breeds with lasso weight >= floor (prior pi, renormalised)
  states      ordered ancestry pairs (i, j) for the two haplotypes -> K^2
  emission    P(dosage d | f_i, f_j) with genotype error eps, at every site
  transition  each haplotype switches ancestry between sites with prob
              r = 1 - exp(-gen * dM) (dM in Morgans, constant 1 cM/Mb);
              on a switch the new ancestry ~ pi. The K^2 transition is the
              Kronecker product of T = (1-r) I + r 1 pi^T, so a forward step
              is alpha' = (T^T A T) o E with A the K x K forward matrix.
  outputs     per-site posterior over unordered pairs; per-site ancestry
              marginals; genome-wide mean = final proportions; painting =
              posterior-max pair per site, run-length encoded into segments
              with their mean posterior; gen chosen by likelihood over a grid.

Usage
  python3 lai_hmm.py <query.vcf.gz> <breed_result.json> <out_prefix>
        [--panel DIR] [--floor 0.02] [--eps 0.05] [--gen auto|N]
Query: phased or unphased GT at the panel sites (Dog10K REF/ALT orientation
or any orientation - alleles are matched to the panel's a1/a2 per site).
"""
import gzip, json, sys, math
import numpy as np

def load_panel(panel):
    P = np.load(f'{panel}/phat.npy')                       # sites x breeds, freq of allele a1
    breeds = [l.strip() for l in open(f'{panel}/breeds.txt') if l.strip()]
    sites = {}
    with open(f'{panel}/sites.tsv') as f:
        next(f)
        for i, l in enumerate(f):
            c, p, a1, a2 = l.split()[:4]; sites[(c, p)] = (i, a1, a2)
    return P, breeds, sites

def load_query(path, sites):
    """dosage of panel allele a1 (0/1/2, -1 missing) per panel site, plus chrom/pos arrays"""
    n = len(sites); dos = np.full(n, -1, dtype=np.int8)
    with gzip.open(path, 'rt') as f:
        for l in f:
            if l[0] == '#': continue
            p = l.rstrip('\n').split('\t', 10)
            hit = sites.get((p[0], p[1]))
            if not hit: continue
            i, a1, a2 = hit; ref, alt = p[3], p[4]
            g = p[9].split(':')[0]
            if '.' in g: continue
            alt_count = (g[0] == '1') + (g[2] == '1') if len(g) >= 3 else int(g) if g.isdigit() else -1
            if alt_count < 0: continue
            if a1 == alt and a2 == ref: dos[i] = alt_count
            elif a1 == ref and a2 == alt: dos[i] = 2 - alt_count
    return dos

def site_coords(panel):
    chroms, pos = [], []
    with open(f'{panel}/sites.tsv') as f:
        next(f)
        for l in f:
            c, p = l.split()[:2]; chroms.append(c); pos.append(int(p))
    return np.array(chroms), np.array(pos)

def lasso_prior(breed_json, breeds, floor):
    b = json.load(open(breed_json)); w = {}
    for e in b['breed_composition']:
        comps = e.get('components') or [{'code': e.get('breed'), 'proportion': e.get('proportion', 0)}]
        for c in comps:
            code = c.get('code') or e.get('breed'); w[code] = w.get(code, 0.0) + float(c.get('proportion', 0))
    keep = {k: v for k, v in w.items() if k in set(breeds) and v >= floor}
    tot = sum(keep.values()); keep = {k: v / tot for k, v in keep.items()}
    return keep

def run_hmm(F, dos, chroms, pos, pi, gen, eps, temper=1.0):
    """F: sites x K frequencies of a1 for the K candidates; returns (loglik, gamma_pairs sites x K x K)

    temper: sites are treated as independent by the HMM but neighbouring panel
    SNPs are in LD (one site per ~18 kb; within-breed LD runs to ~1 Mb), so raw
    per-site likelihoods overstate the evidence and the model over-switches.
    Raising emissions to 1/temper deflates the evidence to roughly one
    independent observation per `temper` sites (LAMP-style windowing)."""
    K = F.shape[1]; n = len(dos)
    f = np.clip(F, 0.01, 0.99)
    # emission E[m, i, j] = P(dos_m | f_i, f_j) with error eps (mixture with the genotype's population probability)
    fi = f[:, :, None]; fj = f[:, None, :]
    p2 = fi * fj; p1 = fi * (1 - fj) + fj * (1 - fi); p0 = (1 - fi) * (1 - fj)
    E = np.empty((n, K, K), dtype=np.float64)
    m2, m1, m0, miss = dos == 2, dos == 1, dos == 0, dos < 0
    E[m2] = p2[m2]; E[m1] = p1[m1]; E[m0] = p0[m0]; E[miss] = 1.0
    E = (1 - eps) * E + eps / 3.0
    if temper != 1.0: E = E ** (1.0 / temper)
    E[miss] = 1.0
    pi = np.asarray(pi); I = np.eye(K)
    # per-site haplotype transition T_m = (1-r_m) I + r_m 1 pi^T, r from distance to previous site (0 across chromosomes)
    d = np.zeros(n); same = np.r_[False, chroms[1:] == chroms[:-1]]
    d[same] = (pos[1:] - pos[:-1])[same[1:]] / 1e6 / 100.0          # Morgans at 1 cM/Mb
    r = 1 - np.exp(-gen * d); r[~same] = 1.0                        # new chromosome: redraw from prior
    alpha = np.empty((n, K, K)); c = np.empty(n)
    A = np.outer(pi, pi) * E[0]; c[0] = A.sum(); alpha[0] = A / c[0]
    for m in range(1, n):
        rm = r[m]; T = (1 - rm) * I + rm * np.outer(np.ones(K), pi)
        A = (T.T @ alpha[m - 1] @ T) * E[m]
        c[m] = A.sum(); alpha[m] = A / c[m]
    beta = np.empty((n, K, K)); beta[-1] = 1.0
    for m in range(n - 2, -1, -1):
        rm = r[m + 1]; T = (1 - rm) * I + rm * np.outer(np.ones(K), pi)
        B = T @ (beta[m + 1] * E[m + 1]) @ T.T
        beta[m] = B / c[m + 1]
    gamma = alpha * beta
    gamma /= gamma.sum(axis=(1, 2), keepdims=True)
    return float(np.log(c).sum()), gamma

def paint(gamma, cands, chroms, pos):
    """unordered-pair posterior-max painting -> segments [{chrom,start,end,anc:[a,b],posterior}]"""
    n, K, _ = gamma.shape
    unord = gamma + np.transpose(gamma, (0, 2, 1)); iu = np.triu_indices(K)
    unord[:, iu[0], iu[1]] *= np.where(iu[0] == iu[1], 0.5, 1.0)     # diagonal counted once
    best = np.stack([unord[:, i, j] for i, j in zip(*iu)], axis=1)   # sites x pairs
    top = best.argmax(axis=1); conf = best.max(axis=1)
    pairs = [(cands[i], cands[j]) for i, j in zip(*iu)]
    segs = []
    for m in range(n):
        if segs and segs[-1]['chrom'] == chroms[m] and segs[-1]['_k'] == top[m]:
            segs[-1]['end'] = int(pos[m]); segs[-1]['_c'].append(conf[m])
        else:
            segs.append({'chrom': str(chroms[m]), 'start': int(pos[m]), 'end': int(pos[m]), 'anc': list(pairs[top[m]]), '_k': int(top[m]), '_c': [conf[m]]})
    for s in segs:
        s['posterior'] = round(float(np.mean(s.pop('_c'))), 3); s.pop('_k')
    return segs

def main():
    args = sys.argv[1:]
    q, bj, out = args[:3]
    def opt(name, default):
        return args[args.index(name) + 1] if name in args else default
    panel = opt('--panel', '/Users/matteopellegrini/Downloads/dogs/breed_panel')
    floor = float(opt('--floor', 0.02)); eps = float(opt('--eps', 0.05)); gen_opt = opt('--gen', 'auto')
    temper = float(opt('--temper', 10.0))
    P, breeds, sites = load_panel(panel)
    chroms, pos = site_coords(panel)
    dos = load_query(q, sites)
    keep = lasso_prior(bj, breeds, floor)
    cands = list(keep); idx = [breeds.index(b) for b in cands]; pi = np.array([keep[b] for b in cands])
    F = P[:, idx].astype(np.float64)
    gens = [1, 2, 3, 5, 8, 12] if gen_opt == 'auto' else [float(gen_opt)]
    best = None
    for g in gens:
        ll, gamma = run_hmm(F, dos, chroms, pos, pi, g, eps, temper)
        if best is None or ll > best[0]: best = (ll, g, gamma)
    ll, gen, gamma = best
    marg = (gamma.sum(axis=2) + gamma.sum(axis=1)) / 2.0                # sites x K ancestry marginals
    props = {b: round(float(v), 4) for b, v in zip(cands, marg.mean(axis=0))}
    segs = paint(gamma, cands, chroms, pos)
    res = {'candidates': cands, 'prior': {b: round(float(v), 4) for b, v in zip(cands, pi)}, 'gen': gen, 'loglik': round(ll, 1),
           'eps': eps, 'temper': temper, 'sites_used': int((dos >= 0).sum()), 'proportions': props, 'n_segments': len(segs), 'segments': segs}
    json.dump(res, open(out + '.lai.json', 'w'))
    print(f'{out}: K={len(cands)} gen={gen} sites={res["sites_used"]} segments={len(segs)}; ' +
          ', '.join(f'{b} {v:.3f} (prior {keep[b]:.3f})' for b, v in sorted(props.items(), key=lambda kv: -kv[1])))

if __name__ == '__main__':
    main()
