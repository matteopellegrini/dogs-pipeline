#!/usr/bin/env python3
"""
Local ancestry from breed ALLELE FREQUENCIES: a diploid HMM over the production
breed panel (breed_panel/phat.npy, the same per-breed frequencies the global
lasso uses), with the lasso proportions as the ancestry prior. Stage 9c of
run_dog_pipeline.sh; the genome-wide mean of the local posteriors is the
REPORTED breed composition (user decision 2026-09-11), the lasso estimate is
kept as breed_composition_global.

Why frequencies, not haplotype copying (FLARE): the panel has 5-15 dogs per
breed across 231 breeds. Haplotype-copying models need hundreds of haplotypes
per ancestry and, on this panel, collapse onto one breed with clean data or
onto a diverse "sink" breed with noisy data (simulation, 2026-09-11; see
sim_admixed.py / score_sim_hmm.py). Frequencies are stable at small n (the
lasso is 95% top-1 on them), and the model works on genotype DOSAGES, so
low-pass phase-switch errors are not an error source at all.

Model
  candidates  K breeds with lasso weight >= floor (prior pi, renormalised)
  states      ordered ancestry pairs (i, j) for the two haplotypes -> K^2
  emission    P(dosage d | f_i, f_j) with genotype error eps, at every site,
              tempered (^1/temper) because neighbouring panel SNPs are in LD
  transition  each haplotype switches ancestry between sites with prob
              r = 1 - exp(-gen * dM) (dM in Morgans, constant 1 cM/Mb);
              on a switch the new ancestry ~ pi. The K^2 transition is the
              Kronecker product of T = (1-r) I + r 1 pi^T, so a forward step
              is alpha' = (T^T A T) o E with A the K x K forward matrix.
  outputs     per-site posterior over unordered pairs; ancestry marginals;
              genome-wide mean = proportions; painting = posterior-max pair per
              site, run-length encoded into segments with mean posterior and a
              confident flag (>= 0.5); gen chosen by likelihood over a grid.

Validated (simulated admixed dogs from held-out reference haplotypes, truth at
every site, temper 10): F1 98-99%, 3-generation doodle 84-85%, German x White
Swiss Shepherd 65-69%, six-breed supermutt 77-80% site accuracy; proportions
within 3 points; phase noise has no effect; 5% genotype error costs <= 3 pts.

Usage
  lai_hmm.py <query> <breed_result.json> <out_prefix> [--ds] [--panel DIR]
             [--floor 0.02] [--eps 0.05] [--gen auto|N] [--temper 10]
             [--write-breed-result]
  <query>: a VCF(.gz) with GT at the panel sites (any REF/ALT orientation),
           or with --ds the pipeline's relatives_ds.tsv (chrom pos ref alt DS).
  --write-breed-result: rewrite <breed_result.json> in place with the HMM
           proportions as breed_composition (lasso kept as _global).
Writes <out_prefix>.lai.json (the pipeline names it local_ancestry.json).
"""
import gzip, json, sys, os
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


def _orient(hit, ref, alt, alt_count):
    i, a1, a2 = hit
    if a1 == alt and a2 == ref: return i, alt_count
    if a1 == ref and a2 == alt: return i, 2 - alt_count
    return i, -1


def load_query_vcf(path, sites):
    """dosage of panel allele a1 (0/1/2, -1 missing) per panel site"""
    dos = np.full(len(sites), -1, dtype=np.int8)
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        for l in f:
            if l[0] == '#': continue
            p = l.rstrip('\n').split('\t', 10)
            hit = sites.get((p[0], p[1]))
            if not hit: continue
            g = p[9].split(':')[0]
            if '.' in g: continue
            alt_count = (g[0] == '1') + (g[2] == '1') if len(g) >= 3 else (int(g) if g.isdigit() else -1)
            if alt_count < 0: continue
            i, d = _orient(hit, p[3], p[4], alt_count)
            if d >= 0: dos[i] = d
    return dos


def load_query_ds(path, sites):
    """pipeline relatives_ds.tsv: chrom pos ref alt DS (ALT dosage, GLIMPSE2 posterior mean) -> rounded a1 dosage"""
    dos = np.full(len(sites), -1, dtype=np.int8)
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        for l in f:
            p = l.rstrip('\n').split('\t')
            if len(p) < 5 or p[4] in ('', '.'): continue
            hit = sites.get((p[0], p[1]))
            if not hit: continue
            try: ds = float(p[4])
            except ValueError: continue
            i, d = _orient(hit, p[2], p[3], int(round(min(max(ds, 0.0), 2.0))))
            if d >= 0: dos[i] = d
    return dos


def site_coords(panel):
    chroms, pos = [], []
    with open(f'{panel}/sites.tsv') as f:
        next(f)
        for l in f:
            c, p = l.split()[:2]; chroms.append(c); pos.append(int(p))
    return np.array(chroms), np.array(pos)


def lasso_prior(breed, breeds, floor):
    """per-code lasso weights >= floor, renormalised; plus code -> pretty name and code -> display group"""
    w, pretty, group = {}, {}, {}
    for e in breed.get('breed_composition_raw') or []:
        w[e['breed']] = float(e['proportion']); pretty[e['breed']] = e.get('breed_name', e['breed'])
    for e in breed.get('breed_composition', []):
        comps = e.get('components') or [{'code': e.get('breed'), 'proportion': e.get('proportion', 0)}]
        for c in comps:
            code = c.get('code') or e.get('breed')
            group[code] = e.get('breed_name', code)
            if code not in w: w[code] = float(c.get('proportion', 0)); pretty.setdefault(code, e.get('breed_name', code))
    keep = {k: v for k, v in w.items() if k in set(breeds) and v >= floor}
    tot = sum(keep.values()); keep = {k: v / tot for k, v in keep.items()}
    return keep, pretty, group


def run_hmm(F, dos, chroms, pos, pi, gen, eps, temper=10.0):
    """F: sites x K frequencies of a1 for the K candidates; returns (loglik, gamma sites x K x K)"""
    K = F.shape[1]; n = len(dos)
    f = np.clip(F, 0.01, 0.99)
    fi = f[:, :, None]; fj = f[:, None, :]
    p2 = fi * fj; p1 = fi * (1 - fj) + fj * (1 - fi); p0 = (1 - fi) * (1 - fj)
    E = np.empty((n, K, K), dtype=np.float64)
    m2, m1, m0, miss = dos == 2, dos == 1, dos == 0, dos < 0
    E[m2] = p2[m2]; E[m1] = p1[m1]; E[m0] = p0[m0]; E[miss] = 1.0
    E = (1 - eps) * E + eps / 3.0
    if temper != 1.0: E = E ** (1.0 / temper)
    E[miss] = 1.0
    pi = np.asarray(pi); I = np.eye(K)
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


def paint(gamma, cands, chroms, pos, conf_min=0.5):
    """unordered-pair posterior-max painting -> segments [{chrom,start,end,anc:[a,b],posterior,confident}]"""
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
            segs.append({'chrom': str(chroms[m]), 'start': int(pos[m]), 'end': int(pos[m]),
                         'anc': list(pairs[top[m]]), '_k': int(top[m]), '_c': [conf[m]]})
    for s in segs:
        s['posterior'] = round(float(np.mean(s.pop('_c'))), 3); s.pop('_k')
        s['confident'] = s['posterior'] >= conf_min
    return segs


def rewrite_breed_result(path, breed, props, pretty, group, gen, n_conf_frac):
    """HMM proportions become breed_composition; lasso kept as breed_composition_global."""
    if 'breed_composition_global' not in breed:
        breed['breed_composition_global'] = breed['breed_composition']
    grouped = {}
    for code, p in props.items():
        d = group.get(code, pretty.get(code, code))
        g = grouped.setdefault(d, {'breed': code, 'breed_name': d, 'proportion': 0.0, 'components': []})
        g['proportion'] += p; g['components'].append({'code': code, 'proportion': round(p, 6)})
        if p > max((c['proportion'] for c in g['components'][:-1]), default=-1): g['breed'] = code
    comp = sorted(grouped.values(), key=lambda e: -e['proportion'])
    for e in comp:
        e['proportion'] = round(e['proportion'], 6)
        e['components'].sort(key=lambda c: -c['proportion'])
        if len(e['components']) == 1: e.pop('components')
    breed['breed_composition'] = comp
    breed['ancestry_method'] = ('Two-step: supervised NNLS/lasso on breed allele frequencies for the candidate '
                                'breeds and prior, then a diploid local-ancestry HMM along the genome; the reported '
                                'percentages are the genome-wide mean of the local posteriors (breed_composition), '
                                'the global estimate is kept in breed_composition_global. '
                                f'HMM gen={gen}, {100*n_conf_frac:.0f}% of the genome painted at posterior >= 0.5.')
    json.dump(breed, open(path, 'w'), indent=2)


def main():
    args = sys.argv[1:]
    q, bj, out = args[:3]
    def opt(name, default):
        return args[args.index(name) + 1] if name in args else default
    panel = opt('--panel', os.environ.get('BREED_PANEL', '/Users/matteopellegrini/Downloads/dogs/breed_panel'))
    floor = float(opt('--floor', 0.02)); eps = float(opt('--eps', 0.05)); gen_opt = opt('--gen', 'auto')
    temper = float(opt('--temper', 10.0))
    P, breeds, sites = load_panel(panel)
    chroms, pos = site_coords(panel)
    dos = load_query_ds(q, sites) if '--ds' in args else load_query_vcf(q, sites)
    breed = json.load(open(bj))
    keep, pretty, group = lasso_prior(breed, breeds, floor)
    if not keep:
        print(f'{out}: no candidate breeds at floor {floor}; local ancestry skipped'); return
    cands = list(keep); idx = [breeds.index(b) for b in cands]; pi = np.array([keep[b] for b in cands])
    F = P[:, idx].astype(np.float64)
    gens = [1, 2, 3, 5, 8, 12] if gen_opt == 'auto' else [float(gen_opt)]
    best = None
    for g in gens:
        ll, gamma = run_hmm(F, dos, chroms, pos, pi, g, eps, temper)
        if best is None or ll > best[0]: best = (ll, g, gamma)
    ll, gen, gamma = best
    marg = (gamma.sum(axis=2) + gamma.sum(axis=1)) / 2.0
    props = {b: round(float(v), 4) for b, v in zip(cands, marg.mean(axis=0))}
    segs = paint(gamma, cands, chroms, pos)
    span = sum(s['end'] - s['start'] for s in segs) or 1
    conf_frac = sum(s['end'] - s['start'] for s in segs if s['confident']) / span
    chrom_span = {}
    for s in segs:
        c = chrom_span.setdefault(s['chrom'], [s['start'], s['end']])
        c[0] = min(c[0], s['start']); c[1] = max(c[1], s['end'])
    res = {'method': 'diploid local-ancestry HMM on breed allele frequencies (breed_panel Phat), lasso prior; '
                     'segments = posterior-max unordered ancestry pair per site, run-length encoded; '
                     'confident = mean posterior >= 0.5',
           'candidates': cands, 'names': {b: pretty.get(b, b) for b in cands}, 'groups': {b: group.get(b, pretty.get(b, b)) for b in cands},
           'prior': {b: round(float(v), 4) for b, v in zip(cands, pi)}, 'gen': gen, 'loglik': round(ll, 1),
           'eps': eps, 'temper': temper, 'floor': floor, 'sites_used': int((dos >= 0).sum()),
           'proportions': props, 'confident_fraction': round(conf_frac, 3),
           'chromosomes': {c: {'start': v[0], 'end': v[1]} for c, v in chrom_span.items()},
           'n_segments': len(segs), 'segments': segs}
    json.dump(res, open(out + '.lai.json', 'w'))
    if '--write-breed-result' in args:
        rewrite_breed_result(bj, breed, props, pretty, group, gen, conf_frac)
    print(f'{out}: K={len(cands)} gen={gen} sites={res["sites_used"]} segments={len(segs)} confident={100*conf_frac:.0f}%; ' +
          ', '.join(f'{b} {v:.3f} (prior {keep[b]:.3f})' for b, v in sorted(props.items(), key=lambda kv: -kv[1])))


if __name__ == '__main__':
    main()
