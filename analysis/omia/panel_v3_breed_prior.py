#!/usr/bin/env python3
"""Panel v3, step 2b: genotype calls from a dog's own reads plus a breed-aware prior.

Most catalogue variants are breed founder mutations that no dog in the imputation panel carries, so they cannot be
imputed; at ~1x a dog usually has 0-1 reads at the site. Instead of a flat ">= 5 reads" rule, each dog's genotype
posterior combines

  * its reads at the site (ref/alt counts; per-read error EPS by variant type, higher for homopolymer indels), and
  * a prior from its breed composition: q_dog = sum_b w_dog,b * f_b, with the dog's inbreeding (F_ROH) moving
    mass to homozygotes: P(0) = (1-q)^2 (1-F) + (1-q) F, P(1) = 2q(1-q)(1-F), P(2) = q^2 (1-F) + q F.

The per-breed frequencies f_b are fitted to the whole cohort's reads by EM (alleles attributed to breeds in
proportion to w_b f_b), shrunk towards a prior mean:
  * breeds the literature names for this variant (OMIA / panel breed lists): mean M_DOC, strength K_DOC;
  * every other breed: a site-wide background frequency beta, itself fitted with a Beta(0.5, B_HYPER) hyperprior
    (B_HYPER smaller when the literature says the variant also occurs in mixed breeds), strength K_BG.
So a Labrador with no reads at a Labrador founder mutation stays "uncertain", while a Chihuahua with no reads at
the same site is "not detected" with a stated residual probability.

Calls per dog and site: carrier if P(>=1 copy) >= CALL_HI, clear if P(>=1 copy) <= CALL_LO, otherwise uncertain.

  python3 panel_v3_breed_prior.py fit      WORK_DIR   # cohort fit + calls -> WORK_DIR/bp_*.tsv
  python3 panel_v3_breed_prior.py validate WORK_DIR   # Dog10K sites vs GLIMPSE, and read thinning of deep kits

WORK_DIR holds the step-2a outputs (cohort_reads.tsv.gz, cohort_imputed.tsv.gz, v3_small.norm_ids.tsv,
kit_inputs.json, panel_v23_breeds.txt, d10k_match.pos).
"""
import csv, gzip, json, os, re, subprocess, sys, collections
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from breed_names import BreedMatcher, MIXED
from panel_v3_carriers import load_matrix

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
FASTA = os.path.join(ROOT, 'canFam4.fa')
EPS = {'snv': 0.002, 'indel': 0.005, 'homopolymer': 0.03}
HOMOPOLYMER_RUN = 6        # indel next to >= this many identical reference bases
M_DOC, K_DOC = 0.05, 20.0  # documented breeds: prior mean, pseudo-alleles
K_BG = 100.0               # other breeds: pseudo-alleles around the background frequency
# Beta(0.5, B) on the background frequency. B=4000 under-predicted carriers 2-3x on the Dog10K validation sites
# (P 0.11 -> observed 0.29); B=400 is calibrated (0.11 -> 0.09, 0.32 -> 0.33). Swept 2026-11-03.
B_HYPER = {'none': 400.0, 'documented': 400.0, 'mixed': 400.0}
# CALL_LO 0.002: on Dog10K rare sites 15/268 true carriers fall under it (75 at 0.01); clear-call error 0.017%.
CALL_HI, CALL_LO = 0.95, 0.002
EM_ITERS = 100

# ── inputs ────────────────────────────────────────────────────────────────────────────────────────────────────
def site_breeds():
    """site_id -> list of free-text breed names from the catalogue or the curated candidate tables."""
    kv = os.path.join(ROOT, 'known_variants')
    cat = json.load(open(os.path.join(ROOT, 'reference_json/omia_variants.json')))['variants']
    cand = list(csv.DictReader(open(os.path.join(kv, 'catalogue_expansion_candidates.tsv')), delimiter='\t'))
    by_omia = {r['omia_id']: r['breeds'] for r in cand if r['omia_id']}
    curated = {}
    for f in ('t4_priority_curated.tsv', 't5_curated.tsv'):
        for r in csv.DictReader(open(os.path.join(kv, 'panel_v3', f)), delimiter='\t'):
            curated[(r['gene'], r['pos'])] = r['breeds']
    out = {}
    for s in csv.DictReader(open(os.path.join(kv, 'panel_v3/panel_v3_sites.tsv')), delimiter='\t'):
        sid = s['site_id']
        if sid.startswith('cat'):
            b = cat[int(sid[3:])].get('variant_breed') or ''
        else:
            b = by_omia.get(sid.split('_')[-1]) or curated.get((s['gene'], s['pos'])) or ''
        out[sid] = (s, [x.strip() for x in re.split(r'[;,]\s*', b) if x.strip()])
    return out

def homopolymer_len(chrom, pos, ref, alt):
    """Longest run of one base in the 12 bp after the anchor base (where an indel can slide)."""
    try:
        seq = subprocess.run(['samtools', 'faidx', FASTA, f'{chrom}:{int(pos) + 1}-{int(pos) + 12}'],
                             capture_output=True, text=True).stdout.split('\n', 1)[1].replace('\n', '').upper()
    except Exception:
        return 0
    return max((len(m.group(0)) for m in re.finditer(r'(.)\1*', seq)), default=0)

def load(work):
    samples, sites, R, A = load_matrix(os.path.join(work, 'cohort_reads.tsv.gz'))
    # the rerun lists repeat the 96 DOGS-Gen dogs and COSMO3 under internal names (same result folder): keep one
    # row per folder, preferring the one with a real kit barcode, so no dog is counted twice in the EM
    pub, bc = {}, {}
    for line in open(os.path.join(work, 'manifest.tsv')):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3:
            pub[p[0]], bc[p[0]] = p[2], p[1]
    order = sorted(range(len(samples)), key=lambda j: not bc.get(samples[j], '').isdigit())
    seen, keep = set(), []
    for j in order:
        key = pub.get(samples[j]) or samples[j]
        if key == '-' or key not in seen:
            seen.add(key); keep.append(j)
    keep = sorted(keep)
    samples = [samples[j] for j in keep]; R = R[:, keep]; A = A[:, keep]
    ids = {}
    for line in open(os.path.join(work, 'v3_small.norm_ids.tsv')):
        c, p, ref, alt, sid = line.rstrip('\n').split('\t')[:5]
        ids[(c, p, ref, alt)] = sid
    sid = [ids.get(s) for s in sites]
    ki = json.load(open(os.path.join(work, 'kit_inputs.json')))
    names = dict(ki['breed_names'])
    codes = open(os.path.join(work, 'panel_v23_breeds.txt')).read().split()
    for k in ki['kits'].values():
        for c in (k.get('breeds') or {}):
            if c not in codes:
                codes.append(c)
    for c in codes:
        names.setdefault(c, c.replace('_', ' ').title())
    ci = {c: i for i, c in enumerate(codes)}
    W = np.zeros((len(samples), len(codes)))
    F = np.zeros(len(samples)); depth = np.zeros(len(samples))
    for j, s in enumerate(samples):
        k = ki['kits'].get(s) or {}
        for c, p in (k.get('breeds') or {}).items():
            W[j, ci[c]] = p
        F[j] = min(max(k.get('f_roh') or 0.0, 0.0), 0.5)
        depth[j] = k.get('depth') or np.nan
    tot = W.sum(1, keepdims=True)
    W = np.where(tot > 0, W / np.where(tot > 0, tot, 1), 1.0 / len(codes))  # no breed result -> flat
    return dict(samples=samples, sites=sites, sid=sid, R=R, A=A, W=W, F=F, depth=depth, codes=codes, names=names)

# ── model ─────────────────────────────────────────────────────────────────────────────────────────────────────
def read_lik(r, a, eps):
    p = np.array([eps, 0.5, 1 - eps])
    ll = a[:, None] * np.log(p)[None, :] + r[:, None] * np.log(1 - p)[None, :]
    ll -= ll.max(1, keepdims=True)
    return np.exp(ll)

def geno_prior(q, F):
    q = np.clip(q, 1e-9, 1 - 1e-9)
    return np.stack([(1 - q) ** 2 * (1 - F) + (1 - q) * F, 2 * q * (1 - q) * (1 - F), q ** 2 * (1 - F) + q * F], 1)

def fit_site(r, a, W, F, eps, doc_mask, hyper_b, fit_rows=None):
    """EM for per-breed allele frequencies. Returns f (n_breeds,), background beta."""
    m = (r + a) > 0 if fit_rows is None else ((r + a) > 0) & fit_rows
    r, a, W, F = r[m], a[m], W[m], F[m]
    lik = read_lik(r, a, eps)
    beta = 1e-4
    f = np.where(doc_mask, M_DOC, beta)
    for _ in range(EM_ITERS):
        q = W @ f
        post = lik * geno_prior(q, F)
        post /= post.sum(1, keepdims=True)
        e_alt = post[:, 1] + 2 * post[:, 2]
        qc = np.clip(q, 1e-12, 1 - 1e-12)
        alt_b = (e_alt / qc) @ W * f
        ref_b = ((2 - e_alt) / (1 - qc)) @ W * (1 - f)
        bg = ~doc_mask
        beta = (alt_b[bg].sum() + 0.5) / (alt_b[bg].sum() + ref_b[bg].sum() + 0.5 + hyper_b)
        prior_mean = np.where(doc_mask, M_DOC, beta)
        k = np.where(doc_mask, K_DOC, K_BG)
        f_new = (alt_b + k * prior_mean) / (alt_b + ref_b + k)
        if np.max(np.abs(f_new - f)) < 1e-7:
            f = f_new; break
        f = f_new
    return f, beta

def posterior(r, a, W, F, eps, f):
    post = read_lik(r, a, eps) * geno_prior(W @ f, F)
    return post / post.sum(1, keepdims=True)

def site_setup(D, bm, sb):
    """Per row: (eps, doc_mask, hyper key, documented codes, variant kind)."""
    cidx = {c: i for i, c in enumerate(D['codes'])}
    out = []
    for (c, p, ref, alt), sid in zip(D['sites'], D['sid']):
        kind = 'snv' if len(ref) == 1 and len(alt) == 1 else 'indel'
        eps = EPS[kind]
        if kind == 'indel' and homopolymer_len(c, p, ref, alt) >= HOMOPOLYMER_RUN:
            kind, eps = 'homopolymer', EPS['homopolymer']
        doc = set(); mixed = False
        for name in (sb.get(sid, (None, []))[1] if sid else []):
            hit = bm.match(name)
            mixed |= MIXED in hit
            doc |= {h for h in hit if h != MIXED}
        mask = np.zeros(len(D['codes']), bool)
        for h in doc:
            if h in cidx:
                mask[cidx[h]] = True
        hyper = 'mixed' if mixed else ('documented' if doc else 'none')
        out.append(dict(eps=eps, mask=mask, hyper=hyper, doc=sorted(doc), kind=kind))
    return out

def call_of(p):
    pc = p[:, 1] + p[:, 2]
    return np.where(pc >= CALL_HI, np.where(p[:, 2] >= 0.9, 'two_copies', 'carrier'),
                    np.where(pc <= CALL_LO, 'clear', 'uncertain'))

# ── commands ──────────────────────────────────────────────────────────────────────────────────────────────────
def cmd_fit(work):
    D = load(work)
    sb = site_breeds()
    bm = BreedMatcher(D['names'])
    setup = site_setup(D, bm, sb)
    n = len(D['samples'])
    pc_all = np.zeros((len(D['sites']), n), np.float32)
    with open(os.path.join(work, 'bp_sites.tsv'), 'w') as fo:
        fo.write('site_id\tchrom\tpos\tref\talt\tgene\tcondition\treport_class\tkind\tdocumented_breeds\tprior\t'
                 'background_af\ttop_breed_af\tkits_with_reads\tcarrier\ttwo_copies\tuncertain\tclear\n')
        for i, ((c, p, ref, alt), sid) in enumerate(zip(D['sites'], D['sid'])):
            s = setup[i]
            r, a = D['R'][i], D['A'][i]
            f, beta = fit_site(r, a, D['W'], D['F'], s['eps'], s['mask'], B_HYPER[s['hyper']])
            post = posterior(r, a, D['W'], D['F'], s['eps'], f)
            pc_all[i] = post[:, 1] + post[:, 2]
            calls = collections.Counter(call_of(post))
            top = np.argsort(-f)[:3]
            meta = sb.get(sid, ({}, []))[0] if sid else {}
            fo.write('\t'.join(map(str, [sid or '', c, p, ref, alt, meta.get('gene', ''), meta.get('condition', ''),
                     meta.get('report_class', ''), s['kind'], ','.join(s['doc']), s['hyper'], f'{beta:.2e}',
                     '; '.join(f"{D['names'].get(D['codes'][t], D['codes'][t])} {f[t]:.3f}" for t in top),
                     int(((r + a) > 0).sum()), calls['carrier'], calls['two_copies'], calls['uncertain'],
                     calls['clear']])) + '\n')
    np.savez_compressed(os.path.join(work, 'bp_p_carrier.npz'), p=pc_all, samples=np.array(D['samples']),
                        sid=np.array([x or '' for x in D['sid']]))
    print(f"fit {len(D['sites'])} sites x {n} kits -> {work}/bp_sites.tsv, bp_p_carrier.npz")

def gp_truth(work, samples):
    """(chrom,pos) -> array of GLIMPSE genotypes (0/1/2, -1 when max GP < 0.99) aligned to samples."""
    with gzip.open(os.path.join(work, 'cohort_imputed.tsv.gz'), 'rt') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        col = {s: i for i, s in enumerate(hdr)}
        out = {}
        for line in f:
            p = line.rstrip('\n').split('\t')
            g = np.full(len(samples), -1, int)
            for j, s in enumerate(samples):
                if s not in col: continue
                v = p[col[s]]
                if ':' not in v: continue
                gt, gp = v.split(':')
                try:
                    if max(float(x) for x in gp.split(',')) < 0.99: continue
                except ValueError:
                    continue
                g[j] = sum(x == '1' for x in re.split(r'[|/]', gt))
            out[(p[0], p[1], p[2], p[3])] = g
    return out

def cmd_validate(work):
    D = load(work)
    sb = site_breeds(); bm = BreedMatcher(D['names']); setup = site_setup(D, bm, sb)
    truth = gp_truth(work, D['samples'])
    rng = np.random.default_rng(1)
    deep = np.nan_to_num(D['depth']) >= 3
    rep = collections.defaultdict(collections.Counter)
    cal = []  # (p_carrier, truth carrier) for calibration
    for i, site in enumerate(D['sites']):
        t = truth.get(site)
        if t is None: continue
        s = setup[i]
        r, a = D['R'][i], D['A'][i]
        # A. Dog10K sites vs GLIMPSE: fit on everybody's reads (as in production), score dogs with a confident GLIMPSE call
        f, _ = fit_site(r, a, D['W'], D['F'], s['eps'], s['mask'], B_HYPER[s['hyper']])
        post = posterior(r, a, D['W'], D['F'], s['eps'], f)
        calls = call_of(post); ok = t >= 0
        rare = (t[ok] > 0).mean() < 0.02  # carrier frequency like a founder mutation
        for c, g in zip(calls[ok], t[ok]):
            rep['glimpse'][(c, 'truth_carrier' if g > 0 else 'truth_clear')] += 1
            if rare:
                rep['glimpse_rare_sites'][(c, 'truth_carrier' if g > 0 else 'truth_clear')] += 1
        cal += list(zip(post[ok, 1] + post[ok, 2], t[ok] > 0))
        # B. deep kits thinned to ~1x; frequencies refitted without them so they cannot inform their own prior
        f2, _ = fit_site(r, a, D['W'], D['F'], s['eps'], s['mask'], B_HYPER[s['hyper']], fit_rows=~deep)
        keep = 1.0 / np.maximum(np.nan_to_num(D['depth'], nan=1.0), 1.0)
        rt = rng.binomial(r, keep); at = rng.binomial(a, keep)
        pt = posterior(rt, at, D['W'], D['F'], s['eps'], f2)
        ct = call_of(pt); okd = deep & ok
        for c, g in zip(ct[okd], t[okd]):
            rep['thinned_deep'][(c, 'truth_carrier' if g > 0 else 'truth_clear')] += 1
    for k, cnt in rep.items():
        tc = sum(v for (c, tr), v in cnt.items() if tr == 'truth_carrier')
        tn = sum(v for (c, tr), v in cnt.items() if tr == 'truth_clear')
        print(f'== {k}: {tc} truth-carrier and {tn} truth-clear dog x site pairs')
        for c in ('carrier', 'two_copies', 'uncertain', 'clear'):
            print(f"   called {c:11s} carrier-truth {cnt[(c, 'truth_carrier')]:7d}  clear-truth {cnt[(c, 'truth_clear')]:7d}")
        called_clear = cnt[('clear', 'truth_carrier')] + cnt[('clear', 'truth_clear')]
        called_carr = sum(cnt[(c, tr)] for c in ('carrier', 'two_copies') for tr in ('truth_carrier', 'truth_clear'))
        if called_clear:
            print(f"   'clear' calls that are truly carriers: {cnt[('clear', 'truth_carrier')] / called_clear:.4%}")
        if called_carr:
            fp = sum(cnt[(c, 'truth_clear')] for c in ('carrier', 'two_copies'))
            print(f"   carrier calls that are truly clear:   {fp / called_carr:.4%}")
        if tc:
            print(f"   truth carriers called carrier: {sum(cnt[(c, 'truth_carrier')] for c in ('carrier', 'two_copies')) / tc:.1%}, "
                  f"left uncertain: {cnt[('uncertain', 'truth_carrier')] / tc:.1%}, missed as clear: {cnt[('clear', 'truth_carrier')] / tc:.2%}")
    pc, tr = np.array([x[0] for x in cal]), np.array([x[1] for x in cal])
    print('== calibration (GLIMPSE truth): predicted P(carrier) bin -> observed carrier fraction')
    for lo, hi in [(0, .001), (.001, .01), (.01, .05), (.05, .2), (.2, .5), (.5, .8), (.8, .95), (.95, 1.01)]:
        m = (pc >= lo) & (pc < hi)
        if m.sum():
            print(f'   [{lo:.3f},{hi:.3f})  n={m.sum():7d}  mean predicted {pc[m].mean():.4f}  observed {tr[m].mean():.4f}')

if __name__ == '__main__':
    cmd, work = sys.argv[1], sys.argv[2]
    {'fit': cmd_fit, 'validate': cmd_validate}[cmd](work)
