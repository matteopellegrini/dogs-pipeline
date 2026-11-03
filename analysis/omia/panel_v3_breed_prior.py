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
  python3 panel_v3_breed_prior.py sv       WORK_DIR   # structural variants from sv_counts.tsv (sv_evidence.py)

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
M_DOC, K_DOC = 0.05, float(os.environ.get('PV3_K_DOC', 5.0))  # documented breeds: prior mean, pseudo-alleles (20 held CDDY at 0.15 in Dachshunds whose reads say ~1; 5 validates as well or better)
K_BG = 100.0               # other breeds: pseudo-alleles around the background frequency
# Beta(0.5, B) on the background frequency. B=4000 under-predicted carriers 2-3x on the Dog10K validation sites
# (P 0.11 -> observed 0.29); B=400 is calibrated (0.11 -> 0.09, 0.32 -> 0.33). Swept 2026-11-03.
B_HYPER = {'none': 400.0, 'documented': 400.0, 'mixed': 400.0}
B_HYPER_HOMOPOLYMER = 40000.0
# CALL_LO 0.002: on Dog10K rare sites 15/268 true carriers fall under it (75 at 0.01); clear-call error 0.017%.
CALL_HI, CALL_LO = 0.95, 0.002
EM_ITERS = 100
X_PAR_END = 6_650_000      # canFam4 chrX pseudoautosomal region; males are hemizygous beyond it
MT_MIN_DEPTH, MT_HOMOPLASMIC, MT_CLEAR = 20, 0.8, 0.05  # chrM: allele fraction calls, no prior
# Breed lists for entries with none in the catalogue or OMIA (OMIA 632 lists no breeds)
MANUAL_BREEDS = {
    'NHEJ1': 'Collie, Shetland Sheepdog, Border Collie, Australian Shepherd, Nova Scotia Duck Tolling Retriever, '
             'Lancashire Heeler, Long-Haired Whippet, Boykin Spaniel',  # collie eye anomaly, Parker et al. 2007
}
# sites whose OMIA entry the gene-name join misses ("FGF4 retrogene in CFA12"): OMIA's breed list wins over ours
SITE_OMIA = {'cat093': '855', 'cat469': '694'}
K_EPS = 500.0  # indels: per-site read error estimated by EM, shrunk to EPS[kind] with this many pseudo-reads
# Held out of calls until their coordinates/alleles are re-curated: cohort frequencies impossible for a disease allele
HOLD = {
    'new_ABCB1_442': 'ABCB1 A>C at 38% across breeds: CanFam3 -> canFam4 position or allele under review',
    'cat159': 'BTBD17 insG at 20-28% across breeds: OMIA-resolved allele under review',
    'cat396': 'IDUA 270-bp deletion: both breakpoints sit in a GC-rich tandem repeat, so repeat-length variation in '
              'every breed mimics the junction in short reads; needs a targeted assay',
    'new_CYTB_52': 'm.14474G>A homoplasmic in ~7% of dogs across 27 breeds (pit bulls, chihuahuas, shepherds) and in none '
                   'of 32 cattle dogs and shelties: a maternal-lineage polymorphism, not a reportable disease variant',
}

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
    # OMIA's dog variant table (download of 2026-11-01): breed lists for entries where ours are empty
    omia = list(csv.DictReader(open(os.path.join(kv, 'omia_download_2026-11-01/omia_dog_variants.tsv')), delimiter='\t'))
    omia_by_id = {r['OMIA Variant ID']: r for r in omia}
    omia_by_gene = collections.defaultdict(list)
    for r in omia:
        omia_by_gene[r['Gene']].append(r)
    resolved = {r['site_id']: r['omia_variant_id'] for r in
                csv.DictReader(open(os.path.join(kv, 'panel_v3/catalogue_alleles_resolved.tsv')), delimiter='\t')}
    def omia_type_ok(row, kind):
        t = row['Variant Type']
        return ('substitution' in t) if kind == 'snv' else ('(<=20)' in t) if kind == 'indel' else ('gross' in t or '(>20)' in t)
    def omia_breeds(sid, s):
        oid = SITE_OMIA.get(sid) or resolved.get(sid) or (sid.split('_')[-1] if not sid.startswith('cat') else '')
        if oid in omia_by_id:
            return omia_by_id[oid]['Breed(s)'], f'OMIA variant {oid}'
        gene = {'FVII': 'F7', 'MOP': 'MPO', 'FGF4 retrogene': 'FGF4', 'IGF1-AS (LOC111090066)': 'IGF1'}.get(s['gene'], s['gene'])
        rows = [r for r in omia_by_gene.get(gene, []) if omia_type_ok(r, s['kind'])] or omia_by_gene.get(gene, [])
        if rows:
            return ' '.join(r['Breed(s)'] for r in rows), 'OMIA gene ' + ','.join(r['OMIA Variant ID'] for r in rows)
        return '', ''
    out = {}
    for s in csv.DictReader(open(os.path.join(kv, 'panel_v3/panel_v3_sites.tsv')), delimiter='\t'):
        sid = s['site_id']
        if sid.startswith('cat'):
            b = cat[int(sid[3:])].get('variant_breed') or ''
        else:
            b = by_omia.get(sid.split('_')[-1]) or curated.get((s['gene'], s['pos'])) or ''
        src = 'catalogue' if b.strip() else ''
        if not b.strip() or sid in SITE_OMIA:
            b, src = omia_breeds(sid, s)
        if not b.strip() and s['gene'] in MANUAL_BREEDS:
            b, src = MANUAL_BREEDS[s['gene']], 'manual'
        s = dict(s, breed_source=src)
        out[sid] = (s, [x.strip() for x in re.split(r'\(Dog\)|[;,]', b) if x.strip()])
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
        ids[(c, p, ref.upper(), alt.upper())] = sid  # soft-masked reference gives lowercase alleles (cat083)
    sid = [ids.get((c, p, ref.upper(), alt.upper())) for c, p, ref, alt in sites]
    # chrM sites: counts taken without the MAPQ filter (cluster/panel_v3_chrM.sh), the cohort matrix has none
    mt = os.path.join(work, 'chrM_counts.tsv')
    if os.path.exists(mt):
        col = {x: j for j, x in enumerate(samples)}
        cnt = collections.defaultdict(dict)
        for line in open(mt):
            k, c, p, ref, alt, nr, na = line.rstrip('\n').split('\t')
            cnt[(c, p, ref.upper(), alt.upper())][k] = (int(nr), int(na))
        for key, per in cnt.items():
            sid_mt = ids.get(key)
            if key not in [(c, p, r.upper(), a.upper()) for c, p, r, a in sites]:
                sites.append(key); sid.append(sid_mt)
                R = np.vstack([R, np.zeros(len(samples), R.dtype)]); A = np.vstack([A, np.zeros(len(samples), A.dtype)])
            i = [(c, p, r.upper(), a.upper()) for c, p, r, a in sites].index(key)
            R[i] = 0; A[i] = 0
            for k, (nr, na) in per.items():
                if k in col:
                    R[i, col[k]], A[i, col[k]] = nr, na
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
    male = np.zeros(len(samples), bool)
    sx = os.path.join(work, 'kit_sex.tsv')
    if os.path.exists(sx):
        sexes = dict(line.rstrip('\n').split('\t')[:2] for line in open(sx) if '\t' in line)
        male = np.array([sexes.get(x) == 'M' for x in samples])
    return dict(samples=samples, sites=sites, sid=sid, R=R, A=A, W=W, F=F, depth=depth, codes=codes, names=names,
                male=male)

def hemizygous(D, chrom, pos):
    """Rows (dogs) that carry one copy of this site: males beyond the X pseudoautosomal region."""
    if chrom == 'chrX' and int(pos) > X_PAR_END:
        return D['male']
    return np.zeros(len(D['samples']), bool)

# ── model ─────────────────────────────────────────────────────────────────────────────────────────────────────
def read_lik(r, a, eps):
    p = np.array([eps, 0.5, 1 - eps])
    ll = a[:, None] * np.log(p)[None, :] + r[:, None] * np.log(1 - p)[None, :]
    ll -= ll.max(1, keepdims=True)
    return np.exp(ll)

def geno_prior(q, F, hap=None):
    """Columns: 0, 1, 2 alt copies. Hemizygous rows (hap) use column 2 for their single alt copy."""
    q = np.clip(q, 1e-9, 1 - 1e-9)
    pr = np.stack([(1 - q) ** 2 * (1 - F) + (1 - q) * F, 2 * q * (1 - q) * (1 - F), q ** 2 * (1 - F) + q * F], 1)
    if hap is not None and hap.any():
        pr[hap] = np.stack([1 - q[hap], np.zeros(hap.sum()), q[hap]], 1)
    return pr

def fit_site(r, a, W, F, eps, doc_mask, hyper_b, fit_rows=None, hap=None, fit_eps=False, lik_fn=None, informative=None):
    """EM for per-breed allele frequencies (and, for indels, the site's read error).
    lik_fn(eps, rows) -> (n_rows x 3) replaces the read likelihood (structural variants: breakpoints x depth);
    informative marks dogs with any evidence (default: >= 1 read). Returns f (n_breeds,), background beta, eps."""
    eps0 = eps
    inf = ((r + a) > 0) if informative is None else informative
    m = inf if fit_rows is None else inf & fit_rows
    rows = np.where(m)[0]
    hap = np.zeros(len(r), bool) if hap is None else hap
    r, a, W, F, hap = r[m], a[m], W[m], F[m], hap[m]
    lik_of = (lambda e: lik_fn(e, rows)) if lik_fn is not None else (lambda e: read_lik(r, a, e))
    lik = lik_of(eps)
    n_al = np.where(hap, 1.0, 2.0)
    beta = 1e-4
    f = np.where(doc_mask, M_DOC, beta)
    for _ in range(EM_ITERS):
        q = W @ f
        if fit_eps:
            lik = lik_of(eps)
        post = lik * geno_prior(q, F, hap)
        post /= post.sum(1, keepdims=True)
        e_alt = np.where(hap, post[:, 2], post[:, 1] + 2 * post[:, 2])
        qc = np.clip(q, 1e-12, 1 - 1e-12)
        alt_b = (e_alt / qc) @ W * f
        ref_b = ((n_al - e_alt) / (1 - qc)) @ W * (1 - f)
        bg = ~doc_mask
        beta = (alt_b[bg].sum() + 0.5) / (alt_b[bg].sum() + ref_b[bg].sum() + 0.5 + hyper_b)
        prior_mean = np.where(doc_mask, M_DOC, beta)
        k = np.where(doc_mask, K_DOC, K_BG)
        f_new = (alt_b + k * prior_mean) / (alt_b + ref_b + k)
        if fit_eps:  # reads disagreeing with a homozygous genotype are errors
            err = (post[:, 0] * a).sum() + (post[:, 2] * r).sum()
            tot = ((post[:, 0] + post[:, 2]) * (r + a)).sum()
            eps = float(np.clip((err + K_EPS * eps0) / (tot + K_EPS), 1e-4, 0.2))
        if np.max(np.abs(f_new - f)) < 1e-7:
            f = f_new; break
        f = f_new
    return f, beta, eps

def posterior(r, a, W, F, eps, f, hap=None):
    post = read_lik(r, a, eps) * geno_prior(W @ f, F, hap)
    return post / post.sum(1, keepdims=True)

def mt_posterior(r, a):
    """chrM: one (heteroplasmic) genome per dog. Column 2 = homoplasmic alt, column 1 = heteroplasmic."""
    n = r + a
    frac = np.where(n > 0, a / np.maximum(n, 1), 0.0)
    post = np.zeros((len(r), 3))
    deep = n >= MT_MIN_DEPTH
    post[deep & (frac >= MT_HOMOPLASMIC), 2] = 1
    post[deep & (frac <= MT_CLEAR), 0] = 1
    het = deep & (frac > MT_CLEAR) & (frac < MT_HOMOPLASMIC)
    post[het, 1] = 1
    post[~deep] = [0.5, 0.25, 0.25]  # too few reads: uncertain
    return post

def imputed_gp(work, samples):
    """(chrom,pos,ref,alt) -> (n_kits x 3 GLIMPSE genotype probabilities, has-GP mask)."""
    path = os.path.join(work, 'cohort_imputed.tsv.gz')
    out = {}
    if not os.path.exists(path):
        return out
    with gzip.open(path, 'rt') as f:
        col = {x: i for i, x in enumerate(f.readline().rstrip('\n').split('\t'))}
        for line in f:
            p = line.rstrip('\n').split('\t')
            gp = np.zeros((len(samples), 3)); ok = np.zeros(len(samples), bool)
            for j, x in enumerate(samples):
                v = p[col[x]] if x in col else '.'
                if ':' not in v: continue
                try:
                    gp[j] = [float(y) for y in v.split(':')[1].split(',')][:3]; ok[j] = True
                except ValueError:
                    pass
            out[(p[0], p[1], p[2].upper(), p[3].upper())] = (gp, ok)
    return out

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

def hyper_b(s):
    # homopolymer indels: sporadic alt reads outside the documented breeds are slippage, not a background allele
    # (with B=400 the EM put AP3B1 cyclic neutropenia at 0.5% in every breed and left 96% of dogs uncertain)
    return B_HYPER_HOMOPOLYMER if s['kind'] == 'homopolymer' else B_HYPER[s['hyper']]

def call_of(p, hap=None):
    pc = p[:, 1] + p[:, 2]
    hi = np.where(p[:, 2] >= 0.9, 'two_copies', 'carrier')
    if hap is not None:
        hi = np.where(hap, 'hemizygous', hi)  # a male's single X copy carries it: affected for X-linked recessive
    return np.where(pc >= CALL_HI, hi, np.where(pc <= CALL_LO, 'clear', 'uncertain'))

# ── structural variants ───────────────────────────────────────────────────────────────────────────────────────
EPS_SV, K_EPS_SV = 0.01, 200.0   # breakpoint fragments: junction-read noise in non-carriers, shrunk to EPS_SV
DEPTH_MIN_LEN, FLANK_MAX = 2000, 50000
DEPTH_PCS = 3                    # shared depth components removed (GC / library biases that move many loci at once)
DEPTH_WINSOR = 0.35              # |log ratio| beyond this is left out of the component fit (a real CNV, or noise)
DEPTH_DF = 4                     # Student-t tails: one odd window must not make a carrier on its own

def depth_model(svs, ev, n):
    """Log depth ratio inside/flanks for every long DEL/DUP and every kit, with shared biases removed.

    Per site the raw ratio is centred on its cohort median. The top DEPTH_PCS components of the winsorised
    (site x kit) matrix are biases that move many loci together in one kit (GC, library, batch: the FGF3/4/19
    duplication region read 2-4x in whole batches of non-Ridgebacks). Each kit's loading is fitted by weighted
    least squares over its sites, excluding its own outliers, so a real event in one dog does not cancel itself.
    Returns {site_id: (residual log ratio, variance, expected log ratio of a homozygous deletion)}."""
    rows, sids = [], []
    for sv in svs:
        length = int(sv['end']) - int(sv['start'])
        if sv['type'] not in ('DEL', 'DUP') or length < DEPTH_MIN_LEN:
            continue
        e = ev[sv['site_id']]
        din, dfl = e[:, 3], e[:, 4]
        mu = dfl * length / (2 * min(length, FLANK_MAX))
        ok = np.isfinite(din) & np.isfinite(dfl) & (dfl > 0)
        lr = np.where(ok, np.log((np.nan_to_num(din) + 0.5) / (np.nan_to_num(mu) + 0.5)), np.nan)
        pv = np.where(ok, 1 / (np.nan_to_num(din) + 0.5) + 1 / (np.nan_to_num(dfl) + 0.5), np.nan)
        rows.append((lr, pv, np.where(ok, np.log(0.5 / (np.nan_to_num(mu) + 0.5)), np.nan)))
        sids.append(sv['site_id'])
    if not rows:
        return {}
    LR = np.array([r[0] for r in rows]); PV = np.array([r[1] for r in rows]); HOM = np.array([r[2] for r in rows])
    centre = np.nanmedian(LR, axis=1, keepdims=True)
    LR, HOM = LR - centre, HOM - centre
    ok = np.isfinite(LR)
    W = np.where(ok & (np.abs(LR) <= DEPTH_WINSOR), 1 / np.where(ok, PV + 0.01, 1), 0.0)
    X = np.where(W > 0, LR, 0.0)
    deep = (W > 0).mean(0) > 0.8
    k = min(DEPTH_PCS, len(sids) - 1)
    U = np.linalg.svd(X[:, deep] * np.sqrt(W[:, deep]), full_matrices=False)[0][:, :k] if k > 0 else np.zeros((len(sids), 0))
    R = np.full_like(LR, np.nan)
    for j in range(n):
        w = W[:, j]
        if k and (w > 0).sum() > k + 2:
            A = U * np.sqrt(w)[:, None]
            sc = np.linalg.lstsq(A, X[:, j] * np.sqrt(w), rcond=None)[0]
            R[:, j] = LR[:, j] - U @ sc
        else:
            R[:, j] = LR[:, j]
    mad = lambda x: 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))
    tau_site = np.array([max(0.002, mad(R[i]) ** 2 - np.nanmedian(PV[i])) for i in range(len(sids))])
    base = PV + tau_site[:, None]
    tau_kit = np.maximum(0.0, np.nanmedian(R ** 2, axis=0) / 0.455 - np.nanmedian(base, axis=0))
    V = base + np.nan_to_num(tau_kit)[None, :]
    return {sid: (R[i], V[i], HOM[i]) for i, sid in enumerate(sids)}

def sv_lik_factory(ka, kr, svtype, depth):
    """Per-dog likelihood over 0/1/2 copies of the SV allele from junction fragments and, for long DEL/DUP,
    the bias-corrected depth ratio (Student-t on the log scale)."""
    if depth is not None:
        lr, var, hom = depth
        use_depth = np.isfinite(lr) & np.isfinite(var)
        means = (np.stack([np.zeros_like(lr), np.full_like(lr, np.log(0.5)), np.nan_to_num(hom)], 1) if svtype == 'DEL'
                 else np.tile(np.log([1.0, 1.5, 2.0]), (len(lr), 1)))
    else:
        use_depth = np.zeros(len(ka), bool)
    def lik(eps, rows):
        L = read_lik(kr[rows], ka[rows], eps)
        if use_depth[rows].any():
            z2 = (np.nan_to_num(lr[rows])[:, None] - means[rows]) ** 2 / np.nan_to_num(var[rows], nan=1.0)[:, None]
            ld = (1 + z2 / DEPTH_DF) ** (-(DEPTH_DF + 1) / 2)
            ld = np.where(use_depth[rows][:, None], ld / np.maximum(ld.max(1, keepdims=True), 1e-300), 1.0)
            L = L * ld
        return L / np.maximum(L.max(1, keepdims=True), 1e-300)
    return lik, use_depth

def cmd_sv(work):
    D = load(work)
    sb = site_breeds(); bm = BreedMatcher(D['names'])
    col = {x: j for j, x in enumerate(D['samples'])}
    n = len(D['samples'])
    svs = list(csv.DictReader(open(os.path.join(work, 'sv_sites.tsv')), delimiter='\t'))
    ev = collections.defaultdict(lambda: np.full((n, 5), np.nan))
    for line in open(os.path.join(work, 'sv_counts.tsv')):
        p = line.rstrip('\n').split('\t')
        if len(p) < 8 or p[0] not in col: continue
        ev[p[1]][col[p[0]]] = [float(x) if x != '' else np.nan for x in p[2:7]]
    cidx = {c: i for i, c in enumerate(D['codes'])}
    depth = depth_model(svs, ev, n)
    p_all = np.zeros((len(svs), n), np.float32); calls_all = np.empty((len(svs), n), dtype='U11')
    with open(os.path.join(work, 'bp_sv_sites.tsv'), 'w') as fo:
        fo.write('site_id\tgene\tcondition\treport_class\ttype\tlength\tdocumented_breeds\tbackground_af\tbreakpoint_error\t'
                 'depth_used\ttop_breed_af\tkits_with_evidence\tcarrier\ttwo_copies\themizygous\tuncertain\tclear\theld\n')
        for i, sv in enumerate(svs):
            sid, c, t = sv['site_id'], sv['chrom'], sv['type']
            e = ev[sid]
            ka = np.nan_to_num(e[:, 0]) + np.nan_to_num(e[:, 2]); kr = np.nan_to_num(e[:, 1])
            length = int(sv['end']) - int(sv['start'])
            lik_fn, use_depth = sv_lik_factory(ka, kr, t, depth.get(sid))
            meta, names = sb.get(sid, ({}, []))
            doc, mixed = set(), False
            for nm in names:
                hit = bm.match(nm); mixed |= MIXED in hit; doc |= {h for h in hit if h != MIXED}
            mask = np.zeros(len(D['codes']), bool)
            for h in doc:
                if h in cidx: mask[cidx[h]] = True
            hap = hemizygous(D, c, sv['start'])
            informative = ((ka + kr) > 0) | use_depth
            f, beta, eps = fit_site(kr, ka, D['W'], D['F'], EPS_SV, mask,
                                    B_HYPER['mixed' if mixed else ('documented' if doc else 'none')],
                                    hap=hap, fit_eps=True, lik_fn=lik_fn, informative=informative)
            prior = geno_prior(D['W'] @ f, D['F'], hap)
            post = lik_fn(eps, np.arange(n)) * prior
            post /= post.sum(1, keepdims=True)
            p_all[i] = post[:, 1] + post[:, 2]
            cl = call_of(post, hap)
            if sid in HOLD:
                cl[:] = 'held'
            calls_all[i] = cl
            cnt = collections.Counter(cl); top = np.argsort(-f)[:3]
            fo.write('\t'.join(map(str, [sid, sv['gene'], meta.get('condition', ''), meta.get('report_class', ''), t, length,
                     ','.join(sorted(doc)), f'{beta:.2e}', f'{eps:.4f}', int(use_depth.sum()),
                     '; '.join(f"{D['names'].get(D['codes'][j], D['codes'][j])} {f[j]:.3f}" for j in top),
                     int(informative.sum()), cnt['carrier'], cnt['two_copies'], cnt['hemizygous'], cnt['uncertain'],
                     cnt['clear'], cnt['held']])) + '\n')
    np.savez_compressed(os.path.join(work, 'bp_sv.npz'), p=p_all, call=calls_all, samples=np.array(D['samples']),
                        sid=np.array([s['site_id'] for s in svs]))
    print(f'sv: {len(svs)} sites x {n} kits -> {work}/bp_sv_sites.tsv, bp_sv.npz')

# ── commands ──────────────────────────────────────────────────────────────────────────────────────────────────
def cmd_fit(work):
    D = load(work)
    sb = site_breeds()
    bm = BreedMatcher(D['names'])
    setup = site_setup(D, bm, sb)
    n = len(D['samples'])
    gps = imputed_gp(work, D['samples'])
    pc_all = np.zeros((len(D['sites']), n), np.float32)
    calls_all = np.empty((len(D['sites']), n), dtype='U11')
    with open(os.path.join(work, 'bp_sites.tsv'), 'w') as fo:
        fo.write('site_id\tchrom\tpos\tref\talt\tgene\tcondition\treport_class\tkind\tsource\tdocumented_breeds\t'
                 'prior\tbackground_af\tread_error\ttop_breed_af\tkits_with_reads\tcarrier\ttwo_copies\themizygous\tuncertain\tclear\t'
                 'held\n')
        for i, ((c, p, ref, alt), sid) in enumerate(zip(D['sites'], D['sid'])):
            s = setup[i]
            r, a = D['R'][i], D['A'][i]
            hap = hemizygous(D, c, p)
            f, beta, eps = fit_site(r, a, D['W'], D['F'], s['eps'], s['mask'], hyper_b(s), hap=hap,
                                    fit_eps=s['kind'] != 'snv')
            if c == 'chrM':
                post, source = mt_posterior(r, a), 'chrM allele fraction'
            else:
                post, source = posterior(r, a, D['W'], D['F'], eps, f, hap), 'reads + breed prior'
                g = gps.get((c, p, ref.upper(), alt.upper()))
                if g is not None and g[1].any():  # in the Dog10K panel: GLIMPSE posterior wherever the kit has one
                    post[g[1]] = g[0][g[1]]; source = 'GLIMPSE (Dog10K)'
            pc_all[i] = post[:, 1] + post[:, 2]
            cl = call_of(post, hap)
            if sid in HOLD:
                cl[:] = 'held'
            calls_all[i] = cl
            calls = collections.Counter(cl)
            top = np.argsort(-f)[:3]
            meta = sb.get(sid, ({}, []))[0] if sid else {}
            fo.write('\t'.join(map(str, [sid or '', c, p, ref, alt, meta.get('gene', ''), meta.get('condition', ''),
                     meta.get('report_class', ''), s['kind'], source, ','.join(s['doc']), s['hyper'], f'{beta:.2e}', f'{eps:.4f}',
                     '; '.join(f"{D['names'].get(D['codes'][t], D['codes'][t])} {f[t]:.3f}" for t in top),
                     int(((r + a) > 0).sum()), calls['carrier'], calls['two_copies'], calls['hemizygous'],
                     calls['uncertain'], calls['clear'], calls['held']])) + '\n')
    np.savez_compressed(os.path.join(work, 'bp_p_carrier.npz'), p=pc_all, call=calls_all,
                        samples=np.array(D['samples']), sid=np.array([x or '' for x in D['sid']]),
                        reads=(D['R'] + D['A']).astype(np.int32))
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
        f, _, eps = fit_site(r, a, D['W'], D['F'], s['eps'], s['mask'], hyper_b(s), fit_eps=s['kind'] != 'snv')
        post = posterior(r, a, D['W'], D['F'], eps, f)
        calls = call_of(post); ok = t >= 0
        rare = (t[ok] > 0).mean() < 0.02  # carrier frequency like a founder mutation
        for c, g in zip(calls[ok], t[ok]):
            rep['glimpse'][(c, 'truth_carrier' if g > 0 else 'truth_clear')] += 1
            if rare:
                rep['glimpse_rare_sites'][(c, 'truth_carrier' if g > 0 else 'truth_clear')] += 1
        cal += list(zip(post[ok, 1] + post[ok, 2], t[ok] > 0))
        # B. deep kits thinned to ~1x; frequencies refitted without them so they cannot inform their own prior
        f2, _, eps2 = fit_site(r, a, D['W'], D['F'], s['eps'], s['mask'], hyper_b(s), fit_rows=~deep,
                               fit_eps=s['kind'] != 'snv')
        keep = 1.0 / np.maximum(np.nan_to_num(D['depth'], nan=1.0), 1.0)
        rt = rng.binomial(r, keep); at = rng.binomial(a, keep)
        pt = posterior(rt, at, D['W'], D['F'], eps2, f2)
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
    {'fit': cmd_fit, 'validate': cmd_validate, 'sv': cmd_sv}[cmd](work)
