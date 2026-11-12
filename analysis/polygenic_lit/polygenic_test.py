#!/usr/bin/env python3
"""Test the literature polygenic models (~/Projects/dog literature/results/polygenic_snps_canfam4.csv, scorable rows
in weights.tsv) on the Dog Aging Project genotypes (7,828 dogs, GLIMPSE dosages) against owner-reported outcomes.

usage: polygenic_test.py WEIGHTS DAP_VCF DAP_DIR OUT_DIR

Outcomes available in the 2024 release we hold: follow-up (AFUS) new diagnoses with condition codes (hip dysplasia
1110, elbow dysplasia 1108), osteosarcoma (HLES + AFUS cancer tables). Not available: the HLES baseline condition
table (only the body-system flags), coat colour, hepatic copper, litter size. Baseline orthopedic flag = unknown
hip status, so those dogs are left out of hip-dysplasia controls (controls also need >= 1 follow-up survey).
Score = sum over SNPs of weight x effect-allele dosage (DS for ALT, 2-DS for REF).
"""
import csv, math, os, random, subprocess, sys
from collections import Counter, defaultdict

wts_f, vcf, dap, out = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
BCF = '/u/project/pellegrini/gkislik/dogs/envs/genomics/bin/bcftools'
random.seed(1)

W = list(csv.DictReader(open(wts_f), delimiter='\t'))
with open(f'{out}/regions.tsv', 'w') as f:
    for w in W:
        f.write(f"{w['chrom']}\t{w['pos']}\n")
samples = subprocess.run([BCF, 'query', '-l', vcf], capture_output=True, text=True).stdout.split()
ds = {}   # (chrom,pos,ref,alt) -> list of ALT dosage
q = subprocess.run([BCF, 'query', '-R', f'{out}/regions.tsv', '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n', vcf],
                   capture_output=True, text=True).stdout
for line in q.splitlines():
    c = line.split('\t')
    ds[(c[0], c[1], c[2], c[3])] = [float(x) if x != '.' else math.nan for x in c[4:]]

# ---- phenotypes
def rows(name):
    return csv.DictReader(open(f'{dap}/{name}', newline=''))
ov = {r['dog_id']: r for r in rows('DAP_2024_DogOverview_v1.1.csv')}
ortho0 = {}
for r in rows('DAP_2024_HLES_dog_owner_v1.0.csv'):
    ortho0[r['dog_id']] = r.get('hs_health_conditions_orthopedic', '')
hd, ed, afus_any = set(), set(), set()
for r in rows('DAP_2024_AFUS_health_condition_v1.0.csv'):
    afus_any.add(r['dog_id'])
    if r['afus_hs_new_condition'] == '1110': hd.add(r['dog_id'])
    if r['afus_hs_new_condition'] == '1108': ed.add(r['dog_id'])
osa = set()
for r in rows('DAP_2024_HLES_cancer_conditions_v1.0.csv'):
    if r.get('hs_cancer_types_osteosarcoma') in ('True', 'TRUE', '1'): osa.add(r['dog_id'])
for r in rows('DAP_2024_AFUS_cancer_conditions_v1.0.csv'):
    if r.get('afus_hs_new_cancer_types_osteosarcoma') in ('True', 'TRUE', '1'): osa.add(r['dog_id'])
afus_done = {d for d, r in ov.items() if any('Y' in (r.get(f'AFUS_Year{y}_Complete') or '') for y in range(1, 6))}  # 12-char Y/N per survey part

log = open(f'{out}/report.txt', 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); log.write(s + '\n')
say('genotyped DAP dogs:', len(samples), '| with overview:', sum(s in ov for s in samples))
say('ortho baseline flag values:', dict(Counter(ortho0.get(s, '') for s in samples).most_common(4)))
say('AFUS complete (any year):', sum(s in afus_done for s in samples), '| with any AFUS condition row:', sum(s in afus_any for s in samples))
say('hip dysplasia (AFUS new, 1110):', sum(s in hd for s in samples), '| elbow dysplasia 1108:', sum(s in ed for s in samples),
    '| osteosarcoma (HLES+AFUS):', sum(s in osa for s in samples))

# ---- scores
models = defaultdict(list)
for w in W: models[w['pmid']].append(w)
score = {}
for pm, ws in models.items():
    s = [0.0] * len(samples); used = []
    for w in ws:
        d = ds.get((w['chrom'], w['pos'], w['ref'], w['alt']))
        if d is None:
            say(f'  {pm} {w["snp"]} {w["chrom"]}:{w["pos"]} not in the DAP genotypes'); continue
        used.append(w)
        b = float(w['weight'])
        for i, x in enumerate(d):
            if not math.isnan(x): s[i] += b * (x if w['effect_is'] == 'ALT' else 2 - x)
    score[pm] = s
    say(f'model {pm}: {len(used)}/{len(ws)} SNPs scored')

def breed_of(sid):
    r = ov.get(sid) or {}
    return (r.get('Breed') or '').strip() if (r.get('Breed_Status') or '').lower().startswith('pure') else ''

def eaf(w, idx):
    d = ds.get((w['chrom'], w['pos'], w['ref'], w['alt']))
    if d is None: return math.nan
    v = [d[i] for i in idx if not math.isnan(d[i])]
    f = sum(v) / (2 * len(v)) if v else math.nan
    return f if w['effect_is'] == 'ALT' else 1 - f

def auc(pos, neg):
    allv = sorted([(x, 1) for x in pos] + [(x, 0) for x in neg])
    r, i, ranks = 1, 0, []
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]: j += 1
        for k in range(i, j + 1): ranks.append(((i + j) / 2 + 1, allv[k][1]))
        i = j + 1
    rp = sum(r for r, l in ranks if l)
    return (rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))

def test(label, pm, cases, controls):
    s = score[pm]
    pos = [s[i] for i in cases]; neg = [s[i] for i in controls]
    if len(pos) < 5 or len(neg) < 5:
        say(f'{label}: {len(pos)} cases / {len(neg)} controls, not tested'); return
    a = auc(pos, neg)
    bs = sorted(auc(random.choices(pos, k=len(pos)), random.choices(neg, k=len(neg))) for _ in range(1000))
    allv = pos + neg; ge = 0
    for _ in range(1000):
        random.shuffle(allv)
        ge += auc(allv[:len(pos)], allv[len(pos):]) >= a
    say(f'{label}: cases {len(pos)}, controls {len(neg)}, AUC {a:.3f} (95% CI {bs[25]:.3f}-{bs[974]:.3f}), '
        f'one-sided permutation p {(ge + 1) / 1001:.3f}')

idx = {s: i for i, s in enumerate(samples)}
by_breed = defaultdict(list)
for s in samples: by_breed[breed_of(s)].append(idx[s])
say('\npurebred counts (top):', [(b, len(v)) for b, v in sorted(by_breed.items(), key=lambda x: -len(x[1])) if b][:12])

def hd_sets(ix):
    cases = [i for i in ix if samples[i] in hd]
    ctrl = [i for i in ix if samples[i] not in hd and samples[i] in afus_done and ortho0.get(samples[i]) == '0']   # 0 = no orthopedic disorder at baseline
    return cases, ctrl

lab = [i for b, v in by_breed.items() if 'labrador' in b.lower() for i in v]
gsd = [i for b, v in by_breed.items() if 'german shepherd' in b.lower() for i in v]
say('\n=== Hip dysplasia, Labrador model (PMID 25874693, 7 SNPs, AUC 0.85 dev / 0.80 validation in the paper)')
for w in models['25874693']:
    say(f"  {w['snp']} {w['chrom']}:{w['pos']} effect {w['effect_is']} w {w['weight']}: effect-allele freq Labs {eaf(w, lab):.3f}, all DAP {eaf(w, range(len(samples))):.3f}")
test('  purebred Labradors', '25874693', *hd_sets(lab))
test('  all DAP dogs (breed not adjusted)', '25874693', *hd_sets(range(len(samples))))
say('\n=== Hip dysplasia, German Shepherd single-SNP ORs (PMID 24802516; effect allele assumed minor)')
for w in models['24802516']:
    say(f"  {w['snp']} {w['chrom']}:{w['pos']} effect {w['effect_is']} lnOR {w['weight']}: freq GSD {eaf(w, gsd):.3f}, all {eaf(w, range(len(samples))):.3f}")
test('  purebred German Shepherds', '24802516', *hd_sets(gsd))
test('  all DAP dogs', '24802516', *hd_sets(range(len(samples))))

say('\n=== Does the Labrador hip score track breed hip-dysplasia rates? (breeds with >= 40 genotyped purebreds)')
pts = []
for b, ix in by_breed.items():
    if not b or len(ix) < 40: continue
    c, k = hd_sets(ix)
    if len(c) + len(k) < 30: continue
    m = sum(score['25874693'][i] for i in ix) / len(ix)
    pts.append((b, m, len(c) / (len(c) + len(k)), len(ix)))
def spearman(x, y):
    rk = lambda v: {i: r for r, i in enumerate(sorted(range(len(v)), key=lambda j: v[j]))}
    rx, ry = rk(x), rk(y); n = len(x)
    return 1 - 6 * sum((rx[i] - ry[i]) ** 2 for i in range(n)) / (n * (n * n - 1))
if len(pts) >= 5:
    say(f'  {len(pts)} breeds, Spearman(mean score, hip-dysplasia rate) = {spearman([p[1] for p in pts], [p[2] for p in pts]):+.2f}')
    for p in sorted(pts, key=lambda p: -p[2]): say(f'    {p[0][:28]:28s} n {p[3]:4d}  score {p[1]:.2f}  HD rate {p[2]*100:.1f}%')

say('\n=== Osteosarcoma (PMID 30890123, 21 SNPs, AA-vs-BB ORs from per-breed stepwise models; relative score only)')
big = [idx[s] for s in samples if any(k in (ov.get(s) or {}).get('Breed_Size_Class_at_HLES', '') for k in ('Large', 'Giant'))]
say('  size classes:', dict(Counter((ov.get(s) or {}).get('Breed_Size_Class_at_HLES', '') for s in samples)))
osa_case = lambda ix: [i for i in ix if samples[i] in osa]
osa_ctrl = lambda ix: [i for i in ix if samples[i] not in osa]
test('  large + giant dogs', '30890123', osa_case(big), osa_ctrl(big))
test('  all DAP dogs', '30890123', osa_case(range(len(samples))), osa_ctrl(range(len(samples))))
for b in ('Rottweiler', 'Greyhound', 'Irish Wolfhound'):
    ix = [i for k, v in by_breed.items() if b.lower() in k.lower() for i in v]
    say(f'  {b}: {len(ix)} genotyped, {len(osa_case(ix))} osteosarcoma')

say('\n=== Breed confound check: scores z-scored within breed (purebreds of breeds with >= 20 genotyped), cases pooled')
def within_breed(pm, case_fn, ctrl_fn, label):
    zs = {}
    for b, ix in by_breed.items():
        if not b or len(ix) < 20: continue
        v = [score[pm][i] for i in ix]; m = sum(v) / len(v); sd = (sum((x - m) ** 2 for x in v) / len(v)) ** .5 or 1
        for i in ix: zs[i] = (score[pm][i] - m) / sd
    ix = list(zs)
    keep = score[pm]
    score[pm + 'z'] = [zs.get(i, 0.0) for i in range(len(samples))]
    test(label, pm + 'z', case_fn(ix), ctrl_fn(ix))
within_breed('25874693', lambda ix: hd_sets(ix)[0], lambda ix: hd_sets(ix)[1], '  Labrador hip model, within-breed, all purebreds')
within_breed('24802516', lambda ix: hd_sets(ix)[0], lambda ix: hd_sets(ix)[1], '  German Shepherd hip model, within-breed, all purebreds')
bigset = set(big)
within_breed('30890123', lambda ix: osa_case([i for i in ix if i in bigset]), lambda ix: osa_ctrl([i for i in ix if i in bigset]),
             '  osteosarcoma model, within-breed, large + giant purebreds')

say('\n=== Models with no outcome in the DAP release: effect-allele frequencies only')
for pm in ('34043658', '26747866', '31802524', '28158183', '36780177', '41701742'):
    for w in models.get(pm, []):
        say(f"  {pm} {w['snp'][:16]:16s} {w['gene'][:18]:18s} {w['chrom']}:{w['pos']} effect {w['effect_is']} w {w['weight']}: all DAP {eaf(w, range(len(samples))):.3f}")

models = {k: v for k, v in models.items() if not k.endswith('z')}
with open(f'{out}/scores.tsv', 'w') as f:
    f.write('dog_id\tbreed\t' + '\t'.join(models) + '\thd\tosa\n')
    for i, s in enumerate(samples):
        f.write(f'{s}\t{breed_of(s)}\t' + '\t'.join(f'{score[p][i]:.4f}' for p in models) + f'\t{int(s in hd)}\t{int(s in osa)}\n')
