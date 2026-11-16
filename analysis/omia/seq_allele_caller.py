#!/usr/bin/env python3
"""Count reads carrying a variant's own sequence, for variants a pileup misses (long insertions, complex
deletion-insertions, homopolymer length changes) and for checking SNV sites against a published position.

Every primary, non-duplicate read placed near the site is tested, including unmapped reads placed at their mate's
position and the soft-clipped part of a read: a read that holds the alternate haplotype's junction sequence counts as
alt, one that holds the reference junction (the flanks on both sides of the event, unbroken) as ref. Both strands are
searched. SNV sites capture the base between fixed flanks and keep it only at base quality >= 20.

  python3 seq_allele_caller.py count LIST OUT.tsv      LIST: sample <tab> bam (one per line)
Output: sample, site_id, ref_reads, alt_reads, other_reads (an SNV base that is neither, or an A/C run of another
length), detail (the other run lengths or bases seen).

Genotypes for the report sites (CALLED) combine the reads with a prior from the dog's breed mix: the alternate allele
frequency per breed in Donner et al. 2023 (S3, 1,054,293 dogs), shrunk to the paper's mixed-breed frequency for
breeds with few dogs; breeds the paper does not list and the unassigned share take the mixed-breed frequency.
  python3 seq_allele_caller.py prior S3.tsv PANEL_BREEDS.txt OUT.json   (once; -> reference_panel/panel_v3/seq_prior.json)
  python3 seq_allele_caller.py apply PRIOR.json COUNTS.tsv LIST         LIST: sample <tab> pub_dir; rewrites the CALLED
                                                                         entries of pub_dir/omia_result.json in place
Calls follow panel v3: clear when P(>=1 copy) <= 0.2%, a copy at >= 95%, two copies at P(2) >= 90%, else uncertain.
"""
import csv, json, math, os, re, sys

# Patterns are on the + strand of canFam4 (checked against canFam4_plusY.fa, 2026-11-15).
SITES = {
    # RPGRIP1 cord1/CRD4: chr15:18457641 T > T + GGAAGCAACAGGATG + A(29) (Mellersh 2006); the A run varies in length
    # between dogs, so any run >= 10 counts. Ref junction = left flank, the 15 bp, then AGATC unbroken.
    'cat073': dict(chrom='chr15', pos=18457641, ref=[r'TTTCTTGGAAGCAACAGGATGAGATCA'],
                   alt=[r'TTTCTTGGAAGCAACAGGATGA{10}', r'A{10}GGAAGCAACAGGATGAGATCA']),
    # PNPLA1 ichthyosis: chr12:5599030 ACC > TACTACTA (c.1445_1447delinsTACTACTA)
    'cat147': dict(chrom='chr12', pos=5599030, ref=[r'CAGGGGAAACCCTACTGAAG'], alt=[r'CAGGGGAATACTACTACTACTGAAG']),
    # ABCA4 Stargardt: C7 -> C8 at chr6:55639738-55639744 (c.4176insC)
    'cat010': dict(chrom='chr6', pos=55639737, run=(r'CATTGT(C+)TTTGGTG', 7, 8)),
    # KLKB1 prekallikrein deficiency: chr16:45427398 A>T (catalogue), last base of an A6 run
    'cat446': dict(chrom='chr16', pos=45427398, snv=(r'GTAAGTAAAAA(.)CTGACAGCG', 'A', 'T')),
    # SLC7A9 cystinuria II-B p.A217T at the position Donner 2023 gives (canFam3 chr1:119209134 lifted); the live
    # catalogue entry cat338 sits at chr1:120056493, 2.8 kb away
    'slc7a9_donner': dict(chrom='chr1', pos=120053689, snv=(r'TTTCGAGGGC(.)CGAAGCTGTC', 'G', 'A')),
}
PAD = 400
# Report sites called here (catalogue site_id -> Donner 2023 S3 column = OMIA variant id, per-read error, alt-read
# fraction in a carrier). Error covers sequencing error plus cross-sample reads (index hopping) at a fixed signature;
# carrier allele balance measured in public genomes with >= 4 reads (2026-11-15): the cord1 reads with the A29 run
# are under-sampled (0.38), the others ~0.5.
CALLED = {
    'cat073': dict(omia='699', err=0.005, ab=0.38),
    'cat147': dict(omia='616', err=0.005, ab=0.5),
    'cat010': dict(omia='1050', err=0.01, ab=0.45),
    'cat446': dict(omia='74', err=0.005, ab=0.5),
}
CALL_LO, CALL_HI, HOM_HI = 0.002, 0.95, 0.90
K_SHRINK = 20  # dogs: breeds with few genotyped dogs lean on the mixed-breed frequency
# paper breed names the matcher cannot place, and panel breeds the paper lacks (closest listed breed stands in)
PAPER_NAMES = {'Bulldog (Standard)': {'ENGLISH_BULLDOG'}, 'Bulldog (American)': {'AMERICAN_BULLDOG'}}
PROXY = {'AMERICAN_PIT_BULL_TERRIER': 'AMERICAN_STAFFORDSHIRE_TERRIER'}
COMP = str.maketrans('ACGTN', 'TGCAN')

def rc(s):
    return s.translate(COMP)[::-1]

def classify(site, seq, qual):
    """'ref', 'alt', ('other', detail) or None for one read sequence (+ strand) and its qualities."""
    if 'snv' in site:
        pat, r, a = site['snv']
        m = re.search(pat, seq)
        if not m or qual[m.start(1)] < 20:
            return None
        b = m.group(1)
        return 'ref' if b == r else 'alt' if b == a else ('other', b)
    if 'run' in site:
        pat, nr, na = site['run']
        m = re.search(pat, seq)
        if not m:
            return None
        n = len(m.group(1))
        return 'ref' if n == nr else 'alt' if n == na else ('other', f'C{n}')
    if any(re.search(p, seq) for p in site['alt']):
        return 'alt'
    if any(re.search(p, seq) for p in site['ref']):
        return 'ref'
    return None

def count(bam_path):
    import pysam
    out = {}
    with pysam.AlignmentFile(bam_path) as bam:
        for sid, s in SITES.items():
            n = {'ref': 0, 'alt': 0, 'other': 0}
            other = {}
            seen = set()
            for rd in bam.fetch(s['chrom'], max(0, s['pos'] - PAD), s['pos'] + PAD):
                if rd.flag & 0xD00 or rd.query_sequence is None:  # secondary, duplicate, supplementary
                    continue
                key = (rd.query_name, rd.is_read1)
                if key in seen:
                    continue
                seen.add(key)
                seq, q = rd.query_sequence, rd.query_qualities or [30] * len(rd.query_sequence)
                c = classify(s, seq, q) or classify(s, rc(seq), list(q)[::-1])
                if c is None:
                    continue
                if isinstance(c, tuple):
                    n['other'] += 1
                    other[c[1]] = other.get(c[1], 0) + 1
                else:
                    n[c] += 1
            out[sid] = (n['ref'], n['alt'], n['other'], ','.join(f'{k}:{v}' for k, v in sorted(other.items())))
    return out

def s3_breed_af(rows, i, bm):
    """Donner S3 column i -> (mixed-breed AF, {panel breed code: AF}), each breed's AF shrunk towards the mixed-breed
    AF with K_SHRINK dogs."""
    from breed_names import MIXED
    per, n_row = {}, {}
    for r in rows[4:]:
        if len(r) <= i:
            continue
        if r[0] == 'N dogs genotyped':
            n_row = r
        elif r[0] == 'Alternative allele frequency' and r[i] not in ('', 'NA', 'N/A') and n_row[i] not in ('', 'NA', 'N/A'):
            per[r[1]] = (int(float(n_row[i])), float(r[i]))
    mixed = per.get('Mixed Breed', (0, 0.0))[1]
    # several paper names can map to one panel breed (Dachshund varieties): pool them by dogs genotyped; a name
    # covering several panel breeds ("Poodle") counts only where no specific name does
    pool = {True: {}, False: {}}
    for name, (n, f) in per.items():
        # "Poodle (Toy)" -> "Toy Poodle"
        hit = PAPER_NAMES.get(name) or bm.match(re.sub(r'^(.*?) \((.*)\)$', r'\2 \1', name)) or bm.match(name)
        if not hit or hit == {MIXED} or n == 0:
            continue
        for code in hit:
            a = pool[len(hit) == 1].setdefault(code, [0, 0.0])
            a[0] += n; a[1] += n * f
    af = {}
    for code in set(pool[True]) | set(pool[False]):
        n, nf = pool[True].get(code) or pool[False][code]
        af[code] = round((nf + K_SHRINK * mixed) / (n + K_SHRINK), 6)
    for code, proxy in PROXY.items():
        if code not in af and proxy in af:
            af[code] = af[proxy]
    return mixed, af

def s3_rows(s3_path, breeds_path):
    """(S3 rows, {column header: index}, breed matcher over the panel breed codes)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from breed_names import BreedMatcher
    codes = [l.strip() for l in open(breeds_path) if l.strip()]
    bm = BreedMatcher({c: c.replace('_', ' ').title() for c in codes})
    rows = [l.rstrip('\n').split('\t') for l in open(s3_path)]
    return rows, {rows[2][i].strip(): i for i in range(2, len(rows[2]))}, bm

def build_prior(s3_path, breeds_path, out_path):
    rows, col, bm = s3_rows(s3_path, breeds_path)
    out = {}
    for sid, c in CALLED.items():
        mixed, af = s3_breed_af(rows, col[c['omia']], bm)
        out[sid] = {'omia_variant': c['omia'], 'mixed_af': mixed, 'breed_af': af}
        print(f'{sid}: mixed {mixed:.4f}, {len(af)} breeds, top', sorted(af.items(), key=lambda x: -x[1])[:4])
    json.dump({'source': 'Donner et al. 2023 PLoS Genet 19:e1010651, S3', 'k_shrink': K_SHRINK, 'sites': out},
              open(out_path, 'w'), indent=1)

def dog_af(prior, breeds):
    """Expected alternate allele frequency for a dog's breed mix ({code: proportion})."""
    tot = sum(breeds.values()) if breeds else 0.0
    q = sum(p * prior['breed_af'].get(b, prior['mixed_af']) for b, p in (breeds or {}).items())
    return q + max(0.0, 1.0 - tot) * prior['mixed_af']

def posterior(q, n_ref, n_alt, err, ab):
    """P(1 copy), P(2 copies) given reads, HWE prior at allele frequency q."""
    pri = [(1 - q) ** 2, 2 * q * (1 - q), q * q]
    pa = [err, ab, 1 - err]
    ll = [n_alt * math.log(x) + n_ref * math.log(1 - x) for x in pa]
    m = max(ll)
    w = [pri[k] * math.exp(ll[k] - m) for k in range(3)]
    t = sum(w)
    return w[1] / t, w[2] / t

def call(q, n_ref, n_alt, err, ab):
    p1, p2 = posterior(q, n_ref, n_alt, err, ab)
    p = p1 + p2
    c = ('two_copies' if p2 >= HOM_HI else 'carrier' if p >= CALL_HI else 'clear' if p <= CALL_LO else 'uncertain')
    return c, p, p2

def apply(prior_path, counts_path, lst):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from panel_v3_apply import genotype, gt_key
    prior = json.load(open(prior_path))['sites']
    counts = {}
    for r in csv.DictReader(open(counts_path), delimiter='\t'):
        counts[(r['sample'], r['site_id'])] = (int(r['ref_reads']), int(r['alt_reads']))
    for line in open(lst):
        f = line.rstrip('\n').split('\t')
        if len(f) < 2:
            continue
        sample, pub = f[0], f[1]
        try:
            b = json.load(open(f'{pub}/breed_result.json'))
            breeds = {r['breed']: r['proportion'] for r in b.get('breed_composition_raw') or [] if r.get('proportion')}
        except (OSError, ValueError):
            breeds = {}
        path = f'{pub}/omia_result.json'
        res = json.load(open(path))
        done = []
        for v in res['variants']:
            sid, k = v.get('site_id'), gt_key(v)
            if sid not in CALLED or not k or (sample, sid) not in counts:
                continue
            n_ref, n_alt = counts[(sample, sid)]
            q = dog_af(prior[sid], breeds)
            c, p, p2 = call(q, n_ref, n_alt, CALLED[sid]['err'], CALLED[sid]['ab'])
            g = genotype(c, p, 'sequence' if n_ref + n_alt else 'breed_prior',
                         {'site_reads': n_ref + n_alt, 'alt_reads': n_alt, 'prior_af': round(q, 5), 'p_two': round(p2, 4)})
            v[k] = g
            done.append(f'{sid}:{c}:{p:.3f}:{n_alt}/{n_ref + n_alt}')
        json.dump(res, open(path, 'w'), indent=1)
        print(sample, ' '.join(done))

def main(lst, out_path):
    with open(out_path, 'w') as out:
        out.write('sample\tsite_id\tref_reads\talt_reads\tother_reads\tdetail\n')
        for line in open(lst):
            f = line.rstrip('\n').split('\t')
            if len(f) < 2:
                continue
            try:
                res = count(f[1])
            except (OSError, ValueError) as e:
                print(f'{f[0]}: {e}', file=sys.stderr)
                continue
            for sid, (r, a, o, d) in res.items():
                out.write(f'{f[0]}\t{sid}\t{r}\t{a}\t{o}\t{d}\n')
            out.flush()

if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'count':
        main(*sys.argv[2:4])
    elif cmd == 'prior':
        build_prior(*sys.argv[2:5])
    elif cmd == 'apply':
        apply(*sys.argv[2:5])
    else:
        sys.exit(__doc__)
