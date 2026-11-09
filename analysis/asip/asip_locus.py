"""Read-based A locus (ASIP) caller for low-pass whole genomes.

Dog A-locus patterns are set by two promoter modules (Bannasch et al. 2021):
the hair cycle promoter (HCP) decides yellow vs agouti vs tan points, and the
ventral promoter (VP) decides black saddle vs black back among tan-pointed
dogs. Recessive black (a) is a coding change, ASIP p.Arg96Cys.

  HCP1 -> ay (sable / fawn; dominant yellow with VP1, shaded yellow with VP2)
  HCP2 -> aw (agouti, wolf sable)
  HCP3/4/5 -> at (tan points; loss of hair-cycle expression)
  R96C -> a (recessive black)          dominance: ay > aw > at > a

Reads near the HCP are matched to the five haplotype sequences by k-mers that
occur in exactly one of them; a likelihood over the 15 HCP diplotypes weights
each haplotype by how many read placements can see one of its unique k-mers.
Reads at ASIP p.Ala82Ser/p.Arg83His add evidence for HCP1 (the two missense
changes sit on the HCP1 haplotype in 98% of the paper's genotyped dogs, but
are not causal themselves).
usage: python3 asip_locus.py asip_haplotypes.json sample.bam [...]"""
import itertools, json, math, sys
import pysam

K = 25             # in homopolymer-compressed sequence (~35 bp of real sequence)
EPS = 0.03          # read-assignment error
TAG_EPS = 0.03      # A82S/R83H vs HCP1 discordance
_RC = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(_RC)[::-1]
def hpc(s):
    """Homopolymer-compress: the SINE poly-A tails vary in length between dogs and
    sequencing runs, and would otherwise pass for haplotype differences."""
    return ''.join(c for i, c in enumerate(s) if i == 0 or c != s[i - 1])

def compat(seq, kmers):
    """Haplotypes compatible with a (compressed) sequence: the intersection of the
    signatures of the diagnostic k-mers it contains; None if it has none."""
    sig = None
    for i in range(len(seq) - K + 1):
        h = kmers.get(seq[i:i + K])
        if h is None: continue
        sig = set(h) if sig is None else sig & set(h)
    return frozenset(sig) if sig is not None else None

def kmer_index(kmers):
    idx = dict(kmers)
    idx.update({rc(k): v for k, v in kmers.items()})
    return idx

def placement_sets(haps, kmers, L):
    """For each haplotype, how often a read of compressed length L placed on it
    shows each compatible set: P(set | haplotype) for the likelihood."""
    out = {}
    for h, s in haps.items():
        s = hpc(s); cnt = {}
        for i in range(0, max(len(s) - L, 0) + 1):
            c = compat(s[i:i + L], kmers)
            if c is not None: cnt[c] = cnt.get(c, 0) + 1
        out[h] = cnt
    return out

def read_sets(bam, chrom, region, kmers):
    sets, lens, seen = [], [], set()
    for r in bam.fetch(chrom, region[0], region[1]):
        if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail: continue
        key = (r.query_name, r.is_read1)
        if key in seen: continue
        seen.add(key)
        s = r.query_sequence
        if not s: continue
        lens.append(len(s))
        c = compat(hpc(s), kmers)
        if c: sets.append(c)
    L = int(sorted(lens)[len(lens) // 2]) if lens else 150
    return sets, L

def base_counts(bam, chrom, pos):
    names = {}
    for col in bam.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=20,
                          min_mapping_quality=20, ignore_overlaps=True, ignore_orphans=False):
        for p in col.pileups:
            if p.is_del or p.is_refskip or p.alignment.is_duplicate: continue
            names[p.alignment.query_name] = p.alignment.query_sequence[p.query_position].upper()
    return names

def set_lik(c, h, place, alpha=0.5):
    """P(read shows compatible set c | it came from haplotype h and is informative)."""
    n = place[h]; tot = sum(n.values())
    m = len(place)
    if tot == 0: return 1e-9
    return (n.get(c, 0) + (alpha if h in c else 0)) / (tot + alpha * m)

def diplotype_post(sets, place, tag_alt, tag_ref, classes):
    haps = list(place)
    N = {h: max(sum(place[h].values()), 1) for h in haps}
    post = {}
    for x, y in itertools.combinations_with_replacement(haps, 2):
        prior = 0.04 if x == y else 0.08          # uniform haplotype frequencies (0.2 each)
        ll = 0.0
        for c in sets:
            p = (N[x] * set_lik(c, x, place) + N[y] * set_lik(c, y, place)) / (N[x] + N[y])
            ll += math.log((1 - EPS) * p + EPS / 50)
        f1 = ((classes[x] == 'ay') + (classes[y] == 'ay')) / 2
        p_alt = (1 - TAG_EPS) * f1 + TAG_EPS * (1 - f1)
        ll += tag_alt * math.log(p_alt) + tag_ref * math.log(1 - p_alt)
        post[(x, y)] = math.log(prior) + ll
    m = max(post.values())
    z = sum(math.exp(v - m) for v in post.values())
    return {k: math.exp(v - m) / z for k, v in post.items()}

RANK = {'ay': 0, 'aw': 1, 'at': 2, 'a': 3}

def call_asip(bam_path, hap_json):
    H = json.load(open(hap_json)) if isinstance(hap_json, str) else hap_json
    chrom = H['chrom']
    bam = pysam.AlignmentFile(bam_path, 'rb')
    cls = H['hcp_class']
    hk = kmer_index(H['hcp_kmers'])
    sets, L = read_sets(bam, chrom, H['hcp_region'], hk)
    place = placement_sets(H['hcp'], hk, int(L * 0.8))   # compressed reads are ~80% as long
    # Coding sites
    cod = {k: base_counts(bam, chrom, v['pos']) for k, v in H['coding'].items()}
    tag_reads = {}
    for k in ('A82S', 'R83H'):
        for qn, b in cod[k].items():
            if b == H['coding'][k]['alt']: tag_reads[qn] = 'alt'
            elif b == H['coding'][k]['ref']: tag_reads.setdefault(qn, 'ref')
    tag_alt = sum(1 for v in tag_reads.values() if v == 'alt')
    tag_ref = sum(1 for v in tag_reads.values() if v == 'ref')
    r96 = cod['R96C']
    a_alt = sum(1 for b in r96.values() if b == H['coding']['R96C']['alt'])
    a_ref = sum(1 for b in r96.values() if b == H['coding']['R96C']['ref'])
    # VP
    vk = kmer_index(H['vp_kmers'])
    vsets, _ = read_sets(bam, chrom, H['vp_region'], vk)
    vhits = {v: sum(1 for c in vsets if c == frozenset([v])) for v in H['vp']}
    bam.close()

    post = diplotype_post(sets, place, tag_alt, tag_ref, cls)
    cpost = {}
    for (x, y), p in post.items():
        k = tuple(sorted((cls[x], cls[y]), key=RANK.get))
        cpost[k] = cpost.get(k, 0) + p
    best_dip = max(post, key=post.get)
    best_cls = max(cpost, key=cpost.get)
    p_ay = sum(p for k, p in cpost.items() if 'ay' in k)
    n_inf = len(sets) + tag_alt + tag_ref
    # Sable/fawn is dominant, so one stray read would flip the visible pattern:
    # an ay call needs two independent supporting observations (reads whose
    # k-mers fit only HCP1, or reads with the p.Ala82Ser/p.Arg83His alleles).
    # A Yorkshire Terrier was called ay from a single HCP1 read (2026-11-09).
    ay_support = sum(1 for c in sets if c == frozenset(['HCP1'])) + tag_alt
    if cpost[best_cls] >= 0.90 and n_inf >= 3 and ('ay' not in best_cls or ay_support >= 2):
        alleles, conf = list(best_cls), 'medium'
    elif p_ay >= 0.95 and ay_support >= 2:
        alleles, conf = ['ay', '?'], 'low'
    else:
        alleles, conf = ['?', '?'], 'low'
    # Recessive black: the coding change replaces the lower-ranked allele.
    n96 = a_alt + a_ref
    if n96 >= 3 and a_alt >= 2 and a_ref == 0:
        alleles, conf = ['a', 'a'], 'medium'
    elif a_alt >= 1 and a_ref >= 1:
        alleles = [alleles[0], 'a'] if alleles[0] != '?' else ['?', 'a']
    res = dict(
        alleles=alleles, confidence=conf, informative_reads=n_inf, ay_support=ay_support,
        hcp_reads={'+'.join(sorted(c)): sum(1 for x in sets if x == c) for c in set(sets)}, read_len=L,
        hcp_best=list(best_dip), hcp_best_post=round(post[best_dip], 3),
        class_post={'/'.join(k): round(v, 3) for k, v in sorted(cpost.items(), key=lambda kv: -kv[1])},
        tag_reads={'alt': tag_alt, 'ref': tag_ref},
        r96c_reads={'alt': a_alt, 'ref': a_ref},
        vp_reads=vhits,
    )
    res['interpretation'] = describe(res)
    return res

def describe(r):
    """Owner-facing sentence for coat_color.json, in the style of the other loci."""
    a = r['alleles']; n = r['informative_reads']
    nr = r['r96c_reads']; vp = r['vp_reads']
    if a == ['a', 'a']:
        return (f"Recessive black (a/a) — {nr['alt']} of {nr['alt'] + nr['ref']} reads carry ASIP p.Arg96Cys: "
                'solid black (or liver) wherever the A locus shows.')
    shape = ''
    if 'at' in a and 'ay' not in a and 'aw' not in a:
        if vp['VP2'] >= 2 and vp['VP1'] == 0:
            shape = ' The ventral promoter reads point to the black-back form (tan points kept into adulthood).'
        elif vp['VP1'] >= 2 and vp['VP2'] == 0:
            shape = ' The ventral promoter reads point to the black-saddle form (tan spreads as the dog matures).'
    names = {'ay': 'sable/fawn (ay)', 'aw': 'agouti / wolf sable (aw)', 'at': 'tan points (at)', 'a': 'recessive black (a)'}
    if a == ['?', '?']:
        return (f'Not determined: {n} read{"s" if n != 1 else ""} cross the ASIP promoter and coding variants that '
                'separate sable, agouti and tan points, too few to settle it at this sequencing depth.')
    if a[1] in ('?',):
        return (f'At least one {names[a[0]]} copy, from {n} informative reads across the ASIP hair-cycle promoter; '
                'the other copy cannot be read at this depth. Sable/fawn is dominant, so it is the pattern '
                'wherever the A locus shows.')
    if a[0] == '?':
        return ('One recessive black (a) copy is seen at ASIP p.Arg96Cys; the other copy cannot be read at this depth.')
    plain = {'ay': 'sable/fawn', 'aw': 'agouti / wolf sable', 'at': 'tan points', 'a': 'recessive black'}
    top = plain[a[0]] if a[0] == a[1] else f'{plain[a[0]]} over {plain[a[1]]}'
    return (f'{top[0].upper() + top[1:]} ({a[0]}/{a[1]}) — from {n} informative reads across the ASIP '
            f'hair-cycle promoter haplotypes (Bannasch et al. 2021) and the p.Ala82Ser/p.Arg83His sites.{shape}')

if __name__ == '__main__':
    for b in sys.argv[2:]:
        print(b, json.dumps(call_asip(b, sys.argv[1])))
