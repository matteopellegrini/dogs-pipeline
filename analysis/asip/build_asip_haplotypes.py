"""Build reference_json/asip_haplotypes.json: the ASIP promoter haplotypes of
Bannasch et al. 2021 (Nat Ecol Evol 5:1415, doi 10.1038/s41559-021-01524-x),
placed on canFam4, for the read-based A-locus caller (asip_locus.py).

Hair cycle promoter (HCP), one amplicon-length sequence per haplotype, all
starting at the same primer-side base (canFam4 chr24:23,891,223):
  HCP1  canFam3 (Boxer, dominant yellow) chr24:23,377,633-23,378,892
  HCP2  GenBank MT319114.1 (Saarloos Wolfdog)       agouti
  HCP3  GenBank MT319115.1 (Doberman)               loss of function (tan points)
  HCP4  GenBank MT319116.1 (Beagle; = canFam4, a German Shepherd)
The GenBank sequences end at different points, so 'background' holds the
flanking sequence every haplotype shares (canFam3 and canFam4 either side of
the amplicon, including what follows the HCP4 deletion junction); the caller
drops any k-mer found there, leaving only haplotype-diagnostic k-mers.
  HCP5  GenBank MT319117.1 (Chinook)                loss of function
Ventral promoter (VP): VP1 carries a 185-bp SINE (canFam3 g.23353288_23353472,
present in canFam4 at chr24:23,863,804); VP2 lacks it and one copy of its
target-site duplication.
Inputs: canFam3 chr24:23,290,001-23,420,000 (UCSC), canFam4 chr24:23,830,001-
23,940,000, and the four GenBank records."""
import json, sys
def fa(p):
    s, k = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'): k = l.split()[0][1:]; s[k] = ''
        elif l: s[k] += l.upper()
    return s
d = sys.argv[1] if len(sys.argv) > 1 else '.'
cf3 = next(iter(fa(f'{d}/cf3_asip.fa').values()))   # idx 0 = canFam3 chr24:23,290,001
cf4 = next(iter(fa(f'{d}/cf4_asip.fa').values()))   # idx 0 = canFam4 chr24:23,830,001
g = fa(f'{d}/hcp_genbank.fa')
start = cf4.find(g['MT319116.1'][:40])
assert start == 61222, start
a, b = 23353288 - 23290001, 23353472 - 23290001     # VP SINE, 0-based inclusive in the cf3 slice
F = 150
vp_left = cf4.find(cf3[a - 40:a])
vp_right = cf4.find(cf3[b + 1:b + 41], vp_left)
TSD = 14
assert cf3[a - TSD:a] == cf3[b + 1:b + 1 + TSD] == 'TGACGGCCTGTTTA'
assert vp_left > 0 and 0 < vp_right - vp_left < 400, (vp_left, vp_right)
out = {
    'source': 'Bannasch et al. 2021 Nat Ecol Evol, doi:10.1038/s41559-021-01524-x (Supplementary Tables 1, 5, 7)',
    'chrom': 'chr24',
    'hcp_region': [23830000 + start - 400, 23830000 + start + 1400],   # 0-based fetch window
    'vp_region': [23830000 + vp_left - 300, 23830000 + vp_left + 40 + 185 + 300],
    'hcp': {
        'HCP1': cf3[87632:88892],
        'HCP2': g['MT319114.1'],
        'HCP3': g['MT319115.1'],
        'HCP4': g['MT319116.1'],
        'HCP5': g['MT319117.1'],
    },
    'hcp_background': [cf3[85000:87632], cf3[88892:92500], cf4[59000:start], cf4[start + 690:start + 3500]],
    'hcp_class': {'HCP1': 'ay', 'HCP2': 'aw', 'HCP3': 'at', 'HCP4': 'at', 'HCP5': 'at'},
    # VP2 loses the SINE and one copy of its 14-bp target-site duplication
    # (TGACGGCCTGTTTA), as Cosmo's split reads show
    'vp': {'VP1': cf3[a - F:b + 1 + F], 'VP2': cf3[a - F:a] + cf3[b + 1 + TSD:b + 1 + TSD + F]},
    # the same two alleles built on canFam4 (a VP1 German Shepherd), so k-mers
    # survive the SNPs that separate the two reference genomes near the SINE
    'vp_cf4': {'VP1': cf4[vp_left + 40 - F:vp_right + F], 'VP2': cf4[vp_left + 40 - F:vp_left + 40] + cf4[vp_right + TSD:vp_right + TSD + F]},
    # 1-based canFam4 positions; ASIP is on the plus strand here
    'coding': {
        'A82S': {'pos': 23906172, 'ref': 'G', 'alt': 'T'},   # c.244G>T, dominant-yellow haplotype tag
        'R83H': {'pos': 23906176, 'ref': 'G', 'alt': 'A'},   # c.248G>A, same haplotype
        'R96C': {'pos': 23906214, 'ref': 'C', 'alt': 'T'},   # c.286C>T, recessive black (a)
    },
}
# ── Diagnostic k-mers (homopolymer-compressed; see asip_locus.py) ─────────
# A k-mer is diagnostic for a class only if (1) it occurs in that class's
# haplotypes and in no other class, (2) no other-class or shared-flank k-mer is
# within 2 mismatches of it (one sample per haplotype cannot tell a haplotype
# SNP from a population SNP, so only structural breakpoints count), and (3) it
# occurs nowhere else in canFam4 (SINE bodies have thousands of
# copies whose reads land here at low mapping quality).
sys.path.insert(0, d)
from asip_locus import hpc, rc, K
import numpy as np, subprocess

def kset(s): return {s[i:i + K] for i in range(len(s) - K + 1)}
ENC = np.frombuffer(b'ACGTN', dtype=np.uint8)
def arr(ks): return np.array([np.frombuffer(k.encode(), dtype=np.uint8) for k in ks]) if ks else np.zeros((0, K), np.uint8)

def diagnostic(haps, background, mask):
    """{kmer: [haplotypes containing it]} for k-mers that are absent from at least
    one haplotype, at least 3 mismatches from every k-mer of the haplotypes (and
    shared flanks) that lack it, and not in `mask`."""
    hk = {h: set().union(*(kset(hpc(x)) for x in ([s] if isinstance(s, str) else s))) for h, s in haps.items()}
    bg = set().union(*(kset(hpc(s)) for s in background)) if background else set()
    allk = set().union(*hk.values())
    out = {}
    for k in sorted(allk):
        sig = tuple(h for h in haps if k in hk[h])
        if len(sig) == len(haps): continue
        if k in bg or rc(k) in bg or k in mask or rc(k) in mask: continue
        out[k] = sig
    # SNP-only differences: drop a k-mer if a haplotype outside its signature (or
    # the shared flank) has a k-mer within 2 mismatches of it.
    keep = {}
    arrs = {h: arr(sorted(hk[h] | {rc(x) for x in hk[h]})) for h in haps}
    barr = arr(sorted(bg | {rc(x) for x in bg}))
    for k, sig in out.items():
        a1 = np.frombuffer(k.encode(), dtype=np.uint8)
        near = any(len(arrs[h]) and ((arrs[h] != a1).sum(1) <= 2).any() for h in haps if h not in sig)
        if not near and len(barr) and ((barr != a1).sum(1) <= 2).any(): near = True
        if not near: keep[k] = list(sig)
    return keep

def genome_hits(cands, exclude):
    """Candidate k-mers seen anywhere in the homopolymer-compressed genome outside
    `exclude` (chrom, start, end). One pass: `tr -s` compresses runs, `grep -F`
    matches all k-mers at once."""
    import tempfile
    want = sorted(set(cands) | {rc(k) for k in cands})
    ch, s0, e0 = exclude
    with tempfile.TemporaryDirectory() as td:
        open(f'{td}/regions', 'w').write(''.join(
            (f"{l.split()[0]}\n" if l.split()[0] != ch else f"{ch}:1-{s0}\n{ch}:{e0 + 1}-{l.split()[1]}\n")
            for l in open(FASTA + '.fai')))
        open(f'{td}/kmers', 'w').write('\n'.join(want) + '\n')
        # one line per sequence (headers become line breaks), uppercase, runs squeezed
        cmd = (f"samtools faidx {FASTA} -r {td}/regions | tr 'acgtn' 'ACGTN' | "
               f"awk '/^>/{{if (n++) print \"\"; next}} {{printf \"%s\", $0}} END{{print \"\"}}' | "
               f"tr -s 'ACGTN' | LC_ALL=C grep -F -o -f {td}/kmers | sort -u")
        out = subprocess.run(['bash', '-c', cmd], capture_output=True, text=True).stdout.split()
    seen = set(out)
    print('genome scan:', len(seen), 'candidate k-mers found elsewhere', flush=True)
    return seen

FASTA = sys.argv[3]
pre = diagnostic(out['hcp'], out['hcp_background'], set())
vpre = diagnostic({v: [out['vp'][v], out['vp_cf4'][v]] for v in out['vp']}, [], set())
mask = genome_hits(list(pre) + list(vpre), ('chr24', 23_850_000, 23_920_000))
out['hcp_kmers'] = {k: v for k, v in pre.items() if k not in mask and rc(k) not in mask}
out['vp_kmers'] = {k: v for k, v in vpre.items() if k not in mask and rc(k) not in mask}
from collections import Counter
print('HCP signatures', Counter('+'.join(v) for v in out['hcp_kmers'].values()))
print('VP signatures', Counter('+'.join(v) for v in out['vp_kmers'].values()))
json.dump(out, open(sys.argv[2] if len(sys.argv) > 2 else 'asip_haplotypes.json', 'w'), indent=1)
print({k: len(v) for k, v in out['hcp'].items()}, {k: len(v) for k, v in out['vp'].items()}, out['hcp_region'], out['vp_region'])
