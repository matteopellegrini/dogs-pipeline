#!/usr/bin/env python3
"""
Paternal line: Y-chromosome haplogroup from a whole-genome BAM aligned to a
reference that carries KP081776.1 as contig "chrY" (canFam4_plusY.fa).

  python3 y_haplogroup.py <markdup.bam> <y_tree_oetjens2018.tsv> <out.json>
                          [--auto-depth X] [--sex male|female|unknown]

Method (Oetjens et al. 2018, BMC Genomics 19:350): 795 diagnostic Y-SNVs on
KP081776.1, each tagged with the set of branches that carry the DERIVED
allele. A dog's Y is haploid, so at every covered site the reads should show
one allele: the derived one if his haplogroup lies on that branch, else the
reference one. Each candidate haplogroup H is scored by the log-likelihood
of the observed bases under "expected allele = derived iff H in branch set"
with per-base error eps; the winner's posterior (uniform prior over the six
dog haplogroups) is the confidence. Sites are only informative if the
candidates disagree there, so the JSON also reports how many such sites had
reads. Y coverage relative to autosomal depth is reported as a second sex
signal (male ~0.5 of autosomal on the haploid Y, female ~0).

Haplogroup names follow Oetjens/Ding (HG1-3, HG27, HG6, HG9, HG8, HG23);
Embark uses a different label scheme (e.g. A1a) - mapping is a separate,
documented table once anchored on shared samples.
"""
import json, math, sys
import pysam

DOG_HG = ['HG1-3', 'HG27', 'HG6', 'HG9', 'HG8', 'HG23']
EPS = 0.02
MIN_BQ, MIN_MQ = 20, 30


def load_tree(path):
    sites = []
    for l in open(path):
        if l.startswith('#') or l.startswith('pos\t'): continue
        pos, branch, ref, der, coy = l.rstrip('\n').split('\t')
        carriers = {b.strip() for b in branch.split(';')}
        sites.append((int(pos), carriers, ref, der))
    return sites


def pileup_bases(bam, chrom, pos):
    counts = {}
    for col in bam.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=MIN_BQ,
                          min_mapping_quality=MIN_MQ, ignore_overlaps=True, ignore_orphans=True):
        if col.reference_pos != pos - 1: continue
        for r in col.pileups:
            if r.is_del or r.is_refskip: continue
            b = r.alignment.query_sequence[r.query_position].upper()
            counts[b] = counts.get(b, 0) + 1
    return counts


def main():
    args = sys.argv[1:]; bam_path, tree_path, out_path = args[:3]
    def opt(n, d): return args[args.index(n) + 1] if n in args else d
    auto_depth = float(opt('--auto-depth', 0) or 0); sex = opt('--sex', 'unknown')
    bam = pysam.AlignmentFile(bam_path, 'rb')
    if 'chrY' not in bam.references:
        json.dump({'status': 'no_chrY_in_reference', 'method': 'requires alignment to canFam4_plusY.fa'}, open(out_path, 'w'), indent=2)
        print('y_haplogroup: reference has no chrY; skipped'); return
    ylen = bam.get_reference_length('chrY')

    # Y coverage over the callable, non-amplicon part (the diagnostic-SNV span; the
    # ~1 Mb amplicon region beyond ~1.3 Mb is masked by Oetjens' callability).
    sites = load_tree(tree_path)
    span0, span1 = min(s[0] for s in sites), max(s[0] for s in sites)
    n_reads = bam.count('chrY', span0, span1, read_callback=lambda r: r.mapping_quality >= MIN_MQ and not r.is_duplicate and not r.is_secondary)
    # mean depth ~ reads * read length / span
    rl = 150
    for r in bam.fetch('chrY', span0, span1):
        rl = r.query_length or 150; break
    y_depth = n_reads * rl / max(1, span1 - span0)
    y_auto = (y_depth / auto_depth) if auto_depth else None

    obs = []
    for pos, carriers, ref, der in sites:
        c = pileup_bases(bam, 'chrY', pos)
        n_ref, n_der = c.get(ref, 0), c.get(der, 0)
        if n_ref + n_der == 0: continue
        obs.append((pos, carriers, ref, der, n_ref, n_der))
    ll = {}
    for H in DOG_HG:
        s = 0.0
        for pos, carriers, ref, der, n_ref, n_der in obs:
            exp_der = H in carriers
            good, bad = (n_der, n_ref) if exp_der else (n_ref, n_der)
            s += good * math.log(1 - EPS) + bad * math.log(EPS)
        ll[H] = s
    best = max(ll, key=ll.get)
    m = max(ll.values()); post = {H: math.exp(v - m) for H, v in ll.items()}; z = sum(post.values()); post = {H: v / z for H, v in post.items()}
    ranked = sorted(post.items(), key=lambda kv: -kv[1])
    # informative sites: where the top-2 candidates expect different alleles
    second = ranked[1][0] if len(ranked) > 1 else best
    informative = [o for o in obs if (best in o[1]) != (second in o[1])]
    support = sum(o[5] if best in o[1] else o[4] for o in informative)
    against = sum(o[4] if best in o[1] else o[5] for o in informative)

    call_ok = (sex != 'female') and len(obs) >= 20 and post[best] >= 0.95
    res = {
        'status': 'called' if call_ok else ('not_applicable_female' if sex == 'female' else 'insufficient_data'),
        'haplogroup': best if call_ok else None,
        'posterior': round(post[best], 4),
        'runner_up': second, 'runner_up_posterior': round(post[second], 4),
        'sites_covered': len(obs), 'sites_total': len(sites),
        'informative_sites_vs_runner_up': len(informative), 'reads_supporting': support, 'reads_against': against,
        'y_depth_x': round(y_depth, 3), 'y_to_autosomal': round(y_auto, 3) if y_auto is not None else None,
        'sex_from_x_coverage': sex,
        'sex_from_y': ('male' if (y_auto is not None and y_auto > 0.15) else ('female' if y_auto is not None else None)),
        'nomenclature': 'Oetjens et al. 2018 / Ding et al. 2012 haplogroups',
        'method': ('Haploid genotypes at 795 diagnostic Y-SNVs (KP081776.1 as chrY), scored against the six dog '
                   'haplogroup branch sets with per-base error 0.02; confidence = posterior with a uniform prior.'),
    }
    json.dump(res, open(out_path, 'w'), indent=2)
    print(f"y_haplogroup: {res['status']} {best} (post {post[best]:.3f}; sites {len(obs)}/{len(sites)}; "
          f"informative vs {second}: {len(informative)} sites, {support} reads for / {against} against; "
          f"Y depth {y_depth:.2f}x, Y/auto {y_auto if y_auto is None else round(y_auto,3)})")


if __name__ == '__main__':
    main()
