#!/usr/bin/env python3
"""Per-breed allele frequencies from Donner et al. 2023 (PLoS Genet 19:e1010651, 1,054,293 dogs, S3) for the panel v3
sites they genotyped, as starting frequencies for panel_v3_breed_prior.py (PV3_PUB_PRIOR).

A site takes the paper's column only when the paper's variant is the site's: same position and alleles (either
strand) for a single-base change, an indel within 25 bp for an insertion/deletion. Sites the panel calls another way
(held, sequence caller) are left out.

  python3 published_priors.py S3.tsv catalogue_match.tsv bp_sites.tsv panel_v23_breeds.txt OUT.json
"""
import csv, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seq_allele_caller import CALLED, K_SHRINK, s3_breed_af, s3_rows

COMP = str.maketrans('ACGT', 'TGCA')
SKIP = set(CALLED) | {'cat338'}  # sequence-caller sites; cat338 = SLC7A9 type II-B, the paper has 0 carriers anyway
# checked by hand 2026-11-16: the paper's variant is this site's although position/alleles differ
SAME = {'cat061': 'NIPAL4 c.744del: the paper used the OMIA listing OMIA later corrected (g.52737379 -> g.52737279)',
        'cat281': "L2HGDH: the site types the first base of the paper's GAA>AAG double substitution"}


def same_variant(paper, pos4, site):
    """Whether the paper's variant (CanFam3 string, lifted canFam4 position) is this panel site."""
    d = abs(int(site['pos']) - int(pos4))
    m = re.search(r'([ACGT])>([ACGT])$', paper)
    snv = len(site['ref']) == 1 and len(site['alt']) == 1
    if m:
        return snv and d == 0 and ((m[1], m[2]) == (site['ref'], site['alt']) or
                                    (m[1].translate(COMP), m[2].translate(COMP)) == (site['ref'], site['alt']))
    return not snv and d <= 25


def main(s3, match, bp_sites, breeds, out_path):
    rows, col, bm = s3_rows(s3, breeds)
    sites = {r['site_id']: r for r in csv.DictReader(open(bp_sites), delimiter='\t') if r['site_id']}
    out, rejected = {}, []
    for r in csv.DictReader(open(match), delimiter='\t'):
        if r['match'] == 'none' or r['omia_variant_id'].strip() not in col:
            continue
        for sid in filter(None, r['site_ids'].split(',')):
            if sid in SKIP or sid not in sites:
                continue
            if sid not in SAME and not same_variant(r['paper_variant_canfam3'], r['canfam4_pos'].split(':')[1], sites[sid]):
                rejected.append(f"{sid} {r['gene']} paper {r['paper_variant_canfam3']} -> {r['canfam4_pos']}, "
                                f"site {sites[sid]['chrom']}:{sites[sid]['pos']} {sites[sid]['ref']}>{sites[sid]['alt']}")
                continue
            mixed, af = s3_breed_af(rows, col[r['omia_variant_id'].strip()], bm)
            out[sid] = {'omia_variant': r['omia_variant_id'].strip()[:40], 'gene': r['gene'], 'mixed_af': mixed,
                        'breed_af': af}
    json.dump({'source': 'Donner et al. 2023 PLoS Genet 19:e1010651, S3', 'k_shrink': K_SHRINK, 'sites': out},
              open(out_path, 'w'), indent=1)
    print(f'{len(out)} sites with published frequencies -> {out_path}; {len(rejected)} matched rows rejected:')
    for x in rejected:
        print('  ', x)


if __name__ == '__main__':
    main(*sys.argv[1:6])
