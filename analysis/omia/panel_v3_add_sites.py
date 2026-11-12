#!/usr/bin/env python3
"""Fold sites genotyped after the Panel v3 freeze into the main work dir, without regenotyping the old ones.

  python3 panel_v3_add_sites.py ADD_DIR WORK_DIR

ADD_DIR is a cluster/panel_v3_genotype.sh work dir run with PV3_WORK on the new sites only (2026-11-11 literature
sites: work_panel_v3_lit) holding
  cohort_reads.tsv.gz      CHROM POS REF ALT + one AD column per kit   (bcftools query of the merged reads BCFs)
  cohort_imputed.tsv.gz    CHROM POS REF ALT + one GT:GP column per kit (merged imputed BCFs, Dog10K-matched sites)
  v3_small.norm_ids.tsv    CHROM POS REF ALT site_id of the new left-normalised records
  d10k_match.pos           new sites with a matching Dog10K record
  glimpse_sites_af.tsv     CHROM POS REF ALT AC AN in the Dog10K panel
Rows are appended to the same files in WORK_DIR with WORK_DIR's kit columns (kits missing from ADD_DIR get '.').
A row already present (same CHROM POS REF ALT) is left alone, so a rerun adds nothing. The fit is per site, so a
refit afterwards leaves every old site's frequencies and calls unchanged.
"""
import gzip, os, sys


def key(line):
    return tuple(line.split('\t', 4)[:4])


def merge_matrix(add, work, name):
    src, dst = os.path.join(add, name), os.path.join(work, name)
    if not os.path.exists(src):
        return 0
    with gzip.open(dst, 'rt') as f:
        head = f.readline().rstrip('\n').split('\t')
        have = {key(l) for l in f}
    new = []
    with gzip.open(src, 'rt') as f:
        col = {x: i for i, x in enumerate(f.readline().rstrip('\n').split('\t'))}
        for line in f:
            p = line.rstrip('\n').split('\t')
            if tuple(p[:4]) in have:
                continue
            new.append('\t'.join(p[:4] + [p[col[k]] if k in col else '.' for k in head[4:]]) + '\n')
    if new:
        tmp = dst + '.tmp'
        with gzip.open(dst, 'rt') as f, gzip.open(tmp, 'wt') as o:
            for line in f:
                o.write(line)
            o.writelines(new)
        os.replace(tmp, dst)
    return len(new)


def merge_lines(add, work, name, k=lambda l: l):
    src, dst = os.path.join(add, name), os.path.join(work, name)
    if not os.path.exists(src):
        return 0
    have = {k(l.rstrip('\n')) for l in open(dst)} if os.path.exists(dst) else set()
    new = [l for l in open(src) if l.strip() and k(l.rstrip('\n')) not in have]
    with open(dst, 'a') as o:
        o.writelines(new)
    return len(new)


def main(add, work):
    four = lambda l: tuple(l.split('\t')[:4])
    print('cohort_reads', merge_matrix(add, work, 'cohort_reads.tsv.gz'),
          '| cohort_imputed', merge_matrix(add, work, 'cohort_imputed.tsv.gz'),
          '| norm_ids', merge_lines(add, work, 'v3_small.norm_ids.tsv', four),
          '| d10k_match', merge_lines(add, work, 'd10k_match.pos'),
          '| glimpse_af', merge_lines(add, work, 'glimpse_sites_af.tsv', four), 'rows added')


if __name__ == '__main__':
    main(*sys.argv[1:3])
