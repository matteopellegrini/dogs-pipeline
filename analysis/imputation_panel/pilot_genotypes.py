#!/usr/bin/env python3
"""Unphased genotypes of the pilot's added dogs at the panel v2 sites of one chromosome, as a VCF on stdout
(input to SHAPEIT5 phase_common with panel v2 as --reference).

usage: pilot_genotypes.py CHROM SITES_TSV_GZ DOGS_TSV > out.vcf
DOGS_TSV rows: name <TAB> kind <TAB> path, kind = gvcf (sra_full/<RUN>.g.vcf.gz, bcftools call -m -g 5) or
sites (sra_panel/<BS>.panel_sites.bcf, alleles already forced to the panel's).
A genotype is kept when the dog has >= MIN_DP reads there (gVCF bands: MIN_DP of the band); otherwise it is
missing and SHAPEIT5 imputes it. A gVCF variant whose ALT is not the panel ALT is missing too.
"""
import gzip, os, subprocess, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np

chrom, sites_tsv, dogs_tsv = sys.argv[1:4]
BCFTOOLS = '/u/project/pellegrini/gkislik/dogs/envs/genomics/bin/bcftools'
MIN_DP = int(os.environ.get('MIN_DP', '5'))

POS, REF, ALT = [], [], []
with gzip.open(sites_tsv, 'rt') as f:
    for line in f:
        c, p, al = line.rstrip('\n').split('\t')
        r, a = al.split(',', 1)
        POS.append(int(p)); REF.append(r); ALT.append(a)
POS = np.array(POS)
IDX = {p: i for i, p in enumerate(POS.tolist())}   # first record per position (panel v2 is biallelic, split)


def gt_code(gt):
    a = gt.replace('|', '/').split('/')
    if len(a) != 2 or '.' in a:
        return -1
    return int(a[0] != '0') + int(a[1] != '0') if all(x in '01' for x in a) else -1


def from_sites(path):
    g = np.full(len(POS), -1, np.int8)
    q = subprocess.Popen([BCFTOOLS, 'query', '-r', chrom, '-f', '%POS\t%ALT\t[%GT]\t[%DP]\n', path],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    for line in q.stdout:
        p, alt, gt, dp = line.rstrip('\n').split('\t')
        i = IDX.get(int(p))
        if i is None or (dp in ('.', '') or int(dp) < MIN_DP):
            continue
        if alt not in ('.', ALT[i]):
            continue
        g[i] = gt_code(gt)
    q.wait()
    return g


def from_gvcf(path):
    g = np.full(len(POS), -1, np.int8)
    q = subprocess.Popen([BCFTOOLS, 'query', '-r', chrom,
                          '-f', '%POS\t%REF\t%ALT\t%INFO/END\t%INFO/MIN_DP\t[%DP]\t[%GT]\n', path],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    for line in q.stdout:
        pos, ref, alts, end, mindp, dp, gt = line.rstrip('\n').split('\t')
        pos = int(pos)
        alts = [a for a in alts.split(',') if a not in ('.', '<*>')]
        if alts:   # variant record: map the dog's alleles onto the panel's REF/ALT
            i = IDX.get(pos)
            if i is None or dp in ('.', '') or int(dp) < MIN_DP or ref != REF[i]:
                continue
            a = gt.replace('|', '/').split('/')
            if len(a) != 2 or '.' in a:
                continue
            al = [ref] + alts
            calls = [al[int(x)] for x in a]
            if all(c in (REF[i], ALT[i]) for c in calls):
                g[i] = sum(c == ALT[i] for c in calls)
            continue
        depth = int(mindp) if mindp not in ('.', '') else (int(dp) if dp not in ('.', '') else 0)
        if depth < MIN_DP:
            continue
        last = int(end) if end not in ('.', '') else pos
        lo, hi = np.searchsorted(POS, pos), np.searchsorted(POS, last, side='right')
        sel = g[lo:hi]
        sel[sel == -1] = 0   # hom-ref band; a variant record at the same position (already set) wins
    q.wait()
    return g


def load(row):
    name, kind, path = row
    return name, (from_gvcf if kind == 'gvcf' else from_sites)(path)


dogs = [l.rstrip('\n').split('\t') for l in open(dogs_tsv) if l.strip() and not l.startswith('#')]
with ProcessPoolExecutor(int(os.environ.get('THREADS', '8'))) as ex:
    res = list(ex.map(load, dogs))
names = [n for n, _ in res]
G = np.stack([g for _, g in res], axis=1)
miss = (G == -1).mean(axis=0)
for n, m in zip(names, miss):
    print(f'{n}\tmissing {m:.3f}', file=sys.stderr)
code = {-1: './.', 0: '0/0', 1: '0/1', 2: '1/1'}
out = sys.stdout
out.write('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
out.write(f'##contig=<ID={chrom}>\n')
out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(names) + '\n')
for i in range(len(POS)):
    out.write(f'{chrom}\t{POS[i]}\t.\t{REF[i]}\t{ALT[i]}\t.\t.\t.\tGT\t' + '\t'.join(code[int(x)] for x in G[i]) + '\n')
