#!/usr/bin/env python3
"""Read-count truth for an SRA reference dog at the test-panel sites of one chromosome, from its full-depth gVCF
(sra-full-genotype.sh: bcftools mpileup -a AD,DP -g 5 | call -m -g 5). Writes the cosmo_truth pileup format
(CHROM POS REF ALTS AD) that score_panels.py reads: variant records give their AD; a hom-ref record or band
covering the site gives "DP,0" (MIN_DP for bands), so the >=10-read rule applies the same way.
usage: sra_truth.py GVCF CHROM RAF_TSV_GZ OUT_TSV_GZ"""
import gzip, subprocess, sys
gvcf, chrom, raf, out = sys.argv[1:5]
D = '/u/project/pellegrini/gkislik/dogs'
sites = {}
for line in gzip.open(raf, 'rt'):
    pos, ref, alt, ac, an = line.split()
    sites.setdefault(int(pos), (ref, alt))
q = subprocess.Popen([f'{D}/envs/genomics/bin/bcftools', 'query', '-r', chrom,
                      '-f', '%POS\t%REF\t%ALT\t%INFO/END\t%INFO/MIN_DP\t[%DP]\t[%AD]\n', gvcf],
                     stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
n = 0
with gzip.open(out, 'wt') as fo:
    for line in q.stdout:
        pos, ref, alts, end, mindp, dp, ad = line.rstrip('\n').split('\t')
        pos = int(pos)
        last = int(end) if end not in ('.', '') else pos
        alts = [a for a in alts.split(',') if a not in ('.', '<*>')]
        if alts:   # variant record: AD per allele
            if pos in sites and ad not in ('.', ''):
                fo.write(f'{chrom}\t{pos}\t{ref}\t{",".join(alts)}\t{ad.split(",")[0]},{",".join(ad.split(",")[1:1 + len(alts)])}\n'); n += 1
            continue
        depth = int(mindp) if mindp not in ('.', '') else (int(dp) if dp not in ('.', '') else 0)
        if depth < 10:
            continue
        for p in range(pos, last + 1):   # hom-ref record or band
            s = sites.get(p)
            if s:
                fo.write(f'{chrom}\t{p}\t{s[0]}\t{s[1]}\t{depth},0\n'); n += 1
print(chrom, 'truth rows', n)
