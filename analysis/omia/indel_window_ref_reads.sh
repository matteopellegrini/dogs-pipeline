#!/bin/bash
# Like indel_window_ref.sh, for indels the Dog10K call set under-calls (RPGRIP1 cord1: 5 carriers in 1,987 dogs, every
# American Staffordshire Terrier 0/0, where Donner 2023 finds 18%). The indel genotypes come from our own reads instead:
# deep public genomes genotyped by seq_allele_caller.py (counts TSV) join the reference with their window SNVs from
# sra_panel/<BS>.panel_sites.bcf. Dog10K dogs keep their window SNVs; their indel genotype is the Dog10K call when it is
# a carrier and missing otherwise, and Beagle 4.1 fills it in from the haplotypes it shares with the genotyped dogs.
#   analysis/omia/indel_window_ref_reads.sh NAME CHROM POS REF ALT SITE_ID COUNTS OUT_DIR [FLANK_BP]
# Public-genome genotype from COUNTS: >=2 alt reads and >=15% alt -> carrier (two copies at >=85% alt, >=6 reads);
# no alt read and >=6 reads -> clear; anything else missing.
set -euo pipefail
NAME=$1; C=$2; P=$3; REF=$4; ALT=$5; SID=$6; COUNTS=$7; O=$8; F=${9:-1000000}
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export PATH="$ROOT/envs/genomics/bin:$PATH" OPENBLAS_NUM_THREADS=1
PANEL="$ROOT/dog10k_panel/AutoAndXPAR.Dog10K.phased_plus_disease_rh.bcf"
INDELS="$ROOT/dog10k_indels/AutoAndXPAR.nonSNPs.filter.vcf.gz"
FASTA="$ROOT/canFam4_plusY.fa"
mkdir -p "$O"; T=$(mktemp -d "$O/.$NAME.XXXX")
bcftools query -l "$PANEL" > "$T/samples.txt"
lo=$(( P > F ? P - F : 1 )); hi=$(( P + F ))
bcftools view -r "$C:$lo-$hi" -S "$T/samples.txt" -m2 -M2 -v snps -Ou "$PANEL" | bcftools norm -d exact \
  | bcftools annotate -x INFO,^FORMAT/GT -Oz -o "$T/snv.vcf.gz"
bcftools index -t "$T/snv.vcf.gz"
# Dog10K carriers of the exact catalogue allele (any A-run length is the same insertion for cord1: take every ALT
# that starts with the catalogue ALT's first 17 bases)
bcftools view -r "$C:$P" -S "$T/samples.txt" "$INDELS" -Ou 2>/dev/null | bcftools norm -m -any -f "$FASTA" -Ou 2>/dev/null \
  | bcftools query -i "REF=\"$REF\"" -f '[%SAMPLE\t%ALT\t%GT\n]' | awk -v a="${ALT:0:17}" 'substr($2,1,17)==a && $3 ~ /1/ {print $1"\t0/1"}' \
  | sort -u > "$T/d10k_carriers.tsv"
python3 - "$T" "$COUNTS" "$SID" "$C" "$lo" "$hi" "$P" "$REF" "$ALT" "$ROOT" <<'PY'
import csv, gzip, os, subprocess, sys
T, counts, sid, chrom, lo, hi, pos, ref, alt, root = sys.argv[1:]
truth = {}
for r in csv.DictReader(open(counts), delimiter='\t'):
    s = r['sample']
    if r['site_id'] != sid or s.startswith(('pk-', 'dogs-gen', 'COSMO')) or not os.path.exists(f'{root}/sra_panel/{s}.panel_sites.bcf'):
        continue
    a, n = int(r['alt_reads']), int(r['ref_reads']) + int(r['alt_reads'])
    if a >= 2 and a / n >= 0.15:
        truth[s] = '1/1' if a / n >= 0.85 and n >= 6 else '0/1'
    elif a == 0 and n >= 6:
        truth[s] = '0/0'
pub = sorted(truth)
sites = []  # (pos, ref, alt) of the panel window
for line in subprocess.run(['bcftools', 'query', '-f', '%POS\t%REF\t%ALT\n', f'{T}/snv.vcf.gz'], capture_output=True, text=True).stdout.splitlines():
    p, r, a = line.split('\t'); sites.append((int(p), r, a))
gt = {p: {} for p, _, _ in sites}
for s in pub:
    out = subprocess.run(['bcftools', 'query', '-r', f'{chrom}:{lo}-{hi}', '-i', 'FMT/DP>=3', '-f', '%POS\t%REF\t%ALT\t[%GT]\n',
                          f'{root}/sra_panel/{s}.panel_sites.bcf'], capture_output=True, text=True).stdout
    for line in out.splitlines():
        p, r, a, g = line.split('\t')
        p = int(p)
        if p in gt:
            gt[p][s] = (r, a, g.replace('|', '/'))
d10k = dict(l.split() for l in open(f'{T}/d10k_carriers.tsv'))
panel = open(f'{T}/samples.txt').read().split()
with open(f'{T}/pub.vcf', 'w') as o:
    o.write('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    o.write(f'##contig=<ID={chrom}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(pub) + '\n')
    for p, r, a in sites:
        row = []
        for s in pub:
            x = gt[p].get(s)
            if not x or x[0] != r:
                row.append('./.')
            elif x[2] == '0/0':
                row.append('0/0')
            elif x[1] == a and x[2] in ('0/1', '1/0', '1/1'):
                row.append('0/1' if x[2] != '1/1' else '1/1')
            else:
                row.append('./.')
        o.write(f'{chrom}\t{p}\t.\t{r}\t{a}\t.\t.\t.\tGT\t' + '\t'.join(row) + '\n')
with open(f'{T}/indel.vcf', 'w') as o:
    o.write('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    o.write(f'##contig=<ID={chrom}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(panel + pub) + '\n')
    o.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t' + '\t'.join([d10k.get(s, './.') for s in panel] + [truth[s] for s in pub]) + '\n')
with open(f'{T}/truth.tsv', 'w') as o:
    for s in pub:
        o.write(f'{s}\t{truth[s]}\n')
print(f'public genomes: {len(pub)} ({sum(g != "0/0" for g in truth.values())} carriers); Dog10K carriers: {len(d10k)}', file=sys.stderr)
PY
bcftools view -Oz -o "$T/pub.vcf.gz" "$T/pub.vcf"; bcftools index -t "$T/pub.vcf.gz"
bcftools merge -m none -Oz -o "$T/snv_all.vcf.gz" "$T/snv.vcf.gz" "$T/pub.vcf.gz"; bcftools index -t "$T/snv_all.vcf.gz"
bcftools view -Oz -o "$T/indel.vcf.gz" "$T/indel.vcf"; bcftools index -t "$T/indel.vcf.gz"
bcftools concat -a -Oz -o "$T/merged.vcf.gz" "$T/snv_all.vcf.gz" "$T/indel.vcf.gz"
bcftools index -t "$T/merged.vcf.gz"
"$ROOT/envs/beagle/bin/beagle" -Xmx16g gt="$T/merged.vcf.gz" out="$T/phased" nthreads=${BEAGLE_THREADS:-4} > "$T/beagle.log" 2>&1
bcftools reheader --fai "$FASTA.fai" "$T/phased.vcf.gz" | bcftools +fill-tags -Ob -o "$O/$NAME.ref.bcf" -- -t AC,AN
bcftools index -f "$O/$NAME.ref.bcf"
cp "$T/truth.tsv" "$O/$NAME.truth.tsv"; cp "$T/d10k_carriers.tsv" "$O/$NAME.d10k_carriers.tsv"
echo "$NAME: $(bcftools view -H "$O/$NAME.ref.bcf" | wc -l) sites; indel AC/AN after phasing: $(bcftools query -i "POS=$P && REF=\"$REF\"" -f '%AC/%AN' "$O/$NAME.ref.bcf")"
rm -rf "$T"
