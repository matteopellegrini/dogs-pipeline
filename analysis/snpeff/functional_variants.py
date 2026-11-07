#!/usr/bin/env python3
"""Stage 10: snpEff ANN of a dog's imputed Dog10K genotypes -> functional_variants.json (HIGH rows + MODERATE by gene).
Shared by run_dog_pipeline.sh and the cohort rebuild (cluster/rerun-10-fv-list.sh).
  functional_variants.py <annotated.vcf.gz> <imputed.bcf> <pub dir> <snpEff database name>
ANN format (pipe-delimited per transcript, comma-separated per variant):
  ALT | effect | impact | gene_name | gene_id | feature_type | feature_id |
  biotype | rank | hgvsc | hgvsp | cdna_pos | cds_pos | aa_pos | distance | messages
"""
import gzip, re, subprocess, json, sys

ANN_VCF, BCF, PUB, SNPEFF_DB = sys.argv[1:5]

# Build AF lookup from imputed BCF (chr-prefixed keys); only sites where the dog carries an alt allele get a row,
# so only those are kept (all 26M panel sites would need >6 GB)
print("Building AF lookup...")
af_lookup = {}
result = subprocess.run(
    ['bcftools', 'query', '-i', 'GT="alt"', '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RAF\n', BCF],
    capture_output=True, text=True)
for line in result.stdout.strip().split('\n'):
    if not line: continue
    parts = line.split('\t')
    if len(parts) < 5: continue
    try:
        af_lookup[(parts[0], int(parts[1]), parts[3])] = float(parts[4].split(',')[0])
    except Exception:
        pass
print(f"  AF lookup: {len(af_lookup)} sites")

HIGH = []
HIGH_COUNTS = {}
MOD_BY_GENE = {}
seen = set()   # deduplicate: (chrom, pos, ref, alt, gene, impact)

with gzip.open(ANN_VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        cols = line.strip().split('\t')
        if len(cols) < 8: continue
        chrom, pos_s, _, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
        pos = int(pos_s)

        # Genotype
        gt_str = ''
        if len(cols) >= 10:
            fmt = cols[8].split(':')
            smp = cols[9].split(':')
            if 'GT' in fmt:
                gt_str = smp[fmt.index('GT')]
        alleles = re.split(r'[|/]', gt_str)
        if set(alleles) <= {'0', '.'}: continue
        zyg = ('hom_alt' if set(alleles) == {'1'} else
               'het'     if '0' in alleles and '1' in alleles else 'other')

        af = af_lookup.get((chrom, pos, alt))

        # Parse ANN field
        info = dict(x.split('=', 1) for x in cols[7].split(';') if '=' in x)
        ann_str = info.get('ANN', '')
        if not ann_str: continue

        best = {}   # gene → best (highest impact) annotation for this variant
        for ann in ann_str.split(','):
            fields = ann.split('|')
            if len(fields) < 4: continue
            ann_alt    = fields[0]
            effect     = fields[1]
            impact     = fields[2]
            gene       = fields[3]
            if not gene or impact not in ('HIGH', 'MODERATE'): continue
            if ann_alt != alt: continue   # skip if annotation is for a different ALT
            # keep highest-impact annotation per gene
            rank = {'HIGH': 0, 'MODERATE': 1}
            if gene not in best or rank[impact] < rank[best[gene]['impact']]:
                best[gene] = {'impact': impact, 'effect': effect}

        for gene, ann_data in best.items():
            impact = ann_data['impact']
            effect = ann_data['effect']
            key = (chrom, pos, ref, alt, gene, impact)
            if key in seen: continue
            seen.add(key)

            row = {'impact': impact, 'gene': gene,
                   'chr': chrom.replace('chr', ''), 'pos': str(pos),
                   'ref': ref, 'alt': alt, 'effect': effect, 'zygosity': zyg,
                   'af_dog10k': round(af, 6) if af is not None else None}

            if impact == 'HIGH':
                HIGH.append(row)
                base_effect = effect.split('&')[0]
                HIGH_COUNTS[base_effect] = HIGH_COUNTS.get(base_effect, 0) + 1
            else:
                MOD_BY_GENE.setdefault(gene, []).append(row)

def rare(thresh):
    return sum(1 for r in HIGH if r['zygosity'] == 'hom_alt'
               and r['af_dog10k'] is not None and r['af_dog10k'] < thresh)

h_hom     = [r for r in HIGH if r['zygosity'] == 'hom_alt']
h_het     = [r for r in HIGH if r['zygosity'] == 'het']
mod_total = sum(len(v) for v in MOD_BY_GENE.values())
mod_hom   = sum(sum(1 for r in v if r['zygosity'] == 'hom_alt') for v in MOD_BY_GENE.values())
mod_het   = sum(sum(1 for r in v if r['zygosity'] == 'het')     for v in MOD_BY_GENE.values())

high_sorted = sorted(HIGH, key=lambda r: (
    0 if r['zygosity'] == 'hom_alt' else 1,
    r['af_dog10k'] if r['af_dog10k'] is not None else 1.0))

mod_gene_list = []
for g, rows in sorted(MOD_BY_GENE.items(), key=lambda kv: -len(kv[1])):
    hom_rows = [r for r in rows if r['zygosity'] == 'hom_alt']
    het_rows = [r for r in rows if r['zygosity'] == 'het']
    effects = list(dict.fromkeys(r['effect'].split('&')[0] for r in rows))
    hom_afs = [r['af_dog10k'] for r in hom_rows if r['af_dog10k'] is not None]
    min_af = min(hom_afs) if hom_afs else None
    mod_gene_list.append({
        'gene': g,
        'n_moderate': len(rows),
        'hom_alt': len(hom_rows),
        'het': len(het_rows),
        'effects': effects,
        'min_af': round(min_af, 6) if min_af is not None else None,
    })

fv = {
    'summary': {
        'high_total':          len(HIGH),
        'high_hom_alt':        len(h_hom),
        'high_het':            len(h_het),
        'high_hom_rare_1pct':  rare(0.01),
        'high_hom_rare_5pct':  rare(0.05),
        'high_hom_rare_10pct': rare(0.10),
        'moderate_total':      mod_total,
        'moderate_hom_alt':    mod_hom,
        'moderate_het':        mod_het,
    },
    'high_effect_counts': HIGH_COUNTS,
    'high_variants':      high_sorted,
    'moderate_by_gene':   mod_gene_list,
    'source':   f'GLIMPSE2 Dog10K imputation + SnpEff ({SNPEFF_DB}), canonical transcripts',
    'af_note':  'AF = allele frequency in Dog10K reference panel',
}
with open(f'{PUB}/functional_variants.json', 'w') as f:
    json.dump(fv, f)
print(f"functional_variants.json: {len(HIGH)} HIGH, {mod_total} MODERATE variants")
