#!/usr/bin/env python3
"""dogs-app lib/panelHomLof.json: rare (<1% of Dog10K dogs) predicted loss-of-function variants (snpEff HIGH, first
ANN = canonical transcript, highest impact) in disease genes (GenCC asset dogs-app lib/geneDiseases.json plus the
Panel v3 genes), and the ones seen with two copies in a Dog10K reference dog (AC_Hom > 0). The app lists a two-copy
predicted disruption only when it is NOT in `hom` (lib/predictedVariants.ts).
  build_panel_hom_lof.py <work_dog10k_ann dir> <geneDiseases.json> <panel_v3_sites.tsv> <out.json> <snpEff db>
Reads <dir>/chr*.damaging.tsv from cluster/dog10k_annotate.sh (CHROM POS REF ALT AF GENE EFFECT IMPACT HGVS_P AC_HOM).
"""
import datetime, glob, json, sys

work, gene_diseases, panel_sites, out, db = sys.argv[1:6]
genes = set(json.load(open(gene_diseases))['genes'])
for line in open(panel_sites).read().splitlines()[1:]:
    genes.add(line.split('\t')[3].strip())
sites, hom = 0, []
for path in sorted(glob.glob(f'{work}/chr*.damaging.tsv')):
    for line in open(path):
        c = line.rstrip('\n').split('\t')
        if len(c) < 10 or c[7] != 'HIGH' or c[5] not in genes: continue
        try:
            af = float(c[4])
        except ValueError:
            continue
        if not 0 < af < 0.01: continue
        sites += 1
        if c[9] not in ('', '.', '0'): hom.append(f'{c[0]}:{c[1]}:{c[3]}')
json.dump({'meta': {'source': 'Dog10K panel (AutoAndXPAR.Dog10K.phased_plus_disease_rh.bcf), AC_Hom via bcftools +fill-tags',
                    'built': datetime.date.today().isoformat(),
                    'scope': f'rare (<1%) predicted loss-of-function variants (snpEff {db}, -canon) in GenCC/OMIA disease genes',
                    'n_sites': sites},
           'hom': sorted(hom)}, open(out, 'w'))
print(f'{out}: {sites} rare LoF sites in {len(genes)} disease genes, {len(hom)} seen homozygous in Dog10K')
