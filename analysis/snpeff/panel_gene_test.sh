#!/bin/bash
# Does a snpEff database put Panel v3 variants in their own gene? Panel v3 small sites (kind snv/indel, status
# live*/verified*) have canFam4-native positions and a curated gene. Take the Dog10K panel records within 1 kb of each
# site (reference_panel/dog10k_ann/dog10k_sites.ann.bcf) plus the panel variants themselves, annotate them with the
# database under test, and count the sites where a genic ANN (not upstream/downstream/intergenic) names the panel gene.
# The ANN already in dog10k_sites.ann.bcf is scored too ("before").
#   sbatch -A pi_matteop -p campus -c 2 --mem=12G -t 1:00:00 -J pv3gene analysis/snpeff/panel_gene_test.sh [DB]
# Output: analysis/snpeff/panel_gene_test/<DB>.sites.tsv (one row per site) + a summary on stdout.
#SBATCH -o logs/%x.%j.log
set -euo pipefail
D="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$D"
export OPENBLAS_NUM_THREADS=1
E="$D/envs/genomics"
JAVA="$(find "$E/lib/jvm/bin" -name java 2>/dev/null | head -1 || true)"; [[ -n "$JAVA" ]] || JAVA="$E/bin/java"
export PATH="$E/bin:$(dirname "$JAVA"):$PATH"
DB="${1:-canFam4_NCBI_AR106}"
ANN_REF="${ANN_REF:-$D/reference_panel/dog10k_ann/dog10k_sites.ann.bcf}"
SITES="$D/known_variants/panel_v3/panel_v3_sites.tsv"
O="$D/analysis/snpeff/panel_gene_test"; mkdir -p "$O"
T="$O/$DB.tmp"; mkdir -p "$T"

awk -F'\t' 'NR>1 && ($11=="snv"||$11=="indel") && $13~/^(live|verified)/ {s=$7-1000; if(s<1)s=1; print $6"\t"s"\t"$8+1000}' "$SITES" \
  | sort -k1,1V -k2,2n > "$T/regions.tsv"
bcftools view -R "$T/regions.tsv" -Ov "$ANN_REF" | bcftools norm -d exact -Oz -o "$T/before.vcf.gz"
bcftools annotate -x INFO/ANN,INFO/LOF,INFO/NMD -Ov "$T/before.vcf.gz" \
  | snpEff ann -canon -noStats -noLog "$DB" 2> "$T/near.snpeff.log" | bgzip -c > "$T/after.vcf.gz"
awk -F'\t' 'NR>1 && ($11=="snv"||$11=="indel") && $13~/^(live|verified)/ && $9!="" && $10!="" && $9!~/[^ACGTacgt]/ && $10!~/[^ACGTacgt]/' "$SITES" \
  | awk -F'\t' 'BEGIN{OFS="\t"; print "##fileformat=VCFv4.2"; print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print $6,$7,$1,toupper($9),toupper($10),".",".","."}' \
  | snpEff ann -canon -noStats -noLog "$DB" 2> "$T/panel.snpeff.log" > "$T/panel.ann.vcf"

"$E/bin/python" - "$SITES" "$T" "$O/$DB.sites.tsv" << 'PYEOF'
import gzip, sys, collections
sites_tsv, T, out = sys.argv[1:4]
NON_GENIC = {'upstream_gene_variant', 'downstream_gene_variant', 'intergenic_region', 'intragenic_variant'}

def records(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        for line in f:
            if line[0] == '#': continue
            c = line.rstrip('\n').split('\t')
            info = dict(x.split('=', 1) for x in c[7].split(';') if '=' in x)
            anns = [a.split('|') for a in info.get('ANN', '').split(',') if a]
            yield c[0], int(c[1]), c[2], anns

def index(path):
    by_chrom = collections.defaultdict(list)
    for chrom, pos, _id, anns in records(path):
        by_chrom[chrom].append((pos, anns))
    return by_chrom

def genic_genes(anns):
    return {a[3].upper() for a in anns if len(a) > 3 and a[3] and not set(a[1].split('&')) <= NON_GENIC}

before, after = index(f'{T}/before.vcf.gz'), index(f'{T}/after.vcf.gz')
panel = {sid: anns for _c, _p, sid, anns in records(f'{T}/panel.ann.vcf')}

rows, n = [], collections.Counter()
for line in open(sites_tsv).read().splitlines()[1:]:
    c = line.split('\t')
    if c[10] not in ('snv', 'indel') or not (c[12].startswith('live') or c[12].startswith('verified')): continue
    sid, gene, chrom, pos, end = c[0], c[3].upper(), c[5], int(c[6]), int(c[7])
    def near(idx):
        return [anns for p, anns in idx.get(chrom, []) if pos - 1000 <= p <= end + 1000]
    nb, na = near(before), near(after)
    gb = set().union(*[genic_genes(a) for a in nb]) if nb else set()
    ga = set().union(*[genic_genes(a) for a in na]) if na else set()
    pa = panel.get(sid)
    pg = genic_genes(pa) if pa else set()
    peff = (pa[0][1] + '|' + pa[0][2] + '|' + pa[0][3]) if pa else ''
    n['sites'] += 1
    if nb: n['with_dog10k'] += 1; n['before_named'] += gene in gb; n['after_named'] += gene in ga
    if pa is not None: n['panel_annotated'] += 1; n['panel_named'] += gene in pg
    rows.append([sid, gene, chrom, str(pos), str(len(nb)), str(gene in gb), str(gene in ga),
                 ','.join(sorted(ga))[:200], str(gene in pg) if pa is not None else 'NA', peff])
with open(out, 'w') as f:
    f.write('site_id\tpanel_gene\tchrom\tpos\tn_dog10k_1kb\tbefore_named\tafter_named\tafter_genic_genes\tpanel_variant_named\tpanel_variant_top_ann\n')
    f.writelines('\t'.join(r) + '\n' for r in rows)
print(f"sites {n['sites']}; with Dog10K records within 1 kb {n['with_dog10k']}: "
      f"panel gene named before {n['before_named']}, after {n['after_named']}; "
      f"panel variants annotated {n['panel_annotated']}, own gene named {n['panel_named']}")
PYEOF
grep -c WARNING_REF_DOES_NOT_MATCH_GENOME <(zcat "$T/after.vcf.gz") || true
echo "PV3GENE-DONE $DB"
