#!/bin/bash
# Panel v3: one kit's inputs for `panel_v3_breed_prior.py score` (a kit sequenced after the cohort fit). The same
# extraction the cohort run did with cluster/panel_v3_genotype.sh, panel_v3_chrM.sh and panel_v3_sv.sh, for one BAM:
#
#   KIT_DIR/reads.tsv.gz     ref,alt read counts (AD) at the frozen small-variant sites
#   KIT_DIR/imputed.tsv.gz   GLIMPSE GT:GP at the v3 sites that are in the Dog10K panel
#   KIT_DIR/chrM_counts.tsv  allele counts at the mitochondrial sites, MAPQ 0 allowed (NUMTs make every read MAPQ 0)
#   KIT_DIR/sv.tsv           structural-variant read evidence (analysis/omia/sv_evidence.py count)
#   KIT_DIR/kit.json         breed mix (breed_composition_raw), F_ROH, depth, inferred sex from the kit's report files
#
#   analysis/omia/panel_v3_kit.sh SAMPLE BAM IMPUTED_BCF PUB_DIR MODEL_DIR KIT_DIR
# MODEL_DIR is reference_panel/panel_v3 (the frozen fit and its site lists). Needs bcftools on PATH, FASTA (default
# canFam4_plusY.fa next to the repo) and DATA_PYTHON with pysam.
set -euo pipefail
SAMPLE="${1:?sample}"; BAM="${2:?bam}"; IMP="${3:?imputed bcf (may be missing)}"; PUB="${4:?pub dir}"
M="${5:?model dir}"; K="${6:?kit dir}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FASTA="${FASTA:-$ROOT/canFam4_plusY.fa}"
PY="${DATA_PYTHON:-python3}"
export OPENBLAS_NUM_THREADS=1
mkdir -p "$K"

# -B off (no BAQ) keeps single low-coverage reads usable; -m 1: one gapped read is enough to consider an indel;
# call -A keeps the catalogue ALT on hom-ref records; norm puts indels in the catalogue's left-normalised form.
{ printf 'CHROM\tPOS\tREF\tALT\t%s\n' "$SAMPLE"
  bcftools mpileup -f "$FASTA" -R "$M/v3_small.regions.tsv" -q 20 -Q 20 -B -m 1 -d 1000 -a FORMAT/AD,FORMAT/DP -Ou "$BAM" 2>/dev/null \
    | bcftools call -m -A -C alleles -T "$M/v3_small.alleles.tsv.gz" -Ou 2>/dev/null \
    | bcftools norm -f "$FASTA" -Ou 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n'
} | gzip > "$K/reads.tsv.gz.tmp"
mv "$K/reads.tsv.gz.tmp" "$K/reads.tsv.gz"

{ printf 'CHROM\tPOS\tREF\tALT\t%s\n' "$SAMPLE"
  if [[ -s "$IMP" ]]; then
    bcftools view -T "$M/d10k_match.pos" "$IMP" -Ou | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%GP]\n'
  fi
} | gzip > "$K/imputed.tsv.gz"

: > "$K/chrM_counts.tsv"
zcat "$M/v3_small.alleles.tsv.gz" | awk '$1=="chrM"' | while IFS=$'\t' read -r c p al; do
  ref=${al%,*}; alt=${al#*,}
  bcftools mpileup -f "$FASTA" -r "$c:$p-$p" -q 0 -Q 20 -B -m 1 -d 100000 -a FORMAT/AD "$BAM" 2>/dev/null \
    | bcftools query -f '%REF\t%ALT\t[%AD]\n' 2>/dev/null \
    | awk -v k="$SAMPLE" -v c="$c" -v p="$p" -v r="$ref" -v a="$alt" -F'\t' '
        $1==r { n=split($2,al,","); split($3,ad,","); x=0; for(i=1;i<=n;i++) if(al[i]==a) x=ad[i+1]; print k"\t"c"\t"p"\t"r"\t"a"\t"ad[1]"\t"x; found=1 }
        # an indel site with no indel record means no read carried it: every read at the base record is reference
        length($1)==1 && length(r)>1 { m=split($3,ad,","); d=0; for(i=1;i<=m;i++) d+=ad[i] }
        END { if(!found) print k"\t"c"\t"p"\t"r"\t"a"\t"d+0"\t0" }' >> "$K/chrM_counts.tsv"
done

"$PY" "$ROOT/analysis/omia/sv_evidence.py" count "$M/sv_sites.tsv" "$SAMPLE" "$BAM" "$FASTA" > "$K/sv.tsv.tmp"
mv "$K/sv.tsv.tmp" "$K/sv.tsv"

"$PY" - "$SAMPLE" "$PUB" "$K/kit.json" <<'PYEOF'
import json, os, sys
sample, pub, out = sys.argv[1:4]
def load(name):
    try:
        return json.load(open(os.path.join(pub, name)))
    except (OSError, ValueError):
        return {}
b = load('breed_result.json')
breeds = {r['breed']: r['proportion'] for r in b.get('breed_composition_raw') or [] if r.get('proportion')}
cnv = load('cnv50_result.json') or {}
json.dump({'sample': sample, 'breeds': breeds or None, 'f_roh': load('inbreeding_result.json').get('f_roh'),
           'depth': load('qc_result.json').get('genome_mean_depth'), 'sex': cnv.get('sex_inferred')},
          open(out, 'w'))
PYEOF
echo "panel v3 inputs for $SAMPLE -> $K"
