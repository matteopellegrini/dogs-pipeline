#!/bin/bash
# Re-run stage 8 (known-variant genotyping: Dog10K imputed + direct reads) on every published kit after the
# variant catalogue (reference_json/omia_variants.json) gains entries. No alignment, no imputation: reads the
# kept imputed BCF and markdup.bam (else sites.bam) in the kit's work dir. ~1-2 min per kit.
#
#   qsub -t 1-N -tc 50 cluster/rerun-8.sh <list> <tag>
#
# Rows: name, barcode, pub dir, work dir (same lists as the stage 9 rerun). The previous omia_result.json is kept
# once per tag as omia_result.pre-<tag>.json. Nothing is uploaded here; republish from a login node afterwards.
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_data=4G,h_rt=1:00:00
#$ -pe shared 2
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONIOENCODING=utf-8
D="${SGE_O_WORKDIR:?run via qsub from the repo root}"
LIST="${1:?usage: qsub -t 1-N cluster/rerun-8.sh <list> <tag>}"
TAG="${2:?tag required, e.g. cys3}"
cd "$D"
row=$(sed -n "${SGE_TASK_ID}p" "$LIST"); [ -n "$row" ] || exit 0
pk=$(echo "$row" | cut -f1); PUB=$(echo "$row" | cut -f3); WORK=$(echo "$row" | cut -f4)
bcf=$(compgen -G "$WORK/glimpse2/*_imputed_dog10k.bcf" | head -1)
[ -n "$bcf" ] || { echo "RERUN8-SKIP $pk (no bcf)"; exit 0; }
# The pipeline finds the BCF by sample name, case-sensitively (UCLA kits are DOGS-Gen-NNN on disk but dogs-gen-NNN
# in the lists), so run under the name the BCF was written with.
name=$(basename "$bcf" _imputed_dog10k.bcf)
[ -f "$PUB/omia_result.pre-$TAG.json" ] || cp "$PUB/omia_result.json" "$PUB/omia_result.pre-$TAG.json" 2>/dev/null || true
sheet=$(mktemp "$D/logs/sheet8.XXXXXX.tsv")
printf 'sample_id\tfastq_dir\tage\toutput_name\twork_dir\tpub_dir\tfrom_stage\tsex\tnotes\n' > "$sheet"
printf '%s\t%s\t\t%s\t%s\t%s\t8\t\trerun stage 8 (%s)\n' "$name" "$D/fastq/none" "$name" "$WORK" "$PUB" "$TAG" >> "$sheet"
# Snapshot in the repo root (the pipeline resolves its paths from its own location; edits mid-task poison bash).
snap=$(mktemp "$D/.rdp_snapshot8.XXXXXX.sh"); cp run_dog_pipeline.sh "$snap"
TO_STAGE=8 PUBLISH_RESULTS=0 bash "$snap" "$sheet" 2 8; rc=$?
rm -f "$sheet" "$snap"
# ADD-ONLY: keep every published call exactly as it was and append just the catalogue entries the previous result
# did not have. A full stage-8 rerun reproduces ProsperKits exactly but shifts some calls on older kits whose results
# came from earlier caller versions (seen on COSMO3 and DOGS-Gen-102, 2026-10-28), so no call is replaced unless
# named. Entries are matched by site (chrom, pos, gene); a kept entry takes the new catalogue's labels (breed, HGVS,
# reference) with its published call. Sites dropped from the catalogue are removed. REPLACE="chrom:pos:GENE ..."
# names sites whose call must be recomputed (an entry whose alleles were re-encoded).
if [ $rc -eq 0 ] && [ -f "$PUB/omia_result.pre-$TAG.json" ]; then
  REPLACE="${REPLACE:-}" "$D/envs/genomics/bin/python" - "$PUB/omia_result.pre-$TAG.json" "$PUB/omia_result.json" <<'PY' || rc=3
import json, os, sys
old = json.load(open(sys.argv[1])); new = json.load(open(sys.argv[2]))
site = lambda v: f"{v.get('chrom')}:{v.get('pos')}:{str(v.get('gene') or '').strip().upper()}"
full = lambda v: (site(v), v.get('hgvs_c') or '')
replace = set(os.environ.get('REPLACE', '').split())
is_call = lambda x: isinstance(x, dict) and 'zygosity' in x
old_by = {}
for v in old['variants']: old_by.setdefault(site(v), []).append(v)
out, kept, added, replaced = [], 0, 0, 0
for v in new['variants']:
    prev = old_by.get(site(v), [])
    match = next((o for o in prev if full(o) == full(v)), None) or (prev[0] if len(prev) == 1 else None)
    if match is None or site(v) in replace:
        out.append(v); added += match is None; replaced += match is not None
        continue
    calls = {k: x for k, x in match.items() if is_call(x)}
    out.append({**{k: x for k, x in v.items() if not is_call(x)}, **calls}); kept += 1
dropped = len(old['variants']) - kept - replaced
old['variants'] = out
if isinstance(old.get('summary'), dict): old['summary']['total_screened'] = len(out)
json.dump(old, open(sys.argv[2], 'w'))
print(f"ADD-ONLY: kept {kept} published calls, added {added}, replaced {replaced}, dropped {dropped}")
PY
fi
echo "RERUN8-DONE $pk rc=$rc"
