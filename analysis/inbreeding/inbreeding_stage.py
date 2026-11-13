"""Stage 12: inbreeding from the imputed BCF, corrected for sequencing depth.

  inbreeding_result.json             F_ROH: bcftools roh (HMM, --GTs-only 30, panel AF) runs >= 500 kb over 2,200 Mb
  inbreeding_froh_dog10k_result.json genotype F (observed vs expected homozygosity at Dog10K panel AF), placed in
                                     the Dog10K reference distribution (what the Inbreeding tab shows)

Low-pass imputation inflates genotype F and distorts the ROH length spectrum; both are corrected with
reference_json/inbreeding_depth_calibration.json (analysis/inbreeding/fit_depth_calibration.py) using the kit's mean
depth from qc_result.json. Raw values are kept as *_raw.

  python3 inbreeding_stage.py BCF PUB --af-file panel_af.tsv.gz --het-file dog10k_het.het --calibration cal.json
          [--depth X] [--bcftools path]"""
import argparse, gzip, json, math, os, subprocess
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("bcf"); ap.add_argument("pub")
ap.add_argument("--af-file", required=True); ap.add_argument("--het-file", required=True)
ap.add_argument("--calibration", required=True)
ap.add_argument("--depth", type=float, default=None)
ap.add_argument("--bcftools", default="bcftools")
a = ap.parse_args()
AUTOSOMES = [f"chr{i}" for i in range(1, 39)]
AUTO_MB = 2200.0
MIN_ROH = 500e3

depth = a.depth
if depth is None:
    try: depth = float(json.load(open(f"{a.pub}/qc_result.json"))["genome_mean_depth"])
    except Exception: depth = None
cal = json.load(open(a.calibration))

def depth_corrected(raw, predict):
    """predict(level_str) -> full-depth estimate from the level's fit; interpolated in log10(depth)."""
    if depth is None: return raw
    levels = cal["levels"]; top = cal["identity_at"]
    key = lambda l: f"{l:g}"
    if depth >= top: return raw
    if depth <= levels[0]: return predict(key(levels[0]))
    if depth >= levels[-1]:
        w = (math.log10(depth) - math.log10(levels[-1])) / (math.log10(top) - math.log10(levels[-1]))
        return (1 - w) * predict(key(levels[-1])) + w * raw
    for lo, hi in zip(levels, levels[1:]):
        if lo <= depth <= hi:
            w = (math.log10(depth) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
            return (1 - w) * predict(key(lo)) + w * predict(key(hi))
    return raw

def lin(coef, feats): return coef[0] + sum(c * f for c, f in zip(coef[1:], feats))

# ── F_ROH: bcftools roh ──────────────────────────────────────
print("Computing F_ROH with bcftools roh...")
r = subprocess.run(f"{a.bcftools} roh --GTs-only 30 --AF-tag RAF -O r -r {','.join(AUTOSOMES)} {a.bcf}",
                   shell=True, capture_output=True, text=True)
if r.returncode != 0: raise SystemExit(f"bcftools roh failed: {r.stderr[-500:]}")
segs = []
for l in r.stdout.splitlines():
    x = l.split("\t")
    if x[0] == "RG" and x[2] in AUTOSOMES and int(x[5]) >= MIN_ROH:
        segs.append({"chrom": x[2], "start": int(x[3]), "end": int(x[4]), "length_mb": round(int(x[5]) / 1e6, 3)})
segs.sort(key=lambda s: (int(s["chrom"][3:]), s["start"]))
mb = [s["length_mb"] for s in segs]
roh_total_mb = sum(mb)
froh_raw = roh_total_mb / AUTO_MB
classes = [sum(m for m in mb if m < 1) / AUTO_MB, sum(m for m in mb if 1 <= m < 5) / AUTO_MB, sum(m for m in mb if m >= 5) / AUTO_MB]
froh = depth_corrected(froh_raw, lambda lev: lin(cal["f_roh"]["fit"][lev]["coef"], classes))
froh = min(max(froh, 0.0), 1.0)
level = ('Very Low' if froh < 0.03125 else 'Low' if froh < 0.0625 else
         'Moderate' if froh < 0.125 else 'High' if froh < 0.25 else 'Very High')
print(f"F_ROH={froh:.4f} (raw {froh_raw:.4f}, depth {depth}x, {level}), {roh_total_mb:.1f} Mb in {len(segs)} runs")
inbreeding = {
    'f_roh': round(froh, 4), 'f_roh_pct': round(froh * 100, 2), 'f_roh_raw': round(froh_raw, 4),
    'roh_total_mb': round(roh_total_mb, 2), 'roh_n_segments': len(segs),
    'level': level.replace(' ', '_'),
    'autosomal_genome_mb': AUTO_MB,
    'depth_x': depth,
    'roh_segments': segs,
    'method': 'Runs of homozygosity >= 500 kb (bcftools roh HMM on imputed genotypes, reference-panel allele '
              'frequencies), corrected for sequencing depth with a calibration from deeply sequenced dogs downsampled '
              'to 0.1-2x. Segments listed are the uncorrected runs.'
}
json.dump(inbreeding, open(f'{a.pub}/inbreeding_result.json', 'w'), indent=2)
print("inbreeding_result.json written")

# ── Genotype F vs Dog10K distribution ────────────────────────
print("Computing genotype F vs Dog10K distribution...")
f_vals = []
with open(a.het_file) as fh:
    next(fh)
    for line in fh:
        parts = line.split()
        if len(parts) > 4: f_vals.append(float(parts[4]))
f_arr = np.array(f_vals)
proc = subprocess.Popen([a.bcftools, 'query', '-r', ','.join(AUTOSOMES), '-f', '%CHROM\t%POS\t[%GT]\n', a.bcf],
                        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
panel = gzip.open(a.af_file, 'rt')
def next_panel():
    line = panel.readline()
    if not line: return None, None, None
    p = line.rstrip().split('\t'); v = p[2].split(',')[0]
    return p[0], int(p[1]), None if v == '.' else float(v)
o_hom = 0; e_hom = 0.0; n_obs = 0
pc, pp, pa = next_panel()
for line in proc.stdout:
    c, pos, gt = line.rstrip().split('\t'); pos = int(pos); cn = int(c[3:])
    while pc is not None:
        if pc == c and pp == pos: break
        if int(pc[3:]) < cn or (pc == c and pp < pos): pc, pp, pa = next_panel()
        else: break
    if pc != c or pp != pos or pa is None: continue
    o_hom += gt in ('0|0', '1|1', '0/0', '1/1'); e_hom += pa * pa + (1 - pa) ** 2; n_obs += 1
    pc, pp, pa = next_panel()
proc.wait(); panel.close()
if n_obs == 0 or (n_obs - e_hom) == 0:
    print("WARNING: no sites matched for genotype F")
else:
    f_raw = (o_hom - e_hom) / (n_obs - e_hom)
    dog_F = depth_corrected(f_raw, lambda lev: lin(cal["genotype_f"]["fit"][lev]["coef"], [f_raw]))
    pct = float(np.mean(f_arr < dog_F) * 100)
    hist_counts, hist_edges = np.histogram(f_arr, bins=40, range=(float(f_arr.min()), float(f_arr.max()) + 0.001))
    corr = '' if depth is None or depth >= cal["identity_at"] else f', corrected for {depth:g}x sequencing depth'
    json.dump({
        'sample_froh': round(dog_F, 4), 'sample_froh_raw': round(f_raw, 4),
        'sample_percentile': round(pct, 1),
        'n_samples': len(f_arr),
        'ref_froh_mean': round(float(f_arr.mean()), 4),
        'ref_froh_p25': round(float(np.percentile(f_arr, 25)), 4),
        'ref_froh_p50': round(float(np.percentile(f_arr, 50)), 4),
        'ref_froh_p75': round(float(np.percentile(f_arr, 75)), 4),
        'hist_counts': hist_counts.tolist(), 'hist_edges': hist_edges.tolist(),
        'metric': 'genotype_F', 'depth_x': depth,
        'note': f'Genotype F from {n_obs/1e6:.1f}M autosomal SNPs vs Dog10K panel ({len(f_arr)} dogs){corr}.'
    }, open(f'{a.pub}/inbreeding_froh_dog10k_result.json', 'w'), indent=2)
    print(f"inbreeding_froh_dog10k_result.json: F={dog_F:.4f} (raw {f_raw:.4f}, {pct:.1f}th pct)")
