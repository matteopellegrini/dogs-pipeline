"""Depth calibration for the two inbreeding numbers (stage 12), from 19 dogs sequenced at 4-18x, each downsampled
to 0.1/0.25/0.5/1/2x and re-imputed with the v2 panel (analysis/imputation_panel/froh_calib.sbatch).

  genotype F : full-depth F ~ a + b * F_obs, per depth level
  F_ROH      : full-depth F_ROH(>=500 kb) ~ c0 + c1*short + c2*mid + c3*long, where short/mid/long are the
               fractions of the genome in bcftools-roh runs of 0.5-1 / 1-5 / >=5 Mb at that depth

Stage 12 interpolates the corrected value between levels in log10(depth), clamps below 0.1x, and fades the
correction out between 2x and IDENTITY_AT (no correction from there up).
  python3 analysis/inbreeding/fit_depth_calibration.py > reference_json/inbreeding_depth_calibration.json"""
import csv, json, os, sys
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
LEVELS = ["0.1", "0.25", "0.5", "1", "2"]
G = 2200.0
IDENTITY_AT = 4.0

def loo(X, y):
    e = []
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        c = np.linalg.lstsq(X[m], y[m], rcond=None)[0]
        e.append(X[i] @ c - y[i])
    return np.abs(e)

def fit(rows, key, feats):
    full = {r["dog"]: float(r[key]) for r in rows if r["level"].startswith("full")}
    out = {}
    for lev in LEVELS:
        R = [r for r in rows if r["level"] == lev and r["dog"] in full]
        y = np.array([full[r["dog"]] for r in R]); raw = np.array([float(r[key]) for r in R])
        X = np.array([[1.0] + feats(r) for r in R])
        c = np.linalg.lstsq(X, y, rcond=None)[0]; e = loo(X, y)
        out[lev] = {"coef": [round(float(v), 5) for v in c], "n": len(y),
                    "raw_mae": round(float(np.mean(np.abs(raw - y))), 4), "raw_max": round(float(np.max(np.abs(raw - y))), 4),
                    "loo_mae": round(float(e.mean()), 4), "loo_max": round(float(e.max()), 4)}
    return out

froh = list(csv.DictReader(open(f"{H}/froh_calib_scores.tsv"), delimiter="\t"))
gf = list(csv.DictReader(open(f"{H}/gf_calib_scores.tsv"), delimiter="\t"))
res = {
    "description": "Sequencing-depth correction for stage 12 inbreeding, fitted on 19 dogs (4-18x) downsampled to 0.1-2x",
    "levels": [float(l) for l in LEVELS], "identity_at": IDENTITY_AT,
    "genotype_f": {"features": ["1", "genotype_f"], "fit": fit(gf, "genotype_f", lambda r: [float(r["genotype_f"])])},
    "f_roh": {"features": ["1", "frac_0.5_1Mb", "frac_1_5Mb", "frac_5Mb_plus"],
              "fit": fit(froh, "f500", lambda r: [float(r["mb_0.5_1"]) / G, float(r["mb_1_5"]) / G, float(r["mb_5plus"]) / G])},
}
json.dump(res, sys.stdout, indent=1); print()
for k in ("genotype_f", "f_roh"):
    for lev, v in res[k]["fit"].items():
        print(f"{k:10} {lev:>4}x  MAE {v['raw_mae']:.4f} -> {v['loo_mae']:.4f}   max {v['raw_max']:.4f} -> {v['loo_max']:.4f}", file=sys.stderr)
