"""Per-sample table for the ProsperKits genotype matrix.

usage: pk_geno_samples.py samples.txt RESULTS_DIR SHEET_BATCH1 SHEET_BATCH2
Writes TSV to stdout: sample, kit, batch, mean_depth, qc_status, sex.
"""
import json
import os
import sys


def sheet_ids(path):
    with open(path) as f:
        next(f)
        return {line.split("\t", 1)[0] for line in f if line.strip()}


def load(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


samples_txt, results, sheet1, sheet2 = sys.argv[1:5]
b1, b2 = sheet_ids(sheet1), sheet_ids(sheet2)
print("sample\tkit\tbatch\tmean_depth\tqc_status\tsex")
for s in open(samples_txt).read().split():
    qc = load(os.path.join(results, s, "qc_result.json"))
    y = load(os.path.join(results, s, "y_result.json"))
    batch = "1" if s in b1 else "2" if s in b2 else "NA"
    sex = y.get("sex_from_x_coverage") or "NA"
    print("\t".join(str(v) for v in (
        s, s[3:], batch, qc.get("genome_mean_depth", "NA"),
        qc.get("qc_status", "NA"), sex)))
