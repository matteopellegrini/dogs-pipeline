#!/usr/bin/env python3
"""Build reference_panel/oral_function_panel.json from the HUMAnN cohort scores.

The panel holds, per functional module, the sorted CPM values of the reference
dogs (clean kits: >=1M microbial reads, Enterobacterales <20%; for gingipains,
clean kits with >=2M reads, because detection of that enzyme depends on depth)
plus the module definitions and the rules the pipeline applies. A new sample's
percentile = share of reference values below its CPM (ties count half).

Usage (Hoffman): python3 build_oral_function_panel.py humann_cohort/function_scores.tsv out.json
"""
import csv, json, sys, datetime

scores, out = sys.argv[1], sys.argv[2]
MODULES = [
    {"key": "sulfur", "col": "sulfur_vsc", "name": "Breath sulfur compounds",
     "description": "Bacteria that release smelly sulfur gases.", "direction": "bad",
     "ecs": ["4.4.1.11", "4.4.1.1", "4.4.1.28", "4.4.1.15"], "min_reads": 1_000_000},
    {"key": "ging", "col": "gingipain", "name": "Gum-damaging enzymes",
     "description": "Gingipains, enzymes that break down gum tissue.", "direction": "bad",
     "ecs": ["3.4.22.37", "3.4.22.47"], "min_reads": 2_000_000},
    {"key": "nit", "col": "nitrate", "name": "Nitrate recycling",
     "description": "Bacteria linked to healthier mouths.", "direction": "good",
     "ecs": ["1.7.5.1", "1.9.6.1", "1.7.99.4"], "min_reads": 1_000_000},
    {"key": "amm", "col": "ammonia", "name": "Ammonia production",
     "description": "Neutralises acids in plaque.", "direction": "neutral",
     "ecs": ["3.5.1.5", "3.5.3.6"], "min_reads": 1_000_000},
]
rows = list(csv.DictReader(open(scores), delimiter="\t"))
def reads(r):
    try: return float(r["microbial_reads"])
    except ValueError: return float("nan")
clean = [r for r in rows if r["overgrowth_flag"] == "0" and float(r["enterobacterales_pct"]) < 5 and reads(r) >= 1e6]
panel = {
    "version": 1, "built": datetime.date.today().isoformat(),
    "source": "HUMAnN 3.9 cohort run (humann_cohort/), unmapped reads capped at 5M, EC regroup, CPM",
    "rules": {"min_microbial_reads": 1_000_000, "read_cap": 5_000_000,
              "overgrowth_enterobacterales_pct": 20.0, "high_percentile": 75, "low_percentile": 25},
    "reference_description": "Percentile compared with dogs we have tested",
    "modules": [],
}
for m in MODULES:
    ref = [r for r in clean if reads(r) >= m["min_reads"]]
    vals = sorted(round(float(r[f"{m['col']}_cpm"]), 3) for r in ref)
    entry = {k: m[k] for k in ("key", "name", "description", "direction", "ecs", "min_reads")}
    entry.update({"n_ref": len(vals), "ref_cpm_sorted": vals})   # Hoffman's login python is 3.6: no dict union
    panel["modules"].append(entry)
    print(f"{m['key']:7s} n_ref={len(vals)} median={vals[len(vals)//2]:.1f}")
json.dump(panel, open(out, "w"))
print("wrote", out)
