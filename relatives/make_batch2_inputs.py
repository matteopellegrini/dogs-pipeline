#!/usr/bin/env python3
# Build relatives/batch2_eligible.txt + batch2_meta.json for the v3 reference.
# Eligibility mirrors v2: depth >= 0.3x, imputed BCF present, not bisulfite.
import json, os, glob
D = "/u/project/pellegrini/gkislik/dogs"
elig, meta, skip = [], {}, {"quar": 0, "depth": 0, "nobcf": 0, "noqc": 0}
with open(f"{D}/sample_sheet.batch2.tsv") as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        name, work, pub = p[0], p[4], p[5]
        if os.path.isdir(f"{D}/results_quarantine_bisulfite/{os.path.basename(pub)}"):
            skip["quar"] += 1; continue
        try:
            q = json.load(open(f"{pub}/qc_result.json"))
        except Exception:
            skip["noqc"] += 1; continue
        if float(q.get("genome_mean_depth") or 0) < 0.3:
            skip["depth"] += 1; continue
        if not glob.glob(f"{work}/glimpse2/*_imputed_dog10k.bcf"):
            skip["nobcf"] += 1; continue
        m = {"id": name, "source": "prosperkits"}
        try:
            comp = json.load(open(f"{pub}/breed_result.json"))["breed_composition"]
            if comp:
                top = comp[0]
                m["label"] = top["breed_name"] + ("" if top["proportion"] >= 0.85 else " mix")
        except Exception: pass
        try:
            w = json.load(open(f"{pub}/prs_result.json"))["physical_traits"]["weight_kg"]["pred_kg"]
            m["weight_kg"] = round(float(w), 1)
        except Exception: pass
        try:
            s = json.load(open(f"{pub}/coat_color.json"))["summary"]
            base = str(s.get("predicted_base_color", ""))
            if s.get("overall_confidence") != "low" and 0 < len(base) < 40:
                m["coat"] = base.lower()
        except Exception: pass
        elig.append(name); meta[name] = m
with open(f"{D}/relatives/batch2_eligible.txt", "w") as f:
    f.write("\n".join(elig) + "\n")
json.dump(meta, open(f"{D}/relatives/batch2_meta.json", "w"))
print(f"eligible: {len(elig)}  skipped: {skip}")
labels = sum(1 for m in meta.values() if "label" in m)
coats = sum(1 for m in meta.values() if "coat" in m)
print(f"with breed label: {labels}, with coat: {coats}, with weight: {sum(1 for m in meta.values() if 'weight_kg' in m)}")
