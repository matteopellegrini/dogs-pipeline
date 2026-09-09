#!/usr/bin/env python3
# Build cluster/publish_batch2.txt (name<TAB>barcode) from the batch-2 sheet,
# excluding: quarantined bisulfite, depth < 0.2x, qc_status FAIL.
import json, os, sys
D = "/u/project/pellegrini/gkislik/dogs"
rows_ok, excl = [], {"quarantine": 0, "low_depth": [], "qc_fail": [], "no_qc": []}
with open(f"{D}/sample_sheet.batch2.tsv") as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        name, pub = p[0], p[5]
        if os.path.isdir(f"{D}/results_quarantine_bisulfite/{os.path.basename(pub)}"):
            excl["quarantine"] += 1; continue
        try:
            q = json.load(open(f"{pub}/qc_result.json"))
        except Exception:
            excl["no_qc"].append(name); continue
        depth = float(q.get("genome_mean_depth") or 0)
        status = str(q.get("qc_status", ""))
        if depth < 0.2: excl["low_depth"].append((name, depth)); continue
        if status == "FAIL": excl["qc_fail"].append((name, status, depth)); continue
        barcode = name[3:] if name.startswith("pk-") else name
        rows_ok.append((name, barcode))
with open(f"{D}/cluster/publish_batch2.txt", "w") as f:
    for n, b in rows_ok: f.write(f"{n}\t{b}\n")
print(f"publishable: {len(rows_ok)}")
print(f"excluded quarantine: {excl['quarantine']}")
print(f"excluded depth<0.2x: {len(excl['low_depth'])}")
print(f"excluded QC FAIL: {len(excl['qc_fail'])}")
print(f"excluded no qc_result: {len(excl['no_qc'])}")
for n, d in excl["low_depth"][:15]: print(f"  low {n} {d}x")
for t in excl["qc_fail"][:15]: print("  fail", *t)
for n in excl["no_qc"][:5]: print("  noqc", n)
