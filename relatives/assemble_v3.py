#!/usr/bin/env python3
# Append batch-2 columns to the relatives reference (v2 -> v3).
import json, gzip, numpy as np, os
os.chdir("/u/project/pellegrini/gkislik/dogs")
sites = {}
with open("relatives_ref/keys.tsv") as f:
    for i, l in enumerate(f):
        c = l.rstrip("\n").split("\t")
        sites[(c[0], c[1], c[2], c[3])] = i
n_sites = len(sites)
meta = json.load(gzip.open("relatives_ref/meta.json.gz", "rt"))
prior = len(meta["samples"])
have = {m["id"] for m in meta["samples"] if isinstance(m, dict) and "id" in m}
g = np.load("relatives_ref/geno.npy")
assert g.shape[1] == prior, f"geno {g.shape[1]} != meta {prior}"
elig = [l.strip() for l in open("relatives/batch2_eligible.txt") if l.strip()]
pm = json.load(open("relatives/batch2_meta.json"))
cols, added = [], []
for s in elig:
    if s in have:
        print("SKIP already-present:", s); continue
    p = f"relatives/batch2/{s}.ds.tsv.gz"
    if not os.path.exists(p):
        print("SKIP no extract:", s); continue
    col = np.full(n_sites, -1, dtype=np.int8)
    with gzip.open(p, "rt") as f:
        for l in f:
            c = l.rstrip("\n").split("\t")
            idx = sites.get((c[0], c[1], c[2], c[3]))
            if idx is None: continue
            try: col[idx] = int(round(float(c[4])))
            except Exception: pass
    if (col >= 0).sum() < 50000:
        print("SKIP sparse:", s, int((col >= 0).sum())); continue
    cols.append(col); added.append(pm.get(s, {"id": s, "source": "prosperkits"}))
G = np.concatenate([g, np.stack(cols, axis=1)], axis=1)
np.save("relatives_ref/geno_v3.npy", G)
meta["samples"] = meta["samples"] + added
meta["built_v3"] = "2026-09-08"
meta["note_v3"] = f"+{len(added)} ProsperKits batch-2 dogs (depth>=0.3x; predicted weights)"
json.dump(meta, gzip.open("relatives_ref/meta_v3.json.gz", "wt"))
os.replace("relatives_ref/geno_v3.npy", "relatives_ref/geno.npy")
os.replace("relatives_ref/meta_v3.json.gz", "relatives_ref/meta.json.gz")
print(f"assembled v3: {G.shape[1]} dogs x {n_sites} sites | was {prior}, added {len(added)}")
