#!/usr/bin/env python3
"""Rewrite CNV gene labels in published result dirs from the lifted canFam4 gene table.

genes_canFam4.json used to be ROS_Cfam_1.0 coordinates (see
analysis/genes/build_genes_canfam4.py). This recomputes, per kit, exactly what
stages 6e and 10 now write, without rerunning the pipeline:
  cnv50_result.json  segments[].genes (first 100, table order) + n_genes
  cnv_homdel.json    regions[].disrupted_genes / disrupted_gene_details,
                     top-level disrupted_genes, summary.unique_genes
  cnv_genes.json     the genome-wide gene map (biotype by gene name, taken from
                     the kits' old snpEff-derived maps: a gene's biotype does not
                     depend on where the ROS models put it)
Each original is kept once as pre-genes/<name>.json: a subdirectory, because
publish-results.mjs uploads every top-level *.json and the old cnv_genes.json is
~5 MB. Idempotent: reruns read the backup, never a file this script rewrote.

Usage: backfill_gene_labels.py <genes_canFam4.json> <result_dir> [<result_dir> ...]
"""
import json
import os
import shutil
import sys
from collections import Counter

TAG = "pre-genes"
genes = json.load(open(sys.argv[1]))
dirs = sys.argv[2:]


def original(d, name):
    """Path of the untouched file, making the pre-genes/ backup on first sight."""
    cur = os.path.join(d, name + ".json")
    bak = os.path.join(d, TAG, name + ".json")
    if not os.path.exists(bak):
        os.makedirs(os.path.dirname(bak), exist_ok=True)
        if not os.path.exists(cur):
            return None
        shutil.copy2(cur, bak)
    return bak


def overlap(chrom, start, end):
    """Genes with start < end_seg and end > start_seg (stage 6e rule)."""
    return [g["name"] for g in genes.get(chrom.replace("chr", ""), [])
            if g["start"] < end and g["end"] > start]


biotype = {}
for d in dirs[:50]:
    p = os.path.join(d, TAG, "cnv_genes.json")
    p = p if os.path.exists(p) else os.path.join(d, "cnv_genes.json")
    try:
        for items in json.load(open(p)).values():
            for g in items:
                if g.get("name") and g.get("biotype"):
                    biotype.setdefault(g["name"], g["biotype"])
    except Exception:
        pass
gene_map = {c: [{"gene_id": g["name"], "name": g["name"], "start": g["start"], "end": g["end"],
                 "strand": "+", "biotype": biotype.get(g["name"], "unknown"), "exons": [], "cds": []}
                for g in items] for c, items in genes.items()}
print(f"biotypes known for {sum(1 for c in gene_map.values() for g in c if g['biotype'] != 'unknown')} "
      f"of {sum(len(v) for v in gene_map.values())} genes", flush=True)
by_name = {g["name"]: g for c in gene_map.values() for g in c}

st = Counter()
for d in dirs:
    st["dirs"] += 1
    src = original(d, "cnv50_result")
    if src:
        r = json.load(open(src))
        if r:
            for s in r.get("segments", []):
                hits = overlap(s["chrom"], s["start"], s["end"])
                st["seg_changed"] += set(hits[:100]) != set(s.get("genes", []))
                s["genes"] = hits[:100]
                s["n_genes"] = len(hits)
                st["segments"] += 1
        json.dump(r, open(os.path.join(d, "cnv50_result.json"), "w"), indent=2)
        st["cnv50"] += 1

    src = original(d, "cnv_homdel")
    if src:
        cnv = json.load(open(src))
        all_d = {}
        for r in cnv.get("regions", []):
            dets = []
            for c_g in genes.get(r["chrom"].replace("chr", ""), []):
                if c_g["end"] < r["start"] or c_g["start"] > r["end"]:
                    continue
                dets.append({"gene": c_g["name"], "biotype": by_name[c_g["name"]]["biotype"],
                             "chrom": r["chrom"], "start": c_g["start"], "end": c_g["end"],
                             "overlap": "full" if c_g["start"] >= r["start"] and c_g["end"] <= r["end"] else "partial"})
            st["homdel_changed"] += {x["gene"] for x in dets} != set(r.get("disrupted_genes", []))
            r["disrupted_genes"] = [x["gene"] for x in dets]
            r["disrupted_gene_details"] = dets
            for x in dets:
                all_d[x["gene"]] = x
        cnv["disrupted_genes"] = list(all_d.values())
        cnv.setdefault("summary", {})["unique_genes"] = len(all_d)
        json.dump(cnv, open(os.path.join(d, "cnv_homdel.json"), "w"), indent=2)
        st["homdel"] += 1

    if original(d, "cnv_genes"):
        json.dump(gene_map, open(os.path.join(d, "cnv_genes.json"), "w"))
        st["cnv_genes"] += 1
print(dict(st))
