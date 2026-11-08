#!/usr/bin/env python3
"""Lift the Ensembl ROS_Cfam_1.0 gene table to canFam4 and rebuild its derived JSONs.

The table the pipeline used as genes_canFam4.json was never canFam4: it is
Ensembl's ROS_Cfam_1.0 (Labrador) gene set (ENSCAFG00845 ids), and only ~9% of
Panel v3 small variants fall inside their gene with it as-is. The original is
kept as reference_json/genes_ROS_Cfam_1.0.json; this script lifts every span
with UCSC GenArk's GCF_014441545.1ToCanFam4 chain (source names NC_051805.1..
NC_051842.1 = chr1..38, NC_051843.1 = chrX).

Each gene is lifted by both ends, stepping inward up to 10% of its span when an
end falls in a chain gap; a gene is dropped when an end cannot be placed, when
the ends land on different chromosomes, or when the lifted length differs from
the original by more than 50% (same rules as analysis/cnv/build_cnv_catalogue.py).

Outputs (reference_json/):
  genes_canFam4.json   {chrom without 'chr': [{name, start, end}]}, sorted by start
  genes_1mb.json       {chrN: {Mb bin: [gene names overlapping the bin]}}
  cnv_genes.json       the legacy stage-6 fallback table, exons lifted too
  genes_canFam4.build_report.txt   counts, dropped genes, the panel-variant check

Run from the repo root:
  python3 analysis/genes/build_genes_canfam4.py [rosCfam1ToCanFam4.over.chain.gz] [panel_v3_sites.tsv]
"""
import json
import os
import sys
from collections import Counter, defaultdict

from pyliftover import LiftOver

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
REF = os.path.join(ROOT, "reference_json")
D = os.environ.get("D", "/Users/matteopellegrini/Downloads/dogs")
CHAIN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    D, "known_variants/omia_download_2026-11-01/chains/rosCfam1ToCanFam4.over.chain.gz")
PANEL = sys.argv[2] if len(sys.argv) > 2 else os.path.join(D, "known_variants/panel_v3/panel_v3_sites.tsv")

ROS_ACC = {f"chr{i}": f"NC_0518{4 + i:02d}.1" for i in range(1, 39)}
ROS_ACC["chrX"] = "NC_051843.1"
LO = LiftOver(CHAIN)


def lift(chrom, start, end):
    """Lift a 1-based closed ROS interval; returns (chrom, start, end) or (None, reason)."""
    src = ROS_ACC[chrom]
    span = end - start
    step = max(1, span // 100)

    def one(pos, direction):
        for k in range(11):
            p = pos + direction * k * step
            if not start <= p <= end:
                break
            hit = LO.convert_coordinate(src, p - 1)
            if hit:
                return hit[0][0], hit[0][1] + 1
        return None, None

    c1, s = one(start, +1)
    c2, e = one(end, -1)
    if c1 is None or c2 is None:
        return None, "unmapped"
    if c1 != c2:
        return None, "split"
    if s > e:
        s, e = e, s
    if span > 0 and abs((e - s) / span - 1) > 0.5:
        return None, "length"
    return (c1, s, e), None


def unnamed(name):
    return name.startswith("ENSCAF")


report = []
ros = json.load(open(os.path.join(REF, "genes_ROS_Cfam_1.0.json")))
out = defaultdict(list)
why = Counter()
dropped = []
moved = []
for c, items in ros.items():
    for g in items:
        hit, reason = lift("chr" + c, g["start"], max(g["end"], g["start"] + 1))
        if hit is None:
            why[reason] += 1
            dropped.append(f"{g['name']}\tchr{c}:{g['start']}-{g['end']}\t{reason}")
            continue
        lc, ls, le = hit
        if lc.startswith("chrUn"):   # the pipeline reports chr1-38 and X only
            why["unplaced"] += 1
            dropped.append(f"{g['name']}\tchr{c}:{g['start']}-{g['end']}\tunplaced {lc}")
            continue
        if lc != "chr" + c:
            moved.append(f"{g['name']}\tchr{c} -> {lc}:{ls}-{le}")
        out[lc.replace("chr", "")].append({"name": g["name"], "start": ls, "end": le})
genes = {c: sorted(v, key=lambda g: (g["start"], g["end"])) for c, v in sorted(out.items())}
n_kept = sum(len(v) for v in genes.values())
report.append(f"genes: {sum(len(v) for v in ros.values())} ROS_Cfam_1.0 -> {n_kept} canFam4, "
              f"dropped {sum(why.values())} ({dict(why)}), {len(moved)} moved to another chromosome")
json.dump(genes, open(os.path.join(REF, "genes_canFam4.json"), "w"), separators=(",", ":"))

# genes_1mb.json: every gene overlapping each 1 Mb bin (the rule the old file was built with)
bins = {}
for c, items in genes.items():
    by_bin = defaultdict(list)
    for g in items:
        for b in range((g["start"] - 1) // 1_000_000, (g["end"] - 1) // 1_000_000 + 1):
            by_bin[b].append(g["name"])
    bins["chr" + c] = {str(b): sorted(set(n), key=lambda x: (unnamed(x), x)) for b, n in sorted(by_bin.items())}
order = sorted(bins, key=lambda k: (k == "chrX", int(k[3:]) if k[3:].isdigit() else 99))
json.dump({k: bins[k] for k in order}, open(os.path.join(REF, "genes_1mb.json"), "w"), separators=(",", ":"))
report.append(f"genes_1mb: {sum(len(v) for v in bins.values())} bins")

# cnv_genes.json: legacy fallback with exons/CDS; lift every interval, keep the gene if its span lifts
legacy = json.load(open(os.path.join(REF, "cnv_genes_ROS_Cfam_1.0.json")))
cg = defaultdict(list)
n_leg = 0
for c, items in legacy.items():
    for g in items:
        hit, _ = lift("chr" + c, g["start"], g["end"])
        if hit is None:
            continue
        n_leg += 1
        h = dict(g, start=hit[1], end=hit[2])
        for key in ("exons", "cds"):
            iv = []
            for s, e in g.get(key, []):
                x, _ = lift("chr" + c, s, max(e, s + 1))
                if x and x[0] == hit[0]:
                    iv.append(sorted([x[1], x[2]]))
            h[key] = sorted(iv)
        cg[hit[0].replace("chr", "")].append(h)
json.dump(cg, open(os.path.join(REF, "cnv_genes.json"), "w"))
report.append(f"cnv_genes: {sum(len(v) for v in legacy.values())} -> {n_leg}")


# Validation: Panel v3 small variants inside their own gene, before vs after
def inside(table, chrom, pos, name):
    return any(g["name"] == name and g["start"] <= pos <= g["end"] for g in table.get(chrom, []))


n = before = after = 0
for line in open(PANEL):
    f = line.rstrip("\n").split("\t")
    if f[0] == "site_id" or f[10] not in ("snv", "indel") or not f[12].startswith(("live", "verified")):
        continue
    chrom, pos, gene = f[5].replace("chr", ""), int(f[6]), f[3]
    if not any(g["name"] == gene for c in genes for g in genes[c]):
        continue
    n += 1
    before += inside(ros, chrom, pos, gene)
    after += inside(genes, chrom, pos, gene)
report.append(f"panel v3 small variants whose gene is in the table: {n}; inside their gene: "
              f"ROS table {before} ({before / n:.1%}), lifted {after} ({after / n:.1%})")

report.append("\ndropped genes:")
report += sorted(dropped)
report.append("\ngenes on another chromosome in canFam4:")
report += sorted(moved)
open(os.path.join(REF, "genes_canFam4.build_report.txt"), "w").write("\n".join(report) + "\n")
print("\n".join(report[:5]))
