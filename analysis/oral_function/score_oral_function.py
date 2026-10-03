#!/usr/bin/env python3
"""Score one sample's bacterial-activity card (oral_function_result.json) from its
HUMAnN EC table and MetaPhlAn profile against reference_panel/oral_function_panel.json.

Usage:
  score_oral_function.py --ec <sample>_ec_cpm.tsv --profile <sample>_metaphlan.txt \
      --reads <microbial read count> --panel oral_function_panel.json --out oral_function_result.json

Rules (same as the cohort build, write_function_results.py):
  - reads < panel min_microbial_reads -> no card (exit 0, nothing written, prints SKIP)
  - Enterobacterales >= 20% of the community -> status 'overgrowth', no scores
  - each module: CPM = sum of its ECs (community rows), percentile vs the module's
    reference (ties count half); 'high' >= 75th, 'low' <= 25th; tone from direction
  - gingipain row needs >= its min_reads (2M), else 'not measurable'
  - note: vet (high sulfur or gingipain) > watch (any other concern) > typical
"""
import argparse, bisect, collections, json, re, sys, datetime

ap = argparse.ArgumentParser()
ap.add_argument("--ec", required=True); ap.add_argument("--profile", required=True)
ap.add_argument("--reads", type=int, required=True); ap.add_argument("--panel", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
panel = json.load(open(a.panel)); rules = panel["rules"]

if a.reads < rules["min_microbial_reads"]:
    print(f"SKIP: {a.reads} microbial reads < {rules['min_microbial_reads']}"); sys.exit(0)

# Enterobacterales share from the MetaPhlAn profile (order-level row).
entero = 0.0
for line in open(a.profile):
    if line.startswith("#"): continue
    p = line.split("\t")
    if p[0].split("|")[-1] == "o__Enterobacterales": entero = float(p[2])

ec2mod = {e: m["key"] for m in panel["modules"] for e in m["ecs"]}
total = collections.Counter(); carriers = collections.defaultdict(collections.Counter)
for line in open(a.ec):
    if line.startswith("#"): continue
    k, v = line.rstrip("\n").split("\t")[:2]
    if "|" in k:
        e, sp = k.split("|", 1)
        if e in ec2mod:
            name = re.sub(r"^g__[^.]*\.s__", "", sp).replace("_", " ")
            if name.lower() != "unclassified": carriers[ec2mod[e]][name] += float(v)
    elif k in ec2mod:
        total[ec2mod[k]] += float(v)

base = {"version": 1, "generated": datetime.date.today().isoformat(),
        "source": "HUMAnN 3.9 on unmapped reads (capped at 5M)", "microbial_reads": a.reads,
        "reference": panel["reference_description"],
        "caveat": "These results reflect the bacteria on the swab on the day it was taken and can change between swabs."}

if entero >= rules["overgrowth_enterobacterales_pct"]:
    doc = dict(base, status="overgrowth", enterobacterales_pct=round(entero, 1),
               message="This sample was dominated by environmental bacteria that grew after collection, so bacterial-activity scores are not shown.")
    json.dump(doc, open(a.out, "w"), indent=1); print(f"overgrowth: Enterobacterales {entero:.1f}%"); sys.exit(0)

def percentile(v, ref):
    return 100.0 * (bisect.bisect_left(ref, v) + 0.5 * (bisect.bisect_right(ref, v) - bisect.bisect_left(ref, v))) / len(ref)

mods = []
for m in panel["modules"]:
    key = m["key"]
    if a.reads < m["min_reads"]:
        mods.append({"key": key, "name": m["name"], "description": m["description"], "direction": m["direction"],
                     "measurable": False, "message": "Not measurable at this sample's sequencing depth.", "carriers": []})
        continue
    cpm = total[key]; p = int(round(percentile(round(cpm, 3), m["ref_cpm_sorted"])))
    hi, lo = p >= rules["high_percentile"], p <= rules["low_percentile"]
    if m["direction"] == "neutral" or not (hi or lo): level, tone = "typical", "neutral"
    else:
        level = "high" if hi else "low"
        bad = (m["direction"] == "bad" and hi) or (m["direction"] == "good" and lo)
        tone = "concern" if bad else "good"
    tot = sum(carriers[key].values()) or 1
    car = [{"species": s, "share": round(100 * v / tot)} for s, v in carriers[key].most_common(3) if v / tot >= 0.08][:2]
    mods.append({"key": key, "name": m["name"], "description": m["description"], "direction": m["direction"],
                 "measurable": True, "percentile": p, "cpm": round(cpm, 2), "level": level, "tone": tone, "carriers": car})

concern = [m for m in mods if m.get("tone") == "concern"]
if any(m["key"] in ("sulfur", "ging") for m in concern):
    note = ("Worth raising at your next vet visit: the bacteria behind gum disease are more active than in most dogs. "
            "Ask about a dental check. At home, daily tooth brushing helps most, and dental chews or a dental diet with the "
            "VOHC seal (vohc.org) are shown to reduce plaque and tartar.")
    kind = "vet"
elif concern:
    note = ("Fewer of the bacteria linked to healthier mouths than in most dogs. Daily tooth brushing and VOHC-approved "
            "dental chews (vohc.org) help keep the balance healthy.")
    kind = "watch"
else:
    note = "The mouth bacteria look typical. Keep up what you're doing; daily brushing and VOHC-approved chews help keep it that way."
    kind = "typical"
doc = dict(base, status="ok", enterobacterales_pct=round(entero, 1), modules=mods, note=note, note_kind=kind)
json.dump(doc, open(a.out, "w"), indent=1)
print("ok " + kind + " " + " ".join(f"{m['key']}={m.get('percentile', 'NA')}" for m in mods))
