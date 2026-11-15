"""Score the F_ROH depth-calibration BCFs with bcftools roh (--GTs-only 30, --AF-tag RAF, autosomes).
Rows: dog, level (0.1..2 or full), F at >=500 kb and >=1 Mb, segment count, Mb in 0.5-1 / 1-5 / >5 Mb classes.
  envs/genomics/bin/python3 analysis/imputation_panel/froh_calib_score.py > scratch/froh_calib/scores.tsv"""
import os, subprocess, sys
from multiprocessing import Pool
D = "/u/project/pellegrini/gkislik/dogs"
B = f"{D}/envs/genomics/bin/bcftools"
KEEP = "/u/project/pellegrini_archive/data/dogs_imputed_v2"
AUTO = {f"chr{i}" for i in range(1, 39)}
G = 2.2e9

def roh(bcf):
    r = subprocess.run(f"{B} roh --GTs-only 30 --AF-tag RAF -O r {bcf} 2>/dev/null", shell=True,
                       capture_output=True, text=True).stdout
    L = [int(x[5]) for x in (l.split("\t") for l in r.splitlines()) if x[0] == "RG" and x[2] in AUTO]
    big = [l for l in L if l >= 500e3]
    return (sum(big) / G, sum(l for l in L if l >= 1e6) / G, len(big),
            sum(l for l in big if l < 1e6) / 1e6, sum(l for l in big if 1e6 <= l < 5e6) / 1e6,
            sum(l for l in big if l >= 5e6) / 1e6)

def job(t):
    dog, lev, bcf = t
    return (dog, lev, *roh(bcf))

todo = []
for row in open(f"{D}/scratch/froh_calib_dogs.tsv"):
    dog, bam, depth = row.rstrip("\n").split("\t")
    for lev in ("0.1", "0.25", "0.5", "1", "2"):
        f = f"{D}/scratch/froh_calib/{dog}/{lev}.bcf"
        if os.path.exists(f): todo.append((dog, lev, f))
    f = f"{KEEP}/{dog}.bcf"
    if os.path.exists(f): todo.append((dog, f"full:{depth}", f))
print("dog\tlevel\tf500\tf1000\tn500\tmb_0.5_1\tmb_1_5\tmb_5plus", flush=True)
with Pool(int(os.environ.get("NPROC", "8"))) as p:
    for r in p.imap(job, todo):
        print(r[0], r[1], *[f"{v:.4f}" if isinstance(v, float) else v for v in r[2:]], sep="\t", flush=True)
