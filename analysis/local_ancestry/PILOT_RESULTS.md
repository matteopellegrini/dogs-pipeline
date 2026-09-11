# Local ancestry pilot — results (2026-09-11)

Two-step design (user): global lasso proportions -> FLARE prior; FLARE's genome-wide
mean of the local posteriors = final breed proportions; painting from AN1/AN2.
Reference: merged Parker + Dog10K (+42 Darwin's Ark APBT) phased at the 131k panel
sites, all 231 production breeds (min-n 6), 2,937 dogs. Query: GLIMPSE2 haplotypes
at the same sites (113,316 matched). Runs on the Mac (`run_dog_lai.py`).

## Runtime: the reference must be lasso-restricted
| reference | chr1 | genome |
|---|---|---|
| all 231 breeds / 5,874 haplotypes, 1 thread (cluster) | — | >75 min, never finished |
| 27 breeds (lasso > 0) / 624 haplotypes, 8 threads | 9 s | ~3 min (est.) |
| 7 breeds (lasso >= 2%) / 162 haplotypes, 8 threads | <1 s | **1–2 s** |
Memory is ~3 GB regardless (the 15 GB on Hoffman was JVM reservation). EM on/off makes
no difference at small reference size; with the prior, `em=false` is used.

## Genome-wide proportions (cosmo3, 2.2x, 2% floor, gen=3) vs lasso vs Embark
| breed | lasso | FLARE | Embark |
|---|---|---|---|
| Standard Poodle | 45.5 | 46.3 | 40.2 |
| Small poodles (Mini+Toy) | 30.0 | 34.4 | 34.9 |
| Labrador Retriever | 8.0 | 8.8 | 10.4 |
| Cocker Spaniel | 7.5 | 5.2 | 9.4 |
| English Cocker Spaniel | 4.6 | 3.4 | 5.1 |
| Barbet | 4.4 | 1.9 | — |
Stable across depth after the sink guard (Standard Poodle 0.45–0.49 at 2.2/1.0/0.5/0.25/0.1x)
and across platforms (MGI replicates ucla-20001/6164: 0.435/0.410; they shift Toy->Mini
and Cocker up to 13–15%).

## Sink breeds
A small lasso component (2–4%) from a tiny reference breed (ZAGAR n=7) absorbs
related-breed segments in the local step (9.5% -> 17.5% of the genome at 0.1x).
Guard in `run_dog_lai.py`: lasso < 5% and painted share > 3x lasso -> prune, re-run.
Barbet (n=6, Poodle relative) is a *local* sink (paints whole chr27/chr32) without
ballooning globally — not caught by the guard; removing it did not improve Embark
agreement (62.2 -> 59.7% family level), so the disagreement is broader.

## Painting vs Embark (compare_flare_embark.py; positions as chromosome fractions,
unordered haplotype pairs; Embark's painting is itself an estimate, not truth)
- exact 5-class pair: ~40%; at least one haplotype: ~90%
- family level (Poodle / Cocker / Lab / other): 60–66% (cosmo3 and MGI replicates)
- segments: FLARE 134–280 (gen=3) vs Embark 526 — Embark paints finer/more switches
- gen=3 (Embark's tree: purebred great-grandparents) halves segment count vs gen=10
  with the same proportions.

## Verdict
Two-step pipeline is fast and its proportions are consistent (depth, platform,
Embark totals). The PAINTING is not validated: no ground truth, and Embark disagrees
at the segment level. Next: (1) simulated admixed dogs spliced from held-out
reference haplotypes (known truth) -> segment accuracy vs depth via the cosmo3
downsample recipe; (2) an F1 cross if one exists in the cohort (one Poodle + one
Cocker haplotype on every chromosome is an exact test); (3) posterior-based
"uncertain" class + minimum segment length before anything is shown.
