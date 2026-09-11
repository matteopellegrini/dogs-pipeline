# Cosmo: Embark (2020, SNP array) vs ProsperK9 (cosmo3, 2.2x WGS)

Embark report ingested 2026-09-10 (37 pages; structured copy in
`embark_cosmo_report.json`, chromosome painting in `embark_cosmo_segments.json`).
ProsperK9 side = `results/cosmo3` on Hoffman (DOGS-Gen-50, Illumina 2.2x,
current 230-breed lasso, coat/OMIA from the Aug-2026 pipeline).

## Breed

| | Embark | ProsperK9 (lasso) |
|---|---|---|
| Poodle (all) | **75.1** (Std 40.2 + Small 34.9) | **61.2** (Std 36.9, Mini 15.2, Toy 9.1) |
| Labrador Retriever | 10.4 | 6.5 |
| Cocker Spaniel | 9.4 | 6.1 |
| English Cocker Spaniel | 5.1 | 3.7 |
| everything else | 0 | 22.5 (Barbet 3.6, Zagar 1.9, Pont-Audemer Spaniel 1.8, Cavalier 1.6, … ; Gray Wolf 1.1) |
| wolf | "wolfiness 2.8% HIGH" (ancient-allele metric) | 1.1% Gray Wolf ancestry component |

Same five breeds, same rank order (Poodle > Lab > Cocker > English Cocker), and
the Standard-vs-small Poodle split agrees (Embark 40/35, ours 37/24). The
difference is entirely the tail: Embark assigns 100% to five breeds; our lasso
leaves ~19% spread over water-dog/spaniel relatives of the same five (Barbet,
Pont-Audemer, Zagar are all Poodle/Cocker-adjacent). Embark's family tree
(4 grandparents: Std Poodle×Cocker, Small Poodle×English Cocker, Std Poodle mix,
Small Poodle×Lab) is a story consistent with both. **The FLARE painting is the
arbiter: if the tail breeds do not form segments, they are lasso smear and the
prune-and-renormalise presentation is right.**

## Lineages, inbreeding, size

| | Embark | ProsperK9 |
|---|---|---|
| Maternal (mtDNA) | B1 / B84 | haplogroup **B** (no subtyping) — concordant |
| Paternal (Y) | A1a / H1a.59 | not analysed (no Y pipeline) |
| Inbreeding | COI 10% | F_ROH 7.85% ("moderate"; 59 ROH segments, 173 Mb); Dog10K-relative F_ROH 0.169, 34th percentile |
| Predicted adult weight | 33 lb | **40.6 lb** (PRS blend) — check against Cosmo's actual weight |
| MHC (DRB1, DQA1/DQB1) diversity | High / High | not analysed |

## Coat and body traits

| Locus | Embark | ProsperK9 | |
|---|---|---|---|
| E (MC1R) | ee → cream/red, no mask | **e/e** (low conf, proxy SNP; e1 stop uncovered) | agree |
| I (intensity) | intermediate (yellow/tan) | not analysed | — |
| B (TYRP1) | Bb | **B/b** (medium) | agree |
| D (MLPH) | DD | **D/D** (high) | agree |
| K (CBD103) | KBKB | **KB/KB** (high) | agree |
| A (ASIP) | atat | undetermined (needs SV analysis) | gap |
| S (MITF) | SS little white | not typed | gap |
| M (PMEL) | mm | m/m (low; 2 reads at junction) | agree |
| H (PSMB7) | hh | (H locus added to pipeline 2026-09-08; cosmo3 predates) | — |
| Furnishings (RSPO2) | FF | not in catalogue | gap |
| Coat length (FGF5) | LhLh long | one FGF5 variant **alt/alt** (high, imputed) → long | agree |
| Texture (KRT71) | CC (wavy via furnishings) | c1/c2 indels not callable; SNP ref/ref | partial |
| Shedding (MC5R) | CT | no_call (insufficient reads) | gap |
| Muzzle (BMP3) | AC | **het** (high, imputed) | agree |
| Size genes (IGF1, IGFR1, STC2, GHR×2) | NN/GA/TT/GA/CC | expressed via PRS, not per-gene | — |
| POMC appetite | NN | indel_no_call | gap |

## Health

Embark: 197 conditions tested, **all Clear** (25 breed-relevant + 172 other).
ProsperK9 catalogue: 449 variants. Cosmo3 calls: 225 ref/ref, 7 het, 3 alt/alt,
**212 not called** (114 insufficient reads / not in panel, 98 indels, 2 not
callable) — i.e. 47% of the catalogue is uncallable at 2.2x with this pipeline
version.

Of the 46 catalogue variants in Embark-tested disease genes: **29 called ref/ref
(all concordant with Embark "Clear"), 17 not called** (5 indel: ADAMTS17×2,
HEXB×2, VWF; 12 insufficient reads: COL4A4×2, PFKM, ATF2, CNGA3, COLQ×2,
HCRTR2, MTM1×3, SUV39H2). No discordant calls.

Our non-reference disease-gene calls: BMP3 het (trait, matches Embark AC);
**TPO het (low confidence, SV-type read call)** — Embark reports both TPO
congenital-hypothyroidism variants Clear → likely a low-depth false positive on
our side, worth a read-level review; ATP7B/ATP1B2/ATP2A2/RELN het (not in
Embark's list; ATP1B2/ATP2A2/RELN low confidence).

## Take-aways

1. Breed identity and ranking agree exactly; only the tail differs. FLARE will
   settle whether the tail is real.
2. Every coat locus both tests call agrees (E, B, D, K, M, FGF5, BMP3). Our gaps
   are the untyped loci (A, S, furnishings, MC5R, intensity) — catalogue/method
   additions, not errors.
3. Health: perfect concordance where we call, but a **coverage gap** — a third
   of Embark's breed-relevant list is uncallable for a 2.2x dog (indels + low
   reads). The chip calls every one of its 197 sites every time. This is the
   real competitive weakness of low-pass WGS for the disease panel and argues for
   (a) indel calling from reads, (b) a "not assessable" count shown honestly.
4. Weight: 33 vs 40.6 lb — needs Cosmo's real weight to judge.
5. Wolf: Embark's "wolfiness" is not ancestry; our 1.1% Gray Wolf component is
   below the ~5% display threshold noted in the breed-panel README.
