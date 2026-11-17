# Imputation panel v3: Dog10K + DAP + streamed public genomes + all catalogue variants

Plan, 2026-11-15. Goal: one imputation panel that carries every catalogue variant we have seen in a deep genome, so stage 12
(GLIMPSE2) imputes each of them for every kit from its haplotype, as it already does for the 105 catalogue SNVs that are
Dog10K sites. Replaces the per-variant windows (MDR1, GH1, PDK4, STRN, FGF5, cord1) with one build.

## Inputs

| Source | Dogs | Panel-site genotypes | Catalogue-site genotypes |
|---|---|---|---|
| Panel v2 (Dog10K) | 1,929 | phased (ref_panel_v2) | SNVs: the 105 panel sites; indels: Dog10K non-SNP calls where trustworthy (44 catalogue indels with carriers), else missing |
| Panel v2 (DAP) | 7,828 | phased (imputed from ~0.8x) | missing (never sequenced deep enough) |
| Streamed public genomes | ~9,500 at <=10x (613 pilot/existing, 1,378 breed-gap, 7,548 all) | sra_panel/<BS>.panel_sites.bcf (unphased) | read from sra_slices/<BS>.bam: mpileup/call at the 548 small sites + seq_allele_caller.py for sites a pileup misses |

Not in this build: chrX non-PAR (62 catalogue sites; panel v2 has no non-PAR X), the 121 structural variants (phase 2: sv_evidence.py on the slices).

## Steps

1. **Catalogue genotypes for every streamed dog** (cluster job over the slices, like catalogue_public_carriers.sh).
   - GT only at DP >= 6 (else missing). Hard sites (long/complex insertions, homopolymer runs, cord1, PNPLA1, ABCA4) from seq_allele_caller.py, with a signature added per site.
   - The 426 slices cut before 2026-11-15 lack 24 small sites: those dogs stay missing there.
   - QC per site: drop HOLD/artefact sites (BTBD17, HSF4 homopolymers, ABCB1_442, IDUA, CYTB, MT-TY) and any site whose streamed-dog frequency is impossible for its disease (vs Donner 2023 where matched).
2. **Relatedness and duplicate check**: KING on a thinned SNV set across panel v2 + streamed dogs; drop duplicates of panel dogs, keep one of each identical pair, drop contaminated genomes (foreign-allele rate, as for cosmo3).
3. **Phase the streamed dogs onto panel v2** per chromosome (Beagle 5.5, ref = panel v2, impute=false; pilot recipe 4c9a236/3618fc9), then merge into the panel ("pall" in the pilot). 38 autosomes + X PAR. Drop AC==0 / AC==AN sites and exact duplicates (GLIMPSE2 aborts on them).
4. **Add the catalogue sites, without forcing them into dogs where they were never observed.** 2026-11-15 cord1 test: letting Beagle fill the indel into 1,929 Dog10K dogs (missing in ~90% of samples, 16 observed carriers) gave 12% allele frequency and put cord1 in 25% of wolves and 83% of Polish Lowland Sheepdogs. So the catalogue sites live in a second reference made only of dogs whose catalogue genotypes were read directly (the ~9,500 streamed genomes, phased with their catalogue sites), and kits are imputed against it by window from their own stage-12 genotypes (indel_window_target_impute.sh). SHAPEIT5_phase_rare stays a candidate for phasing the observed genotypes onto the step-3 scaffold, but only for dogs with an observed genotype.
5. **Which catalogue sites go in**: those with >= 3 carrier dogs genotyped directly (not filled in). Today ~34 qualify (43 minus ~9 known artefacts); the 37 sites with 1-2 carriers join as streaming finds more. The rest stay on reads + breed prior.
6. **Validate** (before any customer change)
   - Genome-wide accuracy: the pilot scoring (Cosmo 0.1-2.2x + held-out SRA dogs at 0.1-1x) on all chromosomes, pv2 vs v3.
   - Catalogue sites: hold out ~10% of streamed carriers per site from the build, impute their 1x down-sample, score; and compare with customers' own reads at the site (as for MDR1 and cord1).
   - Calibration: cohort imputed frequency per site vs Donner 2023 mixed-breed frequency.
7. **Wire into production**
   - Stage 12 chunks rebuilt for the new panel; panel_v3_breed_prior.py: catalogue sites in the panel get basis 'haplotype' (imputed GP floored at 0.1% x reads at the site), replacing the per-window indel_imp path; refit the frozen model.
   - Recalibrate what depends on the panel: inbreeding depth calibration, MHC vs_inbreeding, breed-panel scoring check (imp2 lesson).
8. **Cohort re-impute in staging** (results_imp3, pi_shared), before/after comparison (condition calls, carriers, breeds, F_ROH), Matteo's go, promote + silent republish.

## Timing and resources

- Streaming: breed-gap list ~3-4 days at 80 concurrent; the 7,548 list ~1-2 weeks depending on FAST mode (testing now).
- Build: Beagle phasing hours per chromosome at ~19k dogs (pilot: 13 min for chr38 with 613 added dogs, 42 GB); phase_rare per 5-Mb chunk, modest memory.
- Cohort re-impute: ~1 day on pi_shared (imp2 precedent).
- Disk: panel_sites.bcf ~190 MB/dog -> ~1.8 TB for 9,500 dogs; panel v2 is 48 GB, v3 ~100 GB.

## Proposed sequence

- Dry run of steps 1-6 once the breed-gap list finishes (~2,000 added genomes): proves the build and gives a first catalogue-site validation.
- Customer rollout once, after the full 7,548 list: one re-impute and one republish instead of two.

## Catalogue fixes found against Donner et al. 2023 (2026-11-15)
- IMPG2 PRA4 (new_IMPG2_1253): was the gene span; now the LINE-1 insertion site canFam4 chr33:8,197,950-8,197,966
  (17-bp TSD, Hitti-Malin et al. 2020). In the SV pass of this build since 2026-11-16 (`sv_evidence.py sites`): short
  reads can exclude it from spanning reads but cannot confirm a copy (sv_confirmable), so dogs read clear or uncertain
  only. `sites` now reads an insertion's pos..end span as its target-site duplication (breakpoint = its middle, size
  unknown unless SVLEN): this also stops FAM161A_1706's 15-bp TSD passing as a 15-bp insertion that reads could confirm.
- SLC7A9 cat338 (cystinuria type II-B, Miniature Pinscher) is at the right position; the paper finds 0 carriers in
  1,054,232 dogs. The 2.2%-frequency SLC7A9 variant in our reads is p.A217T (Donner "type I-B", recessive, incomplete
  penetrance, chr1:120,053,689), which the catalogue lacks; adding it is a separate decision.
- Paper describes a different change from our catalogue entry (published_priors.py rejected 7). Re-curated 2026-11-16
  (build_panel_v3_sites.py FIXES, checked against OMIA HGVS and the canFam4 sequence):
  - GALC cat006 -> 78-bp insertion (Irish Setter, SVLEN=78) and DLX6 cat167 -> 2,056-bp LINE-1 insertion (NSDTR,
    SVLEN=2056): now SV-pass sites (clear/uncertain only); the catalogue typed an invented SNV at the insertion point.
  - SLC19A3 cat329 -> C>TTGCA (c.624delinsTGCAA); COL1A2 cat232 -> CTGA>TGTCATTGG (c.3856_3859delins; the paper's
    "CTAG" is CTGA in canFam4). Both were the first base of the delins as an SNV.
  - Ours was right: NIPAL4 cat061 (the paper used OMIA's old g.52737379del, OMIA corrected it to g.52737279del = our
    site), L2HGDH cat281 (our G>A is the first base of the GAA>AAG double substitution; the second base is A>G at
    26,947,411), ITGA2B cat401 (= the mixed-breed R453*; the paper row is the Great Pyrenees 14-bp dup, OMIA 369, which
    the catalogue lacks; breed prior now from OMIA 1568, not the Otterhound allele cat400).
  - published_priors.py: SAME accepts cat061/cat281; with the re-curated alleles 213 sites get published priors (+cat061
    American Bulldog 8.8%, cat281 Staffordshire Bull Terrier 1.4%, cat232, cat329).
  - Live reports still use the frozen model (old alleles at cat006/cat167/cat232/cat329, all 0-1 per 10k dogs in the
    paper): they change at the panel v3 rollout, or earlier by per-kit rescore on Matteo's go.
