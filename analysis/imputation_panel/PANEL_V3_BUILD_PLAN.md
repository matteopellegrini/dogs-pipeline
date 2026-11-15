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
4. **Add the catalogue sites onto that scaffold**: SHAPEIT5_phase_rare (scaffold = step 3 panel; input = all dogs at panel sites + catalogue sites, catalogue GT missing for DAP and most Dog10K dogs). It phases each observed allele onto its haplotype and fills the missing dogs from the haplotypes they share. Test first on chr15 against the cord1 window build (job 1748821).
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
