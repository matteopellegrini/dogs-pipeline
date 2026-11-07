# snpEff on canFam4 (2026-11-07)

`ROS_Cfam_1.115` (the stock snpEff dog database, Ensembl 115) is **ROS_Cfam_1.0**, not canFam4 (UU_Cfam_GSD_1.0).
The pipeline annotated canFam4 positions with it, so gene names and effects in `functional_variants.json` were
computed against gene models at the wrong place.

## Database: `canFam4_NCBI_AR106`
`cluster/build_snpeff_canfam4.sh` builds it from NCBI RefSeq annotation release 106 of GCF_011100685.1 (= canFam4).
The RefSeq accessions are renamed to UCSC names from the assembly report, and the sequence comes from `canFam4.fa`.
It is installed in `envs/genomics/share/snpeff-5.4.0c-0/data/` on Hoffman2 (485 MB).
- Size: 37,987 genes (22,083 protein-coding), 98,085 transcripts.
- Protein check: 62,257 of 63,650 NCBI proteins match snpEff's own translation of canFam4.fa (97.8%). There are 1,393 errors (2.2%), the expected rate for RefSeq models with genome exceptions or partial models.
- chrM has no RefSeq annotation in this assembly, so it gets no genes. Stage 10 has no chrM sites.
- Known gaps in the annotation:
  - NCBI merges COLQ and HACL1 into a single gene, LOC608697.
  - SLC19A3 and CHST6 are named LOC486151 and LOC489707. Their products are thiamine transporter 2 and carbohydrate sulfotransferase 6.
  - About 3,780 protein-coding genes are `LOC` IDs. Only 116 of them have a human ortholog in NCBI `gene_orthologs`, so none were renamed.

## Validation: `panel_gene_test.sh` (Panel v3 small sites, kind snv/indel, status live*/verified*, 523 sites)
| | ROS_Cfam_1.115 | canFam4_NCBI_AR106 |
|---|---|---|
| Sites with Dog10K records within 1 kb (480): panel gene named in a genic ANN | 33 | **442** |
| Panel variants themselves (523): own gene named | 31 | **475** |
| WARNING_REF_DOES_NOT_MATCH_GENOME (nearby records) | 610 | **0** |

The 48 panel variants that still miss break down as follows:
- Synonyms (17): FVIII/F8, FIX/F9, FVII/F7, G6PC/G6PC1, ACPT/ACP4, FAM134B/RETREG1, LARGE/LARGE1, TYPR1/TYRP1, MOP/MPO, T-BOX/TBXT.
- `LOC` names (5): COLQ, SLC19A3, CHST6, LOC477558.
- Regulatory or intergenic sites by design: SHH ZRS, RSPO2, MSTN, IGF1-AS, NEBL1-6, chrM.
- The panel-side questions already noted in `gene_table_ros_fix.md`: FAM20C, MICOS13, CYP27B1.

Per-site detail is in `panel_gene_test/*.sites.tsv`.

## Dog10K panel re-annotation
`cluster/dog10k_annotate.sh` (39 tasks plus a finalize job) gives 30,452,768 sites, the same count as before.
- HIGH-impact sites: 25,372 to **6,430**.
- Rare (<1%) loss-of-function sites in disease genes: 5,441 to **692**. Of those, 99 are seen homozygous in Dog10K (was 1,458).
- `cluster/dog10k_indel_annotate.sh`: HIGH-impact indel records 291,993 to **53,387**. This probably explains why the Dog10K "damaging" indels looked like artefacts during the indel-imputation work.

The new panel annotation is staged on Hoffman2 as `reference_panel/dog10k_ann/dog10k_sites.ann.canFam4_NCBI_AR106.bcf`, with per-chromosome files in `work_dog10k_ann_cf4/` and `work_dog10k_indel_cf4/`. The live `dog10k_sites.ann.bcf` (ROS) is **not swapped yet**.

## App assets (dogs-app `lib/`, not committed there)
`build_panel_hom_lof.py` and `build_cohort_hom_counts.py` reproduce the shipped files from the old annotation:
- panelHomLof: 5,452 sites against 5,441 shipped, with all 1,458 shipped homozygous keys recovered.
- cohortHomCounts: identical, 567 of 567 genes, over the 96 DOGS-Gen dogs.

The canFam4 versions are in `app_assets_canFam4/`:
- `panelHomLof.json`: 692 sites, 99 seen homozygous.
- `cohortHomCounts.json`: 128 genes, none hit in 5 or more dogs.

## Cohort rebuild, sized by running it into scratch (`snpeff_build/fv_new/`, nothing published)
`cluster/rerun-10-fv-list.sh` with `CHUNK=10` processed the 1,978 published kits (`logs/coat_b1101_list.tsv`).
- Each array task copies ANN onto the kit's imputed BCF and runs `functional_variants.py`.
- About 4 minutes and under 8 GB per kit. With 100 tasks running it took about 1.5 hours wall time. All 1,978 kits finished, with 0 failures.
- Parity check: run with the old ROS panel, the refactored parser reproduces a published `functional_variants.json` byte for byte.
- HIGH rows: 4,191,130 to **613,497**. Homozygous HIGH rows: 1,435,334 to 202,302. Every report's list changes.
- Rows listed by default (OmiaTable PredictedVariants, `listed_rows.py`): **75 rows on 66 reports, 14 inside the named gene, become 7 rows on 7 reports, 7 inside**. The 7 are FICD, MASP2 x2, NFATC1 x2, ANAPC1 and CA12 stop-gain or splice-site changes. See `listed_rows/`.

## Rollout (approved and run 2026-11-07)
- Merged to `dogs` main (53d0a97). Scripts synced to Hoffman2.
- Panel swapped: `dog10k_sites.ann.bcf` is now canFam4_NCBI_AR106; the old one is kept as `dog10k_sites.ann.ROS_Cfam_1.115.bcf` (+ .csi).
- All 1,978 published `functional_variants.json` replaced from `snpeff_build/fv_new/` (identical panel + parser), each with `functional_variants.pre-snp1107.json` beside it.
- dogs-app 6dfba8e: new `lib/panelHomLof.json` + `lib/cohortHomCounts.json`, deployed to production.
- Files-only republish tag `snp1107` (log: `logs/republish_snp1107.log` on Hoffman2).
- Cached AI text: the republish itself deletes every kit's `ai_summary_<tab>.json` and bumps `report_uploaded_at` (so dashboards go stale and regenerate on view). For the 60 published reports that listed rows by default (`ai_cache_cleared_barcodes.txt`), both were also cleared outright via `DELETE /api/admin/summaries {barcodes, dashboard:true}` (dogs-app 8928808): 300 summary keys, 1 dashboard row (the only one cached).

Rollback: move the ROS BCF back, copy each `functional_variants.pre-snp1107.json` over `functional_variants.json`, revert dogs-app 6dfba8e, republish.

## Rollout plan (each step needs approval)
1. Merge this branch. Copy `run_dog_pipeline.sh`, `analysis/snpeff/functional_variants.py` and the cluster scripts to Hoffman2. Then move the canFam4 panel BCF into place as `dog10k_sites.ann.bcf`, keeping the old one as `dog10k_sites.ann.ROS_Cfam_1.115.bcf`. New kits are correct from then on.
2. Rebuild the published kits in place:
   `CHUNK=10 sbatch -A pi_matteop -J rerun10fv -a 1-198%100 -p campus -c 1 --mem=8G -t 4:00:00 cluster/rerun-10-fv-list.sh logs/coat_b1101_list.tsv snp1107`
   This keeps `functional_variants.pre-snp1107.json`. Alternatively, copy the finished scratch files from `snpeff_build/fv_new/<pk>/` after the swap, which takes minutes.
3. In dogs-app, copy `app_assets_canFam4/*.json` over `lib/panelHomLof.json` and `lib/cohortHomCounts.json`, and change the AR106 source wording if it is shown anywhere. Deploy.
4. Silent files-only republish: `FILES="functional_variants.json" bash cluster/republish_files.sh snp1107 logs/coat_b1101_list.tsv`. Do steps 3 and 4 together, so the app assets and the files match.
5. Clear cached AI dashboard and summary text that cites predicted variants. Today the Advisor gets the top 30 HIGH rows, so this affects up to every report. At minimum clear the 60 published reports that listed rows by default. The other 44 `functional_variants.json` files in `results/` are unpublished test or downsample directories.
6. Optional follow-ups: a symbol alias map for the panel synonyms; a curated name for the COLQ, SLC19A3 and CHST6 `LOC` genes; then re-check the indel imputation that rejected Dog10K damaging indels.
