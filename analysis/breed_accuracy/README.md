# Breed panel: validation and the Parker + Dog10K merge

Everything here evaluates the breed reference behind Stage 9 of
`run_dog_pipeline.sh`. State as of 2026-08-09.

## The one thing to know first

**Phat is the per-breed empirical allele frequency.** Verified against the
genotypes at corr 0.9997+, and against Matteo's own PLINK output
(`scope_freq.frq.strat`) at 0.9993-0.9998. So a reference panel for any set of
individuals is just their per-breed allele frequencies at the shared SNPs —
which means panels can be built and compared in minutes without re-running
SCOPE.

**But SCOPE's Qhat is not the same as the pipeline's NNLS projection.** SCOPE
estimates a latent subspace by randomized PCA first, then solves Q against that
rank-k reconstruction. That smoothing regularizes it. Measured on the Parker
panel, in-sample:

    SCOPE Qhat top-1 == own breed    96.83%
    NNLS       top-1 == own breed   100.00%   <- circular
    the two agree with each other    96.83%

So `Qhat` on the full reference is an honest-ish estimate, while NNLS in-sample
is not. Practically: **SCOPE Qhat in-sample ~ NNLS leave-one-out** (96.8% vs
97.2%). Stage 9 uses NNLS for customer dogs, so the two estimators differ on
~3% of dogs — that is a real property of production, not an artifact.

## Evaluation methods, weakest to strongest

| method | Parker panel | note |
|---|---|---|
| NNLS in-sample | 100.0% | worthless — every dog defines its own target |
| SCOPE Qhat in-sample | 96.8% | LSE regularizes; roughly honest |
| NNLS leave-one-out | 97.2% | contribution removed analytically |
| **held out on Dog10K** | **96.2%** | nothing shared between fit and test |

**Always prefer held-out.** Leave-one-out hid a real regression: pooling the
regional Salukis looked like an improvement under LOO but cost three Anatolian
Shepherd calls on held-out data, because a pooled `SALU` (n=19) becomes a broad
attractor. Reverted in `9f3be59`; they are summed for display instead.

**General rule that came out of this: merge for display, not in the model,**
unless held-out data says otherwise. Display grouping cannot change what the
model competes over.

## Scripts

| script | what it does |
|---|---|
| `check_reference.py` | leave-one-out on the Parker panel; flags for `--no-wolves`, `--merge-wolves`, `--merge-saluki`, `--akc-display` |
| `project_dog10k.py` | projects every Dog10K sample onto Parker — the held-out test |
| `merged_panel.py` | builds Parker / Dog10K / merged panels and evaluates by NNLS; `--cross` for fully held-out |
| `make_merged_plink.py` | writes a merged `.bed/.bim/.fam` + matching `.frq.strat` for supervised SCOPE |
| `compare_breeds.py` | predictions vs owner-reported breeds |
| `build_scope.sh` | clones, patches and builds SCOPE (three macOS fixes documented inside) |

`compare_breeds.py` takes the owner data as an argument and never reads it from
a fixed path: `dog_details.xlsx` contains owner names and email addresses and
**must not enter this repo, which is public.**

## Reproducing the merged-panel run

```bash
# 1. Parker sites out of the Dog10K panel (~45 s)
awk 'BEGIN{OFS="\t"} {c=$1; if(c !~ /^chr/) c="chr"c; print c,$4-1,$4}' \
    COSMO/analysis/cosmo_parker_full.bim | sort -k1,1 -k2,2n > sites.bed
bcftools query -R sites.bed -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
    dog10k_panel/AutoAndXPAR.Dog10K.phased_plus_disease_rh.bcf > dog10k_gt.txt
bcftools query -l dog10k_panel/...bcf > dog10k_samples.txt

# 2. merged PLINK + frequency file (~95 s for all shared SNPs)
python3 analysis/breed_accuracy/make_merged_plink.py \
    COSMO/analysis/cosmo_parker_full COSMO/analysis/scope_clust.txt \
    dog10k_gt.txt dog10k_samples.txt merged --step 1 --min-n 2

# 3. supervised SCOPE (~5 min at 131k SNPs, k=381)
scope -g merged -freq merged.frq.strat -k 381 -o out_ -nt 4 -seed 1
```

SNP counts, so nobody re-derives them: Parker has 143,933; 143,385 (99.6%) are
present in the Dog10K panel; 136,502 survive allele matching; **131,353 are
typed in both**.

## Results

Supervised SCOPE Qhat, merged panel, 3,238 dogs, 381 breeds:

| | 32,839 SNPs | 131,353 SNPs |
|---|---|---|
| overall | 90.70% | 90.27% |
| n=2-3 | 40.1% | 38.9% |
| n=4-7 | 84.8% | 83.9% |
| **n>=8** | **96.4%** | **96.2%** |
| Parker dogs | 96.1% | 96.1% |
| Dog10K dogs | 86.9% | 86.1% |

**4x more SNPs changed nothing.** The panel is SNP-saturated at ~33k for 381-way
discrimination — accuracy is limited by reference composition, not marker
density. Future effort belongs in sample counts and label harmonisation, not
more markers.

**The headline understates the merge.** Parker dogs score 96.1% in the merged
panel against 96.83% in Parker-only SCOPE: doubling the breed count costs the
already-well-covered dogs 0.7 points, while adding 250 breeds, 280 village dogs
and 57 wolves.

Platform worry was unfounded and is now measured: for the 119 breeds in both
sources, Parker-derived vs Dog10K-derived frequencies differ by 0.0347 against
0.1162 between different breeds. Array and WGS agree 3.3x more closely than
breeds differ, so pooling injects no meaningful batch structure.

## What Dog10K adds (1,929 samples, 369 prefixes)

| | Parker | Dog10K |
|---|---|---|
| breeds shared | 119 | 119 |
| breeds Parker lacks | — | 250 |
| wolves | 7 populations, n=1 each | 11 (`CLUP*`), **57 samples** |
| village dogs | **none** | 25 countries (`VILL*`), **280 samples** |

Village dogs are the structural gap. Parker is built from registered pedigree
breeds and cannot represent unbred dogs at all, so that ancestry is forced onto
whichever pedigree breeds fit least badly — which likely explains the Italian
landraces (`CPAT_Italy`, `MAAB_Italy`) predicting as German Shepherd, and why
only 50% of mixed-breed owner reports had all named components in the top 3
against 95% top-1 for purebreds.

## To finalise: three cleanups, none touching the model

Of the 315 errors in the full run:

1. **Drop populations with n<8.** 68 breeds at n=2-3 score 38.9% and can only
   absorb ancestry they cannot support. This single filter takes the panel from
   90.3% to 96.2%.
2. **Harmonise duplicate codes across panels.** `TERV`/`TURV` (Belgian
   Tervuren), `CIRN`/`CIRN_Italy`, `POMR`/`GSPK`/`POM`, `FXTE`/`WFOX`/`SMFX`.
   These count as errors while being pure bookkeeping. Needs the Dog10K
   breed-code table — the same table that would settle `SSHP`.
3. **Pool village dogs and wolves regionally.** 12% of errors are
   village->village and 3% wolf->wolf (`VILLTH`->`VILLCN`, `CLUPCN`->`CLUPRU`).
   These are geographic continua, not breeds. Pooling converts the errors into
   correct calls and gives a more robust "village dog ancestry" signal for the
   mixed customers where the model is weakest.

## Open questions

- **`SSHP`** — Smooth Collie or Shetland Sheepdog? The code is consistent across
  both panels, but Dog10K's `SSHP` dogs show heavy `COLL` secondary (0.21-0.24
  in half of them), leaning Smooth Collie; against that, neither panel has a
  `SHED` code and a 1,929-dog panel should contain Shelties. **The
  `COLL`+`SSHP` display merge stays backed out until this is settled.**
- **`GDJK`** ("Grand Basset Griffon Vendéen", n=2) — the code matches no obvious
  breed name and Petit Basset uses the sensible `PBGV`. Label unverified.
- `NELK` is **settled**: Norwegian Elkhound, not Norrbottenspets. Confirmed by
  Matteo's own sample and by all 7 Dog10K `NELK` dogs projecting onto Parker
  `NELK` at 0.51-0.62. Fixed in `a22041d` — it had been mislabelling live
  reports.
