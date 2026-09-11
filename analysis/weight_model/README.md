# Adult weight predictor: breed composition x Darwin's Ark size PRS

Fitted and validated by `fit_weight_model.py` on 2026-09-10; model shipped as `reference_json/weight_breed_model.json`, consumed by stage 11 of `run_dog_pipeline.sh`.

## Data

- 935 ProsperK9 customer-reported weights (`$D/analysis_weights.tsv` on Hoffman, from `merged_microbiome_age_weight_3.18_final.csv`); 879 have `breed_result.json` + `prs_result.json` in `$D/results_prosper/pk-<kit>/`.
- **302 of 879 are under 1 year old** and were dropped: a puppy weighed at swab time is not an adult-weight label. This alone accounts for most of the "miscalibration" originally reported (all ages: r 0.63, MAE 8.2 kg, bias +4.1 kg; adults: r 0.73, MAE 6.5 kg, bias +0.4 kg).
- The remaining heavy-dog shortfall is label noise: the ten heaviest "adults" are 65-77 kg pit-bull mixes whose ancestry predicts ~28 kg.
- 471 of the 577 adults have a top breed below 50%; only 20 are >=80% one breed. Purebred behaviour is therefore barely tested by CV and was checked separately (table below).

## Cross-validation (5 repeats x 5 folds on 577 adults, never in-sample; shrinkage by nested CV)

Bias columns are mean(pred - actual) in kg and % of bin mean, binned by **actual** weight. Regression to the mean is expected in these columns for any r~0.73 predictor; the "by predicted" table below is the calibration check.

| model | r | MAE kg | bias kg | <8 | 8-15 | 15-25 | 25-40 | >40 | dMAE vs current [95% CI] |
|---|---|---|---|---|---|---|---|---|---|
| current pipeline (dense LMM x DA blend) | 0.729 | 6.47 | +0.44 | +0.4 (+7%) | +2.0 (+18%) | +5.1 (+24%) | +0.2 (+1%) | -17.4 (-34%) |  |
| AKC breed-standard prior only (no fit) | 0.649 | 7.49 | +0.28 | +9.5 (+170%) | +6.9 (+61%) | +2.5 (+12%) | -4.8 (-15%) | -22.1 (-44%) | +1.03 [+0.61, +1.47] |
| linear recal  y ~ a + b*pred | 0.727 | 6.09 | +0.00 | +2.7 (+49%) | +3.2 (+28%) | +4.1 (+19%) | -1.7 (-5%) | -19.7 (-39%) | -0.38 [-0.53, -0.22] |
| quadratic recal | 0.730 | 6.15 | +0.01 | +2.0 (+36%) | +3.3 (+29%) | +4.6 (+22%) | -1.7 (-5%) | -20.1 (-39%) | -0.32 [-0.49, -0.15] |
| linear recal, Huber | 0.727 | 5.94 | -1.22 | +2.2 (+40%) | +2.3 (+21%) | +2.7 (+13%) | -3.3 (-11%) | -21.3 (-42%) | -0.52 [-0.76, -0.29] |
| log-linear recal | 0.722 | 6.20 | -2.09 | +3.0 (+54%) | +3.2 (+28%) | +1.8 (+8%) | -5.5 (-18%) | -24.2 (-48%) | -0.27 [-0.63, +0.08] |
| NNLS breed model, prior = global mean | 0.709 | 6.40 | +0.24 | +4.1 (+73%) | +4.6 (+41%) | +4.7 (+22%) | -2.7 (-9%) | -21.3 (-42%) | -0.09 [-0.44, +0.26] |
| NNLS breed model, AKC prior | 0.718 | 6.27 | +0.28 | +3.4 (+61%) | +3.8 (+34%) | +4.6 (+22%) | -2.1 (-7%) | -19.7 (-39%) | -0.22 [-0.52, +0.08] |
| NNLS breed, AKC prior, Huber | 0.722 | 6.06 | -0.85 | +3.2 (+58%) | +3.2 (+29%) | +3.1 (+14%) | -3.6 (-12%) | -21.0 (-41%) | -0.42 [-0.74, -0.09] |
| NNLS breed, AKC prior, log space | 0.726 | 5.76 | -1.43 | +1.6 (+29%) | +1.7 (+15%) | +2.8 (+13%) | -3.6 (-12%) | -20.9 (-41%) | -0.72 [-1.05, -0.38] |
| M2_akc + LMM prs_z | 0.723 | 6.15 | -0.02 | +2.6 (+46%) | +3.2 (+28%) | +4.2 (+20%) | -1.9 (-6%) | -19.6 (-39%) | -0.33 [-0.56, -0.12] |
| M2_akc + DA size PRS | 0.736 | 5.99 | -0.02 | +1.9 (+33%) | +3.0 (+27%) | +4.6 (+22%) | -1.6 (-5%) | -19.7 (-39%) | -0.48 [-0.72, -0.26] |
| M2_akc + dense pred | 0.723 | 6.15 | -0.02 | +2.6 (+46%) | +3.1 (+28%) | +4.2 (+20%) | -1.9 (-6%) | -19.6 (-39%) | -0.34 [-0.56, -0.12] |
| M2_akc + height pred | 0.727 | 6.23 | -0.03 | +2.3 (+41%) | +3.0 (+27%) | +4.5 (+21%) | -1.8 (-6%) | -19.9 (-39%) | -0.24 [-0.47, -0.01] |
| M2_akc + all four | 0.731 | 6.03 | -0.03 | +1.8 (+32%) | +3.0 (+27%) | +4.5 (+21%) | -1.5 (-5%) | -19.7 (-39%) | -0.45 [-0.65, -0.25] |
| M2_akc + current pred | 0.732 | 6.02 | -0.01 | +2.3 (+41%) | +2.9 (+26%) | +4.3 (+20%) | -1.6 (-5%) | -19.5 (-38%) | -0.46 [-0.64, -0.28] |
| M2_akc + current pred, Huber | 0.733 | 5.87 | -1.26 | +1.8 (+32%) | +2.1 (+19%) | +2.9 (+14%) | -3.3 (-11%) | -21.2 (-42%) | -0.61 [-0.86, -0.37] |
| bounded log breed model (box 1.5x prior) | 0.727 | 5.88 | -1.17 | +2.8 (+50%) | +2.4 (+22%) | +2.8 (+13%) | -3.7 (-12%) | -21.2 (-42%) | -0.61 [-0.94, -0.29] |
| M2B + current pred | 0.735 | 5.93 | -0.01 | +2.3 (+41%) | +2.8 (+25%) | +4.3 (+20%) | -1.6 (-5%) | -19.2 (-38%) | -0.56 [-0.76, -0.36] |
| M2B + current pred, Huber | 0.735 | 5.79 | -1.24 | +1.8 (+32%) | +1.9 (+17%) | +2.9 (+14%) | -3.2 (-10%) | -20.8 (-41%) | -0.69 [-0.96, -0.43] |
| M2B + DA size PRS, Huber (kg space) | 0.739 | 5.72 | -1.27 | +1.6 (+28%) | +2.0 (+18%) | +3.0 (+14%) | -3.3 (-11%) | -21.0 (-41%) | -0.76 [-1.03, -0.47] |
| **M2B + DA size PRS, Huber, log space (SHIPPED)** | 0.729 | 5.72 | -1.44 | +1.7 (+31%) | +1.2 (+11%) | +2.2 (+10%) | -3.0 (-10%) | -20.3 (-40%) | -0.76 [-1.04, -0.49] |
| M2B + log current pred, log space | 0.724 | 5.84 | -1.55 | +1.8 (+31%) | +1.5 (+14%) | +2.5 (+12%) | -3.7 (-12%) | -20.9 (-41%) | -0.65 [-0.97, -0.33] |
| M2B + log current pred + DA PRS, log space | 0.726 | 5.74 | -1.43 | +1.7 (+31%) | +1.1 (+9%) | +2.2 (+10%) | -2.9 (-9%) | -20.0 (-39%) | -0.74 [-1.04, -0.46] |

### Calibration by predicted weight (bias kg (%), n per bin)

```
M0            -2.8( -42%) n=106   +1.1( +11%) n=57    +1.1(  +5%) n=114   +0.5(  +2%) n=279  +11.3( +32%) n=21  
P0            -0.6(  -8%) n=6     +5.0( +69%) n=69    +1.3(  +6%) n=323   -3.9( -12%) n=171  +10.7( +31%) n=7   
M1            -0.9( -14%) n=74    +2.6( +32%) n=77    +0.3(  +2%) n=147   -1.1(  -4%) n=265   +9.5( +27%) n=12  
M1Q           -1.8( -27%) n=83    +2.8( +34%) n=68    +1.3(  +6%) n=129   -0.8(  -3%) n=290   +5.9( +16%) n=6   
M1H           -1.0( -15%) n=82    +1.7( +19%) n=81    -1.3(  -6%) n=176   -2.6(  -8%) n=229   +7.6( +21%) n=7   
M1L           -1.8( -28%) n=59    +3.4( +44%) n=91    -2.7( -11%) n=250   -4.1( -13%) n=175         n/a        
M2_mean       -0.1(  -2%) n=45    +3.5( +46%) n=101   -0.3(  -1%) n=131   -0.6(  -2%) n=295         n/a        
M2_akc        -0.7( -11%) n=57    +3.2( +40%) n=96    -0.4(  -2%) n=123   -0.4(  -1%) n=291   +6.5( +17%) n=9   
M2H_akc       -0.6( -10%) n=58    +2.8( +34%) n=101   -1.6(  -7%) n=163   -2.1(  -7%) n=246   +6.3( +16%) n=6   
M2L_akc       -0.4(  -6%) n=106   +0.3(  +2%) n=71    -2.9( -12%) n=156   -1.7(  -6%) n=235   +7.9( +21%) n=7   
M3_prs_z      -1.3( -20%) n=65    +2.5( +30%) n=90    -0.3(  -1%) n=124   -0.8(  -3%) n=285  +10.3( +29%) n=10  
M3_da         -1.2( -18%) n=87    +2.4( +28%) n=70    +0.6(  +3%) n=122   -0.6(  -2%) n=292   +4.7( +12%) n=4   
M3_dense      -1.4( -21%) n=66    +2.5( +31%) n=90    -0.3(  -1%) n=124   -0.8(  -3%) n=284  +10.1( +28%) n=11  
M3_height     -1.9( -29%) n=65    +2.5( +31%) n=88    +0.4(  +2%) n=123   -0.7(  -2%) n=290   +9.1( +26%) n=6   
M3_all        -1.1( -16%) n=86    +2.3( +27%) n=70    +0.4(  +2%) n=130   -0.6(  -2%) n=282   +8.4( +24%) n=6   
M3_blend      -1.1( -17%) n=77    +2.2( +26%) n=79    +0.4(  +2%) n=130   -0.9(  -3%) n=278   +9.5( +27%) n=10  
M3H_blend     -1.2( -18%) n=90    +1.4( +15%) n=76    -1.4(  -6%) n=159   -2.3(  -7%) n=245   +8.5( +24%) n=6   
M2B           -0.2(  -3%) n=68    +1.7( +18%) n=101   -2.5( -10%) n=156   -2.0(  -6%) n=243   +6.6( +17%) n=6   
M3B_blend     -0.7( -10%) n=83    +1.6( +17%) n=80    +0.3(  +1%) n=115   -0.8(  -2%) n=286   +9.1( +25%) n=11  
M3BH_blend    -1.0( -15%) n=90    +1.1( +11%) n=80    -1.7(  -7%) n=144   -2.0(  -7%) n=255   +8.0( +21%) n=6   
M3BH_da       -1.0( -15%) n=95    +1.0( +10%) n=75    -1.5(  -7%) n=154   -2.0(  -7%) n=248   +5.0( +13%) n=4   
M3BL_da       -0.4(  -6%) n=97    +0.0(  +0%) n=82    -3.3( -13%) n=178   -1.4(  -4%) n=207   +7.2( +19%) n=11  
M3BL_blend    -0.6( -10%) n=94    +0.4(  +4%) n=81    -3.2( -13%) n=163   -1.7(  -6%) n=230   +7.8( +20%) n=7   
M3BL_both     -0.4(  -6%) n=100   -0.3(  -2%) n=80    -3.4( -14%) n=176   -1.2(  -4%) n=206   +7.6( +20%) n=12  
```

## Decision

- The current blend is already unbiased on adults; the gain available from any model is ~0.6-0.8 kg MAE (10-12%), all of it from (a) robust fitting that stops 65-77 kg outliers pulling the fit and (b) better handling of mixed-breed dogs. No candidate raises r beyond ~0.74: within-breed body condition is invisible to genetics.
- The LMM `prs_z`, dense prediction and height prediction add no skill on top of the breed model (dMAE vs the recalibrated blend within +-0.1 kg). The Darwin's Ark size PRS adds a small but consistent gain (-0.1 to -0.15 kg) and is kept.
- Plain NNLS breed models (M2*) win on mixed-breed MAE but drive toy breeds to ~0-1 kg (Chihuahua 1.0, Shih Tzu 2.2, Toy Poodle 4.4 kg for a 100% dog) because large "noise" components in small dogs push small-breed weights down. Boxing each breed to 1.5x its prior (M2B family) fixes that at a cost of ~0.1 kg MAE.
- The kg-space stack (M3BH_da, lowest MAE 5.72, r 0.739) still floors toy dogs (a 93% Maltese at 1.45 kg; the current pipeline gives 1.0 kg) because the PRS term is additive in kg. The log-space stack (M3BL_da) has the same MAE 5.72, r 0.729, and predicts that Maltese at 2.6 kg, a 50% Chihuahua mix at 4.3 kg. **M3BL_da is shipped.**
- Cosmo (13.6 kg): current 18.4 kg -> new 17.8 kg. His composition is 37% Standard Poodle, so the model cannot get much lower.

## Shipped model

- `breed_pred = exp(sum_b p_b * log W_b)` over `breed_composition_raw` (unknown label -> 17.8 kg), W_b from box-constrained ridge (lambda 1.0, box 1.5x) in log space; prior = AKC breed-standard weight for 98/231 labels (from `BREED_WEIGHT_KG` in the pipeline), proportion-weighted marginal mean for the rest.
- `pred_kg = exp(0.757 + 0.727 * log(breed_pred) + 0.00025 * da_size_prs)`, Huber c=0.3 in log space.
- Fields kept: `pred_kg`, `pred_lbs`, `prs_z`, `percentile` (still the LMM z); added `pred_kg_breed`, `pred_kg_blend` (old value), `validation`; `method_note` now carries the CV figures. Falls back to the old blend if the model JSON or `breed_result.json` is missing.

## Implied purebred weights (100% of one breed), 30 best-supported labels

```
label                               W_b kg   prior
AMERICAN_PIT_BULL_TERRIER             28.0    27.1
GERMAN_SHEPHERD                       31.9    32.0
SIBERIAN_HUSKY                        24.2    22.0
BOXER                                 33.7    30.0
LABRADOR_RETRIEVER                    42.9    30.0
CHIHUAHUA                              1.7     2.5
WHITE_SWISS_SHEPHERD                  36.1    28.3
VILLAGE_EastAsia                      20.2    24.1
AMERICAN_STAFFORDSHIRE_TERRIER        23.1    26.0
TOY_POODLE                             6.7    10.1
AUSTRALIAN_CATTLE                     20.3    20.0
GOLDEN_RETRIEVER                      33.8    30.0
GRAY_WOLF                             20.0    23.3
YORKSHIRE_TERRIER                      3.8     3.0
SHIH_TZU                               4.1     6.0
BORDER_COLLIE                         21.3    22.6
ROTTWEILER                            43.2    48.0
STANDARD_POODLE                       30.8    25.5
VILLAGE_Oceania                       22.5    24.4
GREAT_PYRENEES                        46.3    45.0
MINIATURE_POODLE                       7.4    11.1
DACHSHUND                              8.2    12.2
MALTESE                                3.2     3.0
COCKER_SPANIEL                        10.3    12.0
ALASKAN_MALAMUTE                      43.9    38.0
AUSTRALIAN_SHEPHERD                   28.7    25.0
POMERANIAN                             3.3     2.5
CONTINENTAL_BULLDOG                   28.4    25.4
CHOW_CHOW                             28.3    28.0
VILLAGE_Americas                      12.9    19.3
```

Labrador (42.9 kg, prior 30) and Chihuahua (1.7 kg, prior 2.5) sit on the box edge: customer Labs and Lab mixes really are heavy, and the additive model wants Chihuahua lower still. The stack's exponent 0.73 and intercept pull both back toward the middle.

## Re-running

```bash
# on a machine with the kit JSONs (Hoffman: $D/results_prosper), weights TSV and the age/weight TSV
python3 analysis/weight_model/fit_weight_model.py <kits_dir> analysis_weights.tsv customer_age_weight.tsv run_dog_pipeline.sh <out_dir> M3BL_da
cp <out_dir>/weight_breed_model.json reference_json/
```

Open follow-ups: 56 weighted kits had no results yet; refit once the ProsperKits batch (1,100 dogs) has weights; a neuter/sex field would likely add more than any further genomic covariate.
