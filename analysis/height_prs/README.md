# Height PRS vs actual breed height (Parker panel)

Does a genetic height score predict real breed height? Answer: **r = 0.732** per dog
(1,326 dogs, 162 breed codes), **r = 0.750** on breed means.

## Run order

```bash
python3 map_breeds.py     # Parker breed codes -> AKC breed names + heights -> mapping.json
python3 height_prs.py     # PRS with leave-one-breed-out CV -> results.npz, stats.json
python3 plot_height.py    # -> height_prs_scatter.png
```

Requires the interpreter with the data-science stack (`/usr/bin/python3`, not the
genomics micromamba env — that one lacks pandas/scipy).

`map_breeds.py` downloads AKC heights from `tmfilho/akcdata` and scrapes the
`PARKER_TO_AKC` table straight out of `run_dog_pipeline.sh`, so it stays in sync
with the pipeline rather than duplicating 150 breed-name mappings.

## Method

Same as pipeline Stage 11 — ridge-regularised marginal (per-SNP) effect sizes over
28,787 LD-thinned Parker SNPs — with the target swapped from AKC behavioural
ratings to AKC breed height in cm.

**Leave-one-breed-out cross-validation**: a breed's own dogs are excluded when
fitting the weights that score them. This matters a lot:

| Scoring | r |
|---|---|
| Leave-one-breed-out (honest) | 0.732 |
| In-sample (circular) | 0.912 |

Computation is float64 throughout — in float32 the effect-size division overflows
on near-monomorphic SNPs and silently corrupts the PRS.

## Caveats

- Height is the **AKC breed-standard midpoint, sexes combined**, not measured per
  dog. All dogs of a breed share one x value, so within-breed variation is invisible.
- 162 of 178 Parker codes used. Excluded: wolves, Italian/Chinese landraces, COSMO
  itself, and Manchester Terrier (code doesn't distinguish Standard from Toy, which
  differ ~10 cm). `FOXH` assumed American Foxhound (vs English, ~2 cm apart).
- SPOO/MPOO/TPOO are mapped to Standard/Miniature/Toy Poodle **separately**. The
  pipeline's own table collapses all three to "Poodles" — fine for temperament,
  wrong for height.
- Marginal effects on a panel where breed *is* the population structure partly
  capture ancestry, so this is a breed-discriminating score as much as a causal
  height score. Dog body size being unusually oligogenic (*IGF1*, *IGFBP4*,
  *HMGA2*, *SMAD2*) is likely why a simple model does this well.

## Notable outliers

- **Bull Terrier** — 55 cm but PRS −2.04, scores like a toy breed.
- **Boxer, Siberian Husky, Chow Chow, Shiba Inu** — score well above their height.
