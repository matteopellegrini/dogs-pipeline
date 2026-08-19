# Darwin's Ark individual-level PRS

Training/build scripts for the `reference_json/darwins_ark/` artifact.
They run from a data directory (untracked `darwins_ark/`) holding the
Dryad download (doi:10.5061/dryad.g4f4qrfr0, CC0): the 8.5M-SNP PLINK set
(canFam3.1), survey CSVs, the UCSC canFam3ToCanFam4 chain, and GEMMA
outputs. Order: liftover_bim -> build_gemma_inputs -> run_queue (GEMMA
sweep) -> evaluate_cv / dump_cv_preds (held-out validation; p<=0.1 chosen)
-> filter_and_rescore (platform-site filter + honest references) ->
build_artifact (manifest + weights). score_dog.py / score_cohort.sh score
a dog or cohort from its GLIMPSE BCF. breed_benchmark.py validates the
230-breed panel against Darwin's Ark published ancestry.

Key numbers (held-out, dog-level CV): size r=0.77; biddability 0.24,
human sociability 0.21, down to agonistic threshold 0.08 — tracking h2.
Weight blend (dense breed-level + DA size PRS): r=0.92, MAE 5.1 kg on 94
owner-weighed dogs. z_scale corrects in-sample reference inflation
(fold-based 0.796 for size == 0.795 from 96 external dogs).
