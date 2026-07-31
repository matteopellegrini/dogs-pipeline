"""
Height PRS vs actual breed height across the Parker panel.

Method mirrors the pipeline's Stage 11 PRS (ridge-regularised marginal effect
sizes on breed-level trait values), with the trait swapped from AKC behavioural
scores to AKC breed height in cm.

Critically this uses LEAVE-ONE-BREED-OUT cross-validation: a breed's own samples
are excluded when fitting the SNP weights used to score them. Scoring a breed
with a model that saw its height is circular and inflates the correlation; the
in-sample version is computed too, purely to show the size of that inflation.
"""
import json, numpy as np
from scipy import stats

BASE = "/private/tmp/claude-501/-Users-matteopellegrini-Downloads-dogs/0fb26461-2afd-4142-936c-a58295d35959/scratchpad/height"
ROOT = "/Users/matteopellegrini/Downloads/dogs/COSMO/analysis"
FAM, BIM, BED = f"{ROOT}/cosmo_parker_full.fam", f"{ROOT}/cosmo_parker_full.bim", f"{ROOT}/cosmo_parker_full.bed"
PRUNE = 5          # keep every 5th SNP (crude LD thinning, as in Stage 11)
LAMBDA_FRAC = 0.1  # ridge strength, as in Stage 11

m = json.load(open(f"{BASE}/mapping.json"))
code_to_breed, height_cm = m['code_to_breed'], m['height_cm']

codes = np.array([l.split()[0] for l in open(FAM)])
n_samples = len(codes)
n_snps = sum(1 for _ in open(BIM))
print(f"Parker panel: {n_samples} samples x {n_snps} SNPs")

# ── vectorised PLINK .bed reader (SNP-major, 2 bits/sample) ──────────────
bps = (n_samples + 3) // 4
raw = np.fromfile(BED, dtype=np.uint8, offset=3).reshape(n_snps, bps)
raw = raw[::PRUNE]                                   # thin before expanding: saves memory
kept = raw.shape[0]
codes2 = np.empty((kept, bps * 4), dtype=np.uint8)
for k in range(4):
    codes2[:, k::4] = (raw >> (2 * k)) & 3
codes2 = codes2[:, :n_samples]
lut = np.array([0.0, np.nan, 1.0, 2.0], dtype=np.float32)
G = lut[codes2]                                      # (snps, samples) copies of A2
del raw, codes2
print(f"Loaded genotypes, LD-thinned to {kept} SNPs")

# ── restrict to samples whose breed has a known AKC height ──────────────
y_all = np.array([height_cm.get(code_to_breed.get(c, ''), np.nan) for c in codes])
keep = ~np.isnan(y_all)
G, y, breeds = G[:, keep], y_all[keep], codes[keep]
print(f"Samples with breed height: {keep.sum()} across {len(set(breeds))} breed codes")

# mean-impute missing genotypes per SNP, then drop monomorphic SNPs
row_mean = np.nanmean(G, axis=1, keepdims=True)
G = np.where(np.isnan(G), row_mean, G)
# float64 throughout: in float32 the effect-size division overflows for
# near-monomorphic SNPs and silently poisons the PRS.
G = G.astype(np.float64)
poly = G.std(axis=1) > 1e-8
G = G[poly]
y = y.astype(np.float64)
print(f"Polymorphic SNPs used: {G.shape[0]}")


def marginal_beta(Gtr, ytr):
    """Ridge-regularised marginal (per-SNP) effect sizes — same form as Stage 11."""
    Gc = Gtr - Gtr.mean(axis=1, keepdims=True)
    yc = ytr - ytr.mean()
    var_j = np.einsum('ij,ij->i', Gc, Gc)
    denom = var_j + LAMBDA_FRAC * var_j.mean()
    beta = np.zeros_like(var_j)
    ok = denom > 0                      # SNP monomorphic within this training split
    beta[ok] = (Gc[ok] @ yc) / denom[ok]
    return beta


# ── in-sample (circular) ────────────────────────────────────────────────
beta_full = marginal_beta(G, y)
prs_insample = G.T @ beta_full

# ── leave-one-breed-out ─────────────────────────────────────────────────
prs_lobo = np.full(len(y), np.nan)
uniq = np.unique(breeds)
for i, b in enumerate(uniq):
    test = breeds == b
    beta = marginal_beta(G[:, ~test], y[~test])
    prs_lobo[test] = G[:, test].T @ beta
    if (i + 1) % 25 == 0:
        print(f"  LOBO {i+1}/{len(uniq)} breeds")

# standardise for readability (z across panel)
z = lambda v: (v - np.nanmean(v)) / np.nanstd(v)
prs_lobo_z, prs_in_z = z(prs_lobo), z(prs_insample)

r_lobo = stats.pearsonr(prs_lobo_z, y)
r_in = stats.pearsonr(prs_in_z, y)
rho_lobo = stats.spearmanr(prs_lobo_z, y)

# breed-level means (one point per breed)
bmean = {}
for b in uniq:
    s = breeds == b
    bmean[b] = (float(np.mean(prs_lobo_z[s])), float(y[s][0]), int(s.sum()))
bx = np.array([v[0] for v in bmean.values()]); by = np.array([v[1] for v in bmean.values()])
r_breed = stats.pearsonr(bx, by)

print("\n─── Height PRS vs actual breed height ───")
print(f"  per-sample, leave-one-breed-out : r = {r_lobo[0]:.3f}  (p = {r_lobo[1]:.2e}, n = {len(y)})")
print(f"  per-sample, Spearman (LOBO)     : rho = {rho_lobo[0]:.3f}")
print(f"  per-breed mean, LOBO            : r = {r_breed[0]:.3f}  (n = {len(bx)} breeds)")
print(f"  per-sample, IN-SAMPLE (circular): r = {r_in[0]:.3f}   <- inflated, for contrast")

np.savez(f"{BASE}/results.npz", prs_lobo_z=prs_lobo_z, prs_in_z=prs_in_z,
         height=y, breeds=breeds)
json.dump({'r_lobo': r_lobo[0], 'p_lobo': r_lobo[1], 'rho_lobo': rho_lobo[0],
           'r_breed': r_breed[0], 'r_insample': r_in[0],
           'n_samples': int(len(y)), 'n_breeds': int(len(uniq)), 'n_snps': int(G.shape[0]),
           'breed_means': {k: v for k, v in bmean.items()}},
          open(f"{BASE}/stats.json", 'w'), indent=1)
print(f"\nSaved -> {BASE}/results.npz, stats.json")
