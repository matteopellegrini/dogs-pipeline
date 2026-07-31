"""Scatter: height PRS vs actual AKC breed height, Parker panel."""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

BASE = "/private/tmp/claude-501/-Users-matteopellegrini-Downloads-dogs/0fb26461-2afd-4142-936c-a58295d35959/scratchpad/height"
d = np.load(f"{BASE}/results.npz", allow_pickle=True)
S = json.load(open(f"{BASE}/stats.json"))
prs, height, breeds = d['prs_lobo_z'], d['height'], d['breeds']

BLUE, ORANGE = '#2a78d6', '#eb6834'
INK, INK2, MUTED, GRID = '#1a1a19', '#4a4a46', '#77766e', '#e6e5e0'
SURFACE = '#fcfcfb'

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans'],
    'axes.edgecolor': GRID, 'axes.linewidth': 1.0,
    'xtick.color': MUTED, 'ytick.color': MUTED,
    'xtick.labelsize': 9.5, 'ytick.labelsize': 9.5,
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
})

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.2, 5.9))

# ── Panel A: every dog ───────────────────────────────────────────────────
axA.grid(True, color=GRID, lw=0.8, zorder=0)
axA.set_axisbelow(True)
axA.scatter(height, prs, s=26, c=BLUE, alpha=0.30, linewidths=0, zorder=2)

sl, ic, r, p, _ = stats.linregress(height, prs)
xs = np.linspace(height.min() - 1, height.max() + 1, 100)
# 2px surface ring under the fit so it reads over dense points
axA.plot(xs, sl * xs + ic, color=SURFACE, lw=5.0, zorder=3, solid_capstyle='round')
axA.plot(xs, sl * xs + ic, color=ORANGE, lw=2.4, zorder=4, solid_capstyle='round')

axA.set_title('Every dog in the panel', fontsize=12.5, color=INK, fontweight='bold',
              loc='left', pad=10)
axA.set_xlabel('Actual breed height  (cm, AKC standard)', fontsize=10.5, color=INK2, labelpad=8)
axA.set_ylabel('Height PRS  (z, leave-one-breed-out)', fontsize=10.5, color=INK2, labelpad=8)
axA.text(0.035, 0.955, f"r = {S['r_lobo']:.3f}\nn = {S['n_samples']} dogs",
         transform=axA.transAxes, va='top', ha='left', fontsize=11.5, color=INK,
         fontweight='bold', linespacing=1.5,
         bbox=dict(boxstyle='round,pad=0.55', fc=SURFACE, ec=GRID, lw=1.0))

# ── Panel B: breed means, selectively labelled ───────────────────────────
bm = S['breed_means']                       # code -> [prs_mean, height, n]
codes = list(bm)
bx = np.array([bm[c][1] for c in codes])    # height
by = np.array([bm[c][0] for c in codes])    # mean PRS
bn = np.array([bm[c][2] for c in codes])

axB.grid(True, color=GRID, lw=0.8, zorder=0)
axB.set_axisbelow(True)
axB.scatter(bx, by, s=np.clip(bn * 5, 22, 150), c=BLUE, alpha=0.62,
            linewidths=1.2, edgecolors=SURFACE, zorder=2)

slb, icb, *_ = stats.linregress(bx, by)
xs2 = np.linspace(bx.min() - 1, bx.max() + 1, 100)
axB.plot(xs2, slb * xs2 + icb, color=SURFACE, lw=5.0, zorder=3, solid_capstyle='round')
axB.plot(xs2, slb * xs2 + icb, color=ORANGE, lw=2.4, zorder=4, solid_capstyle='round')

# Label only the height extremes and the largest residuals — never every point.
# Offsets are hand-set because the labelled set is deterministic; this is the
# only way to guarantee no collisions in the tight upper-right cluster.
resid = by - (slb * bx + icb)
show = set(np.argsort(bx)[:4]) | set(np.argsort(bx)[-4:]) | \
       set(np.argsort(np.abs(resid))[-5:])
LABEL = {   # code: (display name, dx, dy, ha)
    'CHIH': ('Chihuahua',          7, -13, 'left'),
    'POM':  ('Pomeranian',         7,   6, 'left'),
    'DACH': ('Dachshund',          7,   3, 'left'),
    'YORK': ('Yorkshire Terrier',  7,   3, 'left'),
    'SHIB': ('Shiba Inu',          7,   3, 'left'),
    'CHOW': ('Chow Chow',          7,   3, 'left'),
    'BULT': ('Bull Terrier',       7,   3, 'left'),
    'HUSK': ('Siberian Husky',    -7,   4, 'right'),
    'BOX':  ('Boxer',              7,   2, 'left'),
    'MAST': ('Mastiff',           -7, -12, 'right'),
    'DANE': ('Great Dane',         7, -12, 'left'),
    'DEER': ('Scottish Deerhound', 7,  -2, 'left'),
    'IWOF': ('Irish Wolfhound',    7,   9, 'left'),
}
for i in sorted(show):
    name, dx, dy, ha = LABEL.get(codes[i], (codes[i], 7, 3, 'left'))
    axB.annotate(name, (bx[i], by[i]), textcoords='offset points',
                 xytext=(dx, dy), ha=ha, fontsize=8.4, color=INK2, zorder=5)

axB.set_xlim(12, 92)   # room for the right-hand labels

axB.set_title('Averaged per breed', fontsize=12.5, color=INK, fontweight='bold',
              loc='left', pad=10)
axB.set_xlabel('Actual breed height  (cm, AKC standard)', fontsize=10.5, color=INK2, labelpad=8)
axB.set_ylabel('Mean height PRS  (z)', fontsize=10.5, color=INK2, labelpad=8)
axB.text(0.035, 0.955, f"r = {S['r_breed']:.3f}\n{S['n_breeds']} breeds",
         transform=axB.transAxes, va='top', ha='left', fontsize=11.5, color=INK,
         fontweight='bold', linespacing=1.5,
         bbox=dict(boxstyle='round,pad=0.55', fc=SURFACE, ec=GRID, lw=1.0))

for ax in (axA, axB):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.suptitle('Genetic height score tracks actual breed height across the Parker panel',
             fontsize=15, color=INK, fontweight='bold', x=0.062, ha='left', y=0.975)
fig.text(0.062, 0.918,
         f"Ridge marginal-effect PRS over {S['n_snps']:,} LD-thinned Parker SNPs, scored under "
         f"leave-one-breed-out cross-validation\n(a breed's own dogs are excluded when fitting the "
         f"weights that score it). Dot size in the right panel = dogs per breed.",
         fontsize=9.8, color=MUTED, ha='left', va='top', linespacing=1.55)

fig.tight_layout(rect=[0.008, 0.01, 0.995, 0.875])
out = f"{BASE}/height_prs_scatter.png"
fig.savefig(out, dpi=200, facecolor=SURFACE)
print("wrote", out)
