"""Step 3: figures for Models/PARADOXES.md, from _cache/sweeps.json.

fig1_sign_table   -- each variant x each test, at the repo's parameters
fig2_sigma        -- every asymmetry against normalization strength
fig3_pool_size    -- how the context effect depends on spatial pool size
fig4_mechanism    -- (A) why q > p is needed to invert, (B) why a broad pool
                     sees a polar grating as orientation-balanced
"""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import paradox_lib as L

OUT = L.os.path.join(L.HERE, "figures")
L.os.makedirs(OUT, exist_ok=True)
d = json.load(open(L.os.path.join(L.CACHE, "sweeps.json")))
S, R = d["sigmas"], d["results"]
REPO_GAIN = "repo (card 2, rad 1.5) | both"

# V1, matched 7 observers (Reproduction/RESULTS.md), % signal change
DATA = {"dg": {"co": -0.193, "rt": 0.091, "hv": -0.514},
        "da": {"co": -0.044, "rt": 0.198, "hv": -0.212}}


def rel(a, s, k):
    """Asymmetry as % of the model's mean response over both stimulus sets."""
    m = (a["dg"]["mean"] + a["da"]["mean"]) / 2
    return 100 * a[s][k] / m


def tests(a):
    """The six numbers the three paradoxes turn on, in data-sign convention."""
    return {
        "dg card-obl": rel(a, "dg", "co"),
        "da card-obl": rel(a, "da", "co"),
        "dg rad-tang": rel(a, "dg", "rt"),
        "da rad-tang": rel(a, "da", "rt"),
        "context card-obl\n(dg - da)": rel(a, "dg", "co") - rel(a, "da", "co"),
        "context rad-tang\n(dg - da)": rel(a, "dg", "rt") - rel(a, "da", "rt"),
    }


DATA_T = {"dg card-obl": DATA["dg"]["co"], "da card-obl": DATA["da"]["co"],
          "dg rad-tang": DATA["dg"]["rt"], "da rad-tang": DATA["da"]["rt"],
          "context card-obl\n(dg - da)": DATA["dg"]["co"] - DATA["da"]["co"],
          "context rad-tang\n(dg - da)": DATA["dg"]["rt"] - DATA["da"]["rt"]}

# ---- fig 1 --------------------------------------------------------------------------
ROWS = [
    ("linear (no normalization)", "linear"),
    ("M1 untuned, p=q=1", "M1 untuned"),
    ("M2 untuned, q=2", "M2 untuned q=2 [input]"),
    ("M2, q=2, 12-deg pool", "M2 untuned q=2 pool 12deg [input]"),
    ("M3 tuned, p=q=1", "M3 tuned"),
    ("M4 anisotropy, global", "M4 anisotropy global"),
    ("M4 anisotropy, 3-deg window", "M4 anisotropy window 3deg"),
    ("M5 homogeneity (repo)", "M5 homogeneity [input]"),
    ("M5 homogeneity, circular ori", "M5 homogeneity surround 10deg pool wrap [input]"),
    ("M2 + radial on output", "M2 untuned q=2 [input] [radial on output]"),
    ("M4 global + radial on output", "M4 anisotropy global [radial on output]"),
    ("M4 3-deg window + radial on output", "M4 anisotropy window 3deg [radial on output]"),
    ("M5 circular + radial on output", "M5 homogeneity surround 10deg pool wrap [input] [radial on output]"),
]
REPO_SIGMA = {"M5": 0.1}            # every other model uses 0.01 in the repo


def repo_sigma(key):
    return REPO_SIGMA.get(key[:2], 0.01)


cols = list(DATA_T)
M = np.array([[tests(R[f"{k} || {REPO_GAIN}"][S.index(repo_sigma(k))])[c] for c in cols]
              for _, k in ROWS])
fig, ax = plt.subplots(figsize=(10.5, 6.8))
lim = np.nanmax(np.abs(M))
ax.imshow(M, cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
for r in range(M.shape[0]):
    for c in range(M.shape[1]):
        ok = np.sign(M[r, c]) == np.sign(DATA_T[cols[c]]) and abs(M[r, c]) > 0.5
        ax.text(c, r, f"{M[r, c]:+.0f}" + (" ✓" if ok else ""), ha="center", va="center",
                fontsize=9, fontweight="bold" if ok else "normal")
ax.set_xticks(range(len(cols)))
ax.set_xticklabels([f"{c}\ndata {DATA_T[c]:+.3f}" for c in cols], fontsize=8)
ax.set_yticks(range(len(ROWS)))
ax.set_yticklabels([r for r, _ in ROWS], fontsize=9)
ax.set_title("Model asymmetries (% of mean response) at the repo's parameters\n"
             "cardinal gain 2, radial gain 1.5, 100% contrast.  ✓ = same sign as the data",
             fontsize=10)
fig.tight_layout()
fig.savefig(L.os.path.join(OUT, "fig1_sign_table.png"), dpi=150)

# ---- fig 2: sigma ---------------------------------------------------------------------
PANEL = [("M1 untuned", "M1 untuned, p=q=1"), ("M2 untuned q=2 [input]", "M2 untuned, q=2"),
         ("M3 tuned", "M3 tuned"), ("M4 anisotropy global", "M4 anisotropy, global"),
         ("M5 homogeneity surround 10deg pool wrap [input]", "M5 homogeneity, circular ori"),
         ("M4 anisotropy window 3deg [radial on output]", "M4 3-deg window + radial on output")]
fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
for ax, (k, title) in zip(axes.flat, PANEL):
    lst = R[f"{k} || {REPO_GAIN}"]
    for s, key, col, ls in [("dg", "co", "C0", "-"), ("da", "co", "C0", "--"),
                            ("dg", "rt", "C3", "-"), ("da", "rt", "C3", "--")]:
        ax.plot(S, [rel(a, s, key) for a in lst], ls, color=col, marker="o", ms=3,
                label=f"{s} {'card-obl' if key == 'co' else 'rad-tang'}")
    rs = repo_sigma(k)
    q = 2 if "q=2" in k else (1.5 if k.startswith("M5") else 1)
    ax.axvline(rs, color="k", lw=0.8)
    # energy scales as contrast^2, so at contrast c the pool term scales as c^(2q):
    # sigma at 50% contrast is equivalent to sigma / 0.5^(2q) at 100%
    ax.axvline(rs * 0.5 ** -(2 * q), color="k", lw=0.8, ls=":")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xscale("log"); ax.set_title(title, fontsize=10)
    ax.set_ylabel("% of mean response")
for ax in axes[1]:
    ax.set_xlabel("sigma (semi-saturation constant), 100%-contrast units")
axes[0, 0].legend(fontsize=8)
fig.suptitle("Signs depend on normalization strength. Solid line: repo sigma at 100% contrast. "
             "Dotted: the same sigma at the experiment's 50% contrast.\n"
             "Large sigma = weak normalization; every model then converges on the linear model.",
             fontsize=10)
fig.tight_layout()
fig.savefig(L.os.path.join(OUT, "fig2_sigma.png"), dpi=150)

# ---- fig 3: pool size ----------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
fam = {"M2 untuned, q=2": [(1, "M2 untuned q=2 pool 1deg [input]"),
                           (3, "M2 untuned q=2 pool 3deg [input]"),
                           (12, "M2 untuned q=2 pool 12deg [input]")],
       "M4 anisotropy": [(1, "M4 anisotropy window 1deg"), (3, "M4 anisotropy window 3deg"),
                         (6, "M4 anisotropy window 6deg"), (12, "M4 anisotropy window 12deg"),
                         (30, "M4 anisotropy global")]}
for ax, (title, items) in zip(axes, fam.items()):
    xs = [x for x, _ in items]
    for s, key, col, ls in [("dg", "co", "C0", "-"), ("da", "co", "C0", "--"),
                            ("dg", "rt", "C3", "-"), ("da", "rt", "C3", "--")]:
        ax.plot(xs, [rel(R[f"{k} || {REPO_GAIN}"][S.index(0.01)], s, key) for _, k in items],
                ls, color=col, marker="o", label=f"{s} {'card-obl' if key == 'co' else 'rad-tang'}")
    ax.axhline(0, color="0.6", lw=0.8); ax.set_xscale("log")
    ax.set_xticks(xs); ax.set_xticklabels([("global" if x == 30 else f"{x}") for x in xs])
    ax.set_xlabel("spatial pool / window size (deg, Gaussian sigma)")
    ax.set_ylabel("% of mean response"); ax.set_title(title + "  (sigma = 0.01)", fontsize=10)
axes[0].legend(fontsize=8)
fig.suptitle("The cardinal context effect needs a pool wider than the distance over which a polar "
             "grating's orientation changes", fontsize=10)
fig.tight_layout()
fig.savefig(L.os.path.join(OUT, "fig3_pool_size.png"), dpi=150)

# ---- fig 4: mechanism ----------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
g = np.linspace(1, 3, 200)
k, s = 1.0, 1.0
for q, lab, col in [(1, "p = q = 1", "C0"),
                    (2, "p = 1, q = 2, gain before exponent", "C3")]:
    r = g / (s + k * g ** q)
    axes[0].plot(g, r / r[0], color=col, label=lab)
axes[0].axhline(1, color="0.6", lw=0.8)
axes[0].set_xlabel("overrepresentation gain g of the driven channel")
axes[0].set_ylabel("response relative to g = 1")
axes[0].set_title("A. Only suppression that grows faster than\nexcitation can turn more into less",
                  fontsize=10)
axes[0].legend(fontsize=8)

E = L.stimulus_energy()
theta, ecc = L.angle_maps()
pt = (ecc > 5.8) & (ecc < 6.2) & (np.abs((theta + 180) % 360 - 180) < 2)   # 6 deg, right HM
sizes = [0.5, 1, 2, 3, 6, 12]
for (i, j), lab, col in [((0, 2), "Cartesian horizontal grating", "C0"),
                         ((1, 0), "pinwheel (locally horizontal here)", "C3")]:
    e = E[i, j].mean(axis=0)                       # (4, H, W)
    frac = []
    for sd in sizes:
        b = L.blur_xy(e, sd)
        frac.append(float((b[2][pt] / b.sum(axis=0)[pt]).mean()))
    axes[1].plot(sizes, frac, "o-", color=col, label=lab)
axes[1].axhline(0.25, color="0.6", lw=0.8, ls=":")
axes[1].text(0.6, 0.27, "orientation-balanced (1/4)", fontsize=8, color="0.4")
axes[1].set_xscale("log"); axes[1].set_xticks(sizes); axes[1].set_xticklabels(sizes)
axes[1].set_xlabel("pool size (deg, Gaussian sigma)")
axes[1].set_ylabel("share of pooled energy in the horizontal channel")
axes[1].set_title("B. A broad pool at 6 deg on the horizontal meridian\n"
                  "sees a polar grating as nearly orientation-balanced", fontsize=10)
axes[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(L.os.path.join(OUT, "fig4_mechanism.png"), dpi=150)
print("figures written to", OUT)
