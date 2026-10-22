"""Step 2: how robust is each result? Sweeps sigma, gain sizes, pool size,
exponent, M4's spatial window, M5's construction, and how 'overrepresentation'
enters the model:

  input  -- the repo's choice: the gain multiplies energy before any exponent,
            so the pool sees (g*E)**q. Overrepresentation acts like extra
            effective contrast for those channels.
  count  -- more neurons, each unchanged: numerator g*E**p, pool built from
            g*E**q. Overrepresentation adds units, not contrast.
  (Identical when p = q = 1, so only M2 and M5 are run both ways.)

Usage: python 02_sweeps.py [quick]   (quick = one gain condition, for timing)
Writes _cache/sweeps.json and _cache/02_sweeps.log (progress + time remaining)."""
import sys, json, time, itertools, numpy as np
import paradox_lib as L

LOG = L.os.path.join(L.CACHE, "02_sweeps.log")
open(LOG, "w").close()
QUICK = len(sys.argv) > 1 and sys.argv[1] == "quick"

E = L.stimulus_energy()
masks = L.wedge_masks()
band = L.band_index()
masks_band = masks[:, band]
bar = np.load(L.os.path.join(L.CACHE, "bar_ori.npy"))
SIGMAS = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]
H_cache = {}

def H_of(key, e, s_deg):
    k = (key, s_deg, e.ctypes.data)
    if k not in H_cache:
        H_cache[k] = L.homogeneity_map(e, sigma_s_deg=s_deg)
    return H_cache[k]

# --- each model: returns (num, pool) on band pixels for one stimulus ----------------
def make(kind, mode="input", **kw):
    def f(key, E0, g):
        Eg = E0 * g[None]                       # gained energy
        q = kw.get("q", 1.0)
        if mode == "input":
            num, src = Eg, Eg
        else:   # count: numerator g*E, pool from g*E**q (passed as (g*E**q)**(1/q))
            num, src = Eg, (g[None] * E0 ** q) ** (1.0 / q)
        if kind == "none":
            pool = np.zeros_like(num)
        elif kind == "untuned":
            pool = L.pool_untuned(src, q=q, sigma_xy_deg=kw.get("xy", 3.0))
        elif kind == "tuned":
            pool = L.pool_tuned(src, q=q, sigma_xy_deg=kw.get("xy", 3.0))
        elif kind == "aniso":
            pool = L.pool_anisotropy_local(src, kw.get("window"))
        elif kind == "homog":
            # homogeneity is read from the gained energy, as in the repo
            H = H_of(key, Eg, kw.get("surround", 10.0))
            if kw.get("pool", "local") == "wrap":
                pool = L.pool_homogeneity_wrap(src if mode == "input" else Eg, H)
            else:
                pool = L.pool_homogeneity(src if mode == "input" else Eg, H,
                                          pool=kw.get("pool", "local"))
            if mode == "count":          # pool of g*E**q with the per-pixel q
                qmap = 1 + H
                X = g[None] * E0 ** qmap[None, None]
                from scipy.ndimage import gaussian_filter
                if kw.get("pool", "local") == "local":
                    pool = gaussian_filter(X, sigma=(1, 1, 1, 1))
                elif kw.get("pool") == "wrap":
                    pool = gaussian_filter(X, sigma=(1, 1, 1, 1),
                                           mode=["reflect", "wrap", "reflect", "reflect"])
                else:
                    pool = np.broadcast_to(L.blur_xy(X.mean(axis=(0, 1), keepdims=True), 3.0), X.shape)
        if mode == "count" and kind in ("untuned", "tuned", "aniso"):
            pass                             # p = 1: numerator already g*E
        return num[..., band], np.asarray(pool)[..., band]
    return f

BASE = {
    "M1 untuned":                 make("untuned"),
    "M2 untuned q=2 [input]":     make("untuned", q=2),
    "M2 untuned q=2 [count]":     make("untuned", "count", q=2),
    "M3 tuned":                   make("tuned"),
    "M4 anisotropy global":       make("aniso"),
    "M5 homogeneity [input]":     make("homog"),
    "M5 homogeneity [count]":     make("homog", "count"),
    "linear":                     make("none"),
}
EXTRA = {}
for xy in [1, 6, 12]:
    EXTRA[f"M1 untuned pool {xy}deg"] = make("untuned", xy=xy)
for q, xy in itertools.product([1.25, 1.5, 2], [1, 3, 12]):
    EXTRA[f"M2 untuned q={q} pool {xy}deg [input]"] = make("untuned", q=q, xy=xy)
for w in [1, 3, 6, 12]:
    EXTRA[f"M4 anisotropy window {w}deg"] = make("aniso", window=w)
for s, p in itertools.product([3, 10], ["local", "wrap", "untuned"]):
    EXTRA[f"M5 homogeneity surround {s}deg pool {p} [input]"] = make("homog", surround=s, pool=p)

GAINSETS = {"repo (card 2, rad 1.5)": (2, 1.5),
            "behaviour (card 1.3, rad 1.1)": (1.3, 1.1),
            "macaque V1 (card 1.11, rad 1.15)": (1.11, 1.15)}
GAINS = {"balanced": dict()}
for gname, (c, r) in GAINSETS.items():
    GAINS[f"{gname} | cardinal only"] = dict(g_h=c, g_v=c)
    GAINS[f"{gname} | radial only"] = dict(g_rad=r)
    GAINS[f"{gname} | both"] = dict(g_h=c, g_v=c, g_rad=r)

# radial bias outside the pool: cardinal gain enters before normalization,
# radial gain multiplies the normalized output
RADOUT = {
    "M1 untuned [radial on output]": "M1 untuned",
    "M2 untuned q=2 [input] [radial on output]": "M2 untuned q=2 [input]",
    "M2 untuned q=2 pool 1deg [input] [radial on output]": "M2 untuned q=2 pool 1deg [input]",
    "M4 anisotropy global [radial on output]": "M4 anisotropy global",
    "M4 anisotropy window 3deg [radial on output]": "M4 anisotropy window 3deg",
    "M5 homogeneity surround 10deg pool wrap [input] [radial on output]":
        "M5 homogeneity surround 10deg pool wrap [input]",
}

jobs = [(m, gk) for gk in GAINS for m in BASE]
jobs += [(m, gk) for gk in GAINS if gk == "balanced" or gk.startswith("repo") for m in EXTRA]
if QUICK:
    jobs = [(m, gk) for (m, gk) in jobs if gk == "repo (card 2, rad 1.5) | both"]
jobs += [(m, gk) for gk in GAINS if gk.endswith("| both") for m in RADOUT]
MODELS = {**BASE, **EXTRA}
gain_maps = {}
results = {}
pool_scale = {}     # median pool value in the 4-8 deg band, balanced population
t0 = time.time()
for n, (mname, gk) in enumerate(jobs, 1):
    kw = dict(GAINS[gk])
    out_gain = None
    if mname in RADOUT:
        g_rad = kw.pop("g_rad", 1.0)
        out_gain = L.channel_gain(g_rad=g_rad)[:, band]
    key = (gk, mname in RADOUT)
    if key not in gain_maps:
        gain_maps = {key: L.channel_gain(**kw)}           # keep one map in memory
    g = gain_maps[key]
    resp = np.zeros((len(SIGMAS), 2, 4, 8))
    for (i, j), e in E.items():
        num, pool = MODELS[RADOUT.get(mname, mname)]((i, j), e, g)
        if mname == "linear":
            resp[:, i, j] = L.readout_band(num, pool + 1.0, [0.0], masks_band)[0]
        else:
            resp[:, i, j] = L.readout_band(num, pool, SIGMAS, masks_band, out_gain)
        if gk == "balanced" and mname not in RADOUT:
            pool_scale.setdefault(mname, []).append(float(np.median(pool)))
    results[f"{mname} || {gk}"] = [L.asymmetries(resp[s], bar) for s in range(len(SIGMAS))]
    el = time.time() - t0
    L.log(f"{n}/{len(jobs)} ({100*n/len(jobs):.0f}%) elapsed {el/60:.1f} min, "
          f"remaining ~{el/n*(len(jobs)-n)/60:.1f} min :: {mname} || {gk}", LOG)
    if n % 20 == 0 or n == len(jobs):
        json.dump({"sigmas": SIGMAS, "results": results, "pool_scale": pool_scale},
                  open(L.os.path.join(L.CACHE, "sweeps_quick.json" if QUICK else "sweeps.json"), "w"))
