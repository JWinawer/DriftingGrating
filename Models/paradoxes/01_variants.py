"""Step 1: every model variant x every imbalance, at the repo's own parameters.

Output: _cache/variants.json -- the four asymmetries (and mean response) per
stimulus set, for each (normalization, imbalance) pair; plus a progress log."""
import json, time, numpy as np
import paradox_lib as L

LOG = L.os.path.join(L.CACHE, "01_variants.log")
open(LOG, "w").close()

E = L.stimulus_energy()
masks = L.wedge_masks()
bar = np.load(L.os.path.join(L.CACHE, "bar_ori.npy"))

# gain conditions: the repo's defaults are cardinal=2, radial=1.5
IMB = {
    "balanced": dict(),
    "cardinal x2": dict(g_h=2, g_v=2),
    "radial x1.5": dict(g_rad=1.5),
    "cardinal x2 + radial x1.5": dict(g_h=2, g_v=2, g_rad=1.5),
}
# the repo's six conditions, at the repo's parameter values
def homog(Eg):
    return {k: e / (0.1 + L.pool_homogeneity(e, L.homogeneity_map(e))) for k, e in Eg.items()}
MODELS = {
    "none (linear)":            lambda Eg: Eg,
    "M1 untuned, p=q=1":        lambda Eg: L.normalize(Eg, L.pool_untuned, 0.01, q=1),
    "M2 untuned, q=2":          lambda Eg: L.normalize(Eg, L.pool_untuned, 0.01, q=2),
    "M3 tuned, p=q=1":          lambda Eg: L.normalize(Eg, L.pool_tuned, 0.01, q=1),
    "M4 anisotropy (global)":   lambda Eg: L.normalize(Eg, L.pool_anisotropy, 0.01),
    "M5 homogeneity-gated q":   homog,
}

# the gain map must reproduce utils_image.imposeImbalance
g = L.channel_gain(g_h=2, g_v=2, g_rad=1.5)
ref = L.utils_image.imposeImbalance({(0, 0): E[0, 0]}, cardinal=2, radial=1.5)[0, 0]
assert np.allclose(ref, E[0, 0] * g[None], rtol=1e-5), "gain map differs from imposeImbalance"

results, n, t0 = {}, len(IMB) * len(MODELS), time.time()
done = 0
for iname, ikw in IMB.items():
    Eg = L.apply_gain(E, L.channel_gain(**ikw))
    for mname, f in MODELS.items():
        R = f(Eg)
        a = L.asymmetries(L.readout(R, masks), bar)
        results[f"{mname} | {iname}"] = a
        done += 1
        el = time.time() - t0
        L.log(f"{done}/{n} ({100*done/n:.0f}%) {mname} | {iname}  elapsed {el:.0f}s  "
              f"remaining ~{el/done*(n-done):.0f}s", LOG)
json.dump(results, open(L.os.path.join(L.CACHE, "variants.json"), "w"), indent=1)
