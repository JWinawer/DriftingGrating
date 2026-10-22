"""Shared machinery for the three-paradox analysis (see Models/PARADOXES.md).

Builds the notebook's stimuli and steerable-pyramid energies exactly as
model_DriftingGratings.ipynb does, then reads every model variant out the way
the fMRI analysis reads V1: one number per stimulus per 45-degree polar-angle
wedge, 4-8 deg eccentricity, and the four asymmetries formed from those.

Nothing in helper_functions/ is modified. The existing normalization functions
are called as they are; the only new model code here is (a) an imbalance that
lets horizontal and vertical gains differ, and (b) a parameterised version of
the untuned / tuned / homogeneity-gated pools, so sigma, pool size and the
exponents can be swept. Each re-implementation is checked against the
original function at the notebook's parameter values (check_against_originals).
"""
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from helper_functions import pyrmodel, utils_image  # noqa: E402

# ---- geometry, copied from the notebook (cell 5) ---------------------------------
IMG_SIZE = 896
EXP_STIM_DIAM_DEG = 24.4
_x = int(np.round(IMG_SIZE * 0.85)); _x = _x if _x % 2 == 0 else _x + 1
PIXPDEG = _x / EXP_STIM_DIAM_DEG
PADDED_DEG = IMG_SIZE / PIXPDEG
N_SF, N_ORI = 6, 4

# Pyramid orientation channel index -> preferred BAR orientation, conventional
# degrees (0 = horizontal, counter-clockwise). Measured, not assumed: see
# check_channel_orientations().
CH_BAR_ORI = np.array([90, 45, 0, 135])

# Notebook stimulus order: oris = [0, 135, 90, 45] (the grating's *frequency*
# direction). Cartesian set: bar orientation is frequency + 90.
CART_NAMES = ["V", "D45", "H", "D135"]           # D45 = '/', D135 = '\'
CART_BAR_ORI = np.array([90, 45, 0, 135])
POLAR_NAMES = ["pin", "spiralA", "ann", "spiralB"]  # handedness measured below

WEDGES = np.arange(0, 360, 45)                    # conventional polar angle

CACHE = os.path.join(HERE, "_cache")


def log(msg, logfile=None):
    line = time.strftime("%H:%M:%S ") + msg
    print(line, flush=True)
    if logfile:
        with open(logfile, "a") as f:
            f.write(line + "\n")


# ---- stimuli and energies ----------------------------------------------------------
def make_model():
    return pyrmodel.SteerablePyramidSF(n_ori=N_ORI, n_sf=N_SF, img_size=IMG_SIZE,
                                       img_degrees=PADDED_DEG, is_complex=True,
                                       downsample=False, tight_frame=True, device="cpu")


def make_stimuli(contrast=1.0):
    """(2, 4, H, W): set 0 = Cartesian, set 1 = polar, in the notebook's order."""
    import contextlib, io
    out = []
    for ref in ["cartesian", "polar"]:
        with contextlib.redirect_stdout(io.StringIO()):   # the class prints SFs
            s = pyrmodel.InputStimuli(sfs=np.array([1]), oris=np.array([0, 135, 90, 45]),
                                      img_size=IMG_SIZE, img_degrees=PADDED_DEG,
                                      referenceframe=ref, mask=True,
                                      paddedBondaryPerc=EXP_STIM_DIAM_DEG / PADDED_DEG)
        out.append(np.stack([s.input[(0, j)] for j in range(4)]))
    return contrast * np.stack(out)


def stimulus_energy(stims=None, use_cache=True):
    """dict (set, stim) -> (n_sf, n_ori, H, W) energy, as notebook cell 15."""
    path = os.path.join(CACHE, "energy.npz")
    if use_cache and stims is None and os.path.exists(path):
        z = np.load(path)
        return {(i, j): z[f"e{i}{j}"] for i in range(2) for j in range(4)}
    if stims is None:
        stims = make_stimuli()
    m = make_model()
    E = {}
    for i in range(2):
        for j in range(4):
            E[i, j] = m.get_energy_of_filters(m.forward(stims[i, j]), to_magnitude=False,
                                              to_numpy=True, reshape=True, normalize=False)
    if use_cache:
        os.makedirs(CACHE, exist_ok=True)
        np.savez(path, **{f"e{i}{j}": E[i, j].astype(np.float32) for (i, j) in E})
    return E


# ---- polar-angle geometry -----------------------------------------------------------
def angle_maps():
    """Conventional polar angle (deg, y up) and eccentricity (deg) per pixel.
    Same centre convention as utils_image.get_stats_within_eccentricity_by_angle."""
    y, x = np.ogrid[:IMG_SIZE, :IMG_SIZE]
    c = IMG_SIZE / 2
    theta = (np.degrees(np.arctan2(-(y - c), (x - c))) + 360) % 360
    ecc = np.sqrt((x - c) ** 2 + (y - c) ** 2) / PIXPDEG
    return np.broadcast_to(theta, (IMG_SIZE, IMG_SIZE)), np.broadcast_to(ecc, (IMG_SIZE, IMG_SIZE))


def wedge_masks(ecc_lo=4, ecc_hi=8):
    theta, ecc = angle_maps()
    band = (ecc >= ecc_lo) & (ecc <= ecc_hi)
    masks = []
    for w in WEDGES:
        d = np.abs((theta - w + 180) % 360 - 180)
        masks.append(band & (d < 22.5))
    return np.stack(masks)


def local_bar_orientation(E, masks):
    """Dominant-channel bar orientation per stimulus x wedge, from the SF channel
    closest to 1 cpd (index 3). Used to label the polar stimuli empirically."""
    out = np.zeros((2, 4, len(WEDGES)))
    for (i, j), e in E.items():
        for k, m in enumerate(masks):
            out[i, j, k] = CH_BAR_ORI[np.argmax(e[3][:, m].mean(axis=1))]
    return out


# ---- readout: what a V1 voxel at each wedge reports ---------------------------------
def readout(R, masks):
    """(2, 4, 8): mean over SF, sum over orientation channels, mean over the wedge.
    This is the notebook's readout (condense_energy mean over SF, then
    plot_energy_per_stim sums over orientation)."""
    out = np.zeros((2, 4, len(WEDGES)))
    for (i, j), r in R.items():
        img = r.mean(axis=0).sum(axis=0)
        for k, m in enumerate(masks):
            out[i, j, k] = img[m].mean()
    return out


def classify(bar_ori, wedge):
    """Local-orientation labels for a bar orientation at a polar angle."""
    rel = (bar_ori - wedge) % 180           # 0 = radial, 90 = tangential
    return {
        "H": bar_ori % 180 == 0, "V": bar_ori % 180 == 90,
        "card": bar_ori % 90 == 0, "obl": bar_ori % 90 == 45,
        "rad": rel == 0, "tang": rel == 90,
        "polc": rel % 90 == 0, "polo": rel % 90 == 45,
    }


def asymmetries(resp, bar_ori):
    """Four asymmetries per stimulus set, equal weight over the 8 wedges.

    resp, bar_ori: (2, 4, 8). Within a wedge each asymmetry is the mean of the
    stimuli in one class minus the mean of those in the other; then averaged
    over wedges -- the equal-coverage definition used in Reproduction/.
    Also returns the mean response (the 'intercept')."""
    out = {}
    for i, setname in enumerate(["dg", "da"]):
        acc = {k: [] for k in ["hv", "co", "rt", "pp"]}
        for k, w in enumerate(WEDGES):
            labs = [classify(bar_ori[i, j, k], w) for j in range(4)]
            r = resp[i, :, k]
            def m(key):
                sel = [r[j] for j in range(4) if labs[j][key]]
                return np.mean(sel) if sel else np.nan
            acc["hv"].append(m("H") - m("V"))
            acc["co"].append(m("card") - m("obl"))
            acc["rt"].append(m("rad") - m("tang"))
            acc["pp"].append(m("polc") - m("polo"))
        out[setname] = {k: float(np.nanmean(v)) for k, v in acc.items()}
        out[setname]["mean"] = float(resp[i].mean())
    return out


# ---- orientation imbalance ----------------------------------------------------------
def radial_weight():
    """(4, H, W) weight in [0,1]: 1 where a channel's bar orientation is radial,
    0 where tangential, linear in between -- utils_image.orientation_weight's
    ramp, written in the conventional frame. Centre at (N-1)/2, as in
    imposeImbalance (the readout masks use N/2, as get_stats_* does)."""
    y, x = np.ogrid[:IMG_SIZE, :IMG_SIZE]
    c = (IMG_SIZE - 1) / 2
    theta = (np.degrees(np.arctan2(-(y - c), (x - c))) + 360) % 360
    w = np.zeros((N_ORI, IMG_SIZE, IMG_SIZE))
    for c, ori in enumerate(CH_BAR_ORI):
        d = np.abs((theta - ori + 90) % 180 - 90)      # 0..90, distance from radial
        w[c] = 1 - d / 90
    return w


def channel_gain(g_h=1.0, g_v=1.0, g_rad=1.0, radial_profile="linear"):
    """(4, H, W) multiplicative gain per orientation channel.
    g_h, g_v: horizontal and vertical channel gains (utils_image.imposeImbalance
    uses one 'cardinal' value for both). g_rad: peak gain where the channel is
    radial, falling to 1 at tangential."""
    gain = np.ones((N_ORI, IMG_SIZE, IMG_SIZE))
    gain[CH_BAR_ORI == 0] *= g_h
    gain[CH_BAR_ORI == 90] *= g_v
    if g_rad != 1:
        w = radial_weight()
        if radial_profile == "vonmises":
            # same peak and trough, smoother: w' = (exp(k cos 2d) - exp(-k)) / (exp(k) - exp(-k))
            k = 1.0
            d = np.arccos(np.clip(2 * w - 1, -1, 1)) / 2   # recover distance, radians
            w = (np.exp(k * np.cos(2 * d)) - np.exp(-k)) / (np.exp(k) - np.exp(-k))
        gain *= 1 + (g_rad - 1) * w
    return gain


def apply_gain(E, gain):
    return {k: e * gain[None] for k, e in E.items()}


# ---- normalization pools (parameterised) --------------------------------------------
def blur_xy(a, sigma_deg):
    """Spatial Gaussian on the last two axes, same implementation as the originals."""
    s = sigma_deg * PIXPDEG
    nd = a.ndim
    return utils_image.gaussian_filter_gpu(a, sigma=(0,) * (nd - 2) + (s, s), device="cpu")


def pool_untuned(e, q=1.0, sigma_xy_deg=3.0):
    """Notebook 'untuned': mean of e**q over all 24 SF x orientation channels,
    blurred in space. (The original's ~90-octave Gaussian over the 6- and
    4-element axes with reflect boundaries IS a flat mean; checked below.)"""
    z = (e ** q).mean(axis=(0, 1), keepdims=True)
    return np.broadcast_to(blur_xy(z, sigma_xy_deg), e.shape)


def pool_tuned(e, q=1.0, sigma_xy_deg=3.0, cross_weight=1.0):
    """Notebook 'tuned' (div_normalization tuned=True): each orientation
    channel's own e**q, blurred over SF (the same huge sigma -> SF mean) and
    space, plus a local (unblurred) cross-orientation term summing the other
    three channels' energy (to the power 1, as in the original)."""
    own = (e ** q).mean(axis=0, keepdims=True)
    own = np.broadcast_to(blur_xy(own, sigma_xy_deg), e.shape)
    cross = e.sum(axis=1, keepdims=True) - e
    return own + cross_weight * cross


def homogeneity_map(e, sigma_s_deg=10.0, sigma_c_px=1.0):
    """Cosine similarity, per pixel, of the 24-channel (SF x orientation) energy
    vector in a 1-pixel 'centre' with the same vector in a Gaussian 'surround'
    (normalization_byStimHomogeneity: sigma_c = 1 px, sigma_s = 10 deg)."""
    from scipy.ndimage import gaussian_filter
    c = gaussian_filter(e, sigma=(0, 0, sigma_c_px, sigma_c_px)).reshape(-1, *e.shape[2:])
    s = blur_xy(e, sigma_s_deg).reshape(c.shape)
    return (c * s).sum(0) / (np.sqrt((c ** 2).sum(0)) * np.sqrt((s ** 2).sum(0)) + 1e-9)


def pool_homogeneity(e, H, q_lo=1.0, q_hi=2.0, pool="local", sigma_xy_deg=3.0):
    """Pool e**q with a per-pixel exponent q = q_lo + (q_hi - q_lo) * H.
    pool='local' is the current repo version: a scipy Gaussian of sigma = 1
    index step on all four axes (SF, orientation, x, y) -- spatially local,
    broad across SF and orientation. pool='untuned' uses the broad spatial
    pool of models 1-2 instead (the earlier version of this function)."""
    q = q_lo + (q_hi - q_lo) * H
    eq = e ** q[None, None]
    if pool == "local":
        from scipy.ndimage import gaussian_filter
        return gaussian_filter(eq, sigma=(1, 1, 1, 1))
    z = eq.mean(axis=(0, 1), keepdims=True)
    return np.broadcast_to(blur_xy(z, sigma_xy_deg), e.shape)


def pool_anisotropy(e):
    """normalization_byAnisotropy_NOA: one scalar per stimulus, the std across
    orientation of energy averaged over SF and the whole image."""
    eori = e.mean(axis=(0, 2, 3))
    return np.full(e.shape, eori.std())


def normalize(E, pool_fn, sigma, p=1.0, **kw):
    return {k: e ** p / (sigma + pool_fn(e, **kw)) for k, e in E.items()}


# ---- checks -------------------------------------------------------------------------
def check_channel_orientations(m=None):
    m = m or make_model()
    found = []
    for freq_dir in [0, 45, 90, 135]:   # frequency direction, conventional
        g = m.generate_sinusoidal_grating(IMG_SIZE / 2 ** 5, np.deg2rad(freq_dir))
        e = m.get_energy_of_filters(m.forward(g), normalize=False).sum(axis=(2, 3))[3]
        found.append((freq_dir + 90) % 180)
        assert CH_BAR_ORI[np.argmax(e)] == (freq_dir + 90) % 180, (freq_dir, e)
    return found


def check_against_originals(E, tol=1e-4):
    """Re-implementations must match the repo's functions at the notebook's
    parameters (relative error on the 4-8 deg band)."""
    _, ecc = angle_maps()
    band = (ecc >= 4) & (ecc <= 8)
    sub = {k: E[k] for k in [(0, 0), (1, 0)]}
    report = {}
    pairs = {
        "untuned q=1": (utils_image.div_normalization(sub, PIXPDEG, tuned=False),
                        normalize(sub, pool_untuned, 0.01, q=1)),
        "untuned q=2": (utils_image.div_normalization(sub, PIXPDEG, q_exp=2, tuned=False),
                        normalize(sub, pool_untuned, 0.01, q=2)),
        "tuned q=1": (utils_image.div_normalization(sub, PIXPDEG, tuned=True),
                      normalize(sub, pool_tuned, 0.01, q=1)),
        # the repo's NOA averages over SF before dividing; compare SF means
        "anisotropy": (utils_image.normalization_byAnisotropy_NOA(sub, PIXPDEG),
                       {k: v.mean(axis=0, keepdims=True)
                        for k, v in normalize(sub, pool_anisotropy, 0.01).items()}),
        "homogeneity": (utils_image.normalization_byStimHomogeneity(sub, PIXPDEG),
                        {k: sub[k] / (0.1 + pool_homogeneity(sub[k], homogeneity_map(sub[k])))
                         for k in sub}),
    }
    for name, (orig, mine) in pairs.items():
        err = max(np.abs(orig[k][..., band] - mine[k][..., band]).max()
                  / np.abs(orig[k][..., band]).max() for k in sub)
        report[name] = err
    return report


# ---- extensions used only by the sweeps ---------------------------------------------
def pool_anisotropy_local(e, window_deg=None):
    """M4 with a spatial window: std across orientation of the SF-averaged
    energy, after Gaussian blurring each orientation channel by window_deg.
    window_deg=None is the repo's global version."""
    if window_deg is None:
        return pool_anisotropy(e)
    m = blur_xy(e.mean(axis=0), window_deg)          # (n_ori, H, W)
    return np.broadcast_to(m.std(axis=0)[None, None], e.shape)


def pool_homogeneity_wrap(e, H, q_lo=1.0, q_hi=2.0):
    """M5's local pool, but with the orientation axis treated as circular.
    The channel order [90, 45, 0, 135] steps by -45 deg, so wrapping index 3
    onto index 0 is the correct circular neighbour. The repo version reflects
    at the ends, which treats vertical and 135-deg channels differently from
    the 45-deg and horizontal channels."""
    from scipy.ndimage import gaussian_filter
    q = q_lo + (q_hi - q_lo) * H
    return gaussian_filter(e ** q[None, None], sigma=(1, 1, 1, 1),
                           mode=["reflect", "wrap", "reflect", "reflect"])


def band_index(ecc_lo=4, ecc_hi=8):
    _, ecc = angle_maps()
    return (ecc >= ecc_lo) & (ecc <= ecc_hi)


def readout_band(num, pool, sigmas, masks_band, out_gain=None):
    """Readout for several sigma at once, using only 4-8 deg pixels.
    num, pool: (n_sf, n_ori, n_band_pixels). masks_band: (8, n_band_pixels).
    out_gain: optional (n_ori, n_band_pixels) gain applied after normalization."""
    out = np.zeros((len(sigmas), masks_band.shape[0]))
    for s_i, s in enumerate(sigmas):
        r = num / (s + pool)
        if out_gain is not None:
            r = r * out_gain[None]
        img = r.mean(axis=0).sum(axis=0)
        out[s_i] = [img[m].mean() for m in masks_band]
    return out
