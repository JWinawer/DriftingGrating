"""Steerable-pyramid front end: energy of the model stimuli, 8 orientation channels, 50% contrast.

Stimuli and pyramid are the repo's own (helper_functions/pyrmodel.py, unmodified), with three
parameter changes from the notebook, each to match the experiment:

  * 8 orientation channels (pyramid order 7) instead of 4. The energy tuning then has a half-width
    at half-height of about 18 deg (cos^14), against 27 deg with 4 channels, and a tuned pool has
    neighbouring orientations to weight.
  * 50% Michelson contrast (ExperimentCode/Config/constConfig.m, const.contrast = .5).
  * Polar gratings at 1 cpd at 6 deg eccentricity. The notebook's InputStimuli matches 1 cpd at
    a quarter of the padded image (7.2 deg). The experiment matches it at half the stimulus radius
    (ExperimentCode/Stim/my_stim.m), 6.1 deg, which gives round(2 pi 6.1) = 38 cycles around the
    circle; passing sf = 6 / 7.2 to InputStimuli gives 38 as well.

Energies are computed at full resolution (896 px, 31.2 px/deg), then averaged in DS x DS blocks
(default 4: 224 px, 7.8 px/deg) for the factorial. Energy here is the squared magnitude of a
complex (quadrature) filter, so it is smooth on the scale of a block; check_backend.py measures
what the block averaging does to the asymmetries.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from helper_functions import pyrmodel, utils_image  # noqa: E402

from backend import Geometry, WEDGES, circdiff  # noqa: E402

# geometry, as in model_DriftingGratings.ipynb cell 5 and Models/paradoxes/paradox_lib.py
IMG_SIZE = 896
EXP_STIM_DIAM_DEG = 24.4
_x = int(np.round(IMG_SIZE * 0.85)); _x = _x if _x % 2 == 0 else _x + 1
PIXPDEG = _x / EXP_STIM_DIAM_DEG
PADDED_DEG = IMG_SIZE / PIXPDEG
N_SF = 6
CONTRAST = 0.5
POLAR_MATCH_DEG = 6.0
CACHE = os.path.join(HERE, '_cache')


class PyramidGeometry(Geometry):
    """Image grid, optionally block-averaged by ds. theta: bar orientation of each channel."""

    def __init__(self, theta, ds=1, ecc_lo=4, ecc_hi=8):
        self.theta = np.asarray(theta, float)
        self.ds = ds
        n = IMG_SIZE // ds
        self.ppd = PIXPDEG / ds
        c = (IMG_SIZE - 1) / 2
        pix = ds * np.arange(n) + (ds - 1) / 2                     # block centres, full-res pixels
        x = (pix[None, :] - c) / PIXPDEG
        y = -(pix[:, None] - c) / PIXPDEG
        self.alpha = (np.degrees(np.arctan2(y, x)) + 360) % 360
        self.ecc = np.sqrt(x ** 2 + y ** 2)
        band = (self.ecc >= ecc_lo) & (self.ecc <= ecc_hi)
        self.masks = np.stack([band & (np.abs(circdiff(self.alpha, w, 360)) < 22.5) for w in WEDGES])

    def blur(self, X, s_deg):
        if s_deg == 0:
            return X
        if np.isinf(s_deg):
            return np.broadcast_to(self.spatial_mean(X), X.shape)
        s = s_deg * self.ppd
        return utils_image.gaussian_filter_gpu(np.ascontiguousarray(X, dtype=np.float64),
                                               sigma=(0,) * (X.ndim - 2) + (s, s), device='cpu')


def make_pyramid(n_ori):
    return pyrmodel.SteerablePyramidSF(n_ori=n_ori, n_sf=N_SF, img_size=IMG_SIZE,
                                       img_degrees=PADDED_DEG, is_complex=True,
                                       downsample=False, tight_frame=True, device='cpu')


def channel_orientations(m, n_ori):
    """Bar orientation of each orientation channel, measured with gratings at ~1 cpd."""
    out = np.zeros(n_ori)
    for freq_dir in np.arange(n_ori) * 180 / n_ori:
        g = m.generate_sinusoidal_grating(IMG_SIZE / 2 ** 5, np.deg2rad(freq_dir))
        e = m.get_energy_of_filters(m.forward(g), normalize=False).sum(axis=(2, 3))[3]
        out[np.argmax(e)] = (freq_dir + 90) % 180
    assert len(set(out)) == n_ori, out
    return out


def make_stimuli(contrast=CONTRAST, polar_match_deg=POLAR_MATCH_DEG):
    """(name, image) for the 8 stimuli: Cartesian then polar, in the notebook's order."""
    import contextlib
    import io
    r_quarter_deg = (IMG_SIZE / 4) / PIXPDEG                       # where InputStimuli matches
    out = []
    for ref, sf in (('cartesian', 1.0), ('polar', polar_match_deg / r_quarter_deg)):
        with contextlib.redirect_stdout(io.StringIO()):
            s = pyrmodel.InputStimuli(sfs=np.array([sf]), oris=np.array([0, 135, 90, 45]),
                                      img_size=IMG_SIZE, img_degrees=PADDED_DEG,
                                      referenceframe=ref, mask=True,
                                      paddedBondaryPerc=EXP_STIM_DIAM_DEG / PADDED_DEG)
        out += [(f'{ref}{j}', contrast * s.input[(0, j)]) for j in range(4)]
    return out


def block_mean(a, ds):
    if ds == 1:
        return a
    n = a.shape[-1] // ds
    return a.reshape(*a.shape[:-2], n, ds, n, ds).mean(axis=(-3, -1))


def dominant_bar(E, geom, sf_index=3):
    """Energy-weighted mean bar orientation (doubled-angle) in each wedge, rounded to 45 deg."""
    out = np.zeros(len(geom.masks))
    for k, m in enumerate(geom.masks):
        w = E[sf_index][:, m].mean(1)
        z = (w * np.exp(2j * np.deg2rad(geom.theta))).sum()
        out[k] = (np.round(np.degrees(np.angle(z)) / 2 / 45) * 45) % 180
    return out


def energies(n_ori=8, ds=4, contrast=CONTRAST, log=print):
    """Energies (F, K, H/ds, W/ds) and labels for the 8 stimuli, cached in _cache/."""
    path = os.path.join(CACHE, f'pyr_K{n_ori}_ds{ds}_c{contrast:g}.npz')
    if os.path.exists(path):
        z = np.load(path, allow_pickle=True)
        return list(z['E']), list(z['stims']), z['theta']
    t0 = time.time()
    m = make_pyramid(n_ori)
    theta = channel_orientations(m, n_ori)
    log(f'channel bar orientations: {theta}')
    geom = PyramidGeometry(theta, ds=ds)
    Es, stims = [], []
    imgs = make_stimuli(contrast)
    for i, (name, img) in enumerate(imgs):
        e = m.get_energy_of_filters(m.forward(img), to_magnitude=False, to_numpy=True,
                                    reshape=True, normalize=False)
        Es.append(block_mean(e, ds).astype(np.float64))
        bar = dominant_bar(Es[-1], geom)
        stims.append(dict(exp='dg' if name.startswith('cart') else 'da', name=name, bar=bar))
        log(f'  {i + 1}/{len(imgs)} {name}: bar at wedges {bar.astype(int)}  '
            f'({time.time() - t0:.0f}s elapsed)')
    os.makedirs(CACHE, exist_ok=True)
    np.savez(path, E=np.stack(Es), stims=np.array(stims, dtype=object), theta=theta)
    return Es, stims, theta
