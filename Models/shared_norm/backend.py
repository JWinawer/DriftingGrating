"""Shared divisive-normalization back end for the orientation-paradox models.

One normalization formula, with switches, so that the ring toy (Models/ring_toy/) and the
steerable-pyramid model (Models/paradoxes/) can be run with every component defined the same
way. A front end supplies only three things:

  E      energy per stimulus, array (F, K, *S): F spatial-frequency channels, K orientation
         channels, S spatial locations (a ring (L,) or an image (H, W)).
  geom   a Geometry: channel orientations theta (K,), polar angle of every location (*S),
         a spatial blur in visual degrees, and the eight 45-deg wedge masks at 4-8 deg.
  stims  one label per stimulus: experiment ('dg' Cartesian, 'da' polar) and the bar
         orientation at each wedge centre.

The response of channel (f, k) at location x is

    R = W * (g E)^p / (sigma + P)                  W = readout weight (density, output gain)
    P = Pool[ D * (g E)^q ]        (q 'inside')     g = input gain (applied before exponents)
    P = Pool[ D * (g E)^p ]^q      (q 'outside')    D = unit density (feeds the pool)

and BOLD at a location is the mean of R over SF and orientation channels. Three
over-representations (horizontal > vertical, cardinal > oblique, radial > tangential) each enter
at one of three sites:

  'density'  more units with unchanged gain: they raise the readout and feed the pool
             (the ring toy's density; the 'count' mode of Models/paradoxes/02_sweeps.py)
  'gain'     a gain on energy before any exponent (the repo's imposeImbalance)
  'output'   a gain on the normalized response; it never reaches the pool
             (the ring toy's rad_site='output'; "radial gain on the output" in PARADOXES.md)

With p = q = 1 'density' and 'gain' are the same model; they differ once q != p.

Pool switches: orientation weighting (von Mises in kappa on the doubled angle, rows sum to 1;
0 = untuned), spatial s.d. in visual degrees (0 = local, inf = global mean), SF weighting (flat
mean by default), pool type ('sum': weighted mean of the drive, the repo's M1-M3 and the toy;
'spread': weighted s.d. of the drive across orientation, the repo's M4), and homogeneity gating
(the repo's M5: the pool exponent runs from q_lo to q_hi with the cosine similarity of the local
channel vector to its spatial surround). Two legacy terms exist only to reproduce the repo's M3
and M5 exactly (cross_local; explicit Wo, Wf); the factorial does not use them.

Angles in degrees; orientation 0 = horizontal, counter-clockwise; polar angle 0 = right
horizontal meridian, counter-clockwise. Radial orientation at polar angle a is a (mod 180).
"""
import numpy as np

COMPONENTS = ('hv', 'card', 'rad')
SITES = ('density', 'gain', 'output')
WEDGES = np.arange(8) * 45

DEFAULTS = dict(
    # over-representation: factor 1 + a * u per channel and location
    a_hv=0.10, a_card=0.15, a_rad=0.15,          # the ring toy's amplitudes
    site_hv='density', site_card='density', site_rad='density',
    form='sym',            # 'sym': u in [-1, 1] (toy); 'up': u in [0, 1] (repo: obliques at 1)
    rad_profile='cos',     # 'cos': cos 2(theta - alpha) (toy); 'linear': repo's ramp in angle
    combine='mul',         # components at a site multiply; 'add' sums deviations (toy)
    # exponents
    p=1.0, q=1.0, q_where='inside',
    gate=False, q_lo=1.0, q_hi=2.0, s_surround_deg=10.0, s_center_deg=0.0,
    # pool
    pool='sum',            # 'sum' | 'spread' | 'none'
    kappa=0.0,             # orientation tuning of the pool; 0 untuned, np.inf own channel only
    s_deg=3.0,             # spatial s.d. of the pool, visual degrees; 0 local, np.inf global
    sigma=0.01, sigma_mode='abs',   # 'rel': sigma is a multiple of the mean pool (see run)
    # legacy, for exact reproduction of the repo's M3 and M5 only
    cross_local=0.0, Wo=None, Wf=None,
)


def cosd(x):
    return np.cos(np.deg2rad(x))


def circdiff(a, b, period):
    return (a - b + period / 2) % period - period / 2


class Geometry:
    """Base class. Front ends set theta (K,), alpha (*S), masks (8, *S) and implement blur."""
    theta = alpha = masks = None

    def blur(self, X, s_deg):
        """Gaussian of s.d. s_deg (visual degrees) over the trailing spatial axes of X.
        s_deg = 0 returns X; s_deg = inf returns the mean over all locations."""
        raise NotImplementedError

    def spatial_mean(self, X):
        nd = self.alpha.ndim
        return X.mean(axis=tuple(range(-nd, 0)), keepdims=True)


# ---- over-representation --------------------------------------------------------------
def u_component(name, geom, prm):
    """Signed profile in [-1, 1], shape (K, *S) or broadcastable to it."""
    th = geom.theta.reshape((-1,) + (1,) * geom.alpha.ndim)
    if name == 'hv':
        u = cosd(2 * th)
    elif name == 'card':
        u = cosd(4 * th)
    else:
        d = circdiff(th, geom.alpha[None], 180)                 # 0 radial, +-90 tangential
        u = cosd(2 * d) if prm['rad_profile'] == 'cos' else 1 - 2 * np.abs(d) / 90
    return u if prm['form'] == 'sym' else (1 + u) / 2


def site_factors(geom, prm):
    """(gain, pool density, readout weight), each broadcastable to (K, *S)."""
    f = {c: 1 + prm['a_' + c] * u_component(c, geom, prm) for c in COMPONENTS}

    def comb(names):
        if prm['combine'] == 'add':
            return 1 + sum((f[c] - 1 for c in names), 0.0)
        out = 1.0
        for c in names:
            out = out * f[c]
        return out

    at = {s: [c for c in COMPONENTS if prm['site_' + c] == s] for s in SITES}
    return comb(at['gain']), comb(at['density']), comb(at['density'] + at['output'])


# ---- pool weights -----------------------------------------------------------------------
def von_mises_weights(theta, kappa):
    if np.isinf(kappa):
        return np.eye(len(theta))
    w = np.exp(kappa * (cosd(2 * circdiff(theta[:, None], theta[None, :], 180)) - 1))
    return w / w.sum(1, keepdims=True)


def index_gauss_weights(n, sigma):
    """scipy gaussian_filter1d (mode 'reflect') along an index axis, as a matrix: the repo's M5
    pools SF and orientation channels this way."""
    from scipy.ndimage import gaussian_filter1d
    return gaussian_filter1d(np.eye(n), sigma, axis=0, mode='reflect')


def mix(W, X, axis):
    """Apply an (n, n) weight matrix along one axis of X."""
    return np.moveaxis(np.tensordot(W, np.moveaxis(X, axis, 0), axes=(1, 0)), 0, axis)


def homogeneity(V, geom, prm):
    """Cosine similarity of each location's (F x K) channel vector with its surround's."""
    c = geom.blur(V, prm['s_center_deg']).reshape(-1, *V.shape[2:])
    s = geom.blur(V, prm['s_surround_deg']).reshape(c.shape)
    return (c * s).sum(0) / (np.sqrt((c ** 2).sum(0)) * np.sqrt((s ** 2).sum(0)) + 1e-9)


def pool_drive(E, geom, prm, fac=None, H=None):
    """Numerator, pool and readout weight for one stimulus: R = W * num / (sigma + P).
    H: precomputed homogeneity map, optional (it depends only on which biases reach the pool)."""
    g, D, W = site_factors(geom, prm) if fac is None else fac
    gE = g * E
    num = gE ** prm['p']
    if prm['pool'] == 'none':
        return num, None, W
    if prm['gate']:
        if H is None:
            H = homogeneity(D * gE, geom, prm)
        q = prm['q_lo'] + (prm['q_hi'] - prm['q_lo']) * H
    else:
        q = prm['q']
    X = D * gE ** q if prm['q_where'] == 'inside' else D * num
    Wo = von_mises_weights(geom.theta, prm['kappa']) if prm['Wo'] is None else prm['Wo']
    if prm['pool'] == 'sum':
        Y = X.mean(0, keepdims=True) if prm['Wf'] is None else mix(prm['Wf'], X, 0)
        P = geom.blur(mix(Wo, Y, 1), prm['s_deg'])
    elif prm['pool'] == 'spread':
        Y = geom.blur(X.mean(0, keepdims=True), prm['s_deg'])       # (1, K, *S)
        mu = mix(Wo, Y, 1)
        # weighted variance around each channel's own weighted mean
        Y2 = mix(Wo, Y ** 2, 1)
        P = np.sqrt(np.maximum(Y2 - mu ** 2, 0))
    else:
        raise ValueError(prm['pool'])
    if prm['cross_local']:
        P = P + prm['cross_local'] * D * (gE.sum(1, keepdims=True) - gE)
    if prm['q_where'] == 'outside':
        P = P ** q
    return num, P, W


# ---- run and read out -------------------------------------------------------------------
def readout_map(R):
    """BOLD map: mean over SF and orientation channels."""
    return R.mean(axis=(0, 1))


def wedge_means(img, geom):
    return np.array([img[m].mean() for m in geom.masks])


def run(Es, geom, prm, sigmas=None, hcache=None):
    """Wedge responses (n_sigma, n_stim, 8) for every stimulus energy in Es.

    sigma_mode 'abs': sigma as given. 'rel': sigma = value x the mean pool P over the 4-8 deg
    band and all stimuli, so that one number means equally strong normalization on every front
    end, exponent and pool type. Also returns the absolute sigmas used."""
    prm = dict(DEFAULTS, **prm)
    sigmas = [prm['sigma']] if sigmas is None else list(sigmas)
    fac = site_factors(geom, prm)
    band = geom.masks.any(0)
    Hs = [None] * len(Es)
    if prm['gate'] and hcache is not None:
        key = (prm['combine'], prm['form'], prm['rad_profile'], prm['s_center_deg'], prm['s_surround_deg'],
               *[(prm['a_' + c], prm['site_' + c]) for c in COMPONENTS])
        if key not in hcache:
            hcache[key] = [homogeneity(fac[1] * fac[0] * E, geom, prm) for E in Es]
        Hs = hcache[key]
    parts = [pool_drive(E, geom, prm, fac, H) for E, H in zip(Es, Hs)]
    if prm['pool'] == 'none':
        resp = np.array([wedge_means(readout_map(W * num), geom) for num, _, W in parts])
        return np.repeat(resp[None], len(sigmas), 0), np.full(len(sigmas), np.nan)
    if prm['sigma_mode'] == 'rel':
        ref = np.mean([np.broadcast_to(P, num.shape)[..., band].mean() for num, P, _ in parts])
        abs_s = [s * ref for s in sigmas]
    else:
        abs_s = sigmas
    # only the 4-8 deg band is read out
    mb = geom.masks[:, band]
    out = np.zeros((len(sigmas), len(Es), len(geom.masks)))
    for j, (num, P, W) in enumerate(parts):
        nb = num[..., band]
        Pb = np.broadcast_to(P, num.shape)[..., band]
        Wb = np.broadcast_to(W, num.shape[1:])[..., band]
        for i, s in enumerate(abs_s):
            img = (Wb * nb / (s + Pb)).mean(axis=(0, 1))
            out[i, j] = [img[m].mean() for m in mb]
    return out, np.array(abs_s)


# ---- asymmetries and scoring ------------------------------------------------------------
def labels(bar):
    """Classes of a bar orientation at each wedge centre."""
    rel = (bar - WEDGES) % 180
    return dict(H=bar % 180 == 0, V=bar % 180 == 90, card=bar % 90 == 0, obl=bar % 90 == 45,
                rad=rel == 0, tang=rel == 90, polc=rel % 90 == 0, polo=rel % 90 == 45)


def asymmetries(resp, stims):
    """The four asymmetries per experiment, as fractions of that experiment's mean response.

    resp (n_stim, 8). stims: list of dicts with 'exp' and 'bar' (bar orientation at each of the
    8 wedge centres). Within each wedge an asymmetry is the mean of one class minus the mean of
    the other; wedges are then weighted equally (the equal-coverage rule of Reproduction/).
    Raw (un-normalized) values and the mean are returned too."""
    out = {}
    for exp in ('dg', 'da'):
        idx = [j for j, s in enumerate(stims) if s['exp'] == exp]
        lab = [labels(np.asarray(stims[j]['bar'])) for j in idx]
        r = resp[idx]
        acc = {}
        for name, (a, b) in dict(hv=('H', 'V'), co=('card', 'obl'), rt=('rad', 'tang'),
                                 pp=('polc', 'polo')).items():
            d = []
            for k in range(r.shape[1]):
                ia = [i for i in range(len(idx)) if lab[i][a][k]]
                ib = [i for i in range(len(idx)) if lab[i][b][k]]
                if ia and ib:
                    d.append(r[ia, k].mean() - r[ib, k].mean())
            acc[name] = float(np.mean(d))
        m = float(r.mean())
        out[exp] = {**{k: v / m for k, v in acc.items()},
                    **{'raw_' + k: v for k, v in acc.items()}, 'mean': m}
    return out


# V1, 4-8 deg, matched 7 observers (Reproduction/RESULTS.md section 2; ring_toy.DATA)
DATA = {'dg': dict(hv=-0.514, co=-0.193, rt=0.091, pp=0.056),
        'da': dict(hv=-0.212, co=-0.044, rt=0.198, pp=0.050)}

TARGETS = {
    # the cardinal / radial paradoxes; the last two are figures.py's beyond-geometry targets
    'P1 dg card<obl': lambda c: c['dg']['co'] < 0,
    'P2 rad>tang both': lambda c: c['dg']['rt'] > 0 and c['da']['rt'] > 0,
    'P3 card context dg<da': lambda c: c['dg']['co'] < c['da']['co'],
    'P3 rad context da>dg': lambda c: c['da']['rt'] > c['dg']['rt'],
    'P3 strict card: da > dg/2': lambda c: c['da']['co'] > 0.5 * c['dg']['co'],
    'P3 strict rad: da > 1.5 dg': lambda c: c['da']['rt'] > 1.5 * c['dg']['rt'],
    # horizontal - vertical, scored separately
    'HV dg h<v': lambda c: c['dg']['hv'] < 0,
    'HV da h<v': lambda c: c['da']['hv'] < 0,
    'HV context dg<da': lambda c: c['dg']['hv'] < c['da']['hv'],
}
CORE = [k for k in TARGETS if not k.startswith('HV')]
HV = [k for k in TARGETS if k.startswith('HV')]


def score(c):
    t = {k: bool(f(c)) for k, f in TARGETS.items()}
    t['core all'] = all(t[k] for k in CORE)
    t['HV all'] = all(t[k] for k in HV)
    t['all'] = t['core all'] and t['HV all']
    # sign agreement with the data: 8 asymmetries, and the 3 context effects that are
    # detectable in the data (polc-polo's context effect is +0.006, p = .81, so not scored)
    t['n signs (of 8)'] = sum(np.sign(c[e][a]) == np.sign(DATA[e][a]) for e in DATA for a in DATA[e])
    t['n context signs (of 3)'] = sum(
        np.sign(c['dg'][a] - c['da'][a]) == np.sign(DATA['dg'][a] - DATA['da'][a])
        for a in ('hv', 'co', 'rt'))
    return t
