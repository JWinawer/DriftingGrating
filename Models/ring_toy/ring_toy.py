"""Ring-by-orientation toy model of the cardinal / radial BOLD paradoxes.

No images. V1 is a ring of polar-angle locations (the stimulus annulus) x a bank of
orientation channels. Each location has a density of units per preferred orientation, with
optional horizontal, cardinal and radial over-representation. A stimulus is one orientation
per location: Cartesian gratings (dg) have the same orientation everywhere, polar gratings (da)
have an orientation fixed relative to the radial direction. Units are divisively normalized by a
pool that can be tuned or untuned in orientation, local or global in space, gated by stimulus
homogeneity, and raised to a power q. BOLD at a location is the density-weighted sum of unit
responses; wedge responses are averaged over the 45-deg wedge, as in the data.

Angles are in degrees. Polar angle 0 = right horizontal meridian, counter-clockwise. Orientation
0 = horizontal. The radial orientation at polar angle a is a (mod 180).
"""
import numpy as np

L, K = 72, 36
ALPHA = np.arange(L) * 360 / L          # ring locations (polar angle)
THETA = np.arange(K) * 180 / K          # preferred orientations
WEDGES = np.arange(8) * 45              # wedge centres, as in meanBOLDpa

DEFAULTS = dict(
    a_hv=0.10,     # density: horizontal > vertical (cos 2theta)
    a_card=0.15,   # density: cardinal > oblique (cos 4theta)
    a_rad=0.15,    # density: radial > tangential (cos 2(theta - alpha))
    b_card=0.0,    # tuning: kappa multiplied by (1 + b_card cos 4theta), i.e. narrower at cardinals
    kappa=2.0,     # tuning concentration (von Mises on doubled angle); 2 ~ 35 deg half-width
    n=2.0,         # numerator exponent (energy)
    c=1.0,         # contrast
    norm=True,
    sigma=0.1,     # semi-saturation; smaller = stronger normalization
    q=1.0,         # exponent on the pooled drive; q = 1 is the standard (Heeger) form
    kappa_pool=0.0,  # orientation tuning of the pool; 0 = untuned
    s_sp=np.inf,   # spatial s.d. of the pool along the ring (deg of polar angle); inf = global
    gate=0.0,      # homogeneity gate exponent; 0 = off
    a_pool=0.0,    # extra cardinal weighting of the suppressive pool only
    rad_site='density',  # 'density': radial units feed the pool too; 'output': radial gain after normalization
)


def cosd(x):
    return np.cos(np.deg2rad(x))


def circdiff(a, b, period):
    return (a - b + period / 2) % period - period / 2


def density(p, with_rad=True):
    th, al = THETA[None, :], ALPHA[:, None]
    d = 1 + p['a_hv'] * cosd(2 * th) + p['a_card'] * cosd(4 * th)
    d = d + (p['a_rad'] * cosd(2 * (th - al)) if with_rad else 0 * al)
    assert d.min() > 0, 'density must stay positive'
    return d


def tuning(phi, p):
    kap = p['kappa'] * (1 + p['b_card'] * cosd(4 * THETA))[None, :]
    return np.exp(kap * (cosd(2 * (phi[:, None] - THETA[None, :])) - 1))


def von_mises_rows(delta, kappa, period):
    """Row-normalized weights over a circular difference; kappa = 0 gives uniform weights."""
    w = np.exp(kappa * (np.cos(2 * np.pi * delta / period) - 1))
    return w / w.sum(1, keepdims=True)


def spatial_weights(s_sp):
    dl = circdiff(ALPHA[:, None], ALPHA[None, :], 360)
    w = np.ones((L, L)) if np.isinf(s_sp) else np.exp(-0.5 * (dl / s_sp) ** 2)
    return w / w.sum(1, keepdims=True)


WO_CACHE, WS_CACHE = {}, {}


def respond(phi, p):
    """BOLD at each ring location for stimulus orientations phi (length L)."""
    d = density(p)
    En = (p['c'] * tuning(phi, p)) ** p['n']
    if not p['norm']:
        R = En
    else:
        if p['kappa_pool'] not in WO_CACHE:
            WO_CACHE[p['kappa_pool']] = von_mises_rows(
                circdiff(THETA[:, None], THETA[None, :], 180), p['kappa_pool'], 180)
        if p['s_sp'] not in WS_CACHE:
            WS_CACHE[p['s_sp']] = spatial_weights(p['s_sp'])
        Wo, Ws = WO_CACHE[p['kappa_pool']], WS_CACHE[p['s_sp']]
        if p['gate'] > 0:
            # surround counts only as far as its stimulus orientation matches the centre's
            h = ((1 + cosd(2 * (phi[:, None] - phi[None, :]))) / 2) ** p['gate']
            Ws = Ws * h
        # with rad_site='output' the radial bias multiplies responses but does not feed the pool
        d_pool = density(p, with_rad=p['rad_site'] == 'density')
        X = d_pool * En * (1 + p['a_pool'] * cosd(4 * THETA))[None, :]
        pool = Ws @ X @ Wo.T
        R = En / (p['sigma'] ** p['n'] + pool ** p['q'])
    return (d * R).mean(1)


def stimuli():
    """(experiment, stimulus orientation parameter) -> orientation at each ring location."""
    out = {}
    for s in (0, 45, 90, 135):
        out[('dg', s)] = np.full(L, float(s))                 # Cartesian: same everywhere
        out[('da', s)] = (ALPHA + s) % 180                    # polar: s relative to radial
    return out


def wedge_table(p):
    """One row per (experiment, stimulus, wedge): response and its labels at the wedge centre."""
    rows = []
    for (exp, s), phi in stimuli().items():
        b = respond(phi, p)
        for w in WEDGES:
            inw = np.abs(circdiff(ALPHA, w, 360)) <= 22.5
            loc = phi[np.argmin(np.abs(circdiff(ALPHA, w, 360)))] % 180   # orientation at centre
            rel = (loc - w) % 180
            rows.append(dict(exp=exp, stim=s, wedge=w, resp=b[inw].mean(),
                             cart='card' if loc % 90 == 0 else 'obl',
                             hv={0: 'h', 90: 'v'}.get(loc, ''),
                             polar={0: 'rad', 90: 'tang'}.get(rel, 'spiral')))
    return rows


def contrasts(p):
    """The paper's four asymmetries per experiment, as (A - B) / mean response of that experiment."""
    rows = wedge_table(p)
    out = {}
    for exp in ('dg', 'da'):
        r = [x for x in rows if x['exp'] == exp]
        m = np.mean([x['resp'] for x in r])

        def mean_of(key, val):
            return np.mean([x['resp'] for x in r if x[key] == val])

        out[exp] = {
            'horiz-vert': (mean_of('hv', 'h') - mean_of('hv', 'v')) / m,
            'card-obl': (mean_of('cart', 'card') - mean_of('cart', 'obl')) / m,
            'rad-tang': (mean_of('polar', 'rad') - mean_of('polar', 'tang')) / m,
            'polc-polo': ((mean_of('polar', 'rad') + mean_of('polar', 'tang')) / 2
                          - mean_of('polar', 'spiral')) / m,
        }
    return out


# Signs and context effects in V1, 4-8 deg (Reproduction/RESULTS.md, section 2).
DATA = {'dg': {'horiz-vert': -0.514, 'card-obl': -0.193, 'rad-tang': 0.091, 'polc-polo': 0.056},
        'da': {'horiz-vert': -0.212, 'card-obl': -0.044, 'rad-tang': 0.198, 'polc-polo': 0.050}}


def targets(c):
    """The qualitative patterns the model should reproduce, each True/False."""
    return {
        'P1 card<obl (dg)': c['dg']['card-obl'] < 0,
        'P1 horiz<vert (dg)': c['dg']['horiz-vert'] < 0,
        'P1 horiz<vert (da)': c['da']['horiz-vert'] < 0,
        'P2 rad>tang (both)': c['dg']['rad-tang'] > 0 and c['da']['rad-tang'] > 0,
        'P3 card-obl more neg in dg': c['dg']['card-obl'] < c['da']['card-obl'],
        'P3 rad-tang larger in da': c['da']['rad-tang'] > c['dg']['rad-tang'],
    }
