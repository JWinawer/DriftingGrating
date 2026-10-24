"""Step 1: the back end reproduces the repo's five pools, the repo's imbalance, and the ring toy.

  A. Rania's M1-M5 (helper_functions/utils_image.py, called unmodified) on the notebook's own
     4-channel, 100%-contrast energies (Models/paradoxes/_cache/energy.npz, built by
     paradoxes/00_checks.py), against the back end with the switch settings in REPO_MODELS.
     Max relative error over the 4-8 deg band, every channel. paradox_lib.py matches the same
     functions to ~1e-5; the originals blur in float32.
  B. The repo's imposeImbalance (cardinal 2, radial 1.5) against the back end's 'gain' site.
  C. ring_toy.respond against the back end on the ring front end, BOLD at every ring location,
     for the toy's default, its best model and seven variations.
  D. p = q = 1 makes 'density' and 'gain' the same model (a consistency check).
  E. What 4 x 4 block averaging of the pyramid energy (used by the factorial) does to the
     asymmetries, on the 4-channel energies, for six pool settings.

Run from Models/shared_norm/. About 3 minutes; writes _cache/check_backend.log.
"""
import copy
import os
import sys
import time

import numpy as np

import backend as B
import frontend_pyramid as FP
import frontend_ring as FR

sys.path.insert(0, os.path.join(os.path.dirname(FP.HERE), 'ring_toy'))
import ring_toy as T  # noqa: E402
from helper_functions import utils_image as U  # noqa: E402

LOG = os.path.join(FP.CACHE, 'check_backend.log')
os.makedirs(FP.CACHE, exist_ok=True)
open(LOG, 'w').close()


def log(msg):
    line = time.strftime('%H:%M:%S ') + msg
    print(line, flush=True)
    with open(LOG, 'a') as f:
        f.write(line + '\n')


CH4 = np.array([90, 45, 0, 135])     # bar orientation of the notebook's 4 channels (measured)
OFF = dict(a_hv=0.0, a_card=0.0, a_rad=0.0)
REPO_MODELS = {
    'M1 untuned, q=1': dict(kappa=0.0, s_deg=3.0, q=1.0, sigma=0.01),
    'M2 untuned, q=2': dict(kappa=0.0, s_deg=3.0, q=2.0, sigma=0.01),
    'M3 tuned': dict(kappa=np.inf, s_deg=3.0, q=1.0, sigma=0.01, cross_local=1.0),
    'M4 anisotropy': dict(pool='spread', kappa=0.0, s_deg=np.inf, sigma=0.01),
    'M5 homogeneity': dict(gate=True, q_lo=1.0, q_hi=2.0, s_center_deg=1 / FP.PIXPDEG,
                           s_surround_deg=10.0, Wf=B.index_gauss_weights(6, 1.0),
                           Wo=B.index_gauss_weights(4, 1.0), s_deg=1 / FP.PIXPDEG, sigma=0.1),
}


def originals(sub, name):
    f = {'M1 untuned, q=1': lambda: U.div_normalization(sub, FP.PIXPDEG, tuned=False),
         'M2 untuned, q=2': lambda: U.div_normalization(sub, FP.PIXPDEG, q_exp=2, tuned=False),
         'M3 tuned': lambda: U.div_normalization(sub, FP.PIXPDEG, tuned=True),
         'M4 anisotropy': lambda: U.normalization_byAnisotropy_NOA(sub, FP.PIXPDEG),
         'M5 homogeneity': lambda: U.normalization_byStimHomogeneity(sub, FP.PIXPDEG, device='cpu')}
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        return f[name]()


def response(E, geom, prm):
    prm = dict(B.DEFAULTS, **prm)
    num, P, W = B.pool_drive(E, geom, prm)
    return W * num if P is None else W * num / (prm['sigma'] + P)


def part_a(E4):
    geom = FP.PyramidGeometry(CH4, ds=1)
    band = geom.ecc
    band = (band >= 4) & (band <= 8)
    sub = {k: E4[k] for k in [(0, 0), (1, 0)]}
    log('A. repo pools vs back end (max relative error on the 4-8 deg band):')
    worst = 0
    for name, prm in REPO_MODELS.items():
        t0 = time.time()
        orig = originals(sub, name)
        err = 0
        for k, e in sub.items():
            mine = response(e, geom, dict(OFF, **prm))
            if name == 'M4 anisotropy':          # the repo's M4 averages over SF before dividing
                mine = mine.mean(0, keepdims=True)
            err = max(err, np.abs(orig[k][..., band] - mine[..., band]).max()
                      / np.abs(orig[k][..., band]).max())
        worst = max(worst, err)
        log(f'   {name:18s} {err:.1e}   ({time.time() - t0:.0f}s)')
    return worst


def part_b(E4):
    geom = FP.PyramidGeometry(CH4, ds=1)
    sub = {(0, 0): E4[1, 0]}       # one key: imposeImbalance loops over range(n_set)
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        orig = U.imposeImbalance(copy.deepcopy(sub), cardinal=2, radial=1.5)[0, 0]
    prm = dict(B.DEFAULTS, form='up', a_hv=0.0, a_card=1.0, a_rad=0.5, rad_profile='linear',
               site_card='gain', site_rad='gain')
    g, _, _ = B.site_factors(geom, prm)
    mine = g * E4[1, 0]
    err = np.abs(orig - mine).max() / np.abs(orig).max()
    log(f'B. repo imposeImbalance(cardinal 2, radial 1.5) vs back end gain site: {err:.1e}')
    return err


def toy_to_backend(p, R_ecc=FR.ECC):
    """ring_toy parameters -> back-end switches. Requires gate = a_pool = b_card = 0."""
    assert p['gate'] == 0 and p['a_pool'] == 0 and p['b_card'] == 0
    return dict(a_hv=p['a_hv'], a_card=p['a_card'], a_rad=p['a_rad'], form='sym',
                rad_profile='cos', combine='add', site_hv='density', site_card='density',
                site_rad='density' if p['rad_site'] == 'density' else 'output',
                p=p['n'] / 2, q=p['q'], q_where='outside', kappa=p['kappa_pool'],
                s_deg=np.inf if np.isinf(p['s_sp']) else np.deg2rad(p['s_sp']) * R_ecc,
                sigma=p['sigma'] ** p['n'], pool='sum' if p['norm'] else 'none')


def part_c():
    geom = FR.RingGeometry(T.L, T.K)
    best = dict(T.DEFAULTS, rad_site='output', kappa_pool=4.0, s_sp=15.0, q=2.0, sigma=0.03)
    cases = {'toy default': T.DEFAULTS, 'toy best': best,
             'no normalization': dict(best, norm=False), 'untuned global': dict(best, kappa_pool=0.0, s_sp=np.inf),
             'local 5 deg': dict(best, s_sp=5.0), 'radial in pool': dict(best, rad_site='density'),
             'q = 1': dict(best, q=1.0), 'sigma 1': dict(best, sigma=1.0),
             'contrast 0.5, n 1': dict(best, c=0.5, n=1.0)}
    log('C. ring_toy.respond vs back end on the ring front end (max relative error, BOLD per location):')
    worst = 0
    for name, p in cases.items():
        Es, stims = FR.stimuli(geom, contrast=p['c'])
        prm = toy_to_backend(p)
        err = 0
        for (exp, s), phi in T.stimuli().items():
            j = [i for i, st in enumerate(stims) if st['name'] == f'{exp}{s}'][0]
            toy = T.respond(phi, p)
            mine = B.readout_map(response(Es[j], geom, prm))
            err = max(err, np.abs(toy - mine).max() / np.abs(toy).max())
        worst = max(worst, err)
        log(f'   {name:18s} {err:.1e}')
    return worst


def part_d(Es, geom, stims):
    base = dict(B.DEFAULTS, kappa=1.0, s_deg=3.0, sigma=0.1, sigma_mode='rel')
    a = B.run(Es, geom, dict(base, site_hv='density', site_card='density', site_rad='density'))[0]
    b = B.run(Es, geom, dict(base, site_hv='gain', site_card='gain', site_rad='gain'))[0]
    err = np.abs(a - b).max() / np.abs(a).max()
    log(f'D. p = q = 1: density vs gain site, max relative difference in wedge responses: {err:.1e}')
    return err


def part_e(E4):
    g1, g4 = FP.PyramidGeometry(CH4, ds=1), FP.PyramidGeometry(CH4, ds=4)
    keys = [(i, j) for i in range(2) for j in range(4)]
    E1 = [FP.CONTRAST ** 2 * E4[k] for k in keys]
    E4b = [FP.block_mean(e, 4) for e in E1]
    stims = [dict(exp='dg' if i == 0 else 'da', bar=FP.dominant_bar(E4b[n], g4))
             for n, (i, j) in enumerate(keys)]
    cases = {'sum, untuned, 3 deg, q 1': dict(kappa=0.0, s_deg=3.0),
             'sum, tuned, 1.5 deg, q 2 inside': dict(kappa=4.0, s_deg=1.5, q=2.0),
             'sum, local, q 2 outside': dict(s_deg=0.0, q=2.0, q_where='outside'),
             'spread, 3 deg': dict(pool='spread', s_deg=3.0),
             'sum, gated, 3 deg': dict(gate=True, s_deg=3.0),
             'spread, global, q 2 outside': dict(pool='spread', s_deg=np.inf, q=2.0, q_where='outside')}
    log('E. block averaging 4 x 4 (factorial grid) vs full resolution, 4-channel energies, 50% '
        'contrast, sigma 0.1 x mean pool, toy amplitudes; asymmetries in % of mean response:')
    worst = 0
    for name, prm in cases.items():
        prm = dict(prm, sigma=0.1, sigma_mode='rel', site_rad='output')
        a = B.asymmetries(B.run(E1, g1, prm)[0][0], stims)
        b = B.asymmetries(B.run(E4b, g4, prm)[0][0], stims)
        d = max(abs(a[e][k] - b[e][k]) for e in a for k in ('hv', 'co', 'rt', 'pp'))
        worst = max(worst, d)
        log(f'   {name:34s} max |difference| {100 * d:.2f} points; dg card-obl {100 * a["dg"]["co"]:+.1f} '
            f'vs {100 * b["dg"]["co"]:+.1f}, da rad-tang {100 * a["da"]["rt"]:+.1f} vs {100 * b["da"]["rt"]:+.1f}')
    return worst, stims


if __name__ == '__main__':
    t0 = time.time()
    z = np.load(os.path.join(FP.HERE, '..', 'paradoxes', '_cache', 'energy.npz'))
    E4 = {(i, j): z[f'e{i}{j}'].astype(np.float64) for i in range(2) for j in range(4)}
    res = dict(A=part_a(E4), B=part_b(E4), C=part_c())
    Es, stims, theta = FP.energies(8, 4)
    res['D'] = part_d(Es, FP.PyramidGeometry(theta, ds=4), stims)
    res['E'], st4 = part_e(E4)
    log(f'summary: {res}   total {time.time() - t0:.0f}s')
