"""Step 3: one factorial over the back end's switches, run identically on each front end.

Front ends:
  ring36  the ring toy as written (36 orientation channels), on a ring at 6 deg
  ring8   the same ring with 8 channels, to separate orientation resolution from everything else
  pyr8    steerable-pyramid energy of the model stimuli, 8 channels, 50% contrast (4x4 blocks)

Switches (every combination):
  site_hv, site_card, site_rad   density | gain | output
  kappa                          0 (untuned) | 1 | 4           pool orientation weighting
  s_deg                          0 (local) | 1.5 | 3 | 6 | global
  exponent                       p1q1 | p1q2in | p1q2out | p2q2 | gate_in | gate_out
  pool                           sum | spread
  sigma_rel                      0.03 | 0.1 | 0.3 | 1 | 3  (x the mean pool over the band)
plus no normalization at every site combination. Bias amplitudes are the ring toy's: horizontal
1 +- 0.10, cardinal 1 +- 0.15, radial 1 +- 0.15; components at one site multiply.

sigma is set relative to the mean pool so that "strong normalization" means the same thing on
every front end, exponent and pool type. Because of that, stimulus contrast only rescales sigma
(exactly, except for the gated exponents).

Usage: python factorial.py FRONTEND [--quick] [--workers N]
Writes _cache/factorial_FRONTEND.csv and logs progress to _cache/factorial_FRONTEND.log.
"""
import argparse
import csv
import itertools
import multiprocessing as mp
import os
import time

import numpy as np

import backend as B

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, '_cache')

EXPONENTS = {
    'p1q1': dict(p=1.0, q=1.0),
    'p1q2in': dict(p=1.0, q=2.0, q_where='inside'),
    'p1q2out': dict(p=1.0, q=2.0, q_where='outside'),
    'p2q2': dict(p=2.0, q=2.0, q_where='inside'),
    'gate_in': dict(p=1.0, gate=True, q_where='inside'),
    'gate_out': dict(p=1.0, gate=True, q_where='outside'),
}
GRID = dict(site_hv=list(B.SITES), site_card=list(B.SITES), site_rad=list(B.SITES),
            kappa=[0.0, 1.0, 4.0], s_deg=[0.0, 1.5, 3.0, 6.0, np.inf],
            exponent=list(EXPONENTS), pool=['sum', 'spread'])
SIGMAS = [0.03, 0.1, 0.3, 1.0, 3.0]
SWITCHES = list(GRID)
ASYM = ('hv', 'co', 'rt', 'pp')


def configs(quick=False):
    out = [dict(zip(SWITCHES, v)) for v in itertools.product(*GRID.values())]
    nonorm = [dict(site_hv=a, site_card=b, site_rad=c, kappa=np.nan, s_deg=np.nan,
                   exponent='none', pool='none')
              for a, b, c in itertools.product(B.SITES, repeat=3)]
    out = nonorm + out
    if quick:
        rng = np.random.default_rng(0)
        out = [out[i] for i in rng.choice(len(out), 40, replace=False)]
    return out


def to_params(cfg):
    prm = {k: cfg[k] for k in ('site_hv', 'site_card', 'site_rad')}
    if cfg['pool'] == 'none':
        return dict(prm, pool='none')
    prm.update(kappa=cfg['kappa'], s_deg=cfg['s_deg'], pool=cfg['pool'], sigma_mode='rel')
    prm.update(EXPONENTS[cfg['exponent']])
    return prm


def load(frontend):
    if frontend.startswith('ring'):
        import frontend_ring as FR
        geom = FR.RingGeometry(72, int(frontend[4:]))
        Es, stims = FR.stimuli(geom, contrast=0.5)
    else:
        import frontend_pyramid as FP
        Es, stims, theta = FP.energies(int(frontend[3:]), 4)
        geom = FP.PyramidGeometry(theta, ds=4)
    return Es, stims, geom


_STATE = {}


def _init(frontend):
    try:
        import torch
        torch.set_num_threads(1)
    except ImportError:
        pass
    _STATE['data'] = load(frontend)
    _STATE['hcache'] = {}


def evaluate(cfg):
    Es, stims, geom = _STATE['data']
    resp, abs_s = B.run(Es, geom, to_params(cfg), sigmas=SIGMAS, hcache=_STATE['hcache'])
    rows = []
    for i, s in enumerate(SIGMAS):
        if cfg['pool'] == 'none' and i > 0:
            break
        c = B.asymmetries(resp[i], stims)
        row = dict(cfg, sigma_rel=np.nan if cfg['pool'] == 'none' else s, sigma_abs=abs_s[i])
        for e in ('dg', 'da'):
            row.update({f'{e}_{a}': c[e][a] for a in ASYM})
            row.update({f'{e}_raw_{a}': c[e]['raw_' + a] for a in ASYM})
            row[f'{e}_mean'] = c[e]['mean']
        row.update(B.score(c))
        rows.append(row)
    return rows


def log(msg, path):
    line = time.strftime('%H:%M:%S ') + msg
    print(line, flush=True)
    with open(path, 'a') as f:
        f.write(line + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('frontend', choices=['ring36', 'ring8', 'pyr8'])
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--workers', type=int, default=max(1, os.cpu_count() - 2))
    a = ap.parse_args()
    os.makedirs(CACHE, exist_ok=True)
    tag = a.frontend + ('_quick' if a.quick else '')
    logpath = os.path.join(CACHE, f'factorial_{tag}.log')
    open(logpath, 'w').close()
    cfgs = configs(a.quick)
    log(f'{a.frontend}: {len(cfgs)} configurations x {len(SIGMAS)} sigma, {a.workers} workers', logpath)
    rows, t0 = [], time.time()
    step = max(1, len(cfgs) // 50)
    with mp.get_context('spawn').Pool(a.workers, initializer=_init, initargs=(a.frontend,)) as pool:
        for i, r in enumerate(pool.imap_unordered(evaluate, cfgs, chunksize=4)):
            rows += r
            if (i + 1) % step == 0 or i + 1 == len(cfgs):
                el = time.time() - t0
                log(f'{i + 1}/{len(cfgs)} ({100 * (i + 1) / len(cfgs):.0f}%), {el:.0f}s elapsed, '
                    f'~{el / (i + 1) * (len(cfgs) - i - 1):.0f}s left', logpath)
    with open(os.path.join(CACHE, f'factorial_{tag}.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    log(f'done: {len(rows)} rows in {time.time() - t0:.0f}s', logpath)


if __name__ == '__main__':
    main()
