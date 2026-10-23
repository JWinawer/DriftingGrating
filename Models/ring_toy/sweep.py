"""Parameter sweep of the ring toy: which components reproduce each paradox?

Writes sweep_results.csv, prints pass rates per component, and runs ablations from the
simplest configuration that reproduces every target. Takes about a minute.
"""
import csv
import itertools
import time

import numpy as np

from ring_toy import DEFAULTS, DATA, contrasts, targets

GRID = dict(
    rad_site=['density', 'output'],
    b_card=[0.0, 0.3],
    kappa_pool=[0.0, 1.0, 4.0],
    s_sp=[5.0, 15.0, 45.0, np.inf],
    q=[1.0, 1.5, 2.0],
    gate=[0.0, 2.0],
    a_pool=[0.0, 0.3],
    sigma=[1.0, 0.3, 0.1, 0.03, 0.01],
)
KEYS = list(GRID)
ASYM = ['horiz-vert', 'card-obl', 'rad-tang', 'polc-polo']


def run(p):
    c = contrasts(p)
    row = {k: p[k] for k in KEYS}
    row.update({f'{e}:{a}': c[e][a] for e in ('dg', 'da') for a in ASYM})
    t = targets(c)
    row.update(t)
    row['all'] = all(t.values())
    return row


def main():
    combos = list(itertools.product(*GRID.values()))
    rows, t0 = [], time.time()
    for i, vals in enumerate(combos):
        rows.append(run(dict(DEFAULTS, **dict(zip(KEYS, vals)))))
        if (i + 1) % 500 == 0 or i + 1 == len(combos):
            el = time.time() - t0
            print(f'{i + 1}/{len(combos)} ({100 * (i + 1) / len(combos):.0f}%), '
                  f'{el:.0f}s elapsed, ~{el / (i + 1) * (len(combos) - i - 1):.0f}s left', flush=True)

    with open('sweep_results.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    tnames = [k for k in rows[0] if k.startswith('P')] + ['all']
    print(f'\n{len(rows)} runs; pass rate per target, overall:')
    for t in tnames:
        print(f'  {t:30s} {np.mean([r[t] for r in rows]):.2f}')

    print('\nPass rate of "all targets" by component level:')
    for k in KEYS:
        print(f'  {k:11s}', '  '.join(f'{v}: {np.mean([r["all"] for r in rows if r[k] == v]):.3f}'
                                     for v in GRID[k]))

    passing = [r for r in rows if r['all']]
    print(f'\n{len(passing)} runs reproduce every target.')
    if not passing:
        return
    print('Levels never seen among passing runs (i.e. they preclude the full pattern):')
    for k in KEYS:
        absent = [v for v in GRID[k] if not any(r[k] == v for r in passing)]
        if absent:
            print(f'  {k}: {absent}')

    # simplest passing run: fewest optional components switched on, then closest to the data
    def complexity(r):
        return (r['kappa_pool'] > 0) + (r['gate'] > 0) + (r['a_pool'] > 0) + (r['b_card'] > 0) \
            + (r['rad_site'] == 'output') + (r['q'] > 1) + (not np.isinf(r['s_sp']))

    def dist(r):
        return sum((r[f'{e}:{a}'] - DATA[e][a]) ** 2 for e in DATA for a in ASYM)

    best = min(passing, key=lambda r: (complexity(r), dist(r)))
    print('\nSimplest passing configuration:', {k: best[k] for k in KEYS})
    report(best, 'simplest passing')

    base = dict(DEFAULTS, **{k: best[k] for k in KEYS})
    print('\nAblations from it (which targets break):')
    off = dict(q=1.0, gate=0.0, a_pool=0.0, b_card=0.0, kappa_pool=0.0, s_sp=np.inf,
               rad_site='density', a_card=0.0, a_rad=0.0, a_hv=0.0)
    for k, v in off.items():
        if base[k] == v:
            continue
        r = run(dict(base, **{k: v}))
        broken = [t for t in tnames[:-1] if not r[t]]
        print(f'  {k} -> {v}: breaks {broken or "nothing"}')
    for s in GRID['sigma']:
        r = run(dict(base, sigma=s))
        print(f'  sigma -> {s}: breaks {[t for t in tnames[:-1] if not r[t]] or "nothing"}')


def report(r, label):
    print(f'  {label}:')
    for e in ('dg', 'da'):
        print(f'    {e}  model ' + '  '.join(f'{a}={r[f"{e}:{a}"]:+.3f}' for a in ASYM))
        print(f'    {e}  data  ' + '  '.join(f'{a}={DATA[e][a]:+.3f}' for a in ASYM))


if __name__ == '__main__':
    main()
