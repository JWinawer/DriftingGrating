"""Figures and ablations for the ring toy, using the stricter (beyond-geometry) targets.

Run after sweep.py. Writes figures/*.png and prints the ablation table.
"""
import csv
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from ring_toy import DEFAULTS, DATA, contrasts, targets

ASYM = ['horiz-vert', 'card-obl', 'rad-tang', 'polc-polo']
KEYS = ['rad_site', 'b_card', 'kappa_pool', 's_sp', 'q', 'gate', 'a_pool', 'sigma']
os.makedirs('figures', exist_ok=True)


def strict(c):
    """All qualitative targets, plus context effects larger than wedge geometry alone gives."""
    t = targets(c)
    t['P3 strict: da card-obl < half of dg'] = c['da']['card-obl'] > 0.5 * c['dg']['card-obl']
    t['P3 strict: da rad-tang > 1.5x dg'] = c['da']['rad-tang'] > 1.5 * c['dg']['rad-tang']
    return t


def parse(v):
    try:
        return float(v)
    except ValueError:
        return v


# simplest strict-passing configuration from the sweep: standard exponent, no gate, no extra pool
# anisotropy, no tuning-width anisotropy; among those, closest to the data
rows = list(csv.DictReader(open('sweep_results.csv')))
cands = []
for x in rows:
    p = dict(DEFAULTS, **{k: parse(x[k]) for k in KEYS})
    if x['all'] == 'True':
        c = contrasts(p)
        if all(strict(c).values()):
            extras = (p['q'] > 1) + (p['gate'] > 0) + (p['a_pool'] > 0) + (p['b_card'] > 0)
            d = sum((c[e][a] - DATA[e][a]) ** 2 for e in DATA for a in ASYM)
            cands.append((extras, d, p, c))
print(f'{len(cands)} strict-passing runs; optional extras used by the simplest: '
      f'{min(t[0] for t in cands)}')
_, _, best, cbest = min(cands, key=lambda t: (t[0], t[1]))
print('Simplest strict-passing model:', {k: best[k] for k in KEYS})

print('\nAblations (targets that break):')
abl = [('no normalization', dict(norm=False)),
       ('untuned pool (kappa_pool 0)', dict(kappa_pool=0.0)),
       ('weakly tuned pool (kappa_pool 1)', dict(kappa_pool=1.0)),
       ('local pool (s_sp 5 deg)', dict(s_sp=5.0)),
       ('radial bias feeds the pool', dict(rad_site='density')),
       ('weak normalization (sigma 1)', dict(sigma=1.0)),
       ('no cardinal density bias', dict(a_card=0.0)),
       ('no horizontal density bias', dict(a_hv=0.0)),
       ('no radial bias', dict(a_rad=0.0)),
       ('standard exponent (q 1)', dict(q=1.0)),
       ('add homogeneity gate', dict(gate=2.0))]
ablated = {}
for name, change in abl:
    c = contrasts(dict(best, **change))
    ablated[name] = c
    broken = [k for k, v in strict(c).items() if not v]
    print(f'  {name:34s} breaks: {broken or "nothing"}')

# Figure 1: model vs data, for geometry only, the standard untuned pool, and the best model
models = [('Data (V1, 4-8 deg)', DATA),
          ('No normalization\n(geometry + density)', ablated['no normalization']),
          ('Untuned pool', ablated['untuned pool (kappa_pool 0)']),
          ('Radial bias feeds pool', ablated['radial bias feeds the pool']),
          ('Best: tuned, broad pool;\nradial gain after norm.', cbest)]
fig, axs = plt.subplots(1, len(models), figsize=(17, 3.6), sharey=False)
x = np.arange(len(ASYM))
for ax, (title, c) in zip(axs, models):
    ax.bar(x - 0.2, [c['dg'][a] for a in ASYM], 0.4, label='dg (Cartesian)', color='#4C72B0')
    ax.bar(x + 0.2, [c['da'][a] for a in ASYM], 0.4, label='da (polar)', color='#DD8452')
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xticks(x, ASYM, rotation=30, fontsize=8)
    ax.set_title(title, fontsize=9)
axs[0].legend(fontsize=7, frameon=False)
axs[0].set_ylabel('asymmetry')
fig.text(0.5, -0.04, 'Data are BOLD estimates; model values are (A - B) / mean response. Compare signs '
         'and dg-vs-da pattern, not magnitudes.', ha='center', fontsize=8)
fig.tight_layout()
fig.savefig('figures/fig1_model_vs_data.png', dpi=150, bbox_inches='tight')

# Figure 2: pool orientation tuning x spatial extent, at the best sigma, radial gain after norm.
kps = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
sps = [5.0, 10.0, 15.0, 30.0, 45.0, 90.0, np.inf]
grids = {k: np.zeros((len(kps), len(sps))) for k in ('dg card-obl', 'da/dg card-obl', 'da/dg rad-tang', 'pass')}
for i, kp in enumerate(kps):
    for j, sp in enumerate(sps):
        c = contrasts(dict(best, kappa_pool=kp, s_sp=sp))
        grids['dg card-obl'][i, j] = c['dg']['card-obl']
        grids['da/dg card-obl'][i, j] = c['da']['card-obl'] / c['dg']['card-obl']
        grids['da/dg rad-tang'][i, j] = c['da']['rad-tang'] / c['dg']['rad-tang']
        grids['pass'][i, j] = all(strict(c).values())
fig, axs = plt.subplots(1, 4, figsize=(17, 3.8))
spec = [('dg card-obl', 'RdBu_r', 'dg cardinal - oblique\n(data: negative)'),
        ('da/dg card-obl', 'viridis', 'da / dg, cardinal - oblique\n(data 0.23; geometry alone 0.64)'),
        ('da/dg rad-tang', 'viridis', 'da / dg, radial - tangential\n(data 2.2; geometry alone 1.1)'),
        ('pass', 'Greys', 'all targets met (black)')]
for ax, (k, cmap, title) in zip(axs, spec):
    g = grids[k]
    if k == 'dg card-obl':
        v = np.abs(g).max()
        im = ax.imshow(g, cmap=cmap, vmin=-v, vmax=v, origin='lower', aspect='auto')
    elif k == 'pass':
        im = ax.imshow(g, cmap=cmap, vmin=0, vmax=1, origin='lower', aspect='auto')
    else:
        im = ax.imshow(np.clip(g, -1, 3), cmap=cmap, origin='lower', aspect='auto')
    if k != 'pass':
        fig.colorbar(im, ax=ax, fraction=0.046)
    ax.set_xticks(range(len(sps)), ['global' if np.isinf(s) else f'{s:g}' for s in sps], fontsize=8)
    ax.set_yticks(range(len(kps)), [f'{k:g}' for k in kps], fontsize=8)
    ax.set_xlabel('pool spatial s.d. (deg polar angle)')
    ax.set_title(title, fontsize=9)
axs[0].set_ylabel('pool orientation tuning (kappa; 0 = untuned)')
fig.tight_layout()
fig.savefig('figures/fig2_pool_tuning_by_extent.png', dpi=150, bbox_inches='tight')

# Figure 3: normalization strength, for the best model with the radial bias after vs. in the pool
sigmas = np.logspace(-2.5, 0.5, 25)
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
for ax, site in zip(axs, ('output', 'density')):
    cs = [contrasts(dict(best, sigma=s, rad_site=site)) for s in sigmas]
    for e, ls in (('dg', '-'), ('da', '--')):
        ax.semilogx(sigmas, [c[e]['card-obl'] for c in cs], ls, color='#4C72B0', label=f'{e} card-obl')
        ax.semilogx(sigmas, [c[e]['rad-tang'] for c in cs], ls, color='#C44E52', label=f'{e} rad-tang')
    ax.axhline(0, color='k', lw=0.8)
    ax.axvline(best['sigma'], color='grey', lw=0.8, ls=':')
    ax.set_xlabel('sigma (smaller = stronger normalization)')
    ax.set_title('radial bias ' + ('applied after normalization' if site == 'output' else 'feeds the pool'),
                 fontsize=9)
axs[0].set_ylabel('asymmetry')
axs[0].legend(fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig('figures/fig3_normalization_strength.png', dpi=150, bbox_inches='tight')
print('\nFigures written to figures/.')
print('Best model values:')
for e in ('dg', 'da'):
    print(f'  {e} model ' + '  '.join(f'{a}={cbest[e][a]:+.3f}' for a in ASYM))
    print(f'  {e} data  ' + '  '.join(f'{a}={DATA[e][a]:+.3f}' for a in ASYM))
