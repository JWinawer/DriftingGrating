"""Sign table for the ring toy, in the same layout as Models/paradoxes fig1_sign_table.

Rows are model variants, columns the six numbers the three paradoxes turn on, with the V1 data
in its own row on top. Writes figures/fig4_sign_table.png.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from ring_toy import DEFAULTS, DATA, contrasts

BEST = dict(DEFAULTS, rad_site='output', kappa_pool=4.0, s_sp=15.0, q=2.0, sigma=0.03)
ROWS = [
    ('no normalization', dict(norm=False)),
    ('untuned global pool, q=1', dict(kappa_pool=0.0, s_sp=np.inf, q=1.0, rad_site='density')),
    ('untuned global pool, q=2', dict(kappa_pool=0.0, s_sp=np.inf, rad_site='density')),
    ('tuned pool, local (5 deg), q=2', dict(s_sp=5.0, rad_site='density')),
    ('tuned pool, broad, q=2', dict(rad_site='density')),
    ('tuned, broad, q=1 + radial on output', dict(q=1.0)),
    ('untuned, broad, q=2 + radial on output', dict(kappa_pool=0.0)),
    ('tuned, local, q=2 + radial on output', dict(s_sp=5.0)),
    ('best, weak normalization (sigma 1)', dict(sigma=1.0)),
    ('best + homogeneity gate', dict(gate=2.0)),
    ('best: tuned, broad, q=2 + radial on output', dict()),
]


def tests(c):
    """The six numbers, in % of mean response (data: % signal change)."""
    return {
        'dg card-obl': c['dg']['card-obl'],
        'da card-obl': c['da']['card-obl'],
        'dg rad-tang': c['dg']['rad-tang'],
        'da rad-tang': c['da']['rad-tang'],
        'context card-obl\n(dg - da)': c['dg']['card-obl'] - c['da']['card-obl'],
        'context rad-tang\n(dg - da)': c['dg']['rad-tang'] - c['da']['rad-tang'],
    }


data_t = tests(DATA)
cols = list(data_t)
M = np.array([[100 * tests(contrasts(dict(BEST, **ch)))[k] for k in cols] for _, ch in ROWS])
D = np.array([[data_t[k] for k in cols]])

fig, (axd, ax) = plt.subplots(2, 1, figsize=(10.5, 6.6), sharex=True,
                              gridspec_kw={'height_ratios': [1, len(ROWS)], 'hspace': 0.08})
# data row on its own colour scale: % signal change and % of model mean are not comparable
dlim = np.abs(D).max()
axd.imshow(D, cmap='RdBu_r', norm=TwoSlopeNorm(0, -dlim, dlim), aspect='auto')
for c in range(D.shape[1]):
    axd.text(c, 0, f'{D[0, c]:+.3f}', ha='center', va='center', fontsize=9, fontweight='bold',
             color='white' if abs(D[0, c]) > 0.6 * dlim else 'black')
axd.set_yticks([0])
axd.set_yticklabels(['Data: V1, 4-8 deg\n(% signal change)'], fontsize=9)
axd.tick_params(axis='x', bottom=False)

lim = np.abs(M).max()
ax.imshow(M, cmap='RdBu_r', norm=TwoSlopeNorm(0, -lim, lim), aspect='auto')
for r in range(M.shape[0]):
    for c in range(M.shape[1]):
        ok = np.sign(M[r, c]) == np.sign(D[0, c]) and abs(M[r, c]) > 0.5
        label = '0' if round(M[r, c]) == 0 else f'{M[r, c]:+.0f}'   # no "-0"
        ax.text(c, r, label + (' ✓' if ok else ''), ha='center', va='center',
                fontsize=9, fontweight='bold' if ok else 'normal',
                color='white' if abs(M[r, c]) > 0.75 * lim else 'black')
ax.axhline(len(ROWS) - 1.5, color='k', lw=1.5)          # set the best model apart
ax.set_xticks(range(len(cols)))
ax.set_xticklabels(cols, fontsize=8)
ax.set_yticks(range(len(ROWS)))
ax.set_yticklabels([r for r, _ in ROWS], fontsize=9)
ax.get_yticklabels()[-1].set_fontweight('bold')
axd.set_title('Data (top, % signal change) and ring-toy asymmetries (% of mean response).\n'
              'Unlisted settings as in the best model: kappa_pool 4, pool s.d. 15 deg, q 2, sigma 0.03.'
              '  ✓ = same sign as the data', fontsize=10)
fig.savefig('figures/fig4_sign_table.png', dpi=150, bbox_inches='tight')
print('wrote figures/fig4_sign_table.png')
