"""Do the conclusions depend on the size of the density biases? Strict pass rates over the
mechanism grid of sweep.py (restricted to the components that matter), for several bias sizes."""
import itertools
import numpy as np
from ring_toy import DEFAULTS, contrasts, targets


def strict_pass(c):
    return all(targets(c).values()) and c['da']['card-obl'] > 0.5 * c['dg']['card-obl'] \
        and c['da']['rad-tang'] > 1.5 * c['dg']['rad-tang']


mech = dict(rad_site=['density', 'output'], kappa_pool=[0.0, 4.0], s_sp=[5.0, 45.0],
            q=[1.0, 2.0], sigma=[0.3, 0.1, 0.03, 0.01])
keys = list(mech)
for a_card, a_rad, a_hv in itertools.product([0.05, 0.15, 0.3], [0.05, 0.15, 0.3], [0.05, 0.1, 0.2]):
    if a_card + a_rad + a_hv >= 0.95:
        continue
    res = {}
    for vals in itertools.product(*mech.values()):
        p = dict(DEFAULTS, a_card=a_card, a_rad=a_rad, a_hv=a_hv, **dict(zip(keys, vals)))
        res[vals] = strict_pass(contrasts(p))
    ok = [v for v, r in res.items() if r]
    needed = {k: sorted({v[i] for v in ok}) for i, k in enumerate(keys)} if ok else {}
    print(f'a_card={a_card:.2f} a_rad={a_rad:.2f} a_hv={a_hv:.2f}: {len(ok):2d} strict passes;',
          'levels present: ' + '; '.join(f'{k}={v}' for k, v in needed.items() if k != 'sigma'))
