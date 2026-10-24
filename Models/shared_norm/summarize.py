"""Step 4: the matched table (component x front end -> paradoxes produced or ruled out).

Reads _cache/factorial_{ring36,ring8,pyr8}.csv. Writes results/tables.md (generated; the
plain-language reading is in RESULTS.md) and figures/matched_table.png.

Paradox groups, each a conjunction of the shared targets in backend.TARGETS:
  P1      dg cardinal < oblique
  P2      radial > tangential in both experiments
  P3card  P1, and the dg - da context effect has the data's sign, and da is less than half of dg
  P3rad   P2, and da > dg, and da > 1.5 x dg
  HV      horizontal < vertical in dg and in da, and more negative in dg (scored separately)
  core    P1, P2, P3card and P3rad together (figures.py's strict targets without HV)
  all     core and HV
"""
import os

import numpy as np
import pandas as pd

import backend as B
from factorial import GRID, SIGMAS

HERE = os.path.dirname(os.path.abspath(__file__))
FRONTENDS = ['ring36', 'ring8', 'pyr8']
SWITCHES = list(GRID) + ['sigma_rel']
GROUPS = ['P1', 'P2', 'P3card', 'P3rad', 'HV', 'core', 'all']
LEVEL_NAMES = {('s_deg', 'inf'): 'global', ('s_deg', '0.0'): 'local (0)', ('kappa', '0.0'): 'untuned (0)'}


def add_groups(d, prefix=''):
    t = lambda k: d[prefix + k]
    d[prefix + 'P1'] = t('P1 dg card<obl')
    d[prefix + 'P2'] = t('P2 rad>tang both')
    d[prefix + 'P3card'] = d[prefix + 'P1'] & t('P3 card context dg<da') & t('P3 strict card: da > dg/2')
    d[prefix + 'P3rad'] = d[prefix + 'P2'] & t('P3 rad context da>dg') & t('P3 strict rad: da > 1.5 dg')
    d[prefix + 'HV'] = t('HV dg h<v') & t('HV da h<v') & t('HV context dg<da')
    d[prefix + 'core'] = d[prefix + 'P3card'] & d[prefix + 'P3rad']
    d[prefix + 'all'] = d[prefix + 'core'] & d[prefix + 'HV']
    return d


def common_mean_targets(d):
    """The same targets with every asymmetry divided by the mean over BOTH experiments, so that
    a difference in overall response level between dg and da is kept rather than divided out."""
    m = (d['dg_mean'] + d['da_mean']) / 2
    rows = []
    for _, r in d.iterrows():
        c = {e: {a: r[f'{e}_raw_{a}'] / m[_] for a in ('hv', 'co', 'rt', 'pp')} for e in ('dg', 'da')}
        rows.append({'cm_' + k: v for k, v in B.score(c).items()})
    return add_groups(pd.DataFrame(rows, index=d.index), 'cm_')


def load():
    out = {}
    for fe in FRONTENDS:
        d = pd.read_csv(os.path.join(HERE, '_cache', f'factorial_{fe}.csv'))
        for k in B.TARGETS:
            d[k] = d[k].astype(bool)
        d = add_groups(d)
        d = pd.concat([d, common_mean_targets(d)], axis=1)
        for s in SWITCHES:
            d[s] = d[s].astype(str)
        out[fe] = d
    return out


def level_table(D):
    """For every switch level: % of normalized configurations that pass each group."""
    rows = []
    for s in SWITCHES:
        levels = GRID.get(s, SIGMAS)
        for lv in levels:
            lv = str(lv)
            row = dict(switch=s, level=LEVEL_NAMES.get((s, lv), lv))
            for fe, d in D.items():
                n = d[d.pool != 'none']
                sub = n[n[s] == lv]
                for g in GROUPS:
                    p = sub[g].mean() * 100
                    others = n[(n[s] != lv) & n[g]]
                    required = sub[g].any() and others.empty
                    row[f'{fe}:{g}'] = 'req' if required else ('✗' if p == 0 else f'{p:.0f}' if p >= 1 else '<1')
            rows.append(row)
    return pd.DataFrame(rows)


def agreement(D):
    """Do two front ends give the same verdict for the same switch settings?"""
    key = SWITCHES
    lines = []
    for a, b in (('ring36', 'ring8'), ('ring8', 'pyr8'), ('ring36', 'pyr8')):
        m = D[a].merge(D[b], on=key, suffixes=('_a', '_b'))
        m = m[m.pool != 'none']
        g = '  '.join(f'{x} {100 * (m[x + "_a"] == m[x + "_b"]).mean():.0f}%' for x in GROUPS)
        r = '  '.join(f'{e} {a_}: {np.corrcoef(m[f"{e}_{a_}_a"], m[f"{e}_{a_}_b"])[0, 1]:+.2f}'
                      for e in ('dg', 'da') for a_ in ('hv', 'co', 'rt'))
        both = (m['core_a'] & m['core_b']).sum()
        lines.append(f'| {a} vs {b} | {g} | {r} | {both} |')
    return lines


def complexity(r):
    return ((r['kappa'] != '0.0') + (r['exponent'] != 'p1q1') + (r['pool'] != 'sum')
            + sum(r[f'site_{c}'] != 'density' for c in ('hv', 'card', 'rad'))
            + (r['s_deg'] == '0.0'))


def shared_passers(D, group='core'):
    m = None
    for fe in FRONTENDS:
        d = D[fe][D[fe][group]][SWITCHES]
        m = d if m is None else m.merge(d, on=SWITCHES)
    return m


def broken(r):
    return [g for g in GROUPS[:5] if not r[g]]


def ablations(D, ref):
    """Change one switch at a time from ref; list the groups that fail on each front end."""
    out = []
    def row(cfg):
        res = []
        for fe in FRONTENDS:
            d = D[fe]
            sel = np.logical_and.reduce([d[s] == cfg[s] for s in SWITCHES])
            r = d[sel]
            res.append(', '.join(broken(r.iloc[0])) or 'none' if len(r) else 'n/a')
        return res
    out.append(('reference', *row(ref)))
    for s in SWITCHES:
        for lv in GRID.get(s, SIGMAS):
            lv = str(lv)
            if lv == ref[s]:
                continue
            out.append((f'{s} = {LEVEL_NAMES.get((s, lv), lv)}', *row(dict(ref, **{s: lv}))))
    nn = dict(ref, kappa='nan', s_deg='nan', exponent='none', pool='none', sigma_rel='nan')
    out.append(('no normalization', *row(nn)))
    return out


def md_table(df):
    cols = list(df.columns)
    lines = ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    lines += ['| ' + ' | '.join(str(v) for v in r) + ' |' for r in df.itertuples(index=False)]
    return lines


def figure(T, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    groups = GROUPS
    cols = [f'{fe}:{g}' for g in groups for fe in FRONTENDS]
    def val(v):
        return 100.0 if v == 'req' else 0.0 if v == '✗' else 0.5 if v == '<1' else float(v)
    M = np.array([[val(v) for v in T[c]] for c in cols]).T
    fig, ax = plt.subplots(figsize=(14, 0.32 * len(T) + 1.5))
    im = ax.imshow(np.log10(M + 0.3), cmap='viridis', aspect='auto', vmin=np.log10(0.3), vmax=2)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = T[cols[j]].iloc[i]
            ax.text(j, i, v, ha='center', va='center', fontsize=6,
                    color='white' if M[i, j] < 20 else 'black')
    ax.set_yticks(range(len(T)), [f'{s}: {l}' for s, l in zip(T.switch, T.level)], fontsize=7)
    ax.set_xticks(range(len(cols)), [c.split(':')[0] for c in cols], rotation=90, fontsize=6)
    for k, g in enumerate(groups):
        ax.text(3 * k + 1, -1.2, g, ha='center', fontsize=9, fontweight='bold')
        if k:
            ax.axvline(3 * k - 0.5, color='w', lw=2)
    prev = None
    for i, s in enumerate(T.switch):
        if prev is not None and s != prev:
            ax.axhline(i - 0.5, color='w', lw=1)
        prev = s
    ax.set_title('% of configurations with this setting that reproduce each paradox group\n'
                 '(✗ = never: the setting rules the group out; req = only setting of its switch that ever does)',
                 fontsize=9, pad=24)
    fig.colorbar(im, ax=ax, fraction=0.02, label='log10(% + 0.3)')
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches='tight')


def agreement_figure(D, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    m = D['ring8'].merge(D['pyr8'], on=SWITCHES, suffixes=('_r', '_p'))
    m = m[m.pool != 'none'].sample(4000, random_state=0)
    names = dict(hv='horizontal - vertical', co='cardinal - oblique', rt='radial - tangential')
    fig, axs = plt.subplots(2, 3, figsize=(11, 7))
    for i, e in enumerate(('dg', 'da')):
        for j, a in enumerate(('hv', 'co', 'rt')):
            ax = axs[i, j]
            x, y = 100 * m[f'{e}_{a}_r'], 100 * m[f'{e}_{a}_p']
            ax.scatter(x, y, s=2, alpha=0.3, color='#4C72B0')
            lim = [min(x.min(), y.min()), max(x.max(), y.max())]
            ax.plot(lim, lim, 'k', lw=0.7)
            ax.axhline(0, color='grey', lw=0.5); ax.axvline(0, color='grey', lw=0.5)
            ax.set_title(f'{e}: {names[a]}  (r = {np.corrcoef(x, y)[0, 1]:.2f})', fontsize=9)
            ax.set_xlabel('ring, 8 channels (% of mean)', fontsize=8)
            ax.set_ylabel('pyramid, 8 channels (% of mean)', fontsize=8)
    fig.suptitle('Same switch settings, two front ends: each dot is one configuration (4,000 of 24,300)', fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches='tight')


def main():
    D = load()
    T = level_table(D)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
    agreement_figure(D, os.path.join(HERE, 'figures', 'agreement.png'))
    figure(T, os.path.join(HERE, 'figures', 'matched_table.png'))
    L = ['# Shared-normalization factorial: generated tables', '',
         'Generated by `summarize.py` from `_cache/factorial_*.csv`. Read `RESULTS.md` (this folder) for '
         'what these mean. Percentages are over all configurations with normalization on.', '']
    L += ['## Counts', '', '| front end | configurations | ' + ' | '.join(GROUPS)
          + ' | core, common-mean readout | all, common-mean readout |', '|' + '---|' * (len(GROUPS) + 4)]
    for fe, d in D.items():
        n = d[d.pool != 'none']
        L.append(f'| {fe} | {len(n)} | ' + ' | '.join(str(int(n[g].sum())) for g in GROUPS)
                 + f' | {int(n["cm_core"].sum())} | {int(n["cm_all"].sum())} |')
    L += ['', '## Matched table: % passing, by switch level', '']
    for g in GROUPS:
        sub = T[['switch', 'level'] + [f'{fe}:{g}' for fe in FRONTENDS]]
        sub.columns = ['switch', 'level'] + FRONTENDS
        L += [f'### {g}', ''] + md_table(sub) + ['']
    L += ['## Agreement between front ends, same switch settings', '',
          '| pair | same verdict, per group | correlation of asymmetry values across configurations | both pass core |',
          '|---|---|---|---|'] + agreement(D) + ['']
    sp = shared_passers(D, 'core')
    L += [f'## Configurations that pass core on all three front ends: {len(sp)}', '']
    if len(sp):
        sp = sp.assign(cx=sp.apply(complexity, axis=1)).sort_values('cx')
        L += md_table(sp.head(25)) + ['']
        ref = sp.iloc[0][SWITCHES].to_dict()
    else:
        ref = None
    spa = shared_passers(D, 'all')
    L += [f'Configurations that pass all (core and HV) on all three: {len(spa)}', '']
    for name, r in [('simplest configuration passing core on all three front ends', ref),
                    ('the ring toy\'s best model (sites density / density / output, kappa 4, 1.5 deg, '
                     'q 2 outside, sum)', dict(site_hv='density', site_card='density', site_rad='output',
                                                kappa='4.0', s_deg='1.5', exponent='p1q2out', pool='sum',
                                                sigma_rel='0.1')),
                    ('PARADOXES.md\'s best model (cardinal and horizontal as input gain, radial on the '
                     'output, untuned spread pool, 3 deg, q 1)',
                     dict(site_hv='gain', site_card='gain', site_rad='output', kappa='0.0', s_deg='3.0',
                          exponent='p1q1', pool='spread', sigma_rel='0.1'))]:
        if r is None:
            continue
        L += [f'## Ablations from the {name}', '', '`' + str(r) + '`', '',
              '| change | ' + ' | '.join(FRONTENDS) + ' |', '|---|---|---|---|']
        L += [f'| {a} | {b} | {c} | {d} |' for a, b, c, d in ablations(D, r)] + ['']
    with open(os.path.join(HERE, 'results', 'tables.md'), 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('\n'.join(L[:12]))


if __name__ == '__main__':
    main()
