"""
Is the green channel reporting tagged receptor, or is it reporting the tissue?

This exists because the three-arm argument in v2_ish_arms.py rests entirely on
one premise -- that ex vivo the SEP channel reports TOTAL GluA1 -- and the arm
built on it behaved strangely: `SEP / auto` correlates with Gria1 expression at
-0.11, while plain nano manages +0.62. Either the premise is wrong or the arm
is, and that has to be settled on the channels themselves, before any ratio is
taken.

So this looks at the three registered channels with no denominators anywhere:
one mean per structure per adult, in log2, straight from the per-mouse files.
Three things decide it.

  dynamic range   how much a channel varies across the brain. Total receptor
                  is not flat -- hippocampus against thalamus is a large,
                  well known difference -- so a channel reporting it should
                  have a range comparable to the nanobody's.
  who it tracks   if the green channel is tag, its spatial profile should look
                  like the nanobody's and like Gria1's. If it is tissue, it
                  should look like the autofluorescence channel's.
  what is left    log2(SEP) regressed on log2(auto), and the residual asked
                  whether it still carries a receptor map. This is the
                  constructive part: if the tag survives underneath the
                  autofluorescence, an unmixed SEP could still serve as the
                  total-receptor arm, and this says how much is there.

Read the output as a property of THIS tissue, fixed, cleared and mounted --
not of SEP-GluA1 mice. Superecliptic pHluorin is quenched in acidic
compartments in a living cell, which is what makes it a surface reporter in
vivo; after fixation and clearing the pH gradients are gone, and whatever
green emission survives the protocol is what these volumes contain.

Outputs, under data\\adult_v2\\arms:
  sep_channel_check.csv     the per-mouse numbers behind the figure
  sep_channel_check.png     the figure

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe mapping\\v2_sep_channel_check.py
"""

import csv
import math
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from v2_per_mouse import annotation_20, MICE, CSV_MAP, OUT as PER_MOUSE, DATA
from v2_cohort import NAIVE, RWS

OUT = os.path.join(DATA, 'adult_v2', 'arms')
GENES = os.path.join(DATA, 'adult_v2', 'ish', 'gene_region_table.csv')

ADULTS = NAIVE + RWS
CHANNELS = ('sig', 'auto', 'sep')
NICE = {'sig': 'nano', 'auto': 'autofluo', 'sep': 'SEP (green)'}
MIN_VOX = 250
MIN_ISH_VOXELS = 10
CONTROL = 'Gria1'


def structure_names():
    names = {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['parcellation_term_set_name'] == 'structure':
                names[int(row['parcellation_index'])] = row['parcellation_term_name']
    return names


def mouse_channels(mouse, names):
    """{channel: {structure: log2 mean}} for one adult, raw, no denominators."""
    z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
    tissue = z['tissue']
    ann = annotation_20(MICE[mouse][1])
    lab = ann[tissue]
    n = np.bincount(lab, minlength=int(ann.max()) + 1)

    out = {}
    for key in CHANNELS:
        s = np.bincount(lab, weights=z[key].astype(np.float32)[tissue], minlength=len(n))
        acc = defaultdict(lambda: [0, 0.0])
        for i in np.nonzero(n)[0]:
            if i == 0:
                continue
            a = acc[names.get(int(i), f'id{i}')]
            a[0] += int(n[i]); a[1] += s[i]
        out[key] = {k: math.log2(t / c) for k, (c, t) in acc.items()
                    if c >= MIN_VOX and t > 0}
    return out


def gria1_profile():
    prof = {}
    with open(GENES, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            if r['symbol'] == CONTROL and int(r['n_voxels']) >= MIN_ISH_VOXELS:
                prof[r['structure']] = float(r['ish_mean'])
    return prof


def residual(y, x):
    """log2(SEP) with the part predicted by log2(auto) taken out, linearly."""
    a = np.column_stack([x, np.ones_like(x)])
    return y - a @ np.linalg.lstsq(a, y, rcond=None)[0]


def main():
    os.makedirs(OUT, exist_ok=True)
    names = structure_names()
    per = {}
    for mouse in ADULTS:
        per[mouse] = mouse_channels(mouse, names)
        print(f'{mouse:20s} {len(per[mouse]["sep"])} structures', flush=True)

    gria = gria1_profile()
    rows = []
    for mouse in ADULTS:
        ch = per[mouse]
        common = sorted(set.intersection(*[set(ch[k]) for k in CHANNELS]))
        v = {k: np.array([ch[k][s] for s in common]) for k in CHANNELS}
        withg = [s for s in common if s in gria]
        g = np.array([gria[s] for s in withg])
        idx = [common.index(s) for s in withg]
        res = residual(v['sep'], v['auto'])
        row = dict(mouse=mouse, n_structures=len(common),
                   range_nano=float(np.diff(np.percentile(v['sig'], [10, 90]))[0]),
                   range_auto=float(np.diff(np.percentile(v['auto'], [10, 90]))[0]),
                   range_sep=float(np.diff(np.percentile(v['sep'], [10, 90]))[0]),
                   rho_sep_auto=float(spearmanr(v['sep'], v['auto']).statistic),
                   rho_sep_nano=float(spearmanr(v['sep'], v['sig']).statistic),
                   rho_nano_auto=float(spearmanr(v['sig'], v['auto']).statistic),
                   rho_nano_gria=float(spearmanr(v['sig'][idx], g).statistic),
                   rho_auto_gria=float(spearmanr(v['auto'][idx], g).statistic),
                   rho_sep_gria=float(spearmanr(v['sep'][idx], g).statistic),
                   rho_sepresid_gria=float(spearmanr(res[idx], g).statistic),
                   rho_sepresid_nano=float(spearmanr(res, v['sig']).statistic))
        rows.append(row)

    path = os.path.join(OUT, 'sep_channel_check.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()}
                    for r in rows)
    print(f'\n{len(rows)} mice -> {path}')

    col = lambda k: np.array([r[k] for r in rows])
    say = lambda k: f'{np.mean(col(k)):+.3f} +- {np.std(col(k)):.3f}'
    print(f'\ndynamic range across structures, p90-p10 of log2, mean over {len(rows)} adults')
    for k, label in (('range_nano', 'nano'), ('range_auto', 'autofluo'), ('range_sep', 'SEP')):
        print(f'  {label:10s} {np.mean(col(k)):.2f} +- {np.std(col(k)):.2f}')
    print('\nwhat the green channel tracks (per mouse, over structures)')
    print(f'  SEP  ~ autofluo   {say("rho_sep_auto")}')
    print(f'  SEP  ~ nano       {say("rho_sep_nano")}')
    print(f'  nano ~ autofluo   {say("rho_nano_auto")}   <- the nano channel is its own thing')
    print(f'\nagainst {CONTROL} expression')
    for k, label in (('rho_nano_gria', 'nano'), ('rho_auto_gria', 'autofluo'),
                     ('rho_sep_gria', 'SEP'), ('rho_sepresid_gria', 'SEP minus autofluo')):
        print(f'  {label:20s} {say(k)}')
    print(f'\n  SEP minus autofluo, against nano: {say("rho_sepresid_nano")}')

    # ------------------------------------------------------------------ figure
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.0))
    rng = np.random.default_rng(0)

    ax = axes[0]
    for i, k in enumerate(('range_nano', 'range_auto', 'range_sep')):
        v = col(k)
        ax.scatter(np.full(len(v), i) + rng.uniform(-.1, .1, len(v)), v, s=18,
                   facecolor='#c0392b' if k == 'range_sep' else '0.6',
                   edgecolor='0.25', linewidth=0.4, zorder=2)
        ax.plot([i - .28, i + .28], [np.median(v)] * 2, color='0.15', lw=1.7, zorder=3)
    ax.set_xticks(range(3)); ax.set_xticklabels(['nano', 'autofluo', 'SEP'], fontsize=8)
    ax.set_ylabel('p90 - p10 across structures (log2)', fontsize=8)
    ax.set_ylim(bottom=0)
    ax.set_title('how much each channel varies\nacross the brain', fontsize=9)

    ax = axes[1]
    mouse = ADULTS[0]
    ch = per[mouse]
    common = sorted(set.intersection(*[set(ch[k]) for k in CHANNELS]))
    x = np.array([ch['auto'][s] for s in common]); y = np.array([ch['sep'][s] for s in common])
    ax.scatter(x, y, s=10, facecolor='0.55', edgecolor='0.25', linewidth=0.3)
    ax.set_xlabel('log2 autofluorescence', fontsize=8)
    ax.set_ylabel('log2 SEP', fontsize=8)
    ax.set_title(f'{mouse}\nSEP against autofluo, rho = '
                 f'{spearmanr(x, y).statistic:+.2f}', fontsize=9)

    ax = axes[2]
    y2 = np.array([ch['sig'][s] for s in common])
    ax.scatter(y2, y, s=10, facecolor='0.55', edgecolor='0.25', linewidth=0.3)
    ax.set_xlabel('log2 nano', fontsize=8)
    ax.set_ylabel('log2 SEP', fontsize=8)
    ax.set_title(f'SEP against nano, rho = {spearmanr(y2, y).statistic:+.2f}', fontsize=9)

    ax = axes[3]
    keys = ('rho_nano_gria', 'rho_auto_gria', 'rho_sep_gria', 'rho_sepresid_gria')
    labels = ['nano', 'autofluo', 'SEP', 'SEP minus\nautofluo']
    for i, k in enumerate(keys):
        v = col(k)
        ax.scatter(np.full(len(v), i) + rng.uniform(-.1, .1, len(v)), v, s=18,
                   facecolor='#c0392b' if 'sep' in k else '0.6',
                   edgecolor='0.25', linewidth=0.4, zorder=2)
        ax.plot([i - .28, i + .28], [np.median(v)] * 2, color='0.15', lw=1.7, zorder=3)
    ax.axhline(0, color='0.8', lw=0.7, zorder=0)
    ax.set_xticks(range(len(keys))); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(f'Spearman with {CONTROL} expression', fontsize=8)
    ax.set_title(f'if the green channel were total receptor,\n'
                 f'it would beat nano here', fontsize=9)

    for ax in axes:
        ax.tick_params(labelsize=7)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    fig.suptitle('The green channel in fixed, cleared tissue: one dot per adult, '
                 'structure means with no denominator anywhere', fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path = os.path.join(OUT, 'sep_channel_check.png')
    fig.savefig(path, dpi=200); plt.close(fig)
    print(f'\n{path}')


if __name__ == '__main__':
    main()
