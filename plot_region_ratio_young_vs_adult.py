"""
SUPERSEDED with region_means_raw_per_mouse.py; the plot it makes is now
v2_region_plot.py, which also marks the P16 brain and applies FDR.

Per-mouse dot plot of the P20-vs-adult region comparison, for the cortical
areas that matter to the critical-period question plus the structures that
dominate the brain-wide picture.

Three rows, one per way of reading the same raw data (see
region_ratio_young_vs_adult.py for what each means):
  nano/auto, no reference     absolute-ish: nanobody per unit autofluorescence
  nano-bg  / isocortex        share of the cortex
  nano-bg  / subcortex-HPF-STR share of the subcortex excluding HPF and STR

Every dot is one mouse; bars are cohort means. The naive and rws adults are
drawn in two greys so the reader sees the size of a difference that carries no
developmental meaning.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe plot_region_ratio_young_vs_adult.py
"""

import csv
import math
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA = r'D:\sep_histology\data'
OUT = os.path.join(DATA, 'comparisons', 'young_P20_vs_adult_nano')
MIN_VOX = 2000

AREAS = ['VISp', 'VISl', 'VISal', 'VISrl', 'VISpm', 'VISam', 'SSp-bfd', 'SSp-ul', 'SSp-ll', 'SSp-m', 'SSp-n', 'SSs',
         'AUDp', 'AUDd', 'MOp', 'MOs', 'RSPd', 'RSPv', 'ACAd', 'ACAv', 'PL', 'ILA', 'ORBl',
         '|', 'VPM', 'VPL', 'LGd', 'LP', 'CP', 'ACB', 'CA1', 'CA3', 'DG', 'GPe', 'PVH', 'ZI']
ROWS = [('nano/auto', 'none', 'nanobody / autofluorescence, both background-subtracted  (log2)'),
        ('nano-bg', 'isocortex', 'background-subtracted nanobody, relative to the mouse\'s own isocortex  (log2)'),
        ('nano-bg', 'subcortex-HPF-STR', 'background-subtracted nanobody, relative to subcortex excluding HPF and STR  (log2)')]
COL = {'young_P20': '#c0392b', 'naive': '#555555', 'rws': '#9a9a9a'}
LABEL = {'young_P20': 'P20 (n = 3)', 'naive': 'adult naive (n = 5)', 'rws': 'adult rws (n = 5)'}


def load(path):
    per_mouse = defaultdict(lambda: defaultdict(lambda: [0.0, 0]))
    bg, cohort_of = {}, {}
    with open(path, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            m = r['mouse']; cohort_of[m] = r['cohort']
            if r['acronym'] == 'bg':
                bg[m] = float(r['mean_raw']); continue
            cell = per_mouse[m][r['acronym']]
            cell[0] += float(r['mean_raw']) * int(r['n_vox']); cell[1] += int(r['n_vox'])
    return per_mouse, bg, cohort_of


def main():
    nano, bg_n, cohort_of = load(os.path.join(OUT, 'region_means_raw_per_mouse.csv'))
    auto, bg_a, _ = load(os.path.join(OUT, 'region_means_raw_per_mouse_auto.csv'))
    # division membership for the references, from the ratio CSV (already has it per acronym)
    divi = {}
    with open(os.path.join(OUT, 'region_ratio_young_vs_adult.csv'), newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            divi[r['acronym']] = r['division']
    not_sub = {'Isocortex', 'HPF', 'STR', 'OLF', 'CTXsp', 'fiber tracts', 'VS', 'CB', ''}
    refs = {'none': None, 'isocortex': lambda a: divi.get(a) == 'Isocortex',
            'subcortex-HPF-STR': lambda a: divi.get(a, '') not in not_sub}

    def val(mode, m, a):
        sv, cv = nano[m].get(a, (0.0, 0))
        if cv < MIN_VOX:
            return None
        v = sv / cv - bg_n[m]
        if mode == 'nano-bg':
            return v if v > 0 else None
        sa, ca = auto[m].get(a, (0.0, 0))
        if ca < MIN_VOX:
            return None
        au = sa / ca - bg_a[m]
        return v / au if (au > 0 and v > 0) else None

    def ref_level(mode, m, pred):
        if pred is None:
            return 1.0
        s = c = 0.0
        for a, (sv, cv) in nano[m].items():
            v = val(mode, m, a)
            if v is not None and pred(a):
                s += v * cv; c += cv
        return s / c if c else float('nan')

    mice = sorted(nano, key=lambda m: (cohort_of[m], m))
    fig, axes = plt.subplots(len(ROWS), 1, figsize=(15, 3.6 * len(ROWS)), sharex=True)
    xs = [i for i, a in enumerate(AREAS) if a != '|']
    for ax, (mode, ref, title) in zip(axes, ROWS):
        rl = {m: ref_level(mode, m, refs[ref]) for m in mice}
        for cohort in ('naive', 'rws', 'young_P20'):
            ms = [m for m in mice if cohort_of[m] == cohort]
            jit = np.linspace(-0.22, 0.22, len(ms)) if cohort != 'young_P20' else np.linspace(-0.12, 0.12, len(ms))
            means = []
            for i, a in enumerate(AREAS):
                if a == '|':
                    means.append(np.nan); continue
                ys = []
                for k, m in enumerate(ms):
                    v = val(mode, m, a)
                    if v is not None and rl[m] > 0:
                        y = math.log2(v / rl[m]); ys.append(y)
                        ax.plot(i + jit[k] + (0.28 if cohort == 'young_P20' else -0.1), y, 'o', ms=4.5,
                                color=COL[cohort], alpha=0.85, mec='none')
                means.append(np.mean(ys) if ys else np.nan)
            xo = 0.28 if cohort == 'young_P20' else -0.1
            ax.plot(np.array(xs) + xo, [means[i] for i in xs], '_', ms=14, mew=2.2, color=COL[cohort], label=LABEL[cohort])
        ax.axhline(0, color='k', lw=0.6)
        sep = AREAS.index('|')
        ax.axvline(sep, color='k', lw=0.6, ls=':')
        ax.text(sep - 0.5, ax.get_ylim()[1], 'cortex', ha='right', va='top', fontsize=9, color='#333')
        ax.text(sep + 0.5, ax.get_ylim()[1], 'subcortex', ha='left', va='top', fontsize=9, color='#333')
        ax.set_title(title, fontsize=10.5, loc='left')
        ax.grid(axis='y', lw=0.3, alpha=0.6)
        ax.set_xlim(-0.8, len(AREAS) - 0.2)
    axes[0].legend(loc='lower left', fontsize=9, frameon=False, ncol=3)
    axes[-1].set_xticks(range(len(AREAS)))
    axes[-1].set_xticklabels(['' if a == '|' else a for a in AREAS], rotation=60, ha='right', fontsize=9)
    fig.suptitle('P20 vs adult, nano channel, per mouse, from the raw registered stacks on each cohort\'s own atlas', fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'region_ratio_young_vs_adult.png'), dpi=110)
    print('wrote', os.path.join(OUT, 'region_ratio_young_vs_adult.png'))


if __name__ == '__main__':
    main()
