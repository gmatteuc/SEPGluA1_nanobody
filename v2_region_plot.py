"""
v2, per-mouse region statistics and dot plot, computed on each brain's OWN
atlas -- no warping anywhere.

This is where the numbers to quote come from. A region-wise comparison only
needs every voxel's label, and the DeMBA annotations were remapped to Allen
parcellation_index when they were built, so VISp in a P16 brain is measured
against the P16 annotation, VISp in a P20 brain against the P20 one, and only
the resulting per-mouse means are compared. Pooling ages is therefore clean
here in a way it is not for a voxelwise map, where the young brains have to be
carried into the CCF first (v2_to_ccf).

Per mouse and structure: the mean of ratio (sig/auto) and of sig over the
tissue voxels, then per mouse
  ratio                       nano per unit autofluorescence, no reference.
                              Not an absolute measure: the young cortex is
                              2.0 log2 below the adult in nano and 1.0 log2
                              below it in auto, so the denominator carries its
                              own age effect (see v2_cohort).
  sig / isocortex mean        share of the cortex
  sig / subcortex-HPF-STR     share of the subcortex without the two
                              structures that dominate the scale
  zref                        range-matched: the same cortex-relative values,
                              minus that brain's median over structures and
                              divided by its own p90-p10 spread. Every brain
                              then has the same level AND the same dynamic
                              range, so the question becomes where a region
                              sits inside its own brain's range.
                              Why it is needed: the pup brain is genuinely
                              flatter, p90-p10 = 0.89 +- 0.25 log2 against
                              1.79 +- 0.24 in adults, and no single-number
                              reference can touch that -- dividing by cortex,
                              by subcortex or by the hippocampus shifts every
                              point equally and only moves where zero sits.
                              What it costs: that compression is defined away,
                              so this reading can show re-ordering but says
                              nothing about amplitude.
and per structure a Welch test of young against the ten adults, with the
naive-vs-rws difference printed beside it as the size of a difference that
carries no developmental meaning. The young group is P20 + P16 pooled and the
P16 brain is one young mouse like any other -- same marker, counted in the
median and in the tests. The P20-only contrast stays in the CSV, so what the
P16 brain does to the answer can still be checked.

The figure marks each region with the Mann-Whitney (rank-sum) test of the six
young against the ten adults, uncorrected: * p<0.05, ** p<0.01. Ranks rather
than means because six against ten is small and log ratios are not guaranteed
normal. region_stats.csv carries the Welch p, the Mann-Whitney p and the
Benjamini-Hochberg q of each, so a corrected reading is one column away.

Writes region_means_per_mouse.csv, region_stats.csv and region_plot.png into
data/comparisons_v2/young_vs_adult/.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_region_plot.py
"""

import csv
import math
import os
from collections import defaultdict

import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from v2_per_mouse import annotation_20, MICE, CSV_MAP, OUT as PER_MOUSE, DATA
from v2_cohort import RATIO_CLIP, YOUNG_P20, YOUNG_P16, YOUNG_P22, NAIVE, RWS

OUT = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult')
MIN_VOX = 250      # 20 um voxels = 2 nl, the same volume as the earlier tables
AREAS = ['VISp', 'VISl', 'VISal', 'VISrl', 'VISpm', 'VISam', 'SSp-bfd', 'SSp-ul', 'SSp-ll', 'SSp-m', 'SSp-n', 'SSs',
         'AUDp', 'AUDd', 'MOp', 'MOs', 'RSPd', 'RSPv', 'ACAd', 'ACAv', 'PL', 'ILA', 'ORBl',
         '|', 'VPM', 'VPL', 'LGd', 'LP', 'CP', 'ACB', 'CA1', 'CA3', 'DG', 'GPe', 'PVH', 'ZI']
NOT_SUBCORTEX = {'Isocortex', 'HPF', 'STR', 'OLF', 'CTXsp', 'fiber tracts', 'VS', 'CB', ''}
READINGS = [('ratio', 'nanobody / autofluorescence, both background-subtracted  (log2)'),
            ('cref', 'background-subtracted nanobody, relative to the mouse\'s own isocortex  (log2)'),
            ('subref', 'background-subtracted nanobody, relative to subcortex excluding HPF and STR  (log2)'),
            ('zref', "range-matched: cortex-relative, then centred and scaled by each brain's own spread")]
GROUPS = {'young': YOUNG_P20 + YOUNG_P16 + YOUNG_P22, 'naive': NAIVE, 'rws': RWS}
ADULTS = NAIVE + RWS
COL = {'young': '#c0392b', 'naive': '#555555', 'rws': '#9a9a9a'}
LABEL = {'young': f'young P16-P22 (n = {len(YOUNG_P20) + len(YOUNG_P16) + len(YOUNG_P22)})',
         'naive': f'adult naive (n = {len(NAIVE)})', 'rws': f'adult rws (n = {len(RWS)})'}


def save_figure(fig, path):
    """Save, and if the file is open in a viewer say so instead of dying.

    Windows refuses to overwrite a PNG that an image viewer holds open, and a
    run that writes several figures should not lose the rest because one of
    them was being looked at.
    """
    try:
        fig.savefig(path, dpi=105)
    except OSError:
        alt = path.replace('.png', '_new.png')
        fig.savefig(alt, dpi=105)
        print(f'  NOTE: {os.path.basename(path)} is open elsewhere; wrote {os.path.basename(alt)} instead', flush=True)


def bh_fdr(p):
    """Benjamini-Hochberg q-values for one family of tests.

    A few hundred structures are tested per reading, so a handful of p < 0.05
    is expected from noise alone. The q-value is what should be quoted for
    anything other than the regions named in advance (SS, VIS, prefrontal).
    """
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if ok.sum() == 0:
        return q
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    n = len(ranked)
    adj = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def mannwhitney(a, b):
    """Two-sided rank-sum p for two independent samples, or NaN if too small."""
    if len(a) < 2 or len(b) < 2:
        return float('nan')
    from scipy.stats import mannwhitneyu
    return float(mannwhitneyu(a, b, alternative='two-sided').pvalue)


def welch(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return float('nan')
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    if va + vb == 0:
        return float('nan')
    t = (a.mean() - b.mean()) / math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    from scipy.stats import t as tdist
    return float(2 * tdist.sf(abs(t), df))


def main():
    names, acro, divi = {}, {}, {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            idx = int(row['parcellation_index'])
            if row['parcellation_term_set_name'] == 'structure':
                names[idx] = row['parcellation_term_name']; acro[idx] = row['parcellation_term_acronym']
            elif row['parcellation_term_set_name'] == 'division':
                divi[idx] = row['parcellation_term_acronym']

    anns, per = {}, {}
    mice = [m for g in GROUPS.values() for m in g]
    for mouse in mice:
        cohort, atlas_key = MICE[mouse][:2]
        if atlas_key not in anns:
            anns[atlas_key] = annotation_20(atlas_key)
        ann = anns[atlas_key]
        z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
        sig = z['sig'].astype(np.float32); auto = z['auto'].astype(np.float32); tissue = z['tissue']
        auto_s = gaussian_filter(np.where(tissue, auto, 0), 1.0) / np.maximum(gaussian_filter(tissue.astype(np.float32), 1.0), 1e-3)
        ratio = np.clip(np.where(auto_s > 0, sig / np.maximum(auto_s, 1e-3), 0), -RATIO_CLIP, RATIO_CLIP)
        lab = ann[tissue]; nlab = int(ann.max()) + 1
        n = np.bincount(lab, minlength=nlab)
        s_sig = np.bincount(lab, weights=sig[tissue], minlength=nlab)
        s_rat = np.bincount(lab, weights=ratio[tissue], minlength=nlab)
        d = defaultdict(lambda: [0, 0.0, 0.0])
        for idx in np.nonzero(n)[0]:
            if idx == 0:
                continue
            key = names.get(int(idx), f'id{idx}')
            d[key][0] += int(n[idx]); d[key][1] += s_sig[idx]; d[key][2] += s_rat[idx]
        per[mouse] = {k: (v[0], v[1] / v[0], v[2] / v[0]) for k, v in d.items() if v[0] >= MIN_VOX}
        print(f'{mouse:20s} {atlas_key:10s} {len(per[mouse])} structures', flush=True)

    meta = {nm: (acro[idx], divi.get(idx, '')) for idx, nm in names.items()}
    by_acro = {meta[k][0]: k for k in meta}
    group_of = {m: g for g, ms in GROUPS.items() for m in ms}

    def ref(m, pred):
        s = c = 0.0
        for k, (n, ms, _) in per[m].items():
            if pred(meta.get(k, ('', ''))[1]):
                s += ms * n; c += n
        return s / c if c else float('nan')
    # Each reference is a single number per brain, so every reading is a pure
    # scale and region ratios inside a brain survive it exactly. What changes
    # between them is only the question being asked -- see the header.
    refs = {m: {'cref': ref(m, lambda d: d == 'Isocortex'),
                'subref': ref(m, lambda d: d not in NOT_SUBCORTEX)} for m in mice}

    # The range match needs two numbers per brain rather than one, and they have
    # to come from the same set of structures in every brain or the spread would
    # depend on which regions a section happened to cover.
    common = set.intersection(*[set(per[m]) for m in mice])
    norm = {}
    for m in mice:
        v = np.array([math.log2(per[m][k][1] / refs[m]['cref']) for k in sorted(common)
                      if per[m][k][1] > 0])
        p10, med, p90 = np.percentile(v, [10, 50, 90])
        norm[m] = (med, max(p90 - p10, 1e-6))
    print('dynamic range per brain (p90-p10 of log2 over %d shared structures):' % len(common))
    for m in mice:
        print(f'  {m:20s} median {norm[m][0]:+.2f}   spread {norm[m][1]:.2f}')

    def value(reading, m, k):
        if k is None or k not in per[m]:
            return None
        n, ms, mr = per[m][k]
        if reading == 'zref':
            if ms <= 0:
                return None
            med, spread = norm[m]
            return (math.log2(ms / refs[m]['cref']) - med) / spread
        v = mr if reading == 'ratio' else ms / refs[m][reading]
        return math.log2(v) if v > 0 else None

    # ------------------------------------------------------- tables
    rows_pm, rows_st = [], []
    for k in sorted(meta, key=lambda k: (meta[k][1], meta[k][0])):
        for reading, _ in READINGS:
            v = {m: value(reading, m, k) for m in mice}
            v = {m: x for m, x in v.items() if x is not None}
            for m, x in v.items():
                rows_pm.append((reading, group_of[m], MICE[m][0], m, k, meta[k][0], meta[k][1], per[m][k][0], x))
            yo = [v[m] for m in GROUPS['young'] if m in v]
            y20 = [v[m] for m in YOUNG_P20 if m in v]
            p16 = [v[m] for m in YOUNG_P16 if m in v]
            ad = [v[m] for m in ADULTS if m in v]
            nv = [v[m] for m in NAIVE if m in v]
            rw = [v[m] for m in RWS if m in v]
            if len(yo) < 2 or len(ad) < 4:
                continue
            rows_st.append((reading, k, meta[k][0], meta[k][1], len(yo), len(ad),
                            np.mean(yo), np.mean(ad), np.mean(yo) - np.mean(ad),
                            np.median(yo) - np.median(ad), welch(yo, ad), mannwhitney(yo, ad),
                            (np.mean(y20) - np.mean(ad)) if len(y20) >= 2 else float('nan'),
                            welch(y20, ad) if len(y20) >= 2 else float('nan'),
                            (p16[0] - np.mean(ad)) if p16 else float('nan'),
                            (np.mean(nv) - np.mean(rw)) if nv and rw else float('nan')))
    # q-values within each reading, so the brain-wide lists can be read honestly
    q_welch, q_mw = {}, {}
    for reading, _ in READINGS:
        idx = [i for i, r in enumerate(rows_st) if r[0] == reading]
        for store, col in ((q_welch, 10), (q_mw, 11)):
            for i, qi in zip(idx, bh_fdr([rows_st[i][col] for i in idx])):
                store[i] = qi
    rows_st = [r[:12] + (q_welch[i], q_mw[i]) + r[12:] for i, r in enumerate(rows_st)]

    with open(os.path.join(OUT, 'region_means_per_mouse.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['reading', 'group', 'cohort', 'mouse', 'structure', 'acronym', 'division', 'n_vox20', 'log2_value'])
        w.writerows([r[:8] + (f'{r[8]:.4f}',) for r in rows_pm])
    with open(os.path.join(OUT, 'region_stats.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['reading', 'structure', 'acronym', 'division', 'n_young', 'n_adult',
                    'young_mean_log2', 'adult_mean_log2', 'diff_log2', 'diff_median_log2',
                    'welch_p', 'mannwhitney_p', 'welch_q_BH', 'mannwhitney_q_BH',
                    'diff_log2_P20only', 'welch_p_P20only', 'diff_log2_P16_single', 'naive_minus_rws_log2'])
        w.writerows([r[:6] + tuple(f'{x:.4f}' for x in r[6:]) for r in rows_st])

    st = {(r[0], r[2]): r for r in rows_st}
    print(f'\nCORTEX  log2(young / adult), young = {len(GROUPS["young"])} mice (P20 + P16) vs {len(ADULTS)} adults '
          '(* p<0.05, ** p<0.01, Mann-Whitney, uncorrected; q in the CSV). '
          'P20only = without the P16 brain; P16 = that brain alone; naive-rws = the null scale.')
    print(f'  {"area":9s} ' + ' '.join(f'{r:>10s}' for r, _ in READINGS) + f' {"P20only":>9s} {"P16":>7s} {"naive-rws":>10s}')
    for a in AREAS:
        if a == '|':
            print('  ' + '-' * 60); continue
        cells = []
        for reading, _ in READINGS:
            r = st.get((reading, a))
            if r is None:
                cells.append(f'{"--":>10s}'); continue
            star = '**' if r[11] < 0.01 else ('*' if r[11] < 0.05 else '')      # rank-sum p
            cells.append(f'{r[8]:+7.2f}{star:3s}')
        r = st.get(('cref', a))
        tail = f'{r[14]:+9.2f} {r[16]:+7.2f} {r[17]:+10.2f}' if r else ''
        print(f'  {a:9s} ' + ' '.join(cells) + ' ' + tail)

    # ------------------------------------------------------- figure
    # Every mouse is a dot of the same size, the P16 brain included: it is one
    # young animal among six. The bar is the group MEDIAN, to match the rank-sum
    # test that puts the stars on.
    star_of = {(r[0], r[2]): ('**' if r[11] < 0.01 else ('*' if r[11] < 0.05 else '')) for r in rows_st}
    ylab = {'ratio': 'log2  nano / auto', 'cref': 'log2  relative to own isocortex',
            'subref': 'log2  relative to subcortex', 'zref': 'range-matched (median 0, spread 1)'}
    fig, axes = plt.subplots(len(READINGS), 1, figsize=(15, 3.8 * len(READINGS)), sharex=True)
    xs = [i for i, a in enumerate(AREAS) if a != '|']
    for ax, (reading, title) in zip(axes, READINGS):
        for g in ('naive', 'rws', 'young'):
            ms = GROUPS[g]
            jit = np.linspace(-0.22, 0.22, len(ms))
            xo = 0.28 if g == 'young' else -0.1
            mids = []
            for i, a in enumerate(AREAS):
                if a == '|':
                    mids.append(np.nan); continue
                k = by_acro.get(a); ys = []
                for j, m in enumerate(ms):
                    y = value(reading, m, k)
                    if y is None:
                        continue
                    ys.append(y)
                    ax.plot(i + jit[j] + xo, y, 'o', ms=4.5, color=COL[g], alpha=0.9, mec='none')
                mids.append(np.median(ys) if ys else np.nan)
            ax.plot(np.array(xs) + xo, [mids[i] for i in xs], '_', ms=14, mew=2.2, color=COL[g], label=LABEL[g])
        # stars for the young-vs-adult rank-sum test, just under the top of the panel
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.12 * (hi - lo))
        lo, hi = ax.get_ylim()
        for i, a in enumerate(AREAS):
            st = star_of.get((reading, a), '')
            if st:
                ax.text(i, hi - 0.04 * (hi - lo), st, ha='center', va='top', fontsize=11, color=COL['young'])
        ax.axhline(0, color='k', lw=0.6)
        sep = AREAS.index('|'); ax.axvline(sep, color='k', lw=0.6, ls=':')
        box = dict(facecolor='w', edgecolor='none', alpha=0.85, pad=1.5)
        ax.text(sep - 0.5, lo + 0.02 * (hi - lo), 'cortex', ha='right', va='bottom', fontsize=9, color='#333', bbox=box)
        ax.text(sep + 0.5, lo + 0.02 * (hi - lo), 'subcortex', ha='left', va='bottom', fontsize=9, color='#333', bbox=box)
        ax.set_title(title, fontsize=10.5, loc='left')
        ax.set_ylabel(ylab[reading], fontsize=10)
        ax.grid(axis='y', lw=0.3, alpha=0.6)
        ax.set_xlim(-0.8, len(AREAS) - 0.2)
    axes[0].legend(loc='lower left', fontsize=9, frameon=True, framealpha=0.9, edgecolor='none', ncol=3)
    axes[-1].set_xticks(range(len(AREAS)))
    axes[-1].set_xticklabels(['' if a == '|' else a for a in AREAS], rotation=60, ha='right', fontsize=9)
    fig.suptitle('Young vs adult, nano channel: one dot per mouse, each brain measured on the atlas of its own age\n'
                 f'Bars are group medians.  * p<0.05, ** p<0.01, Mann-Whitney '
                 f'{len(GROUPS["young"])} vs {len(ADULTS)}, uncorrected',
                 fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    save_figure(fig, os.path.join(OUT, 'region_plot.png'))
    print('\nwrote', os.path.join(OUT, 'region_plot.png'))


if __name__ == '__main__':
    main()
