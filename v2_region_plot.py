"""
v2, per-mouse region plot: the same figure as plot_region_ratio_young_vs_adult
but computed from the v2 per-mouse volumes, so it rests on exactly the tissue
mask, background subtraction and readings that the v2 maps use.

Per mouse and structure: mean of ratio (sig/auto) and of sig over the tissue
voxels, on the cohort's own 20 um atlas labels. Then, per mouse, log2 of
  ratio                       nano per unit autofluorescence, no reference
  sig / isocortex mean        share of the cortex
  sig / subcortex-HPF-STR     share of the subcortex without the two
                              structures that dominate the scale

Writes region_means_per_mouse.csv (both readings, every structure with at
least MIN_VOX tissue voxels), region_stats.csv (cohort means, Welch p vs the
ten adults, naive-vs-rws null) and region_plot.png.

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

from v2_per_mouse import annotation_20, MICE, CSV_MAP
from v2_cohort import PER_MOUSE, V2, RATIO_CLIP

OUT = os.path.join(V2, 'young_P20_vs_adult')
MIN_VOX = 250      # 20 um voxels = 2 nl, the same volume as before
AREAS = ['VISp', 'VISl', 'VISal', 'VISrl', 'VISpm', 'VISam', 'SSp-bfd', 'SSp-ul', 'SSp-ll', 'SSp-m', 'SSp-n', 'SSs',
         'AUDp', 'AUDd', 'MOp', 'MOs', 'RSPd', 'RSPv', 'ACAd', 'ACAv', 'PL', 'ILA', 'ORBl',
         '|', 'VPM', 'VPL', 'LGd', 'LP', 'CP', 'ACB', 'CA1', 'CA3', 'DG', 'GPe', 'PVH', 'ZI']
NOT_SUBCORTEX = {'Isocortex', 'HPF', 'STR', 'OLF', 'CTXsp', 'fiber tracts', 'VS', 'CB', ''}
READINGS = [('ratio', 'nanobody / autofluorescence, both background-subtracted  (log2)'),
            ('cref', 'background-subtracted nanobody, relative to the mouse\'s own isocortex  (log2)'),
            ('subref', 'background-subtracted nanobody, relative to subcortex excluding HPF and STR  (log2)')]
COL = {'young_P20': '#c0392b', 'naive': '#555555', 'rws': '#9a9a9a'}
LABEL = {'young_P20': 'P20 (n = 3)', 'naive': 'adult naive (n = 5)', 'rws': 'adult rws (n = 5)'}


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

    # per mouse, per structure NAME: voxel-weighted mean of sig and of ratio
    anns = {}
    per = {}                                   # per[mouse][name] = (n, mean_sig, mean_ratio)
    for mouse, (cohort, atlas_key, *_r) in MICE.items():
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
        print(f'{mouse:20s} {len(per[mouse])} structures', flush=True)
    meta = {}
    for idx, nm in names.items():
        meta[nm] = (acro[idx], divi.get(idx, ''))
    mice = sorted(per, key=lambda m: (MICE[m][0], m))
    cohort_of = {m: MICE[m][0] for m in mice}
    by_acro = {meta[k][0]: k for k in meta}

    def ref(m, pred):
        s = c = 0.0
        for k, (n, ms, _) in per[m].items():
            if pred(meta.get(k, ('', ''))[1]):
                s += ms * n; c += n
        return s / c if c else float('nan')
    refs = {m: {'cref': ref(m, lambda d: d == 'Isocortex'), 'subref': ref(m, lambda d: d not in NOT_SUBCORTEX)} for m in mice}

    def value(reading, m, k):
        if k not in per[m]:
            return None
        n, ms, mr = per[m][k]
        v = mr if reading == 'ratio' else ms / refs[m][reading]
        return math.log2(v) if v > 0 else None

    # ------------------------------------------------------- tables
    rows_pm, rows_st = [], []
    for k in sorted(meta, key=lambda k: (meta[k][1], meta[k][0])):
        for reading, _ in READINGS:
            v = {m: value(reading, m, k) for m in mice}
            v = {m: x for m, x in v.items() if x is not None}
            for m, x in v.items():
                rows_pm.append((reading, cohort_of[m], m, k, meta[k][0], meta[k][1], per[m][k][0], x))
            yp = [v[m] for m in mice if cohort_of[m] == 'young_P20' and m in v]
            ad = [v[m] for m in mice if cohort_of[m] != 'young_P20' and m in v]
            nv = [v[m] for m in mice if cohort_of[m] == 'naive' and m in v]
            rw = [v[m] for m in mice if cohort_of[m] == 'rws' and m in v]
            if len(yp) >= 2 and len(ad) >= 4:
                rows_st.append((reading, k, meta[k][0], meta[k][1], len(yp), len(ad), np.mean(yp), np.mean(ad),
                                np.mean(yp) - np.mean(ad), welch(yp, ad),
                                (np.mean(nv) - np.mean(rw)) if nv and rw else float('nan')))
    with open(os.path.join(OUT, 'region_means_per_mouse.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh); w.writerow(['reading', 'cohort', 'mouse', 'structure', 'acronym', 'division', 'n_vox20', 'log2_value'])
        w.writerows([r[:7] + (f'{r[7]:.4f}',) for r in rows_pm])
    with open(os.path.join(OUT, 'region_stats.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['reading', 'structure', 'acronym', 'division', 'n_P20', 'n_adult', 'P20_mean_log2', 'adult_mean_log2',
                    'diff_log2', 'welch_p', 'naive_minus_rws_log2'])
        w.writerows([r[:6] + tuple(f'{x:.4f}' for x in r[6:]) for r in rows_st])
    st = {(r[0], r[2]): r for r in rows_st}
    print('\nCORTEX  log2(P20/adult) per reading (* p<0.05, ** p<0.01), last column naive-rws under cref:')
    print(f'  {"area":9s} ' + ' '.join(f'{r:>10s}' for r, _ in READINGS) + f' {"naive-rws":>10s}')
    for a in AREAS:
        if a == '|':
            print('  ' + '-' * 50); continue
        k = by_acro.get(a)
        cells = []
        for reading, _ in READINGS:
            r = st.get((reading, a))
            if r is None:
                cells.append(f'{"--":>10s}'); continue
            star = '**' if r[9] < 0.01 else ('*' if r[9] < 0.05 else '')
            cells.append(f'{r[8]:+7.2f}{star:3s}')
        r = st.get(('cref', a))
        print(f'  {a:9s} ' + ' '.join(cells) + (f' {r[10]:+10.2f}' if r else ''))

    # ------------------------------------------------------- figure
    fig, axes = plt.subplots(len(READINGS), 1, figsize=(15, 3.6 * len(READINGS)), sharex=True)
    xs = [i for i, a in enumerate(AREAS) if a != '|']
    for ax, (reading, title) in zip(axes, READINGS):
        for cohort in ('naive', 'rws', 'young_P20'):
            ms = [m for m in mice if cohort_of[m] == cohort]
            jit = np.linspace(-0.22, 0.22, len(ms)) if cohort != 'young_P20' else np.linspace(-0.12, 0.12, len(ms))
            xo = 0.28 if cohort == 'young_P20' else -0.1
            means = []
            for i, a in enumerate(AREAS):
                if a == '|':
                    means.append(np.nan); continue
                k = by_acro.get(a); ys = []
                for j, m in enumerate(ms):
                    y = value(reading, m, k) if k else None
                    if y is not None:
                        ys.append(y); ax.plot(i + jit[j] + xo, y, 'o', ms=4.5, color=COL[cohort], alpha=0.85, mec='none')
                means.append(np.mean(ys) if ys else np.nan)
            ax.plot(np.array(xs) + xo, [means[i] for i in xs], '_', ms=14, mew=2.2, color=COL[cohort], label=LABEL[cohort])
        ax.axhline(0, color='k', lw=0.6)
        sep = AREAS.index('|'); ax.axvline(sep, color='k', lw=0.6, ls=':')
        ax.text(sep - 0.5, ax.get_ylim()[1], 'cortex', ha='right', va='top', fontsize=9, color='#333')
        ax.text(sep + 0.5, ax.get_ylim()[1], 'subcortex', ha='left', va='top', fontsize=9, color='#333')
        ax.set_title(title, fontsize=10.5, loc='left'); ax.grid(axis='y', lw=0.3, alpha=0.6)
        ax.set_xlim(-0.8, len(AREAS) - 0.2)
    axes[0].legend(loc='lower left', fontsize=9, frameon=False, ncol=3)
    axes[-1].set_xticks(range(len(AREAS)))
    axes[-1].set_xticklabels(['' if a == '|' else a for a in AREAS], rotation=60, ha='right', fontsize=9)
    fig.suptitle('P20 vs adult, nano channel, per mouse -- v2 tissue masks (auto channel) and background subtraction, '
                 'each cohort on its own atlas', fontsize=11.5)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'region_plot.png'), dpi=110)
    print('wrote', os.path.join(OUT, 'region_plot.png'))


if __name__ == '__main__':
    main()
