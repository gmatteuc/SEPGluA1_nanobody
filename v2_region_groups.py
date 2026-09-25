"""
The same per-mouse measurements as v2_region_plot, aggregated two ways that
ask coarser questions than "one area at a time":

  by system    primary and higher-order visual, somatosensory and auditory,
               then frontal, motor, retrosplenial, lateral/insular and the
               subcortical divisions -- so a shift spread thinly over many
               small areas is visible instead of being split into
               non-significant pieces, and so the primary/higher-order
               distinction the critical-period argument rests on is explicit.
  by layer     within each cortical system, supragranular (1, 2/3), granular
               (4) and infragranular (5, 6a, 6b), read from the ontology's
               substructure names. This is the laminar pattern Sami noticed
               by eye in the videos, measured.

Everything is computed on each brain's OWN atlas, no warping: a group mean is
the voxel-weighted mean over its member labels in that brain. The three
readings and the references are the ones v2_region_plot uses (ratio = nano per
autofluorescence; cref = relative to that brain's isocortex; subref = relative
to the subcortex excluding HPF and STR; sepratio = nano per unit SEP, i.e.
surface receptor per unit receptor expressed; zref = range-matched to each brain's
own spread), and the test is the same Mann-Whitney
of the young group against the adults, uncorrected in the figure, with BH
q-values in the CSV.

Writes into data/comparisons_v2/young_vs_adult/:
  group_stats.csv        one row per (reading, grouping, group), with medians,
                         both tests, the P20-only contrast and the naive-rws null
  group_plot.png         systems, one dot per mouse
  laminar_plot.png       layers within each cortical system

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_region_groups.py
"""

import csv
import math
import os
import re
from collections import defaultdict

import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from v2_per_mouse import annotation_20, MICE, CSV_MAP, DATA, OUT as PER_MOUSE
from v2_cohort import RATIO_CLIP, YOUNG_P20, YOUNG_P16, YOUNG_P22, NAIVE, RWS, SIGNED_READINGS
from v2_region_plot import bh_fdr, mannwhitney, welch, READINGS, NOT_SUBCORTEX, COL

OUT = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult')
GROUPS = {'young': YOUNG_P20 + YOUNG_P16 + YOUNG_P22, 'naive': NAIVE, 'rws': RWS}
ADULTS = NAIVE + RWS
LABEL = {'young': f'young P16-P22 (n = {len(YOUNG_P20) + len(YOUNG_P16) + len(YOUNG_P22)})',
         'naive': f'adult naive (n = {len(NAIVE)})', 'rws': f'adult rws (n = {len(RWS)})'}

# Cortical systems, split primary vs higher order, which is the distinction the
# critical-period argument rests on: a thalamorecipient primary area and its
# higher-order neighbours mature on different schedules.
SYSTEMS = [
    ('primary visual (VISp)', lambda a: a == 'VISp'),
    # RL and AL stand apart from the rest of the higher visual belt: they sit on
    # the visual-somatosensory border and are the ones that look modulated in
    # the maps, so Sami asked for them as a group of their own rather than
    # averaged into the belt that surrounds them.
    ('associative VT (RL+AL)', lambda a: a in ('VISrl', 'VISal')),
    ('higher visual (all other)', lambda a: a.startswith('VIS') and a not in ('VISp', 'VISC', 'VISrl', 'VISal')),
    ('primary somatosensory (SSp)', lambda a: a.startswith('SSp')),
    ('higher somatosensory (SSs)', lambda a: a == 'SSs'),
    ('primary auditory (AUDp)', lambda a: a == 'AUDp'),
    ('higher auditory (AUDd/po/v)', lambda a: a.startswith('AUD') and a != 'AUDp'),
    ('frontal', lambda a: a.startswith(('ACA', 'ORB')) or a in ('PL', 'ILA', 'FRP', 'DP')),
    ('motor', lambda a: a in ('MOp', 'MOs')),
    ('retrosplenial', lambda a: a.startswith('RSP')),
    ('lateral/insular', lambda a: a.startswith('AI') or a in ('GU', 'VISC', 'TEa', 'PERI', 'ECT')),
]
# the laminar figure would be unreadable with every system on it; these are the
# ones the layer question is actually about
# Naming: every sensory system says in brackets what it contains, so a reader
# never has to guess whether RL and AL are also inside 'higher visual'.
LAMINAR_SYSTEMS = ('primary visual (VISp)', 'associative VT (RL+AL)', 'higher visual (all other)',
                   'primary somatosensory (SSp)', 'higher somatosensory (SSs)',
                   'primary auditory (AUDp)', 'higher auditory (AUDd/po/v)', 'frontal', 'retrosplenial')
# and the subcortical side, taken straight from the ontology's division
DIVISIONS = [('thalamus', 'TH'), ('striatum', 'STR'), ('pallidum', 'PAL'), ('hippocampus', 'HPF'),
             ('hypothalamus', 'HY'), ('midbrain', 'MB'), ('olfactory', 'OLF'), ('cortical subplate', 'CTXsp')]
LAYERS = [('supragranular', ('1', '2', '3', '2/3')), ('granular', ('4',)), ('infragranular', ('5', '6', '6a', '6b'))]


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


def layer_of(substructure_name):
    """'Primary visual area, layer 2/3' -> '2/3'; anything without a layer -> None."""
    m = re.search(r'layer\s*([0-9]+(?:/[0-9]+)?[ab]?)', substructure_name, re.I)
    return m.group(1) if m else None


def main():
    # index -> (structure acronym, division acronym, layer or None)
    stru, divi, sub = {}, {}, {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            i = int(r['parcellation_index'])
            if r['parcellation_term_set_name'] == 'structure':
                stru[i] = r['parcellation_term_acronym']
            elif r['parcellation_term_set_name'] == 'division':
                divi[i] = r['parcellation_term_acronym']
            elif r['parcellation_term_set_name'] == 'substructure':
                sub[i] = r['parcellation_term_name']
    layer = {i: layer_of(nm) for i, nm in sub.items()}

    # every group is a set of parcellation indices
    groups = {}
    for name, pred in SYSTEMS:
        groups[('system', name)] = {i for i, a in stru.items() if divi.get(i) == 'Isocortex' and pred(a)}
    for name, div in DIVISIONS:
        groups[('system', name)] = {i for i in stru if divi.get(i) == div}
    for name, pred in SYSTEMS:
        if name not in LAMINAR_SYSTEMS:
            continue
        for lname, tokens in LAYERS:
            groups[('layer', f'{name} {lname}')] = {i for i, a in stru.items()
                                                    if divi.get(i) == 'Isocortex' and pred(a) and layer.get(i) in tokens}
    groups = {k: v for k, v in groups.items() if v}

    # per mouse: the voxel-weighted mean of sig and of ratio in every group,
    # plus the two references, all on that brain's own atlas
    anns, per, refs, struct_mean = {}, {}, {}, {}
    mice = [m for g in GROUPS.values() for m in g]
    for mouse in mice:
        atlas_key = MICE[mouse][1]
        if atlas_key not in anns:
            anns[atlas_key] = annotation_20(atlas_key)
        ann = anns[atlas_key]
        z = np.load(os.path.join(PER_MOUSE, mouse + '.npz'))
        sig = z['sig'].astype(np.float32); auto = z['auto'].astype(np.float32); tissue = z['tissue']
        if any(r == 'sepratio' for r, _ in READINGS) and 'sep' not in z.files:
            raise SystemExit(f'{mouse}: no SEP channel in its per-mouse file. Run\n'
                             f'  P4bis_add_sep_channel.m for this brain, then v2_per_mouse.py,\n'
                             f'  or drop the reading with V2_READINGS.')

        # sig per unit of a reference CHANNEL, the denominator smoothed by one
        # 20 um voxel and mask-normalised, as in v2_cohort and v2_region_plot
        def per_unit(ref):
            ref_s = gaussian_filter(np.where(tissue, ref, 0), 1.0) / np.maximum(
                gaussian_filter(tissue.astype(np.float32), 1.0), 1e-3)
            return np.clip(np.where(ref_s > 0, sig / np.maximum(ref_s, 1e-3), 0), -RATIO_CLIP, RATIO_CLIP)

        ratio = per_unit(auto)
        sepratio = per_unit(z['sep'].astype(np.float32)) if 'sep' in z.files else np.zeros_like(sig)
        lab = ann[tissue]; nlab = int(ann.max()) + 1
        n = np.bincount(lab, minlength=nlab)
        s_sig = np.bincount(lab, weights=sig[tissue], minlength=nlab)
        s_rat = np.bincount(lab, weights=ratio[tissue], minlength=nlab)
        s_sep = np.bincount(lab, weights=sepratio[tissue], minlength=nlab)
        per[mouse] = {}
        for key, ids in groups.items():
            ids = [i for i in ids if i < nlab]
            c = n[ids].sum()
            per[mouse][key] = (int(c), s_sig[ids].sum() / c, s_rat[ids].sum() / c,
                               s_sep[ids].sum() / c) if c >= 250 else None
        iso = [i for i in stru if divi.get(i) == 'Isocortex' and i < nlab]
        sub_ids = [i for i in stru if divi.get(i, '') not in NOT_SUBCORTEX and i < nlab]
        refs[mouse] = {'cref': s_sig[iso].sum() / n[iso].sum(),
                       'subref': s_sig[sub_ids].sum() / n[sub_ids].sum()}
        # the brain's own distribution over structures, for the range match --
        # computed at structure level so it does not depend on the grouping
        by_struct = defaultdict(lambda: [0, 0.0])
        for i in np.nonzero(n)[0]:
            if i == 0 or i not in stru:
                continue
            by_struct[stru[i]][0] += int(n[i]); by_struct[stru[i]][1] += s_sig[i]
        struct_mean[mouse] = {k: v[1] / v[0] for k, v in by_struct.items() if v[0] >= 250}
        print(f'{mouse:20s} {sum(v is not None for v in per[mouse].values())}/{len(groups)} groups', flush=True)

    common = set.intersection(*[set(struct_mean[m]) for m in mice])
    norm = {}
    for m in mice:
        v = np.array([math.log2(struct_mean[m][k] / refs[m]['cref']) for k in sorted(common)
                      if struct_mean[m][k] > 0])
        p10, med, p90 = np.percentile(v, [10, 50, 90])
        norm[m] = (med, max(p90 - p10, 1e-6))

    def value(reading, mouse, key):
        cell = per[mouse].get(key)
        if cell is None:
            return None
        _, m_sig, m_rat, m_sep = cell
        if reading in SIGNED_READINGS:
            if m_sig <= 0:
                return None
            med, spread = norm[mouse]
            return (math.log2(m_sig / refs[mouse]['cref']) - med) / spread
        if reading in ('ratio', 'sepratio'):
            v = m_rat if reading == 'ratio' else m_sep
        else:
            v = m_sig / refs[mouse][reading]
        return math.log2(v) if v > 0 else None

    # ------------------------------------------------------------- stats
    rows = []
    for key in groups:
        for reading, _ in READINGS:
            v = {m: value(reading, m, key) for m in mice}
            v = {m: x for m, x in v.items() if x is not None}
            yo = [v[m] for m in GROUPS['young'] if m in v]
            y20 = [v[m] for m in YOUNG_P20 if m in v]
            ad = [v[m] for m in ADULTS if m in v]
            nv = [v[m] for m in NAIVE if m in v]
            rw = [v[m] for m in RWS if m in v]
            if len(yo) < 3 or len(ad) < 5:
                continue
            rows.append(dict(reading=reading, grouping=key[0], group=key[1], n_young=len(yo), n_adult=len(ad),
                             young_median=np.median(yo), adult_median=np.median(ad),
                             diff_median=np.median(yo) - np.median(ad), diff_mean=np.mean(yo) - np.mean(ad),
                             mannwhitney_p=mannwhitney(yo, ad), welch_p=welch(yo, ad),
                             diff_P20only=(np.median(y20) - np.median(ad)) if len(y20) >= 3 else float('nan'),
                             naive_minus_rws=(np.median(nv) - np.median(rw)) if nv and rw else float('nan')))
    for reading, _ in READINGS:
        for grouping in ('system', 'layer'):
            idx = [i for i, r in enumerate(rows) if r['reading'] == reading and r['grouping'] == grouping]
            for i, q in zip(idx, bh_fdr([rows[i]['mannwhitney_p'] for i in idx])):
                rows[i]['mannwhitney_q_BH'] = q
    with open(os.path.join(OUT, 'group_stats.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow({k: (f'{x:.4f}' if isinstance(x, float) else x) for k, x in r.items()})
    star = {(r['reading'], r['grouping'], r['group']):
            ('**' if r['mannwhitney_p'] < 0.01 else ('*' if r['mannwhitney_p'] < 0.05 else '')) for r in rows}

    print(f'\n{"group":22s} ' + ' '.join(f'{r:>12s}' for r, _ in READINGS) + '   (log2 young - adult, medians)')
    for grouping in ('system', 'layer'):
        print(f'--- by {grouping}')
        for key in [k for k in groups if k[0] == grouping]:
            cells = []
            for reading, _ in READINGS:
                r = next((r for r in rows if r['reading'] == reading and r['group'] == key[1]
                          and r['grouping'] == grouping), None)
                cells.append(f'{r["diff_median"]:+8.2f}{star[(reading, grouping, key[1])]:<3s}' if r else f'{"--":>11s}')
            print(f'  {key[1]:22s} ' + ' '.join(cells))

    # ------------------------------------------------------------ figures
    def dotplot(keys, labels, fname, title):
        fig, axes = plt.subplots(len(READINGS), 1, figsize=(max(9, 0.62 * len(keys) + 4), 3.6 * len(READINGS)), sharex=True)
        for ax, (reading, rtitle) in zip(axes, READINGS):
            for g in ('naive', 'rws', 'young'):
                ms = GROUPS[g]
                jit = np.linspace(-0.2, 0.2, len(ms))
                xo = 0.26 if g == 'young' else -0.1
                mids = []
                for i, key in enumerate(keys):
                    ys = [value(reading, m, key) for m in ms]
                    ys = [y for y in ys if y is not None]
                    for j, y in enumerate(ys):
                        ax.plot(i + jit[j] + xo, y, 'o', ms=4.5, color=COL[g], alpha=0.9, mec='none')
                    mids.append(np.median(ys) if ys else np.nan)
                ax.plot(np.arange(len(keys)) + xo, mids, '_', ms=16, mew=2.4, color=COL[g], label=LABEL[g])
            lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.14 * (hi - lo)); lo, hi = ax.get_ylim()
            for i, key in enumerate(keys):
                st = star.get((reading, key[0], key[1]), '')
                if st:
                    ax.text(i, hi - 0.04 * (hi - lo), st, ha='center', va='top', fontsize=12, color=COL['young'])
            ax.axhline(0, color='k', lw=0.6)
            ax.set_title(rtitle, fontsize=10.5, loc='left')
            ax.set_ylabel({'ratio': 'log2  nano / auto', 'sepratio': 'log2  nano / SEP',
                           'cref': 'log2  vs own isocortex', 'subref': 'log2  vs subcortex',
                           'zref': 'range-matched'}[reading], fontsize=10)
            ax.grid(axis='y', lw=0.3, alpha=0.6); ax.set_xlim(-0.7, len(keys) - 0.3)
        axes[0].legend(loc='lower left', fontsize=9, frameon=True, framealpha=0.9, edgecolor='none', ncol=3)
        axes[-1].set_xticks(range(len(keys)))
        axes[-1].set_xticklabels(labels, rotation=55, ha='right', fontsize=9)
        fig.suptitle(title, fontsize=11.5)
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        save_figure(fig, os.path.join(OUT, fname)); plt.close(fig)

    sys_keys = [k for k in groups if k[0] == 'system']
    dotplot(sys_keys, [k[1] for k in sys_keys], 'group_plot.png',
            'Young vs adult by system, one dot per mouse, each brain on the atlas of its own age\n'
            f'Bars are group medians.  * p<0.05, ** p<0.01 '
            f'(Mann-Whitney {len(GROUPS["young"])} vs {len(ADULTS)}, uncorrected)')
    lay_keys = [k for k in groups if k[0] == 'layer']
    dotplot(lay_keys, [k[1].replace(' supragranular', ' L1-3').replace(' granular', ' L4').replace(' infragranular', ' L5-6')
                       for k in lay_keys], 'laminar_plot.png',
            'Young vs adult by cortical layer within each system (layers from the ontology)\n'
            f'Bars are group medians.  * p<0.05, ** p<0.01 '
            f'(Mann-Whitney {len(GROUPS["young"])} vs {len(ADULTS)}, uncorrected)')
    print('\nwrote group_plot.png, laminar_plot.png, group_stats.csv')


if __name__ == '__main__':
    main()
