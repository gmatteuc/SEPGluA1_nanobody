"""
v2, cortical flatmaps: the whole isocortex on one picture, and by layer.

A coronal plane shows a sliver of each area and nothing of its neighbours along
the anterior-posterior axis, which is why "RL and AL specifically" is hard to
see in a slice and easy to see here. The Allen streamline flatmap unrolls the
cortex along the streamlines that run from pia to white matter, so every area
becomes a patch on a 2D map and the depth becomes an axis of its own.

Two figures per reading:

  <reading>_flatmap.png         young | adult | difference, averaged through
                                the full cortical depth
  <reading>_flatmap_layers.png  the same three, separately for supragranular
                                (L1 + L2/3), granular (L4) and infragranular
                                (L5 + L6), with the depth normalised per
                                streamline so a layer lands at the same depth
                                everywhere despite cortex being thicker in some
                                areas than others

What is projected. Cohort volumes come from v2_cohort at 20 um in the adult CCF;
the Allen streamline assets are defined on the 10 um CCF, so each volume is
repeated 2x along each axis, which is exactly the inverse of the block-mean that
produced it -- no interpolation is introduced. Hemispheres are folded first (as
everywhere else in v2) and the fold is mirrored back into both halves, so the
butterfly view shows the hemisphere-averaged data on both wings and coverage is
as complete as the data allow.

Depth averaging is done as a ratio of two projections, sum(value) / sum(mask),
so a voxel with too few brains behind it contributes to neither. Projecting a
volume with NaN in it would poison whole streamlines instead.

Assets (once, ~0.6 GB, see data/atlas_flatmap): surface_paths_10_v3.h5,
flatmap_butterfly.h5, flatmap_butterfly.nrrd, labelDescription_ITKSNAPColor.txt,
avg_layer_depths.json, cortical_layers_10_v2.h5, all from
download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/
cortical_coordinates/ccf_2017/ccf_streamlines_assets/

This one runs in its own environment, tools\\venv_flat, because ccf_streamlines
pulls in its own numpy and scikit-image and the analysis venv is not worth
disturbing:

  D:\\sep_histology\\code\\tools\\venv_flat\\Scripts\\python.exe v2_flatmap.py [reading ...]
"""

import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as path_effects

DATA = r'D:\sep_histology\data'
ASSETS = os.path.join(DATA, 'atlas_flatmap')
CCF_ROOT = os.path.join(DATA, 'comparisons_v2', 'ccf')
OUT = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult')

YOUNG, ADULT = 'young', 'adult'
MIN_N_YOUNG, MIN_N_ADULT = 2, 5          # as in v2_compare, so the maps agree
MEAN_VMAX = {'zref': 1.0, 'cref': 2.0, 'subref': 2.0, 'ratio': 2.0, 'sepratio': 0.6}
DIFF_LIM = {'zref': 0.5, 'cref': 1.0, 'subref': 1.0, 'ratio': 1.0, 'sepratio': 1.0}
# the areas worth naming on a map that already has every border drawn
LABEL_AREAS = ['VISp', 'VISrl', 'VISal', 'VISl', 'VISam', 'SSp-bfd', 'SSp-m', 'SSs',
               'AUDp', 'MOp', 'MOs', 'RSPd', 'ACAd', 'ORBl', 'TEa']
# drawn in white and named in full: the areas the argument is about
HIGHLIGHT = ('VISrl', 'VISal', 'VISp', 'SSp-bfd')
# RL and AL are neighbours and their centroids are 100 flatmap pixels apart, so
# at three-panels-wide their labels collide. Nudged apart by hand; the outlines
# stay where they are.
LABEL_NUDGE = {'VISrl': (55, 0), 'VISal': (-55, 0)}
# Streamlines sample the volume one path at a time, so a map straight out of the
# projector carries fine radial streaks that are sampling, not signal. Two
# flatmap pixels of blur removes them and is far below any area's size.
SMOOTH_PX = 2.0
BANDS = [('supragranular  L1 + L2/3', ('Isocortex layer 1', 'Isocortex layer 2/3')),
         ('granular  L4', ('Isocortex layer 4',)),
         ('infragranular  L5 + L6', ('Isocortex layer 5', 'Isocortex layer 6a', 'Isocortex layer 6b'))]


def fold(v):
    """Average the two hemispheres, as every other v2 figure does."""
    h = v.shape[2] // 2
    return np.nanmean(np.stack([v[:, :, :h], v[:, :, v.shape[2] - h:][:, :, ::-1]]), axis=0)


def fold_n(n):
    h = n.shape[2] // 2
    return np.maximum(n[:, :, :h], n[:, :, n.shape[2] - h:][:, :, ::-1])


def cohort_volume(cohort, reading, min_n):
    """(value, mask) for one cohort on the 10 um CCF grid.

    Hemispheres folded and mirrored back into both halves; value is zero where
    the mask is zero, so a sum/sum ratio downstream ignores it.
    """
    mean = fold(np.load(os.path.join(CCF_ROOT, cohort, f'{reading}_mean.npy')))
    n = fold_n(np.load(os.path.join(CCF_ROOT, cohort, f'{reading}_n.npy')))
    ok = (n >= min_n) & np.isfinite(mean)
    half_v = np.where(ok, mean, 0).astype(np.float32)
    half_m = ok.astype(np.float32)

    full_v = np.concatenate([half_v, half_v[:, :, ::-1]], axis=2)
    full_m = np.concatenate([half_m, half_m[:, :, ::-1]], axis=2)
    up = lambda a: np.repeat(np.repeat(np.repeat(a, 2, 0), 2, 1), 2, 2)   # 20 um -> 10 um, exactly
    return up(full_v), up(full_m)


def layer_bands():
    """(name, first depth index, last depth index) per band, from the Allen depths.

    avg_layer_depths.json gives the cumulative depth of each layer's lower
    border in um; with thickness_type='normalized_layers' the slab is built to
    exactly those thicknesses, so the band edges are the same numbers over 10.
    """
    d = json.load(open(os.path.join(ASSETS, 'avg_layer_depths.json')))
    order = ['2/3', '4', '5', '6a', '6b', 'wm']
    lower = {k: d[k] for k in order}
    thick = {'Isocortex layer 1': d['2/3']}
    prev = d['2/3']
    for key, name in (('4', 'Isocortex layer 2/3'), ('5', 'Isocortex layer 4'), ('6a', 'Isocortex layer 5'),
                      ('6b', 'Isocortex layer 6a'), ('wm', 'Isocortex layer 6b')):
        thick[name] = lower[key] - prev
        prev = lower[key]
    # the slab's depth axis is the layers stacked in order, 10 um per bin
    edges, at = {}, 0
    for name in ['Isocortex layer 1', 'Isocortex layer 2/3', 'Isocortex layer 4',
                 'Isocortex layer 5', 'Isocortex layer 6a', 'Isocortex layer 6b']:
        n = int(round(thick[name] / 10))
        edges[name] = (at, at + n)
        at += n
    return thick, edges


def smooth(img, sigma=SMOOTH_PX):
    """Blur within the mapped cortex only, so the outline stays crisp."""
    from scipy.ndimage import gaussian_filter
    ok = np.isfinite(img).astype(np.float32)
    num = gaussian_filter(np.where(ok > 0, img, 0).astype(np.float32), sigma)
    den = gaussian_filter(ok, sigma)
    return np.where(ok > 0, num / np.maximum(den, 1e-6), np.nan)


def draw(ax, img, cmap, lim, title, border_sets, label_xy):
    # border_sets is a LIST, not a merged dict: the two hemispheres carry the
    # same acronyms as keys, so merging them silently keeps one hemisphere only.
    m = np.ma.masked_invalid(img)
    ax.imshow(np.zeros(img.shape + (4,)), origin='upper')
    h = ax.imshow(m, cmap=cmap, vmin=lim[0], vmax=lim[1], origin='upper', interpolation='nearest')
    for borders in border_sets:
        for acro, coords in borders.items():
            hi = acro in HIGHLIGHT
            ax.plot(coords[:, 0], coords[:, 1], c='w' if hi else '#bbbbbb', lw=1.1 if hi else 0.35)
    for acro, (x, y) in label_xy.items():
        hi = acro in HIGHLIGHT
        dx, dy = LABEL_NUDGE.get(acro, (0, 0))
        t = ax.text(x + dx, y + dy, acro, color='w', fontsize=8 if hi else 6, ha='center', va='center',
                    fontweight='bold' if hi else 'normal')
        t.set_path_effects([path_effects.withStroke(linewidth=1.6, foreground='k')])
    ax.set_title(title, color='w', fontsize=10)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor('k')
    for s in ax.spines.values():
        s.set_visible(False)
    return h


def main(readings):
    from ccf_streamlines.projection import Isocortex2dProjector, Isocortex3dProjector, BoundaryFinder

    proj_file = os.path.join(ASSETS, 'flatmap_butterfly.h5')
    path_file = os.path.join(ASSETS, 'surface_paths_10_v3.h5')

    hot = plt.get_cmap('hot')
    hot_cut = LinearSegmentedColormap.from_list('hot_cut', hot(np.linspace(0, 0.82, 256)))
    for cm in (hot_cut,):
        cm.set_bad((0, 0, 0, 0))
    puor = plt.get_cmap('PuOr_r').copy(); puor.set_bad((0, 0, 0, 0))
    rdbu = plt.get_cmap('RdBu_r').copy(); rdbu.set_bad((0, 0, 0, 0))

    print('loading borders...', flush=True)
    bf = BoundaryFinder(projected_atlas_file=os.path.join(ASSETS, 'flatmap_butterfly.nrrd'),
                        labels_file=os.path.join(ASSETS, 'labelDescription_ITKSNAPColor.txt'))
    borders = {k: v for k, v in bf.region_boundaries().items() if len(v)}
    borders_r = {k: v for k, v in bf.region_boundaries(
        hemisphere='right_for_both', view_space_for_other_hemisphere='flatmap_butterfly').items() if len(v)}
    label_xy = {a: borders[a].mean(axis=0) for a in LABEL_AREAS if a in borders}

    print('loading streamlines (0.5 GB, takes a minute)...', flush=True)
    t0 = time.time()
    p2 = Isocortex2dProjector(proj_file, path_file, hemisphere='both',
                              view_space_for_other_hemisphere='flatmap_butterfly')
    thick, edges = layer_bands()
    p3 = Isocortex3dProjector(proj_file, path_file, thickness_type='normalized_layers',
                              layer_thicknesses=thick,
                              streamline_layer_thickness_file=os.path.join(ASSETS, 'cortical_layers_10_v2.h5'),
                              hemisphere='both', view_space_for_other_hemisphere='flatmap_butterfly')
    print(f'  projectors ready, {time.time() - t0:.0f} s', flush=True)

    for reading in readings:
        t0 = time.time()
        signed = reading == 'zref'
        cmap_mean = puor if signed else hot_cut
        vm, dl = MEAN_VMAX[reading], DIFF_LIM[reading]
        lim_mean = (-vm, vm) if signed else (0, vm)

        flat, slab = {}, {}
        for cohort, min_n in ((YOUNG, MIN_N_YOUNG), (ADULT, MIN_N_ADULT)):
            vol, msk = cohort_volume(cohort, reading, min_n)
            # depth average as a ratio of sums: a voxel nobody covered adds nothing
            #
            # project_volume indexes its output [x, y], while region_boundaries
            # hands back (x, y) pairs to plot; imshow wants [row, col] = [y, x].
            # Everything below is therefore transposed once, here, and the
            # borders can then be drawn exactly as they come.
            s_v = p2.project_volume(vol, kind='sum').T
            s_m = p2.project_volume(msk, kind='sum').T
            flat[cohort] = smooth(np.where(s_m > 0, s_v / np.maximum(s_m, 1e-6), np.nan))
            slab[cohort] = (p3.project_volume(vol).swapaxes(0, 1),
                            p3.project_volume(msk).swapaxes(0, 1))
            del vol, msk
            print(f'  {reading} {cohort} projected   {time.time() - t0:.0f} s', flush=True)

        diff = flat[YOUNG] - flat[ADULT] if signed else np.log2(
            np.maximum(flat[YOUNG], 0.02) / np.maximum(flat[ADULT], 0.02))

        fig, axes = plt.subplots(1, 3, figsize=(19, 5.2), facecolor='k')
        panels = ((flat[YOUNG], cmap_mean, lim_mean, f'young (n = 7)   {reading}'),
                  (flat[ADULT], cmap_mean, lim_mean, f'adult (n = 10)   {reading}'),
                  (diff, rdbu, (-dl, dl), 'young - adult' if signed else 'log2( young / adult )'))
        for ax, (im, cm, lim, ttl) in zip(axes, panels):
            h = draw(ax, im, cm, lim, ttl, (borders, borders_r), label_xy)
            cb = fig.colorbar(h, ax=ax, fraction=0.03, pad=0.01)
            cb.ax.yaxis.set_tick_params(color='w', labelcolor='w')
        fig.suptitle(f'Cortical flatmap (Allen streamlines), averaged through the full depth, {SMOOTH_PX:g} px blur  -  {reading}',
                     color='w', fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(os.path.join(OUT, f'{reading}_flatmap.png'), dpi=110, facecolor='k')
        plt.close(fig)

        fig, axes = plt.subplots(3, 3, figsize=(19, 13), facecolor='k')
        for row, (bname, keys) in enumerate(BANDS):
            lo, hi = edges[keys[0]][0], edges[keys[-1]][1]
            band = {}
            for cohort in (YOUNG, ADULT):
                sv, sm = slab[cohort]
                num = sv[:, :, lo:hi].sum(axis=2)
                den = sm[:, :, lo:hi].sum(axis=2)
                band[cohort] = smooth(np.where(den > 0, num / np.maximum(den, 1e-6), np.nan))
            d = band[YOUNG] - band[ADULT] if signed else np.log2(
                np.maximum(band[YOUNG], 0.02) / np.maximum(band[ADULT], 0.02))
            panels = ((band[YOUNG], cmap_mean, lim_mean, f'young   {bname}'),
                      (band[ADULT], cmap_mean, lim_mean, f'adult   {bname}'),
                      (d, rdbu, (-dl, dl), f'young - adult   {bname}'))
            for ax, (im, cm, lim, ttl) in zip(axes[row], panels):
                h = draw(ax, im, cm, lim, ttl, (borders, borders_r), label_xy)
                cb = fig.colorbar(h, ax=ax, fraction=0.03, pad=0.01)
                cb.ax.yaxis.set_tick_params(color='w', labelcolor='w')
        fig.suptitle(f'Cortical flatmap by layer, depth normalised per streamline, {SMOOTH_PX:g} px blur  -  {reading}',
                     color='w', fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(os.path.join(OUT, f'{reading}_flatmap_layers.png'), dpi=110, facecolor='k')
        plt.close(fig)
        print(f'{reading}: wrote {reading}_flatmap.png and {reading}_flatmap_layers.png   '
              f'{time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    main(sys.argv[1:] or ['zref'])
