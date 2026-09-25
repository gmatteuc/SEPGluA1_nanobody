"""
v2, close inspection of one reading: four views of exactly the same numbers.

The cohort volumes are prepared ONCE here and every output below is drawn from
them, so a coronal frame, the video it comes from and the flatmaps cannot
disagree with each other about contrast, coverage or smoothing:

  detail_plane<P>_<reading>.png        one coronal plane, young | adult | difference
  detail_video_<reading>.mp4           the same three panels, plane by plane
  detail_flatmap_<reading>.png         the isocortex unrolled, full cortical depth
  detail_flatmap_layers_<reading>.png  the same, supragranular / granular / infragranular

Why it exists beside the standard figures. The standard videos and slice figures
run every reading at one fixed scale so cohorts stay comparable. These are for
looking closely at one reading: the colour range is tightened until cortex is
readable (the hippocampus then clips, which is the point), and a light 3D
smoothing is applied to the volumes before anything is drawn.

The smoothing is cosmetic and it is declared in every title. It is mask
normalised, so tissue is never averaged with the black outside it, and it never
touches the statistics -- the region tables and tests read the unsmoothed
per-mouse volumes and know nothing about this script. What it removes is
sampling, not signal: sections sit ~150 um apart and the registered volumes
interpolate between them, which leaves faint coronal banding, and the flatmap
samples one streamline at a time, which leaves fine radial streaks. The default
sigma is 3 x 1 x 1 voxels, i.e. 60 um along AP and 20 um across -- anisotropic
because the banding is, and far below the size of any area.

  tools\\venv_flat\\Scripts\\python.exe v2_inspect.py [reading ...]
        [--plane 790] [--vmax 1.0] [--dlim 0.5] [--smooth 3,1,1] [--no-video] [--no-flatmap]

`--smooth` takes one number or three, comma separated, as sigma in 20 um voxels
along (AP, DV, ML); 0 turns it off.

The flatmaps need ccf_streamlines, which brings its own numpy and scikit-image,
so this script runs in its own environment, tools\\venv_flat, rather than the
analysis one:

  py -m venv tools\\venv_flat
  tools\\venv_flat\\Scripts\\python.exe -m pip install ccf-streamlines matplotlib nibabel imageio-ffmpeg

Its assets go in data\\atlas_flatmap, about 0.6 GB, fetched once from
https://download.alleninstitute.org/informatics-archive/current-release/
mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/ :

  streamlines/surface_paths_10_v3.h5        where each streamline runs (0.5 GB)
  view_lookup/flatmap_butterfly.h5          which streamline lands on which pixel
  master_updated/flatmap_butterfly.nrrd     the areas, for their borders
  master_updated/labelDescription_ITKSNAPColor.txt   and their names
  cortical_metrics/avg_layer_depths.json    where the layers sit
  cortical_metrics/cortical_layers_10_v2.h5 and where they sit per streamline
"""

import csv
import json
import os
import sys
import time

import numpy as np
import nibabel as nib
import imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as path_effects
from scipy.ndimage import center_of_mass, gaussian_filter

DATA = r'D:\sep_histology\data'
ASSETS = os.path.join(DATA, 'atlas_flatmap')
CCF_ROOT = os.path.join(DATA, 'comparisons_v2', 'ccf')
OUT = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult')
CSV_MAP = os.path.join(DATA, 'atlas', 'parcellation_to_parcellation_term_membership.csv')

YOUNG, ADULT = 'young', 'adult'
MIN_N_YOUNG, MIN_N_ADULT = 2, 5          # as in v2_compare, so the maps agree
# The floor that keeps log2 finite where a reading is near zero, the same value
# v2_cohort uses. It bites only the readings compared as a ratio, never zref.
LOG2_FLOOR = 0.02
# Kept in step with v2_cohort.SIGNED_READINGS by hand, because this script runs
# in the flatmap environment and cannot import the analysis chain. A signed
# reading is a position within a brain's own range: compared by difference,
# drawn diverging, scaled -v..+v.
SIGNED_READINGS = ('zref',)
VMAX = {'zref': 1.0, 'cref': 2.0, 'subref': 2.0, 'ratio': 2.0, 'sepratio': 0.6}
DLIM = {'zref': 0.5, 'cref': 1.0, 'subref': 1.0, 'ratio': 1.0, 'sepratio': 1.0}
PLANE = 790                              # CCF plane at 10 um, where RL and AL are cut
# Sigma in 20 um voxels, along (AP, DV, ML). Anisotropic on purpose: the
# artefact being removed is anisotropic. Sections sit ~150 um apart -- 7.5
# voxels -- and the registered volumes interpolate between them, so the banding
# runs along AP and nothing else. 3 voxels of AP blur is 60 um, well under a
# section spacing and two orders below the distance between V1 and RL.
SMOOTH = (3.0, 1.0, 1.0)
FPS = 12
MIN_LABEL_AREA = 150                     # 20 um voxels in a coronal plane

# the areas the argument is about: outlined in white on the flatmap and named
HIGHLIGHT = ('VISp', 'VISrl', 'VISal', 'SSp-bfd')
LABEL_AREAS = ['VISp', 'VISrl', 'VISal', 'VISam', 'SSp-bfd', 'SSs', 'AUDp', 'MOp', 'RSPd', 'ACAd', 'TEa']
# RL and AL are neighbours whose centroids are 100 flatmap pixels apart, which
# is not enough at three panels wide
LABEL_NUDGE = {'VISrl': (60, 0), 'VISal': (-60, 0)}
BANDS = [('supragranular  L1 + L2/3', ('Isocortex layer 1', 'Isocortex layer 2/3')),
         ('granular  L4', ('Isocortex layer 4',)),
         ('infragranular  L5 + L6', ('Isocortex layer 5', 'Isocortex layer 6a', 'Isocortex layer 6b'))]


def save_figure(fig, path, dpi):
    """Save, and if the file is open in a viewer say so instead of dying.

    Windows refuses to overwrite a PNG an image viewer holds open, and these
    figures are made to be looked at while the next one is being drawn.
    """
    try:
        fig.savefig(path, dpi=dpi, facecolor='k')
    except OSError:
        alt = path.replace('.png', '_new.png')
        fig.savefig(alt, dpi=dpi, facecolor='k')
        print(f'  NOTE: {os.path.basename(path)} is open elsewhere; wrote {os.path.basename(alt)}', flush=True)


def cohort_size(cohort):
    """How many brains are behind a cohort, read from the list v2_cohort wrote.

    Never a literal: a hardcoded 'n = 6' once survived into figures built from
    seven brains, and a caption that cannot go stale is worth four lines.
    """
    with open(os.path.join(CCF_ROOT, cohort, 'mice.txt'), encoding='utf-8') as fh:
        return sum(1 for line in fh if line.strip())


def fold(v):
    """Average the two hemispheres of an (AP, DV, ML) volume, NaN-aware."""
    h = v.shape[2] // 2
    return np.nanmean(np.stack([v[:, :, :h], v[:, :, v.shape[2] - h:][:, :, ::-1]]), axis=0)


def fold_n(n):
    """The same fold for an n map: a voxel counts if either hemisphere had it."""
    h = n.shape[2] // 2
    return np.maximum(n[:, :, :h], n[:, :, n.shape[2] - h:][:, :, ::-1])


def prepare(reading, sigma):
    """(value, mask) per cohort, hemispheres folded, at 20 um. One source for all four figures.

    The mask is left alone by the smoothing: blurring a value inside the tissue
    is cosmetic, growing the tissue would not be.
    """
    out = {}
    for cohort, min_n in ((YOUNG, MIN_N_YOUNG), (ADULT, MIN_N_ADULT)):
        mean = fold(np.load(os.path.join(CCF_ROOT, cohort, f'{reading}_mean.npy')))
        n = fold_n(np.load(os.path.join(CCF_ROOT, cohort, f'{reading}_n.npy')))
        m = ((n >= min_n) & np.isfinite(mean)).astype(np.float32)
        v = np.where(m > 0, mean, 0).astype(np.float32)
        if np.any(sigma):
            num = gaussian_filter(v, sigma)
            den = gaussian_filter(m, sigma)
            v = np.where(m > 0, num / np.maximum(den, 1e-6), 0).astype(np.float32)
        out[cohort] = (v, m)
    return out


def difference(vals, signed):
    y, a = vals[YOUNG], vals[ADULT]
    both = (y[1] > 0) & (a[1] > 0)
    d = (y[0] - a[0]) if signed else np.log2(np.maximum(y[0], LOG2_FLOOR) / np.maximum(a[0], LOG2_FLOOR))
    return np.where(both, d, np.nan), both


def shown(value, mask):
    return np.where(mask > 0, value, np.nan)


# ----------------------------------------------------------------- coronal

def annotation_half():
    ann = np.asarray(nib.load(os.path.join(DATA, 'atlas', 'annotation_10.nii.gz')).dataobj)[::2, ::2, ::2]
    return ann[:, :, :ann.shape[2] // 2]


def boundaries(lab):
    b = np.zeros(lab.shape, bool)
    b[1:, :] |= lab[1:, :] != lab[:-1, :]
    b[:, 1:] |= lab[:, 1:] != lab[:, :-1]
    return b & (lab > 0)


def coronal_frame(fig, axes, caxes, k, panels, ann_h, acro, header):
    lab = ann_h[k]; inside = lab > 0; bnd = boundaries(lab)
    for ax, cax, (im, cmap, lim, ttl) in zip(axes, caxes, panels):
        ax.clear(); cax.clear()
        bg = np.zeros(lab.shape + (4,)); bg[inside] = (0.23, 0.23, 0.23, 1.0)
        ax.imshow(bg, origin='upper', interpolation='nearest', aspect='equal')
        h = ax.imshow(np.ma.masked_invalid(np.where(inside, im, np.nan)), cmap=cmap,
                      vmin=lim[0], vmax=lim[1], origin='upper', interpolation='nearest', aspect='equal')
        ov = np.zeros(lab.shape + (4,)); ov[bnd] = (0.75, 0.75, 0.75, 0.9)
        ax.imshow(ov, origin='upper', interpolation='nearest', aspect='equal')
        for idx in np.unique(lab):
            if idx == 0:
                continue
            m = lab == idx
            if m.sum() < MIN_LABEL_AREA:
                continue
            cy, cx = center_of_mass(m)
            ax.text(cx, cy, acro.get(int(idx), ''), color='w', fontsize=5.5, ha='center', va='center')
        ax.set_facecolor('k'); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(ttl, color='w', fontsize=12)
        cb = fig.colorbar(h, cax=cax); cb.ax.yaxis.set_tick_params(color='w', labelcolor='w')
    fig.texts.clear()
    fig.text(0.5, 0.93, header, color='w', fontsize=12, ha='center')


def coronal(reading, vals, diff, both, plane, lim_mean, cmaps, sigma_txt, want_video, n, signed):
    ann_h = annotation_half()
    acro = {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['parcellation_term_set_name'] == 'structure':
                acro[int(row['parcellation_index'])] = row['parcellation_term_acronym']
    cmap_mean, rdbu = cmaps

    fig = plt.figure(figsize=(19.2, 7.6), dpi=100, facecolor='k')
    axes = [fig.add_axes([0.02 + i * 0.325, 0.05, 0.27, 0.82]) for i in range(3)]
    caxes = [fig.add_axes([0.295 + i * 0.325, 0.12, 0.009, 0.68]) for i in range(3)]

    def panels_at(k):
        return ((shown(vals[YOUNG][0][k], vals[YOUNG][1][k]), cmap_mean, lim_mean[0], f'young (n = {n[YOUNG]})   {reading}'),
                (shown(vals[ADULT][0][k], vals[ADULT][1][k]), cmap_mean, lim_mean[0], f'adult (n = {n[ADULT]})   {reading}'),
                (diff[k], rdbu, lim_mean[1], 'young - adult' if signed else 'log2( young / adult )'))

    def header_at(k):
        return (f'CCF plane {2 * k} / 10 um    young {n[YOUNG]} brains, adult {n[ADULT]} brains    '
                f'{sigma_txt}, colour range tightened for cortex')

    k = plane // 2
    coronal_frame(fig, axes, caxes, k, panels_at(k), ann_h, acro, header_at(k))
    out = os.path.join(OUT, f'detail_plane{plane}_{reading}.png')
    save_figure(fig, out, 100)
    print(f'  wrote {os.path.basename(out)}', flush=True)

    if want_video:
        t0 = time.time()
        frames = [i for i in range(ann_h.shape[0]) if both[i].sum() > 200]
        out = os.path.join(OUT, f'detail_video_{reading}.mp4')
        writer = imageio_ffmpeg.write_frames(out, (1920, 760), fps=FPS, quality=7, macro_block_size=8)
        writer.send(None)
        for i in frames:
            coronal_frame(fig, axes, caxes, i, panels_at(i), ann_h, acro, header_at(i))
            fig.canvas.draw()
            writer.send(np.ascontiguousarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3]))
        writer.close()
        print(f'  wrote {os.path.basename(out)}, {len(frames)} frames, {time.time() - t0:.0f} s', flush=True)
    plt.close(fig)


# ---------------------------------------------------------------- flatmap

def layer_thicknesses():
    """Thickness of each cortical layer in um, from the Allen's average depths.

    The file gives the depth of each layer's lower border below the pia, so the
    thicknesses are the differences between them.
    """
    d = json.load(open(os.path.join(ASSETS, 'avg_layer_depths.json')))
    names = ['Isocortex layer 1', 'Isocortex layer 2/3', 'Isocortex layer 4',
             'Isocortex layer 5', 'Isocortex layer 6a', 'Isocortex layer 6b']
    thick, prev = {}, 0.0
    for name, low in zip(names, [d['2/3'], d['4'], d['5'], d['6a'], d['6b'], d['wm']]):
        thick[name] = low - prev
        prev = low
    return thick


def band_edges(p3, slab_depth):
    """First and last depth bin of each layer, as the projector actually built them.

    Ask the projector, never recompute. The slab does NOT have one bin per 10 um:
    it has one per sample along a streamline -- 200 of them -- and the layers are
    given bins in proportion to their thickness. Deriving the edges here from the
    layer depths instead put 96 bins where there are 200, which silently labelled
    the middle of layer 2/3 as layer 4 and layer 4 as layer 6.
    """
    bins = p3.reference_layer_thicknesses_in_voxels()
    edges, at = {}, 0
    for name in p3.ISOCORTEX_LAYER_KEYS:
        edges[name] = (at, at + bins[name])
        at += bins[name]
    if at != slab_depth:
        raise SystemExit(f'layer bins sum to {at} but the slab is {slab_depth} deep; '
                         f'the band edges cannot be trusted.')
    return edges


def draw_flat(ax, img, cmap, lim, title, border_sets, label_xy):
    # border_sets is a LIST: the two hemispheres carry the same acronyms as keys,
    # so merging them into one dict silently keeps a single hemisphere.
    ax.imshow(np.zeros(img.shape + (4,)), origin='upper')
    h = ax.imshow(np.ma.masked_invalid(img), cmap=cmap, vmin=lim[0], vmax=lim[1],
                  origin='upper', interpolation='nearest')
    for borders in border_sets:
        for acro, coords in borders.items():
            hi = acro in HIGHLIGHT
            ax.plot(coords[:, 0], coords[:, 1], c='#e6e6e6' if hi else '#9a9a9a',
                    lw=0.7 if hi else 0.3, alpha=0.9 if hi else 0.65)
    for acro, (x, y) in label_xy.items():
        dx, dy = LABEL_NUDGE.get(acro, (0, 0))
        t = ax.text(x + dx, y + dy, acro, color='w', fontsize=7, ha='center', va='center')
        # the map runs from near-white to dark red, so plain text loses either
        # end of it; a soft dark halo keeps it readable without shouting
        t.set_path_effects([path_effects.withStroke(linewidth=1.3, foreground=(0, 0, 0, 0.55))])
    ax.set_title(title, color='w', fontsize=10)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor('k')
    for s in ax.spines.values():
        s.set_visible(False)
    return h


def flatmaps(reading, vals, signed, lim_mean, cmaps, sigma_txt, n):
    from ccf_streamlines.projection import Isocortex2dProjector, Isocortex3dProjector, BoundaryFinder
    cmap_mean, rdbu = cmaps
    proj_file = os.path.join(ASSETS, 'flatmap_butterfly.h5')
    path_file = os.path.join(ASSETS, 'surface_paths_10_v3.h5')

    bf = BoundaryFinder(projected_atlas_file=os.path.join(ASSETS, 'flatmap_butterfly.nrrd'),
                        labels_file=os.path.join(ASSETS, 'labelDescription_ITKSNAPColor.txt'))
    left = {k: v for k, v in bf.region_boundaries().items() if len(v)}
    right = {k: v for k, v in bf.region_boundaries(
        hemisphere='right_for_both', view_space_for_other_hemisphere='flatmap_butterfly').items() if len(v)}
    label_xy = {a: left[a].mean(axis=0) for a in LABEL_AREAS if a in left}

    p2 = Isocortex2dProjector(proj_file, path_file, hemisphere='both',
                              view_space_for_other_hemisphere='flatmap_butterfly')
    p3 = Isocortex3dProjector(proj_file, path_file, thickness_type='normalized_layers',
                              layer_thicknesses=layer_thicknesses(),
                              streamline_layer_thickness_file=os.path.join(ASSETS, 'cortical_layers_10_v2.h5'),
                              hemisphere='both', view_space_for_other_hemisphere='flatmap_butterfly')

    def to_10um(half):
        full = np.concatenate([half, half[:, :, ::-1]], axis=2)
        return np.repeat(np.repeat(np.repeat(full, 2, 0), 2, 1), 2, 2)   # exact inverse of the block mean

    flat, slab = {}, {}
    for cohort in (YOUNG, ADULT):
        v10, m10 = to_10um(vals[cohort][0]), to_10um(vals[cohort][1])
        # project_volume indexes [x, y] while region_boundaries returns (x, y) to
        # plot, and imshow wants [row, col] = [y, x]: transpose once, here
        s_v = p2.project_volume(v10, kind='sum').T
        s_m = p2.project_volume(m10, kind='sum').T
        flat[cohort] = np.where(s_m > 0, s_v / np.maximum(s_m, 1e-6), np.nan)
        slab[cohort] = (p3.project_volume(v10).swapaxes(0, 1), p3.project_volume(m10).swapaxes(0, 1))
        del v10, m10
    edges = band_edges(p3, slab[YOUNG][0].shape[2])

    def diff_of(y, a):
        return (y - a) if signed else np.log2(np.maximum(y, LOG2_FLOOR) / np.maximum(a, LOG2_FLOOR))

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.2), facecolor='k')
    panels = ((flat[YOUNG], cmap_mean, lim_mean[0], f'young (n = {n[YOUNG]})   {reading}'),
              (flat[ADULT], cmap_mean, lim_mean[0], f'adult (n = {n[ADULT]})   {reading}'),
              (diff_of(flat[YOUNG], flat[ADULT]), rdbu, lim_mean[1], 'young - adult' if signed else 'log2( young / adult )'))
    for ax, (im, cm, lim, ttl) in zip(axes, panels):
        h = draw_flat(ax, im, cm, lim, ttl, (left, right), label_xy)
        cb = fig.colorbar(h, ax=ax, fraction=0.03, pad=0.01)
        cb.ax.yaxis.set_tick_params(color='w', labelcolor='w')
    fig.suptitle(f'Cortical flatmap, averaged through the full depth  -  {reading}, {sigma_txt}',
                 color='w', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(fig, os.path.join(OUT, f'detail_flatmap_{reading}.png'), 110)
    plt.close(fig)

    fig, axes = plt.subplots(3, 3, figsize=(19, 13), facecolor='k')
    for row, (bname, keys) in enumerate(BANDS):
        lo, hi = edges[keys[0]][0], edges[keys[-1]][1]
        band = {}
        for cohort in (YOUNG, ADULT):
            sv, sm = slab[cohort]
            num = sv[:, :, lo:hi].sum(axis=2)
            den = sm[:, :, lo:hi].sum(axis=2)
            band[cohort] = np.where(den > 0, num / np.maximum(den, 1e-6), np.nan)
        panels = ((band[YOUNG], cmap_mean, lim_mean[0], f'young   {bname}'),
                  (band[ADULT], cmap_mean, lim_mean[0], f'adult   {bname}'),
                  (diff_of(band[YOUNG], band[ADULT]), rdbu, lim_mean[1], f'young - adult   {bname}'))
        for ax, (im, cm, lim, ttl) in zip(axes[row], panels):
            h = draw_flat(ax, im, cm, lim, ttl, (left, right), label_xy)
            cb = fig.colorbar(h, ax=ax, fraction=0.03, pad=0.01)
            cb.ax.yaxis.set_tick_params(color='w', labelcolor='w')
    fig.suptitle(f'Cortical flatmap by layer, depth normalised per streamline  -  {reading}, {sigma_txt}',
                 color='w', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, os.path.join(OUT, f'detail_flatmap_layers_{reading}.png'), 110)
    plt.close(fig)
    print(f'  wrote detail_flatmap_{reading}.png and detail_flatmap_layers_{reading}.png', flush=True)


def main(readings, plane, vmax, dlim, sigma, want_video, want_flatmap):
    hot = plt.get_cmap('hot')
    hot_cut = LinearSegmentedColormap.from_list('hot_cut', hot(np.linspace(0, 0.82, 256)))
    hot_cut.set_bad((0, 0, 0, 0))
    puor = plt.get_cmap('PuOr_r').copy(); puor.set_bad((0, 0, 0, 0))
    rdbu = plt.get_cmap('RdBu_r').copy(); rdbu.set_bad((0, 0, 0, 0))

    sig = np.atleast_1d(sigma).astype(float)
    sigma_txt = ('no smoothing' if not np.any(sig)
                 else 'smoothed sigma ' + ' x '.join(f'{v * 20:.0f}' for v in (sig if sig.size > 1 else np.repeat(sig, 3)))
                      + ' um (AP x DV x ML)')
    n = {c: cohort_size(c) for c in (YOUNG, ADULT)}
    for reading in readings:
        t0 = time.time()
        signed = reading in SIGNED_READINGS
        vm = VMAX[reading] if vmax is None else vmax
        dl = DLIM[reading] if dlim is None else dlim
        lim_mean = ((-vm, vm) if signed else (0, vm), (-dl, dl))
        cmaps = (puor if signed else hot_cut, rdbu)

        print(f'{reading}: preparing volumes ({sigma_txt})', flush=True)
        vals = prepare(reading, sigma)
        diff, both = difference(vals, signed)

        coronal(reading, vals, diff, both, plane, lim_mean, cmaps, sigma_txt, want_video, n, signed)
        if want_flatmap:
            flatmaps(reading, vals, signed, lim_mean, cmaps, sigma_txt, n)
        print(f'{reading}: done in {time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    argv, readings, opts, flags = sys.argv[1:], [], {}, set()
    i = 0
    while i < len(argv):
        if argv[i] in ('--plane', '--vmax', '--dlim', '--smooth'):
            opts[argv[i][2:]] = argv[i + 1]; i += 2
        elif argv[i] in ('--no-video', '--no-flatmap'):
            flags.add(argv[i]); i += 1
        else:
            readings.append(argv[i]); i += 1
    s = ([float(x) for x in opts['smooth'].split(',')] if 'smooth' in opts else list(SMOOTH))
    main(readings or ['zref'],
         plane=int(opts.get('plane', PLANE)),
         vmax=float(opts['vmax']) if 'vmax' in opts else None,
         dlim=float(opts['dlim']) if 'dlim' in opts else None,
         sigma=(s[0] if len(s) == 1 else s),
         want_video='--no-video' not in flags,
         want_flatmap='--no-flatmap' not in flags)
