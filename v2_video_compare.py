"""
Young and adult side by side, plane by plane, in the adult CCF.

Three panels per frame: the young cohort mean, the adult cohort mean on the
same colour scale, and log2(young / adult) beside them, so a pattern can be
compared where it sits rather than by flicking between two windows. Atlas
outlines and acronyms on all three; the header says which CCF plane and how
many brains are behind each side there.

Everything comes from the cohort volumes v2_cohort built in CCF, the same ones
v2_compare makes its figures from. The colour scale of the two means is fixed
and shared (linear, not log); only the third panel is a log2 ratio.

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_video_compare.py [reading ...]
"""

import csv
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
from scipy.ndimage import center_of_mass

from v2_per_mouse import CSV_MAP, DATA
from v2_cohort import COHORTS, OUT_ROOT as CCF_ROOT
from v2_compare import MIN_N_YOUNG, MIN_N_ADULT, YOUNG, fold, fold_n

OUT = os.path.join(DATA, 'comparisons_v2', 'young_vs_adult')
MEAN_VMAX = {'cref': 2.0, 'ratio': 2.0}
LOG2_LIM = 1.5
FPS = 12
MIN_LABEL_AREA = 150


def boundaries(lab):
    b = np.zeros(lab.shape, bool)
    b[1:, :] |= lab[1:, :] != lab[:-1, :]
    b[:, 1:] |= lab[:, 1:] != lab[:, :-1]
    return b & (lab > 0)


def main(readings):
    acro = {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['parcellation_term_set_name'] == 'structure':
                acro[int(row['parcellation_index'])] = row['parcellation_term_acronym']
    hot = plt.get_cmap('hot')
    hot_cut = LinearSegmentedColormap.from_list('hot_cut', hot(np.linspace(0, 0.82, 256)))
    hot_cut.set_bad((0, 0, 0, 0))
    rdbu = plt.get_cmap('RdBu_r').copy(); rdbu.set_bad((0, 0, 0, 0))

    ann = np.asarray(nib.load(os.path.join(DATA, 'atlas', 'annotation_10.nii.gz')).dataobj)[::2, ::2, ::2]
    ann_h = ann[:, :, :ann.shape[2] // 2]
    y_n = fold_n(np.load(os.path.join(CCF_ROOT, YOUNG, 'cref_n.npy')))
    a_n = fold_n(np.load(os.path.join(CCF_ROOT, 'adult', 'cref_n.npy')))

    for reading in readings:
        t0 = time.time()
        y = fold(np.load(os.path.join(CCF_ROOT, YOUNG, f'{reading}_mean.npy')))
        a = fold(np.load(os.path.join(CCF_ROOT, 'adult', f'{reading}_mean.npy')))
        ok_y = (y_n >= MIN_N_YOUNG) & np.isfinite(y)
        ok_a = (a_n >= MIN_N_ADULT) & np.isfinite(a)
        both = ok_y & ok_a & (ann_h > 0)
        log2 = np.where(both, np.log2(np.maximum(y, 0.02) / np.maximum(a, 0.02)), np.nan)
        frames = [k for k in range(ann_h.shape[0]) if both[k].sum() > 200]
        out = os.path.join(OUT, f'video_side_by_side_{reading}.mp4')
        writer = imageio_ffmpeg.write_frames(out, (1920, 760), fps=FPS, quality=7, macro_block_size=8)
        writer.send(None)
        fig = plt.figure(figsize=(19.2, 7.6), dpi=100, facecolor='k')
        axes = [fig.add_axes([0.02 + i * 0.325, 0.05, 0.27, 0.82]) for i in range(3)]
        caxes = [fig.add_axes([0.295 + i * 0.325, 0.12, 0.009, 0.68]) for i in range(3)]
        for k in frames:
            lab = ann_h[k]; inside = lab > 0; bnd = boundaries(lab)
            panels = ((np.where(ok_y[k], y[k], np.nan), hot_cut, (0, MEAN_VMAX[reading]),
                       f'young (n = {len(COHORTS[YOUNG])})   {reading}'),
                      (np.where(ok_a[k], a[k], np.nan), hot_cut, (0, MEAN_VMAX[reading]),
                       f'adult (n = {len(COHORTS["adult"])})   {reading}'),
                      (log2[k], rdbu, (-LOG2_LIM, LOG2_LIM), 'log2( young / adult )'))
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
            fig.text(0.5, 0.93, f'CCF plane {2 * k} / 10 um    young: {int(np.nanmax(np.where(ok_y[k], y_n[k], 0)))} brains   '
                                f'adult: {int(np.nanmax(np.where(ok_a[k], a_n[k], 0)))} brains   '
                                f'(means on a linear scale, shared range; only the right panel is log2)',
                     color='w', fontsize=12, ha='center')
            fig.canvas.draw()
            writer.send(np.ascontiguousarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3]))
        writer.close(); plt.close(fig)
        print(f'{reading:6s} {len(frames)} frames -> {out}   {time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    main(sys.argv[1:] or ['cref', 'ratio'])
