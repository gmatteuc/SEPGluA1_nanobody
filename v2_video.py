"""
v2 videos, one per cohort and reading, in the layout of the P8 videos:
left, the hemisphere-averaged cohort mean; right, reliability t = mean / SEM
over the mice with tissue at that voxel; atlas outlines and acronyms on both;
one frame per 20 um plane, front to back.

Everything is drawn in the adult CCF, because that is where v2_cohort now
builds the cohort volumes -- the young brains having been carried there one by
one (v2_to_ccf), which is what lets the pooled young group mix P20 and P16.

  cohorts   young (P20 + P16 pooled), young_P20, adult, naive, rws
  readings  the same four the region tables carry: ratio (per unit
            autofluorescence), cref (relative to the brain's own isocortex),
            subref (relative to the subcortex without HPF and STR) and zref
            (range-matched; a signed position, not an intensity)

Colour range of the mean panel is fixed per reading so cohorts can be compared
by eye; the t panel runs to that cohort's 95th percentile. Grey = fewer than
MIN_N mice with tissue.

Output: data/comparisons_v2/ccf/<cohort>/video_<reading>_<cohort>.mp4

  D:\\sep_histology\\code\\tools\\venv_atlas\\Scripts\\python.exe v2_video.py [cohort ...]
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
from v2_cohort import OUT_ROOT as CCF_ROOT, COHORTS, MODES

MEAN_VMAX = {'cref': 2.0, 'ratio': 2.0, 'subref': 2.0, 'zref': 2.0}   # cortex sits near 1 in both readings; HPF saturates by design
T_PCT = 95.0                              # t panel range: 0 .. this percentile of t over the cohort's voxels
MIN_N = {'young': 2, 'young_P20': 2, 'young_P16': 1, 'naive': 3, 'rws': 3, 'adult': 5}
FPS = 12
MIN_LABEL_AREA = 150     # 20 um voxels in the plane, below which no acronym is drawn


def fold(v):
    h = v.shape[2] // 2
    return np.nanmean(np.stack([v[:, :, :h], v[:, :, v.shape[2] - h:][:, :, ::-1]]), axis=0)


def fold_n(n):
    h = n.shape[2] // 2
    return np.maximum(n[:, :, :h], n[:, :, n.shape[2] - h:][:, :, ::-1])


def boundaries(lab):
    b = np.zeros(lab.shape, bool)
    b[1:, :] |= lab[1:, :] != lab[:-1, :]
    b[:, 1:] |= lab[:, 1:] != lab[:, :-1]
    return b & (lab > 0)


def annotation_ccf20():
    """The full CCF annotation at 20 um, the grid the cohort volumes live on."""
    return np.asarray(nib.load(os.path.join(DATA, 'atlas', 'annotation_10.nii.gz')).dataobj)[::2, ::2, ::2]


def main(cohorts):
    acro = {}
    with open(CSV_MAP, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['parcellation_term_set_name'] == 'structure':
                acro[int(row['parcellation_index'])] = row['parcellation_term_acronym']
    hot = plt.get_cmap('hot')
    hot_cut = LinearSegmentedColormap.from_list('hot_cut', hot(np.linspace(0, 0.82, 256)))
    hot_cut.set_bad((0, 0, 0, 0))            # masked = transparent, so the grey/black ground shows through
    ann = annotation_ccf20()
    ann_h = ann[:, :, :ann.shape[2] // 2]
    for cohort in cohorts:
        n_h = fold_n(np.load(os.path.join(CCF_ROOT, cohort, 'cref_n.npy')))
        for reading in MODES:
            t0 = time.time()
            mean = fold(np.load(os.path.join(CCF_ROOT, cohort, f'{reading}_mean.npy')))
            sd = fold(np.load(os.path.join(CCF_ROOT, cohort, f'{reading}_sd.npy')))
            ok = (n_h >= MIN_N[cohort]) & np.isfinite(mean)
            with np.errstate(divide='ignore', invalid='ignore'):
                tval = np.where(ok & (sd > 0), mean / (sd / np.sqrt(np.maximum(n_h, 1))), np.nan)
            t_vmax = float(np.nanpercentile(tval[ok & (sd > 0)], T_PCT))
            frames = [k for k in range(ann_h.shape[0]) if ok[k].sum() > 200]
            out = os.path.join(CCF_ROOT, cohort, f'video_{reading}_{cohort}.mp4')
            writer = imageio_ffmpeg.write_frames(out, (1600, 800), fps=FPS, quality=7, macro_block_size=8)
            writer.send(None)
            fig = plt.figure(figsize=(16, 8), dpi=100, facecolor='k')
            axes = [fig.add_axes([0.04, 0.06, 0.40, 0.82]), fig.add_axes([0.53, 0.06, 0.40, 0.82])]
            caxes = [fig.add_axes([0.445, 0.12, 0.012, 0.70]), fig.add_axes([0.935, 0.12, 0.012, 0.70])]
            title = f'{cohort.replace("_", " ")} nano, v2 {reading} (n = {len(COHORTS[cohort])})'
            for k in frames:
                lab = ann_h[k]; inside = lab > 0
                bnd = boundaries(lab)
                panels = ((np.where(ok[k], mean[k], np.nan), hot_cut, (0, MEAN_VMAX[reading]), f'{title} - mean (hemispheres averaged)'),
                          (np.where(ok[k], tval[k], np.nan), hot_cut, (0, t_vmax), f'{title} - reliability t = mean/SEM  (range = {T_PCT:.0f}th pct)'))
                for ax, cax, (im, cmap, lim, ttl) in zip(axes, caxes, panels):
                    ax.clear(); cax.clear()
                    # black outside the atlas, grey inside it where there is no data, colour where there is
                    bg = np.zeros(lab.shape + (4,)); bg[inside] = (0.23, 0.23, 0.23, 1.0)
                    ax.imshow(bg, origin='upper', interpolation='nearest', aspect='equal')
                    im_full = np.ma.masked_invalid(np.where(inside, im, np.nan))
                    h = ax.imshow(im_full, cmap=cmap, vmin=lim[0], vmax=lim[1], origin='upper', interpolation='nearest', aspect='equal')
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
                fig.text(0.5, 0.93, f'CCF plane {2 * k} / 10 um   n = {int(np.nanmax(np.where(ok[k], n_h[k], 0)))} mice at this plane',
                         color='w', fontsize=12, ha='center')
                fig.canvas.draw()
                rgba = np.asarray(fig.canvas.buffer_rgba())
                writer.send(np.ascontiguousarray(rgba[:, :, :3]))
            writer.close(); plt.close(fig)
            print(f'{cohort:10s} {reading:6s} {len(frames)} frames -> {out}   {time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    main(sys.argv[1:] or ['young', 'adult', 'young_P20', 'naive', 'rws'])
