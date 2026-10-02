"""Cohort videos: the mean and its reliability t, plane by plane, per reading.

One video per cohort and reading, in the adult CCF, where volumes.cohort builds the
cohort volumes; the young brains were carried there one by one (volumes.to_ccf),
which is what lets the pooled young group mix ages. Each frame is one 20 um plane,
front to back: on the left the hemisphere-averaged cohort mean, on the right the
reliability t = mean / SEM over the mice with tissue at that voxel, with the atlas
outlines and acronyms on both.

    cohorts   young (every registered young brain), young_P20, adult, naive, rws
    readings  those of the region tables: ratio (per unit autofluorescence),
              sepratio (per unit SEP), cref (relative to the brain's own
              isocortex), subref (relative to the subcortex without HPF and STR)
              and zref (range-matched; a signed position, not an intensity)

The colour range of the mean is fixed per reading, so that cohorts can be compared
by eye; the t panel runs to that cohort's 95th percentile of t. Grey is fewer than
MIN_N mice with tissue.

Writes comparisons_v2/ccf/<cohort>/video_<reading>_<cohort>.mp4.

Run by run_video.py.
"""

import csv
import os
import time

import imageio_ffmpeg
import matplotlib
import nibabel as nib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import center_of_mass

from sepmap.config import SETTINGS
from sepmap.volumes.cohort import COHORTS, MODES, SIGNED_READINGS
from sepmap.volumes.cohort import OUT_ROOT as CCF_ROOT
from sepmap.volumes.per_mouse import CSV_MAP, DATA

# colour range of the mean per reading, set so that cortex sits near half of it: cref
# is 1 by construction, adult cortex is 1.01 in ratio and 0.29 in sepratio; the
# hippocampus saturates by design
MEAN_VMAX = {"ratio": 2.0, "sepratio": 0.6, "cref": 2.0, "subref": 2.0, "zref": 2.0}

# range of the t panel: 0 to this percentile of t over the cohort's voxels
T_PCT = 95.0

# brains with tissue that a voxel needs, per cohort
MIN_N = {"young": 2, "young_P20": 2, "young_P16": 1, "naive": 3, "rws": 3, "adult": 5}

# frames per second, and the 20 um voxels a structure needs in the plane to get its
# acronym drawn
FPS = SETTINGS["videos"]["fps"]
MIN_LABEL_AREA = SETTINGS["videos"]["min_label_area"]


def fold(v):
    """Average the two hemispheres of an (AP, DV, ML) volume, ignoring NaN."""
    h = v.shape[2] // 2
    return np.nanmean(
        np.stack([v[:, :, :h], v[:, :, v.shape[2] - h :][:, :, ::-1]]), axis=0
    )


def fold_count(n):
    """Fold a count map as fold does, keeping the larger count of the two sides."""
    h = n.shape[2] // 2
    return np.maximum(n[:, :, :h], n[:, :, n.shape[2] - h :][:, :, ::-1])


def boundaries(lab):
    """Pixels of the label image `lab` that border another label, inside the atlas."""
    b = np.zeros(lab.shape, bool)
    b[1:, :] |= lab[1:, :] != lab[:-1, :]
    b[:, 1:] |= lab[:, 1:] != lab[:, :-1]
    return b & (lab > 0)


def annotation_ccf20():
    """The full CCF annotation at 20 um, the grid the cohort volumes live on."""
    return np.asarray(
        nib.load(os.path.join(DATA, "atlas", "annotation_10.nii.gz")).dataobj
    )[::2, ::2, ::2]


def main(cohorts):
    """Write the videos of each of `cohorts`, one per reading in force."""
    # acronyms by parcellation index
    acro = {}
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["parcellation_term_set_name"] == "structure":
                acro[int(row["parcellation_index"])] = row["parcellation_term_acronym"]

    # hot up to 0.82 of its range, and purple-orange for zref, a position rather than
    # an intensity; masked voxels transparent, so the grey or black ground shows
    hot = plt.get_cmap("hot")
    hot_cut = LinearSegmentedColormap.from_list("hot_cut", hot(np.linspace(0, 0.82, 256)))
    hot_cut.set_bad((0, 0, 0, 0))
    puor = plt.get_cmap("PuOr_r").copy()
    puor.set_bad((0, 0, 0, 0))

    # the annotation, the half that the folded volumes cover
    ann = annotation_ccf20()
    ann_h = ann[:, :, : ann.shape[2] // 2]
    for cohort in cohorts:
        n_h = fold_count(np.load(os.path.join(CCF_ROOT, cohort, "cref_n.npy")))
        for reading in MODES:
            t0 = time.time()

            # folded mean and SD, and t where there are enough brains and an SD
            mean = fold(np.load(os.path.join(CCF_ROOT, cohort, f"{reading}_mean.npy")))
            sd = fold(np.load(os.path.join(CCF_ROOT, cohort, f"{reading}_sd.npy")))
            signed = reading in SIGNED_READINGS
            ok = (n_h >= MIN_N[cohort]) & np.isfinite(mean)
            with np.errstate(divide="ignore", invalid="ignore"):
                tval = np.where(
                    ok & (sd > 0), mean / (sd / np.sqrt(np.maximum(n_h, 1))), np.nan
                )

            # colour limits: fixed for the mean; for the t panel the T_PCT percentile
            # of t, of |t| for a signed reading
            t_vmax = float(
                np.nanpercentile(
                    np.abs(tval[ok & (sd > 0)]) if signed else tval[ok & (sd > 0)], T_PCT
                )
            )
            cmap_use = puor if signed else hot_cut
            if signed:
                lim_mean = (-MEAN_VMAX[reading], MEAN_VMAX[reading])
            else:
                lim_mean = (0, MEAN_VMAX[reading])
            lim_t = (-t_vmax, t_vmax) if signed else (0, t_vmax)

            # a video of the planes with more than 200 voxels with data, drawn into
            # one 1600 x 800 figure plane by plane
            frames = [k for k in range(ann_h.shape[0]) if ok[k].sum() > 200]
            out = os.path.join(CCF_ROOT, cohort, f"video_{reading}_{cohort}.mp4")
            writer = imageio_ffmpeg.write_frames(
                out, (1600, 800), fps=FPS, quality=7, macro_block_size=8
            )
            writer.send(None)
            fig = plt.figure(figsize=(16, 8), dpi=100, facecolor="k")
            axes = [
                fig.add_axes([0.04, 0.06, 0.40, 0.82]),
                fig.add_axes([0.53, 0.06, 0.40, 0.82]),
            ]
            caxes = [
                fig.add_axes([0.445, 0.12, 0.012, 0.70]),
                fig.add_axes([0.935, 0.12, 0.012, 0.70]),
            ]
            title = (
                f"{cohort.replace('_', ' ')} nano, v2 {reading} "
                f"(n = {len(COHORTS[cohort])})"
            )
            for k in frames:
                lab = ann_h[k]
                inside = lab > 0
                bnd = boundaries(lab)
                panels = (
                    (
                        np.where(ok[k], mean[k], np.nan),
                        cmap_use,
                        lim_mean,
                        f"{title} - mean (hemispheres averaged)",
                    ),
                    (
                        np.where(ok[k], tval[k], np.nan),
                        cmap_use,
                        lim_t,
                        f"{title} - reliability t = mean/SEM  "
                        f"(range = {T_PCT:.0f}th pct)",
                    ),
                )
                for ax, cax, (im, cmap, lim, ttl) in zip(axes, caxes, panels):
                    ax.clear()
                    cax.clear()

                    # black outside the atlas, grey inside it where there is no data,
                    # colour where there is
                    bg = np.zeros(lab.shape + (4,))
                    bg[inside] = (0.23, 0.23, 0.23, 1.0)
                    ax.imshow(bg, origin="upper", interpolation="nearest", aspect="equal")
                    im_full = np.ma.masked_invalid(np.where(inside, im, np.nan))
                    h = ax.imshow(
                        im_full,
                        cmap=cmap,
                        vmin=lim[0],
                        vmax=lim[1],
                        origin="upper",
                        interpolation="nearest",
                        aspect="equal",
                    )

                    # area borders, and the acronyms of the structures large enough
                    ov = np.zeros(lab.shape + (4,))
                    ov[bnd] = (0.75, 0.75, 0.75, 0.9)
                    ax.imshow(ov, origin="upper", interpolation="nearest", aspect="equal")
                    for idx in np.unique(lab):
                        if idx == 0:
                            continue
                        m = lab == idx
                        if m.sum() < MIN_LABEL_AREA:
                            continue
                        cy, cx = center_of_mass(m)
                        ax.text(
                            cx,
                            cy,
                            acro.get(int(idx), ""),
                            color="w",
                            fontsize=5.5,
                            ha="center",
                            va="center",
                        )

                    # a black panel without ticks or frame, white title and labels
                    ax.set_facecolor("k")
                    ax.set_xticks([])
                    ax.set_yticks([])
                    for s in ax.spines.values():
                        s.set_visible(False)
                    ax.set_title(ttl, color="w", fontsize=12)
                    cb = fig.colorbar(h, cax=cax)
                    cb.ax.yaxis.set_tick_params(color="w", labelcolor="w")

                # header: the plane, and the most brains behind any voxel of it
                fig.texts.clear()
                fig.text(
                    0.5,
                    0.93,
                    f"CCF plane {2 * k} / 10 um   "
                    f"n = {int(np.nanmax(np.where(ok[k], n_h[k], 0)))} "
                    "mice at this plane",
                    color="w",
                    fontsize=12,
                    ha="center",
                )
                fig.canvas.draw()
                rgba = np.asarray(fig.canvas.buffer_rgba())
                writer.send(np.ascontiguousarray(rgba[:, :, :3]))
            writer.close()
            plt.close(fig)
            print(
                f"{cohort:10s} {reading:6s} {len(frames)} frames -> {out}   "
                f"{time.time() - t0:.0f} s",
                flush=True,
            )
