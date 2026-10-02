"""Young and adult side by side, plane by plane, in the adult CCF.

Three panels per frame: the young cohort mean, the adult cohort mean on the same
colour scale, and their comparison, so a pattern can be compared where it sits
rather than by flicking between two windows. Atlas outlines and acronyms on all
three; the header gives the CCF plane and how many brains are behind each side
there.

Everything comes from the cohort volumes volumes.cohort built in the CCF, the ones
young_vs_adult.compare draws its figures from, and every reading is available. The
colour scale of the two means is fixed and shared. For ratio, sepratio, cref and
subref it is linear and the third panel is a log2 ratio of the means, each floored
at 0.02; for zref the means are a signed position within each brain's range and
the third panel is their difference.

One plane can be drawn on its own instead, as a still for a talk or a message, and
the colour range can be tightened for that run only: the cortex occupies a
fraction of the range the hippocampus needs, so saturating the extremes is often
the only way to see it. A tightened range of the means goes into the file name.

Writes video_side_by_side_<reading>.mp4, or plane<P>_side_by_side_<reading>.png
for one plane, into comparisons_v2/young_vs_adult/.

Run by run_video_compare.py.
"""

import csv
import os
import time

import numpy as np
import nibabel as nib
import imageio_ffmpeg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import center_of_mass

from sepmap.config import SETTINGS
from sepmap.volumes.per_mouse import CSV_MAP, DATA
from sepmap.volumes.cohort import COHORTS, OUT_ROOT as CCF_ROOT, SIGNED_READINGS
from sepmap.young_vs_adult.compare import MIN_N_YOUNG, MIN_N_ADULT, YOUNG, fold, fold_n

OUT = os.path.join(DATA, "comparisons_v2", "young_vs_adult")

# colour range of the two means per reading, fixed and shared
MEAN_VMAX = {"ratio": 2.0, "sepratio": 0.6, "cref": 2.0, "subref": 2.0, "zref": 2.0}

# colour range of the comparison, symmetric about zero
LOG2_LIM = 1.5

# frames per second, and the 20 um voxels a structure needs in the plane to get its
# acronym drawn
FPS = SETTINGS["videos"]["fps"]
MIN_LABEL_AREA = SETTINGS["videos"]["min_label_area"]


def boundaries(lab):
    """Pixels of the label image `lab` that border another label, inside the atlas."""
    b = np.zeros(lab.shape, bool)
    b[1:, :] |= lab[1:, :] != lab[:-1, :]
    b[:, 1:] |= lab[:, 1:] != lab[:, :-1]
    return b & (lab > 0)


def main(readings, plane=None, vmax=None, dlim=None):
    """Young beside adult for each of `readings`: a video, or a still of CCF `plane`.

    `plane` is numbered at 10 um; `vmax` and `dlim` replace MEAN_VMAX and LOG2_LIM
    for this run.
    """
    # acronyms by parcellation index
    acro = {}
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["parcellation_term_set_name"] == "structure":
                acro[int(row["parcellation_index"])] = row["parcellation_term_acronym"]

    # hot up to 0.82 of its range, and red-blue for the comparison; masked voxels
    # transparent, so the grey or black ground shows
    hot = plt.get_cmap("hot")
    hot_cut = LinearSegmentedColormap.from_list("hot_cut", hot(np.linspace(0, 0.82, 256)))
    hot_cut.set_bad((0, 0, 0, 0))
    rdbu = plt.get_cmap("RdBu_r").copy()
    rdbu.set_bad((0, 0, 0, 0))

    # purple-orange for zref, which is a position rather than an intensity, so
    # red-blue means one thing only: a young-minus-adult difference
    puor = plt.get_cmap("PuOr_r").copy()
    puor.set_bad((0, 0, 0, 0))

    # the annotation, the half that the folded volumes cover, and the brains with
    # tissue behind each voxel
    ann = np.asarray(
        nib.load(os.path.join(DATA, "atlas", "annotation_10.nii.gz")).dataobj
    )[::2, ::2, ::2]
    ann_h = ann[:, :, : ann.shape[2] // 2]
    y_n = fold_n(np.load(os.path.join(CCF_ROOT, YOUNG, "cref_n.npy")))
    a_n = fold_n(np.load(os.path.join(CCF_ROOT, "adult", "cref_n.npy")))

    for reading in readings:
        t0 = time.time()

        # folded means, the voxels where both have enough brains, and the comparison
        y = fold(np.load(os.path.join(CCF_ROOT, YOUNG, f"{reading}_mean.npy")))
        a = fold(np.load(os.path.join(CCF_ROOT, "adult", f"{reading}_mean.npy")))
        ok_y = (y_n >= MIN_N_YOUNG) & np.isfinite(y)
        ok_a = (a_n >= MIN_N_ADULT) & np.isfinite(a)
        both = ok_y & ok_a & (ann_h > 0)
        signed = reading in SIGNED_READINGS
        log2 = np.where(
            both,
            (y - a) if signed else np.log2(np.maximum(y, 0.02) / np.maximum(a, 0.02)),
            np.nan,
        )

        # colour ranges of this run
        v_mean = MEAN_VMAX[reading] if vmax is None else vmax
        v_diff = LOG2_LIM if dlim is None else dlim

        # a still of one plane, or a video of the planes with more than 200 compared
        # voxels; a CCF plane is quoted at 10 um in every caption, the volumes are 20 um
        frames = (
            [plane // 2]
            if plane is not None
            else [k for k in range(ann_h.shape[0]) if both[k].sum() > 200]
        )
        tag = "" if vmax is None else f"_vmax{v_mean:g}"
        if plane is None:
            out = os.path.join(OUT, f"video_side_by_side_{reading}{tag}.mp4")
            writer = imageio_ffmpeg.write_frames(
                out, (1920, 760), fps=FPS, quality=7, macro_block_size=8
            )
            writer.send(None)
        else:
            out = os.path.join(OUT, f"plane{plane}_side_by_side_{reading}{tag}.png")
            writer = None

        # one 1920 x 760 figure, three panels with their colour bars
        fig = plt.figure(figsize=(19.2, 7.6), dpi=100, facecolor="k")
        axes = [fig.add_axes([0.02 + i * 0.325, 0.05, 0.27, 0.82]) for i in range(3)]
        caxes = [fig.add_axes([0.295 + i * 0.325, 0.12, 0.009, 0.68]) for i in range(3)]
        for k in frames:
            lab = ann_h[k]
            inside = lab > 0
            bnd = boundaries(lab)
            cmap_mean = puor if signed else hot_cut
            lim_mean = (-v_mean, v_mean) if signed else (0, v_mean)
            panels = (
                (
                    np.where(ok_y[k], y[k], np.nan),
                    cmap_mean,
                    lim_mean,
                    f"young (n = {len(COHORTS[YOUNG])})   {reading}",
                ),
                (
                    np.where(ok_a[k], a[k], np.nan),
                    cmap_mean,
                    lim_mean,
                    f"adult (n = {len(COHORTS['adult'])})   {reading}",
                ),
                (
                    log2[k],
                    rdbu,
                    (-v_diff, v_diff),
                    "young - adult" if signed else "log2( young / adult )",
                ),
            )
            for ax, cax, (im, cmap, lim, ttl) in zip(axes, caxes, panels):
                ax.clear()
                cax.clear()

                # the atlas in dark grey under the data, black outside it
                bg = np.zeros(lab.shape + (4,))
                bg[inside] = (0.23, 0.23, 0.23, 1.0)
                ax.imshow(bg, origin="upper", interpolation="nearest", aspect="equal")
                h = ax.imshow(
                    np.ma.masked_invalid(np.where(inside, im, np.nan)),
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

            # header: the plane, the most brains behind any voxel on each side, and
            # what the panels show
            fig.texts.clear()
            fig.text(
                0.5,
                0.93,
                f"CCF plane {2 * k} / 10 um    "
                f"young: {int(np.nanmax(np.where(ok_y[k], y_n[k], 0)))} brains   "
                f"adult: {int(np.nanmax(np.where(ok_a[k], a_n[k], 0)))} brains   "
                + (
                    "(means are a position within each brain's own range; "
                    "right panel is their difference)"
                    if signed
                    else "(means on a linear scale, shared range; "
                    "only the right panel is log2)"
                ),
                color="w",
                fontsize=12,
                ha="center",
            )

            # into the video, or the still
            fig.canvas.draw()
            if writer is not None:
                writer.send(
                    np.ascontiguousarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3])
                )
            else:
                fig.savefig(out, dpi=100, facecolor="k")
        if writer is not None:
            writer.close()
        plt.close(fig)
        print(
            f"{reading:8s} {len(frames)} frame(s) -> {out}   {time.time() - t0:.0f} s",
            flush=True,
        )
