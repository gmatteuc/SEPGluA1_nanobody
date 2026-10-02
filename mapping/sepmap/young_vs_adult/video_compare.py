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

import os
import time

import imageio_ffmpeg
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

from sepmap.config import SETTINGS
from sepmap.plotting import coronal_figure, coronal_frame, hot_cut, transparent_bad
from sepmap.volumes.cohort import COHORTS, SIGNED_READINGS
from sepmap.volumes.cohort import OUT_ROOT as CCF_ROOT
from sepmap.volumes.per_mouse import DATA, structure_terms
from sepmap.young_vs_adult.compare import MIN_N_ADULT, MIN_N_YOUNG, YOUNG
from sepmap.young_vs_adult.hemispheres import fold, fold_count

OUT = os.path.join(DATA, "comparisons_v2", "young_vs_adult")

# colour range of the two means per reading, fixed and shared
MEAN_VMAX = {"ratio": 2.0, "sepratio": 0.6, "cref": 2.0, "subref": 2.0, "zref": 2.0}

# colour range of the comparison, symmetric about zero
LOG2_LIM = 1.5

# frames per second
FPS = SETTINGS["videos"]["fps"]


def main(
    readings: list[str],
    plane: int | None = None,
    vmax: float | None = None,
    dlim: float | None = None,
) -> None:
    """Young beside adult for each of `readings`: a video, or a still of CCF `plane`.

    `plane` is numbered at 10 um; `vmax` and `dlim` replace MEAN_VMAX and LOG2_LIM
    for this run.
    """
    # acronyms by parcellation index
    _, acro, _ = structure_terms()

    # hot up to 0.82 of its range, and red-blue for the comparison; masked voxels
    # transparent, so the grey or black ground shows
    hot = hot_cut()
    rdbu = transparent_bad("RdBu_r")

    # purple-orange for zref, which is a position rather than an intensity, so
    # red-blue means one thing only: the young-adult difference
    puor = transparent_bad("PuOr_r")

    # the annotation, the half that the folded volumes cover, and the brains with
    # tissue behind each voxel
    ann = np.asarray(
        nib.load(os.path.join(DATA, "atlas", "annotation_10.nii.gz")).dataobj
    )[::2, ::2, ::2]
    ann_h = ann[:, :, : ann.shape[2] // 2]
    y_n = fold_count(np.load(os.path.join(CCF_ROOT, YOUNG, "cref_n.npy")))
    a_n = fold_count(np.load(os.path.join(CCF_ROOT, "adult", "cref_n.npy")))

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
        if plane is not None:
            frames = [plane // 2]
        else:
            frames = [k for k in range(ann_h.shape[0]) if both[k].sum() > 200]
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
        fig, axes, caxes = coronal_figure()
        for k in frames:
            cmap_mean = puor if signed else hot
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

            # header: the plane, the most brains behind any voxel on each side, and
            # what the panels show
            if signed:
                note = (
                    "(means are a position within each brain's own range; "
                    "right panel is their difference)"
                )
            else:
                note = (
                    "(means on a linear scale, shared range; "
                    "only the right panel is log2)"
                )
            header = (
                f"CCF plane {2 * k} / 10 um    "
                f"young: {int(np.nanmax(np.where(ok_y[k], y_n[k], 0)))} brains   "
                f"adult: {int(np.nanmax(np.where(ok_a[k], a_n[k], 0)))} brains   " + note
            )

            # the atlas in dark grey under the data, black outside it, the area
            # borders and the acronyms of the structures large enough
            coronal_frame(fig, axes, caxes, k, panels, ann_h, acro, header)

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
