"""One reading looked at closely: a coronal plane, its video and the cortical flatmaps.

The cohort volumes are prepared once and every output is drawn from them, so a
coronal frame, the video it comes from and the flatmaps cannot disagree about
contrast, coverage or smoothing:

    detail_plane<P>_<reading>.png        one coronal plane: young, adult, difference
    detail_video_<reading>.mp4           the same three panels, plane by plane
    detail_flatmap_<reading>.png         the isocortex unrolled, full cortical depth
    detail_flatmap_layers_<reading>.png  the same by depth band: supragranular,
                                         granular, infragranular

The standard videos and slice figures show every reading at one fixed scale, so
that cohorts stay comparable. These look closely at one reading: the colour range
is tightened until cortex is readable (the hippocampus then clips, on purpose),
and the volumes are smoothed lightly in 3D before anything is drawn.

The smoothing is cosmetic and every title declares it. It is mask normalised, so
tissue is never averaged with the black outside it, and it never reaches the
statistics: the region tables and tests read the unsmoothed per-mouse volumes. It
removes sampling, not signal. Sections sit about 150 um apart and the registered
volumes interpolate between them, which leaves a faint coronal banding; the
flatmap samples one streamline at a time, which leaves fine radial streaks. The
default sigma, 3 x 1 x 1 voxels, is 60 um along AP and 20 um across: anisotropic
because the banding is, and far below the size of any area.

A colormap named with --cmap draws the mean panels and sends the whole set into a
subfolder of that name, so the default set stays. The difference panel stays
red-blue whatever the colormap: it is a signed quantity read against zero, and a
rainbow has no neutral middle to read zero against.

The flatmaps need ccf_streamlines, which brings its own numpy and scikit-image, so
this module runs in its own environment, tools\\venv_flat (made from
tools\\requirements_flat.txt), and imports only config from the package. Its
assets, about 0.6 GB, are fetched once into atlas_flatmap/ under the data root
from the Allen Institute's ccf_streamlines_assets folder,
    https://download.alleninstitute.org/informatics-archive/current-release/
    mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/

    streamlines/surface_paths_10_v3.h5                where each streamline runs
    view_lookup/flatmap_butterfly.h5                  the streamline of each pixel
    master_updated/flatmap_butterfly.nrrd             the areas, for their borders
    master_updated/labelDescription_ITKSNAPColor.txt  and their names
    cortical_metrics/avg_layer_depths.json            where the layers sit
    cortical_metrics/cortical_layers_10_v2.h5         the same per streamline

surface_paths_10_v3.h5 alone is 0.5 GB.

Run by run_closeup.py.
"""

import csv
import json
import os
import time

import imageio_ffmpeg
import matplotlib
import nibabel as nib
import numpy as np

matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import center_of_mass, gaussian_filter

from sepmap.config import DATA, SETTINGS

ASSETS = os.path.join(DATA, "atlas_flatmap")
CCF_ROOT = os.path.join(DATA, "comparisons_v2", "ccf")
OUT = os.path.join(DATA, "comparisons_v2", "young_vs_adult")
CSV_MAP = os.path.join(DATA, "atlas", "parcellation_to_parcellation_term_membership.csv")

YOUNG, ADULT = "young", "adult"

# this module runs in the flatmap environment and cannot import the analysis chain,
# so what it shares with the chain comes from settings.toml, which both read

# brains with tissue that a voxel needs, as in young_vs_adult.compare, so the maps
# agree
MIN_N_YOUNG = SETTINGS["young_vs_adult"]["min_n_young"]
MIN_N_ADULT = SETTINGS["young_vs_adult"]["min_n_adult"]

# floor that keeps log2 finite where a reading is near zero, the value volumes.cohort
# uses; only the readings compared as a ratio meet it, never zref
LOG2_FLOOR = SETTINGS["readings"]["log2_floor"]

# readings that are a position within a brain's own range, as volumes.cohort's
# SIGNED_READINGS: compared by difference, drawn diverging from -v to +v
SIGNED_READINGS = tuple(SETTINGS["readings"]["signed"])

# colour ranges of the means and of the difference, per reading, tightened for cortex
VMAX = {"zref": 0.9, "cref": 2.0, "subref": 2.0, "ratio": 2.0, "sepratio": 0.6}
DLIM = {"zref": 0.5, "cref": 1.0, "subref": 1.0, "ratio": 1.0, "sepratio": 1.0}

# CCF plane of the still, at 10 um, where RL and AL are cut
PLANE = 790

# Gaussian sigma in 20 um voxels along (AP, DV, ML), anisotropic because the artefact
# is: sections sit about 150 um (7.5 voxels) apart and the registered volumes
# interpolate between them, so the banding runs along AP only. 3 voxels of AP blur
# is 60 um, well under a section spacing and two orders below the distance from V1
# to RL.
SMOOTH = (3.0, 1.0, 1.0)

# frames per second of the video, and the 20 um voxels a structure needs in a
# coronal plane to get its acronym drawn
FPS = SETTINGS["videos"]["fps"]
MIN_LABEL_AREA = SETTINGS["videos"]["min_label_area"]

# the areas the argument is about, outlined brighter on the flatmaps
HIGHLIGHT = ("VISp", "VISrl", "VISal", "SSp-bfd")

# the areas named on the flatmaps
LABEL_AREAS = [
    "VISp",
    "VISrl",
    "VISal",
    "VISam",
    "SSp-bfd",
    "SSs",
    "AUDp",
    "MOp",
    "RSPd",
    "ACAd",
    "TEa",
]

# RL and AL are neighbours whose centroids are 100 flatmap pixels apart, which is
# not enough at three panels wide: their labels move apart, in flatmap pixels
LABEL_NUDGE = {"VISrl": (60, 0), "VISal": (-60, 0)}

# depth bands of the layer flatmap: the title, and the layers each pools
BANDS = [
    ("supragranular  L1 + L2/3", ("Isocortex layer 1", "Isocortex layer 2/3")),
    ("granular  L4", ("Isocortex layer 4",)),
    (
        "infragranular  L5 + L6",
        ("Isocortex layer 5", "Isocortex layer 6a", "Isocortex layer 6b"),
    ),
]


def save_figure(fig, path, dpi):
    """Save `fig` as a PNG at `path` on black, and as an EPS beside it.

    Windows refuses to overwrite a PNG that an image viewer holds open, and these
    figures are made to be looked at while the next one is drawn. The figure then
    goes to <name>_new.png, with a note.

    The EPS is what goes into a figure for a paper. PostScript has no
    transparency, so the image layers are rasterised and composited by Agg
    first; otherwise a no-data region, transparent here, would come out opaque
    black instead of showing the ground beneath it. Text, lines and axes stay
    vector, the part that has to be editable.
    """
    try:
        fig.savefig(path, dpi=dpi, facecolor="k")
    except OSError:
        alt = path.replace(".png", "_new.png")
        fig.savefig(alt, dpi=dpi, facecolor="k")
        print(
            f"  NOTE: {os.path.basename(path)} is open elsewhere; "
            f"wrote {os.path.basename(alt)}",
            flush=True,
        )

    # the EPS, with the image layers rasterised
    eps = os.path.splitext(path)[0] + ".eps"
    for ax in fig.axes:
        for im in ax.images:
            im.set_rasterized(True)
    try:
        fig.savefig(eps, dpi=dpi, facecolor=fig.get_facecolor(), format="eps")
    except OSError:
        print(
            f"  NOTE: {os.path.basename(eps)} is open elsewhere; "
            "the PNG was still written",
            flush=True,
        )


def cohort_size(cohort):
    """Number of brains in `cohort`, counted in the mice list volumes.cohort wrote.

    Counted rather than typed, so a caption cannot go stale when a brain is added.
    """
    with open(os.path.join(CCF_ROOT, cohort, "mice.txt"), encoding="utf-8") as fh:
        return sum(1 for line in fh if line.strip())


def fold(v):
    """Average the two hemispheres of an (AP, DV, ML) volume, ignoring NaN."""
    h = v.shape[2] // 2
    return np.nanmean(
        np.stack([v[:, :, :h], v[:, :, v.shape[2] - h :][:, :, ::-1]]), axis=0
    )


def fold_n(n):
    """Fold a count map as fold does, keeping the larger count of the two sides."""
    h = n.shape[2] // 2
    return np.maximum(n[:, :, :h], n[:, :, n.shape[2] - h :][:, :, ::-1])


def prepare(reading, sigma):
    """Value and mask of each cohort for `reading`, hemispheres folded, at 20 um.

    Returns {cohort: (value, mask)}, the one source of all four figures. The mask
    holds the voxels with enough brains and a value; the value is 0 outside it and,
    with `sigma`, smoothed inside it. The mask itself is not smoothed: blurring a
    value inside the tissue is cosmetic, growing the tissue would not be.
    """
    out = {}
    for cohort, min_n in ((YOUNG, MIN_N_YOUNG), (ADULT, MIN_N_ADULT)):
        mean = fold(np.load(os.path.join(CCF_ROOT, cohort, f"{reading}_mean.npy")))
        n = fold_n(np.load(os.path.join(CCF_ROOT, cohort, f"{reading}_n.npy")))
        m = ((n >= min_n) & np.isfinite(mean)).astype(np.float32)
        v = np.where(m > 0, mean, 0).astype(np.float32)

        # mask-normalised smoothing: tissue is never averaged with the zeros around it
        if np.any(sigma):
            num = gaussian_filter(v, sigma)
            den = gaussian_filter(m, sigma)
            v = np.where(m > 0, num / np.maximum(den, 1e-6), 0).astype(np.float32)
        out[cohort] = (v, m)
    return out


def difference(vals, signed):
    """Young against adult where both cohorts have a value (NaN elsewhere), and that mask.

    A signed reading is compared by difference, the others by log2 ratio, each value
    floored at LOG2_FLOOR.
    """
    y, a = vals[YOUNG], vals[ADULT]
    both = (y[1] > 0) & (a[1] > 0)
    d = (
        (y[0] - a[0])
        if signed
        else np.log2(np.maximum(y[0], LOG2_FLOOR) / np.maximum(a[0], LOG2_FLOOR))
    )
    return np.where(both, d, np.nan), both


def shown(value, mask):
    """`value` where `mask` is set, NaN elsewhere."""
    return np.where(mask > 0, value, np.nan)


# ===== Coronal plane and video =====


def annotation_half():
    """The CCF annotation at 20 um, the half that the folded volumes cover."""
    ann = np.asarray(
        nib.load(os.path.join(DATA, "atlas", "annotation_10.nii.gz")).dataobj
    )[::2, ::2, ::2]
    return ann[:, :, : ann.shape[2] // 2]


def boundaries(lab):
    """Pixels of the label image `lab` that border another label, inside the atlas."""
    b = np.zeros(lab.shape, bool)
    b[1:, :] |= lab[1:, :] != lab[:-1, :]
    b[:, 1:] |= lab[:, 1:] != lab[:, :-1]
    return b & (lab > 0)


def coronal_frame(fig, axes, caxes, k, panels, ann_h, acro, header, vector_outline=False):
    """Draw plane `k` into the three panels of `fig`, and the header above them.

    `panels` holds (image, colormap, limits, title) per panel. The atlas is dark
    grey under the data and its borders lie on top: as lines with
    `vector_outline`, otherwise as a pixel overlay. Structures with at least
    MIN_LABEL_AREA voxels in the plane get their acronym.
    """
    lab = ann_h[k]
    inside = lab > 0
    bnd = boundaries(lab)
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

        # area borders
        if vector_outline:
            # lines keep the atlas an editable layer of its own in the EPS; at about a
            # second per panel they suit a still, not the hundreds of frames of a video
            ax.contour(bnd.astype(float), levels=[0.5], colors="#bfbfbf", linewidths=0.3)
        else:
            ov = np.zeros(lab.shape + (4,))
            ov[bnd] = (0.75, 0.75, 0.75, 0.9)
            ax.imshow(ov, origin="upper", interpolation="nearest", aspect="equal")

        # acronyms of the structures large enough to name
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

        # a black panel without ticks or frame, white title and colour bar labels
        ax.set_facecolor("k")
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(ttl, color="w", fontsize=12)
        cb = fig.colorbar(h, cax=cax)
        cb.ax.yaxis.set_tick_params(color="w", labelcolor="w")
    fig.texts.clear()
    fig.text(0.5, 0.93, header, color="w", fontsize=12, ha="center")


def coronal(
    reading,
    vals,
    diff,
    both,
    plane,
    lim_mean,
    cmaps,
    sigma_txt,
    want_video,
    n,
    signed,
    out_dir,
):
    """Draw the still at CCF `plane` and, with `want_video`, the video of every plane.

    Young, adult and their difference side by side. `lim_mean` holds the limits of
    the means and of the difference, `cmaps` their colormaps, `n` the brains of
    each cohort. The video runs over the planes with more than 200 compared voxels.
    """
    ann_h = annotation_half()
    acro = {}
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["parcellation_term_set_name"] == "structure":
                acro[int(row["parcellation_index"])] = row["parcellation_term_acronym"]
    cmap_mean, rdbu = cmaps

    # one 1920 x 760 figure, three panels with their colour bars
    fig = plt.figure(figsize=(19.2, 7.6), dpi=100, facecolor="k")
    axes = [fig.add_axes([0.02 + i * 0.325, 0.05, 0.27, 0.82]) for i in range(3)]
    caxes = [fig.add_axes([0.295 + i * 0.325, 0.12, 0.009, 0.68]) for i in range(3)]

    def panels_at(k):
        """The three panels of plane `k`: young, adult and their difference."""
        return (
            (
                shown(vals[YOUNG][0][k], vals[YOUNG][1][k]),
                cmap_mean,
                lim_mean[0],
                f"young (n = {n[YOUNG]})   {reading}",
            ),
            (
                shown(vals[ADULT][0][k], vals[ADULT][1][k]),
                cmap_mean,
                lim_mean[0],
                f"adult (n = {n[ADULT]})   {reading}",
            ),
            (
                diff[k],
                rdbu,
                lim_mean[1],
                "young - adult" if signed else "log2( young / adult )",
            ),
        )

    def header_at(k):
        """Header of plane `k`, numbered on the 10 um CCF grid."""
        return (
            f"CCF plane {2 * k} / 10 um    young {n[YOUNG]} brains, "
            f"adult {n[ADULT]} brains    "
            f"{sigma_txt}, colour range tightened for cortex"
        )

    # the still; the volumes are at 20 um, the plane is numbered at 10 um
    k = plane // 2
    coronal_frame(
        fig, axes, caxes, k, panels_at(k), ann_h, acro, header_at(k), vector_outline=True
    )
    out = os.path.join(out_dir, f"detail_plane{plane}_{reading}.png")
    save_figure(fig, out, 100)
    print(f"  wrote {os.path.basename(out)}", flush=True)

    # the video, redrawing the same figure plane by plane
    if want_video:
        t0 = time.time()
        frames = [i for i in range(ann_h.shape[0]) if both[i].sum() > 200]
        out = os.path.join(out_dir, f"detail_video_{reading}.mp4")
        writer = imageio_ffmpeg.write_frames(
            out, (1920, 760), fps=FPS, quality=7, macro_block_size=8
        )
        writer.send(None)
        for i in frames:
            coronal_frame(fig, axes, caxes, i, panels_at(i), ann_h, acro, header_at(i))
            fig.canvas.draw()
            writer.send(
                np.ascontiguousarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3])
            )
        writer.close()
        print(
            f"  wrote {os.path.basename(out)}, {len(frames)} frames, "
            f"{time.time() - t0:.0f} s",
            flush=True,
        )
    plt.close(fig)


# ===== Flatmaps =====


def layer_thicknesses():
    """Thickness of each cortical layer in um, from the Allen's average depths.

    The file gives the depth of each layer's lower border below the pia, so the
    thicknesses are the differences between them.
    """
    d = json.load(open(os.path.join(ASSETS, "avg_layer_depths.json")))
    names = [
        "Isocortex layer 1",
        "Isocortex layer 2/3",
        "Isocortex layer 4",
        "Isocortex layer 5",
        "Isocortex layer 6a",
        "Isocortex layer 6b",
    ]
    thick, prev = {}, 0.0
    for name, low in zip(names, [d["2/3"], d["4"], d["5"], d["6a"], d["6b"], d["wm"]]):
        thick[name] = low - prev
        prev = low
    return thick


def band_edges(p3, slab_depth):
    """First and last depth bin of each layer, as the projector built them.

    The edges come from the projector, never recomputed: the slab has one bin per
    sample along a streamline (200 of them), not one per 10 um, and the layers get
    bins in proportion to their thickness. Edges derived here from the layer depths
    would put 96 bins where there are 200, and label the middle of layer 2/3 as
    layer 4 and layer 4 as layer 6. Stops if the bins do not add up to the slab.
    """
    bins = p3.reference_layer_thicknesses_in_voxels()
    edges, at = {}, 0
    for name in p3.ISOCORTEX_LAYER_KEYS:
        edges[name] = (at, at + bins[name])
        at += bins[name]
    if at != slab_depth:
        raise SystemExit(
            f"layer bins sum to {at} but the slab is {slab_depth} deep; "
            f"the band edges cannot be trusted."
        )
    return edges


def draw_flat(ax, img, cmap, lim, title, border_sets, label_xy):
    """Draw one flatmap panel with the area borders and labels; returns the image.

    `border_sets` holds the borders of each hemisphere separately: the two carry
    the same acronyms as keys, so merging them into one dict would keep one
    hemisphere only.
    """
    ax.imshow(np.zeros(img.shape + (4,)), origin="upper")
    h = ax.imshow(
        np.ma.masked_invalid(img),
        cmap=cmap,
        vmin=lim[0],
        vmax=lim[1],
        origin="upper",
        interpolation="nearest",
    )

    # area borders, brighter and wider for the areas the argument is about
    for borders in border_sets:
        for acro, coords in borders.items():
            hi = acro in HIGHLIGHT
            ax.plot(
                coords[:, 0],
                coords[:, 1],
                c="#e6e6e6" if hi else "#9a9a9a",
                lw=0.7 if hi else 0.3,
                alpha=0.9 if hi else 0.65,
            )

    # area names
    for acro, (x, y) in label_xy.items():
        dx, dy = LABEL_NUDGE.get(acro, (0, 0))
        t = ax.text(x + dx, y + dy, acro, color="w", fontsize=7, ha="center", va="center")

        # the map runs from near-white to dark red, so plain text loses either end
        # of it; a soft dark halo keeps it readable without shouting
        t.set_path_effects(
            [path_effects.withStroke(linewidth=1.3, foreground=(0, 0, 0, 0.55))]
        )
    ax.set_title(title, color="w", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("k")
    for s in ax.spines.values():
        s.set_visible(False)
    return h


def flatmaps(reading, vals, signed, lim_mean, cmaps, sigma_txt, n, out_dir):
    """Draw the flatmaps of `reading`: through the full depth, and by depth band.

    Each cohort's folded volumes are mirrored back to both hemispheres at 10 um and
    projected along the streamlines, value and mask apart, so that their ratio
    averages over tissue only: summed through the full depth for the first figure,
    and over the depth bins of each of BANDS for the second.
    """
    # ccf_streamlines only for the flatmaps, so a run with --no-flatmap does without it
    from ccf_streamlines.projection import (
        BoundaryFinder,
        Isocortex2dProjector,
        Isocortex3dProjector,
    )

    cmap_mean, rdbu = cmaps
    proj_file = os.path.join(ASSETS, "flatmap_butterfly.h5")
    path_file = os.path.join(ASSETS, "surface_paths_10_v3.h5")

    # area borders of both hemispheres, and where each name goes
    bf = BoundaryFinder(
        projected_atlas_file=os.path.join(ASSETS, "flatmap_butterfly.nrrd"),
        labels_file=os.path.join(ASSETS, "labelDescription_ITKSNAPColor.txt"),
    )
    left = {k: v for k, v in bf.region_boundaries().items() if len(v)}
    right = {
        k: v
        for k, v in bf.region_boundaries(
            hemisphere="right_for_both",
            view_space_for_other_hemisphere="flatmap_butterfly",
        ).items()
        if len(v)
    }
    label_xy = {a: left[a].mean(axis=0) for a in LABEL_AREAS if a in left}

    # projectors: onto the flat surface, and into a slab of depth bins by layer
    p2 = Isocortex2dProjector(
        proj_file,
        path_file,
        hemisphere="both",
        view_space_for_other_hemisphere="flatmap_butterfly",
    )
    p3 = Isocortex3dProjector(
        proj_file,
        path_file,
        thickness_type="normalized_layers",
        layer_thicknesses=layer_thicknesses(),
        streamline_layer_thickness_file=os.path.join(ASSETS, "cortical_layers_10_v2.h5"),
        hemisphere="both",
        view_space_for_other_hemisphere="flatmap_butterfly",
    )

    def to_10um(half):
        """Mirror a folded half back onto both hemispheres and repeat it to 10 um."""
        full = np.concatenate([half, half[:, :, ::-1]], axis=2)

        # each voxel repeated 2 x 2 x 2: the exact inverse of the block mean
        return np.repeat(np.repeat(np.repeat(full, 2, 0), 2, 1), 2, 2)

    # each cohort through the full depth, and its slab
    flat, slab = {}, {}
    for cohort in (YOUNG, ADULT):
        v10, m10 = to_10um(vals[cohort][0]), to_10um(vals[cohort][1])

        # project_volume indexes [x, y] while region_boundaries returns (x, y) to
        # plot, and imshow wants [row, col] = [y, x]: transpose once, here
        s_v = p2.project_volume(v10, kind="sum").T
        s_m = p2.project_volume(m10, kind="sum").T
        flat[cohort] = np.where(s_m > 0, s_v / np.maximum(s_m, 1e-6), np.nan)
        slab[cohort] = (
            p3.project_volume(v10).swapaxes(0, 1),
            p3.project_volume(m10).swapaxes(0, 1),
        )
        del v10, m10
    edges = band_edges(p3, slab[YOUNG][0].shape[2])

    def diff_of(y, a):
        """Young against adult: a difference if signed, else a floored log2 ratio."""
        return (
            (y - a)
            if signed
            else np.log2(np.maximum(y, LOG2_FLOOR) / np.maximum(a, LOG2_FLOOR))
        )

    # through the full depth: young, adult, difference
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.2), facecolor="k")
    panels = (
        (flat[YOUNG], cmap_mean, lim_mean[0], f"young (n = {n[YOUNG]})   {reading}"),
        (flat[ADULT], cmap_mean, lim_mean[0], f"adult (n = {n[ADULT]})   {reading}"),
        (
            diff_of(flat[YOUNG], flat[ADULT]),
            rdbu,
            lim_mean[1],
            "young - adult" if signed else "log2( young / adult )",
        ),
    )
    for ax, (im, cm, lim, ttl) in zip(axes, panels):
        h = draw_flat(ax, im, cm, lim, ttl, (left, right), label_xy)
        cb = fig.colorbar(h, ax=ax, fraction=0.03, pad=0.01)
        cb.ax.yaxis.set_tick_params(color="w", labelcolor="w")
    fig.suptitle(
        f"Cortical flatmap, averaged through the full depth  -  {reading}, {sigma_txt}",
        color="w",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(fig, os.path.join(out_dir, f"detail_flatmap_{reading}.png"), 110)
    plt.close(fig)

    # by depth band, one row per band
    fig, axes = plt.subplots(3, 3, figsize=(19, 13), facecolor="k")
    for row, (bname, keys) in enumerate(BANDS):
        lo, hi = edges[keys[0]][0], edges[keys[-1]][1]
        band = {}
        for cohort in (YOUNG, ADULT):
            sv, sm = slab[cohort]
            num = sv[:, :, lo:hi].sum(axis=2)
            den = sm[:, :, lo:hi].sum(axis=2)
            band[cohort] = np.where(den > 0, num / np.maximum(den, 1e-6), np.nan)
        panels = (
            (band[YOUNG], cmap_mean, lim_mean[0], f"young   {bname}"),
            (band[ADULT], cmap_mean, lim_mean[0], f"adult   {bname}"),
            (
                diff_of(band[YOUNG], band[ADULT]),
                rdbu,
                lim_mean[1],
                f"young - adult   {bname}",
            ),
        )
        for ax, (im, cm, lim, ttl) in zip(axes[row], panels):
            h = draw_flat(ax, im, cm, lim, ttl, (left, right), label_xy)
            cb = fig.colorbar(h, ax=ax, fraction=0.03, pad=0.01)
            cb.ax.yaxis.set_tick_params(color="w", labelcolor="w")
    fig.suptitle(
        "Cortical flatmap by layer, depth normalised per streamline  -  "
        f"{reading}, {sigma_txt}",
        color="w",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, os.path.join(out_dir, f"detail_flatmap_layers_{reading}.png"), 110)
    plt.close(fig)
    print(
        f"  wrote detail_flatmap_{reading}.png and detail_flatmap_layers_{reading}.png",
        flush=True,
    )


def main(readings, plane, vmax, dlim, sigma, want_video, want_flatmap, cmap_name=None):
    """Draw the close-up views of each of `readings`.

    `vmax` and `dlim` replace VMAX and DLIM when given; `sigma` is one sigma or
    three, in 20 um voxels, 0 for none. With `cmap_name` the mean panels use that
    colormap and the set goes into a subfolder of that name.
    """
    # a named colormap sends the whole set to its own subfolder, so the default
    # figures are never overwritten by an experiment with the colours
    out_dir = OUT if cmap_name is None else os.path.join(OUT, cmap_name)
    os.makedirs(out_dir, exist_ok=True)

    # hot up to 0.82 of its range; every colormap transparent where there is no value
    hot = plt.get_cmap("hot")
    hot_cut = LinearSegmentedColormap.from_list("hot_cut", hot(np.linspace(0, 0.82, 256)))
    hot_cut.set_bad((0, 0, 0, 0))
    puor = plt.get_cmap("PuOr_r").copy()
    puor.set_bad((0, 0, 0, 0))
    rdbu = plt.get_cmap("RdBu_r").copy()
    rdbu.set_bad((0, 0, 0, 0))

    # what every title says about the smoothing and the colormap
    sig = np.atleast_1d(sigma).astype(float)
    cmap_txt = "" if cmap_name is None else f", {cmap_name} colormap"
    sigma_txt = (
        "no smoothing"
        if not np.any(sig)
        else "smoothed sigma "
        + " x ".join(
            f"{v * 20:.0f}" for v in (sig if sig.size > 1 else np.repeat(sig, 3))
        )
        + " um (AP x DV x ML)"
    ) + cmap_txt
    n = {c: cohort_size(c) for c in (YOUNG, ADULT)}
    for reading in readings:
        t0 = time.time()

        # colour limits and colormaps of this reading
        signed = reading in SIGNED_READINGS
        vm = VMAX[reading] if vmax is None else vmax
        dl = DLIM[reading] if dlim is None else dlim
        lim_mean = ((-vm, vm) if signed else (0, vm), (-dl, dl))
        if cmap_name is None:
            cmap_mean = puor if signed else hot_cut
        else:
            cmap_mean = plt.get_cmap(cmap_name).copy()
            cmap_mean.set_bad((0, 0, 0, 0))
        cmaps = (cmap_mean, rdbu)

        # volumes, then the figures
        print(f"{reading}: preparing volumes ({sigma_txt})", flush=True)
        vals = prepare(reading, sigma)
        diff, both = difference(vals, signed)
        coronal(
            reading,
            vals,
            diff,
            both,
            plane,
            lim_mean,
            cmaps,
            sigma_txt,
            want_video,
            n,
            signed,
            out_dir,
        )
        if want_flatmap:
            flatmaps(reading, vals, signed, lim_mean, cmaps, sigma_txt, n, out_dir)
        print(f"{reading}: done in {time.time() - t0:.0f} s", flush=True)
