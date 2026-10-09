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
tools\\requirements_flat.txt), and imports only config, hemispheres and plotting
from the package. Its assets, about 0.6 GB, are fetched once into atlas_flatmap/
under the data root from the Allen Institute's ccf_streamlines_assets folder,
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
import time
from pathlib import Path

import imageio_ffmpeg
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.colors import Colormap
from matplotlib.image import AxesImage
from scipy.ndimage import gaussian_filter

from sepmap.config import DATA, SETTINGS
from sepmap.hemispheres import fold, fold_count
from sepmap.plotting import (
    coronal_figure,
    coronal_frame,
    hot_cut,
    save_figure,
    transparent_bad,
)

CLOSEUP = SETTINGS["closeup"]
READINGS = SETTINGS["readings"]
VIDEOS = SETTINGS["videos"]
YOUNG_VS_ADULT = SETTINGS["young_vs_adult"]

ASSETS = DATA / "atlas_flatmap"
CCF_ROOT = DATA / "comparisons_v2" / "ccf"
OUT = DATA / "comparisons_v2" / "young_vs_adult"
CSV_MAP = DATA / "atlas" / "parcellation_to_parcellation_term_membership.csv"

YOUNG, ADULT = "young", "adult"

# this module runs in the flatmap environment and cannot import the analysis chain,
# so what it shares with the chain comes from settings.toml, which both read: the
# brains a voxel needs (young_vs_adult, as in compare, so the maps agree) and the
# log2 floor (readings, the value volumes.cohort uses)

# readings that are a position within a brain's own range, as volumes.cohort's
# SIGNED_READINGS: compared by difference, drawn diverging from -v to +v
SIGNED_READINGS = tuple(READINGS["signed"])

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


def cohort_size(cohort: str) -> int:
    """Number of brains in `cohort`, counted in the mice list volumes.cohort wrote.

    Counted rather than typed, so a caption cannot go stale when a brain is added.
    """
    with open(CCF_ROOT / cohort / "mice.txt", encoding="utf-8") as fh:
        return sum(1 for line in fh if line.strip())


def smooth_within(v: np.ndarray, m: np.ndarray, sigma: float | list[float]) -> np.ndarray:
    """`v` smoothed inside the mask `m` with Gaussian `sigma`, 0 outside it.

    Mask normalised: tissue is never averaged with the zeros around it. Also used by
    adult.layers, which prepares each adult the way a cohort is prepared here.
    """
    num = gaussian_filter(v, sigma)
    den = gaussian_filter(m, sigma)
    return np.where(m > 0, num / np.maximum(den, 1e-6), 0).astype(np.float32)


def prepare(
    reading: str, sigma: float | list[float]
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Value and mask of each cohort for `reading`, hemispheres folded, at 20 um.

    Returns {cohort: (value, mask)}, the one source of all four figures. The mask
    holds the voxels with enough brains and a value; the value is 0 outside it and,
    with `sigma`, smoothed inside it. The mask itself is not smoothed: blurring a
    value inside the tissue is cosmetic, growing the tissue would not be.
    """
    out = {}
    for cohort, min_n in (
        (YOUNG, YOUNG_VS_ADULT["min_n_young"]),
        (ADULT, YOUNG_VS_ADULT["min_n_adult"]),
    ):
        mean = fold(np.load(CCF_ROOT / cohort / f"{reading}_mean.npy"))
        n = fold_count(np.load(CCF_ROOT / cohort / f"{reading}_n.npy"))
        m = ((n >= min_n) & np.isfinite(mean)).astype(np.float32)
        v = np.where(m > 0, mean, 0).astype(np.float32)

        # mask-normalised smoothing: tissue is never averaged with the zeros around it
        if np.any(sigma):
            v = smooth_within(v, m, sigma)
        out[cohort] = (v, m)
    return out


def difference(
    vals: dict[str, tuple[np.ndarray, np.ndarray]], signed: bool
) -> tuple[np.ndarray, np.ndarray]:
    """Young against adult where both cohorts have a value (NaN elsewhere), and that mask.

    A signed reading is compared by difference, the others by log2 ratio, each value
    floored at readings.log2_floor.
    """
    y, a = vals[YOUNG], vals[ADULT]
    both = (y[1] > 0) & (a[1] > 0)
    if signed:
        d = y[0] - a[0]
    else:
        floor = READINGS["log2_floor"]
        d = np.log2(np.maximum(y[0], floor) / np.maximum(a[0], floor))
    return np.where(both, d, np.nan), both


def shown(value: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """`value` where `mask` is set, NaN elsewhere."""
    return np.where(mask > 0, value, np.nan)


# ===== Coronal plane and video =====


def annotation_half() -> np.ndarray:
    """The CCF annotation at 20 um, the half that the folded volumes cover."""
    ann = np.asarray(nib.load(DATA / "atlas" / "annotation_10.nii.gz").dataobj)[
        ::2, ::2, ::2
    ]
    return ann[:, :, : ann.shape[2] // 2]


def structure_acronyms() -> dict[int, str]:
    """Structure acronyms by parcellation index, read from CSV_MAP.

    The package reads them through volumes.per_mouse, which this module cannot
    import in the flatmap environment.
    """
    acro = {}
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["parcellation_term_set_name"] == "structure":
                acro[int(row["parcellation_index"])] = row["parcellation_term_acronym"]
    return acro


def write_detail_video(
    reading: str,
    fig: plt.Figure,
    axes: list[plt.Axes],
    caxes: list[plt.Axes],
    both: np.ndarray,
    ann_h: np.ndarray,
    acro: dict[int, str],
    panels_at,
    header_at,
    out_dir: Path,
) -> None:
    """Redraw the still's figure plane by plane into detail_video_<reading>.mp4.

    `panels_at` and `header_at` give the panels and the header of a plane; the
    video runs over the planes with more than 200 compared voxels.
    """
    t0 = time.time()
    frames = [i for i in range(ann_h.shape[0]) if both[i].sum() > 200]
    out = out_dir / f"detail_video_{reading}.mp4"
    writer = imageio_ffmpeg.write_frames(
        out, (1920, 760), fps=VIDEOS["fps"], quality=7, macro_block_size=8
    )
    writer.send(None)
    for i in frames:
        coronal_frame(fig, axes, caxes, i, panels_at(i), ann_h, acro, header_at(i))
        fig.canvas.draw()
        writer.send(np.ascontiguousarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3]))
    writer.close()
    print(
        f"  wrote {out.name}, {len(frames)} frames, {time.time() - t0:.0f} s",
        flush=True,
    )


def coronal(
    reading: str,
    vals: dict[str, tuple[np.ndarray, np.ndarray]],
    diff: np.ndarray,
    both: np.ndarray,
    plane: int,
    lim_mean: tuple[tuple[float, float], tuple[float, float]],
    cmaps: tuple[Colormap, Colormap],
    sigma_txt: str,
    want_video: bool,
    n: dict[str, int],
    signed: bool,
    out_dir: Path,
) -> None:
    """Draw the still at CCF `plane` and, with `want_video`, the video of every plane.

    Young, adult and their difference side by side. `lim_mean` holds the limits of
    the means and of the difference, `cmaps` their colormaps, `n` the brains of
    each cohort. The video runs over the planes with more than 200 compared voxels.
    """
    ann_h = annotation_half()
    acro = structure_acronyms()
    cmap_mean, rdbu = cmaps

    # one 1920 x 760 figure, three panels with their colour bars
    fig, axes, caxes = coronal_figure()

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
    out = out_dir / f"detail_plane{plane}_{reading}.png"
    save_figure(fig, out, dpi=100, facecolor="k")
    print(f"  wrote {out.name}", flush=True)

    # the video, redrawing the same figure plane by plane
    if want_video:
        write_detail_video(
            reading, fig, axes, caxes, both, ann_h, acro, panels_at, header_at, out_dir
        )
    plt.close(fig)


# ===== Flatmaps =====


def layer_thicknesses() -> dict[str, float]:
    """Thickness of each cortical layer in um, from the Allen's average depths.

    The file gives the depth of each layer's lower border below the pia, so the
    thicknesses are the differences between them.
    """
    with open(ASSETS / "avg_layer_depths.json") as fh:
        d = json.load(fh)
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


def band_edges(p3, slab_depth: int) -> dict[str, tuple[int, int]]:
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
        raise RuntimeError(
            f"layer bins sum to {at} but the slab is {slab_depth} deep; "
            "the band edges cannot be trusted."
        )
    return edges


def draw_flat(
    ax: plt.Axes,
    img: np.ndarray,
    cmap: Colormap,
    lim: tuple[float, float],
    title: str,
    border_sets: tuple[dict[str, np.ndarray], ...],
    label_xy: dict[str, np.ndarray],
) -> AxesImage:
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


def flatmap_borders() -> tuple[
    dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray]
]:
    """Area borders of both hemispheres on the flatmap, and where each name goes."""
    # ccf_streamlines only for the flatmaps, so a run with --no-flatmap does without it
    from ccf_streamlines.projection import BoundaryFinder

    bf = BoundaryFinder(
        projected_atlas_file=str(ASSETS / "flatmap_butterfly.nrrd"),
        labels_file=str(ASSETS / "labelDescription_ITKSNAPColor.txt"),
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
    return left, right, label_xy


def flatmap_projectors():
    """The projectors onto the flat surface, and into a slab of depth bins by layer."""
    from ccf_streamlines.projection import Isocortex2dProjector, Isocortex3dProjector

    proj_file = str(ASSETS / "flatmap_butterfly.h5")
    path_file = str(ASSETS / "surface_paths_10_v3.h5")
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
        streamline_layer_thickness_file=str(ASSETS / "cortical_layers_10_v2.h5"),
        hemisphere="both",
        view_space_for_other_hemisphere="flatmap_butterfly",
    )
    return p2, p3


def to_10um(half: np.ndarray) -> np.ndarray:
    """Mirror a folded half back onto both hemispheres and repeat it to 10 um."""
    full = np.concatenate([half, half[:, :, ::-1]], axis=2)

    # each voxel repeated 2 x 2 x 2: the exact inverse of the block mean
    return np.repeat(np.repeat(np.repeat(full, 2, 0), 2, 1), 2, 2)


def project_slab(v10: np.ndarray, m10: np.ndarray, p3) -> tuple[np.ndarray, np.ndarray]:
    """The value and mask slabs of a volume at 10 um, as [row, col, depth bin].

    Projected apart, so that their ratio averages over tissue only. project_volume
    indexes [x, y], and imshow wants [row, col] = [y, x]: swapped once, here.
    """
    return p3.project_volume(v10).swapaxes(0, 1), p3.project_volume(m10).swapaxes(0, 1)


def band_average(sv: np.ndarray, sm: np.ndarray, lo: int, hi: int) -> np.ndarray:
    """The mean of a value slab over depth bins `lo` to `hi`, over tissue only.

    `sm` is the mask slab; NaN where the band holds no tissue.
    """
    num = sv[:, :, lo:hi].sum(axis=2)
    den = sm[:, :, lo:hi].sum(axis=2)
    return np.where(den > 0, num / np.maximum(den, 1e-6), np.nan)


def project_cohorts(
    vals: dict[str, tuple[np.ndarray, np.ndarray]], p2, p3
) -> tuple[dict[str, np.ndarray], dict[str, tuple[np.ndarray, np.ndarray]]]:
    """Each cohort through the full depth, and its slab of depth bins.

    Value and mask are projected apart, so that their ratio averages over tissue
    only. Returns ({cohort: flatmap}, {cohort: (value slab, mask slab)}).
    """
    flat, slab = {}, {}
    for cohort in (YOUNG, ADULT):
        v10, m10 = to_10um(vals[cohort][0]), to_10um(vals[cohort][1])

        # project_volume indexes [x, y] while region_boundaries returns (x, y) to
        # plot, and imshow wants [row, col] = [y, x]: transpose once, here
        s_v = p2.project_volume(v10, kind="sum").T
        s_m = p2.project_volume(m10, kind="sum").T
        flat[cohort] = np.where(s_m > 0, s_v / np.maximum(s_m, 1e-6), np.nan)
        slab[cohort] = project_slab(v10, m10, p3)
        del v10, m10
    return flat, slab


def diff_of(y: np.ndarray, a: np.ndarray, signed: bool) -> np.ndarray:
    """Young against adult: a difference if signed, else a floored log2 ratio."""
    if signed:
        return y - a
    else:
        floor = READINGS["log2_floor"]
        return np.log2(np.maximum(y, floor) / np.maximum(a, floor))


def draw_full_depth(
    reading: str,
    flat: dict[str, np.ndarray],
    signed: bool,
    lim_mean: tuple[tuple[float, float], tuple[float, float]],
    cmaps: tuple[Colormap, Colormap],
    sigma_txt: str,
    n: dict[str, int],
    borders: tuple[dict[str, np.ndarray], dict[str, np.ndarray]],
    label_xy: dict[str, np.ndarray],
    out_dir: Path,
) -> None:
    """Draw detail_flatmap_<reading>.png: young, adult, difference, full depth."""
    cmap_mean, rdbu = cmaps
    left, right = borders
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.2), facecolor="k")
    panels = (
        (flat[YOUNG], cmap_mean, lim_mean[0], f"young (n = {n[YOUNG]})   {reading}"),
        (flat[ADULT], cmap_mean, lim_mean[0], f"adult (n = {n[ADULT]})   {reading}"),
        (
            diff_of(flat[YOUNG], flat[ADULT], signed),
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
    save_figure(
        fig,
        out_dir / f"detail_flatmap_{reading}.png",
        dpi=110,
        facecolor="k",
    )
    plt.close(fig)


def draw_bands(
    reading: str,
    slab: dict[str, tuple[np.ndarray, np.ndarray]],
    edges: dict[str, tuple[int, int]],
    signed: bool,
    lim_mean: tuple[tuple[float, float], tuple[float, float]],
    cmaps: tuple[Colormap, Colormap],
    sigma_txt: str,
    borders: tuple[dict[str, np.ndarray], dict[str, np.ndarray]],
    label_xy: dict[str, np.ndarray],
    out_dir: Path,
) -> None:
    """Draw detail_flatmap_layers_<reading>.png: one row per depth band of BANDS."""
    cmap_mean, rdbu = cmaps
    left, right = borders
    fig, axes = plt.subplots(3, 3, figsize=(19, 13), facecolor="k")
    for row, (bname, keys) in enumerate(BANDS):
        lo, hi = edges[keys[0]][0], edges[keys[-1]][1]
        band = {}
        for cohort in (YOUNG, ADULT):
            sv, sm = slab[cohort]
            band[cohort] = band_average(sv, sm, lo, hi)
        panels = (
            (band[YOUNG], cmap_mean, lim_mean[0], f"young   {bname}"),
            (band[ADULT], cmap_mean, lim_mean[0], f"adult   {bname}"),
            (
                diff_of(band[YOUNG], band[ADULT], signed),
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
    save_figure(
        fig,
        out_dir / f"detail_flatmap_layers_{reading}.png",
        dpi=110,
        facecolor="k",
    )
    plt.close(fig)


def flatmaps(
    reading: str,
    vals: dict[str, tuple[np.ndarray, np.ndarray]],
    signed: bool,
    lim_mean: tuple[tuple[float, float], tuple[float, float]],
    cmaps: tuple[Colormap, Colormap],
    sigma_txt: str,
    n: dict[str, int],
    out_dir: Path,
) -> None:
    """Draw the flatmaps of `reading`: through the full depth, and by depth band.

    Each cohort's folded volumes are mirrored back to both hemispheres at 10 um and
    projected along the streamlines, value and mask apart, so that their ratio
    averages over tissue only: summed through the full depth for the first figure,
    and over the depth bins of each of BANDS for the second.
    """
    # area borders of both hemispheres and where each name goes; the projectors
    left, right, label_xy = flatmap_borders()
    p2, p3 = flatmap_projectors()

    # each cohort through the full depth, and its slab
    flat, slab = project_cohorts(vals, p2, p3)
    edges = band_edges(p3, slab[YOUNG][0].shape[2])

    # through the full depth, then by depth band
    draw_full_depth(
        reading,
        flat,
        signed,
        lim_mean,
        cmaps,
        sigma_txt,
        n,
        (left, right),
        label_xy,
        out_dir,
    )
    draw_bands(
        reading,
        slab,
        edges,
        signed,
        lim_mean,
        cmaps,
        sigma_txt,
        (left, right),
        label_xy,
        out_dir,
    )
    print(
        f"  wrote detail_flatmap_{reading}.png and detail_flatmap_layers_{reading}.png",
        flush=True,
    )


def main(
    readings: list[str],
    plane: int,
    vmax: float | None,
    dlim: float | None,
    sigma: float | list[float],
    want_video: bool,
    want_flatmap: bool,
    cmap_name: str | None = None,
) -> None:
    """Draw the close-up views of each of `readings`.

    `vmax` and `dlim` replace closeup.vmax and closeup.dlim when given; `sigma` is
    one sigma or three, in 20 um voxels, 0 for none. With `cmap_name` the mean
    panels use that colormap and the set goes into a subfolder of that name.
    """
    # a named colormap sends the whole set to its own subfolder, so the default
    # figures are never overwritten by an experiment with the colours
    out_dir = OUT if cmap_name is None else OUT / cmap_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # hot up to 0.82 of its range; every colormap transparent where there is no value
    hot = hot_cut()
    puor = transparent_bad("PuOr_r")
    rdbu = transparent_bad("RdBu_r")

    # what every title says about the smoothing and the colormap
    sig = np.atleast_1d(sigma).astype(float)
    cmap_txt = "" if cmap_name is None else f", {cmap_name} colormap"
    if not np.any(sig):
        sigma_txt = "no smoothing"
    else:
        sigmas_3d = sig if sig.size > 1 else np.repeat(sig, 3)
        sigma_txt = (
            "smoothed sigma "
            + " x ".join(f"{v * 20:.0f}" for v in sigmas_3d)
            + " um (AP x DV x ML)"
        )
    sigma_txt += cmap_txt
    n = {c: cohort_size(c) for c in (YOUNG, ADULT)}
    for reading in readings:
        t0 = time.time()

        # colour limits and colormaps of this reading
        signed = reading in SIGNED_READINGS
        vm = CLOSEUP["vmax"][reading] if vmax is None else vmax
        dl = CLOSEUP["dlim"][reading] if dlim is None else dlim
        lim_mean = ((-vm, vm) if signed else (0, vm), (-dl, dl))
        if cmap_name is None:
            cmap_mean = puor if signed else hot
        else:
            cmap_mean = transparent_bad(cmap_name)
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
