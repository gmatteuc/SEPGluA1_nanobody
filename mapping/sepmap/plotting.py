"""The palette, the colormaps, the save function and the drawing the figures share.

The colours are those of docs/STYLE.md (Figures), so the young against adult
figures, the adult figures and the ISH figures read as one set: red for the
young group and for what a figure is about, dark and mid grey for the two adult
groups and for comparison and context, dark blue for the receptor subunits, and a
flat grey for no data, which no data colormap produces.

Intensity maps use hot cut at 0.82 of its range, so the brightest values read as
yellow and never as the white of an empty page; masked values are transparent, so
the ground under an image (the atlas in grey, the black of the video frames)
shows where there is no value.

The coronal frames of the videos and of the close-up still share one drawing
function: the atlas in dark grey under the data, the area borders on top, and
the acronym of every structure large enough to name. draw_plane is the same
drawing for one panel of a figure, and paint fills a plane with one value per
structure.

The ISH figures add the light greys of context and of the reading notes, the two
channels and their per-mouse dots (orange for nano, yellow for autofluorescence, as
sep_palette in MATLAB), four colours for the groups of divisions in a scatter of
structures (cortex, hippocampal formation, thalamus, the rest), the colours of the
gene sets, a pale blue for the band of a null distribution, the blues of the green
(SEP) channel, and bars_grey, the grey of a bar whose darkness says how reliable
its value is.

Only numpy, scipy, matplotlib and config are imported, so the flatmap
environment (young_vs_adult.closeup) can import this module too. Imported by
the modules that draw.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Colormap, LinearSegmentedColormap
from matplotlib.image import AxesImage
from scipy.ndimage import center_of_mass

from sepmap.config import SETTINGS

# the 20 um voxels a structure needs in a coronal plane to get its acronym drawn
VIDEOS = SETTINGS["videos"]

# the palette: red for the young group and for what a figure shows, dark and mid
# grey for the naive and rws adults and for comparison and context, dark blue for
# the receptor subunits
RED = "#c0392b"
DARK_GREY = "#555555"
MID_GREY = "#9a9a9a"
DARK_BLUE = "#1f3b73"

# no data, drawn flat under the data
NO_DATA_GREY = "#bfbfbf"

# the light greys of context: what a figure is not about, and a fill behind the data;
# and the grey of the reading notes at the foot of a figure
LIGHT_GREY = "#c8c8c8"
PALE_GREY = "#e6e6e6"
NOTE_GREY = "0.4"

# two blues beside the dark blue of the subunits: the thalamus and the green (SEP)
# channel, and the paler density and GABAergic markers
MID_BLUE = "#3a6db5"
PALE_BLUE = "#7f9cc9"

# the groups of the young against adult figures
GROUP_COLOURS = {"young": RED, "naive": DARK_GREY, "rws": MID_GREY}

# the two channels, their per-mouse dots, and the line joining the values of one mouse
NANO = (0.95, 0.55, 0.10)
AUTO = (0.95, 0.85, 0.20)
NANO_DOT = (0.65, 0.30, 0.00)
AUTO_DOT = (0.70, 0.60, 0.00)
PAIR_LINE = (0.6, 0.6, 0.6)

# the groups of divisions in a scatter of structures: few enough to tell apart, and
# the three that carry most of a whole-brain correlation (cortex and hippocampus high,
# thalamus low) each with a colour of its own
DIVISION_GROUP = {
    "Isocortex": "cortex",
    "OLF": "cortex",
    "CTXsp": "cortex",
    "HPF": "hippocampal formation",
    "TH": "thalamus",
    "STR": "other grey matter",
    "PAL": "other grey matter",
    "HY": "other grey matter",
    "MB": "other grey matter",
    "P": "other grey matter",
    "MY": "other grey matter",
    "CB": "other grey matter",
}
DIVISION_GROUP_COLOURS = {
    "cortex": "#e07b00",
    "hippocampal formation": RED,
    "thalamus": MID_BLUE,
    "other grey matter": DARK_GREY,
}

# the gene sets of the ISH analysis; the subunits keep the dark blue they have everywhere
SET_COLOURS = {
    "subunits": DARK_BLUE,
    "localisation": RED,
    "other postsynaptic": DARK_GREY,
    "presynaptic": MID_GREY,
    "GABAergic markers": PALE_BLUE,
    "glia": LIGHT_GREY,
}

# a gene in no set, paler than any set
NO_SET_GREY = "#ececec"

# the 95% band of a null distribution, drawn behind the data
NULL_BAND = "#c9d6ea"

# the density step of the variance budget, beside the orange of abundance
DENSITY_BLUE = PALE_BLUE

# the third channel, the tag's own green (SEP) fluorescence, in blue so the three
# channels stay apart: its bars and boxes, its per-mouse dots, and what is left of it
# once autofluorescence is regressed out
SEP = MID_BLUE
SEP_DOT = "#24427f"
SEP_REMAINDER = "#8fb3e0"


def set_style() -> None:
    """The figures' style: white ground, fonts of 7 to 10 points, no top or right spine.

    Called once by every run script that draws, and by a notebook before it draws.
    Type 42 fonts keep the text of an EPS editable as text.
    """
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "figure.titlesize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "ps.fonttype": 42,
            "pdf.fonttype": 42,
        }
    )


def bars_grey(t: float | np.ndarray, t_max: float) -> np.ndarray:
    """The grey of a bar whose value has reliability `t`: 0.78 at 0, black at `t_max`.

    Darker is more reliable, and never white, so a bar of low reliability still
    shows on a white page (sep_palette('bars') in MATLAB). Returns one RGB row per
    value, shape (n, 3); t is clamped to 0..t_max, and NaN reads as 0.
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    share = np.clip(np.nan_to_num(t, nan=0.0) / t_max, 0.0, 1.0)
    level = 0.78 * (1.0 - share)
    return np.repeat(level[:, None], 3, axis=1)


def hot_cut() -> LinearSegmentedColormap:
    """hot up to 0.82 of its range, transparent where there is no value."""
    hot = plt.get_cmap("hot")
    cmap = LinearSegmentedColormap.from_list("hot_cut", hot(np.linspace(0, 0.82, 256)))
    cmap.set_bad((0, 0, 0, 0))
    return cmap


def transparent_bad(name: str) -> Colormap:
    """A copy of the colormap `name`, transparent where there is no value."""
    cmap = plt.get_cmap(name).copy()
    cmap.set_bad((0, 0, 0, 0))
    return cmap


def save_figure(
    fig: plt.Figure,
    path: Path,
    dpi: int,
    facecolor: str | None = None,
    eps: bool = True,
) -> None:
    """Save `fig` as a PNG at `path`, and with `eps` as an EPS beside it.

    `facecolor` is the PNG's background, the figure's default when None. The
    figure is not closed: a caller may draw into it again (a video after its
    still).

    Windows refuses to overwrite a PNG that an image viewer holds open. The
    figure then goes to <name>_new.png, with a note, so a run that writes
    several figures does not lose the rest because one of them was being
    looked at.

    The EPS is what goes into a figure for a paper. PostScript has no
    transparency, so the image layers are rasterised and composited by Agg
    first; otherwise a no-data region, transparent here, would come out opaque
    black instead of showing the ground beneath it. Text, lines and axes stay
    vector, the part that has to be editable.
    """
    try:
        fig.savefig(path, dpi=dpi, facecolor=facecolor)
    except OSError:
        alt = path.with_name(path.name.replace(".png", "_new.png"))
        fig.savefig(alt, dpi=dpi, facecolor=facecolor)
        print(f"  {path.name} is open elsewhere; wrote {alt.name} instead", flush=True)
    if not eps:
        return

    # the EPS, with the image layers rasterised
    eps_path = path.with_suffix(".eps")
    for ax in fig.axes:
        for im in ax.images:
            im.set_rasterized(True)
    try:
        fig.savefig(eps_path, dpi=dpi, facecolor=fig.get_facecolor(), format="eps")
    except OSError:
        print(
            f"  {eps_path.name} is open elsewhere; the PNG was still written", flush=True
        )


def tidy(ax: plt.Axes) -> None:
    """Small tick labels, no top or right spine."""
    ax.tick_params(labelsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


# ===== Coronal frames =====


def boundaries(lab: np.ndarray) -> np.ndarray:
    """Pixels of the label image `lab` that border another label, inside the atlas."""
    b = np.zeros(lab.shape, bool)
    b[1:, :] |= lab[1:, :] != lab[:-1, :]
    b[:, 1:] |= lab[:, 1:] != lab[:, :-1]
    return b & (lab > 0)


def paint(lab: np.ndarray, value: dict[str, float], names: dict[int, str]) -> np.ndarray:
    """A plane of structure values: each pixel its structure's value, NaN if none.

    `lab` holds parcellation indices, `names` the structure of each index and `value`
    a value per structure name.
    """
    out = np.full(lab.shape, np.nan, dtype=np.float32)
    for idx in np.unique(lab):
        if idx == 0:
            continue
        v = value.get(names.get(int(idx), ""))
        if v is not None:
            out[lab == idx] = v
    return out


def draw_borders(ax: plt.Axes, lab: np.ndarray) -> None:
    """The area borders of a plane on top of what is drawn, light and half transparent."""
    ov = np.zeros(lab.shape + (4,))
    ov[boundaries(lab)] = (0.8, 0.8, 0.8, 0.55)
    ax.imshow(ov, origin="upper", interpolation="nearest", aspect="equal")


def draw_plane(
    ax: plt.Axes,
    img: np.ndarray,
    lab: np.ndarray,
    cmap: Colormap,
    vmin: float,
    vmax: float,
) -> AxesImage:
    """One coronal plane in a panel: the atlas dark grey under `img`, borders on top.

    `img` and `lab` are (DV, ML), dorsal up; `img` is NaN where there is no value,
    so the dark grey atlas shows there. The panel is black, without ticks or frame.
    Returns the image of the data, for a colour bar.
    """
    inside = lab > 0

    # the atlas in dark grey under the data, black outside it
    bg = np.zeros(lab.shape + (4,))
    bg[inside] = (0.23, 0.23, 0.23, 1.0)
    ax.imshow(bg, origin="upper", interpolation="nearest", aspect="equal")
    h = ax.imshow(
        np.ma.masked_invalid(np.where(inside, img, np.nan)),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        origin="upper",
        interpolation="nearest",
        aspect="equal",
    )

    draw_borders(ax, lab)
    ax.set_facecolor("k")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return h


def coronal_figure() -> tuple[plt.Figure, list[plt.Axes], list[plt.Axes]]:
    """One 1920 x 760 figure on black: three panels and their colour bars."""
    fig = plt.figure(figsize=(19.2, 7.6), dpi=100, facecolor="k")
    axes = [fig.add_axes([0.02 + i * 0.325, 0.05, 0.27, 0.82]) for i in range(3)]
    caxes = [fig.add_axes([0.295 + i * 0.325, 0.12, 0.009, 0.68]) for i in range(3)]
    return fig, axes, caxes


def coronal_frame(
    fig: plt.Figure,
    axes: list[plt.Axes],
    caxes: list[plt.Axes],
    k: int,
    panels: tuple,
    ann_h: np.ndarray,
    acro: dict[int, str],
    header: str,
    vector_outline: bool = False,
) -> None:
    """Draw plane `k` into the panels of `fig`, and the header above them.

    `panels` holds (image, colormap, limits, title) per panel. The atlas is dark
    grey under the data and its borders lie on top: as lines with
    `vector_outline`, otherwise as a pixel overlay. Structures with at least
    videos.min_label_area voxels in the plane get their acronym.
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
            if m.sum() < VIDEOS["min_label_area"]:
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
