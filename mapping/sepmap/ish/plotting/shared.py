"""What every figure of the ISH analysis shares: the page, the pieces, the save.

The guided figures are numbered in the order of the argument (03_beyond.png,
04_beyond_where.png, ...). A main figure answers one question in two to four
panels: its title is the number and the question, under it one line that says what
to take from it, computed from the tables, and at its foot a line or two in grey on
how to read it. Its detailed version, a supplementary figure, carries the same
number with an s (05s_one_comparison_detail.png; 03s1 and 03s2 where a main figure
has two), with every panel and number of the analysis. A takeaway that is a verdict
follows its numbers, so a rerun cannot leave a figure saying what its numbers
contradict.

The drawing follows docs/STYLE.md: maps of intensities and ranks in hot on black,
the atlas dark grey beneath and the area borders on top, structures left out flat
grey; scatters of structures with 35-point dots coloured by group of divisions;
bars whose grey says how reliable a value is; the subunits in dark blue. The
modules that compute draw nothing: every figure function takes tables and returns
the figure, saved as PNG and EPS at 150 dpi when `save` is given.

Here: the title, the line under it and the notes at the foot, panel titles, the
save; ranks and a plane painted with them; the scatter of structures by group of
divisions and its legend; names spread apart beside their dots; a spatial p as the
figures write it; one gene set as a column of dots against its null band.

Imported by the other modules of ish.plotting and by adult.plotting, whose figures of
part 1 are drawn with the same pieces.
"""

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.stats import rankdata

from sepmap.config import SETTINGS
from sepmap.ish.figure_index import FIGURE_NUMBERS, QUESTIONS
from sepmap.plotting import (
    DARK_BLUE,
    DARK_GREY,
    DIVISION_GROUP,
    DIVISION_GROUP_COLOURS,
    LIGHT_GREY,
    MID_GREY,
    NO_DATA_GREY,
    NOTE_GREY,
    NULL_BAND,
    draw_borders,
    draw_plane,
    hot_cut,
    paint,
    save_figure,
)

# the reliability below which a gene's Allen map does not reproduce, marked on the
# figures
ISH_PANEL_TEST = SETTINGS["ish_panel_test"]
LOW_RELIABILITY = ISH_PANEL_TEST["min_reliability"]

# dpi of every figure of the ISH analysis (docs/STYLE.md: new figures use 150)
DPI = 150

# characters of the line under a figure's title per inch of the figure's width,
# beyond which it wraps
HEADING_CHARS_PER_INCH = 13

# the lower end of hot for a map of ranks: below 0, so the lowest structure is dark
# red, not the black of the ground
RANK_FLOOR = -0.1

# the grey-matter divisions in the order the figures list them, front to back
DIVISION_ORDER = ["Isocortex", "OLF", "CTXsp", "HPF", "STR", "PAL", "TH", "HY", "MB"]
DIVISION_ORDER += ["P", "MY", "CB"]


# ===== Every figure =====


def heading(fig: plt.Figure, key: str, line: str) -> None:
    """The figure's number and question, and one line under it computed from the tables.

    What to take from a main figure, the numbers of a detailed one.

    Placed in inches from the top, so a short main figure and a tall detailed one
    leave the same room between the two lines.
    """
    height = fig.get_figheight()
    title = f"{FIGURE_NUMBERS[key]}.  {QUESTIONS[key]}"
    fig.text(0.5, 1 - 0.16 / height, title, ha="center", va="top", fontsize=12)

    # about 13 characters of 9.5-point text to an inch of the figure's width; a line
    # too long for it is cut into lines of about equal length, not one long and one
    # short
    width = int(HEADING_CHARS_PER_INCH * fig.get_figwidth())
    n_lines = max(1, -(-len(line) // width))
    fig.text(
        0.5,
        1 - 0.47 / height,
        textwrap.fill(line, -(-len(line) // n_lines) + 8),
        ha="center",
        va="top",
        fontsize=9.5,
        color=DARK_GREY,
    )


def footer(fig: plt.Figure, lines: list[str]) -> None:
    """The grey lines at the foot of a figure: how to read it, what would mean what."""
    fig.text(
        0.04,
        0.012,
        "\n".join(lines),
        ha="left",
        va="bottom",
        fontsize=8,
        color=NOTE_GREY,
        linespacing=1.5,
    )


def panel_title(ax: plt.Axes, letter: str, text: str, numbers: str = "") -> None:
    """A panel's title: its letter and what it shows, and a line of numbers under it."""
    title = f"{letter}.  {text}"
    if numbers:
        title += f"\n{numbers}"
    ax.set_title(title, loc="left", fontsize=9)


def ranks01(values: np.ndarray) -> np.ndarray:
    """Ranks scaled to 0 (lowest) .. 1 (highest), ties averaged."""
    r = rankdata(values)
    return (r - 1) / max(len(r) - 1, 1)


def overlay(ax: plt.Axes, mask: np.ndarray, colour: str, lab: np.ndarray) -> None:
    """Paint `mask` flat in `colour` over a plane, and the area borders on top again."""
    rgba = np.zeros(mask.shape + (4,))
    rgba[mask] = plt.matplotlib.colors.to_rgba(colour)
    ax.imshow(rgba, origin="upper", interpolation="nearest", aspect="equal")
    draw_borders(ax, lab)


def colour_bar(
    fig: plt.Figure, ax: plt.Axes, image, label: str, limits: tuple | None = None
) -> None:
    """A slim colour bar beside an image panel, cut to `limits` when given.

    Rank maps start hot a little below 0, so the lowest structure is dark red
    rather than the black of the ground; their bar is cut at 0 and 1.
    """
    cb = fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
    cb.ax.tick_params(labelsize=7)
    cb.set_label(label, fontsize=8)
    if limits is not None:
        cb.ax.set_ylim(*limits)


def saved(fig: plt.Figure, path: Path | None, eps: bool = True) -> plt.Figure:
    """Save the figure as PNG, and with `eps` as EPS, when a path is given; return it.

    The QC sheets are PNG only: they are for review, not for a paper, and 750 of
    them as EPS would take about 2 GB.
    """
    if path is not None:
        save_figure(fig, Path(path), DPI, eps=eps)
    return fig


# ===== Pieces of several figures =====


def p_text(p: float, n_surrogates: int) -> str:
    """A spatial p as the figures write it: the floor when no surrogate reached it."""
    floor = 1.0 / (n_surrogates + 1)
    if p <= floor * 1.5:
        return f"p ≤ {floor:.4f}"
    if p < 0.001:
        return f"p = {p:.4f}"
    return f"p = {p:.3f}"


def gene_colour(symbol: str, subunits: set[str]) -> str:
    """Dark blue for the AMPA receptor subunits, near black for every other gene."""
    return DARK_BLUE if symbol in subunits else "0.1"


def scatter_groups(ax: plt.Axes, x: np.ndarray, y: np.ndarray, groups: list[str]) -> None:
    """Structures as dots coloured by group of divisions: 35 points, no edge."""
    groups = np.asarray(groups)
    for group, colour in DIVISION_GROUP_COLOURS.items():
        mine = groups == group
        if mine.any():
            ax.scatter(
                x[mine],
                y[mine],
                s=35,
                color=colour,
                linewidths=0,
                alpha=0.85,
                zorder=2,
            )


def spread_positions(desired: np.ndarray, gap: float) -> np.ndarray:
    """Heights for a column of names: as near their dots as can be, `gap` apart.

    Neighbours closer than `gap` are pushed apart, half the shortfall each, until
    none is; the column stays centred on its dots rather than running off one end.
    """
    order = np.argsort(desired)
    y = np.asarray(desired, float)[order].copy()
    for _ in range(500):
        moved = False
        for k in range(len(y) - 1):
            short = gap - (y[k + 1] - y[k])
            if short > 1e-9:
                y[k] -= short / 2
                y[k + 1] += short / 2
                moved = True
        if not moved:
            break
    out = np.empty_like(y)
    out[order] = y
    return out


def spread_labels(
    ax: plt.Axes,
    points: list[tuple[float, float, str, str]],
    fontsize: float = 7,
    dx: float = 0.012,
) -> None:
    """Names beside their dots, moved apart vertically where they would overlap.

    `points` holds (x, y, text, colour) in data coordinates. Names go right of their
    dot; one that would sit on an earlier one (closer than a line of text, and within
    a name's width) is pushed away from the nearer edge of the axes, down in the upper
    half and up in the lower half, so no name leaves the axes, and is joined to its
    dot by a thin line when moved.
    """
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    gap = abs(y1 - y0) * 0.032 * fontsize / 7
    width = abs(x1 - x0) * 0.14
    middle = (y0 + y1) / 2
    upper = sorted((p for p in points if p[1] >= middle), key=lambda t: -t[1])
    lower = sorted((p for p in points if p[1] < middle), key=lambda t: t[1])
    placed: list[tuple[float, float]] = []
    for x, y, text, colour in upper + lower:
        ty = y
        step = -gap * 0.6 if y >= middle else gap * 0.6
        while any(abs(px - x) < width and abs(py - ty) < gap for px, py in placed):
            ty += step
        placed.append((x, ty))
        moved = abs(ty - y) > 1e-9
        ax.annotate(
            text,
            (x, y),
            xytext=(x + dx * abs(x1 - x0), ty),
            textcoords="data",
            fontsize=fontsize,
            color=colour,
            va="center",
            arrowprops=dict(arrowstyle="-", color=MID_GREY, lw=0.5) if moved else None,
            zorder=6,
        )


def group_handles() -> list:
    """Legend entries of the groups of divisions."""
    return [
        plt.Line2D([], [], ls="", marker="o", ms=6, mfc=c, mec="none", label=g)
        for g, c in DIVISION_GROUP_COLOURS.items()
    ]


def group_of(set_table: pd.DataFrame) -> dict[str, str]:
    """{structure: group of divisions} of the structure table."""
    return {
        s: DIVISION_GROUP.get(d, "other grey matter")
        for s, d in zip(set_table["structure"], set_table["division"])
    }


def acronym_of(set_table: pd.DataFrame) -> dict[str, str]:
    """{structure: acronym} of the structure table."""
    return dict(zip(set_table["structure"], set_table["acronym"]))


def rank_plane(
    fig: plt.Figure,
    ax: plt.Axes,
    ranks: dict[str, float],
    lab: np.ndarray,
    names: dict[int, str],
) -> None:
    """A plane painted with one rank per structure; structures without one flat grey."""
    image = draw_plane(ax, paint(lab, ranks, names), lab, hot_cut(), RANK_FLOOR, 1.0)
    shown = {names.get(int(i)) for i in np.unique(lab) if i}
    mask = np.isin(lab, [i for i, n in names.items() if n in shown - set(ranks)])
    overlay(ax, mask, NO_DATA_GREY, lab)
    colour_bar(fig, ax, image, "rank among structures (0 low, 1 high)", (0, 1))


# sets of at least this many genes get a violin outline; smaller ones only dots
VIOLIN_MIN = 15

# sets of at most this many genes have every gene named; larger ones their top three
# and their lowest
NAME_ALL_MAX = 8


def set_column(
    ax: plt.Axes,
    x: float,
    rows: pd.DataFrame,
    test: pd.Series | None,
    colour,
    rng: np.random.Generator,
    hollow: bool = False,
) -> None:
    """One set: its null band, violin, a dot per gene, the median, the named genes.

    `rows` holds the set's genes (gene_sets.csv), `test` its row of set_tests.csv
    for the nano map (None for the context group, which has no test); a set too
    small to test gets no band.
    """
    rho = rows["rho_nano"].to_numpy(float)
    if test is not None and bool(test["tested"]):
        ax.add_patch(
            Rectangle(
                (x - 0.4, test["null_lo"]),
                0.8,
                test["null_hi"] - test["null_lo"],
                color=NULL_BAND,
                lw=0,
                zorder=0,
            )
        )
    if len(rho) >= VIOLIN_MIN:
        parts = ax.violinplot(
            rho, positions=[x], widths=0.8, showextrema=False, showmedians=False
        )
        for body in parts["bodies"]:
            body.set_facecolor("none")
            body.set_edgecolor(MID_GREY)
            body.set_linewidth(0.8)
            body.set_alpha(1)
    width = 0.3 if len(rho) >= VIOLIN_MIN else 0.12
    jitter = rng.uniform(-width, width, len(rho))
    p9 = rows["p9_gene"].to_numpy(bool)
    face = "white" if hollow else colour
    # every dot has a thin dark edge so pale sets show on the pale band; P9's genes a
    # black ring
    ax.scatter(
        x + jitter,
        rho,
        s=18,
        facecolors=face,
        edgecolors=["0.05" if p else "0.35" for p in p9],
        linewidths=[1.0 if p else 0.4 for p in p9],
        zorder=3,
    )
    ax.plot([x - 0.32, x + 0.32], [np.median(rho)] * 2, color="0.1", lw=2, zorder=4)

    # the genes named: all of a small set, the top three and the lowest of a large one,
    # right of the column and moved apart where they would overlap
    order = np.argsort(rho)[::-1]
    if len(rho) <= NAME_ALL_MAX:
        named = order
    else:
        named = np.concatenate([order[:3], order[-1:]])
    heights = spread_positions(rho[named], 0.045)
    for i, y in zip(named, heights):
        ax.annotate(
            rows["symbol"].iloc[i],
            (x + jitter[i], rho[i]),
            xytext=(x + 0.36, y),
            textcoords="data",
            fontsize=6,
            va="center",
            color=DARK_BLUE if colour == DARK_BLUE else "0.15",
            arrowprops=dict(arrowstyle="-", color=LIGHT_GREY, lw=0.4),
            zorder=5,
        )
