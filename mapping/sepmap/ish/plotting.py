"""The figures of the ISH analysis: the numbered guide, the QC sheets, the gene sheets.

Each guided figure answers one question and is numbered in reading order
(01_structures.png, 02_genes.png, ...). Its title is the number and the question;
under it, a line of numbers computed from the tables; at its foot, in grey, how to
read it and what each outcome would mean. The verdict in words is in
docs/ISH_ANALYSIS.md, so a figure never states a conclusion that its numbers could
contradict after a rerun.

The drawing follows docs/STYLE.md: maps of intensities and ranks in hot on black,
the atlas dark grey beneath and the area borders on top, structures left out flat
grey; scatters of structures with 35-point dots coloured by group of divisions;
bars whose grey says how reliable a value is; the subunits in dark blue. The
modules that compute draw nothing: every figure function takes tables and returns
the figure, saved as PNG and EPS at 150 dpi when `save` is given.

Called by the run scripts of the ISH line.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import rankdata

from sepmap.plotting import (
    DARK_GREY,
    DIVISION_GROUP,
    DIVISION_GROUP_COLOURS,
    GROUP_COLOURS,
    MID_GREY,
    NO_DATA_GREY,
    PAIR_LINE,
    boundaries,
    draw_plane,
    hot_cut,
    save_figure,
    tidy,
)

# dpi of every figure of the ISH analysis (docs/STYLE.md: new figures use 150)
DPI = 150

# the grey of the reading notes at the foot of a figure, and of light context
NOTE_GREY = "0.4"
LIGHT_GREY = "#c8c8c8"

# the lower end of hot for a map of ranks: below 0, so the lowest structure is dark
# red, not the black of the ground
RANK_FLOOR = -0.1

# the grey-matter divisions in the order the figures list them, front to back
DIVISION_ORDER = ["Isocortex", "OLF", "CTXsp", "HPF", "STR", "PAL", "TH", "HY", "MB"]
DIVISION_ORDER += ["P", "MY", "CB"]


# ===== Shared pieces =====


def heading(fig: plt.Figure, number: int, question: str, numbers: str) -> None:
    """The figure's number and question, and under it the line of computed numbers."""
    fig.text(0.5, 0.985, f"{number}.  {question}", ha="center", va="top", fontsize=12)
    fig.text(0.5, 0.955, numbers, ha="center", va="top", fontsize=9, color=DARK_GREY)


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


def paint(lab: np.ndarray, value: dict[str, float], names: dict[int, str]) -> np.ndarray:
    """A plane of structure values: each pixel its structure's value, NaN if none."""
    out = np.full(lab.shape, np.nan, dtype=np.float32)
    for idx in np.unique(lab):
        if idx == 0:
            continue
        v = value.get(names.get(int(idx), ""))
        if v is not None:
            out[lab == idx] = v
    return out


def overlay(ax: plt.Axes, mask: np.ndarray, colour: str, lab: np.ndarray) -> None:
    """Paint `mask` flat in `colour` over a plane, and the area borders on top again."""
    rgba = np.zeros(mask.shape + (4,))
    rgba[mask] = plt.matplotlib.colors.to_rgba(colour)
    ax.imshow(rgba, origin="upper", interpolation="nearest", aspect="equal")
    ov = np.zeros(lab.shape + (4,))
    ov[boundaries(lab)] = (0.8, 0.8, 0.8, 0.55)
    ax.imshow(ov, origin="upper", interpolation="nearest", aspect="equal")


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


def saved(fig: plt.Figure, path: Path | None) -> plt.Figure:
    """Save the figure as PNG and EPS when a path is given; return it."""
    if path is not None:
        save_figure(fig, Path(path), DPI)
    return fig


# ===== 01 Which structures, and why =====


def funnel_panel(ax: plt.Axes, set_table: pd.DataFrame, min_adults: int) -> None:
    """A: the rule as a funnel, from the adult table to the declared set."""
    n_table = len(set_table)
    n_all = int((set_table["n_adults"] >= min_adults).sum())
    n_set = int(set_table["in_set"].sum())
    steps = [
        ("in the adult table\n(measured in at least one adult)", n_table, LIGHT_GREY),
        (f"measured in all {min_adults} adults", n_all, MID_GREY),
        ("and grey matter: the declared set", n_set, DARK_GREY),
    ]
    for i, (label, n, colour) in enumerate(steps):
        ax.barh(i, n, color=colour, height=0.6)
        ax.text(n + 4, i, f"{n}", va="center", fontsize=9)
    ax.set_yticks(range(len(steps)))
    ax.set_yticklabels([s[0] for s in steps], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, n_table * 1.15)
    ax.set_xlabel("structures")
    panel_title(ax, "A", "The rule, applied in order")
    tidy(ax)


def division_rows(set_table: pd.DataFrame) -> list[tuple[str, pd.DataFrame]]:
    """The structures of each grey division, then every other division pooled."""
    rows = []
    for division in DIVISION_ORDER:
        rows.append((division, set_table[set_table["division"] == division]))
    other = set_table[~set_table["division"].isin(DIVISION_ORDER)]
    rows.append(("fibre tracts,\nventricles, other", other))
    return rows


def division_panel(ax: plt.Axes, set_table: pd.DataFrame, min_adults: int) -> None:
    """B: kept and dropped structures per division, dropped split by reason."""
    few = f"measured in fewer than {min_adults} adults"
    for i, (division, sub) in enumerate(division_rows(set_table)):
        kept = int(sub["in_set"].sum())
        n_few = int((sub["reason"] == few).sum())
        n_catch = int(sub["reason"].str.startswith("catch-all").sum())
        n_white = len(sub) - kept - n_few - n_catch
        group = DIVISION_GROUP.get(division)
        colour = DIVISION_GROUP_COLOURS[group] if group else LIGHT_GREY
        left = 0
        segments = [
            (kept, colour, None),
            (n_few, MID_GREY, None),
            (n_white, LIGHT_GREY, None),
            (n_catch, "white", "////"),
        ]
        for n, c, hatch in segments:
            if n:
                ax.barh(
                    i,
                    n,
                    left=left,
                    color=c,
                    height=0.65,
                    hatch=hatch,
                    edgecolor=MID_GREY if hatch else c,
                    linewidth=0.5 if hatch else 0,
                )
            left += n
        ax.text(left + 0.6, i, f"{kept} of {len(sub)}", va="center", fontsize=7.5)
    labels = [d for d, _ in division_rows(set_table)]
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("structures")
    handles = [
        Patch(color=DIVISION_GROUP_COLOURS["cortex"], label="kept (coloured by group)"),
        Patch(color=MID_GREY, label=few),
        Patch(color=LIGHT_GREY, label="not grey matter"),
        Patch(
            facecolor="white", edgecolor=MID_GREY, hatch="////", label="catch-all label"
        ),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=7.5)
    panel_title(ax, "B", "Kept and dropped, by division")
    tidy(ax)


def plane_panel(
    fig: plt.Figure,
    ax: plt.Axes,
    set_table: pd.DataFrame,
    profile: pd.DataFrame,
    lab: np.ndarray,
    names: dict[int, str],
    plane: int,
) -> None:
    """C: the adult map as ranks on one plane, the declared structures only."""
    in_set = profile[profile["in_set"]]
    rank = dict(zip(in_set.index, ranks01(in_set["zref_nano"].to_numpy())))
    image = draw_plane(ax, paint(lab, rank, names), lab, hot_cut(), RANK_FLOOR, 1.0)

    # the structures left out, flat grey
    dropped = set(set_table.loc[~set_table["in_set"], "structure"])
    mask = np.isin(lab, [i for i, n in names.items() if n in dropped])
    overlay(ax, mask, NO_DATA_GREY, lab)
    colour_bar(
        fig, ax, image, "rank among the declared structures (0 low, 1 high)", (0, 1)
    )
    n_shown = len({names.get(int(i)) for i in np.unique(lab)} & set(rank))
    panel_title(
        ax,
        "C",
        f"The adult nano map, CCF plane {plane}: one rank per declared structure",
        f"{n_shown} declared structures in this plane; flat grey: left out",
    )


def reference_panel(ax: plt.Axes, reference: pd.DataFrame, column: str, label: str):
    """One half of D: each adult's zero or spread, before and after, joined per mouse."""
    nano = reference[reference["channel"] == "nano"]
    for _, r in nano.iterrows():
        before, after = r[f"{column}_stored"], r[column]
        colour = GROUP_COLOURS[r["group"]]
        ax.plot([0, 1], [before, after], color=PAIR_LINE, lw=0.8, zorder=1)
        ax.scatter([0, 1], [before, after], s=28, color=colour, zorder=2, linewidths=0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["17 brains'\nshared set", "declared\nset"], fontsize=7.5)
    ax.set_xlim(-0.4, 1.4)
    ax.set_ylabel(label)
    tidy(ax)


def plot_structures(
    set_table: pd.DataFrame,
    profile: pd.DataFrame,
    reference: pd.DataFrame,
    stored_agreement: float,
    lab: np.ndarray,
    names: dict[int, str],
    plane: int,
    min_adults: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 01: which structures enter every comparison, and what A1 does to zref.

    `profile` is indexed by structure (adult_profile.csv), `reference` holds one
    row per adult and channel (zref_reference.csv), `stored_agreement` is the
    Spearman between the cohort map before and after over the declared set, `lab`
    the CCF plane (DV, ML) of parcellation indices and `names` their structures.
    """
    fig = plt.figure(figsize=(15.5, 10.5))
    n_table = len(set_table)
    n_all = int((set_table["n_adults"] >= min_adults).sum())
    n_set = int(set_table["in_set"].sum())
    heading(
        fig,
        1,
        "Which structures enter every comparison, and how much does the declared "
        "reference move zref?",
        f"{n_table} structures in the adult table, {n_all} measured in all "
        f"{min_adults} adults, {n_set} of those grey matter: the declared set",
    )

    # A and B on the top row, C and D below
    ax_a = fig.add_axes([0.13, 0.66, 0.24, 0.22])
    ax_b = fig.add_axes([0.52, 0.53, 0.42, 0.35])
    ax_c = fig.add_axes([0.03, 0.12, 0.47, 0.42])
    ax_d0 = fig.add_axes([0.60, 0.15, 0.13, 0.26])
    ax_d1 = fig.add_axes([0.82, 0.15, 0.13, 0.26])
    funnel_panel(ax_a, set_table, min_adults)
    division_panel(ax_b, set_table, min_adults)
    plane_panel(fig, ax_c, set_table, profile, lab, names, plane)

    # D: the zero and the spread of each adult, before and after
    nano = reference[reference["channel"] == "nano"]
    shift = (nano["median"] - nano["median_stored"]).abs()
    ratio = nano["spread"] / nano["spread_stored"]
    reference_panel(ax_d0, reference, "median", "zero of zref (log2 over own isocortex)")
    reference_panel(ax_d1, reference, "spread", "spread of zref (p90 - p10, log2)")
    ax_d0.set_title(
        "D.  Each adult's zero and spread of zref, before and after\n"
        f"zero moves by {shift.min():.2f} to {shift.max():.2f}, spread x "
        f"{ratio.min():.2f} to {ratio.max():.2f}; the cohort map's order "
        f"agrees at rho {stored_agreement:.4f}",
        loc="left",
        fontsize=9,
    )
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=GROUP_COLOURS[g], label=g)
        for g in ("naive", "rws")
    ]
    ax_d1.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0))

    footer(
        fig,
        [
            f"How to read: a structure enters when all {min_adults} adults measure it "
            "(at least 250 tissue voxels of 20 um and a nano mean above background) and "
            "it is grey matter; the rule is applied before any correlation. zref's zero "
            "(the median structure) and spread (p90 - p10)",
            "are taken over the declared structures only, so they depend on that brain "
            "and that list, not on which other brains are in a run (before: the "
            "structures all 17 brains share).",
            "What would mean what: a moved zero shifts every structure of a brain "
            "alike and keeps its order, so a cohort map whose order agrees near 1 "
            "means the comparisons downstream change through the set, not through "
            "the reference.",
        ],
    )
    return saved(fig, save)
