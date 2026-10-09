"""The figures of the adult map's analyses: part 1 of the ISH line and the green channel.

The guided figures of part 1, how much of the map Gria1 and synapse density leave
(03, 03s1, 03s2, 03s3, 04, 11s1, 14, 14s1), and of analysis 5, what the green channel
reports (15, 15s), drawn as the ISH line's: number and question at the top, one
line to take from a main figure, grey notes on how to read it at the foot, a
takeaway that follows its numbers; their shared pieces (heading, footer, panel
titles, the scatter of structures by group of divisions) are those of ish.plotting,
so the walk reads as one set, and they are saved as PNG and EPS at its dpi. Beside
them, the working figures of the steps (fig0 to fig6, E_regression, F_maps,
sep_channel_check.png), PNG only at 200 dpi, which show each step's numbers as it
runs. The modules that compute draw nothing: every function here takes tables and
returns the figure.

Called by the run scripts of steps 21 to 29 and by adult.beyond_figures.
"""

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle
from scipy.stats import rankdata, spearmanr

from sepmap.adult import synaptome
from sepmap.adult.beyond_calibration import FLOOR_MAP
from sepmap.ish.figure_index import figure_ref
from sepmap.ish.gene_sets import AMPA_FAMILY, LEFTOVER_GENE
from sepmap.ish.plotting import (
    DIVISION_ORDER,
    RANK_FLOOR,
    colour_bar,
    footer,
    gene_colour,
    group_handles,
    heading,
    overlay,
    p_text,
    panel_title,
    ranks01,
    saved,
    scatter_groups,
    set_column,
    spread_labels,
)
from sepmap.plotting import (
    AUTO,
    AUTO_DOT,
    DARK_BLUE,
    DARK_GREY,
    DENSITY_BLUE,
    DIVISION_GROUP,
    LIGHT_GREY,
    MID_GREY,
    NANO,
    NANO_DOT,
    NO_DATA_GREY,
    NULL_BAND,
    PAIR_LINE,
    RED,
    SEP,
    SEP_DOT,
    SEP_REMAINDER,
    SET_COLOURS,
    bars_grey,
    draw_plane,
    hot_cut,
    paint,
    save_figure,
    tidy,
)

# dpi of the working figures, which only the run's log points to
WORKING_DPI = 200


# ===== Shared pieces =====


def saved_working(fig: plt.Figure, path: Path | None) -> plt.Figure:
    """Save a working figure as PNG at WORKING_DPI when a path is given; return it."""
    if path is not None:
        save_figure(fig, Path(path), WORKING_DPI, eps=False)
    return fig


# ===== 14 and 14s1 The measured synapse density (adult.synaptome) =====


# how each density is drawn in the agreement panel: the one the model uses filled, in
# the density colour of the budget; the variants open, in greys
DENSITY_MARKERS = {
    "psd95": ("o", DENSITY_BLUE, 46, "PSD95 puncta (the one used)"),
    "psd95_only": ("s", DARK_GREY, 22, "PSD95 alone"),
    "sap102": ("^", DARK_GREY, 24, "SAP102 puncta"),
    "all_puncta": ("D", MID_GREY, 18, "every punctum"),
}

# the terms of the agreement panel, top to bottom, and how they are named
AGREEMENT_TERMS = {
    "density": "density term",
    "Gria1": "Gria1 mRNA",
    "nano": "nano map",
    "autofluorescence": "autofluorescence map",
}


def coverage_panel(ax: plt.Axes, coverage: pd.DataFrame) -> None:
    """A: per division, the structures of the fit with a measured density and without."""
    rows = coverage[coverage["fit"] > 0].set_index("division")
    order = [d for d in DIVISION_ORDER if d in rows.index]
    y = np.arange(len(order))
    measured = rows.loc[order, "fit_measured"].to_numpy()
    total = rows.loc[order, "fit"].to_numpy()
    ax.barh(y, measured, color=DENSITY_BLUE, height=0.7, label="measured")
    ax.barh(
        y,
        total - measured,
        left=measured,
        color=LIGHT_GREY,
        height=0.7,
        label="not sampled",
    )
    for k in range(len(order)):
        ax.text(
            total[k] + 0.4,
            k,
            f"{measured[k]} of {total[k]}",
            va="center",
            fontsize=7,
            color=DARK_GREY,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("structures of the fit", fontsize=8)
    ax.set_xlim(0, total.max() * 1.3)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    tidy(ax)


def rank_scatter(
    ax: plt.Axes,
    x: pd.Series,
    y: pd.Series,
    groups: dict[str, str],
    labels: tuple[str, str],
) -> tuple[int, float]:
    """Two maps over the structures both have, as ranks; returns (n, Spearman)."""
    pair = pd.DataFrame(dict(x=x, y=y)).dropna()
    group = [groups.get(s, "other grey matter") for s in pair.index]
    scatter_groups(ax, rankdata(pair["x"]), rankdata(pair["y"]), group)
    ax.set_xlabel(labels[0], fontsize=8)
    ax.set_ylabel(labels[1], fontsize=8)
    tidy(ax)
    return len(pair), float(spearmanr(pair["x"], pair["y"]).statistic)


def agreement_dots(ax: plt.Axes, agreement: pd.DataFrame) -> None:
    """D: each density's Spearman with each term, over the fit's measured structures."""
    fit = agreement[agreement["structures"] == "fit"]
    terms = [t for t in AGREEMENT_TERMS if t in set(fit["term"])]
    for k, term in enumerate(terms):
        ax.axhline(k, color="0.92", lw=0.6, zorder=0)
        for density, (marker, colour, size, _) in DENSITY_MARKERS.items():
            row = fit[(fit["term"] == term) & (fit["density"] == density)]
            filled = density == "psd95"
            ax.scatter(
                row["rho"],
                [k] * len(row),
                marker=marker,
                s=size,
                facecolors=colour if filled else "none",
                edgecolors=colour,
                linewidths=0 if filled else 1.0,
                zorder=3 if filled else 2,
            )
    ax.axvline(0, color=MID_GREY, lw=0.6)
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels([AGREEMENT_TERMS[t] for t in terms], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(-0.2, 1.0)
    ax.set_xlabel("Spearman rho over the measured structures of the fit", fontsize=8)
    handles = []
    for density, (marker, colour, size, label) in DENSITY_MARKERS.items():
        face = colour if density == "psd95" else "none"
        handles.append(
            plt.Line2D(
                [],
                [],
                ls="",
                marker=marker,
                ms=np.sqrt(size),
                mfc=face,
                mec=colour,
                label=label,
            )
        )
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="upper left")
    tidy(ax)


# the rows of figure 14 C, top to bottom: the main model on all its structures, then
# the two check rows that meet on the structures PSD95 covers (check_rows.csv)
PSD95_ROWS = (
    ("main", "the main model,\nall its structures"),
    ("psd95_main", "the main model,\nwhere PSD95 is measured"),
    ("psd95", "PSD95 puncta as the\ndensity term, there"),
)


def synaptome_groups(table: pd.DataFrame) -> dict[str, str]:
    """{structure: group of divisions} of the synaptome's density table."""
    return {
        s: DIVISION_GROUP.get(d, "other grey matter")
        for s, d in table["division"].items()
    }


def coverage_counts(density: pd.DataFrame) -> tuple[int, int]:
    """The structures of the fit, and those of them with a measured density."""
    fit = density[density["in_fit"]]
    return len(fit), int(fit["measured"].sum())


def points_interval(contrast: tuple[float, tuple[float, float]]) -> str:
    """A difference of shares in points with its 95% interval: '+7 (-11 to +25)'."""
    point, (lo, hi) = contrast
    return f"{100 * point:+.0f} ({100 * lo:+.0f} to {100 * hi:+.0f})"


def density_rows_panel(ax: plt.Axes, rows: pd.DataFrame) -> None:
    """C of figure 14: what each density term leaves, and that minus its own floor.

    `rows` holds the rows of PSD95_ROWS by key (the main model's and the two check
    rows of check_rows.csv): the share left as a bar, nano minus the row's floor as
    a dot with its 95% interval.
    """
    keys = [key for key, _ in PSD95_ROWS if key in rows.index]
    y = np.arange(len(keys))
    left = rows.loc[keys, "left"].to_numpy(float)
    margin = rows.loc[keys, "nano_minus_floor"].to_numpy(float)
    lo = rows.loc[keys, "nano_minus_floor_lo"].to_numpy(float)
    hi = rows.loc[keys, "nano_minus_floor_hi"].to_numpy(float)
    ax.barh(y - 0.17, left, height=0.3, color=RED, label="the model leaves")
    ax.errorbar(
        margin,
        y + 0.17,
        xerr=[margin - lo, hi - margin],
        fmt="o",
        ms=5,
        color=DARK_GREY,
        elinewidth=1.2,
        capsize=2,
        label="nano minus its floor (95%)",
    )
    for k in range(len(keys)):
        ax.text(left[k] + 0.01, y[k] - 0.17, f"{left[k]:.0%}", va="center", fontsize=7)
        ax.text(
            hi[k] + 0.01,
            y[k] + 0.17,
            f"{100 * margin[k]:+.0f} points",
            va="center",
            fontsize=7,
            color=DARK_GREY,
        )
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([label for key, label in PSD95_ROWS if key in keys], fontsize=8)
    ax.set_ylim(len(keys) - 0.4, -0.7)
    low = min(0.0, float(lo.min()) - 0.05)
    ax.set_xlim(low, max(0.7, float(np.max(np.maximum(left, hi))) + 0.1))
    ticks = np.arange(np.ceil(low * 10) / 10, ax.get_xlim()[1] + 1e-9, 0.1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    tidy(ax)


def plot_synaptome(
    density: pd.DataFrame,
    coverage: pd.DataFrame,
    density_term: pd.Series,
    rows: pd.DataFrame,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 14: the measured synapse density, against the density genes, in the model.

    `density` and `coverage` are the tables of run_synaptome, `density_term` the
    main model's density term per structure of the fit, `rows` the rows of
    PSD95_ROWS by key: the main model's numbers and those of check_rows.csv, each
    with its share left and nano minus its floor with the 95% interval.
    """
    table = density.set_index("structure")
    n_fit, n_measured = coverage_counts(table)
    measured = table[table["in_fit"] & table["measured"]]
    psd95, there = rows.loc["psd95"], rows.loc["psd95_main"]
    fig = plt.figure(figsize=(16, 6.4))
    heading(
        fig,
        "synaptome",
        f"PSD95 puncta are measured in {n_measured} of the {n_fit} structures of the "
        f"fit; as the density term there the model leaves {psd95['left']:.0%}, with the "
        f"density genes {there['left']:.0%}; nano minus each one's floor "
        f"{100 * psd95['nano_minus_floor']:+.0f} and "
        f"{100 * there['nano_minus_floor']:+.0f} points",
    )
    ax = fig.add_axes([0.07, 0.17, 0.24, 0.57])
    coverage_panel(ax, coverage)
    panel_title(
        ax,
        "A",
        "Structures of the fit with a measured density",
        f"{n_measured} of {n_fit} ({n_measured / n_fit:.0%})",
    )
    ax = fig.add_axes([0.39, 0.17, 0.22, 0.57])
    n, rho = rank_scatter(
        ax,
        measured["psd95"],
        density_term,
        synaptome_groups(table),
        ("PSD95 puncta (rank)", "density term: mean rank of the genes (rank)"),
    )
    ax.legend(handles=group_handles(), fontsize=7, frameon=False, loc="upper left")
    panel_title(
        ax,
        "B",
        "PSD95 puncta against the density term",
        f"rho {rho:+.2f} on {n} structures (the genes were chosen by it,\nso this is "
        f"optimistic; held out: {figure_ref('density_markers')})",
    )
    ax = fig.add_axes([0.74, 0.17, 0.23, 0.57])
    density_rows_panel(ax, rows)
    panel_title(
        ax,
        "C",
        "Each density term in the model, with its floor",
        f"on the {int(psd95['n_structures'])} structures PSD95 covers; nano minus "
        "floor:\n"
        f"PSD95 {points_interval(row_margin(psd95))}, the genes "
        f"{points_interval(row_margin(there))}",
    )
    footer(
        fig,
        [
            "How to read: PSD95 punctum density from one adult mouse (Zhu et al. 2018, "
            "as shared by Hansen et al.), placed in the CCF structures; the main model's "
            "density term is the mean rank of three postsynaptic genes chosen by it.",
            "Each row's floor is what its own model leaves of a map made only of its "
            "terms, predicted from other Allen experiments; PSD95 is measured once, so "
            "its floor misses its mismatch and errs low. In detail: "
            f"{figure_ref('synaptome_detail')}.",
        ],
    )
    return saved(fig, save)


def row_margin(row: pd.Series) -> tuple[float, tuple[float, float]]:
    """A row's nano minus its floor with the 95% interval, as points_interval takes it."""
    return (
        float(row["nano_minus_floor"]),
        (float(row["nano_minus_floor_lo"]), float(row["nano_minus_floor_hi"])),
    )


def plot_synaptome_detail(
    density: pd.DataFrame,
    coverage: pd.DataFrame,
    agreement: pd.DataFrame,
    density_term: pd.Series,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 14s1: the measured synapse density, its coverage, how it compares.

    A: per division, the structures of analysis 4's fit with a measured PSD95 density.
    B: that density against the main model's density term, as ranks. C: the left
    hemisphere against the right, the one check of a one-mouse map. D: each
    density's Spearman with the density term, Gria1, the nano map and
    autofluorescence. `density`, `coverage` and `agreement` are the tables of
    run_synaptome; `density_term` is the density term per structure of the fit.
    """
    table = density.set_index("structure")
    n_fit, n_measured = coverage_counts(table)
    measured = table[table["in_fit"] & table["measured"]]
    groups = synaptome_groups(table)

    fig = plt.figure(figsize=(11.5, 8.8))
    grid = fig.add_gridspec(
        2, 2, hspace=0.42, wspace=0.34, left=0.12, right=0.97, top=0.86, bottom=0.14
    )

    ax = fig.add_subplot(grid[0, 0])
    coverage_panel(ax, coverage)
    panel_title(
        ax,
        "A",
        "Structures of the fit with a measured density",
        f"{n_measured} of {n_fit} ({n_measured / n_fit:.0%})",
    )

    ax = fig.add_subplot(grid[0, 1])
    n, rho = rank_scatter(
        ax,
        measured["psd95"],
        density_term,
        groups,
        ("PSD95 puncta (rank)", "density term (rank)"),
    )
    ax.legend(handles=group_handles(), fontsize=7, frameon=False, loc="upper left")
    panel_title(
        ax,
        "B",
        "PSD95 puncta against the density term",
        f"rho {rho:+.2f}, n = {n}; the genes were chosen by it",
    )

    ax = fig.add_subplot(grid[1, 0])
    n, rho = rank_scatter(
        ax,
        measured["psd95_left"],
        measured["psd95_right"],
        groups,
        ("left hemisphere (rank)", "right hemisphere (rank)"),
    )
    panel_title(
        ax, "C", "One mouse: left hemisphere against right", f"rho {rho:+.2f}, n = {n}"
    )

    ax = fig.add_subplot(grid[1, 1])
    agreement_dots(ax, agreement)
    panel_title(ax, "D", "How each density agrees with the other maps")

    # the title, the numbers and how to read it
    heading(
        fig,
        "synaptome_detail",
        f"PSD95 puncta cover {n_measured} of the {n_fit} structures of the fit "
        f"({n_measured / n_fit:.0%}): the yardstick of the density genes and a check "
        "row of part 1, not its term",
    )
    footer(
        fig,
        [
            "Zhu et al. 2018, one adult male mouse, as Hansen et al. share it. PSD95 "
            f"density: the mean of the {len(synaptome.PSD95)} subtypes whose puncta hold "
            "PSD95, each divided by its largest value as shared;",
            "never a punctum's intensity or size. A structure is measured when it or one "
            "of its parts was sampled; a region above several structures is never spread "
            "onto",
            "them. Every rho is Spearman over the structures both maps have.",
        ],
    )
    return saved(fig, save)


# ===== 03s3 The synapse-density genes (adult.density_markers) =====

# the eligible genes figure 03s3 B draws, highest agreement first, and the genes whose
# share of the half-split choices C draws
N_RANKED = 30
N_SELECTED = 10

# the composites of figure 03s3 D, top to bottom, and how they are named
COMPOSITE_LABELS = {
    "chosen": "the genes chosen",
    "first_proposal": "Dlg4, Homer1, Camk2a",
    "marker_panel": "the 11 synaptic markers",
    "psd_pc1": "psd_pc1",
}


def gene_ranking_dots(ax: plt.Axes, agreement: pd.DataFrame, n_shown: int) -> None:
    """B: the eligible genes highest with PSD95: chosen, in the pool, or excluded."""
    top = agreement.head(n_shown)
    y = np.arange(len(top))
    chosen = top["chosen"].to_numpy(bool)
    excluded = top["excluded"].to_numpy(bool)
    pool = ~chosen & ~excluded
    rho = top["rho"].to_numpy(float)
    for k in y:
        ax.axhline(k, color="0.94", lw=0.6, zorder=0)
    ax.scatter(rho[pool], y[pool], s=30, color=MID_GREY, linewidths=0, zorder=2)
    ax.scatter(rho[chosen], y[chosen], s=46, color=DENSITY_BLUE, linewidths=0, zorder=3)
    ax.scatter(
        rho[excluded],
        y[excluded],
        s=30,
        facecolors="none",
        edgecolors=RED,
        linewidths=1.0,
        zorder=2,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(top["symbol"], fontsize=7)
    for label, out in zip(ax.get_yticklabels(), excluded):
        if out:
            label.set_color(RED)
    ax.set_ylim(len(top) - 0.4, -0.6)
    ax.set_xlabel("Spearman rho with PSD95 punctum density", fontsize=8)
    handles = [
        plt.Line2D([], [], ls="", marker="o", ms=7, mfc=DENSITY_BLUE, mec="none"),
        plt.Line2D([], [], ls="", marker="o", ms=5.5, mfc=MID_GREY, mec="none"),
        plt.Line2D([], [], ls="", marker="o", ms=5.5, mfc="none", mec=RED),
    ]
    ax.legend(
        handles,
        ["chosen", "in the pool", "excluded: AMPA-linked"],
        fontsize=7,
        frameon=False,
        loc="lower right",
    )
    tidy(ax)


def selection_bars(
    ax: plt.Axes, selection: pd.DataFrame, chosen: list[str], n_shown: int
) -> None:
    """C: the share of the random halves on which each gene is chosen."""
    top = selection.head(n_shown)
    y = np.arange(len(top))
    colours = [DENSITY_BLUE if g in chosen else LIGHT_GREY for g in top["symbol"]]
    share = top["share"].to_numpy(float)
    ax.barh(y, share, color=colours, height=0.7)
    for k in y:
        ax.text(share[k] + 0.01, k, f"{share[k]:.0%}", va="center", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels(top["symbol"], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.08)
    ticks = np.arange(0, 1.01, 0.25)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the halves on which it is chosen", fontsize=8)
    tidy(ax)


def held_out_panel(ax: plt.Axes, comparison: pd.DataFrame) -> None:
    """D: each composite's agreement with PSD95 on held-out halves and the full set."""
    rows = comparison.set_index("composite")
    names = [c for c in COMPOSITE_LABELS if c in rows.index]
    for k, name in enumerate(names):
        r = rows.loc[name]
        colour = DENSITY_BLUE if name == "chosen" else DARK_GREY
        ax.plot([r["rho_held_out_lo"], r["rho_held_out_hi"]], [k, k], color=colour, lw=2)
        ax.scatter(r["rho_held_out_median"], k, s=46, color=colour, zorder=3)
        ax.scatter(
            r["rho_full_set"],
            k,
            s=34,
            marker="D",
            facecolors="none",
            edgecolors=colour,
            linewidths=1.0,
            zorder=3,
        )
    labels = []
    for name in names:
        genes = rows.loc[name, "genes"] if name == "chosen" else ""
        labels.append(f"{COMPOSITE_LABELS[name]}\n{genes}".strip())
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(labels, fontsize=7)

    # room under the last row for the legend
    ax.set_ylim(len(names) + 0.4, -0.5)
    ax.set_xlabel("Spearman rho with PSD95 punctum density", fontsize=8)
    handles = [
        plt.Line2D([], [], color=DARK_GREY, lw=2, marker="o", ms=6),
        plt.Line2D([], [], ls="", marker="D", ms=5.5, mfc="none", mec=DARK_GREY),
    ]
    ax.legend(
        handles,
        ["held-out halves: median and 95% range", "full set (optimistic for the chosen)"],
        fontsize=7,
        frameon=False,
        loc="lower left",
    )
    tidy(ax)


def exclusion_bars(ax: plt.Axes, reasons: pd.DataFrame) -> None:
    """E: the genes that met the rule but were left out, by their first reason, named."""
    y = np.arange(len(reasons))
    counts = reasons["n_genes"].to_numpy(int)
    ax.barh(y, counts, color=RED, height=0.62)
    for k, r in enumerate(reasons.itertuples()):
        ax.text(
            r.n_genes + 0.3,
            k,
            f"{r.n_genes}:  {r.genes.replace(' ', ', ')}",
            va="center",
            fontsize=7.5,
        )
    ax.set_yticks(y)
    ax.set_yticklabels([textwrap.fill(r, 58) for r in reasons["reason"]], fontsize=7.5)
    ax.set_ylim(len(reasons) - 0.4, -0.6)
    ax.set_xlim(0, counts.max() * 4.6)
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)
    tidy(ax)


def plot_density_markers(
    agreement: pd.DataFrame,
    composite: pd.Series,
    psd95: pd.Series,
    groups: dict[str, str],
    selection: pd.DataFrame,
    comparison: pd.DataFrame,
    reasons: pd.DataFrame,
    n_halves: int,
    go_release: str,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 03s3: the synapse-density genes, chosen by PSD95 without the map.

    A: PSD95 punctum density against the chosen genes' mean rank, over the declared
    structures both have. B: the eligible genes highest with PSD95, the AMPA-linked
    ones marked. C: how often each gene is chosen on a random half of the structures.
    D: each composite's agreement with PSD95 on the held-out halves and on the full
    set. E: the AMPA-linked genes, by the reason that left each out. The tables are
    those of run_density_markers; `composite` is the chosen genes' mean rank per
    structure, `groups` each structure's group of divisions, `reasons`
    density_markers.exclusion_reasons.
    """
    chosen = list(agreement.loc[agreement["chosen"], "symbol"])
    rows = comparison.set_index("composite")
    mine = rows.loc["chosen"]
    first = rows.loc["first_proposal"]
    fig = plt.figure(figsize=(11.5, 14))
    grid = fig.add_gridspec(
        2, 2, hspace=0.42, wspace=0.42, left=0.12, right=0.97, top=0.89, bottom=0.38
    )

    ax = fig.add_subplot(grid[0, 0])
    n, rho = rank_scatter(
        ax,
        psd95,
        composite,
        groups,
        ("PSD95 punctum density (rank)", "mean rank of the genes chosen (rank)"),
    )
    ax.legend(handles=group_handles(), fontsize=7, frameon=False, loc="upper left")
    panel_title(
        ax,
        "A",
        "PSD95 puncta against the genes chosen",
        f"rho {rho:+.2f}, n = {n} (full set, optimistic)",
    )

    ax = fig.add_subplot(grid[0, 1])
    gene_ranking_dots(ax, agreement, N_RANKED)
    n_pool = int(agreement["in_pool"].sum())
    n_out = int(agreement["excluded"].sum())
    panel_title(
        ax,
        "B",
        "The genes highest with PSD95",
        f"the top {N_RANKED} of {len(agreement)}: {n_pool} in the pool, {n_out} excluded",
    )

    ax = fig.add_subplot(grid[1, 0])
    selection_bars(ax, selection, chosen, N_SELECTED)
    panel_title(
        ax,
        "C",
        "How often each gene is chosen",
        f"the rule on {n_halves} random halves of the structures",
    )

    ax = fig.add_subplot(grid[1, 1])
    held_out_panel(ax, comparison)
    panel_title(
        ax,
        "D",
        "Agreement with PSD95 on the other half",
        f"n = {int(mine['n_held_out'])} structures per half",
    )

    # E reaches further left than the grid, for the names of the GO terms
    ax = fig.add_axes([0.36, 0.08, 0.61, 0.22])
    exclusion_bars(ax, reasons)
    panel_title(
        ax,
        "E",
        "The genes left out as AMPA-linked, and why",
        f"{int(reasons['n_genes'].sum())} genes that met the rule otherwise, each under "
        "its first reason: the subunits, the AMPA receptor\ncomplex family, the "
        "localisation set, then a GO term (any evidence, or a term below it)",
    )

    # the title, the numbers and how to read it
    heading(
        fig,
        "density_markers",
        f"Chosen without the map: {', '.join(chosen[:-1])} and {chosen[-1]}; their mean "
        f"follows PSD95 at rho {mine['rho_held_out_median']:+.2f} on held-out halves "
        f"({mine['rho_held_out_lo']:+.2f} to {mine['rho_held_out_hi']:+.2f}), against "
        f"{first['rho_held_out_median']:+.2f} for Dlg4, Homer1 and Camk2a",
    )
    footer(
        fig,
        [
            f"The pool: genes annotated in mouse GO (release {go_release}) to the "
            "postsynaptic density or specialization, with two or more usable Allen",
            "experiments, measured in 90% of the structures where Gria1 is; the genes "
            "that place or regulate AMPA receptors are excluded (red), as the surface "
            "side rather than density.",
            "PSD95: punctum density of one adult mouse (Zhu et al. 2018). No nano value "
            "is read.",
        ],
    )
    return saved(fig, save)


# ===== Working figures of analysis 4 (adult.beyond_density, beyond_controls) =====

# the readings of control G as the working figure names them: the map's own zref, the
# same with the young-against-adult tables' 17-brain reference, and the stored
# readings
READING_LABELS = {
    "zref": "zref\n(declared reference)",
    "zref_stored": "zref\n(17-brain reference)",
    "cref": "cref",
    "subref": "subref",
    "ratio": "nano / auto",
    "sepratio": "nano / SEP",
}


def plot_structures_used(rows: pd.DataFrame, save: Path | None = None) -> plt.Figure:
    """fig0: the structures used per division, and those left out per reason.

    `rows` is structures_used.csv (beyond_density.structure_rows).
    """
    fig, axes = plt.subplots(
        1, 2, figsize=(11.0, 3.9), gridspec_kw=dict(width_ratios=[1, 1.25])
    )
    counts = rows[rows["used"]].groupby("division").size().sort_values(ascending=False)
    axes[0].bar(
        range(len(counts)), counts.to_numpy(), color="0.65", edgecolor="0.25", lw=0.5
    )
    axes[0].set_xticks(range(len(counts)))
    axes[0].set_xticklabels(counts.index, fontsize=7, rotation=45, ha="right")
    axes[0].set_ylabel("structures used", fontsize=8)
    axes[0].set_title(
        f"what the fit runs on\n{int(rows['used'].sum())} grey-matter structures",
        fontsize=9,
    )
    tidy(axes[0])

    # the reasons, the declared set's own reasons pooled into one
    why = rows.loc[~rows["used"], "reason"].map(
        lambda r: "not in the declared set" if r.startswith("not in the declared") else r
    )
    why_counts = why.value_counts().sort_values()
    axes[1].barh(
        range(len(why_counts)),
        why_counts.to_numpy(),
        color=RED,
        edgecolor="0.25",
        linewidth=0.5,
        alpha=0.85,
    )
    axes[1].set_yticks(range(len(why_counts)))
    axes[1].set_yticklabels([w[:52] for w in why_counts.index], fontsize=7)
    axes[1].set_xlabel("structures left out", fontsize=8)
    axes[1].set_title(
        "and what it leaves out\ndecided by rule, before any fitting", fontsize=9
    )
    tidy(axes[1])
    fig.tight_layout()
    return saved_working(fig, save)


def plot_ceiling(
    agreement: list[float],
    explainable: float,
    reproducible: bool,
    save: Path | None = None,
) -> plt.Figure:
    """fig1: the two half-cohort maps' agreement over every split, and the ceiling.

    `reproducible` says whether the mean agreement replicates
    (beyond_density.replicates), which the title follows.
    """
    half = float(np.mean(agreement))
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    ax.hist(agreement, bins=25, color="0.7", edgecolor="0.35", linewidth=0.4)
    ax.axvline(half, color=RED, lw=1.8)
    ax.set_xlabel("Spearman between the two half-cohort maps", fontsize=8)
    ax.set_ylabel(f"splits of ten adults ({len(agreement)})", fontsize=8)
    verdict = "is reproducible" if reproducible else "does not reproduce"
    ax.set_title(
        f"Step 1. the map {verdict}\n"
        f"half-cohorts agree at {half:.3f}; ceiling (Spearman-Brown) {explainable:.3f}",
        fontsize=9,
    )
    tidy(ax)
    fig.tight_layout()
    return saved_working(fig, save)


def plot_covariates(
    table: pd.DataFrame, explainable: float, save: Path | None = None
) -> plt.Figure:
    """fig2: each model's cross-validated R2 against the ceiling.

    `table` is variance_partition.csv (beyond_density.partition_table).
    """
    fig, ax = plt.subplots(figsize=(7.8, 4.5))
    ax.barh(np.arange(len(table)), table["cv_r2"], color="0.65", edgecolor="0.25", lw=0.5)
    ax.axvline(explainable, color=RED, lw=1.8)
    ax.annotate(
        f"ceiling {explainable:.3f}\n(the map's own reliability)",
        (explainable, len(table) - 0.4),
        color=RED,
        fontsize=7.5,
        ha="right",
        va="top",
        xytext=(-6, 0),
        textcoords="offset points",
    )
    ax.set_yticks(np.arange(len(table)))
    ax.set_yticklabels(table["model"], fontsize=7.5)
    ax.set_xlim(min(0.0, float(table["cv_r2"].min()) - 0.02), 1.02)
    ax.set_xlabel(
        "variance of the map predicted on held-out structures (cross-validated R2)",
        fontsize=8,
    )
    quoted = float(table.set_index("key").loc["model", "cv_r2"])
    ax.set_title(
        f"Step 2. what Gria1 and synapse density predict\n"
        f"the main model: {quoted / explainable:.1%} of the ceiling, "
        f"{1 - quoted / explainable:.1%} left",
        fontsize=9,
    )
    tidy(ax)
    fig.tight_layout()
    return saved_working(fig, save)


def plot_residual(
    map_agreement: list[float],
    agreement: list[float],
    implied: float,
    table: pd.DataFrame,
    save: Path | None = None,
) -> plt.Figure:
    """fig3: the half-cohort agreements, and the structures most off prediction.

    `map_agreement` and `agreement` are the maps' and the leftovers' agreement per
    split, `implied` the replication the ceiling and the fit alone imply, `table`
    residual_by_structure.csv.
    """
    fig, axes = plt.subplots(
        1, 2, figsize=(11.8, 4.5), gridspec_kw=dict(width_ratios=[1, 1.15])
    )
    bins = np.linspace(
        min(min(map_agreement), min(agreement)) - 0.005,
        max(max(map_agreement), max(agreement)) + 0.005,
        30,
    )
    axes[0].hist(
        map_agreement,
        bins=bins,
        color="0.6",
        edgecolor="0.3",
        lw=0.3,
        alpha=0.85,
        label="the map itself",
    )
    axes[0].hist(
        agreement,
        bins=bins,
        color=RED,
        edgecolor="0.3",
        lw=0.3,
        alpha=0.7,
        label="what is left of it",
    )
    axes[0].axvline(implied, color="0.1", lw=1, ls=(0, (3, 2)))
    axes[0].set_xlabel("half-cohort against half-cohort (Spearman)", fontsize=8)
    axes[0].set_ylabel(f"splits of ten adults ({len(agreement)})", fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False, loc="upper left")
    axes[0].set_title(
        f"Step 3. map {np.mean(map_agreement):.3f}, leftover {np.mean(agreement):.3f}\n"
        f"dashed: {implied:.3f}, what the ceiling and the fit alone imply",
        fontsize=9,
    )
    tidy(axes[0])

    show = pd.concat([table.head(8), table.tail(6)])
    axes[1].barh(
        np.arange(len(show)),
        show["residual"],
        color=[RED if v > 0 else DARK_BLUE for v in show["residual"]],
        edgecolor="0.25",
        linewidth=0.4,
    )
    axes[1].set_yticks(np.arange(len(show)))
    axes[1].set_yticklabels([s[:36] for s in show["structure"]], fontsize=6.8)
    axes[1].invert_yaxis()
    axes[1].axvline(0, color="0.3", lw=0.7)
    axes[1].set_xlabel("nano rank minus predicted rank", fontsize=8)
    axes[1].set_title("Step 4. where it is largest", fontsize=9)
    tidy(axes[1])
    fig.tight_layout()
    return saved_working(fig, save)


def plot_artefacts(
    leftover: np.ndarray,
    xyz: np.ndarray,
    sizes: np.ndarray,
    pairs: list[float],
    naive: np.ndarray,
    rws: np.ndarray,
    passed: dict[str, bool],
    save: Path | None = None,
) -> plt.Figure:
    """fig4: controls A to D, position, size, pairs of adults, naive against RWS.

    `xyz` holds each structure's centroid (AP, DV, ML) in mm, `passed` each
    control's verdict by letter; each title says what its verdict says.
    """
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 3.9))

    # A: the leftover against anterior-posterior position
    axes[0].scatter(
        xyz[:, 0], leftover, s=12, facecolor="0.6", edgecolor="0.25", linewidth=0.3
    )
    axes[0].axhline(0, color="0.85", lw=0.7)
    axes[0].set_xlabel("structure centroid, anterior-posterior (mm)", fontsize=8)
    axes[0].set_ylabel("residual (ranks)", fontsize=8)
    if passed["A"]:
        title = "A. not a front-to-back gradient"
    else:
        title = "A. a gradient could explain it"
    axes[0].set_title(title, fontsize=9)
    tidy(axes[0])

    # B: the leftover against structure volume
    good = np.isfinite(sizes)
    axes[1].scatter(
        np.log10(sizes[good]),
        leftover[good],
        s=12,
        facecolor="0.6",
        edgecolor="0.25",
        linewidth=0.3,
    )
    axes[1].axhline(0, color="0.85", lw=0.7)
    axes[1].set_xlabel("log10 structure volume (20 um voxels)", fontsize=8)
    axes[1].set_ylabel("residual (ranks)", fontsize=8)
    title = "B. not small-structure noise" if passed["B"] else "B. size may drive it"
    axes[1].set_title(title, fontsize=9)
    tidy(axes[1])

    # C: the agreement of every pair of adults, the median in red
    axes[2].hist(pairs, bins=20, color="0.7", edgecolor="0.35", linewidth=0.4)
    axes[2].axvline(float(np.median(pairs)), color=RED, lw=1.8)
    axes[2].set_xlabel("leftover of one adult against another (Spearman)", fontsize=8)
    axes[2].set_ylabel("pairs of adults", fontsize=8)
    title = "C. every animal shows it" if passed["C"] else "C. one animal may carry it"
    axes[2].set_title(title, fontsize=9)
    tidy(axes[2])

    # D: the naive leftover against the RWS one, with the identity line
    axes[3].scatter(naive, rws, s=12, facecolor="0.6", edgecolor="0.25", linewidth=0.3)
    lim = [min(naive.min(), rws.min()) - 3, max(naive.max(), rws.max()) + 3]
    axes[3].plot(lim, lim, color="0.75", ls="--", lw=0.8)
    axes[3].set_xlabel("leftover, five naive adults", fontsize=8)
    axes[3].set_ylabel("leftover, five RWS adults", fontsize=8)
    title = (
        "D. not the whisker manipulation" if passed["D"] else "D. naive and RWS disagree"
    )
    axes[3].set_title(f"{title}\nrho {spearmanr(naive, rws).statistic:+.2f}", fontsize=9)
    tidy(axes[3])

    not_ruled_out = [k for k in "ABCD" if not passed[k]]
    if not_ruled_out:
        title = (
            "Four ways the leftover could be an artefact; not ruled out: "
            + ", ".join(not_ruled_out)
        )
    else:
        title = "Four ways the leftover could be an artefact, and is not"
    fig.suptitle(title, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return saved_working(fig, save)


def plot_model_space(
    curve: pd.DataFrame,
    k_most: int,
    explainable: float,
    scores: tuple[float, float, float],
    passed: dict[str, bool],
    save: Path | None = None,
) -> plt.Figure:
    """fig5: controls F and E, the gene-space curve and curving the model.

    `curve` is gene_space.csv, `k_most` the number of components picked most often,
    `scores` the held-out R2 of the model straight, curved to cubes and to fifth
    powers.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ks = curve["n_components"]
    axes[0].plot(ks, curve["r2"], color="0.6", lw=1.5, label="fitted")
    axes[0].plot(ks, curve["cv_r2"], color=RED, lw=1.8, label="cross-validated")
    axes[0].axhline(explainable, color="0.3", ls="--", lw=1.2)
    axes[0].annotate(
        "ceiling", (ks.iloc[-1], explainable), fontsize=7.5, ha="right", va="bottom"
    )
    axes[0].axvline(k_most, color="0.4", ls=":", lw=1.0)
    axes[0].set_xlabel(
        "components of the expression of the genes measured in every structure",
        fontsize=8,
    )
    axes[0].set_ylabel("variance of the map explained", fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False)
    if passed["F"]:
        title = "F. even the whole gene table falls short"
    else:
        title = "F. the whole gene table accounts for the map"
    axes[0].set_title(
        f"{title}\nthe gap between the two lines is overfitting", fontsize=9
    )
    tidy(axes[0])

    axes[1].bar(
        [0, 1, 2], scores, color=[RED, "0.65", "0.8"], edgecolor="0.25", linewidth=0.5
    )
    axes[1].set_xticks([0, 1, 2])
    axes[1].set_xticklabels(
        ["the model\n(straight)", "curved\n(squares, cubes)", "to fifth powers"],
        fontsize=8,
    )
    axes[1].set_ylabel("cross-validated R2", fontsize=8)
    if passed["E"]:
        title = "E. and curving it buys nothing"
    else:
        title = "E. and curving it pays: read the check row curved"
    axes[1].set_title(title, fontsize=9)
    tidy(axes[1])
    fig.tight_layout()
    return saved_working(fig, save)


def plot_readings(
    table: pd.DataFrame, passed: bool, save: Path | None = None
) -> plt.Figure:
    """fig6: control G, per reading the map's and the leftover's replication, and R2.

    `table` is readings.csv.
    """
    fig, ax = plt.subplots(figsize=(8.0, 4.0))
    x = np.arange(len(table))
    ax.bar(
        x - 0.22,
        table["map_replication"],
        width=0.2,
        color="0.55",
        edgecolor="0.25",
        linewidth=0.4,
        label="the map replicates",
    )
    ax.bar(
        x,
        table["leftover_replication"],
        width=0.2,
        color=RED,
        edgecolor="0.25",
        linewidth=0.4,
        label="the leftover replicates",
    )
    ax.bar(
        x + 0.22,
        table["r2"],
        width=0.2,
        color="0.8",
        edgecolor="0.25",
        linewidth=0.4,
        label="the model explains (R2)",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([READING_LABELS[r] for r in table["reading"]], fontsize=7.5)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    if passed:
        title = "G. the same picture under every reading, not just zref"
    else:
        title = "G. not the same picture under every reading"
    ax.set_title(title, fontsize=9)
    tidy(ax)
    fig.tight_layout()
    return saved_working(fig, save)


# ===== Working figures of the regression (adult.beyond_regression) =====


def draw_fit(
    ax: plt.Axes,
    observed: np.ndarray,
    predicted: np.ndarray,
    leftover: np.ndarray,
    structures: list[str],
    fitted_r2: float,
    cv: float,
) -> None:
    """Draw observed against predicted, the five largest residuals named."""
    lim = [
        min(observed.min(), predicted.min()) - 4,
        max(observed.max(), predicted.max()) + 4,
    ]
    ax.plot(lim, lim, color=MID_GREY, ls="--", lw=1.0, zorder=1)
    ax.scatter(
        predicted,
        observed,
        s=16,
        facecolor=DARK_GREY,
        edgecolor="0.2",
        linewidth=0.3,
        zorder=2,
    )
    worst = np.argsort(-np.abs(leftover))[:5]
    for i in worst:
        ax.annotate(
            structures[i][:24],
            (predicted[i], observed[i]),
            fontsize=6.5,
            color=RED,
            xytext=(4, 2),
            textcoords="offset points",
        )
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("predicted from Gria1 and synapse density (rank)", fontsize=8.5)
    ax.set_ylabel("nano map (rank)", fontsize=8.5)
    ax.set_title(f"the fit\nR2 {fitted_r2:.2f} fitted, {cv:.2f} predicted", fontsize=9.5)
    tidy(ax)


def draw_diagnostic(
    ax: plt.Axes,
    predicted: np.ndarray,
    leftover: np.ndarray,
    misspecified: bool,
    p_max: float,
) -> None:
    """Draw the residual against the prediction, the standard diagnostic.

    The title follows `misspecified`, whether a tilt or a fan reached `p_max`.
    """
    ax.axhline(0, color=MID_GREY, lw=1.0)
    ax.scatter(
        predicted,
        leftover,
        s=16,
        facecolor=DARK_GREY,
        edgecolor="0.2",
        linewidth=0.3,
        zorder=2,
    )
    ax.set_xlabel("predicted (rank)", fontsize=8.5)
    ax.set_ylabel("residual (ranks)", fontsize=8.5)
    if misspecified:
        title = (
            f"the diagnostic\na tilt or a fan (p < {p_max:g}), so the model\n"
            "may be mis-specified"
        )
    else:
        title = (
            "the diagnostic\nno tilt and no fan, so the model is\n"
            "incomplete rather than mis-specified"
        )
    ax.set_title(title, fontsize=9.5)
    tidy(ax)


def draw_leftover(ax: plt.Axes, leftover: np.ndarray) -> None:
    """Draw the distribution of the residuals, the leftover."""
    ax.hist(leftover, bins=26, color=MID_GREY, edgecolor="0.3", linewidth=0.4)
    ax.axvline(0, color=DARK_GREY, lw=1.4)
    ax.set_xlabel("residual (ranks)", fontsize=8.5)
    ax.set_ylabel("structures", fontsize=8.5)
    ax.set_title(
        f"the leftover\nspread {leftover.std():.1f} ranks over {len(leftover)} "
        "structures",
        fontsize=9.5,
    )
    tidy(ax)


def plot_regression(
    observed: np.ndarray,
    predicted: np.ndarray,
    leftover: np.ndarray,
    structures: list[str],
    fitted_r2: float,
    cv: float,
    misspecified: bool,
    p_max: float,
    save: Path | None = None,
) -> plt.Figure:
    """E_regression: observed against predicted, the diagnostic, the residuals.

    `misspecified` says whether a tilt or a fan of the residuals reached `p_max`
    (beyond_regression.diagnostic).
    """
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.3))
    draw_fit(axes[0], observed, predicted, leftover, structures, fitted_r2, cv)
    draw_diagnostic(axes[1], predicted, leftover, misspecified, p_max)
    draw_leftover(axes[2], leftover)
    fig.suptitle(
        "E.  The regression of analysis 4: the nano map predicted from receptor "
        "mRNA and synaptic density",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return saved_working(fig, save)


def draw_map(
    fig: plt.Figure,
    ax: plt.Axes,
    img: np.ndarray,
    cmap: str,
    span: float | None,
    floor: float,
    n_structures: int,
    labels: tuple[str | None, str | None, bool],
) -> None:
    """Draw one map of F_maps, with its column's colour bar under the last plane.

    `span` is the symmetric limit of the leftover, None for a map of ranks, whose
    scale starts `floor` of the ranks below rank 1 so that nothing draws as black.
    `labels` holds the column's title (first plane), the plane's label (first
    column), and whether the colour bar goes under it (last plane).
    """
    title, ylabel, with_bar = labels
    if span is None:
        lo = 1 - floor * (n_structures - 1)
        im = ax.imshow(
            img, cmap=cmap, vmin=lo, vmax=n_structures, interpolation="nearest"
        )
    else:
        im = ax.imshow(img, cmap=cmap, vmin=-span, vmax=span, interpolation="nearest")

    # the image rasterised in the EPS, the text kept vector; NaN is
    # transparent, so the black face is the ground
    im.set_rasterized(True)
    ax.set_facecolor("black")
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)
    if title is not None:
        ax.set_title(title, fontsize=9)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=8)
    if with_bar:
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, orientation="horizontal")
        cb.ax.tick_params(labelsize=6.5)
        if span is None:
            label = "rank among structures"
        else:
            label = "observed minus predicted (ranks)"
        cb.set_label(label, fontsize=7)
        if span is None:
            # the floor sits below rank 1 so that nothing draws as black;
            # the ticks still stop at the real range
            cb.set_ticks([1, 25, 50, 75, 100, n_structures])


def plot_regression_maps(
    planes: list[dict],
    observed: np.ndarray,
    predicted: np.ndarray,
    leftover: np.ndarray,
    structures: list[str],
    floor: float,
    save: Path | None = None,
) -> plt.Figure:
    """F_maps: observed, predicted and residual maps on coronal planes.

    `planes` holds per plane its parcellation indices (lab), its number in the full
    CCF at 10 um (ccf_plane) and the structure names of the indices (names).
    """
    # (title, value per structure, colormap, symmetric limit or None for ranks)
    maps = [
        ("the nano map", dict(zip(structures, observed)), "hot", None),
        (
            "predicted from Gria1\nand synapse density",
            dict(zip(structures, predicted)),
            "hot",
            None,
        ),
        (
            "what is left over",
            dict(zip(structures, leftover)),
            "RdBu_r",
            float(np.abs(leftover).max()),
        ),
    ]
    fig, axes = plt.subplots(len(planes), 3, figsize=(10.2, 3.2 * len(planes)))
    axes = np.atleast_2d(axes)
    for r, plane in enumerate(planes):
        # the plane in the full CCF at 10 um, as the route's other figures give it
        ylabel = f"CCF plane {plane['ccf_plane']} / 10 um"
        for c, (title, values, cmap, span) in enumerate(maps):
            img = paint(plane["lab"], values, plane["names"])
            labels = (
                title if r == 0 else None,
                ylabel if c == 0 else None,
                r == len(planes) - 1,
            )
            draw_map(fig, axes[r, c], img, cmap, span, floor, len(structures), labels)
    fig.suptitle(
        "F.  The same three quantities on the brain.  Red in the third column is "
        "a higher nano rank than\nGria1 and synapse density predict, blue is lower.  "
        "Black is outside the brain, or a structure the\nfit does not use.",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return saved_working(fig, save)


# ===== 03, 03s1 and 03s2 How much of the map Gria1 and synapse density predict =====

# the parts of the reproducible map: Gria1 alone in the subunits' dark blue, what the
# two terms share in mid grey, synapse density alone in the pale blue of
# DENSITY_BLUE, what is left in red; control F's components in mid grey
PART_COLOURS = {
    "gria1_only": DARK_BLUE,
    "shared": MID_GREY,
    "density_only": DENSITY_BLUE,
    "left": RED,
    "components": MID_GREY,
}

# how each part is named on a bar
PART_LABELS = {
    "gria1_only": "Gria1 only",
    "shared": "shared",
    "density_only": "density only",
    "left": "left",
}

# the rows of the check-row panel of figure 03s1, by their key in folds.csv (folds)
# or check_rows.csv (checks), short; the number of structures follows a row on other
# structures than the main model's
ROW_LABELS = {
    "main": "the main model",
    "one_shuffling": "one shuffling of the folds",
    "single_shufflings": "single shufflings (95%)",
    "ten_folds": "ten folds",
    "leave_one_out": "leave one out",
    "spatial_blocks": "folds of spatial blocks",
    "curved": "curved (x, x², x³)",
    "autofluorescence": "+ autofluorescence",
    "first_proposal": "density: Dlg4, Homer1, Camk2a",
    "marker_panel": "density: 11 marker genes",
    "psd_pc1": "+ psd_pc1",
    "four_subunits": "abundance: Gria1 to Gria4",
    "psd95": "density: PSD95 puncta",
    "psd95_main": "main model where PSD95 is",
    "large": "structures of 0.4 mm³ or more",
    "allen_grid": "nano on the Allen 200 um grid",
}

# structures named on the scatter of the whole model: this many each way
NAMED_LEFTOVER = 3


def map_scatter(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    groups: list[str],
    xlabel: str,
) -> None:
    """The nano map's ranks against a predictor's, 0 to 1, by group, identity dashed."""
    x, y = ranks01(x), ranks01(y)
    scatter_groups(ax, x, y, groups)
    ax.plot([0, 1], [0, 1], color=MID_GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.04, 1.04)
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("nano map, rank (0 low, 1 high)")
    tidy(ax)


def model_scatter(
    ax: plt.Axes,
    y: np.ndarray,
    predicted: np.ndarray,
    residual: np.ndarray,
    groups: list[str],
    acronyms: list[str],
    xlabel: str,
) -> None:
    """The map against the whole model's held-out prediction.

    The structures furthest above and below it are named.
    """
    map_scatter(ax, predicted, y, groups, xlabel)
    order = np.argsort(residual)
    xr, yr = ranks01(predicted), ranks01(y)
    extremes = list(order[:NAMED_LEFTOVER]) + list(order[-NAMED_LEFTOVER:])
    spread_labels(ax, [(xr[i], yr[i], acronyms[i], "0.15") for i in extremes])


def part_segments(parts: dict[str, float]) -> tuple[list, list]:
    """The segments of the partition bar, and its negative parts.

    The positive parts of PART_LABELS laid end to end from zero, (key, start, end);
    a negative part cannot sit in the stack, so it is returned apart, (key, value),
    to be drawn leftwards from zero.
    """
    segments, negative = [], []
    start = 0.0
    for key in PART_LABELS:
        value = parts[key]
        if value < 0:
            negative.append((key, value))
            continue
        segments.append((key, start, start + value))
        start += value
    return segments, negative


def partition_bar(ax: plt.Axes, y: float, parts: dict[str, float], height: float) -> None:
    """One bar of the four parts of the reproducible map, a negative one hatched.

    A negative part (held out, a term that adds nothing to the other costs a little)
    is drawn from zero leftwards, hatched, so the stack of the others then reaches
    past 100% by as much.
    """
    segments, negative = part_segments(parts)
    for key, a, b in segments:
        w = b - a
        ax.barh(
            y,
            w,
            left=a,
            height=height,
            color=PART_COLOURS[key],
            edgecolor="white",
            linewidth=1.5,
            zorder=2,
        )
        if w > 0.06:
            ink = "white" if key in ("gria1_only", "left") else "0.1"
            ax.text(
                a + w / 2,
                y,
                f"{PART_LABELS[key]}\n{parts[key]:.0%}",
                ha="center",
                va="center",
                fontsize=8 if w > 0.1 else 6.5,
                color=ink,
                zorder=4,
            )
    for key, value in negative:
        ax.barh(
            y,
            value,
            left=0,
            height=height * 0.6,
            color="white",
            edgecolor=PART_COLOURS[key],
            hatch="////",
            linewidth=1.0,
            zorder=2,
        )
        ax.text(
            value - 0.01,
            y,
            f"{PART_LABELS[key]} {value:+.1%}",
            ha="right",
            va="center",
            fontsize=7,
            color=DARK_GREY,
        )


def small_parts_note(parts: dict[str, float]) -> str:
    """The parts too narrow to be labelled on the bar, as text: 'shared +2%'."""
    segments, _ = part_segments(parts)
    narrow = [key for key, a, b in segments if b - a <= 0.06]
    return ", ".join(f"{PART_LABELS[k]} {parts[k]:+.1%}" for k in narrow)


def partition_panel(ax: plt.Axes, n: dict) -> None:
    """B of figure 03: the reproducible map in four parts, with where 'left' begins.

    Gria1 only, shared, density only and left, each a share of what two halves of
    the cohort reproduce; under the bar, the 95% interval over structures of where
    the leftover begins; a part too narrow to label is named above the bar.
    """
    parts = n["parts"]
    partition_bar(ax, 0, parts, height=0.6)
    lo, hi = n["parts_ci"]["left"]
    ax.plot([1 - hi, 1 - lo], [-0.48, -0.48], color="0.1", lw=1.4)
    for x in (1 - hi, 1 - lo):
        ax.plot([x, x], [-0.53, -0.43], color="0.1", lw=1.0)
    ax.text(
        1 - hi - 0.01,
        -0.48,
        f"where 'left' begins, 95% over structures (left {lo:.0%} to {hi:.0%})",
        fontsize=7.5,
        va="center",
        ha="right",
    )
    note = small_parts_note(parts)
    if note:
        ax.text(0, 0.36, note, fontsize=7.5, va="bottom", color=DARK_GREY)
    low = min(0.0, min(parts.values()) - 0.12)
    ax.set_xlim(low, 1)
    ticks = np.arange(0, 1.01, 0.2)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_ylim(-0.7, 0.5)
    ax.set_yticks([])
    ax.set_xlabel(
        "share of the map's reproducible pattern, on structures the fit has not seen"
    )
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)


def floor_panel(ax: plt.Axes, calibration: pd.DataFrame, n: dict) -> None:
    """C of figure 03: nano and the floor on the same structures, and their difference.

    On the calibration's structures: a dot per draw of made-up adults whose map is
    only Gria1 and synapse density, the median as a bar; the nano map read the same
    way, the mean over the two halves of the Allen experiments its predictors come
    from; under them, nano minus the floor in points with its 95% interval over
    structures.
    """
    random = calibration[calibration["folds"] == "random"]
    floor = random.loc[random["map"] == FLOOR_MAP, "left"].to_numpy(float)
    rng = np.random.default_rng(0)
    ax.scatter(
        floor,
        rng.uniform(-0.13, 0.13, len(floor)),
        s=14,
        color=DARK_GREY,
        alpha=0.6,
        linewidths=0,
    )
    median = float(np.median(floor))
    ax.plot([median, median], [-0.28, 0.28], color="0.1", lw=2, zorder=3)
    ax.text(median, -0.34, f"{median:.0%}", ha="center", va="bottom", fontsize=8.5)

    # nano, with the gap from the floor's median drawn as an arrow
    left = n["nano_cal_left"]
    ax.plot([median, median], [0.3, 1], color=MID_GREY, lw=0.8, ls=(0, (2, 2)))
    ax.annotate(
        "",
        xy=(left - 0.012, 1),
        xytext=(median, 1),
        arrowprops=dict(arrowstyle="->", color=MID_GREY, lw=1.2),
    )
    ax.scatter([left], [1], s=90, color=RED, linewidths=0, zorder=3)
    ax.text(left, 0.7, f"{left:.0%}", ha="center", va="bottom", fontsize=8.5)

    # the difference, with its interval over structures
    diff = n["minus_floor"]
    lo, hi = n["minus_floor_ci"]
    ax.errorbar(
        [diff],
        [2],
        xerr=[[diff - lo], [hi - diff]],
        fmt="D",
        ms=7,
        color=DARK_GREY,
        elinewidth=1.8,
        capsize=4,
        zorder=3,
    )
    ax.text(
        diff,
        1.68,
        f"{100 * diff:+.0f} points (95% {100 * lo:+.0f} to {100 * hi:+.0f})",
        ha="center",
        va="bottom",
        fontsize=8.5,
    )
    ax.axvline(0, color="0.3", lw=0.7)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(
        [
            "the floor: a map made only\nof Gria1 and synapse density",
            "the nano map",
            "nano minus the floor",
        ],
        fontsize=8.5,
    )
    ax.set_ylim(2.5, -0.75)
    low = min(0.0, lo - 0.05)
    ax.set_xlim(low, max(0.5, float(floor.max()) + 0.05, left + 0.05, hi + 0.05))
    ticks = np.arange(0, ax.get_xlim()[1] + 1e-9, 0.1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map left (the difference in points)")
    tidy(ax)


def replication_panel(ax: plt.Axes, replication: pd.DataFrame, n: dict) -> None:
    """I of figure 03s1: the half-maps' and the leftovers' agreement, every split."""
    rng = np.random.default_rng(1)
    rows = (
        (1, replication["map_agreement"].to_numpy(), DARK_GREY),
        (0, replication["leftover_agreement"].to_numpy(), RED),
    )
    for y, values, colour in rows:
        ax.scatter(
            values,
            y + rng.uniform(-0.2, 0.2, len(values)),
            s=14,
            color=colour,
            alpha=0.6,
            linewidths=0,
            zorder=2,
        )
        ax.plot([values.mean()] * 2, [y - 0.3, y + 0.3], color="0.1", lw=2, zorder=3)
        ax.text(
            values.mean(),
            y + 0.36,
            f"mean {values.mean():.3f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.axvline(n["implied"], color="0.1", lw=1, ls=(0, (3, 2)), zorder=1)
    ax.text(
        n["implied"] - 0.003,
        -0.62,
        f"{n['implied']:.3f}: what the leftover\nreplicates at anyway, given the\n"
        f"map's reliability and the fit (R² {n['r2']:.2f})",
        ha="right",
        va="bottom",
        fontsize=7.5,
    )
    ax.set_yticks([1, 0])
    ax.set_yticklabels(
        ["the two\nhalf-cohort maps", "their two leftovers\n(map minus fit)"]
    )
    ax.tick_params(axis="y", length=0)
    low = min(replication["leftover_agreement"].min(), n["implied"]) - 0.02
    ax.set_xlim(low, 1.0)
    ax.set_ylim(-0.7, 1.6)
    ax.set_xlabel("Spearman rho between the two halves (5 adults against 5)")
    ax.spines["left"].set_visible(False)
    tidy(ax)


def margin_text(n: dict) -> str:
    """Nano minus the floor with its interval: '+14 points (95% -2 to +30)'."""
    lo, hi = n["minus_floor_ci"]
    return (
        f"{100 * n['minus_floor']:+.0f} points (95% {100 * lo:+.0f} to {100 * hi:+.0f})"
    )


def density_check_panel(ax: plt.Axes, markers: dict, groups: dict[str, str]) -> None:
    """D of figure 03: measured PSD95 punctum density against the density term.

    `markers` is adult.density_markers.load_validation: PSD95 and the chosen genes'
    mean rank over the declared structures where both exist; `groups` each
    structure's group of divisions.
    """
    rank_scatter(
        ax,
        markers["psd95"],
        markers["composite"],
        groups,
        (
            "PSD95 punctum density, measured (rank)",
            f"mean rank of {', '.join(markers['chosen'])} (rank)",
        ),
    )
    ax.legend(handles=group_handles(), fontsize=7, frameon=False, loc="upper left")


def density_check_text(markers: dict) -> str:
    """The numbers of figure 03 D: the held-out rho, and the full set's beside it."""
    chosen = markers["comparison"].loc["chosen"]
    first = markers["comparison"].loc["first_proposal"]
    return (
        f"held out, rho {chosen['rho_held_out_median']:+.2f} (95% "
        f"{chosen['rho_held_out_lo']:+.2f} to {chosen['rho_held_out_hi']:+.2f}) over "
        f"{markers['n_halves']} random halves;\nfull set "
        f"{chosen['rho_full_set']:+.2f}, n = {int(chosen['n_structures'])} "
        "(optimistic); Dlg4, Homer1, Camk2a held out "
        f"{first['rho_held_out_median']:+.2f}"
    )


def beyond_takeaway(n: dict) -> str:
    """The line under figure 03's title: what is left, and how it stands to the floor.

    The words follow the interval of nano minus the floor.
    """
    lo, _ = n["minus_floor_ci"]
    floor = n["floor"]["left_median"]
    if lo > 0:
        against = f"above the floor's {floor:.0%} by {margin_text(n)}"
    elif n["minus_floor"] > 0:
        against = (
            f"above the floor's {floor:.0%} at its point value only, by {margin_text(n)}"
        )
    else:
        against = f"no more than the floor's {floor:.0%}: {margin_text(n)}"
    left_lo, left_hi = n["left_ci"]
    return (
        f"Gria1 and synapse density leave {n['left']:.0%} ({left_lo:.0%} to "
        f"{left_hi:.0%}) of the map's reproducible pattern; on the {n['cal_n']} "
        f"structures of the calibration nano leaves {n['nano_cal_left']:.0%}, {against}"
    )


def checks_above(n: dict) -> str:
    """How many check rows stay above their own floor: '10 of 10 check rows'."""
    checks = n["checks"]
    above = int((checks["nano_minus_floor_lo"] > 0).sum())
    return f"{above} of {len(checks)} check rows"


def plot_beyond(
    y: np.ndarray,
    held_out: np.ndarray,
    residual: np.ndarray,
    groups: list[str],
    acronyms: list[str],
    numbers: dict,
    calibration: pd.DataFrame,
    markers: dict,
    structure_groups: dict[str, str],
    save: Path | None = None,
) -> plt.Figure:
    """Figure 03: what Gria1 expression and synapse density predict of the map.

    What they leave is set against the calibration floor, and the density term
    against the measured synapse density that chose it. `y` is the map's ranks on
    the structures of the fit, `held_out` each structure's prediction by the main
    model from fits that never saw it, `residual` its leftover; `numbers` holds the
    numbers of adult.beyond_figures, `calibration` the table of beyond_calibration,
    `markers` adult.density_markers.load_validation and `structure_groups` every
    declared structure's group of divisions.
    """
    n = numbers
    fig = plt.figure(figsize=(15, 11))
    heading(fig, "beyond", beyond_takeaway(n))

    # A: the map against the main model's held-out prediction
    ax = fig.add_axes([0.07, 0.53, 0.27, 0.34])
    model_scatter(
        ax,
        y,
        held_out,
        residual,
        groups,
        acronyms,
        "what Gria1 and synapse density predict (held out), rank",
    )
    ax.legend(handles=group_handles(), loc="lower right", fontsize=7, frameon=False)
    panel_title(
        ax,
        "A",
        "The map against what the model predicts",
        f"{n['n_structures']} structures, each predicted by fits that never saw it;\n"
        "off the diagonal: the leftover",
    )

    # B: the partition
    ax = fig.add_axes([0.45, 0.62, 0.5, 0.17])
    partition_panel(ax, n)
    weights = n["weights"]
    panel_title(
        ax,
        "B",
        "The reproducible map in four parts",
        f"nano ~ {n['abundance']} + synapse density, two straight terms; weights "
        f"(z-scored) {', '.join(f'{k} {v:+.2f}' for k, v in weights.items())}",
    )

    # C: the leftover beside the floor
    ax = fig.add_axes([0.17, 0.15, 0.3, 0.27])
    floor_panel(ax, calibration, n)
    panel_title(
        ax,
        "C",
        "Is the leftover more than Allen-to-Allen mismatch?",
        f"on the {n['cal_n']} structures where both halves of the Allen experiments "
        f"measure every gene;\n{checks_above(n)} also stay above their own floor "
        f"({figure_ref('beyond_budget')})",
    )

    # D: the density term against the measured synapse density
    ax = fig.add_axes([0.62, 0.15, 0.25, 0.27])
    density_check_panel(ax, markers, structure_groups)
    panel_title(
        ax,
        "D",
        "Does the density term follow measured synapse density?",
        density_check_text(markers),
    )
    footer(
        fig,
        [
            "How to read: the map and both terms are ranks across structures; a share "
            "is the R² on structures the fit has not seen, over what two halves of the "
            "cohort reproduce. Gria1 only, density only and shared come from Gria1 "
            "alone, density alone and both.",
            "The floor is a map made only of Gria1 and synapse density, from half of "
            "each gene's Allen experiments, with ten made-up adults as noisy as ours, "
            "predicted from the other half; it holds Allen-to-Allen mismatch, not that "
            "of Allen's P56 mice with these",
            "brains, so it errs low. Synapse density is the mean rank of three "
            "postsynaptic genes, chosen by their agreement with PSD95 punctum density "
            "(Zhu et al. 2018, one adult mouse) by a rule that reads no nano value, the "
            "AMPA-linked genes left out;",
            "held out, the rule chooses on a random half of the structures and is "
            "tested on the other.",
            f"In detail: the check rows, {figure_ref('beyond_budget')}; the seven "
            f"controls, {figure_ref('beyond_controls')}; the choice of the genes, "
            f"{figure_ref('density_markers')}.",
        ],
    )
    return saved(fig, save)


def predictor_scatters(
    fig: plt.Figure,
    y: np.ndarray,
    terms: dict[str, np.ndarray],
    held_out: np.ndarray,
    residual: np.ndarray,
    groups: list[str],
    acronyms: list[str],
    n: dict,
    bottom: float,
) -> None:
    """C to E of figure 03s1: the map against Gria1, synapse density and the model.

    `terms` holds the ranks of Gria1 and of the density term; the model's prediction
    is held out. `bottom` is the panels' lower edge as a share of the figure's height;
    each panel is square.
    """
    width, height = fig.get_size_inches()
    size = 0.2
    shares, ci = n["shares"], n["shares_ci"]
    ax = fig.add_axes([0.06, bottom, size, size * width / height])
    map_scatter(ax, terms["Gria1"], y, groups, "Gria1 mRNA, rank")
    lo, hi = ci["abundance"]
    panel_title(
        ax,
        "C",
        "The nano map against Gria1 mRNA",
        f"rho {spearmanr(y, terms['Gria1']).statistic:+.2f} on "
        f"{n['n_structures']} structures; alone, held out, it predicts\n"
        f"{shares['abundance']:.0%} of the reproducible map (95% {lo:.0%} to {hi:.0%})",
    )
    ax = fig.add_axes([0.39, bottom, size, size * width / height])
    map_scatter(ax, terms["density"], y, groups, "synapse density term, rank")
    lo, hi = ci["density"]
    panel_title(
        ax,
        "D",
        "The nano map against synapse density",
        f"mean rank of {', '.join(n['density_genes'])}; rho "
        f"{spearmanr(y, terms['density']).statistic:+.2f};\nalone "
        f"{shares['density']:.0%} of the reproducible map (95% {lo:.0%} to {hi:.0%})",
    )
    ax = fig.add_axes([0.72, bottom, size, size * width / height])
    model_scatter(
        ax, y, held_out, residual, groups, acronyms, "what the model predicts (held out)"
    )
    lo, hi = ci["model"]
    panel_title(
        ax,
        "E",
        "The nano map against the whole model",
        f"Gria1 and synapse density, straight: {shares['model']:.0%} of the\n"
        f"reproducible map (95% {lo:.0%} to {hi:.0%}); off the diagonal: the leftover",
    )


def shares_panel(ax: plt.Axes, n: dict) -> None:
    """F of figure 03s1: each model's share held out, with its interval if it has one.

    Gria1 alone, density alone and both (95% over structures), then the model curved
    and autofluorescence alone, and control F's gene space.
    """
    shares, ci = n["shares"], n["shares_ci"]
    rows = [
        ("Gria1 alone", shares["abundance"], ci["abundance"], DARK_BLUE),
        ("synapse density alone", shares["density"], ci["density"], DENSITY_BLUE),
        ("both: the main model", shares["model"], ci["model"], DARK_GREY),
        ("both, curved (check row)", shares["curved"], None, MID_GREY),
        ("autofluorescence alone", shares["autofluorescence"], None, AUTO),
        (f"control F: {n['f_genes']} genes' components", n["f_share"], None, MID_GREY),
    ]
    for k, (_, value, interval, colour) in enumerate(rows):
        ax.barh(k, value, color=colour, height=0.62, zorder=2)
        right = value
        if interval is not None:
            ax.plot(interval, [k, k], color="0.1", lw=1.2, zorder=3)
            right = interval[1]
        ax.text(max(right, 0) + 0.01, k, f"{value:.0%}", va="center", fontsize=7.5)
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlim(min(0.0, min(r[1] for r in rows) - 0.05), 1.0)
    ticks = np.arange(0, 1.01, 0.2)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map, held out")
    tidy(ax)


def parts_panel(ax: plt.Axes, n: dict) -> None:
    """G of figure 03s1: the four parts with their 95% intervals over structures.

    A part below zero is drawn below zero, as computed; the weights are written
    under the panel's title.
    """
    parts, ci = n["parts"], n["parts_ci"]
    keys = list(PART_LABELS)
    for k, key in enumerate(keys):
        ax.barh(k, parts[key], color=PART_COLOURS[key], height=0.62, zorder=2)
        lo, hi = ci[key]
        ax.plot([lo, hi], [k, k], color="0.1", lw=1.2, zorder=3)
        ax.text(
            max(hi, 0) + 0.01,
            k,
            f"{parts[key]:+.0%} ({lo:+.0%} to {hi:+.0%})",
            va="center",
            fontsize=7.5,
        )
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([PART_LABELS[k] for k in keys], fontsize=8.5)
    ax.set_ylim(len(keys) - 0.4, -0.6)
    low = min(0.0, min(lo for lo, _ in ci.values()) - 0.05)
    ax.set_xlim(low, 1.0)

    # ticks every 20%, from the first at or above the axis' left end; adding 0.0
    # turns a rounded -0.0 into 0.0
    ticks = [round(t, 1) + 0.0 for t in np.arange(-1.0, 1.01, 0.2) if t >= low - 1e-9]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map")
    tidy(ax)


def calibration_panel(ax: plt.Axes, calibration: pd.DataFrame, n: dict) -> None:
    """H of figure 03s1: share left against the leftover's replication, floor and nano.

    Random folds only: every draw of the known map, from each half of the
    experiments, and nano with each half's predictors and with both.
    """
    calibration = calibration[calibration["folds"] == "random"]
    floor = calibration[calibration["map"] == FLOOR_MAP]
    for half, face in (("A", DARK_GREY), ("B", "white")):
        mine = floor[floor["truth_from"] == half]
        ax.scatter(
            mine["left"],
            mine["replication"],
            s=20,
            marker="s",
            facecolors=face,
            edgecolors=DARK_GREY,
            alpha=0.7,
            linewidths=0.8,
            zorder=2,
            label=f"the floor: Gria1 + synapse density only, from half {half} "
            f"({len(mine)} draws, median {mine['left'].median():.0%})",
        )
    nano = calibration[calibration["map"] == "nano"]
    ax.scatter(
        nano["left"],
        nano["replication"],
        s=46,
        facecolor="white",
        edgecolor=RED,
        linewidths=1.4,
        zorder=4,
        label=f"nano on the same {n['cal_n']} structures (predictors: half A, half B, "
        "both)",
    )
    lo, hi = n["left_ci"]
    ax.errorbar(
        [n["left"]],
        [n["rep_left"]],
        xerr=[[n["left"] - lo], [hi - n["left"]]],
        fmt="o",
        ms=9,
        color=RED,
        ecolor=RED,
        elinewidth=1.4,
        capsize=3,
        zorder=5,
        label=f"nano, the main model ({n['n_structures']} structures; 95% over "
        "structures)",
    )
    ax.set_xlim(0, max(0.5, calibration["left"].max() + 0.04, hi + 0.04))
    ax.set_ylim(min(0.6, calibration["replication"].min() - 0.03), 1.0)
    ticks = np.arange(0, ax.get_xlim()[1] + 1e-9, 0.1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map left")
    ax.set_ylabel("the leftover's replication (half against half)")
    ax.legend(loc="lower right", fontsize=6.8, frameon=False)
    tidy(ax)


def rows_table(n: dict) -> pd.DataFrame:
    """The rows of figure 03s1 A and B: the main model, its other folds, the check rows.

    One row each with kind, key, n_structures, left (with its 95%: over structures for
    the main model, over shufflings for the single shufflings) and, where a floor
    exists, nano minus it with its interval.
    """
    folds = n["folds"].copy()
    folds["nano_minus_floor"] = np.nan
    main = folds["kind"] == "main"
    folds.loc[main, "lo"] = n["left_ci"][0]
    folds.loc[main, "hi"] = n["left_ci"][1]
    folds.loc[main, "nano_minus_floor"] = n["minus_floor"]
    folds.loc[main, "nano_minus_floor_lo"] = n["minus_floor_ci"][0]
    folds.loc[main, "nano_minus_floor_hi"] = n["minus_floor_ci"][1]
    blocks = folds["key"] == "spatial_blocks"
    folds.loc[blocks, "nano_minus_floor"] = n["blocks"]["nano"] - n["blocks"]["floor"]
    checks = n["checks"].copy()
    checks["kind"] = "check"
    checks["lo"] = checks["hi"] = np.nan
    columns = [
        "kind",
        "key",
        "n_structures",
        "left",
        "lo",
        "hi",
        "nano_minus_floor",
        "nano_minus_floor_lo",
        "nano_minus_floor_hi",
    ]
    out = pd.concat([folds.reindex(columns=columns), checks.reindex(columns=columns)])
    return out.reset_index(drop=True)


def margin_label(r: pd.Series) -> str:
    """A row's nano minus its floor in points, with its interval when it has one."""
    text = f"{100 * r['nano_minus_floor']:+.0f}"
    if np.isfinite(r["nano_minus_floor_lo"]):
        lo, hi = 100 * r["nano_minus_floor_lo"], 100 * r["nano_minus_floor_hi"]
        text += f" ({lo:+.0f} to {hi:+.0f})"
    return text


def check_rows_panels(axes: tuple[plt.Axes, plt.Axes], n: dict) -> None:
    """A and B of figure 03s1: per row, the share left, and nano minus its floor.

    The main model first in red, its other folds, then the check rows; a dotted line
    between the kinds. A row on other structures than the main model's names their
    number; each value is written beside its dot.
    """
    table = rows_table(n)
    y = np.arange(len(table))
    left_ax, margin_ax = axes
    for ax in axes:
        for i in range(1, len(table)):
            if table.loc[i, "kind"] != table.loc[i - 1, "kind"]:
                ax.axhline(i - 0.5, color=MID_GREY, lw=0.6, ls=(0, (1, 2)), zorder=0)
    for i, r in table.iterrows():
        colour = RED if r["kind"] == "main" else "0.15"
        if np.isfinite(r["lo"]):
            left_ax.plot([r["lo"], r["hi"]], [i, i], color=MID_GREY, lw=2, zorder=1)
        left_ax.scatter(r["left"], i, s=34, color=colour, zorder=2)
        left_ax.text(
            max(r["left"], r["hi"] if np.isfinite(r["hi"]) else r["left"]) + 0.02,
            i,
            f"{r['left']:.0%}",
            fontsize=7.5,
            va="center",
            color=colour,
        )
        if not np.isfinite(r["nano_minus_floor"]):
            continue
        right = r["nano_minus_floor"]
        if np.isfinite(r["nano_minus_floor_lo"]):
            right = r["nano_minus_floor_hi"]
            margin_ax.plot(
                [r["nano_minus_floor_lo"], r["nano_minus_floor_hi"]],
                [i, i],
                color=MID_GREY,
                lw=2,
                zorder=1,
            )
        margin_ax.scatter(r["nano_minus_floor"], i, s=34, color=colour, zorder=2)
        margin_ax.text(
            right + 0.015, i, margin_label(r), fontsize=7.5, va="center", color=colour
        )
    labels = []
    for _, r in table.iterrows():
        label = ROW_LABELS.get(r["key"], r["key"])
        if r["n_structures"] != table["n_structures"].iloc[0]:
            label += f", {int(r['n_structures'])}"
        labels.append(label)
    left_ax.set_yticks(y)
    left_ax.set_yticklabels(labels, fontsize=8)
    margin_ax.set_yticks(y)
    margin_ax.set_yticklabels([])
    for ax in axes:
        ax.set_ylim(len(table) - 0.4, -0.8)
        tidy(ax)
    left_ax.set_xlim(0, 0.8)
    left_ax.set_xlabel("share of the reproducible map left")
    margin_ax.axvline(0, color="0.3", lw=0.7)
    lows = table["nano_minus_floor_lo"].dropna()
    highs = table["nano_minus_floor_hi"].dropna()
    margin_ax.set_xlim(min(-0.1, lows.min() - 0.05), max(0.3, highs.max() + 0.2))
    margin_ax.set_xlabel("nano minus its own floor, points (95% over structures)")
    for ax in axes:
        ticks = [t for t in ax.get_xticks() if ax.get_xlim()[0] <= t <= ax.get_xlim()[1]]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t:.0%}" for t in ticks])


def check_rows_line(n: dict) -> str:
    """The line under figure 03s1's title: the check rows against their floors."""
    checks = n["checks"].set_index("key")
    low = checks["nano_minus_floor"].idxmin()
    high = checks["nano_minus_floor"].idxmax()
    return (
        f"{checks_above(n)} stay above their own floor over the whole 95% interval; "
        f"nano minus it runs from {100 * checks.loc[low, 'nano_minus_floor']:+.0f} "
        f"points ({ROW_LABELS[low]}) to "
        f"{100 * checks.loc[high, 'nano_minus_floor']:+.0f} ({ROW_LABELS[high]}); the "
        f"main model leaves {n['left']:.0%}, {margin_text(n)} above its floor"
    )


def axes_at(
    fig: plt.Figure, left: float, bottom_in: float, width: float, height_in: float
) -> plt.Axes:
    """An axes placed by its width as a share of the figure and its height in inches."""
    height = fig.get_figheight()
    return fig.add_axes([left, bottom_in / height, width, height_in / height])


def plot_beyond_budget(
    y: np.ndarray,
    terms: dict[str, np.ndarray],
    held_out: np.ndarray,
    residual: np.ndarray,
    groups: list[str],
    acronyms: list[str],
    numbers: dict,
    calibration: pd.DataFrame,
    replication: pd.DataFrame,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 03s1: part 1's check rows, each against its own floor, and the detail.

    The arguments are those of plot_beyond, `terms`, the ranks of Gria1 and of the
    density term, and `replication`, the half-cohort table of beyond_density. Panels
    are placed in inches from the foot, so their titles and labels keep their room.
    """
    n = numbers
    fig = plt.figure(figsize=(16, 23))
    heading(fig, "beyond_budget", check_rows_line(n))

    # A and B: every check row and fold, what it leaves and nano minus its floor
    axes = (axes_at(fig, 0.22, 16.6, 0.3, 4.55), axes_at(fig, 0.58, 16.6, 0.3, 4.55))
    check_rows_panels(axes, n)
    panel_title(
        axes[0],
        "A",
        "The check rows and other folds: what each leaves",
        "each on its own structures, its ceiling recomputed there;\nthe main model "
        "with its 95% over structures",
    )
    panel_title(
        axes[1],
        "B",
        "Nano minus each one's own floor",
        "a floor misses what is measured once (PSD95, three of the 11 markers);\n"
        f"spatial blocks: nano {n['blocks']['nano']:.0%}, floor "
        f"{n['blocks']['floor']:.0%}",
    )

    # C to E: the map against each term and against the whole model
    fig.legend(
        handles=group_handles(),
        loc="center",
        bbox_to_anchor=(0.5, 15.7 / fig.get_figheight()),
        ncol=4,
        fontsize=8,
        frameon=False,
        title="dots of C to E: structures, by group of divisions",
        title_fontsize=8,
    )
    bottom = 11.7 / fig.get_figheight()
    predictor_scatters(fig, y, terms, held_out, residual, groups, acronyms, n, bottom)

    # F and G: each model's share, and the four parts
    ax = axes_at(fig, 0.15, 7.6, 0.3, 2.8)
    shares_panel(ax, n)
    panel_title(ax, "F", "What each model predicts, held out")
    ax = axes_at(fig, 0.6, 7.6, 0.3, 2.8)
    parts_panel(ax, n)
    weights = ", ".join(
        f"{k} {v:+.2f} ({n['weights_ci'][k][0]:+.2f} to {n['weights_ci'][k][1]:+.2f})"
        for k, v in n["weights"].items()
    )
    panel_title(
        ax, "G", "The four parts, 95% over structures", f"weights (z-scored): {weights}"
    )

    # H: the calibration; I: the replication
    ax = axes_at(fig, 0.07, 2.2, 0.38, 3.9)
    calibration_panel(ax, calibration, n)
    floor = n["floor"]
    panel_title(
        ax,
        "H",
        "The same model on a map whose answer is known",
        f"on {n['cal_n']} structures: the floor {floor['left_median']:.0%} "
        f"({floor['left_lo']:.0%} to {floor['left_hi']:.0%}); nano "
        f"{n['nano_cal_left']:.0%}, minus the floor {margin_text(n)}",
    )
    ax = axes_at(fig, 0.6, 2.2, 0.36, 3.9)
    replication_panel(ax, replication, n)
    panel_title(
        ax,
        "I",
        "Does the leftover replicate across mice?",
        f"{len(replication)} ways to split the {n['n_adults']} adults into two fives; "
        f"two unrelated leftovers\nwould agree between {n['noise'][0]:+.2f} and "
        f"{n['noise'][1]:+.2f}, far left of this axis",
    )
    footer(
        fig,
        [
            "How to read: the map and every predictor are ranks across structures; each "
            "term enters straight, and each share is the held-out R² (20 shufflings of "
            "five folds) over the ceiling; intervals are 95% over structures "
            "(jackknife).",
            "Gria1 only = both - density alone; density only = both - Gria1 alone; "
            "shared = Gria1 alone + density alone - both; left = 1 - both. Held out, a "
            "part can come out negative, and is drawn as computed.",
            "The calibration builds a map made only of Gria1 and synapse density from "
            "half of each gene's Allen experiments, gives it ten made-up adults as noisy "
            "as ours, and predicts it from the other half; every check row has a floor",
            "of its own, built the same way with its own model. A floor holds the "
            "mismatch of one Allen map with another, not that of Allen's P56 mice with "
            "these brains, so it errs low.",
        ],
    )
    return saved(fig, save)


# the seven controls of figure 03s2 by their letter in controls.csv: the question
# each asks, as a panel title
CONTROL_QUESTIONS = {
    "A": "A smooth gradient across the brain?",
    "B": "Small structures, where a mean is noisy?",
    "C": "One or two animals?",
    "D": "The whisker manipulation?",
    "E": "Curvature the straight model misses?",
    "F": "The choice of predictors?",
    "G": "The reading of the map?",
}


def control_title(ax: plt.Axes, letter: str, controls: pd.DataFrame) -> None:
    """A control's panel title: its question, its verdict and its number."""
    row = controls[controls["control"].str[0] == letter].iloc[0]
    verdict = "passes" if row["verdict"] == "pass" else "look closer"
    panel_title(
        ax,
        letter,
        CONTROL_QUESTIONS[letter],
        f"{verdict}: {textwrap.fill(row['number'], 52)}",
    )


def leftover_scatter(ax: plt.Axes, x: np.ndarray, y: np.ndarray, xlabel: str) -> None:
    """The leftover of each structure against one of its properties."""
    ok = np.isfinite(x)
    ax.scatter(x[ok], y[ok], s=35, color=DARK_GREY, linewidths=0, alpha=0.85, zorder=2)
    ax.axhline(0, color=MID_GREY, lw=0.7, zorder=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("leftover (ranks above prediction)")
    tidy(ax)


def curvature_panel(ax: plt.Axes, scores: tuple[float, float, float]) -> None:
    """E: the held-out R2 of the model straight, curved to cubes and to fifth powers."""
    labels = ["straight\n(the model)", "x, x², x³\n(row curved)", "to x⁵"]
    colours = [DARK_GREY, MID_GREY, MID_GREY]
    ax.bar(range(3), scores, color=colours, width=0.6)
    for k, value in enumerate(scores):
        ax.text(k, value + 0.01, f"{value:.3f}", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(range(3))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("R² on structures the fit has not seen")
    tidy(ax)


def gene_space_panel(ax: plt.Axes, calibration: pd.DataFrame) -> None:
    """F: what control F's model leaves of nano and of maps made of the genes."""
    rows = (
        (
            "gene space",
            "a map made of the\ngenes' expression\n(its own floor)",
            DARK_GREY,
        ),
        ("nano", "the nano map", RED),
    )
    rng = np.random.default_rng(0)
    for y, (kind, _, colour) in enumerate(rows):
        left = calibration.loc[calibration["map"] == kind, "left"].to_numpy(float)
        jitter = rng.uniform(-0.12, 0.12, len(left))
        ax.scatter(left, y + jitter, s=24, color=colour, alpha=0.8, linewidths=0)
        ax.text(
            left.max() + 0.012,
            y,
            f"{left.min():.0%} to {left.max():.0%}",
            va="center",
            fontsize=7.5,
            color=colour,
        )
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([label for _, label, _ in rows], fontsize=8)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlim(0, max(0.4, float(calibration["left"].max()) + 0.1))
    ticks = np.arange(0, ax.get_xlim()[1] + 1e-9, 0.1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map left")
    tidy(ax)


def readings_panel(ax: plt.Axes, readings: pd.DataFrame, threshold: float) -> None:
    """G: the leftover's replication under every reading of the map."""
    y = np.arange(len(readings))
    ax.scatter(readings["leftover_replication"], y, s=35, color=DARK_GREY, zorder=2)
    ax.axvline(threshold, color=RED, lw=0.9, ls=(0, (3, 2)))
    ax.text(
        threshold,
        len(readings) - 0.5,
        f" the control's line, {threshold}",
        fontsize=7,
        color=RED,
        va="bottom",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(readings["reading"], fontsize=8)
    ax.set_ylim(len(readings) - 0.3, -0.7)
    ax.set_xlim(min(0.7, threshold - 0.05), 1.0)
    ax.set_xlabel("the leftover's replication, half against half")
    tidy(ax)


def artefact_panels(fig: plt.Figure, grid, data: dict) -> None:
    """A to D of figure 03s2: position, size, pairs of adults, naive against RWS.

    The leftover against position and size, the pairs of adults, the naive group against
    the RWS group.
    """
    controls = data["controls"]
    res = data["residual"]

    # the leftover against each structure's position and size
    ax = fig.add_subplot(grid[0, 0:3])
    leftover_scatter(ax, data["ap"], res, "structure centroid, anterior-posterior (mm)")
    control_title(ax, "A", controls)

    ax = fig.add_subplot(grid[0, 3:6])
    leftover_scatter(ax, np.log10(data["sizes"]), res, "log10 volume (20 um voxels)")
    control_title(ax, "B", controls)

    # every pair of adults, and the naive group against the RWS group
    ax = fig.add_subplot(grid[0, 6:9])
    ax.hist(data["pairs"], bins=15, color=MID_GREY)
    ax.axvline(float(np.median(data["pairs"])), color=RED, lw=1.5)
    ax.axvline(data["thresholds"]["pair_rho"], color=RED, lw=0.9, ls=(0, (3, 2)))
    ax.set_xlabel("one adult's leftover against another's (rho)")
    ax.set_ylabel("pairs of adults")
    tidy(ax)
    control_title(ax, "C", controls)

    ax = fig.add_subplot(grid[0, 9:12])
    naive, rws = data["groups"]["naive"], data["groups"]["rws"]
    ax.scatter(naive, rws, s=35, color=DARK_GREY, linewidths=0, alpha=0.85)
    low = min(naive.min(), rws.min())
    high = max(naive.max(), rws.max())
    ax.plot([low, high], [low, high], color=MID_GREY, lw=0.8, ls=(0, (4, 3)))
    ax.set_xlabel("leftover of the five naive adults")
    ax.set_ylabel("leftover of the five RWS adults")
    tidy(ax)
    control_title(ax, "D", controls)


def plot_beyond_controls(
    data: dict, numbers: dict, save: Path | None = None
) -> plt.Figure:
    """Figure 03s2: the seven attempts to break part 1's leftover, one panel each.

    `data` is adult.beyond_figures.control_data: the leftover per structure with
    its anterior-posterior centroid (mm) and mean volume (20 um voxels), the
    agreement of every pair of adults' leftovers, the naive and the RWS group's
    leftovers, the held-out R2 straight, cubic and to fifth powers, control F's
    own calibration, control G's readings, the verdicts and the thresholds.
    """
    controls = data["controls"]
    passed = controls[controls["verdict"] == "pass"]
    failed = [c[0] for c in controls.loc[controls["verdict"] != "pass", "control"]]
    fig = plt.figure(figsize=(16, 10))
    if failed:
        line = (
            f"{len(passed)} of {len(controls)} controls pass; look closer at "
            + ", ".join(failed)
        )
    else:
        line = (
            f"All {len(controls)} controls pass: the leftover is not a gradient, small "
            "structures, one animal, the whisker manipulation, missing curvature, the "
            "choice of predictors or the reading"
        )
    heading(fig, "beyond_controls", line)
    grid = fig.add_gridspec(
        2, 12, hspace=0.55, wspace=1.6, left=0.06, right=0.98, top=0.83, bottom=0.12
    )

    # A to D: four ways the leftover could be an artefact of the data
    artefact_panels(fig, grid, data)

    # E to G: whether the model, the predictors or the reading made it
    ax = fig.add_subplot(grid[1, 0:4])
    curvature_panel(ax, data["curvature"])
    control_title(ax, "E", controls)

    ax = fig.add_subplot(grid[1, 4:8])
    gene_space_panel(ax, data["f_calibration"])
    control_title(ax, "F", controls)

    ax = fig.add_subplot(grid[1, 8:12])
    readings_panel(ax, data["readings"], data["thresholds"]["readings_replication"])
    control_title(ax, "G", controls)
    footer(
        fig,
        [
            "How to read: each panel tries one way the leftover of the main model "
            f"({figure_ref('beyond')}) could be an artefact; the line each must not "
            "cross is fixed in settings.toml [beyond_controls] (dashed in C and G). In C "
            "the median pair is red.",
            "E and F score on structures the fit has not seen; F gives the model the "
            "components of every gene measured in all the structures, picked inside "
            "each training fold, and sets nano beside its own floor.",
        ],
    )
    return saved(fig, save)


# ===== 04 and 11s1 Where the leftover lives, and every gene against it =====

# the structures shown in the bars, each way
N_LEFTOVER_BARS = 12

# the leftover's colour scale on the planes, in ranks each way: most structures sit
# within 20 ranks of prediction, so a wider scale would wash the planes out
LEFTOVER_SPAN = 30

# the genes drawn against the leftover: the best, and these named in any case
N_LEFTOVER_GENES = 15
LEFTOVER_NAMED = ("Cacng8", "Gria1", "Dlg2")


def leftover_planes(fig: plt.Figure, grid, planes: list[dict], span: float) -> list:
    """A: map, prediction and leftover on each plane; returns the two images to key."""
    images = []
    titles = ("the nano map", "what Gria1 and synapse density predict", "the leftover")
    for r, plane in enumerate(planes):
        lab = plane["lab"]
        no_value = np.isnan(plane["map"]) & (lab > 0)
        for c, key in enumerate(("map", "prediction", "leftover")):
            ax = fig.add_subplot(grid[r, c])
            if key == "leftover":
                image = draw_plane(
                    ax, plane[key], lab, plt.get_cmap("RdBu_r"), -span, span
                )
            else:
                values = np.clip(plane[key], 0, 1)
                image = draw_plane(ax, values, lab, hot_cut(), RANK_FLOOR, 1.0)
            overlay(ax, no_value, NO_DATA_GREY, lab)
            if r == 0:
                ax.set_title(titles[c], fontsize=9)
            if c == 0:
                ax.set_ylabel(f"CCF plane {plane['ccf_plane']}", fontsize=8)
            if r == len(planes) - 1 and c in (0, 2):
                images.append(image)
    return images


def leftover_bars(ax: plt.Axes, residuals: pd.DataFrame, t_max: float) -> None:
    """B: the structures with the largest leftovers each way, grey by t across adults."""
    show = pd.concat(
        [residuals.head(N_LEFTOVER_BARS), residuals.tail(N_LEFTOVER_BARS)]
    ).reset_index(drop=True)
    y = np.arange(len(show))
    ax.barh(
        y,
        show["residual"],
        color=bars_grey(show["t_adults"].abs().to_numpy(), t_max),
        height=0.7,
        zorder=2,
    )
    for i, r in show.iterrows():
        ha = "left" if r["residual"] > 0 else "right"
        ax.text(
            r["residual"] + (1.0 if r["residual"] > 0 else -1.0),
            i,
            f"{r['residual']:+.0f}",
            fontsize=6.5,
            va="center",
            ha=ha,
            color=DARK_GREY,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r['acronym']}  {r['structure']}" for _, r in show.iterrows()],
        fontsize=6.8,
    )
    ax.invert_yaxis()
    ax.axvline(0, color="0.3", lw=0.7)
    span = float(show["residual"].abs().max()) * 1.25
    ax.set_xlim(-span, span)
    ax.set_xlabel("nano rank minus predicted rank")
    grey_key(ax, t_max, "|t| across the adults")
    tidy(ax)


def grey_key(ax: plt.Axes, t_max: float, label: str) -> None:
    """A small key to the grey of a panel's bars, inside it at the lower right.

    From t = 0 to black at `t_max`.
    """
    key = ax.inset_axes([0.66, 0.03, 0.3, 0.025])
    steps = np.linspace(0, t_max, 6)
    colours = bars_grey(steps, t_max)[None, :, :]
    key.imshow(colours, aspect="auto", extent=(0, t_max, 0, 1))
    key.set_yticks([])
    key.set_xticks([0, t_max])
    key.set_xticklabels(["0", f"{t_max:g}"], fontsize=6.5)
    key.tick_params(length=2, pad=1)
    key.set_title(label, fontsize=6.5, pad=2)


def leftover_gene_bars(
    ax: plt.Axes, genes: pd.DataFrame, subunits: set[str], t_max: float, q: float
) -> None:
    """A: the genes closest to the leftover and the named ones, with their null bands.

    Cacng8, named for the leftover in advance, in red; the members of the AMPA
    receptor complex family tagged.
    """
    top = genes.head(N_LEFTOVER_GENES)
    named = genes[
        genes["symbol"].isin(LEFTOVER_NAMED) & ~genes["symbol"].isin(top["symbol"])
    ]
    show = pd.concat([top, named]).reset_index(drop=True)
    y = np.arange(len(show))
    for i, r in show.iterrows():
        ax.add_patch(
            Rectangle(
                (r["null_lo"], i - 0.42),
                r["null_hi"] - r["null_lo"],
                0.84,
                color=NULL_BAND,
                lw=0,
                zorder=0,
            )
        )
    ax.barh(
        y,
        show["rho"],
        color=bars_grey(show["t_boot"].to_numpy(), t_max),
        height=0.62,
        zorder=2,
    )
    ax.axvline(0, color="0.3", lw=0.6, zorder=1)
    if len(named):
        ax.axhline(len(top) - 0.5, color="0.6", lw=0.6, ls=(0, (3, 2)))
    ax.set_yticks(y)
    labels = []
    for _, r in show.iterrows():
        text = f"{r['symbol']}  (rank {int(r['rank_all'])}; p {r['p_spatial']:.4f})"
        if r["symbol"] in AMPA_FAMILY:
            text += "  [AMPA complex]"
        if r["in_model"]:
            text += f"  in the model: {r['in_model']}"
        labels.append(text)
    ax.set_yticklabels(labels, fontsize=7)
    for tick, (_, r) in zip(ax.get_yticklabels(), show.iterrows()):
        tick.set_color(gene_colour(r["symbol"], subunits))
        if r["symbol"] == LEFTOVER_GENE:
            tick.set_color(RED)
        if r["q_all"] < q:
            tick.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlim(-0.6, 0.6)
    ax.set_xlabel("Spearman rho of the gene with the leftover")
    tidy(ax)


def leftover_sets(ax: plt.Axes, genes: pd.DataFrame, sets: pd.DataFrame) -> None:
    """B: the gene sets of analysis 3 against the leftover, each with its null band."""
    rng = np.random.default_rng(0)
    by_set = sets.set_index("gene_set")
    labels = []
    for k, name in enumerate(by_set.index):
        member = [
            name in [t.strip() for t in str(g).split(";")] for g in genes["gene_sets"]
        ]
        rows = genes[member].copy()
        rows["rho_nano"] = rows["rho"]
        test = by_set.loc[name]
        set_column(ax, k, rows, test, SET_COLOURS[name], rng)
        if test["tested"]:
            stats = f"p = {test['p_spatial']:.3f}\nq = {test['q']:.2f}"
        else:
            stats = "not tested:\ntoo few genes"
        labels.append(f"{name}\nn = {int(test['n_genes'])}\n{stats}")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.axhline(0, color="0.6", lw=0.6, zorder=1)
    ax.set_xlim(-0.6, len(labels) - 0.3)
    ax.set_ylabel("Spearman rho with the leftover")
    tidy(ax)


def plot_beyond_where(
    planes: list[dict],
    residuals: pd.DataFrame,
    numbers: dict,
    t_max: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 04: where the leftover lives.

    `planes` holds per plane its labels and the map, prediction and leftover
    painted (adult.beyond_figures.plane_images); `residuals` is
    residual_by_structure.csv; `t_max` the t at which a structure's bar turns black.
    """
    fig = plt.figure(figsize=(16, 12))
    top_up = ", ".join(residuals["acronym"].head(4))
    top_down = ", ".join(residuals["acronym"].tail(4)[::-1])
    clipped = int((residuals["residual"].abs() > LEFTOVER_SPAN).sum())
    heading(
        fig,
        "beyond_where",
        f"{numbers['n_structures']} structures of the fit; most above prediction: "
        f"{top_up}; most below: {top_down}",
    )

    # A: the three planes
    grid = fig.add_gridspec(
        len(planes),
        3,
        left=0.04,
        right=0.58,
        top=0.86,
        bottom=0.16,
        hspace=0.12,
        wspace=0.04,
    )
    images = leftover_planes(fig, grid, planes, LEFTOVER_SPAN)
    cax = fig.add_axes([0.06, 0.12, 0.3, 0.012])
    cb = fig.colorbar(images[0], cax=cax, orientation="horizontal")
    cb.ax.set_xlim(0, 1)
    cb.set_label("rank among the structures of the fit (0 low, 1 high)", fontsize=8)
    cax = fig.add_axes([0.42, 0.12, 0.15, 0.012])
    cb = fig.colorbar(images[1], cax=cax, orientation="horizontal", extend="both")
    cb.set_label(
        f"leftover (ranks): red above prediction;\nclipped at ±{LEFTOVER_SPAN} "
        f"({clipped} structures beyond)",
        fontsize=8,
    )
    fig.text(
        0.04,
        0.885,
        f"A.  The map, the prediction and the leftover on {len(planes)} coronal "
        "planes\nflat grey: structures the fit does not use",
        fontsize=9,
        va="bottom",
    )

    # B: the structures with the largest leftovers
    ax = fig.add_axes([0.8, 0.16, 0.17, 0.7])
    leftover_bars(ax, residuals, t_max)
    panel_title(
        ax,
        "B",
        f"The {N_LEFTOVER_BARS} largest leftovers each way",
        f"grey: how steady across the ten adults'\nown leftovers (t; black at {t_max:g})",
    )
    footer(
        fig,
        [
            "How to read: the leftover is the nano rank minus the rank Gria1 and synapse "
            "density predict (the straight model fitted on every structure), on the "
            "structures of the fit; the bars' grey says how steady it is across the "
            "ten adults' own leftovers.",
            "What would mean what: a structure far from prediction in every adult is "
            "where the departure sits; a claim that one structure stands out needs a "
            "null of its own, which is not drawn here.",
        ],
    )
    return saved(fig, save)


def plot_leftover_genes(
    genes: pd.DataFrame,
    sets: pd.DataFrame,
    numbers: dict,
    n_surrogates: int,
    t_max: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 11s1: does any gene's map, or any gene set, follow the leftover?

    `genes` and `sets` are leftover_genes.csv and leftover_sets.csv; `t_max` the
    rho over its SD across resampled adults at which a gene's bar turns black.
    """
    n = numbers
    q = n["q"]
    fig = plt.figure(figsize=(16, 11))
    tested = sets[sets["tested"]]
    passed = tested[tested["q"] < q]
    heading(
        fig,
        "leftover_genes",
        f"{n['n_genes']} genes against the leftover of the main model "
        f"({figure_ref('beyond_budget')}): {n['n_pass']} past its null at BH q < {q} "
        f"({n['n_p05']} below p 0.05 before correction); gene sets past it: "
        + (", ".join(passed["gene_set"]) if len(passed) else "none"),
    )
    subunits = {
        s for s, t in zip(genes["symbol"], genes["gene_sets"]) if "subunits" in str(t)
    }
    ax = fig.add_axes([0.22, 0.16, 0.24, 0.68])
    leftover_gene_bars(ax, genes, subunits, t_max, q)
    cacng8 = genes.set_index("symbol").loc["Cacng8"]
    panel_title(
        ax,
        "A",
        "The genes closest to the leftover, and the named ones",
        f"pale blue: 95% of rho with the leftover's {n_surrogates} surrogates, each "
        "through\nthe same fit; grey: rho over its SD across resampled adults (black "
        f"at {t_max:g});\nbold: past BH; Cacng8 {cacng8['rho']:+.2f} "
        f"({p_text(cacng8['p_spatial'], n_surrogates)}, rank {int(cacng8['rank_all'])})",
    )
    ax = fig.add_axes([0.56, 0.2, 0.41, 0.62])
    leftover_sets(ax, genes, sets)
    panel_title(
        ax,
        "B",
        "The gene sets of analysis 3 against the leftover",
        "pale blue: 95% of the set's median rho with the surrogates; spatial p, and "
        "q by BH over the tested sets",
    )
    footer(
        fig,
        [
            "How to read: the leftover is what the main model (Gria1 + synapse density, "
            "straight) leaves, so it carries nothing of the model's columns; each "
            "surrogate goes through the same fit before it is correlated (a",
            "Freedman-Lane null), so a gene sharing the model's pattern meets a null of "
            "the right width. Gria1 enters the model and sits near zero by construction; "
            "the density genes enter only through their mean rank.",
            "Named in advance: Cacng8 (red), its uncorrected p the test, and the AMPA "
            f"receptor complex family (tagged; {figure_ref('ampa_family')}). Every other "
            "gene and every set here is exploratory, BH over all of them.",
            "The glia set followed the leftover of the first, four-subunit model of 8 "
            "October, a test not named in advance: an unplanned lead, not pursued.",
        ],
    )
    return saved(fig, save)


# ===== 15 and 15s What the green channel reports =====

# the three channels as figure 15 names and draws them: name, what it records, the
# box colour and the per-mouse dot colour
CHANNEL_BOXES = {
    "nano": ("nanobody against the tag;\nsections not permeabilised", NANO, NANO_DOT),
    "SEP": (
        "the tag's own green fluorescence;\nexpected to show the tagged\nreceptor "
        "wherever it sits",
        SEP,
        SEP_DOT,
    ),
    "auto": ("no label: the tissue's\nautofluorescence", AUTO, AUTO_DOT),
}


def paired_columns(
    ax: plt.Axes,
    rows: pd.DataFrame,
    columns: list[tuple[str, str, str]],
    x0: float = 0.0,
    fmt: str = "+.2f",
) -> list[float]:
    """A column of dots per (label, table column, colour), one dot per adult.

    Lines join the same adult across the columns; a bar and its value mark the
    mean, written with `fmt`. Returns the x of each column.
    """
    xs = [x0 + i for i in range(len(columns))]
    values = [rows[key].to_numpy(float) for _, key, _ in columns]
    for m in range(len(rows)):
        ax.plot(xs, [v[m] for v in values], color=PAIR_LINE, lw=0.7, zorder=1)
    for x, v, (_, _, colour) in zip(xs, values, columns):
        ax.scatter([x] * len(v), v, s=26, color=colour, linewidths=0, zorder=3)
        ax.plot([x - 0.22, x + 0.22], [v.mean()] * 2, color="0.1", lw=2, zorder=4)
        ax.text(
            x + 0.27,
            v.mean(),
            format(v.mean(), fmt),
            fontsize=7.5,
            va="center",
            zorder=5,
            bbox=dict(fc="white", ec="none", pad=0.8),
        )
    return xs


def channel_diagram(ax: plt.Axes, rows: pd.DataFrame) -> None:
    """A: the three channels, and how each pair agrees across structures."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    where = {"nano": (0.17, 0.78), "SEP": (0.83, 0.78), "auto": (0.5, 0.2)}
    edges = (
        ("nano", "SEP", "rho_sep_nano", (0.5, 0.86)),
        ("SEP", "auto", "rho_sep_auto", (0.8, 0.45)),
        ("nano", "auto", "rho_nano_auto", (0.2, 0.45)),
    )
    for a, b, key, (tx, ty) in edges:
        v = rows[key].to_numpy(float)
        strong = key == "rho_sep_auto"
        ax.plot(
            [where[a][0], where[b][0]],
            [where[a][1], where[b][1]],
            color=DARK_GREY if strong else LIGHT_GREY,
            lw=3.0 if strong else 1.4,
            zorder=1,
        )
        ax.text(
            tx,
            ty,
            f"rho {v.mean():+.2f}\n({v.min():+.2f} to {v.max():+.2f})",
            ha="center",
            va="center",
            fontsize=8,
            zorder=5,
            fontweight="bold" if strong else "normal",
            bbox=dict(fc="white", ec="none", pad=1.5),
        )
    for name, (x, y) in where.items():
        body, colour, _ = CHANNEL_BOXES[name]
        ax.add_patch(
            FancyBboxPatch(
                (x - 0.15, y - 0.11),
                0.3,
                0.22,
                boxstyle="round,pad=0.01,rounding_size=0.02",
                fc="white",
                ec=colour,
                lw=2.2,
                zorder=3,
            )
        )
        ax.text(
            x,
            y + 0.06,
            name,
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
            zorder=4,
        )
        ax.text(
            x,
            y - 0.035,
            body,
            ha="center",
            va="center",
            fontsize=7.4,
            color="0.25",
            zorder=4,
            linespacing=1.2,
        )


def raw_planes(
    fig: plt.Figure,
    axes: list[plt.Axes],
    images: dict[str, np.ndarray],
    lab: np.ndarray,
    one: pd.Series,
) -> None:
    """B: one adult's three raw channels on the plane, each on its own scale."""
    titles = {
        "sig": f"nano (with autofluorescence: rho {one['rho_nano_auto']:+.2f} here)",
        "sep": f"SEP (with autofluorescence: rho {one['rho_sep_auto']:+.2f} here)",
        "auto": "autofluorescence",
    }
    for ax, key in zip(axes, ("sig", "sep", "auto")):
        img = images[key]
        inside = np.isfinite(img) & (lab > 0)
        lo, hi = np.percentile(img[inside], [1, 99])
        image = draw_plane(ax, img, lab, hot_cut(), lo, hi)
        colour_bar(fig, ax, image, "counts above background")
        ax.set_title(titles[key], fontsize=9)


def channel_columns(ax: plt.Axes, rows: pd.DataFrame) -> None:
    """Per adult, which channel follows which, joined by adult.

    SEP with autofluorescence, SEP with nano, nano with autofluorescence.
    """
    columns = [
        ("SEP with\nautofluorescence", "rho_sep_auto", SEP_DOT),
        ("SEP with\nnano", "rho_sep_nano", SEP_DOT),
        ("nano with\nautofluorescence", "rho_nano_auto", NANO_DOT),
    ]
    paired_columns(ax, rows, columns)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([label for label, _, _ in columns])
    ax.axhline(0, color="0.8", lw=0.7, zorder=0)
    ax.set_xlim(-0.6, 2.8)
    lowest = float(rows[[key for _, key, _ in columns]].min().min())
    ax.set_ylim(min(-0.05, lowest - 0.05), 1.0)
    ax.set_ylabel("Spearman rho across structures, per adult")
    tidy(ax)


def range_columns(ax: plt.Axes, rows: pd.DataFrame) -> None:
    """Per adult, how much each channel varies across structures (p90 - p10 of log2)."""
    columns = [
        ("nano", "range_nano", NANO_DOT),
        ("autofluorescence", "range_auto", AUTO_DOT),
        ("SEP", "range_sep", SEP_DOT),
    ]
    paired_columns(ax, rows, columns, fmt=".2f")
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([label for label, _, _ in columns])
    ax.set_xlim(-0.6, 2.8)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("p90 - p10 across structures (log2)")
    tidy(ax)


def plot_green_channel(
    rows: pd.DataFrame, n_structures: int, save: Path | None = None
) -> plt.Figure:
    """Figure 15: whether the green channel reports the tag or the tissue.

    `rows` is sep_channel_check.csv, one row per adult.
    """
    sep_auto = rows["rho_sep_auto"]
    if sep_auto.min() > 0.5:
        verdict = (
            "so these brains cannot give total receptor, and the surface fraction "
            "stays an interpretation; a total-GluA1 stain would measure it"
        )
    else:
        verdict = "so it carries more than the tissue"
    fig = plt.figure(figsize=(16, 6.2))
    heading(
        fig,
        "green_channel",
        f"In every adult the green (SEP) channel follows autofluorescence "
        f"({sep_auto.min():.2f} to {sep_auto.max():.2f}), {verdict}",
    )
    ax = fig.add_axes([0.02, 0.12, 0.38, 0.66])
    channel_diagram(ax, rows)
    panel_title(
        ax,
        "A",
        "Three channels of the same sections",
        "Spearman across structures, mean over the adults (range)",
    )
    ax = fig.add_axes([0.48, 0.16, 0.24, 0.6])
    channel_columns(ax, rows)
    panel_title(ax, "B", "Adult by adult", "one dot per adult, lines join one adult")
    ax = fig.add_axes([0.79, 0.16, 0.18, 0.6])
    range_columns(ax, rows)
    panel_title(
        ax,
        "C",
        "How much each channel varies",
        "a receptor channel should vary as nano does",
    )
    footer(
        fig,
        [
            f"How to read: one value per declared structure ({n_structures}) and adult, "
            "the mean of a raw channel in log2; no ratio of channels is taken. One "
            "adult's raw planes, and each channel against Gria1: "
            f"{figure_ref('green_channel_detail')}.",
        ],
    )
    return saved(fig, save)


def plot_green_channel_detail(
    rows: pd.DataFrame,
    images: dict[str, np.ndarray],
    lab: np.ndarray,
    mouse: str,
    plane: int,
    n_structures: int,
    gene: str,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 15s: whether the green channel reports the tag or the tissue, in detail.

    `rows` is sep_channel_check.csv (one row per adult), `images` the raw channels
    (sig, sep, auto) of `mouse` on CCF plane `plane` with its labels `lab`;
    `gene` is the gene the channels are compared with (Gria1).
    """
    fig = plt.figure(figsize=(16, 16))
    sep_auto = rows["rho_sep_auto"]
    nano_auto = rows["rho_nano_auto"]
    heading(
        fig,
        "green_channel_detail",
        f"{len(rows)} adults, {n_structures} declared structures: SEP follows "
        f"autofluorescence at {sep_auto.min():.2f} to {sep_auto.max():.2f} in every "
        f"adult; nano follows autofluorescence at only {nano_auto.min():.2f} to "
        f"{nano_auto.max():.2f}",
    )

    # A: the three channels; C: what each tracks, adult by adult
    ax = fig.add_axes([0.03, 0.6, 0.4, 0.3])
    channel_diagram(ax, rows)
    panel_title(
        ax,
        "A",
        "Three channels of the same sections, and how they agree",
        "Spearman across structures, mean over the adults (range); no ratio of "
        "channels is taken",
    )
    ax = fig.add_axes([0.55, 0.62, 0.42, 0.26])
    channel_columns(ax, rows)
    panel_title(
        ax,
        "C",
        "Whom each channel follows, adult by adult",
        "one dot per adult, lines join the same adult; bar: the mean",
    )

    # B: the raw channels of one adult
    axes = [fig.add_axes([0.03 + i * 0.33, 0.33, 0.27, 0.22]) for i in range(3)]
    one = rows.set_index("mouse").loc[mouse]
    raw_planes(fig, axes, images, lab, one)
    fig.text(
        0.03,
        0.565,
        f"B.  The three raw channels of one adult ({mouse.split('_')[0]}), CCF plane "
        f"{plane}, each on its own scale (1st to 99th percentile in the brain)",
        fontsize=9,
        va="bottom",
    )

    # D: the range of each channel; E: against Gria1, and what is left of SEP
    ax = fig.add_axes([0.06, 0.1, 0.3, 0.18])
    range_columns(ax, rows)
    panel_title(
        ax,
        "D",
        "How much each channel varies across the brain",
        "a channel reporting the receptor should vary about as much as nano",
    )
    ax = fig.add_axes([0.45, 0.1, 0.52, 0.18])
    gria = [
        ("nano", "rho_nano_gria", NANO_DOT),
        ("autofluo-\nrescence", "rho_auto_gria", AUTO_DOT),
        ("SEP", "rho_sep_gria", SEP_DOT),
        ("SEP minus its\nautofluo. part", "rho_sepresid_gria", SEP_REMAINDER),
    ]
    xs = paired_columns(ax, rows, gria)
    xs += paired_columns(
        ax,
        rows,
        [("SEP minus its\nautofluo. part", "rho_sepresid_nano", SEP_REMAINDER)],
        x0=len(gria) + 0.8,
    )
    ax.set_xticks(xs)
    ax.set_xticklabels([g[0] for g in gria] + ["SEP minus its\nautofluo. part"])
    ax.text(1.5, 1.0, f"with {gene} mRNA", ha="center", fontsize=8.5, color=DARK_GREY)
    ax.text(xs[-1], 1.0, "with nano", ha="center", fontsize=8.5, color=DARK_GREY)
    ax.axhline(0, color="0.8", lw=0.7, zorder=0)
    ax.set_xlim(-0.6, xs[-1] + 0.8)
    ax.set_ylim(min(-0.1, float(rows[[g[1] for g in gria]].min().min()) - 0.05), 1.05)
    ax.set_ylabel("Spearman rho across structures, per adult")
    panel_title(
        ax,
        "E",
        f"Against {gene} mRNA, and what is left of SEP once autofluorescence "
        "is taken out",
        "SEP minus its autofluorescence part: the residual of log2 SEP regressed on "
        "log2 autofluorescence, in each adult",
    )
    tidy(ax)
    footer(
        fig,
        [
            "How to read: one value per declared structure and adult, the mean of a "
            "raw channel (counts above its off-tissue background), in log2; no ratio "
            "of channels is taken anywhere.",
            "What would mean what: SEP following nano and Gria1 more than "
            "autofluorescence: the green channel reports the tag. SEP following "
            "autofluorescence in every adult: it reports mostly the tissue,",
            "and whether receptor at the membrane and receptor anywhere differ across "
            "the brain cannot be read from these channels; a total-GluA1 stain on the "
            "same brains would answer it.",
            "What is left of SEP once its autofluorescence part is taken out may be tag "
            "that survived, or nano's fluorescence leaking into the green channel; the "
            "filter sets decide which.",
        ],
    )
    return saved(fig, save)


def channel_dots(
    ax: plt.Axes,
    rows: pd.DataFrame,
    columns: tuple[str, ...],
    red: tuple[str, ...],
    rng: np.random.Generator,
) -> None:
    """A column of jittered dots per table column, one per adult, its median a bar.

    The columns in `red` are drawn red, the others grey.
    """
    for i, k in enumerate(columns):
        v = rows[k].to_numpy(float)
        ax.scatter(
            np.full(len(v), i) + rng.uniform(-0.1, 0.1, len(v)),
            v,
            s=18,
            facecolor=RED if k in red else "0.6",
            edgecolor="0.25",
            linewidth=0.4,
            zorder=2,
        )
        ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color="0.15", lw=1.7, zorder=3)


def plot_channel_check(
    rows: pd.DataFrame,
    first_adult: dict[str, dict[str, float]],
    mouse: str,
    gene: str,
    save: Path | None = None,
) -> plt.Figure:
    """sep_channel_check.png: the ranges, the first adult, the Gria1 correlations.

    `rows` is sep_channel_check.csv, `first_adult` one adult's log2 channel means
    ({channel: {structure: value}}, the channels sig, auto and sep) and `gene` the
    gene the channels are compared with. One jitter generator for the figure,
    drawn from in panel order.
    """
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.0))
    rng = np.random.default_rng(0)

    # how much each channel varies across the brain, a dot per adult
    ax = axes[0]
    ranges = ("range_nano", "range_auto", "range_sep")
    channel_dots(ax, rows, ranges, ("range_sep",), rng)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["nano", "autofluo", "SEP"], fontsize=8)
    ax.set_ylabel("p90 - p10 across structures (log2)", fontsize=8)
    ax.set_ylim(bottom=0)
    ax.set_title("how much each channel varies\nacross the brain", fontsize=9)

    # the green channel against autofluorescence and against nano, the first adult
    common = sorted(set.intersection(*[set(v) for v in first_adult.values()]))
    auto = np.array([first_adult["auto"][s] for s in common])
    sep = np.array([first_adult["sep"][s] for s in common])
    nano = np.array([first_adult["sig"][s] for s in common])
    for ax, x, name in ((axes[1], auto, "autofluorescence"), (axes[2], nano, "nano")):
        ax.scatter(x, sep, s=10, facecolor="0.55", edgecolor="0.25", linewidth=0.3)
        ax.set_xlabel(f"log2 {name}", fontsize=8)
        ax.set_ylabel("log2 SEP", fontsize=8)
    rho_auto = spearmanr(auto, sep).statistic
    axes[1].set_title(f"{mouse}\nSEP against autofluo, rho = {rho_auto:+.2f}", fontsize=9)
    rho_nano = spearmanr(nano, sep).statistic
    axes[2].set_title(f"SEP against nano, rho = {rho_nano:+.2f}", fontsize=9)

    # each channel's correlation with Gria1 mRNA, a dot per adult
    ax = axes[3]
    keys = ("rho_nano_gria", "rho_auto_gria", "rho_sep_gria", "rho_sepresid_gria")
    channel_dots(ax, rows, keys, ("rho_sep_gria", "rho_sepresid_gria"), rng)
    ax.axhline(0, color="0.8", lw=0.7, zorder=0)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(["nano", "autofluo", "SEP", "SEP minus\nautofluo"], fontsize=8)
    ax.set_ylabel(f"Spearman with {gene} mRNA", fontsize=8)
    ax.set_title(
        "a channel that reports the tagged receptor\nshould follow Gria1 at least as "
        "nano does",
        fontsize=9,
    )
    for ax in axes:
        tidy(ax)
    fig.suptitle(
        "The green channel in fixed, mounted tissue: one dot per adult, "
        "structure means with no denominator anywhere",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return saved_working(fig, save)
