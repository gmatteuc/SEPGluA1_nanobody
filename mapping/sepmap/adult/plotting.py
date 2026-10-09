"""The figures of the adult map's analyses: part 1 of the ISH line and the green channel.

The guided figures of part 1, how much of the map Gria1 and synapse density leave
(03, 03s1 to 03s4, 04, 11s1), and of analysis 5, what the green channel reports
(14, 14s), drawn as the ISH line's: number and question at the top, one
line to take from a main figure, grey notes on how to read it at the foot, a
takeaway that follows its numbers; their shared pieces (heading, footer, panel
titles, the scatter of structures by group of divisions) are those of
ish.plotting.shared, so the walk reads as one set, and they are saved as PNG and EPS
at its dpi. Beside them, the working figure of analysis 5, sep_channel_check.png,
PNG only at 200 dpi, which shows the step's numbers as it runs. The modules that
compute draw nothing: every function here takes tables and returns the figure.

Called by run_synaptome.py, run_density_markers.py and run_sep_channel_check.py
(steps 21, 22 and 29) and by adult.beyond_figures (step 27).
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
from sepmap.ish.numbers import points as points_text
from sepmap.ish.plotting.shared import (
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
    MID_BLUE,
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
    save_figure,
    tidy,
)

# dpi of the working figure, which only the run's log points to
WORKING_DPI = 200


# ===== 03s4 The measured synapse density (adult.synaptome) =====


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
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="lower right")
    tidy(ax)


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


def plot_synaptome_detail(
    density: pd.DataFrame,
    coverage: pd.DataFrame,
    agreement: pd.DataFrame,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 03s4: the measured synapse density, its coverage, how it compares.

    A: per division, the structures of analysis 4's fit with a measured PSD95 density.
    B: the left hemisphere against the right, the one check of a one-mouse map. C:
    each density's Spearman with the density term, Gria1, the nano map and
    autofluorescence. `density`, `coverage` and `agreement` are the tables of
    run_synaptome. The density against the density term is figure 03 D.
    """
    table = density.set_index("structure")
    n_fit, n_measured = coverage_counts(table)
    measured = table[table["in_fit"] & table["measured"]]
    groups = synaptome_groups(table)

    fig = plt.figure(figsize=(16, 6.2))
    ax = fig.add_axes([0.07, 0.2, 0.22, 0.56])
    coverage_panel(ax, coverage)
    panel_title(
        ax,
        "A",
        "Structures of the fit with a measured density",
        f"{n_measured} of {n_fit} ({n_measured / n_fit:.0%})",
    )

    ax = fig.add_axes([0.38, 0.2, 0.22, 0.56])
    n, rho = rank_scatter(
        ax,
        measured["psd95_left"],
        measured["psd95_right"],
        groups,
        ("left hemisphere (rank)", "right hemisphere (rank)"),
    )
    ax.legend(handles=group_handles(), fontsize=7, frameon=False, loc="upper left")
    panel_title(
        ax, "B", "One mouse: left hemisphere against right", f"rho {rho:+.2f}, n = {n}"
    )

    ax = fig.add_axes([0.76, 0.2, 0.21, 0.56])
    agreement_dots(ax, agreement)
    panel_title(ax, "C", "How each density agrees with the other maps")

    # the title, the numbers and how to read it
    heading(
        fig,
        "synaptome_detail",
        f"PSD95 puncta cover {n_measured} of the {n_fit} structures of the fit "
        f"({n_measured / n_fit:.0%}): the yardstick of the density genes "
        f"({figure_ref('density_markers')}) and a check row of part 1 "
        f"({figure_ref('beyond_budget')}), not its term",
    )
    footer(
        fig,
        [
            "Zhu et al. 2018, one adult male mouse, as Hansen et al. share it. PSD95 "
            f"density: the mean of the {len(synaptome.PSD95)} subtypes whose puncta hold "
            "PSD95, each divided by its largest value as shared; never a punctum's",
            "intensity or size. A structure is measured when it or one of its parts was "
            "sampled; a region above several structures is never spread onto them. "
            "Every rho is Spearman over the structures both maps have.",
            f"PSD95 against the density term on the declared structures: "
            f"{figure_ref('beyond')} D.",
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
    """D: each composite's agreement with PSD95, whole brain and inside divisions.

    Per composite: on the held-out halves, the median with its 95% range, whole
    brain (dot) and inside divisions (square, below it); on the full set, open
    markers of the same shapes. The labels give each composite's structures, and its
    held-out rho on the structures every composite has.
    """
    rows = comparison.set_index("composite")
    names = [c for c in COMPOSITE_LABELS if c in rows.index]
    for k, name in enumerate(names):
        r = rows.loc[name]
        colour = DENSITY_BLUE if name == "chosen" else DARK_GREY
        kinds = (
            (k - 0.13, "held_out", "rho_full_set", "o"),
            (k + 0.17, "within_held_out", "rho_within_full_set", "s"),
        )
        for y, key, full, marker in kinds:
            lo, hi = r[f"rho_{key}_lo"], r[f"rho_{key}_hi"]
            ax.plot([lo, hi], [y, y], color=colour, lw=2 if marker == "o" else 1.2)
            ax.scatter(
                r[f"rho_{key}_median"], y, s=40, marker=marker, color=colour, zorder=3
            )
            ax.scatter(
                r[full],
                y,
                s=30,
                marker=marker,
                facecolors="none",
                edgecolors=colour,
                linewidths=1.0,
                zorder=3,
            )
    labels = []
    for name in names:
        r = rows.loc[name]
        genes = f"\n{r['genes']}" if name == "chosen" else ""
        labels.append(
            f"{COMPOSITE_LABELS[name]}{genes}\n{int(r['n_structures'])} structures "
            f"({int(r['n_held_out'])} a half);\non the {int(r['n_common'])} all have: "
            f"{r['rho_common_held_out_median']:+.2f}"
        )
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color="0.6", lw=0.6)

    # room under the last row for the legend
    ax.set_ylim(len(names) + 0.4, -0.5)
    ax.set_xlim(-0.6, 1.0)
    ax.set_xlabel("Spearman rho with PSD95 punctum density", fontsize=8)
    handles = [
        plt.Line2D([], [], color=DARK_GREY, lw=2, marker="o", ms=6),
        plt.Line2D([], [], color=DARK_GREY, lw=1.2, marker="s", ms=5.5),
        plt.Line2D([], [], ls="", marker="o", ms=5.5, mfc="none", mec=DARK_GREY),
    ]
    ax.legend(
        handles,
        [
            "whole brain: held-out halves, median and 95% range",
            "inside divisions: the same",
            "full set (open; optimistic for the chosen)",
        ],
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
    # room right of the bars for the longest list of genes
    ax.set_xlim(0, counts.max() * 6.5)
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
    rule: dict,
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
    density_markers.exclusion_reasons, `rule` the numbers the footer and D name
    (min_coverage, n_reference: the declared structures where Gria1 is measured,
    min_division).
    """
    chosen = list(agreement.loc[agreement["chosen"], "symbol"])
    min_coverage, n_reference = rule["min_coverage"], rule["n_reference"]
    min_division = rule["min_division"]
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
    n_held = rows["n_held_out"].astype(int)
    panel_title(
        ax,
        "D",
        "Agreement with PSD95 on the other half",
        f"{n_held.min()} to {n_held.max()} structures a half; inside divisions, "
        f"each with {min_division}\nstructures or more, so their contrast does not "
        "enter it",
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
    third = selection.set_index("symbol").loc[chosen[-1], "share"]
    heading(
        fig,
        "density_markers",
        f"Chosen without the map: {', '.join(chosen[:-1])} and {chosen[-1]} (the third "
        f"on {third:.0%} of the halves). Re-run on random halves, the rule's choice "
        f"follows PSD95 at rho {mine['rho_held_out_median']:+.2f} on the other half "
        f"({mine['rho_held_out_lo']:+.2f} to {mine['rho_held_out_hi']:+.2f}), mostly "
        f"between divisions (inside them {mine['rho_within_held_out_median']:+.2f}); "
        f"Dlg4, Homer1 and Camk2a, excluded as AMPA-linked, "
        f"{first['rho_held_out_median']:+.2f} (inside "
        f"{first['rho_within_held_out_median']:+.2f})",
    )
    footer(
        fig,
        [
            f"The pool: genes annotated in mouse GO (release {go_release}) to the "
            "postsynaptic density or specialization, with two or more usable Allen "
            "experiments, measured in at least",
            f"{min_coverage:.0%} of the {n_reference} declared structures where Gria1 "
            "is measured; the genes that place or regulate AMPA receptors are excluded "
            "(red), as the surface side rather than density.",
            "PSD95: punctum density of one adult mouse (Zhu et al. 2018). No nano value "
            "is read.",
        ],
    )
    return saved(fig, save)


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

    # room under the lowest dots for their names
    ax.set_ylim(-0.12, 1.04)
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


def partition_bar(
    ax: plt.Axes,
    y: float,
    parts: dict[str, float],
    height: float,
    left_ci: tuple[float, float] | None = None,
) -> None:
    """One bar of the four parts of the reproducible map, a negative one hatched.

    A negative part (held out, a term that adds nothing to the other costs a little)
    is drawn from zero leftwards, hatched, so the stack of the others then reaches
    past 100% by as much. With `left_ci`, the share left carries its 95% interval
    inside its own segment.
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
            text = f"{PART_LABELS[key]}\n{parts[key]:.0%}"
            if key == "left" and left_ci is not None:
                text += f"\n(95% {left_ci[0]:.0%} to {left_ci[1]:.0%})"
            ax.text(
                a + w / 2,
                y,
                text,
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
    """B of figure 03: the reproducible map in four parts, what is left with its 95%.

    Gria1 only, shared, density only and left, each a share of what two halves of
    the cohort reproduce; the interval of what is left sits in its own segment, over
    structures, and under the bar the same over spatial blocks; a part too narrow to
    label is named above the bar.
    """
    parts = n["parts"]
    partition_bar(ax, 0, parts, height=0.6, left_ci=n["left_ci"])
    lo, hi = n["left_ci"]
    block_lo, block_hi = n["left_ci_blocks"]
    ax.text(
        1,
        -0.42,
        f"left, 95%: {lo:.0%} to {hi:.0%} over resampled structures, {block_lo:.0%} "
        f"to {block_hi:.0%} leaving out whole spatial blocks",
        fontsize=7.5,
        va="center",
        ha="right",
        color=DARK_GREY,
    )
    note = small_parts_note(parts)
    if note:
        ax.text(0, 0.36, note, fontsize=7.5, va="bottom", color=DARK_GREY)
    low = min(0.0, min(parts.values()) - 0.12)
    ax.set_xlim(low, 1)
    ticks = np.arange(0, 1.01, 0.2)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_ylim(-0.6, 0.5)
    ax.set_yticks([])
    ax.set_xlabel(
        "share of the map's reproducible pattern, on structures the fit has not seen"
    )
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)


def floor_panel(ax: plt.Axes, calibration: pd.DataFrame, n: dict) -> None:
    """C of figure 03, its upper part: nano and the floor on the same structures.

    On the calibration's structures: a dot per draw of made-up adults whose map is
    only Gria1 and synapse density, the median as a bar; the nano map read the same
    way, the mean over the two halves of the Allen experiments its predictors come
    from, an arrow from the floor's median to it.
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
    ax.scatter([left], [1], s=90, color=NANO, linewidths=0, zorder=3)
    ax.text(left, 0.7, f"{left:.0%}", ha="center", va="bottom", fontsize=8.5)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(
        [
            f"the floor: a map made only\nof Gria1 and synapse density\n(dots: "
            f"{len(floor)} made-up cohorts;\nbar: their median)",
            "the nano map",
        ],
        fontsize=8.5,
    )
    ax.set_ylim(1.45, -0.75)
    ax.set_xlim(0, max(0.5, float(floor.max()) + 0.05, left + 0.05))
    ticks = np.arange(0, ax.get_xlim()[1] + 1e-9, 0.1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map left")
    tidy(ax)


def difference_panel(ax: plt.Axes, n: dict) -> None:
    """C of figure 03, its lower part: nano minus the floor in points, on its own axis.

    The point with its 95% over resampled structures (thick) and leaving out whole
    spatial blocks (thin), and the zero line.
    """
    diff = 100 * n["minus_floor"]
    lo, hi = (100 * v for v in n["minus_floor_ci"])
    block_lo, block_hi = (100 * v for v in n["minus_floor_ci_blocks"])
    words = [points_text(v / 100) for v in (diff, lo, hi, block_lo, block_hi)]
    ax.plot([block_lo, block_hi], [0, 0], color=MID_GREY, lw=1.2, zorder=2)
    ax.plot([lo, hi], [0, 0], color="0.1", lw=3, solid_capstyle="butt", zorder=3)
    ax.scatter([diff], [0], s=60, marker="D", color="0.1", zorder=4)
    ax.text(
        diff,
        0.42,
        f"{words[0]} points (95% {words[1]} to {words[2]} over structures, "
        f"{words[3]} to {words[4]} over spatial blocks)",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_yticks([0])
    ax.set_yticklabels(["nano minus the floor"], fontsize=8.5)
    ax.set_ylim(-0.6, 1.1)
    ax.set_xlim(min(-5.0, block_lo - 5), max(50.0, block_hi + 5))
    ax.set_xlabel("difference, points of the reproducible map")
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

        # the mean's value to the left of its bar, clear of the dashed line
        ax.text(
            values.mean() - 0.002,
            y + 0.36,
            f"mean {values.mean():.3f}",
            ha="right",
            va="bottom",
            fontsize=8,
        )
    ax.axvline(n["implied"], color="0.1", lw=1, ls=(0, (3, 2)), zorder=1)
    ax.text(
        n["implied"] - 0.003,
        -0.62,
        f"{n['implied']:.3f}: what the leftover\nreplicates at anyway, given the\n"
        f"map's reliability and the fit (in-sample R² {n['r2']:.2f})",
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
    """Nano minus the floor with both intervals: '+27 points (95% +15 to +40; ...)'."""
    lo, hi = n["minus_floor_ci"]
    block_lo, block_hi = n["minus_floor_ci_blocks"]
    return (
        f"{points_text(n['minus_floor'])} points (95% {points_text(lo)} to "
        f"{points_text(hi)}; over spatial blocks {points_text(block_lo)} to "
        f"{points_text(block_hi)})"
    )


def density_check_panel(ax: plt.Axes, markers: dict, groups: dict[str, str]) -> None:
    """D of figure 03: measured PSD95 punctum density against the density term.

    `markers` is adult.density_markers.load_validation: PSD95 and the chosen genes'
    mean rank over the declared structures where both exist; `groups` each
    structure's group of divisions. Ranks 0 to 1, as in A.
    """
    pair = pd.DataFrame(dict(x=markers["psd95"], y=markers["composite"])).dropna()
    group = [groups.get(s, "other grey matter") for s in pair.index]
    scatter_groups(
        ax, ranks01(pair["x"].to_numpy()), ranks01(pair["y"].to_numpy()), group
    )
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.04, 1.04)
    ax.set_xlabel("PSD95 punctum density, measured, rank (0 low, 1 high)")
    ax.set_ylabel(f"mean rank of {', '.join(markers['chosen'])}, rank")
    ax.legend(
        handles=group_handles(),
        fontsize=7,
        frameon=False,
        loc="upper left",
        title=f"all {len(pair)} structures (the full set)",
        title_fontsize=7,
    )
    tidy(ax)


def density_check_text(markers: dict) -> str:
    """The numbers of figure 03 D: the rule held out, inside divisions, the full set."""
    chosen = markers["comparison"].loc["chosen"]
    first = markers["comparison"].loc["first_proposal"]
    text = (
        f"the rule re-run on {markers['n_halves']} random halves: rho "
        f"{chosen['rho_held_out_median']:+.2f} on the other half (95% "
        f"{chosen['rho_held_out_lo']:+.2f} to {chosen['rho_held_out_hi']:+.2f}), "
        f"inside divisions {chosen['rho_within_held_out_median']:+.2f}; full set "
        f"{chosen['rho_full_set']:+.2f}, n = {int(chosen['n_structures'])} "
        "(optimistic); Dlg4, Homer1, Camk2a (excluded: AMPA-linked) "
        f"{first['rho_held_out_median']:+.2f}, inside divisions "
        f"{first['rho_within_held_out_median']:+.2f}"
    )
    return textwrap.fill(text, 58)


def beyond_takeaway(n: dict) -> str:
    """The line under figure 03's title: what is left, and how it stands to the floor.

    The words follow the interval of nano minus the floor over structures; the one
    over spatial blocks is given beside it.
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
    """How the check rows stand to their own floors, at the point and over each 95%.

    'check rows above their own floor: 10 of 10, 10 over their 95% over structures,
    9 over spatial blocks'.
    """
    checks = n["checks"]
    point = int((checks["nano_minus_floor"] > 0).sum())
    over = int((checks["nano_minus_floor_lo"] > 0).sum())
    blocks = int((checks["nano_minus_floor_blocks_lo"] > 0).sum())
    return (
        f"check rows above their own floor: {point} of {len(checks)}, {over} over their "
        f"95% over structures, {blocks} over spatial blocks"
    )


# what the floor holds and lacks, in the words of figures 03 and 03s1, two lines
FLOOR_ERRORS = (
    "The floor errs both ways: it lacks the mismatch of Allen's P56 mice with these "
    "brains (age, strain, grid, registration), which makes it low,",
    "and it holds the disagreement of two halves of the Allen experiments, where nano, "
    "read with one half, meets one, which makes it high.",
)


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
    fig = plt.figure(figsize=(15, 13))
    heading(fig, "beyond", beyond_takeaway(n))

    # A: the map against the main model's held-out prediction
    ax = fig.add_axes([0.07, 0.57, 0.27, 0.3])
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
    ax = fig.add_axes([0.45, 0.65, 0.5, 0.15])
    partition_panel(ax, n)
    weights = n["weights"]
    panel_title(
        ax,
        "B",
        "The reproducible map in four parts",
        f"nano ~ {n['abundance']} + synapse density, two straight terms; weights "
        f"(z-scored) {', '.join(f'{k} {v:+.2f}' for k, v in weights.items())}",
    )

    # C: the leftover beside the floor, and their difference on an axis of its own
    ax = fig.add_axes([0.2, 0.27, 0.29, 0.19])
    floor_panel(ax, calibration, n)
    numbers_c = (
        f"on the {n['cal_n']} structures both halves of the Allen experiments "
        "measure; the floor lacks the P56 mismatch (low) and counts two Allen halves "
        f"(high); {checks_above(n)} ({figure_ref('beyond_budget')})"
    )
    panel_title(
        ax,
        "C",
        "Is the leftover more than Allen-to-Allen mismatch?",
        textwrap.fill(numbers_c, 64),
    )
    ax = fig.add_axes([0.2, 0.16, 0.29, 0.04])
    difference_panel(ax, n)

    # D: the density term against the measured synapse density
    ax = fig.add_axes([0.64, 0.16, 0.26, 0.29])
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
            "cohort reproduce. Gria1 only, density only",
            "and shared come from Gria1 alone, density alone and both. The floor is a "
            "map made only of Gria1 and synapse density, from half of each gene's Allen "
            "experiments, with ten made-up",
            "adults as noisy as ours, predicted from the other half.",
            *FLOOR_ERRORS,
            "Synapse density is the mean rank of three postsynaptic genes, chosen by "
            "their agreement with PSD95 punctum density (Zhu et al. 2018, one adult "
            "mouse) by a rule that reads no nano value,",
            "the AMPA-linked genes left out; held out, the rule chooses on a random half "
            "of the structures and is tested on the other; inside divisions, the mean "
            "rho within each division",
            f"of {markers['min_division']} structures or more, so the contrast between "
            "divisions does not enter it.",
            f"In detail: the check rows, {figure_ref('beyond_budget')}; the seven "
            f"controls, {figure_ref('beyond_controls')}; the choice of the genes, "
            f"{figure_ref('density_markers')}; the measured density, "
            f"{figure_ref('synaptome_detail')}.",
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
        ("both, curved (a check row)", shares["curved"], None, MID_GREY),
        ("autofluorescence alone", shares["autofluorescence"], None, AUTO),
        (
            f"control F: {n['f_genes']} genes' components\n(a flexible bound; on its "
            f"own floor\nnano leaves {n['f_nano'][0]:.0%} to {n['f_nano'][1]:.0%})",
            n["f_share"],
            None,
            MID_GREY,
        ),
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
    """G of figure 03s1: the four parts with their 95% intervals.

    Over structures (thick) and over spatial blocks (thin, below); a part below zero
    is drawn below zero, as computed; the weights are written under the panel's
    title.
    """
    parts, ci, blocks = n["parts"], n["parts_ci"], n["parts_ci_blocks"]
    keys = list(PART_LABELS)
    for k, key in enumerate(keys):
        ax.barh(k, parts[key], color=PART_COLOURS[key], height=0.62, zorder=2)
        lo, hi = ci[key]
        block_lo, block_hi = blocks[key]
        ax.plot([lo, hi], [k, k], color="0.1", lw=1.6, zorder=3)
        ax.plot([block_lo, block_hi], [k + 0.22] * 2, color="0.1", lw=0.7, zorder=3)
        ax.text(
            max(hi, block_hi, 0) + 0.01,
            k,
            f"{parts[key]:+.0%} ({lo:+.0%} to {hi:+.0%}; {block_lo:+.0%} to "
            f"{block_hi:+.0%})",
            va="center",
            fontsize=7.5,
        )
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([PART_LABELS[k] for k in keys], fontsize=8.5)
    ax.set_ylim(len(keys) - 0.4, -0.6)
    lows = [lo for lo, _ in ci.values()] + [lo for lo, _ in blocks.values()]
    low = min(0.0, min(lows) - 0.05)
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
    experiments, and nano with each half's predictors and with both; the known map
    read with its own half, the animals' noise alone, is in the panel's title.
    """
    calibration = calibration[calibration["folds"] == "random"]
    calibration = calibration[calibration["map"].isin([FLOOR_MAP, "nano"])]
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
    exists, nano minus it with its intervals over structures and over spatial blocks,
    and what its floor misses as measured once (once).
    """
    folds = n["folds"].copy()
    folds["nano_minus_floor"] = np.nan
    folds["once"] = ""
    main = folds["kind"] == "main"
    folds.loc[main, "lo"] = n["left_ci"][0]
    folds.loc[main, "hi"] = n["left_ci"][1]
    folds.loc[main, "nano_minus_floor"] = n["minus_floor"]
    folds.loc[main, "nano_minus_floor_lo"] = n["minus_floor_ci"][0]
    folds.loc[main, "nano_minus_floor_hi"] = n["minus_floor_ci"][1]
    folds.loc[main, "nano_minus_floor_blocks_lo"] = n["minus_floor_ci_blocks"][0]
    folds.loc[main, "nano_minus_floor_blocks_hi"] = n["minus_floor_ci_blocks"][1]
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
        "nano_minus_floor_blocks_lo",
        "nano_minus_floor_blocks_hi",
        "once",
    ]
    out = pd.concat([folds.reindex(columns=columns), checks.reindex(columns=columns)])
    out["once"] = out["once"].fillna("")
    return out.reset_index(drop=True)


def margin_label(r: pd.Series) -> str:
    """A row's nano minus its floor in points, with its intervals when it has them.

    '+27 (+15 to +40; +10 to +45)': over structures, then over spatial blocks.
    """
    text = points_text(r["nano_minus_floor"])
    if np.isfinite(r["nano_minus_floor_lo"]):
        lo, hi = r["nano_minus_floor_lo"], r["nano_minus_floor_hi"]
        text += f" ({points_text(lo)} to {points_text(hi)}"
        if np.isfinite(r["nano_minus_floor_blocks_lo"]):
            block_lo = r["nano_minus_floor_blocks_lo"]
            block_hi = r["nano_minus_floor_blocks_hi"]
            text += f"; {points_text(block_lo)} to {points_text(block_hi)}"
        text += ")"
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
            margin_ax.text(1, i, "no floor", fontsize=7, va="center", color=MID_GREY)
            continue

        # in points: the 95% over structures thick, over spatial blocks thin below
        point = 100 * r["nano_minus_floor"]
        right = point
        if np.isfinite(r["nano_minus_floor_lo"]):
            lo, hi = 100 * r["nano_minus_floor_lo"], 100 * r["nano_minus_floor_hi"]
            margin_ax.plot([lo, hi], [i, i], color=MID_GREY, lw=2, zorder=1)
            right = hi
        if np.isfinite(r["nano_minus_floor_blocks_lo"]):
            block_lo = 100 * r["nano_minus_floor_blocks_lo"]
            block_hi = 100 * r["nano_minus_floor_blocks_hi"]
            margin_ax.plot(
                [block_lo, block_hi], [i + 0.25] * 2, color="0.25", lw=0.8, zorder=1
            )
            right = max(right, block_hi)

        # open where the floor misses what is measured once
        face = "none" if r["once"] else colour
        margin_ax.scatter(
            point, i, s=34, facecolors=face, edgecolors=colour, linewidths=1.2, zorder=2
        )
        margin_ax.text(
            right + 1.5, i, margin_label(r), fontsize=7.5, va="center", color=colour
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
    lows = table[["nano_minus_floor_lo", "nano_minus_floor_blocks_lo"]].min(axis=1)
    highs = table[["nano_minus_floor_hi", "nano_minus_floor_blocks_hi"]].max(axis=1)
    margin_ax.set_xlim(
        min(-10.0, 100 * lows.min() - 5), max(30.0, 100 * highs.max() + 30)
    )
    margin_ax.set_xlabel("nano minus its own floor, points of the reproducible map")
    ticks = [
        t
        for t in left_ax.get_xticks()
        if left_ax.get_xlim()[0] <= t <= left_ax.get_xlim()[1]
    ]
    left_ax.set_xticks(ticks)
    left_ax.set_xticklabels([f"{t:.0%}" for t in ticks])


def check_rows_line(n: dict) -> str:
    """The line under figure 03s1's title: the check rows against their floors."""
    checks = n["checks"].set_index("key")
    low = checks["nano_minus_floor"].idxmin()
    high = checks["nano_minus_floor"].idxmax()
    above = checks_above(n)
    return (
        f"{above[:1].upper()}{above[1:]}; nano minus its floor runs from "
        f"{100 * checks.loc[low, 'nano_minus_floor']:+.0f} points ({ROW_LABELS[low]}) "
        f"to {100 * checks.loc[high, 'nano_minus_floor']:+.0f} ({ROW_LABELS[high]}); "
        f"the main model leaves {n['left']:.0%}, {margin_text(n)} above its floor"
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
        "95% over structures thick, over spatial blocks thin; open: a floor that "
        "misses\nwhat is measured once (PSD95, three of the 11 markers); folds of "
        f"spatial blocks on\nthe {n['cal_n']} structures of the calibration: nano "
        f"{n['blocks']['nano']:.0%}, floor {n['blocks']['floor']:.0%}",
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
    weights = ",\n".join(
        f"{k} {v:+.2f} ({n['weights_ci'][k][0]:+.2f} to {n['weights_ci'][k][1]:+.2f}; "
        f"{n['weights_ci_blocks'][k][0]:+.2f} to {n['weights_ci_blocks'][k][1]:+.2f})"
        for k, v in n["weights"].items()
    )
    panel_title(
        ax,
        "G",
        "The four parts, 95% over structures (thick) and over spatial blocks (thin)",
        f"weights (z-scored): {weights}",
    )

    # H: the calibration; I: the replication
    ax = axes_at(fig, 0.07, 2.2, 0.38, 3.9)
    calibration_panel(ax, calibration, n)
    floor = n["floor"]
    panel_title(
        ax,
        "H",
        "The same model on a map whose answer is known",
        textwrap.fill(
            f"on {n['cal_n']} structures: the floor {floor['left_median']:.0%} "
            f"({floor['left_lo']:.0%} to {floor['left_hi']:.0%}; read with its own half "
            f"{n['own_half']:.0%}); nano {n['nano_cal_left']:.0%}, minus the floor "
            f"{margin_text(n)}",
            80,
        ),
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
            "five folds) over the ceiling; intervals are 95%, by jackknife over "
            "structures and over spatial blocks.",
            "Gria1 only = both - density alone; density only = both - Gria1 alone; "
            "shared = Gria1 alone + density alone - both; left = 1 - both. Held out, a "
            "part can come out negative; it is kept as computed, never clipped.",
            "The calibration builds a map made only of Gria1 and synapse density from "
            "half of each gene's Allen experiments, gives it ten made-up adults as noisy "
            "as ours, and predicts it from the other half; every check row has a floor",
            "of its own, built the same way with its own model.",
            *FLOOR_ERRORS,
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


# the readings of control G as figure 03s2 names them
READING_NAMES = {
    "zref": "zref, the declared reference",
    "zref_stored": "zref, the 17-brain reference",
    "cref": "over the isocortex (cref)",
    "subref": "over the subcortex (subref)",
    "ratio": "nano / autofluorescence",
    "sepratio": "nano / SEP",
}


def control_title(
    ax: plt.Axes, letter: str, controls: pd.DataFrame, number: str | None = None
) -> None:
    """A control's panel title: its question, its verdict and its number.

    The number is controls.csv's, R2 written R², unless `number` replaces it.
    """
    row = controls[controls["control"].str[0] == letter].iloc[0]
    verdict = "passes" if row["verdict"] == "pass" else "does not pass"
    if number is None:
        number = row["number"].replace("R2", "R²")
    panel_title(
        ax,
        letter,
        CONTROL_QUESTIONS[letter],
        f"{verdict}: {textwrap.fill(number, 52)}",
    )


def leftover_scatter(ax: plt.Axes, x: np.ndarray, y: np.ndarray, xlabel: str) -> None:
    """The leftover of each structure against one of its properties."""
    ok = np.isfinite(x)
    ax.scatter(x[ok], y[ok], s=35, color=DARK_GREY, linewidths=0, alpha=0.85, zorder=2)
    ax.axhline(0, color=MID_GREY, lw=0.7, zorder=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("leftover (ranks above prediction)")
    tidy(ax)


def curvature_panel(
    ax: plt.Axes, scores: tuple[float, float, float], gain: float
) -> None:
    """E: the held-out R² of the model straight, curved to cubes and to fifth powers.

    The control's line, the straight model's R² plus `gain`, dashed.
    """
    labels = ["straight\n(the model)", "x, x², x³\n(row curved)", "to x⁵"]
    colours = [DARK_GREY, MID_GREY, MID_GREY]
    ax.bar(range(3), scores, color=colours, width=0.6)
    for k, value in enumerate(scores):
        ax.text(k, value + 0.01, f"{value:.3f}", ha="center", va="bottom", fontsize=7.5)
    line = scores[0] + gain
    ax.axhline(line, color=RED, lw=0.9, ls=(0, (3, 2)))

    # the line's name right of the bars, clear of their values
    ax.text(
        2.4,
        line,
        f"the control's line:\nstraight + {gain}",
        fontsize=7,
        color=RED,
        ha="left",
        va="center",
    )
    ax.set_xlim(-0.5, 3.3)
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
        ("nano", "the nano map\n(predictors from\nhalf A and half B)", RED),
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
    ax.set_yticklabels([READING_NAMES.get(r, r) for r in readings["reading"]], fontsize=8)
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
    gain = data["thresholds"]["curvature_gain"]
    fig = plt.figure(figsize=(16, 10))
    if failed:
        line = (
            f"{len(passed)} of {len(controls)} controls pass; "
            + ", ".join(failed)
            + (" does not" if len(failed) == 1 else " do not")
        )
        if "E" in failed:
            straight, curved, _ = data["curvature"]
            line += (
                f": curving buys {curved - straight:+.3f} of R², more than the {gain} "
                "the control allows; the curved check row "
                f"({figure_ref('beyond_budget')}) shows how much of the leftover is "
                "curvature"
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
    curvature_panel(ax, data["curvature"], gain)
    control_title(ax, "E", controls)

    ax = fig.add_subplot(grid[1, 4:8])
    gene_space_panel(ax, data["f_calibration"])
    control_title(ax, "F", controls)

    ax = fig.add_subplot(grid[1, 8:12])
    readings = data["readings"]
    readings_panel(ax, readings, data["thresholds"]["readings_replication"])
    replication = readings["leftover_replication"]
    control_title(
        ax,
        "G",
        controls,
        f"under every reading the leftover replicates at {replication.min():.2f} to "
        f"{replication.max():.2f}",
    )
    footer(
        fig,
        [
            "How to read: each panel tries one way the leftover of the main model "
            f"({figure_ref('beyond')}) could be an artefact; the line each must not "
            "cross is fixed in settings.toml [beyond_controls] (dashed in C, E and G). "
            "In C the median pair is red.",
            "E and F score on structures the fit has not seen; F gives the model the "
            "components of every gene measured in all the structures, picked inside "
            "each training fold, and sets nano beside its own floor.",
        ],
    )
    return saved(fig, save)


# ===== 04 and 11s1 Where the leftover lives, and every gene against it =====

# the structures shown in the bars, each way
N_LEFTOVER_BARS = 12

# the genes drawn against the leftover: this many of the highest rho and of the lowest,
# and these named in any case
N_LEFTOVER_GENES = 10
LEFTOVER_NAMED = ("Cacng8", "Gria1", "Dlg2")

# the divisions by name, as the line under figure 04's title gives them
DIVISION_NAMES = {
    "Isocortex": "isocortex",
    "OLF": "olfactory areas",
    "HPF": "hippocampal formation",
    "CTXsp": "cortical subplate",
    "STR": "striatum",
    "PAL": "pallidum",
    "TH": "thalamus",
    "HY": "hypothalamus",
    "MB": "midbrain",
    "P": "pons",
    "MY": "medulla",
    "CB": "cerebellum",
}


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


def division_bars(ax: plt.Axes, divisions: pd.DataFrame) -> None:
    """C: the mean leftover of each division, with its 95% over resampled adults.

    `divisions` is residual_by_division.csv, highest mean first; a division whose
    interval lies above zero is red, below zero blue, across it grey.
    """
    y = np.arange(len(divisions))
    colours = [
        RED if r.lo > 0 else (MID_BLUE if r.hi < 0 else MID_GREY)
        for r in divisions.itertuples()
    ]
    ax.barh(y, divisions["mean"], color=colours, height=0.65, zorder=2)
    for k, r in enumerate(divisions.itertuples()):
        ax.plot([r.lo, r.hi], [k, k], color="0.1", lw=1.1, zorder=3)
    ax.axvline(0, color="0.3", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r.division} ({r.n_structures})" for r in divisions.itertuples()], fontsize=8
    )
    ax.invert_yaxis()
    span = float(divisions[["lo", "hi"]].abs().to_numpy().max()) * 1.15
    ax.set_xlim(-span, span)
    ax.set_xlabel("mean leftover of the division's structures (ranks)")
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
    """A: the genes highest and lowest with the leftover, the named ones, null bands.

    Cacng8, named for the leftover in advance, in red; the members of the AMPA
    receptor complex family tagged; a dotted line between the highest, the lowest
    and the named.
    """
    top = genes.head(N_LEFTOVER_GENES)
    bottom = genes.tail(N_LEFTOVER_GENES)
    shown = set(top["symbol"]) | set(bottom["symbol"])
    named = genes[genes["symbol"].isin(LEFTOVER_NAMED) & ~genes["symbol"].isin(shown)]
    show = pd.concat([top, bottom, named]).reset_index(drop=True)
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
    for edge in (len(top), len(top) + len(bottom)):
        if edge < len(show):
            ax.axhline(edge - 0.5, color="0.6", lw=0.6, ls=(0, (3, 2)))
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


def tissue_scatter(ax: plt.Axes, genes: pd.DataFrame, auto: pd.Series, q: float) -> float:
    """C: each gene's rho with the leftover against its rho with autofluorescence.

    `auto` is each gene's rho with the autofluorescence map (gene_ranking.csv); the
    genes past BH on the leftover dark, those above zero red. Returns the Spearman
    over the genes.
    """
    both = genes.set_index("symbol")[["rho", "q_all"]].join(auto.rename("auto"))
    both = both.dropna(subset=["auto"])
    past = both["q_all"] < q
    ax.scatter(both.loc[~past, "auto"], both.loc[~past, "rho"], s=10, color=LIGHT_GREY)
    for above, colour in ((False, DARK_GREY), (True, RED)):
        mine = both[past & ((both["rho"] > 0) == above)]
        ax.scatter(mine["auto"], mine["rho"], s=16, color=colour, linewidths=0)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.axvline(0, color="0.6", lw=0.6)
    ax.set_xlabel("rho with the autofluorescence map")
    ax.set_ylabel("rho with the leftover")
    tidy(ax)
    return float(spearmanr(both["auto"], both["rho"]).statistic)


def where_takeaway(residuals: pd.DataFrame, numbers: dict) -> str:
    """The line under figure 04's title: the divisions as a whole, then the structures."""
    above = [DIVISION_NAMES.get(d, d) for d in numbers["divisions_above"]]
    below = [DIVISION_NAMES.get(d, d) for d in numbers["divisions_below"]]
    top_up = ", ".join(residuals["acronym"].head(4))
    top_down = ", ".join(residuals["acronym"].tail(4)[::-1])
    return (
        f"Above prediction as a whole: {', '.join(above) or 'no division'}; below: "
        f"{', '.join(below) or 'no division'} ({numbers['between']:.0%} of the "
        f"leftover's variance lies between divisions); the structures most above: "
        f"{top_up}; most below: {top_down}"
    )


def plot_beyond_where(
    planes: list[dict],
    residuals: pd.DataFrame,
    divisions: pd.DataFrame,
    numbers: dict,
    t_max: float,
    span: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 04: where the leftover lives.

    `planes` holds per plane its labels and the map, prediction and leftover
    painted (adult.beyond_figures.plane_images); `residuals` is
    residual_by_structure.csv and `divisions` residual_by_division.csv; `t_max` the t
    at which a structure's bar turns black; `span` the leftover's colour scale on the
    planes, in ranks each way.
    """
    fig = plt.figure(figsize=(16, 14))
    clipped = int((residuals["residual"].abs() > span).sum())
    heading(fig, "beyond_where", where_takeaway(residuals, numbers))

    # A: the three planes
    grid = fig.add_gridspec(
        len(planes),
        3,
        left=0.04,
        right=0.58,
        top=0.87,
        bottom=0.37,
        hspace=0.12,
        wspace=0.04,
    )
    images = leftover_planes(fig, grid, planes, span)
    cax = fig.add_axes([0.06, 0.335, 0.3, 0.01])
    cb = fig.colorbar(images[0], cax=cax, orientation="horizontal")
    cb.ax.set_xlim(0, 1)
    cb.set_label("rank among the structures of the fit (0 low, 1 high)", fontsize=8)
    cax = fig.add_axes([0.42, 0.335, 0.15, 0.01])
    cb = fig.colorbar(images[1], cax=cax, orientation="horizontal", extend="both")
    cb.set_label(
        f"leftover (ranks): red above prediction;\nclipped at ±{span:g} "
        f"({clipped} structures beyond)",
        fontsize=8,
    )
    fig.text(
        0.04,
        0.895,
        f"A.  The map, the prediction and the leftover on {len(planes)} coronal "
        "planes\nflat grey: structures the fit does not use",
        fontsize=9,
        va="bottom",
    )

    # B: the structures with the largest leftovers
    ax = fig.add_axes([0.8, 0.08, 0.17, 0.79])
    leftover_bars(ax, residuals, t_max)
    panel_title(
        ax,
        "B",
        f"The {N_LEFTOVER_BARS} largest leftovers each way",
        f"grey: how steady across the ten adults'\nown leftovers (t; black at {t_max:g})",
    )

    # C: the leftover by division
    ax = fig.add_axes([0.12, 0.07, 0.4, 0.17])
    division_bars(ax, divisions)
    panel_title(
        ax,
        "C",
        "The leftover by division",
        f"mean over its structures (in brackets), 95% over resampled adults; "
        f"{numbers['between']:.0%} of the leftover's variance lies between divisions",
    )
    footer(
        fig,
        [
            "How to read: the leftover is the nano rank minus the rank Gria1 and synapse "
            "density predict (the straight model fitted on every structure), on the "
            "structures of the fit; the bars' grey says how steady it is across the "
            "ten adults' own leftovers.",
            "What would mean what: a division above or below as a whole is a contrast "
            "the two terms do not reproduce; a structure far from prediction in every "
            "adult is where the departure sits; a claim that one structure stands out "
            "needs a null of its own, which is not drawn here.",
        ],
    )
    return saved(fig, save)


def plot_leftover_genes(
    genes: pd.DataFrame,
    sets: pd.DataFrame,
    auto: pd.Series,
    numbers: dict,
    n_surrogates: int,
    t_max: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 11s1: whether any gene's map, or any gene set, follows the leftover.

    `genes` and `sets` are leftover_genes.csv and leftover_sets.csv, `auto` each
    gene's rho with the autofluorescence map; `t_max` the rho over its SD across
    resampled adults at which a gene's bar turns black.
    """
    n = numbers
    q = n["q"]
    fig = plt.figure(figsize=(16, 13))
    tested = sets[sets["tested"]]
    passed = tested[tested["q"] < q]
    heading(
        fig,
        "leftover_genes",
        f"{n['n_genes']} genes against the leftover of the main model "
        f"({figure_ref('beyond_budget')}): {n['n_pass']} past its null at BH q < {q}, "
        f"{n['n_pass_negative']} of them below zero (exploratory; {n['n_p05']} below "
        "p 0.05 before correction); gene sets past it: "
        + (", ".join(passed["gene_set"]) if len(passed) else "none"),
    )
    subunits = {
        s for s, t in zip(genes["symbol"], genes["gene_sets"]) if "subunits" in str(t)
    }
    ax = fig.add_axes([0.21, 0.15, 0.23, 0.73])
    leftover_gene_bars(ax, genes, subunits, t_max, q)
    cacng8 = genes.set_index("symbol").loc["Cacng8"]
    panel_title(
        ax,
        "A",
        f"The {N_LEFTOVER_GENES} genes highest and lowest with the leftover, and the "
        "named ones",
        f"pale blue: 95% of rho with the leftover's {n_surrogates} surrogates, each "
        "through\nthe same fit; grey: rho over its SD across resampled adults (black "
        f"at {t_max:g});\nbold: past BH; Cacng8 {cacng8['rho']:+.2f} "
        f"({p_text(cacng8['p_spatial'], n_surrogates)}, rank {int(cacng8['rank_all'])})",
    )
    ax = fig.add_axes([0.55, 0.6, 0.42, 0.27])
    leftover_sets(ax, genes, sets)
    panel_title(
        ax,
        "B",
        "The gene sets of analysis 3 against the leftover",
        "pale blue: 95% of the set's median rho with the surrogates; spatial p, and "
        "q by BH over the tested sets",
    )
    ax = fig.add_axes([0.6, 0.15, 0.27, 0.29])
    rho = tissue_scatter(ax, genes, auto, q)
    panel_title(
        ax,
        "C",
        "Each gene with the leftover and with the tissue",
        f"over the genes rho {rho:+.2f}; past BH on the leftover: dark below zero, red "
        f"above;\n{n['pass_auto_positive']} of the {n['n_pass']} follow "
        f"autofluorescence positively, {n['pass_auto_bh']} past its own null",
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
            "C: the genes that run against the leftover are, in large part, those the "
            "autofluorescence map follows; whether the tissue or neuronal mRNA in "
            "general is behind them is not settled here. The glia set followed the",
            "leftover of an earlier model, a test not named in advance: an unplanned "
            "lead, not pursued.",
        ],
    )
    return saved(fig, save)


# ===== 14 and 14s What the green channel reports =====

# the three channels as figure 14 names and draws them: name, what it records, the
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
    """Figure 14: whether the green channel reports the tag or the tissue.

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
    """Figure 14s: whether the green channel reports the tag or the tissue, in detail.

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


def saved_working(fig: plt.Figure, path: Path | None) -> plt.Figure:
    """Save a working figure as PNG at WORKING_DPI when a path is given; return it."""
    if path is not None:
        save_figure(fig, Path(path), WORKING_DPI, eps=False)
    return fig


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
