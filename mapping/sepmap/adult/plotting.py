"""The figures of the adult map's analyses: part 1 of the ISH line and the green channel.

The guided figures of part 1, how much of the map Gria1 and synapse density leave
(03, 03s1, 03s2, 04, 11s1, 14, 14s), and of analysis 5, what the green channel
reports (15, 15s), drawn as the ISH line's: number and question at the top, one
line to take from a main figure, grey notes on how to read it at the foot, a
takeaway that follows its numbers; their shared pieces (heading, footer, panel
titles, the scatter of structures by group of divisions) are those of ish.plotting,
so the walk reads as one set, and they are saved as PNG and EPS at its dpi. Beside
them, the working figures of the steps (fig0 to fig6, E_regression, F_maps,
sep_channel_check.png), PNG only at 200 dpi, which show each step's numbers as it
runs. The modules that compute draw nothing: every function here takes tables and
returns the figure.

Called by the run scripts of steps 21 to 28 and by adult.beyond_figures.
"""

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle
from scipy.stats import rankdata, spearmanr

from sepmap.adult import synaptome
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
    NOTE_GREY,
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


# ===== 14 and 14s The measured synapse density (adult.synaptome) =====


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
    "markers": "marker mRNA composite",
    "psd_pc1": "psd_pc1 (PSD genes' mRNA)",
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


# the check rows of figure 14 C, top to bottom: the density measure of each, on the
# structures PSD95 covers (beyond.variants.csv)
DENSITY_ROWS = (
    ("panel", "mRNA panel (markers, psd_pc1)"),
    ("psd95", "PSD95 puncta"),
    ("psd95_and_panel", "PSD95 puncta and the mRNA panel"),
    ("sap102", "SAP102 puncta"),
    ("all_puncta", "every punctum"),
)


def synaptome_groups(table: pd.DataFrame) -> dict[str, str]:
    """{structure: group of divisions} of the synaptome's density table."""
    return {
        s: DIVISION_GROUP.get(d, "other grey matter")
        for s, d in table["division"].items()
    }


def coverage_counts(density: pd.DataFrame, min_coverage: float) -> tuple[int, int, int]:
    """Structures of the fit, those with a measured density, and those the rule asks."""
    fit = density[density["in_fit"]]
    n_fit = len(fit)
    return n_fit, int(fit["measured"].sum()), int(np.ceil(min_coverage * n_fit))


def density_rows_panel(ax: plt.Axes, variants: pd.DataFrame) -> None:
    """C of figure 14: per density measure, what it predicts and what is left with it.

    What it predicts alone and what the main model leaves with it, on the structures PSD95
    covers.
    """
    table = variants.set_index("key")
    rows = [(key, label) for key, label in DENSITY_ROWS if key in table.index]
    y = np.arange(len(rows))
    alone = np.array([table.loc[key, "density_alone"] for key, _ in rows])
    left = np.array([table.loc[key, "left"] for key, _ in rows])
    ax.barh(
        y - 0.18, alone, height=0.34, color=DENSITY_BLUE, label="density alone predicts"
    )
    ax.barh(y + 0.18, left, height=0.34, color=RED, label="the main model leaves")
    for k in range(len(rows)):
        ax.text(alone[k] + 0.01, y[k] - 0.18, f"{alone[k]:.0%}", va="center", fontsize=7)
        ax.text(left[k] + 0.01, y[k] + 0.18, f"{left[k]:.0%}", va="center", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels([label for _, label in rows], fontsize=8)

    # room under the last row for the legend
    ax.set_ylim(len(rows) + 0.5, -0.6)
    ax.set_xlim(0, 0.7)
    ticks = np.arange(0, 0.71, 0.1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    tidy(ax)


def points_interval(contrast: tuple[float, tuple[float, float]]) -> str:
    """A difference of shares in points with its 95% interval: '+7 (-11 to +25)'."""
    point, (lo, hi) = contrast
    return f"{100 * point:+.0f} ({100 * lo:+.0f} to {100 * hi:+.0f})"


def compared(contrast: tuple[float, tuple[float, float]], more: str, less: str) -> str:
    """`more` or `less` when the interval of a difference leaves zero.

    'about as much' otherwise.
    """
    point, (lo, hi) = contrast
    if lo > 0:
        return more
    if hi < 0:
        return less
    return "about as much"


def plot_synaptome(
    density: pd.DataFrame,
    coverage: pd.DataFrame,
    markers: pd.Series,
    numbers: dict,
    min_coverage: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 14: whether the measured synapse density covers the fit and predicts better.

    Better, that is, than the mRNA it would replace.

    `density` and `coverage` are the tables of run_synaptome, `markers` the marker
    composite per structure of the fit, `numbers` those of adult.beyond_figures (its
    variants.csv, and the paired differences between the density rows).
    """
    table = density.set_index("structure")
    n_fit, n_measured, needed = coverage_counts(table, min_coverage)
    measured = table[table["in_fit"] & table["measured"]]
    variants = numbers["variants"]
    rows = variants.set_index("key")
    psd95, panel = rows.loc["psd95"], rows.loc["panel"]
    contrasts = numbers["density_contrasts"]
    rule = "not below" if n_measured >= needed else "below"
    alone = compared(contrasts["psd95_minus_panel_alone"], "better", "less well")
    if alone == "about as much":
        alone = "about as well"
    left = compared(contrasts["psd95_minus_panel_left"], "more", "less")
    fig = plt.figure(figsize=(16, 6.4))
    heading(
        fig,
        "synaptome",
        f"PSD95 puncta cover {n_measured / n_fit:.0%} of the fit, {rule} the "
        f"{min_coverage:.0%} the rule asks; where measured they predict the map {alone} "
        f"than the mRNA panel ({psd95['density_alone']:.0%} against "
        f"{panel['density_alone']:.0%}), and the model leaves {left} "
        f"({psd95['left']:.0%} against {panel['left']:.0%})",
    )
    ax = fig.add_axes([0.07, 0.17, 0.24, 0.57])
    coverage_panel(ax, coverage)
    panel_title(
        ax,
        "A",
        "Structures of the fit with a measured density",
        f"{n_measured} of {n_fit}; the rule asks for {needed}",
    )
    ax = fig.add_axes([0.39, 0.17, 0.22, 0.57])
    n, rho = rank_scatter(
        ax,
        measured["psd95"],
        markers,
        synaptome_groups(table),
        ("PSD95 puncta (rank)", "marker mRNA composite (rank)"),
    )
    ax.legend(handles=group_handles(), fontsize=7, frameon=False, loc="upper left")
    panel_title(
        ax, "B", "PSD95 puncta against the mRNA", f"rho {rho:+.2f} on {n} structures"
    )
    ax = fig.add_axes([0.76, 0.17, 0.21, 0.57])
    density_rows_panel(ax, variants)
    panel_title(
        ax,
        "C",
        "Each density measure in the model",
        f"on the {int(psd95['n_structures'])} structures PSD95 covers;\nPSD95 minus "
        "panel, in points: left "
        f"{points_interval(contrasts['psd95_minus_panel_left'])},\nalone "
        f"{points_interval(contrasts['psd95_minus_panel_alone'])}",
    )
    footer(
        fig,
        [
            "How to read: PSD95 punctum density from one adult mouse (Zhu et al. 2018, "
            "as shared by Hansen et al.), placed in the CCF structures; the mRNA panel "
            "is the main model's density.",
            "A difference in C comes with its 95% interval over the same subsamples of "
            "the structures (jackknife). In detail, with the two hemispheres and every "
            f"density's agreement: {figure_ref('synaptome_detail')}.",
        ],
    )
    return saved(fig, save)


def coverage_takeaway(
    n_fit: int, n_measured: int, needed: int, min_coverage: float
) -> str:
    """Figure 14s's line under its title: the coverage against the rule.

    The rule says whether PSD95 enters the model.
    """
    if n_measured >= needed:
        rule = "not below"
        verdict = "it replaces the mRNA density terms in the main model"
    else:
        rule = "below"
        verdict = (
            "the mRNA terms stay in the main model, and PSD95 is a variant on the "
            f"{n_measured} structures it covers"
        )
    return (
        f"PSD95 puncta cover {n_measured / n_fit:.0%} of the fit, {rule} the "
        f"{min_coverage:.0%} the rule asks: {verdict}"
    )


def plot_synaptome_detail(
    density: pd.DataFrame,
    coverage: pd.DataFrame,
    agreement: pd.DataFrame,
    markers: pd.Series,
    min_coverage: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 14s: the measured synapse density, its coverage of the fit, how it compares.

    A: per division, the structures of analysis 4's fit with a measured PSD95 density.
    B: that density against the marker mRNA composite it would replace, as ranks.
    C: the left hemisphere against the right, the one check of a one-mouse map.
    D: each density's Spearman with the mRNA density terms, Gria1, the nano map and
    autofluorescence. `density`, `coverage` and `agreement` are the tables of
    run_synaptome; `markers` is the marker composite per structure of the fit.
    """
    table = density.set_index("structure")
    n_fit, n_measured, needed = coverage_counts(table, min_coverage)
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
        f"{n_measured} of {n_fit} ({n_measured / n_fit:.0%}); the rule asks for {needed} "
        f"({min_coverage:.0%})",
    )

    ax = fig.add_subplot(grid[0, 1])
    n, rho = rank_scatter(
        ax,
        measured["psd95"],
        markers,
        groups,
        ("PSD95 puncta (rank)", "marker mRNA composite (rank)"),
    )
    ax.legend(handles=group_handles(), fontsize=7, frameon=False, loc="upper left")
    panel_title(
        ax,
        "B",
        "PSD95 puncta against the mRNA they would replace",
        f"rho {rho:+.2f}, n = {n}",
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
        coverage_takeaway(n_fit, n_measured, needed, min_coverage),
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
        f"the main model, bent: {quoted / explainable:.1%} of the ceiling, "
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
    cubic: float,
    quintic: float,
    passed: dict[str, bool],
    save: Path | None = None,
) -> plt.Figure:
    """fig5: controls F and E, the gene-space curve and bending further.

    `curve` is gene_space.csv, `k_most` the number of components picked most often,
    `cubic` and `quintic` the held-out R2 bent to cubes and to fifth powers.
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
        [0, 1], [cubic, quintic], color=[RED, "0.65"], edgecolor="0.25", linewidth=0.5
    )
    axes[1].set_xticks([0, 1])
    axes[1].set_xticklabels(
        ["the model\n(squares, cubes)", "bent further\n(to fifth powers)"], fontsize=8
    )
    axes[1].set_ylabel("cross-validated R2", fontsize=8)
    if passed["E"]:
        title = "E. and bending it further buys nothing"
    else:
        title = "E. and bending it further still pays"
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

# the steps of the variance budget: abundance in the subunits' dark blue, density in
# the pale blue of DENSITY_BLUE, autofluorescence in its channel's yellow, what is left
# in red; control F's components in mid grey
BUDGET_COLOURS = {
    "abundance": DARK_BLUE,
    "density": DENSITY_BLUE,
    "autofluorescence": AUTO,
    "components": MID_GREY,
    "left": RED,
}

# the variants of figure 03s1 G by their key in variants.csv, short; the number of
# structures follows a row on other structures than the main model's
VARIANT_LABELS = {
    "main": "the main model",
    "one_shuffling": "one shuffling of the folds",
    "single_shufflings": "single shufflings (95%)",
    "ten_folds": "ten folds",
    "leave_one_out": "leave one out",
    "spatial_blocks": "folds of spatial blocks",
    "psd95": "density: PSD95 puncta",
    "panel": "density: mRNA panel",
    "panel_all": "density: mRNA panel",
    "psd95_and_panel": "density: PSD95 and panel",
    "panel_postsynaptic": "panel, postsynaptic markers",
    "sap102": "density: SAP102 puncta",
    "all_puncta": "density: every punctum",
    "four_subunits": "abundance: Gria1 to Gria4",
    "markers_once": "markers measured once out",
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


def budget_bar(
    ax: plt.Axes,
    y: float,
    parts: list[tuple[str, float, float, str]],
    height: float = 0.56,
) -> None:
    """One bar of segments from 0 to 1: (colour key, start, end, label inside)."""
    for key, a, b, label in parts:
        w = b - a
        if w <= 0:
            continue
        ax.barh(
            y,
            w,
            left=a,
            height=height,
            color=BUDGET_COLOURS[key],
            edgecolor="white",
            linewidth=1.5,
            zorder=2,
        )
        if w > 0.05 and label:
            ink = "white" if key in ("abundance", "left") else "0.1"
            narrow = w < 0.1
            ax.text(
                a + w / 2,
                y,
                label.replace("+ ", "+") if narrow else label,
                ha="center",
                va="center",
                fontsize=6.5 if narrow else 8,
                color=ink,
                fontweight="normal" if narrow else "bold",
                zorder=4,
            )


def hatched(ax: plt.Axes, y: float, start: float, width: float) -> None:
    """A hatched stretch over a bar: the calibration floor inside a leftover."""
    ax.barh(
        y,
        width,
        left=start,
        height=0.56,
        fill=False,
        hatch="////",
        edgecolor="white",
        linewidth=0,
        zorder=3,
    )


def model_parts(
    steps: list[float], abundance: str, left_label: bool
) -> list[tuple[str, float, float, str]]:
    """The segments of one model's budget bar.

    Abundance, + density, + autofluorescence (unlabelled, it adds about nothing), left
    (labelled when `left_label`).
    """
    end = steps[2]
    left = f"left\n{1 - end:.0%}" if left_label else ""
    return [
        ("abundance", 0, steps[0], f"{abundance}\n{steps[0]:.0%}"),
        (
            "density",
            steps[0],
            min(steps[1], end),
            f"+ density\n{steps[1] - steps[0]:.0%}",
        ),
        ("autofluorescence", steps[1], end, ""),
        ("left", end, 1, left),
    ]


def budget_rows(n: dict) -> list[tuple[int, str, list]]:
    """The rows of figure 03s1 D, top to bottom: (y, label, segments).

    Those of the main model, of its four-subunit check row and of control F.
    """
    four = n["variants"].set_index("key").loc["four_subunits"]
    four_steps = [
        four["abundance"],
        four["abundance"] + four["plus_density"],
        1 - four["left"],
    ]
    f_share = n["f_share"]
    return [
        (
            2,
            f"the main model\n({n['abundance']}, synapse density,\nautofluorescence)",
            model_parts(n["steps"], n["abundance"], left_label=False),
        ),
        (
            1,
            "check row: Gria1 to Gria4\nin place of Gria1, four terms",
            model_parts(four_steps, "Gria1-4", left_label=True),
        ),
        (
            0,
            f"control F: components of the {n['f_genes']} genes\nmeasured in every "
            f"structure (most often {n['f_k']})",
            [
                ("components", 0, f_share, f"{f_share:.0%}"),
                ("left", f_share, 1, f"left\n{1 - f_share:.0%}"),
            ],
        ),
    ]


def leftover_marks(ax: plt.Axes, n: dict) -> None:
    """The marks on the main model's bar of figure 03s1 D.

    Its leftover with share and interval, where 'left' begins over the jackknife, the
    benchmark's mark, and what autofluorescence adds.
    """
    steps = n["steps"]
    end_model = steps[2]
    floor = n["floor"]
    left = 1 - end_model
    lo, hi = n["left_ci"]

    # a bracket over the leftover with its share and what the hatch is
    ax.plot([end_model, 1], [2.36, 2.36], color="0.1", lw=1.2, zorder=5)
    for x in (end_model, 1):
        ax.plot([x, x], [2.3, 2.36], color="0.1", lw=1.2, zorder=5)
    ax.text(
        1.0,
        2.4,
        f"left {left:.0%} (95% {lo:.0%} to {hi:.0%}, jackknife over structures);\n"
        f"hatched: the floor, {floor['left_median']:.0%} on the calibration's "
        f"{n['cal_n']} structures (nano there {n['nano_cal_left']:.0%})",
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="0.1",
        zorder=6,
    )

    # under the bar, the interval of where the leftover begins
    ax.plot([1 - hi, 1 - lo], [1.64, 1.64], color="0.1", lw=1.4, zorder=5)
    for x in (1 - hi, 1 - lo):
        ax.plot([x, x], [1.6, 1.68], color="0.1", lw=1.0, zorder=5)
    ax.text(
        1 - lo + 0.005,
        1.64,
        "where 'left' begins, 95%",
        fontsize=7,
        va="center",
        ha="left",
        color="0.2",
    )

    # where the leftover of a map that is one Allen Gria1 experiment would begin
    mark = 1 - n["gria1"]["left_median"]
    ax.plot([mark, mark], [1.5, 1.72], color="0.1", lw=1.4, zorder=6)
    ax.plot([mark, mark], [1.72, 2.28], color="white", lw=1.6, ls=(0, (2, 1)), zorder=6)
    ax.text(
        mark - 0.006,
        1.45,
        f"a map that is one Allen Gria1 experiment leaves "
        f"{n['gria1']['left_median']:.0%}: its leftover would begin here",
        fontsize=7,
        color="0.1",
        ha="right",
        va="center",
        zorder=6,
    )
    ax.annotate(
        f"+ autofluorescence {end_model - steps[1]:+.1%} (alone "
        f"{n['partition']['autofluorescence']:+.0%})",
        xy=(min(steps[1], end_model), 2.28),
        xytext=(0.4, 2.5),
        fontsize=7.4,
        color="0.2",
        ha="right",
        va="bottom",
        arrowprops=dict(arrowstyle="-", color=NOTE_GREY, lw=0.7),
        zorder=6,
    )


def budget_panel(ax: plt.Axes, n: dict) -> None:
    """D: the variance budget of the main model, its four-subunit variant and control F.

    The main model's leftover carries a bracket with its share and jackknife
    interval, the calibration floor hatched at its start, and a dark blue mark where
    the leftover would begin if the map left as much as a map of one Allen Gria1
    experiment does; control F's leftover carries its own floor.
    """
    rows = budget_rows(n)
    for y, _, parts in rows:
        budget_bar(ax, y, parts)

    # the floors hatched at the start of the leftovers (the model's, control F's own)
    hatched(ax, 2, n["steps"][2], n["floor"]["left_median"])
    f_floor = float(np.mean(n["f_floor"]))
    hatched(ax, 0, n["f_share"], f_floor)
    leftover_marks(ax, n)
    ax.text(
        n["f_share"] + f_floor,
        -0.42,
        f"hatched: control F's own floor, {n['f_floor'][0]:.0%} to {n['f_floor'][1]:.0%}",
        fontsize=7,
        ha="right",
        va="top",
        color="0.2",
    )
    ax.set_yticks([r[0] for r in rows])
    ax.set_yticklabels([r[1] for r in rows], fontsize=8.5)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, 3.05)
    ticks = np.linspace(0, 1, 6)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel(
        f"share of the map's reproducible variance (the ceiling: {n['ceiling']:.1%} of "
        "the total, from the agreement of two halves of the cohort)"
    )
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)


def calibration_panel(ax: plt.Axes, calibration: pd.DataFrame, n: dict) -> None:
    """E: the share left against the leftover's replication, known maps and nano.

    Random folds only; the Gria1 map from each half of the experiments apart.
    """
    calibration = calibration[calibration["folds"] == "random"]
    kinds = (
        (
            "abundance and density",
            "",
            "s",
            DARK_GREY,
            DARK_GREY,
            "known answer: Gria1 + synapse density (the floor)",
        ),
        (
            "Gria1 mRNA",
            "A",
            "^",
            DARK_BLUE,
            DARK_BLUE,
            "known answer: one Gria1 experiment, from half A",
        ),
        (
            "Gria1 mRNA",
            "B",
            "^",
            "white",
            DARK_BLUE,
            "known answer: one Gria1 experiment, from half B",
        ),
    )
    for kind, half, marker, face, edge, label in kinds:
        mine = calibration[calibration["map"] == kind]
        if half:
            mine = mine[mine["truth_from"] == half]
        ax.scatter(
            mine["left"],
            mine["replication"],
            s=20,
            marker=marker,
            facecolors=face,
            edgecolors=edge,
            alpha=0.7,
            linewidths=0.8,
            zorder=2,
            label=f"{label} ({len(mine)} draws, median {mine['left'].median():.0%})",
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
    left = 1 - n["steps"][2]
    lo, hi = n["left_ci"]
    ax.errorbar(
        [left],
        [n["rep_left"]],
        xerr=[[left - lo], [hi - left]],
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
    ax.set_xlim(0, max(0.5, calibration["left"].max() + 0.04))
    ax.set_ylim(min(0.6, calibration["replication"].min() - 0.03), 1.0)
    ticks = np.arange(0, ax.get_xlim()[1] + 1e-9, 0.1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map left")
    ax.set_ylabel("the leftover's replication (half against half)")
    ax.legend(loc="lower right", fontsize=6.8, frameon=False)
    tidy(ax)


def replication_panel(ax: plt.Axes, replication: pd.DataFrame, n: dict) -> None:
    """F: the two half-maps' and the two leftovers' agreement over every split."""
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


def variants_panel(ax: plt.Axes, n: dict) -> None:
    """G: the leftover under other folds, models and structures, the main one red.

    A row on structures other than the main model's names their number; a dotted
    line separates the kinds of variant.
    """
    table = n["variants"].reset_index(drop=True)
    y = np.arange(len(table))
    for i in range(1, len(table)):
        if table.loc[i, "kind"] != table.loc[i - 1, "kind"]:
            ax.axhline(i - 0.5, color=MID_GREY, lw=0.6, ls=(0, (1, 2)), zorder=0)
    for i, r in table.iterrows():
        colour = RED if r["kind"] == "main" else "0.15"
        if np.isfinite(r["lo"]):
            ax.plot([r["lo"], r["hi"]], [i, i], color=MID_GREY, lw=2, zorder=1)
        ax.scatter(r["left"], i, s=30, color=colour, zorder=2)
        right = r["hi"] if np.isfinite(r["hi"]) else r["left"]
        ax.text(
            right + 0.025,
            i,
            f"{r['left']:.0%}",
            fontsize=7,
            ha="left",
            va="center",
            color=colour,
        )
    labels = []
    for _, r in table.iterrows():
        label = VARIANT_LABELS.get(r["key"], r["variant"])
        if r["n_structures"] != table["n_structures"].iloc[0]:
            label += f", {int(r['n_structures'])}"
        labels.append(label)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_ylim(len(table) - 0.4, -0.8)
    ax.set_xlim(0, 0.6)
    ticks = np.arange(0, 0.61, 0.2)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map left")
    tidy(ax)


def main_budget_panel(ax: plt.Axes, n: dict) -> None:
    """B of figure 03: the main model's budget as one bar, and where 'left' begins.

    Autofluorescence adds nothing, so it is said under the bar rather than drawn.
    """
    steps = n["steps"]
    left = 1 - steps[2]
    lo, hi = n["left_ci"]
    parts = [
        ("abundance", 0, steps[0], f"{n['abundance']}\n{steps[0]:.0%}"),
        (
            "density",
            steps[0],
            min(steps[1], steps[2]),
            f"+ synapse density\n{steps[1] - steps[0]:.0%}",
        ),
        ("autofluorescence", steps[1], steps[2], ""),
        ("left", steps[2], 1, f"left\n{left:.0%}"),
    ]
    budget_bar(ax, 0, parts, height=0.6)

    # the 95% interval of where 'left' begins, under the bar
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
    ax.text(
        0,
        0.36,
        f"autofluorescence adds {steps[2] - steps[1]:+.1%}",
        fontsize=7.5,
        va="bottom",
        color=DARK_GREY,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.7, 0.5)
    ticks = np.linspace(0, 1, 6)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_yticks([])
    ax.set_xlabel(
        "share of the map's reproducible pattern, on structures the fit has not seen"
    )
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)


def benchmarks_panel(ax: plt.Axes, calibration: pd.DataFrame, n: dict) -> None:
    """C of figure 03: what the main model leaves of nano, beside maps of known answer.

    On the calibration's structures.

    A dot per draw of made-up adults, the median as a bar; nano is the mean over the
    two halves of the Allen experiments its predictors come from.
    """
    random = calibration[calibration["folds"] == "random"]
    left = 1 - n["steps"][2]
    nano_label = (
        f"the nano map\n({n['nano_cal_left']:.0%} on these {n['cal_n']}; "
        f"{left:.0%} on all {n['n_structures']})"
    )
    rows = (
        ("abundance and density", "Gria1 + synapse density only\n(the floor)", DARK_GREY),
        ("nano", nano_label, RED),
        ("Gria1 mRNA", "one Allen Gria1 experiment\n(the benchmark)", DARK_BLUE),
    )
    rng = np.random.default_rng(0)
    for y, (kind, _, colour) in enumerate(rows):
        if kind == "nano":
            lo, hi = n["nano_cal_ci"]
            ax.errorbar(
                [n["nano_cal_left"]],
                [y],
                xerr=[[n["nano_cal_left"] - lo], [hi - n["nano_cal_left"]]],
                fmt="o",
                ms=9,
                color=colour,
                elinewidth=1.4,
                capsize=3,
                zorder=3,
            )
            ax.text(
                n["nano_cal_left"],
                y - 0.3,
                f"{n['nano_cal_left']:.0%} (95% {lo:.0%} to {hi:.0%})",
                ha="center",
                fontsize=8,
                color=colour,
            )
            continue
        left = random.loc[random["map"] == kind, "left"].to_numpy(float)
        jitter = rng.uniform(-0.13, 0.13, len(left))
        ax.scatter(left, y + jitter, s=14, color=colour, alpha=0.6, linewidths=0)
        median = float(np.median(left))
        ax.plot([median, median], [y - 0.28, y + 0.28], color="0.1", lw=2, zorder=3)
        ax.text(median, y - 0.3, f"{median:.0%}", ha="center", fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([label for _, label, _ in rows], fontsize=8.5)
    ax.set_ylim(len(rows) - 0.5, -0.7)
    right = max(float(random["left"].max()), n["nano_cal_ci"][1])
    ax.set_xlim(0, max(0.5, right + 0.04))
    ticks = np.arange(0, ax.get_xlim()[1] + 1e-9, 0.1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map left by the same model")
    tidy(ax)


def beyond_takeaway(n: dict) -> str:
    """The line under figure 03's title: what is left, and how it stands.

    Against the floor and the benchmark, each by the interval of its difference from nano.
    """
    left = 1 - n["steps"][2]
    lo, hi = n["minus_floor_ci"]
    floor = n["floor"]["left_median"]
    margin = f"{n['minus_floor']:+.0%}, 95% {lo:+.0%} to {hi:+.0%}"
    if lo > 0:
        against = f"above the floor's {floor:.0%} ({margin})"
    else:
        against = f"the floor {floor:.0%} ({margin})"
    return (
        f"Gria1 and synapse density leave {left:.0%} of the map's reproducible pattern, "
        f"which replicates ({n['rep_left']:.2f}); on {n['cal_n']} structures nano "
        f"leaves {n['nano_cal_left']:.0%}, {against}, "
        f"one Allen Gria1 experiment {n['gria1']['left_median']:.0%}"
    )


def plot_beyond(
    y: np.ndarray,
    held_out: dict[str, np.ndarray],
    residual: np.ndarray,
    groups: list[str],
    acronyms: list[str],
    numbers: dict,
    calibration: pd.DataFrame,
    replication: pd.DataFrame,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 03: what Gria1 expression and synapse density predict of the map.

    What they leave is set against the floor and the benchmark.

    The arguments are those of plot_beyond_budget, its detailed version, but Gria1.
    """
    n = numbers
    fig = plt.figure(figsize=(15, 10.5))
    heading(fig, "beyond", beyond_takeaway(n))

    # A: the map against the main model's held-out prediction
    ax = fig.add_axes([0.07, 0.53, 0.3, 0.36])
    model_scatter(
        ax,
        y,
        held_out["model"],
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
        f"{n['n_structures']} structures; off the diagonal: the leftover",
    )

    # B: the budget
    ax = fig.add_axes([0.45, 0.6, 0.5, 0.2])
    main_budget_panel(ax, n)
    panel_title(
        ax,
        "B",
        "What each term predicts, and what is left",
        f"the main model: {n['abundance']} + synapse density + autofluorescence, each "
        "bent (x, x², x³);\n"
        + textwrap.fill(
            f"synapse density is the mRNA of {n['density_label']} (PSD95 puncta cover "
            f"too few structures: {figure_ref('synaptome')})",
            100,
        ),
    )

    # C: the leftover beside the floor and the benchmark
    ax = fig.add_axes([0.19, 0.09, 0.3, 0.3])
    benchmarks_panel(ax, calibration, n)
    lo, hi = n["minus_floor_ci"]
    glo, ghi = n["minus_gria1_ci"]
    panel_title(
        ax,
        "C",
        "Is the leftover more than Allen-to-Allen mismatch?",
        f"{n['cal_n']} structures; nano minus the floor {n['minus_floor']:+.0%} "
        f"(95% {lo:+.0%} to {hi:+.0%}),\nminus the benchmark {n['minus_gria1']:+.0%} "
        f"({glo:+.0%} to {ghi:+.0%})",
    )

    # D: the replication
    ax = fig.add_axes([0.62, 0.09, 0.34, 0.3])
    replication_panel(ax, replication, n)
    panel_title(
        ax,
        "D",
        "Does the leftover replicate across mice?",
        f"observed {n['rep_left']:.3f}; {n['implied']:.3f} expected from the map's "
        "reliability alone;\nunrelated leftovers "
        f"{n['noise'][0]:+.2f} to {n['noise'][1]:+.2f} (off scale)",
    )
    footer(
        fig,
        [
            "How to read: the map and every predictor are ranks across structures; a "
            "share is the R² on structures the fit has not seen, over what two halves "
            "of the cohort reproduce.",
            "The floor and the benchmark are maps whose answer is known, made from half "
            "of each gene's Allen experiments with ten made-up adults as noisy as ours. "
            f"In detail, with the check rows: {figure_ref('beyond_budget')}; the seven "
            f"controls: {figure_ref('beyond_controls')}.",
        ],
    )
    return saved(fig, save)


def predictor_scatters(
    fig: plt.Figure,
    y: np.ndarray,
    gria1: np.ndarray,
    held_out: dict[str, np.ndarray],
    residual: np.ndarray,
    groups: list[str],
    acronyms: list[str],
    n: dict,
) -> None:
    """A to C of figure 03s1: the map against Gria1, synapse density and the model.

    Against Gria1, what synapse density predicts and what the whole model predicts, the
    latter two held out.
    """
    top = 0.705
    size = 0.2
    rho_g = spearmanr(y, gria1).statistic
    ax = fig.add_axes([0.06, top, size, size * 16 / 19])
    map_scatter(ax, gria1, y, groups, "Gria1 mRNA, rank")
    panel_title(
        ax,
        "A",
        "The nano map against Gria1 mRNA",
        f"rho {rho_g:+.2f} on the fit's {n['n_structures']} structures; Gria1 "
        f"alone,\nbent, predicts {n['partition']['abundance']:.0%} of the reproducible "
        "map",
    )
    ax = fig.add_axes([0.39, top, size, size * 16 / 19])
    map_scatter(
        ax,
        held_out["density"],
        y,
        groups,
        "what synapse density predicts (held out), rank",
    )
    panel_title(
        ax,
        "B",
        "The nano map against synapse density",
        f"{textwrap.fill(n['density_label'], 60)},\nbent: "
        f"{n['partition']['density']:.0%} of the reproducible map",
    )
    ax = fig.add_axes([0.72, top, size, size * 16 / 19])
    model_scatter(
        ax,
        y,
        held_out["model"],
        residual,
        groups,
        acronyms,
        "what the model predicts (held out), rank",
    )
    panel_title(
        ax,
        "C",
        "The nano map against the whole model",
        f"{n['abundance']}, synapse density and autofluorescence, bent: "
        f"{n['steps'][2]:.0%} of\nthe reproducible map; the spread off the diagonal is "
        "the leftover",
    )


def budget_row(fig: plt.Figure, n: dict) -> None:
    """D and G of figure 03s1: the variance budget, and the leftover under other choices.

    Every other choice of folds, model and structures.
    """
    ax = fig.add_axes([0.22, 0.44, 0.5, 0.17])
    budget_panel(ax, n)
    fig.text(
        0.06,
        0.635,
        f"D.  The variance budget, on structures the fit has not seen: what "
        f"{n['abundance']} predicts, what synapse density adds, what is left",
        fontsize=9,
        va="bottom",
    )
    ax = fig.add_axes([0.85, 0.44, 0.12, 0.18])
    variants_panel(ax, n)
    blocks = n["blocks"]
    panel_title(
        ax,
        "G",
        "Other folds, models, structures",
        f"spatial blocks on the {n['cal_n']}\nstructures: nano {blocks['nano']:.0%}, "
        f"floor {blocks['floor']:.0%},\nGria1 map {blocks['gria1']:.0%}",
    )


def calibration_row(
    fig: plt.Figure, calibration: pd.DataFrame, replication: pd.DataFrame, n: dict
) -> None:
    """E and F of figure 03s1: the calibration, and the leftover's replication.

    The same model on maps whose answer is known, and the leftover of one half-cohort
    against the other's.
    """
    ax = fig.add_axes([0.07, 0.1, 0.38, 0.25])
    calibration_panel(ax, calibration, n)
    floor = n["floor"]
    panel_title(
        ax,
        "E",
        "The same model on maps whose answer is known",
        f"on {n['cal_n']} structures: the floor {floor['left_median']:.0%} "
        f"({floor['left_lo']:.0%} to {floor['left_hi']:.0%}); one Gria1 experiment "
        f"{n['gria1_halves']['A']:.0%} from half A, {n['gria1_halves']['B']:.0%} from "
        f"half B;\nnano {n['nano_cal_left']:.0%}: minus the floor "
        f"{n['minus_floor']:+.0%} ({n['minus_floor_ci'][0]:+.0%} to "
        f"{n['minus_floor_ci'][1]:+.0%}), minus the Gria1 map {n['minus_gria1']:+.0%} "
        f"({n['minus_gria1_ci'][0]:+.0%} to {n['minus_gria1_ci'][1]:+.0%})",
    )
    ax = fig.add_axes([0.6, 0.1, 0.36, 0.25])
    replication_panel(ax, replication, n)
    panel_title(
        ax,
        "F",
        "Does the leftover replicate across mice?",
        f"{len(replication)} ways to split the {n['n_adults']} adults into two fives; "
        f"two unrelated leftovers\nwould agree between {n['noise'][0]:+.2f} and "
        f"{n['noise'][1]:+.2f}, far left of this axis",
    )


def plot_beyond_budget(
    y: np.ndarray,
    gria1: np.ndarray,
    held_out: dict[str, np.ndarray],
    residual: np.ndarray,
    groups: list[str],
    acronyms: list[str],
    numbers: dict,
    calibration: pd.DataFrame,
    replication: pd.DataFrame,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 03s1: how much of the map Gria1 and synapse density predict, in detail.

    `y` is the map's ranks on the structures of the fit, `gria1` Gria1's, `held_out`
    each structure's prediction by the density model and by the whole model from
    fits that never saw it, `residual` the model's leftover; `numbers` holds the
    numbers of adult.beyond_figures, `calibration` and `replication` the tables of
    beyond_calibration and beyond_density.
    """
    n = numbers
    left = 1 - n["steps"][2]
    fig = plt.figure(figsize=(16, 19))
    heading(
        fig,
        "beyond_budget",
        f"{n['n_structures']} grey-matter structures, {n['n_adults']} adults: "
        f"{left:.0%} of the reproducible map left ({n['left_ci'][0]:.0%} to "
        f"{n['left_ci'][1]:.0%}); on the {n['cal_n']} structures of the calibration "
        f"nano leaves {n['nano_cal_left']:.0%}, the floor {n['floor']['left_median']:.0%}"
        f" (difference {n['minus_floor']:+.0%}, {n['minus_floor_ci'][0]:+.0%} to "
        f"{n['minus_floor_ci'][1]:+.0%}), one Gria1 experiment "
        f"{n['gria1']['left_median']:.0%}",
    )
    fig.legend(
        handles=group_handles(),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.935),
        ncol=4,
        fontsize=8,
        frameon=False,
        title="dots of A to C: structures, by group of divisions",
        title_fontsize=8,
    )
    predictor_scatters(fig, y, gria1, held_out, residual, groups, acronyms, n)
    budget_row(fig, n)
    calibration_row(fig, calibration, replication, n)
    footer(
        fig,
        [
            "How to read: the map and every predictor are ranks across structures; each "
            "predictor enters as x, x² and x³, and each share is the held-out R² "
            "(20 shufflings of five folds) over the ceiling.",
            "The calibration builds maps whose answer is known from half of each "
            "gene's Allen experiments, gives them ten made-up adults as noisy as ours, "
            "and predicts them from the other half; the differences in E are taken on "
            "the same structures,",
            "resampled together. Control F: the components of the genes measured "
            "everywhere predict most of the map, though no single gene or set follows "
            f"the leftover ({figure_ref('leftover_genes')}).",
            "What would mean what: a leftover well above the floor: part of the map is "
            "set by something Gria1 expression and synapse density do not predict; no "
            "larger than a one-gene map's: its size is within",
            "what one Allen experiment disagreeing with another produces. Either way it "
            "is what the model does not predict, not a measurement of the surface "
            "fraction, and its replication follows from the map's reliability (F).",
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
    "E": "Curvature the model misses?",
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
    """E: the held-out R2 of the model straight, bent to cubes and to fifth powers."""
    labels = ["straight", "x, x², x³\n(the model)", "to x⁵"]
    colours = [MID_GREY, DARK_GREY, MID_GREY]
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
            "How to read: the leftover is the nano rank minus the rank Gria1, synapse "
            "density and autofluorescence predict (in-sample fit), on the "
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
            "How to read: the leftover is what the main model (Gria1 + synapse density "
            "+ autofluorescence) leaves, so it carries nothing of the model's columns; "
            "each surrogate goes through the same fit before it is correlated (a",
            "Freedman-Lane null), so a gene sharing the model's pattern meets a null of "
            "the right width. Gria1 and the markers enter the model and sit near zero by "
            "construction; the PSD genes enter only through their first component.",
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
