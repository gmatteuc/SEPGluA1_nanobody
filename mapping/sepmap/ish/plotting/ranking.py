"""The figures of analysis 1, each gene against the map, and of its two checks.

Figure 05 shows what one comparison is (a gene's ranks against the map's, structure
by structure) and 05s the same step by step; 06 and 06s the spatial null, the rho
that unrelated maps as smooth as the map give; 07s P9's genes against the map and
the Cacng8 - Gria1 gap. The two checks of the ranking: 12 and 12s, whether it is the
label's or the tissue's (each gene against the autofluorescence map), and 13 and
13s, whether it depends on the choices made.

Called by run_ish_gene_ranking.py (05, 06, 07s, 12) and run_ish_robustness.py (13).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle
from scipy.stats import spearmanr

from sepmap.ish import spatial_null
from sepmap.ish.figure_index import figure_ref
from sepmap.ish.gene_ranking import QUOTED_GENES
from sepmap.ish.numbers import ordinal
from sepmap.ish.plotting.shared import (
    RANK_FLOOR,
    acronym_of,
    colour_bar,
    footer,
    gene_colour,
    group_handles,
    group_of,
    heading,
    overlay,
    p_text,
    panel_title,
    rank_plane,
    ranks01,
    saved,
    scatter_groups,
    spread_labels,
)
from sepmap.ish.spatial_null import ALPHA
from sepmap.plotting import (
    AUTO,
    AUTO_DOT,
    DARK_BLUE,
    DARK_GREY,
    LIGHT_GREY,
    MID_GREY,
    NANO,
    NANO_DOT,
    NO_DATA_GREY,
    NULL_BAND,
    PAIR_LINE,
    RED,
    bars_grey,
    draw_plane,
    hot_cut,
    paint,
    tidy,
)

# ===== 05 What one comparison is =====

# the structures named on the scatters of figure 05, by acronym: the top of the map,
# a thalamic relay nucleus, the striatum, the barrel field
NAMED_STRUCTURES = ("CA1", "VPM", "CP", "SSp-bfd")


def comparison_scatter(
    ax: plt.Axes,
    row: dict,
    map_values: pd.Series,
    groups: dict[str, str],
    acronyms: dict[str, str],
    n_surrogates: int,
    letter: str,
) -> None:
    """One gene's scatter of ranks: the map's against the gene's, dots by group."""
    shared = [s for s in map_values.index if s in row["profile"]]
    x = ranks01(map_values[shared].to_numpy())
    y = ranks01(np.array([row["profile"][s] for s in shared]))
    scatter_groups(ax, x, y, [groups[s] for s in shared])
    ax.plot([0, 1], [0, 1], color=MID_GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
    by_acronym = {acronyms[s]: i for i, s in enumerate(shared)}
    for acronym in NAMED_STRUCTURES:
        if acronym in by_acronym:
            i = by_acronym[acronym]
            # above the dot when it sits near the bottom, so no name lands on the axis
            ax.annotate(
                acronym,
                (x[i], y[i]),
                xytext=(4, 4) if y[i] < 0.1 else (4, -9),
                textcoords="offset points",
                fontsize=7,
                color="0.15",
            )
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_aspect("equal")
    ax.set_xlabel("nano map, rank among structures")
    ax.set_ylabel(f"{row['symbol']}, rank among structures")
    r = row["ranking"]
    panel_title(
        ax,
        letter,
        f"Spearman rho {r['rho']:+.2f} on {int(r['n_structures'])} structures",
        f"spatial {p_text(r['p_spatial'], n_surrogates)}; rank {int(r['rank_all'])} "
        f"of {row['n_all']} genes",
    )
    tidy(ax)


def steps_panel(ax: plt.Axes) -> None:
    """The three steps of one comparison, and under them the colours of the scatters."""
    ax.axis("off")
    # the key below the text, reaching into the gap above the next row
    ax.legend(
        handles=group_handles(),
        loc="upper left",
        bbox_to_anchor=(0.0, 0.15),
        title="dots of the scatters: structures by group of divisions",
        title_fontsize=8,
        fontsize=8,
        alignment="left",
    )
    steps = (
        "How one gene is compared with the map\n\n"
        "1  nano: each declared structure's mean in each adult,\n"
        "    as zref (log2 over the brain's isocortex mean,\n"
        "    centred and scaled over the declared set), then\n"
        "    the mean over the 10 adults\n"
        "2  gene: each Allen experiment's mean energy per\n"
        "    structure (200 um grid, failed sections missing),\n"
        "    as ranks, averaged over the gene's experiments\n"
        "3  keep the declared structures both have, rank each\n"
        "    side, correlate the ranks: Spearman rho\n\n"
        f"The spatial p ({figure_ref('spatial_null')}) is how often a map with the\n"
        "nano map's smoothness, and no relation to the gene,\n"
        "correlates with it as strongly."
    )
    ax.text(0, 1.0, steps, va="top", ha="left", fontsize=8.5, linespacing=1.3)


# the two genes of figure 05: one that follows the map, the astrocyte gene that does not
COMPARISON_MAIN = ("Cacng8", "Aqp4")


def plot_one_comparison(
    map_values: pd.Series,
    rows: list[dict],
    set_table: pd.DataFrame,
    lab: np.ndarray,
    names: dict[int, str],
    plane: int,
    n_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 05: a gene's rho with the map, for a gene that follows it and one not.

    `rows` are those of plot_one_comparison_detail; the genes of COMPARISON_MAIN
    are drawn.
    """
    by_symbol = {r["symbol"]: r for r in rows}
    follows = by_symbol[COMPARISON_MAIN[0]]
    unrelated = by_symbol[COMPARISON_MAIN[1]]
    groups = group_of(set_table)
    acronyms = acronym_of(set_table)
    fig = plt.figure(figsize=(17, 5.6))
    heading(
        fig,
        "one_comparison",
        "A gene's rho is how alike its map and the nano map order the structures: "
        f"{follows['symbol']} {follows['ranking']['rho']:+.2f}, "
        f"{unrelated['symbol']} {unrelated['ranking']['rho']:+.2f}",
    )

    # A and B: the nano map and the gene's, each as one rank per structure
    ax = fig.add_axes([0.01, 0.15, 0.22, 0.62])
    rank_plane(
        fig, ax, dict(zip(map_values.index, ranks01(map_values.to_numpy()))), lab, names
    )
    panel_title(ax, "A", "The nano map", f"one rank per structure, CCF plane {plane}")
    ax = fig.add_axes([0.26, 0.15, 0.22, 0.62])
    shared = [s for s in map_values.index if s in follows["profile"]]
    ranks = ranks01(np.array([follows["profile"][s] for s in shared]))
    rank_plane(fig, ax, dict(zip(shared, ranks)), lab, names)
    panel_title(ax, "B", f"{follows['symbol']} mRNA", "the same, from Allen ISH")

    # C and D: the two genes against the map
    for k, row in enumerate((follows, unrelated)):
        ax = fig.add_axes([0.56 + 0.23 * k, 0.15, 0.17, 0.62])
        comparison_scatter(ax, row, map_values, groups, acronyms, n_surrogates, "CD"[k])
    fig.legend(
        handles=group_handles(),
        loc="lower center",
        bbox_to_anchor=(0.78, 0.0),
        ncol=4,
        fontsize=7.5,
        frameon=False,
    )
    footer(
        fig,
        [
            "How to read: grey structures have no rank; in C and D one dot is one "
            "structure, coloured by group of divisions.",
            "Step by step, with the raw maps and Gria1: "
            f"{figure_ref('one_comparison_detail')}.",
        ],
    )
    return saved(fig, save)


def plot_one_comparison_detail(
    map_values: pd.Series,
    nano_image: np.ndarray,
    rows: list[dict],
    set_table: pd.DataFrame,
    lab: np.ndarray,
    names: dict[int, str],
    plane: int,
    n_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 05s: what one gene's rho with the map is, for nano and a few genes.

    `map_values` is the adult map on the declared structures, `nano_image` the
    cohort's cref on the plane; each of `rows` holds a gene's symbol, its ISH
    `image` on the plane and a `caption` naming the experiment shown, its merged
    `profile` ({structure: value}), its `ranking` row (gene_ranking.csv) and
    `n_all`, the genes ranked.
    """
    n_rows = 1 + len(rows)
    fig = plt.figure(figsize=(15, 4.3 * n_rows + 1.2))
    names_line = "; ".join(
        f"{r['symbol']} {r['ranking']['rho']:+.2f} "
        f"({p_text(r['ranking']['p_spatial'], n_surrogates)})"
        for r in rows
    )
    heading(
        fig,
        "one_comparison_detail",
        f"one rank correlation across the declared structures per gene: {names_line}",
    )
    grid = fig.add_gridspec(
        n_rows,
        3,
        width_ratios=[1, 1, 0.78],
        hspace=0.45,
        wspace=0.3,
        left=0.05,
        right=0.98,
        top=1 - 1.45 / fig.get_figheight(),
        bottom=1.25 / fig.get_figheight(),
    )
    groups = group_of(set_table)
    acronyms = acronym_of(set_table)

    # the nano row: the map as measured, as ranks, and the steps
    ax = fig.add_subplot(grid[0, 0])
    image = draw_plane(ax, nano_image, lab, hot_cut(), 0, 2.0)
    colour_bar(fig, ax, image, "nano / own isocortex mean")
    panel_title(ax, "A", f"Adult nano map, 10 adults, CCF plane {plane} (20 um)")
    ax.set_ylabel("nano", fontsize=11)
    ax = fig.add_subplot(grid[0, 1])
    rank_plane(
        fig, ax, dict(zip(map_values.index, ranks01(map_values.to_numpy()))), lab, names
    )
    panel_title(
        ax,
        "B",
        "One value per structure, as its rank",
        f"{len(map_values)} declared structures; flat grey: not declared",
    )
    steps_panel(fig.add_subplot(grid[0, 2]))

    # one row per gene: its ISH section, its ranks, the scatter with the map
    letters = "CDEFGHIJKL"
    for k, row in enumerate(rows, start=1):
        ax = fig.add_subplot(grid[k, 0])
        finite = row["image"][(lab > 0) & np.isfinite(row["image"])]
        vmax = float(np.percentile(finite, 99)) if finite.size else 1.0
        image = draw_plane(ax, row["image"], lab, hot_cut(), 0, vmax)
        colour_bar(fig, ax, image, "expression energy (Allen units)")
        panel_title(ax, letters[3 * (k - 1)], f"{row['symbol']} ISH, {row['caption']}")
        ax.set_ylabel(row["symbol"], fontsize=11)
        ax = fig.add_subplot(grid[k, 1])
        shared = [s for s in map_values.index if s in row["profile"]]
        ranks = dict(zip(shared, ranks01(np.array([row["profile"][s] for s in shared]))))
        rank_plane(fig, ax, ranks, lab, names)
        panel_title(
            ax,
            letters[3 * (k - 1) + 1],
            "Its merged profile, as ranks",
            f"{len(shared)} declared structures with ISH data",
        )
        ax = fig.add_subplot(grid[k, 2])
        comparison_scatter(
            ax, row, map_values, groups, acronyms, n_surrogates, letters[3 * k - 1]
        )
    footer(
        fig,
        [
            "How to read: each row reduces a map to one value per declared structure "
            "(middle) and the scatter puts the two sets of ranks against each other; "
            "rho is 1 when the orders agree and 0 when they are unrelated.",
            "Much of a whole-brain rho is cortex and hippocampus (orange, red) "
            f"against thalamus and the rest (blue, grey): {figure_ref('spatial_null')} "
            "asks how large a rho such a gradient gives by chance, "
            f"{figure_ref('between_within')} what is left inside divisions.",
        ],
    )
    return saved(fig, save)


# ===== 06 The spatial null =====


def null_rho_panel(ax: plt.Axes, calibration: pd.DataFrame) -> None:
    """A: rho between independent smooth maps, against rho between shuffled maps."""
    pairs = calibration[calibration["design"] == "pair"]
    bins = np.linspace(-1, 1, 41)
    ax.hist(pairs["rho"], bins=bins, color=DARK_GREY, label="two smooth maps")
    ax.hist(
        pairs["rho_shuffled"],
        bins=bins,
        histtype="step",
        color=MID_GREY,
        lw=1.5,
        label="one of them shuffled",
    )
    ax.set_xlabel("Spearman rho between two independent maps")
    ax.set_ylabel("pairs of maps")
    ax.legend(loc="upper left")
    panel_title(
        ax,
        "A",
        "Unrelated smooth maps correlate by chance",
        f"{len(pairs)} pairs of random fields; SD of rho {pairs['rho'].std():.2f}, "
        f"shuffled {pairs['rho_shuffled'].std():.2f}",
    )
    tidy(ax)


def variogram_panel(ax: plt.Axes, variogram: pd.DataFrame, map_name: str) -> None:
    """B: the map's variogram, its surrogates' (median and 5 to 95%), a shuffled map's."""
    sub = variogram[variogram["map"] == map_name]
    scale = sub["shuffled_median"].mean()
    h = sub["distance_mm"]
    ax.fill_between(
        h,
        sub["surrogate_p5"] / scale,
        sub["surrogate_p95"] / scale,
        color=NULL_BAND,
        lw=0,
        label="surrogates, 5 to 95%",
    )
    ax.plot(
        h, sub["surrogate_median"] / scale, color=DARK_GREY, lw=1.2, label="surrogates"
    )
    ax.plot(
        h, sub["shuffled_median"] / scale, color=MID_GREY, ls="--", lw=1, label="shuffled"
    )
    ax.plot(
        h,
        sub["variogram"] / scale,
        color=NANO,
        marker="o",
        ms=3,
        lw=1.4,
        label="nano map",
    )
    edge = sub.loc[sub["matched"], "distance_mm"].max()
    ax.axvline(edge, color=DARK_GREY, ls=":", lw=0.9)
    ax.text(edge, 0.05, " matched up to here", fontsize=7, color=DARK_GREY)
    ax.set_ylim(0, None)
    ax.set_xlabel("distance between structures (mm)")
    ax.set_ylabel("semivariance of ranks (shuffled map = 1)")
    ax.legend(loc="upper left", fontsize=7)
    panel_title(
        ax,
        "B",
        "The surrogates keep the map's smoothness",
        "near structures alike (low), far ones not; a shuffled map is flat",
    )
    tidy(ax)


def rates_panel(ax: plt.Axes, calibration: pd.DataFrame, letter: str = "D") -> None:
    """D: false positives at p < 0.05, ordinary against spatial p, for both designs."""
    labels = {"pair": "two random\nsmooth maps", "map": "a random map\nagainst nano"}
    designs = list(dict.fromkeys(calibration["design"]))
    for i, design in enumerate(designs):
        sub = calibration[calibration["design"] == design]
        for j, (column, colour) in enumerate(
            (("p_ordinary", LIGHT_GREY), ("p_spatial", DARK_GREY))
        ):
            rate = float((sub[column] < ALPHA).mean())
            x = i + (j - 0.5) * 0.38
            ax.bar(x, rate, width=0.36, color=colour)
            # above the 5% line when the bar is lower, so the two never overlap
            ax.text(x, max(rate, 0.05) + 0.012, f"{rate:.1%}", ha="center", fontsize=7.5)
    ax.axhline(0.05, color=RED, ls="--", lw=1)
    ax.text(-0.45, 0.062, "5% expected", color=RED, fontsize=7, ha="left")
    ax.set_xticks(range(len(designs)))
    ax.set_xticklabels([labels[d] for d in designs])
    ax.set_ylabel("share of tests with p < 0.05")
    ax.legend(
        handles=[
            Patch(color=LIGHT_GREY, label="ordinary Spearman p"),
            Patch(color=DARK_GREY, label="spatial p"),
        ],
        loc="center",
    )
    n_tests = int(calibration["design"].value_counts().iloc[0])
    panel_title(
        ax,
        letter,
        "False positives on maps with no relation",
        f"{n_tests} tests per design",
    )
    tidy(ax)


def surrogate_planes(
    fig: plt.Figure,
    axes: list[plt.Axes],
    map_ranks: dict[str, float],
    surrogate_ranks: list[dict[str, float]],
    lab: np.ndarray,
    names: dict[int, str],
) -> None:
    """C: the map and three surrogates on one plane, each as ranks."""
    panels = [("nano map", map_ranks)]
    for i, ranks in enumerate(surrogate_ranks):
        panels.append((f"surrogate {i + 1}", ranks))
    shown = {names.get(int(i)) for i in np.unique(lab) if i}
    undeclared = np.isin(
        lab, [i for i, n in names.items() if n in shown - set(map_ranks)]
    )
    for ax, (title, values) in zip(axes, panels):
        image = draw_plane(ax, paint(lab, values, names), lab, hot_cut(), RANK_FLOOR, 1)
        overlay(ax, undeclared, NO_DATA_GREY, lab)
        ax.set_title(title, fontsize=8)

    # one colour bar for the four, beside the last, in an axis of its own so the four
    # planes keep one size
    box = axes[len(panels) - 1].get_position()
    cax = fig.add_axes([box.x1 + 0.005, box.y0 + 0.03, 0.007, box.height - 0.06])
    cb = fig.colorbar(image, cax=cax)
    cb.ax.tick_params(labelsize=7)
    cb.set_label("rank (0 low, 1 high)", fontsize=8)
    cb.ax.set_ylim(0, 1)


def gene_null_panel(ax: plt.Axes, genes: list[tuple]) -> None:
    """E: each gene's null distribution of rho, its observed rho and p."""
    bins = np.linspace(-1, 1, 81)
    for name, observed, null, p in genes:
        colour = DARK_BLUE if name == "Gria1" else RED
        ax.hist(null, bins=bins, histtype="step", color=colour, lw=1.2)
    top = ax.get_ylim()[1]
    # the observed lines stop under the text, which sits above the distributions
    for name, observed, null, p in genes:
        colour = DARK_BLUE if name == "Gria1" else RED
        ax.vlines(observed, 0, top * 1.05, color=colour, lw=1.6)
    ax.set_ylim(0, top * 1.5)
    for k, (name, observed, null, p) in enumerate(genes):
        colour = DARK_BLUE if name == "Gria1" else RED
        lo, hi = np.percentile(null, [2.5, 97.5])
        ax.text(
            -0.97,
            top * (1.4 - 0.18 * k),
            f"{name}: rho {observed:+.2f}, spatial {p_text(p, len(null))}\n"
            f"    95% of its null {lo:+.2f} to {hi:+.2f}",
            color=colour,
            fontsize=7.5,
            va="center",
        )
    ax.set_xlabel("Spearman rho with a surrogate of the nano map")
    ax.set_ylabel("surrogates")
    panel_title(ax, "E", "Two genes against their null", "vertical lines: observed rho")
    tidy(ax)


def rates_by_design(calibration: pd.DataFrame) -> pd.DataFrame:
    """ish.spatial_null.false_positive_rates by design: ordinary and spatial."""
    return spatial_null.false_positive_rates(calibration).set_index("design")


def plot_spatial_null(
    variogram: pd.DataFrame,
    calibration: pd.DataFrame,
    n_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 06: why a spatial null, and that it works.

    `variogram` and `calibration` are the tables of run_ish_spatial_null.
    """
    pairs = calibration[calibration["design"] == "pair"]
    wider = pairs["rho"].std() / pairs["rho_shuffled"].std()
    rates = rates_by_design(calibration)
    ordinary, spatial = rates["ordinary"], rates["spatial"]
    fig = plt.figure(figsize=(16, 5.6))
    heading(
        fig,
        "spatial_null",
        f"Unrelated smooth maps correlate with an SD of rho {pairs['rho'].std():.2f}, "
        f"{wider:.0f} times a shuffled map's; against {n_surrogates} surrogates with the "
        f"map's smoothness, {spatial.get('map', np.nan):.1%} of random maps reach p < "
        f"0.05, where the ordinary p gives {ordinary.get('map', np.nan):.0%}",
    )
    null_rho_panel(fig.add_axes([0.05, 0.15, 0.25, 0.62]), calibration)
    variogram_panel(fig.add_axes([0.38, 0.15, 0.27, 0.62]), variogram, "nano")
    rates_panel(fig.add_axes([0.73, 0.15, 0.24, 0.62]), calibration, letter="C")
    footer(
        fig,
        [
            "How to read: a surrogate is the nano map's ranks shuffled, smoothed over "
            "near structures and rescaled until its variogram matches the map's (Burt "
            "2020); a gene's spatial p is the share of surrogates that correlate with "
            "it as strongly as the map does.",
            f"The surrogates on a plane and two genes against their null: "
            f"{figure_ref('spatial_null_detail')}.",
        ],
    )
    return saved(fig, save)


def plot_spatial_null_detail(
    variogram: pd.DataFrame,
    calibration: pd.DataFrame,
    map_ranks: dict[str, float],
    surrogate_ranks: list[dict[str, float]],
    lab: np.ndarray,
    names: dict[int, str],
    plane: int,
    n_surrogates: int,
    genes: list[tuple] | None = None,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 06s: why a spatial null, and the null itself.

    `map_ranks` and each of `surrogate_ranks` give a 0-1 rank per structure; `genes`,
    when given, holds (name, observed rho, null rhos, p) for panel E, which needs
    the gene ranking and is left out otherwise.
    """
    fig = plt.figure(figsize=(16, 10.5))
    spatial = rates_by_design(calibration)["spatial"]
    heading(
        fig,
        "spatial_null_detail",
        f"{n_surrogates} surrogates of the nano map on {len(map_ranks)} structures; "
        f"spatial p below 0.05 in {spatial.get('pair', np.nan):.1%} of pairs of random "
        f"maps and {spatial.get('map', np.nan):.1%} of random maps against nano "
        "(5% expected)",
    )
    null_rho_panel(fig.add_axes([0.05, 0.6, 0.25, 0.28]), calibration)
    variogram_panel(fig.add_axes([0.37, 0.6, 0.27, 0.28]), variogram, "nano")
    rates_panel(fig.add_axes([0.72, 0.6, 0.25, 0.28]), calibration)
    plane_axes = [fig.add_axes([0.02 + 0.172 * i, 0.12, 0.165, 0.4]) for i in range(4)]
    surrogate_planes(fig, plane_axes, map_ranks, surrogate_ranks, lab, names)
    plane_axes[0].text(
        0,
        1.15,
        f"C.  The nano map and three of its surrogates, CCF plane {plane}",
        transform=plane_axes[0].transAxes,
        fontsize=9,
    )
    if genes:
        gene_null_panel(fig.add_axes([0.78, 0.17, 0.2, 0.28]), genes)
    footer(
        fig,
        [
            "How to read: a surrogate is the nano map's ranks shuffled, smoothed over "
            "near structures and rescaled until its variogram matches the map's (Burt "
            "2020). A gene's spatial p is the share of surrogates that correlate with "
            "it at least as strongly as the map does.",
            "The calibration tests random smooth maps that have no relation to each "
            "other; a calibrated p is below 0.05 in about 5% of them.",
            "What would mean what: spatial p near 5% and ordinary p far above it: the "
            "null is needed and it works; a spatial p far from 5%: no claim rests on it.",
        ],
    )
    return saved(fig, save)


# ===== 07s P9's genes against the map, the gap, the two genes with the tissue =====


def gene_bars(
    ax: plt.Axes,
    table: pd.DataFrame,
    auto: pd.Series,
    subunits: set[str],
    t_max: float,
    q: float,
) -> None:
    """A: one bar per gene in rank order, its null band, autofluorescence's rho.

    `table` holds the genes' rows of gene_ranking.csv for the nano map, best first;
    `auto` the same genes' rho with the autofluorescence map.
    """
    y = np.arange(len(table))
    colours = bars_grey(table["t_boot"].to_numpy(), t_max)
    for i, (_, r) in enumerate(table.iterrows()):
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
    ax.barh(y, table["rho"], color=colours, height=0.62, zorder=2)
    ax.scatter(
        auto.reindex(table["symbol"]).to_numpy(),
        y,
        s=16,
        color=AUTO,
        edgecolors=AUTO_DOT,
        linewidths=0.6,
        zorder=3,
    )
    ax.axvline(0, color="0.3", lw=0.6, zorder=1)

    # past its band before correction: a small open circle at the left margin
    past = (table["p_spatial"] < ALPHA).to_numpy()
    ax.scatter(
        np.full(int(past.sum()), -0.72),
        y[past],
        s=12,
        facecolors="white",
        edgecolors=DARK_GREY,
        linewidths=0.7,
        zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(table["symbol"], fontsize=6.6)
    for tick, (_, r) in zip(ax.get_yticklabels(), table.iterrows()):
        tick.set_color(gene_colour(r["symbol"], subunits))
        if r["q_p9"] < q:
            tick.set_fontweight("bold")
    ax.set_ylim(len(table) - 0.4, -0.6)
    ax.set_xlim(-0.75, 1.0)
    ax.set_xlabel("Spearman rho with the adult map")
    ax.xaxis.set_ticks_position("both")
    ax.tick_params(axis="x", labeltop=True)

    # each gene's Allen reliability in a column at the right
    for i, (_, r) in enumerate(table.iterrows()):
        text = "-" if not np.isfinite(r["reliability"]) else f"{r['reliability']:.2f}"
        ax.text(1.02, i, text, fontsize=6, va="center", ha="left", color=DARK_GREY)
    tidy(ax)


def bar_key(ax: plt.Axes, t_max: float, q: float) -> None:
    """The key to panel A: what length, grey, band, dot and bold name each say."""
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    rows = [
        ("bar", "length: rho with the nano map"),
        ("grey", f"grey: rho over its SD across resampled adults (black at {t_max:g})"),
        ("band", "pale blue: 95% of rho with the map's surrogates"),
        ("circle", "open circle at the left: past its band, p < 0.05 before BH"),
        ("dot", "yellow dot: rho with the autofluorescence map"),
        ("bold", f"bold name: past the null after BH within P9's genes, q < {q}"),
        ("blue", "dark blue name: an AMPA receptor subunit"),
    ]
    for k, (kind, text) in enumerate(rows):
        yy = 0.93 - k * 0.135
        if kind == "bar":
            ax.add_patch(Rectangle((0.0, yy - 0.04), 0.12, 0.08, color="0.25"))
        elif kind == "grey":
            for j, t in enumerate(np.linspace(0, t_max, 5)):
                ax.add_patch(
                    Rectangle(
                        (0.025 * j, yy - 0.04), 0.024, 0.08, color=bars_grey(t, t_max)[0]
                    )
                )
        elif kind == "band":
            ax.add_patch(Rectangle((0.0, yy - 0.06), 0.12, 0.12, color=NULL_BAND))
        elif kind == "circle":
            ax.scatter(
                [0.06], [yy], s=16, facecolors="white", edgecolors=DARK_GREY, lw=0.7
            )
        elif kind == "dot":
            ax.scatter(
                [0.06], [yy], s=22, color=AUTO, edgecolors=AUTO_DOT, linewidths=0.6
            )
        elif kind == "bold":
            ax.text(0.0, yy, "Gene", fontweight="bold", fontsize=8, va="center")
        else:
            ax.text(0.0, yy, "Gria1", color=DARK_BLUE, fontsize=8, va="center")
        ax.text(0.16, yy, text, fontsize=8, va="center")


def gap_panel(
    ax: plt.Axes, gap: pd.DataFrame, null: np.ndarray, n_surrogates: int
) -> None:
    """B: the Cacng8 - Gria1 gap against its two nulls.

    With its adult interval, and the pairings of experiments as rugs beneath.

    `null` holds the null gaps of maps related to both genes alike (first row) and
    of maps unrelated to both (second row).
    """
    merged = gap[gap["kind"] == "merged profiles"].iloc[0]
    pairs = gap[gap["kind"] == "experiment pairing"]
    first, second = merged["first"], merged["second"]
    equal, unrelated = null[0], null[1]
    bins = np.linspace(-0.8, 0.8, 81)
    bound = merged["equal_two_sided"]
    ax.axvspan(-bound, bound, ymin=0.19, color=NULL_BAND, alpha=0.5, lw=0)
    ax.hist(equal, bins=bins, color=NULL_BAND, label="maps related to both genes alike")
    ax.hist(
        unrelated,
        bins=bins,
        histtype="step",
        color=MID_GREY,
        lw=1.1,
        label="maps unrelated to both (a conservative bound)",
    )
    top = ax.get_ylim()[1]
    ax.axvline(merged["gap"], color=RED, lw=1.8, label="observed gap")

    # the rugs: the adult interval and each pairing of the genes' Allen experiments
    ax.plot(
        [merged["boot_lo"], merged["boot_hi"]],
        [-0.08 * top] * 2,
        color=DARK_GREY,
        lw=3,
        solid_capstyle="butt",
        label="95% over resampled adults",
    )
    ax.scatter(
        pairs["gap"],
        [-0.18 * top] * len(pairs),
        marker="|",
        s=80,
        color="0.1",
        linewidths=1.4,
        label=f"each pairing of the two genes' Allen experiments ({len(pairs)})",
    )
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_ylim(-0.25 * top, top * 1.05)
    ax.set_xlim(-0.8, 0.8)
    ax.set_xlabel(f"rho({first}) - rho({second}), on the structures both have")
    ax.set_ylabel("null maps")
    ax.legend(loc="upper left", fontsize=7, frameon=False)
    panel_title(
        ax,
        "B",
        f"The {first} - {second} gap: {merged['gap']:+.3f} on "
        f"{int(merged['n_structures'])} structures ({merged['rho_first']:+.2f} against "
        f"{merged['rho_second']:+.2f})",
        f"maps related to both alike: {p_text(merged['p_equal'], n_surrogates)}, either "
        f"way (shaded: a lead under ±{bound:.2f}, which the test fixed in advance does "
        f"not pass); unrelated maps: {p_text(merged['p_spatial'], n_surrogates)}",
    )
    tidy(ax)


def pair_panel(ax: plt.Axes, ranking: pd.DataFrame, n_surrogates: int) -> None:
    """C: Gria1 and Cacng8 with the nano map and with autofluorescence, null bands."""
    rows = []
    for gene in ("Cacng8", "Gria1"):
        for name, colour in (("nano", NANO), ("auto", AUTO)):
            r = ranking[(ranking["map"] == name) & (ranking["symbol"] == gene)]
            if len(r):
                rows.append((gene, name, colour, r.iloc[0]))
    for i, (gene, name, colour, r) in enumerate(rows):
        ax.add_patch(
            Rectangle(
                (r["null_lo"], i - 0.4),
                r["null_hi"] - r["null_lo"],
                0.8,
                color=NULL_BAND,
                lw=0,
                zorder=0,
            )
        )
        ax.barh(i, r["rho"], color=colour, height=0.55, zorder=2)
        ax.text(
            max(r["rho"], r["null_hi"]) + 0.02,
            i,
            f"{r['rho']:+.2f}, {p_text(r['p_spatial'], n_surrogates)}",
            va="center",
            fontsize=7.5,
        )
    ax.set_yticks(range(len(rows)))
    labels = [
        f"{g}, {'nano' if n == 'nano' else 'autofluorescence'}" for g, n, _, _ in rows
    ]
    ax.set_yticklabels(labels, fontsize=8)
    for tick, (gene, _, _, _) in zip(ax.get_yticklabels(), rows):
        tick.set_color(DARK_BLUE if gene == "Gria1" else "0.1")
    ax.invert_yaxis()
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_xlim(-0.4, 1.15)
    ax.set_xlabel("Spearman rho; pale blue: 95% of each map's own surrogates")
    tidy(ax)


def counts_panel(ax: plt.Axes, ranking: pd.DataFrame, q: float) -> None:
    """D: genes past each map's null at three thresholds, nano against autofluorescence.

    Past the gene's band (p < 0.05 before correction), p of 0.005 or less, and past
    BH within all genes; beside each count, the same among P9's genes.
    """
    thresholds = (
        ("p < 0.05\n(past its band)", lambda t: t["p_spatial"] < ALPHA),
        ("p ≤ 0.005", lambda t: t["p_spatial"] <= 0.005),
        (f"BH q < {q}\n(all genes)", lambda t: t["q_all"] < q),
    )
    top = 0
    for i, (_, rule) in enumerate(thresholds):
        for j, (name, colour) in enumerate((("nano", NANO), ("auto", AUTO))):
            mine = ranking[ranking["map"] == name]
            n_pass = int(rule(mine).sum())
            n_p9 = int(rule(mine[mine["p9_gene"]]).sum())
            x = i + (j - 0.5) * 0.38
            ax.bar(x, n_pass, width=0.36, color=colour)
            ax.text(
                x, n_pass + 2, f"{n_pass}\n({n_p9})", ha="center", va="bottom", fontsize=7
            )
            top = max(top, n_pass)
    ax.set_ylim(0, top * 1.25)
    ax.set_xticks(range(len(thresholds)))
    ax.set_xticklabels([t[0] for t in thresholds], fontsize=8)
    ax.set_ylabel(f"genes of {ranking['symbol'].nunique()} (of P9's in brackets)")
    ax.legend(
        handles=[
            Patch(color=NANO, label="nano map"),
            Patch(color=AUTO, label="autofluorescence map"),
        ],
        loc="upper right",
        frameon=False,
    )
    tidy(ax)


def plot_gene_ranking(
    ranking: pd.DataFrame,
    gap: pd.DataFrame,
    gap_null: np.ndarray,
    subunits: set[str],
    t_max: float,
    q: float,
    n_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 07s: P9's genes against the map and its null, the gap, the tissue.

    The two genes with the tissue close it.

    `ranking` is gene_ranking.csv, `gap` gap.csv and `gap_null` the gap of every
    surrogate.
    """
    nano = ranking[(ranking["map"] == "nano") & ranking["p9_gene"]]
    nano = nano.sort_values("rho", ascending=False)
    auto = ranking[ranking["map"] == "auto"].set_index("symbol")["rho"]
    n_pass = int((nano["q_p9"] < q).sum())
    n_band = int((nano["p_spatial"] < ALPHA).sum())
    fig = plt.figure(figsize=(16, 17))
    heading(
        fig,
        "gene_ranking",
        f"P9's {len(nano)} genes on the declared structures: {n_band} past their band "
        f"(p < 0.05), {n_pass} after BH within them (q < {q}); {n_surrogates} "
        "surrogates per gene",
    )
    ax_a = fig.add_axes([0.08, 0.1, 0.36, 0.785])
    gene_bars(ax_a, nano, auto, subunits, t_max, q)
    ax_a.set_title(
        f"A.  P9's genes, best first: {n_band} past their band, {n_pass} after BH "
        "(bold)\nthe column at the right: the gene's Allen reliability ('-': measured "
        "once)\n\n",
        loc="left",
        fontsize=9,
    )
    bar_key(fig.add_axes([0.56, 0.77, 0.4, 0.16]), t_max, q)
    gap_panel(fig.add_axes([0.58, 0.44, 0.38, 0.26]), gap, gap_null, n_surrogates)
    ax_c = fig.add_axes([0.66, 0.16, 0.24, 0.17])
    pair_panel(ax_c, ranking, n_surrogates)
    panel_title(
        ax_c,
        "C",
        "The two genes with the label and with the tissue",
        "nano against the autofluorescence of the same sections "
        f"({figure_ref('autofluorescence')})",
    )
    footer(
        fig,
        [
            "How to read: a bar past its pale band is a rho that fewer than 5% of maps "
            "with the nano map's smoothness reach with that gene (open circle); BH then "
            f"allows for testing {len(nano)} genes (bold).",
            "The gap asks whether the map follows Cacng8 (a TARP, which brings AMPA "
            "receptors to the surface and holds them at synapses) more closely than "
            "Gria1, the receptor's own mRNA; its null is maps that follow both alike.",
            "What would mean what: Gria1 past its band and autofluorescence not (C): "
            "the map follows receptor expression, not the tissue. The gap past the null "
            "of maps related to both alike: the map follows a regulator of surface",
            "receptor more closely than the receptor's mRNA, as a surface-fraction "
            "reading predicts (a correspondence, not a measurement: Cacng8 is rich in "
            "hippocampus, as the map is). The gap inside it: these maps cannot tell "
            "the two apart.",
        ],
    )
    return saved(fig, save)


# ===== 12 and 12s Is the ranking the nanobody's or the tissue's =====

# genes always named on the scatter of figures 12 and 12s, beside those far from the
# diagonal

# the gene lists of figure 12s E: lines per column, and the width of a column in the
# panel's width
PASSING_LINES = 14
PASSING_WIDTH = 0.19


def auto_scatter(
    ax: plt.Axes,
    ranking: pd.DataFrame,
    subunits: set[str],
    n_named: int = 8,
    letter: str = "A",
    q: float | None = None,
) -> None:
    """A: each gene's rho with autofluorescence against its rho with nano.

    With `q`, the title also counts the genes past each map's null after BH.
    """
    nano = ranking[ranking["map"] == "nano"].set_index("symbol")
    auto = ranking[ranking["map"] == "auto"].set_index("symbol")
    genes = sorted(set(nano.index) & set(auto.index))
    x = nano.loc[genes, "rho"].to_numpy()
    y = auto.loc[genes, "rho"].to_numpy()
    p9 = nano.loc[genes, "p9_gene"].to_numpy(bool)
    sub = np.array([g in subunits for g in genes])
    lim = (-0.75, 0.95)
    ax.plot(lim, lim, color=MID_GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.axhline(0, color="0.88", lw=0.6, zorder=0)
    ax.axvline(0, color="0.88", lw=0.6, zorder=0)
    ax.scatter(x[~p9], y[~p9], s=14, color=LIGHT_GREY, linewidths=0, zorder=2)
    ax.scatter(x[p9 & ~sub], y[p9 & ~sub], s=24, color=DARK_GREY, linewidths=0, zorder=3)
    ax.scatter(x[sub], y[sub], s=30, color=DARK_BLUE, linewidths=0, zorder=4)

    # P9's genes furthest from the diagonal, and the genes always named
    gap = np.where(p9, np.abs(x - y), -1)
    named = [genes[i] for i in np.argsort(gap)[::-1][:n_named]]
    named += [g for g in QUOTED_GENES if g in genes and g not in named]
    agree = spearmanr(x, y).statistic
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    spread_labels(
        ax,
        [
            (x[genes.index(g)], y[genes.index(g)], g, gene_colour(g, subunits))
            for g in named
        ],
    )
    ax.set_aspect("equal")
    ax.set_xlabel("rho with the nano map")
    ax.set_ylabel("rho with the autofluorescence map")
    numbers = (
        f"{len(genes)} genes (P9's dark, subunits blue);\nthe two gene orders agree "
        f"at rho {agree:+.2f}"
    )
    if q is not None:
        past_nano = int((nano["q_all"] < q).sum())
        past_auto = int((auto["q_all"] < q).sum())
        numbers += f";\npast BH: nano {past_nano} genes, autofluorescence {past_auto}"
    panel_title(ax, letter, "Each gene against both maps", numbers)
    tidy(ax)


def rho_distributions(ax: plt.Axes, ranking: pd.DataFrame, q: float) -> None:
    """B: the genes' rho with each map, and how many pass each map's null."""
    bins = np.linspace(-0.8, 0.9, 35)
    lines = []
    for name, colour, label in (
        ("nano", NANO, "nano"),
        ("auto", AUTO, "autofluorescence"),
    ):
        mine = ranking[ranking["map"] == name]
        edge = NANO_DOT if name == "nano" else AUTO_DOT
        ax.hist(mine["rho"], bins=bins, histtype="step", color=edge, lw=1.6, label=label)
        n_all = int((mine["q_all"] < q).sum())
        p9 = mine[mine["p9_gene"]]
        n_p9 = int((p9["q_p9"] < q).sum())
        lines.append(
            f"{label}: median {mine['rho'].median():+.2f}; past BH {n_all} (P9's {n_p9})"
        )
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_xlabel("Spearman rho with the map")
    ax.set_ylabel("genes")
    ax.legend(loc="upper left")
    panel_title(ax, "B", "The genes' rho with each map", "\n".join(lines))
    tidy(ax)


def per_adult_panel(
    ax: plt.Axes, per_adult: pd.DataFrame, genes: tuple[str, ...], letter: str = "C"
) -> None:
    """C: per adult, rho of its nano and its autofluorescence map with each gene."""
    for k, gene in enumerate(genes):
        mine = per_adult[per_adult["symbol"] == gene]
        x0, x1 = 3 * k, 3 * k + 1
        for _, r in mine.iterrows():
            ax.plot([x0, x1], [r["rho_nano"], r["rho_auto"]], color=PAIR_LINE, lw=0.8)
        ax.scatter([x0] * len(mine), mine["rho_nano"], s=28, color=NANO_DOT, zorder=3)
        ax.scatter([x1] * len(mine), mine["rho_auto"], s=28, color=AUTO_DOT, zorder=3)
        ax.text(
            x0 + 0.5,
            1.0,
            f"{gene}\nnano {mine['rho_nano'].min():+.2f} to {mine['rho_nano'].max():+.2f}"
            f"\nauto {mine['rho_auto'].min():+.2f} to {mine['rho_auto'].max():+.2f}",
            ha="center",
            va="bottom",
            fontsize=7.5,
        )
    ticks = [x for k in range(len(genes)) for x in (3 * k, 3 * k + 1)]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["nano", "auto"] * len(genes))
    ax.set_xlim(-0.6, 3 * len(genes) - 1.4)
    ax.set_ylim(-0.3, 1.25)
    ax.axhline(0, color="0.88", lw=0.6, zorder=0)
    ax.set_ylabel("rho of one adult's map with the gene")
    n_mice = per_adult["mouse"].nunique()
    panel_title(ax, letter, "Adult by adult", f"{n_mice} adults, a line joins one brain")
    tidy(ax)


def passing_panel(ax: plt.Axes, ranking: pd.DataFrame, q: float) -> None:
    """D: the genes past each map's null (BH within every gene), with their rho."""
    ax.axis("off")
    column = 0
    for name, label in (("nano", "nano"), ("auto", "autofluorescence")):
        mine = ranking[(ranking["map"] == name) & (ranking["q_all"] < q)]
        mine = mine.sort_values("rho", ascending=False)
        lines = [
            f"{g:9s} {r:+.2f}{' *' if p9 else ''}"
            for g, r, p9 in zip(mine["symbol"], mine["rho"], mine["p9_gene"])
        ]
        ax.text(
            column * PASSING_WIDTH,
            1.0,
            f"{label}: {len(mine)}",
            va="top",
            ha="left",
            fontsize=8,
            fontweight="bold",
        )
        for start in range(0, max(len(lines), 1), PASSING_LINES):
            ax.text(
                column * PASSING_WIDTH,
                0.9,
                "\n".join(lines[start : start + PASSING_LINES]),
                va="top",
                ha="left",
                fontsize=7,
                family="monospace",
            )
            column += 1
    ax.set_title(
        f"E.  The genes past each map's null (BH within every gene, q < {q}); "
        "* P9's genes",
        loc="left",
        fontsize=9,
    )


def plot_autofluorescence(
    ranking: pd.DataFrame,
    per_adult: pd.DataFrame,
    subunits: set[str],
    q: float,
    n_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 12: whether Gria1 and Cacng8 follow the label or the tissue.

    Adult by adult, against each map's null, and gene by gene.

    `ranking` is gene_ranking.csv (both maps), `per_adult` the per-adult rho table.
    """
    genes = ("Gria1", "Cacng8")
    mine = per_adult[per_adult["symbol"].isin(genes)]

    # the adults in which nano's rho is above autofluorescence's for both genes
    nano_above = mine["rho_nano"] > mine["rho_auto"]
    above = int(nano_above.groupby(mine["mouse"]).all().sum())
    auto = ranking[(ranking["map"] == "auto") & ranking["symbol"].isin(genes)]
    inside = "inside" if (auto["p_spatial"] >= ALPHA).all() else "not all inside"
    fig = plt.figure(figsize=(16, 6.2))
    heading(
        fig,
        "autofluorescence",
        f"Gria1 and Cacng8 follow the label, not the tissue: nano's rho with both is "
        f"above autofluorescence's in {above} of {per_adult['mouse'].nunique()} adults, "
        f"and autofluorescence's is {inside} its null",
    )
    per_adult_panel(fig.add_axes([0.05, 0.15, 0.22, 0.6]), per_adult, genes, "A")
    ax = fig.add_axes([0.42, 0.25, 0.2, 0.45])
    pair_panel(ax, ranking, n_surrogates)
    panel_title(
        ax,
        "B",
        "Against each map's own null",
        "pale blue: 95% of the surrogates of that map",
    )
    auto_scatter(
        fig.add_axes([0.7, 0.13, 0.28, 0.64]), ranking, subunits, letter="C", q=q
    )
    footer(
        fig,
        [
            "How to read: the autofluorescence map is the unlabelled channel of the "
            "same sections, read as nano is and tested with surrogates of its own; a "
            "line in A joins one brain.",
            "The tissue has a gene pattern of its own, a different one (C); the genes "
            f"past each map's null: {figure_ref('autofluorescence_detail')}.",
        ],
    )
    return saved(fig, save)


def plot_autofluorescence_detail(
    ranking: pd.DataFrame,
    per_adult: pd.DataFrame,
    subunits: set[str],
    q: float,
    n_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 12s: the same ranking on the autofluorescence map of the same sections.

    `ranking` is gene_ranking.csv (both maps), `per_adult` the per-adult rho table.
    """
    fig = plt.figure(figsize=(16, 12.5))
    nano = ranking[ranking["map"] == "nano"].set_index("symbol")
    auto = ranking[ranking["map"] == "auto"].set_index("symbol")
    parts = [
        f"{g}: nano {nano.loc[g, 'rho']:+.2f} "
        f"({p_text(nano.loc[g, 'p_spatial'], n_surrogates)}), "
        f"autofluorescence {auto.loc[g, 'rho']:+.2f} "
        f"({p_text(auto.loc[g, 'p_spatial'], n_surrogates)})"
        for g in ("Gria1", "Cacng8")
        if g in nano.index and g in auto.index
    ]
    heading(fig, "autofluorescence_detail", "; ".join(parts))
    auto_scatter(fig.add_axes([0.05, 0.45, 0.33, 0.42]), ranking, subunits)
    rho_distributions(fig.add_axes([0.47, 0.6, 0.22, 0.27]), ranking, q)
    per_adult_panel(fig.add_axes([0.76, 0.6, 0.21, 0.27]), per_adult, ("Gria1", "Cacng8"))
    ax_d = fig.add_axes([0.47, 0.15, 0.3, 0.3])
    counts_panel(ax_d, ranking, q)
    widths = {
        name: float(
            (
                ranking.loc[ranking["map"] == name, "null_hi"]
                - ranking.loc[ranking["map"] == name, "null_lo"]
            ).median()
        )
        for name in ("nano", "auto")
    }
    panel_title(
        ax_d,
        "D",
        "How many genes pass, by threshold",
        f"median width of a gene's null band: nano {widths['nano']:.2f}, "
        f"autofluorescence {widths['auto']:.2f};\nBH counts only the smallest p values",
    )
    passing_panel(fig.add_axes([0.05, 0.07, 0.36, 0.27]), ranking, q)
    footer(
        fig,
        [
            "How to read: the autofluorescence map is the unlabelled channel of the "
            "same sections, read as nano is (zref over the declared set) and tested "
            "with surrogates of its own (A8); a line in C joins one brain.",
            "D counts the genes each map passes at three thresholds: the counts can "
            "rank the two maps differently, since BH rewards a few very small p values "
            "and the band counts every gene that leaves it.",
            "What would mean what: a gene past the nano null and not past the "
            "autofluorescence one, and above it in every adult (C), follows the label "
            "and not the tissue; genes past the autofluorescence null share",
            "the tissue's own pattern, and the nano ranking is the label's own as far "
            "as its gene order departs from the autofluorescence one (A).",
        ],
    )
    return saved(fig, save)


# ===== 13 and 13s Does the ranking depend on the choices made =====

# the headers of the kinds of choice in the robustness figure
KIND_HEADERS = {
    "primary": "the primary ranking",
    "statistic": "the statistic",
    "borders": "the borders of a structure",
    "reading": "the reading of the nano map",
    "inputs": "the ISH inputs",
    "structures": "the structure set",
    "5 October": "before this build",
}

# the variants drawn against the primary in panel D, one per kind of choice
SMALL_MULTIPLES = ("pearson_log2", "eroded_both", "ratio", "every_structure")


def variant_heights(summary: pd.DataFrame) -> np.ndarray:
    """The height of each variant's row: one apart, with a gap before each new kind."""
    y, rows, last = 0.0, [], None
    for kind in summary["kind"]:
        if last is not None and kind != last:
            y += 0.9
        rows.append(y)
        y += 1.0
        last = kind
    return np.array(rows)


def kind_headers(ax: plt.Axes, summary: pd.DataFrame, y: np.ndarray) -> None:
    """Each kind of choice named in grey above its rows, left of the row labels."""
    for kind in dict.fromkeys(summary["kind"]):
        first = y[(summary["kind"] == kind).to_numpy()].min()
        ax.text(
            -0.02,
            first - 0.75,
            KIND_HEADERS.get(kind, kind),
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=7.5,
            fontstyle="italic",
            color=DARK_GREY,
        )


def agreement_panel(ax: plt.Axes, summary: pd.DataFrame, y: np.ndarray) -> None:
    """A: each variant's gene order against the primary's, over P9's genes."""
    for yy, (_, r) in zip(y, summary.iterrows()):
        ax.plot([0.5, r["agreement_p9"]], [yy, yy], color=LIGHT_GREY, lw=1, zorder=1)
        ax.scatter(r["agreement_p9"], yy, s=40, color="0.15", zorder=2)
        ax.text(
            r["agreement_p9"] + 0.012,
            yy,
            f"{r['agreement_p9']:.3f}",
            va="center",
            fontsize=7.5,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(summary["label"], fontsize=8)
    ax.set_ylim(y.max() + 0.6, -1.4)
    kind_headers(ax, summary, y)
    ax.set_xlim(0.5, 1.08)
    ax.set_xlabel("agreement with the primary order\n(Spearman over P9's genes)")
    panel_title(ax, "A", "The order of the genes", "")
    tidy(ax)


def ranks_panel(ax: plt.Axes, summary: pd.DataFrame, y: np.ndarray) -> None:
    """B: where Cacng8 and Gria1 sit among P9's genes under each variant."""
    for gene, marker, colour in (("Cacng8", "o", RED), ("Gria1", "s", DARK_BLUE)):
        ax.scatter(
            summary[f"rank_p9_{gene}"],
            y,
            marker=marker,
            s=36,
            color=colour,
            label=gene,
            zorder=2,
        )
    for yy, (_, r) in zip(y, summary.iterrows()):
        ax.text(
            r["rank_p9_Gria1"] + 1,
            yy,
            f"{int(r['rank_p9_Gria1'])}",
            va="center",
            fontsize=7,
            color=DARK_BLUE,
        )
    ax.set_ylim(y.max() + 0.6, -1.4)
    ax.set_yticks([])
    top = float(np.nanmax(summary["rank_p9_Gria1"]))
    ax.set_xlim(0, top + 6)
    ax.set_xlabel("rank among P9's genes (1 = highest rho)")
    ax.legend(loc="lower right", fontsize=7)
    panel_title(ax, "B", "Where the two genes sit", "among P9's genes")
    tidy(ax)


def gap_rows_panel(
    ax: plt.Axes, summary: pd.DataFrame, y: np.ndarray, band: tuple[float, float]
) -> None:
    """C: the Cacng8 - Gria1 gap under each variant, the primary's two-sided test shaded.

    `band` is the lead either way that the primary's test, fixed in advance, does not
    pass (gap.csv equal_two_sided, as minus and plus).
    """
    ax.axvspan(*band, color=NULL_BAND, lw=0, zorder=0)
    ax.axvline(0, color="0.3", lw=0.6, zorder=1)
    ax.scatter(summary["gap"], y, s=30, color="0.15", zorder=3)
    for yy, g in zip(y, summary["gap"]):
        ax.text(g + 0.02, yy, f"{g:+.2f}", va="center", fontsize=7)
    ax.set_ylim(y.max() + 0.6, -1.4)
    ax.set_yticks([])
    ax.set_xlim(min(band[0], 0) - 0.05, max(band[1], summary["gap"].max()) + 0.12)
    ax.set_xlabel("rho(Cacng8) - rho(Gria1)")
    panel_title(
        ax,
        "C",
        "The Cacng8 - Gria1 gap",
        f"pale blue: a lead under ±{band[1]:.2f} either way,\nwhich the primary's "
        "two-sided test does not pass",
    )
    tidy(ax)


def multiples_panel(
    axes: list[plt.Axes],
    robustness: pd.DataFrame,
    summary: pd.DataFrame,
    subunits: set[str],
) -> None:
    """D: each gene's rho under four variants against its rho in the primary."""
    primary = robustness[robustness["variant"] == "primary"].set_index("symbol")
    labels = summary.set_index("variant")["label"]
    lim = (-0.75, 0.95)
    for ax, name in zip(axes, SMALL_MULTIPLES):
        mine = robustness[robustness["variant"] == name].set_index("symbol")
        genes = sorted(set(mine.index) & set(primary.index))
        x = primary.loc[genes, "rho"].to_numpy()
        y = mine.loc[genes, "rho"].to_numpy()
        p9 = primary.loc[genes, "p9_gene"].to_numpy(bool)
        sub = np.array([g in subunits for g in genes])
        ax.plot(lim, lim, color=MID_GREY, lw=0.8, ls=(0, (4, 3)))
        ax.scatter(x[~p9], y[~p9], s=8, color=LIGHT_GREY, linewidths=0)
        ax.scatter(x[p9], y[p9], s=14, color=DARK_GREY, linewidths=0)
        ax.scatter(x[sub], y[sub], s=18, color=DARK_BLUE, linewidths=0)
        for g in ("Cacng8", "Gria1"):
            if g in genes:
                i = genes.index(g)
                ax.annotate(
                    g,
                    (x[i], y[i]),
                    xytext=(-34, 4),
                    textcoords="offset points",
                    fontsize=6.5,
                    color=gene_colour(g, subunits),
                )
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal")
        ax.set_title(
            f"{labels[name]}\nrho {spearmanr(x, y).statistic:.3f}, {len(genes)} genes",
            fontsize=8,
        )
        ax.set_xlabel("rho, primary", fontsize=7.5)
        tidy(ax)
    axes[0].set_ylabel("rho, variant", fontsize=7.5)


def plot_robustness(
    summary: pd.DataFrame,
    gap_band: tuple[float, float],
    save: Path | None = None,
) -> plt.Figure:
    """Figure 13: the gene order, Cacng8's and Gria1's places and their gap, by choice.

    Under other statistics, borders, readings, inputs and structure sets.

    `summary` is robustness_summary.csv in the order of its rows, `gap_band` the 95%
    of the primary gap's null (gap.csv).
    """
    others = summary[summary["variant"] != "primary"]
    first = (
        "first of P9's genes" if others["rank_p9_Cacng8"].max() == 1 else "near the top"
    )
    fig = plt.figure(figsize=(16, 7.4))
    heading(
        fig,
        "robustness",
        f"Under every choice the gene order agrees at {others['agreement_p9'].min():.2f}"
        f" or more and Cacng8 stays {first}; Gria1's place among P9's genes moves "
        f"({ordinal(others['rank_p9_Gria1'].min())} to "
        f"{ordinal(others['rank_p9_Gria1'].max())}), so it is quoted with its null",
    )
    y = variant_heights(summary)
    agreement_panel(fig.add_axes([0.27, 0.15, 0.24, 0.66]), summary, y)
    ranks_panel(fig.add_axes([0.56, 0.15, 0.2, 0.66]), summary, y)
    gap_rows_panel(fig.add_axes([0.8, 0.15, 0.16, 0.66]), summary, y, gap_band)
    footer(
        fig,
        [
            "How to read: each row changes one choice of the primary ranking and keeps "
            "the others; the last row is the route of 5 October as it ran. Every gene "
            f"under four of the choices: {figure_ref('robustness_detail')}.",
        ],
    )
    return saved(fig, save)


def plot_robustness_detail(
    summary: pd.DataFrame,
    robustness: pd.DataFrame,
    subunits: set[str],
    gap_band: tuple[float, float],
    save: Path | None = None,
) -> plt.Figure:
    """Figure 13s: the ranking under other statistics, borders, readings, inputs, sets.

    `summary` is robustness_summary.csv in the order of its rows, `robustness`
    ranking_robustness.csv, `gap_band` the 95% of the primary gap's null (gap.csv).
    """
    fig = plt.figure(figsize=(16, 11.5))
    others = summary[summary["variant"] != "primary"]
    lo, hi = int(others["rank_p9_Cacng8"].min()), int(others["rank_p9_Cacng8"].max())
    cacng8 = f"Cacng8 {lo} to {hi}"
    if hi == 1:
        cacng8 = "Cacng8 first of P9's genes in every variant"
    heading(
        fig,
        "robustness_detail",
        f"{len(others)} variants of the primary ranking; the gene order agrees at "
        f"{others['agreement_p9'].min():.2f} to {others['agreement_p9'].max():.3f} "
        f"over P9's genes; Gria1 ranks {int(others['rank_p9_Gria1'].min())} to "
        f"{int(others['rank_p9_Gria1'].max())}, {cacng8}",
    )
    y = variant_heights(summary)
    agreement_panel(fig.add_axes([0.27, 0.47, 0.24, 0.41]), summary, y)
    ranks_panel(fig.add_axes([0.56, 0.47, 0.2, 0.41]), summary, y)
    gap_rows_panel(fig.add_axes([0.8, 0.47, 0.16, 0.41]), summary, y, gap_band)
    axes = [fig.add_axes([0.06 + 0.235 * k, 0.115, 0.19, 0.23]) for k in range(4)]
    multiples_panel(axes, robustness, summary, subunits)
    axes[0].text(
        0,
        1.25,
        "D.  Gene by gene, four of the variants against the primary (P9's genes "
        "dark, subunits blue)",
        transform=axes[0].transAxes,
        fontsize=9,
    )
    footer(
        fig,
        [
            "How to read: each row changes one choice of the primary ranking (Spearman, "
            "zref with the declared reference, merged profiles after QC, full means, "
            "the declared structures) and keeps the others;",
            "the last row is the route of 5 October as it ran. No row has a p: the "
            "null is drawn on the primary map's structures, and the rows ask about the "
            "order of the genes.",
            "What would mean what: an agreement near 1 everywhere: the conclusions do "
            "not hang on these choices; a row far below the others names the choice "
            "a conclusion depends on. The pale band is the primary's two-sided test "
            "(maps related to both alike); a variant's own null is not drawn, so a gap "
            "near the band's edge is a lead these maps cannot resolve.",
        ],
    )
    return saved(fig, save)
