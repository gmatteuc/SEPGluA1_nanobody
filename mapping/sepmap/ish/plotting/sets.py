"""The figures of analyses 2 and 3: between or within divisions, and kinds of genes.

Figure 09 and 09s show what a whole-brain rho is made of, the contrast between
divisions and the order inside them; a gene sheet shows one gene against the map,
on a plane and division by division. Figure 10 and 10s1 show whether the kinds of
genes that set surface receptor follow the map better than other genes, each gene
set against its null, and 10s2 the localisation genes against expression-matched
postsynaptic controls.

Called by run_ish_divisions.py (09, the gene sheets) and run_ish_gene_sets.py (10).
"""

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import spearmanr

from sepmap.ish.figure_index import figure_ref
from sepmap.ish.gene_sets import detectable
from sepmap.ish.plotting.shared import (
    DIVISION_ORDER,
    acronym_of,
    footer,
    gene_colour,
    group_handles,
    group_of,
    heading,
    p_text,
    panel_title,
    rank_plane,
    ranks01,
    saved,
    scatter_groups,
    set_column,
    spread_labels,
)
from sepmap.ish.spatial_null import ALPHA
from sepmap.plotting import (
    AUTO,
    AUTO_DOT,
    DARK_BLUE,
    DARK_GREY,
    DIVISION_GROUP,
    DIVISION_GROUP_COLOURS,
    LIGHT_GREY,
    MID_GREY,
    NOTE_GREY,
    NULL_BAND,
    PALE_BLUE,
    RED,
    SET_COLOURS,
    tidy,
)

# ===== 09 and 09s What a whole-brain rho is made of =====

# the marker of each detail gene in figure 09s D, all near black but the subunit's dark
# blue, so the genes do not take the colours of the groups of divisions
DETAIL_MARKERS = {
    "Cacng8": ("o", "0.1"),
    "Gria1": ("s", DARK_BLUE),
    "Grm5": ("^", "0.1"),
    "Dlg2": ("D", "0.1"),
    "Aqp4": ("v", MID_GREY),
}


def division_scatter(
    ax: plt.Axes,
    x_values: np.ndarray,
    gene_values: np.ndarray,
    groups: list[str],
    jitter: bool,
) -> None:
    """A map's ranks against a gene's, structures coloured by group of divisions.

    With `jitter`, the many structures that share one x (a division-only map) are
    spread a little so they can be seen.
    """
    x = ranks01(x_values)
    y = ranks01(gene_values)
    if jitter:
        x = x + np.random.default_rng(0).uniform(-0.012, 0.012, len(x))
    scatter_groups(ax, x, y, groups)
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.04, 1.04)
    tidy(ax)


def genes_scatter(
    ax: plt.Axes,
    x: pd.Series,
    y: pd.Series,
    p9: pd.Series,
    subunits: set[str],
    named: list[str],
    filled: pd.Series | None = None,
) -> None:
    """Every gene as a dot, P9's dark, subunits blue; hollow where `filled` is False."""
    lim = (-0.75, 0.95)
    ax.plot(lim, lim, color=MID_GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.axhline(0, color="0.88", lw=0.6, zorder=0)
    ax.axvline(0, color="0.88", lw=0.6, zorder=0)
    if filled is None:
        filled = pd.Series(True, index=x.index)
    sub = pd.Series(x.index.isin(subunits), index=x.index)
    for genes, size, colour in (
        (~p9 & ~sub, 12, LIGHT_GREY),
        (p9 & ~sub, 22, DARK_GREY),
        (sub, 28, DARK_BLUE),
    ):
        on = genes & filled
        off = genes & ~filled
        ax.scatter(x[on], y[on], s=size, color=colour, linewidths=0, zorder=2)
        ax.scatter(
            x[off],
            y[off],
            s=size,
            facecolors="white",
            edgecolors=colour,
            linewidths=0.7,
            zorder=2,
        )
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    spread_labels(
        ax,
        [(x[g], y[g], g, gene_colour(g, subunits)) for g in named if g in x.index],
    )
    ax.set_xlabel("whole-brain rho with the map")
    tidy(ax)


def detail_panel(
    ax: plt.Axes, detail: pd.DataFrame, within: pd.DataFrame, genes: tuple[str, ...]
) -> None:
    """D: each detail gene's rho inside every division, and its weighted mean."""
    divisions = [d for d in DIVISION_ORDER if d in set(detail["division"])]
    rows = divisions + ["weighted mean"]
    offsets = np.linspace(-0.25, 0.25, len(genes))
    means = within.set_index("symbol")["rho_within"]
    for k, gene in enumerate(genes):
        mine = detail[detail["symbol"] == gene].set_index("division")
        marker, colour = DETAIL_MARKERS.get(gene, ("o", DARK_GREY))
        for i, division in enumerate(divisions):
            if division in mine.index:
                ax.scatter(
                    mine.loc[division, "rho"],
                    i + offsets[k],
                    s=6 + 1.2 * mine.loc[division, "n_structures"],
                    marker=marker,
                    color=colour,
                    linewidths=0,
                    alpha=0.85,
                )
        ax.scatter(
            means.get(gene, np.nan),
            len(divisions) + offsets[k],
            s=55,
            marker=marker,
            facecolors="white",
            edgecolors=colour,
            linewidths=1.3,
            label=gene,
        )
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlabel("Spearman rho inside the division")
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=7.5)
    tidy(ax)


def within_rates_panel(ax: plt.Axes, calibration: pd.DataFrame) -> None:
    """E: false positives of the within rho on maps with no relation, by null."""
    nulls = (
        ("p_shuffle", "shuffled inside\ndivisions", LIGHT_GREY),
        ("p_spatial", "surrogates\n(spatial)", DARK_GREY),
    )
    for i, (column, _, colour) in enumerate(nulls):
        rate = float((calibration[column] < ALPHA).mean())
        ax.bar(i, rate, width=0.6, color=colour)
        ax.text(i, max(rate, 0.05) + 0.01, f"{rate:.1%}", ha="center", fontsize=8)
    ax.axhline(0.05, color=RED, ls="--", lw=1)
    ax.text(1.45, 0.055, "5%", color=RED, fontsize=7, ha="right")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([n[1] for n in nulls])
    ax.set_ylabel("share of tests with p < 0.05")
    tidy(ax)


def top_within_panel(
    ax: plt.Axes, table: pd.DataFrame, q: float, n_shown: int = 15
) -> None:
    """F: the genes highest inside divisions, their whole-brain rho and within q."""
    ax.axis("off")
    ranked = table[table["rho_within"].notna()].sort_values("rho_within", ascending=False)
    lines = ["gene       within  whole   q within"]
    for symbol, r in ranked.head(n_shown).iterrows():
        mark = "  past" if r["q_all_spatial"] < q else ""
        lines.append(
            f"{symbol:9s}  {r['rho_within']:+.2f}   {r['rho']:+.2f}   "
            f"{r['q_all_spatial']:.3f}{mark}{' *' if r['p9_gene'] else ''}"
        )
    ax.text(0, 1, "\n".join(lines), va="top", fontsize=7.5, family="monospace")
    ax.set_title(
        "F.  The genes highest inside divisions\nq: BH over all genes (past: "
        f"q < {q}); * P9's genes",
        loc="left",
        fontsize=9,
    )


def division_example(
    fig: plt.Figure,
    rects: tuple[list[float], list[float]],
    map_values: pd.Series,
    coarse: pd.Series,
    example: dict[str, float],
    example_name: str,
    table: pd.DataFrame,
    groups: dict[str, str],
) -> None:
    """A: one gene against the real map and against a map of divisions only.

    The second map knows only each structure's division; `rects` places the two axes.
    """
    shared = [s for s in map_values.index if s in example]
    y = np.array([example[s] for s in shared])
    g = [groups[s] for s in shared]
    ax0 = fig.add_axes(rects[0])
    ax1 = fig.add_axes(rects[1])
    division_scatter(ax0, map_values[shared].to_numpy(), y, g, jitter=False)
    division_scatter(ax1, coarse[shared].to_numpy(), y, g, jitter=True)
    r = table.loc[example_name]
    ax0.set_title(f"the real map: rho {r['rho']:+.2f}", fontsize=9)
    ax1.set_title(f"division only: rho {r['rho_division_only']:+.2f}", fontsize=9)
    ax0.set_xlabel("nano map, rank")
    ax1.set_xlabel("rank of the division's median")
    ax0.set_ylabel(f"{example_name}, rank")
    ax1.set_yticklabels([])
    ax0.text(
        0,
        1.2,
        f"A.  {example_name} against the map, and against a map that knows\nonly "
        "each structure's division (its median over the division)",
        transform=ax0.transAxes,
        fontsize=9,
    )
    ax1.legend(handles=group_handles(), loc="lower right", fontsize=6.5)


def division_only_panel(ax: plt.Axes, table: pd.DataFrame, subunits: set[str]) -> None:
    """B: every gene's whole-brain rho against its rho with the division-only map."""
    agree = spearmanr(table["rho"], table["rho_division_only"]).statistic
    genes_scatter(
        ax,
        table["rho"],
        table["rho_division_only"],
        table["p9_gene"],
        subunits,
        ["Cacng8", "Gria1", "Aqp4"],
    )
    ax.set_ylabel("rho with the division-only map")
    panel_title(
        ax,
        "B",
        "Most of a whole-brain rho is the contrast between divisions",
        f"the two gene orders agree at rho {agree:.2f} ({len(table)} genes)",
    )


def within_handles() -> list:
    """Legend entries of the within panel, each said once.

    The fill says whether a gene passes, the size and grey which panel it comes from.
    """
    entries = (
        (dict(color=MID_GREY, ms=5), "filled: past the within null (BH, all genes)"),
        (dict(mfc="white", mec=MID_GREY, ms=5), "hollow: not past it"),
        (dict(color=DARK_GREY, ms=np.sqrt(22)), "large, dark: P9's genes"),
        (dict(color=LIGHT_GREY, ms=np.sqrt(12)), "small, light: the other genes"),
        (dict(color=DARK_BLUE, ms=np.sqrt(28)), "AMPA receptor subunits"),
    )
    return [
        plt.Line2D([], [], ls="", marker="o", label=label, **style)
        for style, label in entries
    ]


def within_panel(
    ax: plt.Axes,
    table: pd.DataFrame,
    subunits: set[str],
    genes: tuple[str, ...],
    q: float,
    min_structures: int,
) -> None:
    """C: every gene's whole-brain rho against its mean rho inside divisions.

    Filled where the within rho passes its null; the median gene's null band shaded.
    """
    have = table["rho_within"].notna()
    agree = spearmanr(table.loc[have, "rho"], table.loc[have, "rho_within"]).statistic
    band = (table["null_lo"].median(), table["null_hi"].median())
    ax.axhspan(*band, color=NULL_BAND, lw=0, zorder=0)
    top_within = list(table.loc[have, "rho_within"].nlargest(5).index)
    passed = table["q_all_spatial"] < q
    genes_scatter(
        ax,
        table.loc[have, "rho"],
        table.loc[have, "rho_within"],
        table.loc[have, "p9_gene"],
        subunits,
        list(dict.fromkeys(top_within + list(genes))),
        filled=passed[have],
    )
    ax.set_ylabel(f"mean rho inside divisions ({min_structures}+ structures each)")
    panel_title(
        ax,
        "C",
        "Inside divisions the order changes, and rho shrinks",
        f"gene orders agree at rho {agree:.2f};\npale band: 95% of the median gene's "
        "within null",
    )
    ax.legend(handles=within_handles(), loc="lower right", fontsize=6.8, frameon=False)


def plot_between_within(
    within: pd.DataFrame,
    map_values: pd.Series,
    coarse: pd.Series,
    example: dict[str, float],
    example_name: str,
    set_table: pd.DataFrame,
    subunits: set[str],
    genes: tuple[str, ...],
    q: float,
    min_structures: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 09: whether genes follow the map inside divisions, or only between them.

    The arguments are those of plot_between_within_detail, less the tables of its
    panels D and E.
    """
    table = within.set_index("symbol")
    agree_div = spearmanr(table["rho"], table["rho_division_only"]).statistic
    past = table[table["q_all_spatial"] < q].sort_values("rho_within", ascending=False)
    named = ", ".join(f"{g} {r:+.2f}" for g, r in past["rho_within"].head(3).items())
    fig = plt.figure(figsize=(16, 6.4))
    heading(
        fig,
        "between_within",
        "Most of a whole-brain rho is the contrast between divisions (a division-only "
        f"map orders the genes at {agree_div:.2f}); inside divisions {len(past)} genes "
        f"still follow the map past its null, {named} first",
    )
    division_example(
        fig,
        ([0.05, 0.13, 0.15, 0.55], [0.22, 0.13, 0.15, 0.55]),
        map_values,
        coarse,
        example,
        example_name,
        table,
        group_of(set_table),
    )
    division_only_panel(fig.add_axes([0.44, 0.13, 0.22, 0.62]), table, subunits)
    within_panel(
        fig.add_axes([0.75, 0.13, 0.22, 0.62]), table, subunits, genes, q, min_structures
    )
    footer(
        fig,
        [
            "How to read: the division-only map gives each structure its division's "
            "median; the within rho is Spearman inside each division with enough "
            "structures, averaged by their number, so no contrast between divisions "
            "can enter it.",
            "Five genes division by division, and the choice of null: "
            f"{figure_ref('between_within_detail')}.",
        ],
    )
    return saved(fig, save)


def plot_between_within_detail(
    within: pd.DataFrame,
    detail: pd.DataFrame,
    map_values: pd.Series,
    coarse: pd.Series,
    example: dict[str, float],
    example_name: str,
    set_table: pd.DataFrame,
    calibration: pd.DataFrame,
    subunits: set[str],
    genes: tuple[str, ...],
    q: float,
    min_structures: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 09s: what a whole-brain rho is made of, between and within divisions.

    `within` and `detail` are the tables of run_ish_divisions, `map_values` the map
    on the declared structures and `coarse` its division-only version, `example`
    the profile of the gene of panel A, `calibration` the within rho of random maps.
    """
    table = within.set_index("symbol")
    n_pass = int((table["q_all_spatial"] < q).sum())
    n_pass_p9 = int((table.loc[table["p9_gene"], "q_p9_spatial"] < q).sum())
    fig = plt.figure(figsize=(16, 12.5))
    heading(
        fig,
        "between_within_detail",
        f"{len(table)} genes; median rho {table['rho'].median():+.2f} over the whole "
        f"brain, {table['rho_within'].median():+.2f} inside divisions; {n_pass} genes "
        f"past the within null at BH q < {q} within all genes, {n_pass_p9} within "
        "P9's",
    )
    division_example(
        fig,
        ([0.05, 0.58, 0.15, 0.27], [0.22, 0.58, 0.15, 0.27]),
        map_values,
        coarse,
        example,
        example_name,
        table,
        group_of(set_table),
    )
    division_only_panel(fig.add_axes([0.44, 0.55, 0.24, 0.32]), table, subunits)
    within_panel(
        fig.add_axes([0.74, 0.55, 0.24, 0.32]), table, subunits, genes, q, min_structures
    )

    # D: the detail genes division by division; E: the two nulls on random maps;
    # F: the top of the within ranking
    ax_d = fig.add_axes([0.08, 0.1, 0.28, 0.35])
    detail_panel(ax_d, detail, within, genes)
    panel_title(
        ax_d,
        "D",
        "Division by division, five genes",
        "marker size: structures in the division; open marker: the weighted mean",
    )
    ax_e = fig.add_axes([0.47, 0.13, 0.14, 0.29])
    within_rates_panel(ax_e, calibration)
    panel_title(
        ax_e,
        "E",
        "Which null for the within rho",
        f"{len(calibration)} random smooth maps\nwith no relation to the map",
    )
    top_within_panel(fig.add_axes([0.7, 0.1, 0.29, 0.33]), table, q)
    footer(
        fig,
        [
            "How to read: the division-only map gives each structure its division's "
            "median; the within rho is Spearman inside each division with enough "
            "structures, averaged by their number, so no contrast between divisions "
            "can enter it.",
            "It is tested against surrogates of the whole map (spatial), and against "
            "the map's values shuffled inside each division; E shows which of the two "
            "keeps 5% false positives on maps with no relation to the map.",
            "What would mean what: a high rho within divisions: gene and map vary "
            "together at a fine scale, the stronger claim; a whole-brain rho that "
            "vanishes inside divisions: they share only the gross gradient.",
        ],
    )
    return saved(fig, save)


# ===== Gene sheets =====


def division_lines(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    divisions: list[str],
    detail: pd.DataFrame,
) -> None:
    """One least-squares line per division the gene enters, its rho in the legend.

    Divisions of one group share its colour, so each takes its own line style.
    """
    divisions = np.asarray(divisions)
    styles = ("-", "--", ":", "-.")
    used = {}
    for _, r in detail.iterrows():
        mine = divisions == r["division"]
        slope, intercept = np.polyfit(x[mine], y[mine], 1)
        xs = np.array([x[mine].min(), x[mine].max()])
        group = DIVISION_GROUP.get(r["division"], "other grey matter")
        k = used.get(group, 0)
        used[group] = k + 1
        ax.plot(
            xs,
            intercept + slope * xs,
            color=DIVISION_GROUP_COLOURS[group],
            ls=styles[k % len(styles)],
            lw=1.5,
            label=f"{r['division']}: rho {r['rho']:+.2f} ({int(r['n_structures'])})",
        )
        # the division named at the line's end, since a group's lines share a colour
        ax.text(
            xs[1] + 0.01,
            intercept + slope * xs[1],
            r["division"],
            fontsize=6.5,
            color=DIVISION_GROUP_COLOURS[group],
            va="center",
        )


def plot_gene_sheet(
    symbol: str,
    map_values: pd.Series,
    profile: dict[str, float],
    set_table: pd.DataFrame,
    ranking: pd.Series,
    within: pd.Series,
    detail: pd.DataFrame,
    lab: np.ndarray,
    names: dict[int, str],
    plane: int,
    n_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """One gene against the map: both as ranks on a plane, the scatter, the divisions.

    `ranking` and `within` are the gene's rows of gene_ranking.csv (nano map) and
    within_division.csv, `detail` its rows of within_division_detail.csv.
    """
    shared = [s for s in map_values.index if s in profile]
    x = ranks01(map_values[shared].to_numpy())
    y = ranks01(np.array([profile[s] for s in shared]))
    groups = group_of(set_table)
    division = dict(zip(set_table["structure"], set_table["division"]))
    acronyms = acronym_of(set_table)
    fig = plt.figure(figsize=(15, 9.5))
    fig.text(
        0.5,
        0.975,
        f"{symbol} against the adult nano map: whole-brain rho {ranking['rho']:+.2f} "
        f"(spatial {p_text(ranking['p_spatial'], n_surrogates)}); inside divisions "
        f"{within['rho_within']:+.2f} (spatial "
        f"{p_text(within['p_within_spatial'], n_surrogates)}, shuffled inside "
        f"divisions {p_text(within['p_within_shuffle'], n_surrogates)})",
        ha="center",
        va="top",
        fontsize=11,
    )
    ax = fig.add_axes([0.01, 0.5, 0.3, 0.38])
    rank_plane(fig, ax, dict(zip(shared, x)), lab, names)
    panel_title(
        ax,
        "A",
        f"The nano map as ranks, CCF plane {plane}",
        f"among the {len(shared)} declared structures {symbol} has; flat grey: others",
    )
    ax = fig.add_axes([0.01, 0.05, 0.3, 0.38])
    rank_plane(fig, ax, dict(zip(shared, y)), lab, names)
    panel_title(ax, "B", f"{symbol}'s merged profile as ranks", "the same structures")

    # C: the scatter, a line per division, the structures furthest from the diagonal
    ax = fig.add_axes([0.4, 0.08, 0.38, 0.8])
    scatter_groups(ax, x, y, [groups[s] for s in shared])
    division_lines(ax, x, y, [division[s] for s in shared], detail)
    ax.plot([0, 1], [0, 1], color=MID_GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
    for i in np.argsort(np.abs(y - x))[::-1][:6]:
        ax.annotate(
            acronyms[shared[i]],
            (x[i], y[i]),
            xytext=(4, 3),
            textcoords="offset points",
            fontsize=7.5,
        )
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_aspect("equal")
    ax.set_xlabel("nano map, rank among structures")
    ax.set_ylabel(f"{symbol}, rank among structures")
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=7.5,
        title="a line per division (structures)",
        title_fontsize=8,
    )
    panel_title(
        ax,
        "C",
        f"{len(shared)} declared structures; the six furthest from the diagonal named",
    )
    tidy(ax)
    ax.text(
        1.02,
        0.25,
        "\n".join(
            textwrap.wrap(
                "A line rising inside a division means the gene and the map order its "
                "structures alike there; lines that sit apart with no slope inside "
                "carry a rho made of the contrast between divisions "
                f"({figure_ref('between_within')}). "
                "Dots: structures, coloured by group of divisions.",
                38,
            )
        ),
        transform=ax.transAxes,
        fontsize=7.5,
        color=NOTE_GREY,
        va="top",
    )
    return saved(fig, save)


# ===== 10 and 10s1 Which kinds of genes match =====


def set_panel(
    ax: plt.Axes,
    members: pd.DataFrame,
    tests: pd.DataFrame,
    order: list[str],
    context: str,
    q: float,
    n_surrogates: int,
) -> None:
    """A: one column per set, the context group last, each with its null band."""
    rng = np.random.default_rng(0)
    nano = tests[tests["map"] == "nano"].set_index("gene_set")
    auto = tests[tests["map"] == "auto"].set_index("gene_set")
    labels = []
    for k, name in enumerate(order + [context]):
        rows = members[members["gene_set"] == name]
        if name == context:
            set_column(ax, k, rows, None, MID_GREY, rng, hollow=True)
            labels.append(f"{name}\n(context, not tested)\nn = {len(rows)}")
            continue
        test = nano.loc[name]
        set_column(ax, k, rows, test, SET_COLOURS[name], rng)
        ax.scatter(
            k + 0.42,
            auto.loc[name, "median_rho"],
            marker="D",
            s=26,
            color=AUTO,
            edgecolors=AUTO_DOT,
            linewidths=0.6,
            zorder=4,
        )
        if test["tested"]:
            weight = "bold" if test["q"] < q else "normal"
            stats = f"{p_text(test['p_spatial'], n_surrogates)}, q {test['q']:.3f}"
        else:
            weight = "normal"
            stats = "not tested: too few genes,\nno band"
        labels.append((f"{name}\nn = {int(test['n_genes'])}\n{stats}", weight))
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([lab if isinstance(lab, str) else lab[0] for lab in labels])
    for tick, lab in zip(ax.get_xticklabels(), labels):
        if not isinstance(lab, str):
            tick.set_fontweight(lab[1])
    ax.axhline(0, color="0.6", lw=0.6, zorder=1)
    ax.set_xlim(-0.6, len(labels) - 0.3)
    ax.set_ylabel("Spearman rho of the gene with the adult map")
    handles = [
        Patch(color=NULL_BAND, label="95% of the set's median over the surrogates"),
        plt.Line2D([], [], color="0.1", lw=2, label="the set's median"),
        plt.Line2D(
            [],
            [],
            ls="",
            marker="o",
            mfc="white",
            mec="0.05",
            label="a gene of P9's panel (black ring)",
        ),
        plt.Line2D(
            [],
            [],
            ls="",
            marker="D",
            mfc=AUTO,
            mec=AUTO_DOT,
            label="the set's median with the autofluorescence map",
        ),
    ]
    ax.legend(handles=handles, loc="lower left", fontsize=7)
    tidy(ax)


def contrast_panel(
    ax: plt.Axes, row: pd.Series, null: np.ndarray, n_surrogates: int
) -> None:
    """B: one contrast named in advance against its surrogate differences, 95% shaded."""
    bins = np.linspace(-0.8, 0.8, 65)
    ax.axvspan(row["null_lo"], row["null_hi"], color=NULL_BAND, alpha=0.5, lw=0)
    ax.hist(null, bins=bins, color=NULL_BAND)
    ax.axvline(row["difference"], color=RED, lw=1.8)
    ax.axvline(0, color="0.5", lw=0.6)
    ax.set_xlabel("difference of the two sets' median rho")
    ax.set_ylabel("surrogates")
    ax.set_title(
        f"{row['contrast']}\n{row['median_first']:+.2f} ({int(row['n_first'])} genes) "
        f"against {row['median_second']:+.2f} ({int(row['n_second'])}); difference "
        f"{row['difference']:+.2f}\nspatial {p_text(row['p_spatial'], n_surrogates)}; "
        f"labels permuted p = {row['p_labels']:.4f}",
        loc="left",
        fontsize=8.5,
    )
    tidy(ax)


def sets_panel(
    ax: plt.Axes,
    members: pd.DataFrame,
    tests: pd.DataFrame,
    order: list[str],
    n_surrogates: int,
) -> None:
    """A of figure 10: each set fixed in advance, a dot per gene, against its null.

    Its median against the band where the same genes' median falls with 95% of the
    surrogates.
    """
    rng = np.random.default_rng(0)
    nano = tests[tests["map"] == "nano"].set_index("gene_set")
    labels = []
    for k, name in enumerate(order):
        test = nano.loc[name]
        rows = members[members["gene_set"] == name]
        set_column(ax, k, rows, test, SET_COLOURS[name], rng)
        if test["tested"]:
            stats = f"{p_text(test['p_spatial'], n_surrogates)}\nq = {test['q']:.3f}"
        else:
            stats = "too few to test"
        labels.append(f"{name}\nn = {int(test['n_genes'])}\n{stats}")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.axhline(0, color="0.6", lw=0.6, zorder=1)
    ax.set_xlim(-0.6, len(labels) - 0.3)
    ax.set_ylabel("Spearman rho of the gene with the nano map")
    tidy(ax)


def plot_gene_kinds(
    members: pd.DataFrame,
    tests: pd.DataFrame,
    order: list[str],
    table: pd.DataFrame,
    summary: pd.DataFrame,
    pairs: dict[str, str],
    power: pd.DataFrame,
    q: float,
    n_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 10: whether the genes that set surface receptor follow the map better.

    The gene sets against the null, the localisation genes against postsynaptic genes of
    the same expression, and what that test could find.

    `members` and `tests` are gene_sets.csv and set_tests.csv, `order` the sets in
    drawing order; `table`, `summary`, `pairs` and `power` those of
    plot_localisation.
    """
    nano = tests[(tests["map"] == "nano") & tests["tested"]]
    n_pass = int((nano["q"] < q).sum())
    matched = summary.set_index("test").loc["matched controls"]
    found = detectable(power, "power_labels")
    if n_pass:
        sets = f"{n_pass} of {len(nano)} gene sets pass their null"
    else:
        sets = "No gene set passes its null"
    if matched["p_labels"] < ALPHA:
        local = "the localisation genes stand above postsynaptic genes of the same"
    else:
        local = "the localisation genes do no better than postsynaptic genes of the same"
    fig = plt.figure(figsize=(16, 6.4))
    heading(
        fig,
        "gene_kinds",
        f"{sets}, and {local} expression ({matched['difference']:+.3f}, p = "
        f"{matched['p_labels']:.2f}); a difference of {found:+.2f} would have been found",
    )

    # A: the sets fixed in advance
    ax = fig.add_axes([0.05, 0.2, 0.42, 0.56])
    sets_panel(ax, members, tests, order, n_surrogates)
    panel_title(
        ax,
        "A",
        "Kinds of genes, fixed before looking",
        "pale band: where the median falls with 95% of the surrogates",
    )

    # B: localisation genes against their matched controls, the subunits removed
    by = table.set_index("symbol")
    loc = list(table.loc[table["side"] == "localisation", "symbol"])
    paired = sorted(set(pairs.values()))
    ax = fig.add_axes([0.56, 0.2, 0.17, 0.56])
    strips(
        ax,
        [
            ("localisation", by.loc[loc, "rho_partial"].to_numpy(), RED),
            ("matched\ncontrols", by.loc[paired, "rho_partial"].to_numpy(), DARK_GREY),
        ],
        np.random.default_rng(0),
        "partial rho with the map, subunits removed",
    )
    panel_title(
        ax,
        "B",
        "Localisation against matched controls",
        f"{matched['difference']:+.3f} between medians, p = {matched['p_labels']:.2f}",
    )

    # C: what the test could find
    ax = fig.add_axes([0.79, 0.2, 0.19, 0.56])
    power_panel(ax, power, float(matched["difference"]))
    panel_title(
        ax,
        "C",
        "What the test could find",
        f"found in 80% of maps from {found:+.3f}",
    )
    footer(
        fig,
        [
            "How to read: localisation genes (AMPA receptor transport, anchoring, "
            "auxiliary subunits) are paired with the postsynaptic gene closest in "
            "expression; partial rho once the subunit composite is out of both sides.",
            f"In detail: every set and the contrasts named in advance, "
            f"{figure_ref('gene_sets')}; the positive control and every test of the "
            f"design, {figure_ref('localisation')}.",
        ],
    )
    return saved(fig, save)


def plot_gene_sets(
    members: pd.DataFrame,
    tests: pd.DataFrame,
    contrasts: pd.DataFrame,
    contrast_nulls: dict[str, np.ndarray],
    order: list[str],
    context: str,
    rules: dict[str, str],
    q: float,
    n_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 10s1: the gene sets fixed in advance against the map and its null.

    `members`, `tests` and `contrasts` are gene_sets.csv, set_tests.csv and
    contrasts.csv; `contrast_nulls` the surrogate differences of each contrast;
    `order` the sets in drawing order, `context` the group drawn last, untested,
    `rules` where each set comes from.
    """
    nano = tests[(tests["map"] == "nano") & tests["tested"]]
    n_pass = int((nano["q"] < q).sum())
    fig = plt.figure(figsize=(16, 13))
    heading(
        fig,
        "gene_sets",
        f"{len(order)} sets from GO and cited marker lists, fixed before any rho on "
        f"these inputs; {n_pass} of {len(nano)} tested sets past the spatial null at "
        f"BH q < {q}",
    )
    ax = fig.add_axes([0.06, 0.46, 0.9, 0.44])
    set_panel(ax, members, tests, order, context, q, n_surrogates)
    panel_title(
        ax,
        "A",
        "One column per set, a dot per gene",
        "bold set name: its median past the null at the BH level; pale band: where "
        "the median of the same genes falls with 95% of the surrogates",
    )
    for k, (_, row) in enumerate(contrasts.iterrows()):
        ax_b = fig.add_axes([0.08 + 0.3 * k, 0.16, 0.24, 0.135])
        contrast_panel(ax_b, row, contrast_nulls[row["contrast"]], n_surrogates)
    passed = int((contrasts["p_spatial"] < q).sum())
    fig.text(
        0.06,
        0.35,
        "B.  The two contrasts named in advance, against the surrogates (red: observed; "
        f"shaded: 95% of the null): {passed} of {len(contrasts)} past it; the criterion "
        "named in advance needs both",
        fontsize=9,
    )
    pre = members[members["gene_set"] == "presynaptic"].sort_values(
        "rho_nano", ascending=False
    )
    fig.text(
        0.08,
        0.118,
        textwrap.fill(
            f"The presynaptic side ({len(pre)} genes): "
            + ", ".join(f"{g} {r:+.2f}" for g, r in zip(pre["symbol"], pre["rho_nano"]))
            + ". GO annotates the vesicle machinery (Syn1, Vamp2, Stx1a, Bsn, Syt1) "
            "to both sides of the synapse, so the rule leaves it out of this set.",
            150,
        ),
        fontsize=7.5,
        color=DARK_GREY,
        va="top",
    )
    ax_c = fig.add_axes([0.68, 0.14, 0.3, 0.19])
    ax_c.axis("off")
    text = "\n".join(f"{name}: {rules.get(name, '')}" for name in order + [context])
    ax_c.text(
        0,
        1,
        "\n".join(
            textwrap.fill(line, 70, subsequent_indent="   ") for line in text.split("\n")
        ),
        va="top",
        fontsize=7,
    )
    ax_c.set_title("C.  Where each set comes from", loc="left", fontsize=9)
    footer(
        fig,
        [
            "How to read: a set's median is tested against the medians the same genes "
            "give with every surrogate of the map, so co-expressed genes stay together "
            "in the null; the labels-permuted p treats the genes as independent and is "
            "shown for comparison only.",
            "What would mean what: postsynaptic above presynaptic and above glia, both "
            "beyond the null: the map is postsynaptic-like, as any glutamate receptor "
            "label should be (a sanity check). A set past its band follows the map more",
            "than a random smooth map would; for the surface reading, localisation must "
            "stand above other postsynaptic genes of the same expression "
            f"({figure_ref('localisation')}), not only past its own band.",
        ],
    )
    return saved(fig, save)


# ===== 10s2 Localisation genes against matched controls =====


def strips(
    ax: plt.Axes,
    groups: list[tuple[str, np.ndarray, object]],
    rng: np.random.Generator,
    ylabel: str,
) -> None:
    """Dots of each group in a column, jittered, with the group's median as a bar."""
    for i, (label, values, colour) in enumerate(groups):
        ax.scatter(
            i + rng.uniform(-0.15, 0.15, len(values)),
            values,
            s=14,
            color=colour,
            linewidths=0,
            alpha=0.85,
            zorder=2,
        )
        ax.plot([i - 0.3, i + 0.3], [np.median(values)] * 2, color="0.1", lw=2, zorder=3)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([f"{g[0]}\n({len(g[1])})" for g in groups])
    ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax.set_ylabel(ylabel)
    tidy(ax)


def label_null_panel(ax: plt.Axes, null: np.ndarray, row: pd.Series) -> None:
    """The label-permutation null of a difference, the observed one, the critical."""
    ax.hist(null, bins=60, color=LIGHT_GREY)
    ax.axvline(row["difference"], color=RED, lw=1.8)
    for side in (-1, 1):
        ax.axvline(side * row["critical"], color=DARK_GREY, ls="--", lw=0.9)
    ax.set_xlabel("difference of the medians, labels permuted")
    ax.set_ylabel("permutations")
    tidy(ax)


def power_panel(ax: plt.Axes, power: pd.DataFrame, observed: float) -> None:
    """C of figure 10: the share of maps on which each test finds a planted effect."""
    for column, colour, label in (
        ("power_labels", DARK_GREY, "labels permuted (the test)"),
        ("power_spatial", NULL_BAND, "against the remainder's surrogates"),
    ):
        edge = DARK_GREY if column == "power_labels" else PALE_BLUE
        ax.plot(
            power["median_difference"],
            power[column],
            color=edge,
            marker="o",
            ms=4,
            lw=1.4,
            label=label,
        )
    ax.axhline(0.8, color=MID_GREY, ls=":", lw=0.9)
    ax.axhline(0.05, color=RED, ls="--", lw=0.8)
    ax.text(
        power["median_difference"].max(),
        0.06,
        "5%",
        color=RED,
        fontsize=7,
        ha="right",
        va="bottom",
    )
    ax.axvline(observed, color=RED, lw=1.6)
    ax.text(observed, 1.02, " observed", color=RED, fontsize=7.5, va="bottom")
    ax.set_ylim(0, 1.08)
    ax.set_xlabel("matched difference the maps carry (median partial rho)")
    ax.set_ylabel("share of maps with p < 0.05")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.1), fontsize=7, frameon=False)
    tidy(ax)


def plot_localisation(
    table: pd.DataFrame,
    summary: pd.DataFrame,
    nulls: dict[str, np.ndarray],
    pairs: dict[str, str],
    p_spatial: float,
    n_surrogates: int,
    power: pd.DataFrame | None = None,
    panel_table: pd.DataFrame | None = None,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 10s2: localisation genes against expression-matched postsynaptic controls.

    `table` is localisation_test.csv, `summary` localisation_summary.csv, `nulls`
    each test's label null, `pairs` the matching, `p_spatial` the matched
    difference's p against surrogates of the map's remainder; `power` the power
    check (localisation_power.csv); `panel_table`, when given, the same genes with
    the control pool of 5 October, whose positive control panel C draws beside.
    """
    rng = np.random.default_rng(0)
    by = table.set_index("symbol")
    loc = list(table.loc[table["side"] == "localisation", "symbol"])
    ctrl = list(table.loc[table["side"] == "control", "symbol"])
    paired = sorted(set(pairs.values()))
    test = summary.set_index("test")
    matched = test.loc["matched controls"]
    positive = test.loc["positive control"]
    october = "5 October's controls, positive control"
    status = [
        f"positive control, GO pool {positive['difference']:+.3f} "
        f"(p = {positive['p_labels']:.3f})"
    ]
    if october in test.index:
        status.append(
            f"5 October's pool {test.loc[october, 'difference']:+.3f} "
            f"(p = {test.loc[october, 'p_labels']:.3f})"
        )
    found = ""
    if power is not None:
        found = (
            f"; found in 80% of maps from a difference of "
            f"{detectable(power, 'power_labels'):+.3f}"
        )
    fig = plt.figure(figsize=(16, 13))
    heading(
        fig,
        "localisation",
        f"The test: {len(loc)} localisation genes against {len(paired)} "
        f"expression-matched controls, {matched['difference']:+.3f} (p = "
        f"{matched['p_labels']:.3f}; spatial {p_text(p_spatial, n_surrogates)}){found}. "
        "The " + "; ".join(status),
    )
    ax = fig.add_axes([0.06, 0.58, 0.22, 0.29])
    strips(
        ax,
        [
            ("localisation", by.loc[loc, "rho_partial"].to_numpy(), RED),
            ("matched\ncontrols", by.loc[paired, "rho_partial"].to_numpy(), DARK_GREY),
            ("all\ncontrols", by.loc[ctrl, "rho_partial"].to_numpy(), LIGHT_GREY),
        ],
        rng,
        "partial rho with the map, subunit composite removed",
    )
    panel_title(
        ax,
        "A",
        "The test",
        f"medians {matched['median_first']:+.3f} against {matched['median_second']:+.3f}",
    )
    ax = fig.add_axes([0.36, 0.58, 0.26, 0.29])
    label_null_panel(ax, nulls["matched controls"], matched)
    panel_title(
        ax,
        "B",
        "Its null, labels permuted between the two sets",
        f"p = {matched['p_labels']:.3f}; dashed: a difference beyond "
        f"±{matched['critical']:.3f} reaches p < 0.05\nagainst maps whose remainder is a "
        f"surrogate of the real one: spatial {p_text(p_spatial, n_surrogates)}",
    )
    ax = fig.add_axes([0.69, 0.58, 0.28, 0.29])
    groups = []
    pools = [("GO", table, positive)]
    if panel_table is not None and october in test.index:
        pools.append(("5 Oct.", panel_table, test.loc[october]))
    for label, pool, row in pools:
        controls = pool[pool["side"] == "control"]
        cut = row["reliability_cut"]
        high = controls.loc[controls["reliability"] >= cut, "rho"].abs()
        low = controls.loc[controls["reliability"] < cut, "rho"].abs()
        groups.append((f"{label}:\nreproducible", high.to_numpy(), DARK_GREY))
        groups.append((f"{label}:\nunreproducible", low.to_numpy(), LIGHT_GREY))
    strips(ax, groups, rng, "|rho| with the map")
    lines = [
        f"{label} controls split at reliability {row['reliability_cut']:.2f}: "
        f"{row['difference']:+.3f}, p = {row['p_labels']:.3f}"
        for label, _, row in pools
    ]
    panel_title(
        ax,
        "C",
        "The positive control: a difference that must exist (kept from April)",
        "\n".join(lines),
    )
    ax = fig.add_axes([0.06, 0.18, 0.22, 0.26])
    level = by["median_energy"]
    strips(
        ax,
        [
            ("localisation", np.log10(level[loc].to_numpy() + 1e-3), RED),
            ("matched\ncontrols", np.log10(level[paired].to_numpy() + 1e-3), DARK_GREY),
            ("all\ncontrols", np.log10(level[ctrl].to_numpy() + 1e-3), LIGHT_GREY),
        ],
        rng,
        "log10 median expression energy",
    )
    panel_title(ax, "D", "The matching", "a quiet gene correlates with nothing")
    if power is not None:
        ax = fig.add_axes([0.36, 0.18, 0.26, 0.26])
        power_panel(ax, power, float(matched["difference"]))
        zero = power[power["effect"] == 0].iloc[0]
        panel_title(
            ax,
            "E",
            "What the test can find",
            "maps with a planted localisation pattern; found in 80% from "
            f"{detectable(power, 'power_labels'):+.3f};\nwith no pattern the label "
            f"test gives p < 0.05 on {zero['power_labels']:.1%} of maps",
        )
    ax = fig.add_axes([0.69, 0.16, 0.3, 0.28])
    ax.axis("off")
    lines = [
        f"{r['test'][:40]:40s} {r['role'][:9]:9s} {r['difference']:+.3f}  "
        f"{r['p_labels']:.3f}"
        for _, r in summary.iterrows()
    ]
    top = table[table["side"] == "localisation"].nlargest(5, "rho_partial")
    lines += ["", "localisation genes with the highest partial rho:"]
    lines += [f"  {r['symbol']:10s} {r['rho_partial']:+.3f}" for _, r in top.iterrows()]
    header = f"{'test':40s} {'role':9s} {'diff.':6s}  p"
    ax.text(
        0, 1, header + "\n" + "\n".join(lines), va="top", fontsize=7, family="monospace"
    )
    ax.set_title(
        "F.  Every test of the design\n(first side: localisation genes, or the "
        "reproducible controls)",
        loc="left",
        fontsize=9,
    )
    footer(
        fig,
        [
            "How to read: partial rho is the rank correlation of a gene with the map "
            "once the mean rank of the four AMPA subunits has been regressed out of "
            "both (the design of 5 October; F also removes the four as separate terms);",
            "each localisation gene is paired with the control closest in expression "
            f"(controls: the other postsynaptic genes of {figure_ref('gene_sets')}; 5 "
            "Oct.: the ontology panel's postsynaptic-density genes, secondary, added "
            "after the GO pool's positive control failed).",
            "E plants the pattern that sets the localisation genes apart from their "
            "controls into maps that otherwise have the real map's smoothness, and "
            "counts how often the test finds it.",
            "What would mean what: localisation above its matched controls: the genes "
            "that set surface receptor predict the part of the map the subunits do not, "
            "as a surface-fraction reading predicts; a difference inside B,",
            "with E able to find one of the size the reading needs: no evidence for it "
            "from these genes, bounded by the smallest difference E finds.",
        ],
    )
    return saved(fig, save)
