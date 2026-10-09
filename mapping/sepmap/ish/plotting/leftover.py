"""The figures of the genes that follow the map, and of the tests for its leftover.

Figure 07 describes the genes past the map's null on every axis of the ISH line and
what each takes of what Gria1 and synapse density leave; 08 sets the map against
Gria1 and Cacng8 and their gap against maps that follow both alike; 11 holds the
three tiers named in advance for the leftover (Cacng8, the AMPA receptor complex
family, every other gene) and 11s2 the family member by member and on the map. A
top-gene sheet describes one gene on one page.

Called by run_ish_top_genes.py.
"""

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from sepmap.ish.figure_index import figure_ref
from sepmap.ish.gene_sets import LEFTOVER_GENE
from sepmap.ish.plotting.shared import (
    LOW_RELIABILITY,
    acronym_of,
    footer,
    group_handles,
    group_of,
    heading,
    p_text,
    panel_title,
    rank_plane,
    ranks01,
    saved,
    scatter_groups,
    spread_labels,
)
from sepmap.ish.spatial_null import ALPHA
from sepmap.ish.top_genes import MATCHED_TEST, OUTSIDE_TEST, SPATIAL_TEST, WITHOUT_TEST
from sepmap.plotting import (
    AUTO,
    DARK_BLUE,
    DARK_GREY,
    DENSITY_BLUE,
    LIGHT_GREY,
    MID_GREY,
    NANO,
    NULL_BAND,
    PAIR_LINE,
    RED,
    bars_grey,
    tidy,
)

# ===== 07 The genes that follow the map =====

# the band behind the rows of Gria1 and Cacng8, across every panel of figure 07
ROW_TINT = "0.93"

# what a gene is set beside in figure 07 C and on its sheet: column, label, colour
# and marker; a column the main model does not have is skipped
LIKENESS = (
    ("rho_Gria1", "Gria1", DARK_BLUE, "o"),
    ("rho_density", "density term", DENSITY_BLUE, "o"),
    ("rho_autofluo", "autofluorescence", AUTO, "s"),
    ("rho_prediction", "the model's prediction", "0.1", "|"),
)


def named_colour(symbol: str) -> str:
    """Red for Cacng8, the gene named for the leftover; dark blue for Gria1."""
    if symbol == "Cacng8":
        return RED
    if symbol == "Gria1":
        return DARK_BLUE
    return "0.1"


def row_bands(ax: plt.Axes, lo, hi, height: float = 0.84) -> None:
    """A pale band per row from `lo` to `hi`: the 95% of a null behind each value."""
    for i, (a, b) in enumerate(zip(lo, hi)):
        ax.add_patch(
            Rectangle((a, i - height / 2), b - a, height, color=NULL_BAND, lw=0, zorder=0)
        )


def tint_rows(ax: plt.Axes, rows: pd.DataFrame) -> None:
    """The rows of Gria1 and Cacng8 on a light grey band, so they read across panels."""
    for i, symbol in enumerate(rows["symbol"]):
        if symbol in ("Gria1", "Cacng8"):
            ax.axhspan(i - 0.5, i + 0.5, color=ROW_TINT, lw=0, zorder=-1)


def gene_rows_axis(ax: plt.Axes, rows: pd.DataFrame, labels: bool) -> None:
    """One row per gene, top to bottom; the names on the first panel only.

    A gene in the main model has its part of it in brackets after its name.
    """
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_yticks(np.arange(len(rows)))
    if not labels:
        ax.set_yticklabels([])
        return
    names = [f"{s}  ({m})" if m else s for s, m in zip(rows["symbol"], rows["in_model"])]
    ax.set_yticklabels(names, fontsize=8)
    for tick, symbol in zip(ax.get_yticklabels(), rows["symbol"]):
        tick.set_color(named_colour(symbol))
        if symbol in ("Gria1", "Cacng8"):
            tick.set_fontweight("bold")


def top_rows(table: pd.DataFrame) -> pd.DataFrame:
    """The rows of figure 07: the genes past the map's null and the named, best first."""
    rows = table[table["top"] | table["named"]]
    return rows.sort_values("rho", ascending=False).reset_index(drop=True)


def own_term(rows: pd.DataFrame) -> np.ndarray:
    """Whether each gene is a term of the main model itself (Gria1, the abundance).

    Its rho with the leftover is then near zero by construction, and its null so
    narrow that a p says nothing; a gene inside the density term (one of the
    genes averaged into it) is not.
    """
    return (rows["in_model"] == "abundance").to_numpy()


def follows_leftover(rows: pd.DataFrame) -> np.ndarray:
    """Whether each gene that is not a term itself follows the leftover at p < 0.05."""
    return (rows["leftover_p"] < ALPHA).to_numpy() & ~own_term(rows)


def map_bars(ax: plt.Axes, rows: pd.DataFrame, t_max: float) -> None:
    """A: each gene's rho with the map against its null band.

    A red diamond at the right when the gene also follows the leftover (own_term
    genes never do); a gene in the main model has its part of it in brackets.
    """
    y = np.arange(len(rows))
    tint_rows(ax, rows)
    row_bands(ax, rows["null_lo"], rows["null_hi"])
    ax.barh(
        y,
        rows["rho"],
        color=bars_grey(rows["t_boot"].to_numpy(), t_max),
        height=0.62,
        zorder=2,
    )
    follows = follows_leftover(rows)
    ax.scatter(np.full(int(follows.sum()), 1.06), y[follows], marker="D", s=30, color=RED)
    ax.axvline(0, color="0.3", lw=0.6, zorder=1)
    ax.set_xlim(-0.6, 1.14)
    ax.set_xticks([-0.5, 0, 0.5, 1.0])
    ax.set_xlabel("Spearman rho with the nano map")
    gene_rows_axis(ax, rows, labels=True)
    tidy(ax)


def within_dots(ax: plt.Axes, rows: pd.DataFrame, q: float) -> None:
    """B: each gene's mean rho inside divisions against its null; filled past BH."""
    y = np.arange(len(rows))
    tint_rows(ax, rows)
    row_bands(ax, rows["within_null_lo"], rows["within_null_hi"])
    past = (rows["q_within"] < q).to_numpy()
    ax.scatter(rows["rho_within"][past], y[past], s=34, color=DARK_GREY, zorder=3)
    ax.scatter(
        rows["rho_within"][~past],
        y[~past],
        s=34,
        facecolors="white",
        edgecolors=DARK_GREY,
        linewidths=1.0,
        zorder=3,
    )
    ax.axvline(0, color="0.3", lw=0.6, zorder=1)
    ax.set_xlim(-0.4, 0.8)
    ax.set_xlabel("mean rho inside divisions")
    gene_rows_axis(ax, rows, labels=False)
    tidy(ax)


def likeness_dots(ax: plt.Axes, rows: pd.DataFrame) -> None:
    """C: each gene's rho with Gria1, the density terms and the model's prediction."""
    y = np.arange(len(rows))
    tint_rows(ax, rows)
    for column, label, colour, marker in LIKENESS:
        if column not in rows:
            continue
        bar = marker == "|"

        # a gene against itself (Gria1 against the Gria1 term) is left out
        values = rows[column].where(rows["symbol"] != column.removeprefix("rho_"))
        ax.scatter(
            values,
            y,
            marker=marker,
            s=70 if bar else 26,
            color=colour,
            linewidths=1.6 if bar else 0,
            label=label,
            zorder=3,
        )
    ax.axvline(0, color="0.3", lw=0.6, zorder=1)
    ax.set_xlim(-0.4, 1.05)
    ax.set_xlabel("Spearman rho over the structures of the fit")
    ax.legend(loc="lower left", fontsize=7, frameon=False, handletextpad=0.2)
    gene_rows_axis(ax, rows, labels=False)
    tidy(ax)


def points(share: float) -> str:
    """A share of the reproducible map in points, one decimal, never '-0.0'."""
    return f"{round(100 * share, 1) + 0.0:.1f}"


def taken_bars(ax: plt.Axes, rows: pd.DataFrame) -> None:
    """D: the share of the leftover each gene takes, against maps of its smoothness.

    In points of the reproducible map. The pale band runs from 0 to what 95% of the
    plain surrogates of the gene take, the black tick marks what 95% of the maps
    alike to the model take (ish.top_genes). Every bar is one grey, since its
    darkness would otherwise read as reliability; a gene that takes more than the
    band has a filled triangle after its value, more than the tick only an open
    one. A term of the model says so instead of a bar.
    """
    y = np.arange(len(rows))
    tint_rows(ax, rows)
    taken = 100 * rows["taken"].to_numpy(float)
    plain = 100 * rows["taken_null_hi"].to_numpy(float)
    alike = 100 * rows["taken_null_alike_hi"].to_numpy(float)
    own = own_term(rows)
    row_bands(ax, np.zeros(len(rows)), np.where(own, 0, plain))
    ax.barh(y[~own], taken[~own], color=MID_GREY, height=0.62, zorder=2)
    ax.scatter(
        alike[~own], y[~own], marker="|", s=120, color="0.05", linewidths=1.6, zorder=3
    )
    marks = zip(rows["p_taken"], rows["p_taken_alike"], own)
    for i, (p_plain, p_alike, term) in enumerate(marks):
        if term:
            ax.text(0.3, i, "in the model", fontsize=6.5, va="center", color=DARK_GREY)
            continue
        mark = ""
        if p_plain < ALPHA:
            mark = "  ▲"
        elif p_alike < ALPHA:
            mark = "  △"
        ax.text(
            max(taken[i], plain[i], alike[i]) + 0.3,
            i,
            points(taken[i] / 100) + mark,
            fontsize=6.5,
            va="center",
            color=DARK_GREY,
        )
    ax.axvline(0, color="0.3", lw=0.6, zorder=1)
    # room right of the longest bar or band for its value and mark
    right = max(float(np.max(plain)) * 1.7, float(np.max(taken)) * 1.45)
    ax.set_xlim(min(0.0, float(np.min(taken)) - 0.5), right)
    ax.set_xlabel("points of the reproducible map")
    gene_rows_axis(ax, rows, labels=False)
    tidy(ax)


def plot_top_genes(
    table: pd.DataFrame,
    q: float,
    t_max: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 07: the genes that follow the map, what each is, what each takes.

    What kind of map each gene is, and what it takes of the leftover.

    `table` is top_genes.csv; the rows are the genes past the map's null after BH,
    and Gria1 and Cacng8.
    """
    rows = top_rows(table)
    n = len(rows)
    top = rows[rows["top"]]
    follows = rows.loc[follows_leftover(rows), "symbol"]
    height = 3.4 + 0.42 * n
    fig = plt.figure(figsize=(16, height))
    heading(
        fig,
        "top_genes",
        f"{len(top)} genes follow the map past its null, Cacng8 first; they are maps "
        f"much like Gria1 and synapse density (rho {top['rho_prediction'].min():.2f} to "
        f"{top['rho_prediction'].max():.2f} with the model's prediction), and "
        f"{len(follows)} also follow what the model leaves",
    )
    bottom = 1.15 / height
    span = 1 - 1.5 / height - bottom
    ax = fig.add_axes([0.12, bottom, 0.28, span])
    map_bars(ax, rows, t_max)
    panel_title(
        ax,
        "A",
        "With the map",
        "pale: 95% of maps with its smoothness;\nbar darkness: steadiness across "
        f"adults (black at t = {t_max:g});\nred diamond: also follows the leftover",
    )
    ax = fig.add_axes([0.44, bottom, 0.14, span])
    within_dots(ax, rows, q)
    panel_title(ax, "B", "Inside divisions", "filled: past its null")
    ax = fig.add_axes([0.62, bottom, 0.15, span])
    likeness_dots(ax, rows)
    panel_title(ax, "C", "Like Gria1, or like density?", "rho with each term")
    ax = fig.add_axes([0.82, bottom, 0.14, span])
    taken_bars(ax, rows)
    panel_title(
        ax,
        "D",
        "What it adds to the model",
        "pale: 95% of plain surrogates;\ntick: 95% of maps alike to the model;\n"
        "▲ past both, △ past the tick only",
    )
    footer(
        fig,
        [
            "How to read: the leftover is what Gria1 and synapse density leave; D adds "
            "the gene to that model as one more straight term and counts the points of "
            "the reproducible map it takes, against maps of its smoothness in its place.",
            "Maps alike to the model keep the gene's fit on the model's terms and put a "
            "surrogate of the rest in its place; that null was added after the plain "
            "one's numbers were seen, so it is the secondary line.",
            "Only Cacng8's test against the leftover was named in advance. Every gene "
            "in detail: the sheets in top_genes/.",
        ],
    )
    return saved(fig, save)


# ===== 08 Cacng8 against Gria1 =====


def against_gene(
    ax: plt.Axes,
    gene: dict[str, float],
    values: pd.Series,
    groups: dict[str, str],
    labels: tuple[str, str],
    ranks: bool = True,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """`values` (by structure) against a gene's profile, one dot per structure.

    Over the structures both have; the gene as ranks 0 to 1, `values` too when
    `ranks` (identity dashed), else as they are (a leftover in ranks, zero dashed).
    Returns x, y and the structures.
    """
    shared = [s for s in values.index if np.isfinite(gene.get(s, np.nan))]
    x = ranks01(np.array([gene[s] for s in shared]))
    y = values[shared].to_numpy(float)
    if ranks:
        y = ranks01(y)
        ax.plot([0, 1], [0, 1], color=MID_GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax.set_ylim(-0.04, 1.04)
        ax.set_aspect("equal")
    else:
        ax.axhline(0, color=MID_GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
    scatter_groups(ax, x, y, [groups.get(s, "other grey matter") for s in shared])
    ax.set_xlim(-0.04, 1.04)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    tidy(ax)
    return x, y, shared


def gap_verdict(gap: pd.Series) -> str:
    """Figure 08's takeaway on the gap: the test fixed in advance, two-sided.

    Only its p decides the words. The band and the shares each way of the skewed null
    are read in panel C, never in the line under the title.
    """
    p = gap["p_equal"]
    where = "past" if p < ALPHA else "inside"
    return (
        f"their gap, {gap['gap']:+.2f}, is {where} the two-sided test fixed in advance "
        f"(maps that follow both alike, p = {p:.3f})"
    )


def lead_panel(
    ax: plt.Axes, gap: pd.Series, gap_null: np.ndarray, n_map_surrogates: int
) -> None:
    """C of figure 08: the Cacng8 - Gria1 gap against maps that follow both alike.

    The test fixed in advance is two-sided, so the band shaded is the lead either way
    that it does not pass (equal_two_sided), the gap is drawn with its mirror, and the
    share of these maps beyond each is written beside it, as a description.
    """
    bound = gap["equal_two_sided"]
    ax.hist(gap_null, bins=np.linspace(-0.6, 0.6, 61), color=NULL_BAND)
    ax.axvspan(-bound, bound, color=NULL_BAND, alpha=0.4, lw=0)
    ax.axvline(gap["gap"], color=RED, lw=1.8)
    ax.axvline(-gap["gap"], color=RED, lw=1.0, ls=(0, (3, 2)))
    top = ax.get_ylim()[1]
    shares = (
        (
            gap["gap"] + 0.02,
            "left",
            f"Cacng8 leads\nby as much:\n{gap['equal_first_as_large']:.1%}",
        ),
        (
            -gap["gap"] - 0.02,
            "right",
            f"Gria1 leads\nby as much:\n{gap['equal_second_as_large']:.1%}",
        ),
    )
    for x, ha, text in shares:
        ax.text(x, 0.95 * top, text, fontsize=7.5, va="top", ha=ha, color=RED)
    ax.set_xlim(-0.6, 0.6)
    ax.set_xlabel("rho(Cacng8) - rho(Gria1)")
    ax.set_ylabel("maps that follow both alike")
    panel_title(
        ax,
        "C",
        "The Cacng8 - Gria1 gap",
        f"{gap['gap']:+.3f} on {int(gap['n_structures'])} structures (red; its mirror "
        f"dashed);\nshaded: a lead under ±{bound:.3f} either way, which the test\nfixed "
        f"in advance does not pass; {p_text(gap['p_equal'], n_map_surrogates)}, "
        "counting either way",
    )
    tidy(ax)


def plot_cacng8_gria1(
    table: pd.DataFrame,
    map_values: pd.Series,
    profiles: dict[str, dict[str, float]],
    gap: pd.Series,
    gap_null: np.ndarray,
    set_table: pd.DataFrame,
    n_map_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 08: the map against Gria1 and Cacng8, and the gap between them.

    The gap is set against maps that follow both genes alike.

    `table` is top_genes.csv indexed by symbol, `map_values` the map on the declared
    structures, `gap` the merged row of gap.csv and `gap_null` the gaps of maps that
    follow both genes alike; `n_map_surrogates` counts the map's surrogates.
    """
    cacng8, gria1 = table.loc["Cacng8"], table.loc["Gria1"]
    groups = group_of(set_table)
    fig = plt.figure(figsize=(16, 6.4))
    heading(
        fig,
        "cacng8_gria1",
        f"The map follows Cacng8 ({cacng8['rho']:+.2f}) and Gria1 "
        f"({gria1['rho']:+.2f}); {gap_verdict(gap)}",
    )
    width, height, bottom = 0.22, 0.55, 0.22

    # A and B: the map against each gene
    genes = (
        ("Gria1", gria1, "A", "Gria1, the stained protein's mRNA"),
        ("Cacng8", cacng8, "B", "Cacng8, TARP gamma-8"),
    )
    for k, (gene, row, letter, what) in enumerate(genes):
        ax = fig.add_axes([0.05 + 0.31 * k, bottom, width, height])
        against_gene(
            ax,
            profiles[gene],
            map_values,
            groups,
            (f"{gene}, rank among structures", "nano map, rank among structures"),
        )
        panel_title(
            ax,
            letter,
            what,
            f"rho {row['rho']:+.2f} on {int(row['n_structures'])} structures; inside "
            f"divisions {row['rho_within']:+.2f}",
        )
    fig.legend(
        handles=group_handles(),
        loc="lower left",
        bbox_to_anchor=(0.05, 0.075),
        ncol=4,
        frameon=False,
        fontsize=8,
    )

    # C: the lead against maps that follow both genes alike
    ax = fig.add_axes([0.72, bottom, 0.25, height])
    lead_panel(ax, gap, gap_null, n_map_surrogates)
    footer(
        fig,
        [
            "How to read: one dot per structure, coloured by group of divisions. C: maps "
            "made of both genes' patterns alike plus a surrogate of the map; "
            f"{figure_ref('gene_ranking')} B adds the adults' interval and each pairing "
            "of Allen experiments.",
            "The test fixed in advance counts a lead as large either way; the null is "
            "skewed, so the shares beyond the gap and beyond its mirror are written "
            "beside them, as a description, not a test.",
            "Cacng8 against what Gria1 and synapse density leave: "
            f"{figure_ref('leftover')}.",
        ],
    )
    return saved(fig, save)


# ===== 11 and 11s2 The tests named for the leftover, the AMPA receptor complex =====

# the sources of a family member, abbreviated after its name in figure 11s2
SOURCE_TAGS = (("Schwenk", "S"), ("GO:", "GO"), ("partner", "P"))


def source_tag(sources: str) -> str:
    """The abbreviations of a member's sources: 'S GO', 'GO', 'S GO P'."""
    return " ".join(tag for key, tag in SOURCE_TAGS if key in sources)


def family_bars(
    ax: plt.Axes, rows: pd.DataFrame, q: float, t_max: float, n_surrogates: int
) -> None:
    """A: each member's rho with the leftover against its null band, best first."""
    y = np.arange(len(rows))
    row_bands(ax, rows["leftover_null_lo"], rows["leftover_null_hi"])
    ax.barh(
        y,
        rows["leftover_rho"],
        color=bars_grey(rows["leftover_t_boot"].to_numpy(), t_max),
        height=0.62,
        zorder=2,
    )
    ax.axvline(0, color="0.3", lw=0.6, zorder=1)
    labels = []
    for _, r in rows.iterrows():
        text = f"{r['symbol']}  {source_tag(r['family_sources'])}"
        if r["in_model"]:
            text += f"  (in the model: {r['in_model']})"
        labels.append(text)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    for tick, (_, r) in zip(ax.get_yticklabels(), rows.iterrows()):
        colour = DARK_BLUE if "partner" in r["family_sources"] else "0.1"
        tick.set_color(RED if r["symbol"] == "Cacng8" else colour)
        if r["leftover_q_family"] < q:
            tick.set_fontweight("bold")
    for i, r in rows.iterrows():
        ax.text(
            0.62,
            i,
            p_text(r["leftover_p"], n_surrogates),
            fontsize=6.5,
            va="center",
            color=DARK_GREY,
        )
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlim(-0.6, 0.75)
    ax.set_xlabel("Spearman rho with the leftover")
    tidy(ax)


def group_null_panel(ax: plt.Axes, row: pd.Series, null: np.ndarray, colour) -> None:
    """The family's median rho against the same genes' median over the surrogates."""
    ax.hist(null, bins=50, color=NULL_BAND)
    ax.axvspan(row["null_lo"], row["null_hi"], color=NULL_BAND, alpha=0.4, lw=0)
    ax.axvline(row["first"], color=colour, lw=1.8)
    ax.set_xlabel("median rho of the family's genes")
    ax.set_ylabel("surrogates")
    tidy(ax)


def paired_strips(
    ax: plt.Axes, rows: pd.DataFrame, first: str, second: str, rng: np.random.Generator
) -> None:
    """The members and their matched controls, joined pair by pair, with medians.

    `first` and `second` are the columns of the members' and the controls' rho.
    """
    a = rows[first].to_numpy(float)
    b = rows[second].to_numpy(float)
    xa = rng.uniform(-0.12, 0.12, len(a))
    xb = 1 + rng.uniform(-0.12, 0.12, len(b))
    for i in range(len(a)):
        ax.plot([xa[i], xb[i]], [a[i], b[i]], color=PAIR_LINE, lw=0.5, zorder=1)
    ax.scatter(xa, a, s=22, color=RED, linewidths=0, zorder=3)
    ax.scatter(xb, b, s=22, color=DARK_GREY, linewidths=0, zorder=3)
    for x, values in ((0, a), (1, b)):
        ax.plot(
            [x - 0.25, x + 0.25], [np.median(values)] * 2, color="0.1", lw=2, zorder=4
        )
    ax.axhline(0, color="0.6", lw=0.6, zorder=0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"the family\n({len(a)})", f"matched postsynaptic\ncontrols ({len(b)})"]
    )
    ax.set_xlim(-0.5, 1.5)
    tidy(ax)


def bh_panel(ax: plt.Axes, genes: pd.DataFrame, family: list[str], q: float) -> None:
    """D of figure 11: every gene's spatial p against the leftover, smallest first.

    Beside the line Benjamini-Hochberg asks a gene to cross; the family dark, Cacng8 red.
    """
    order = genes.sort_values("p_spatial")
    p = order["p_spatial"].to_numpy(float)
    k = np.arange(1, len(p) + 1)
    symbols = order["symbol"].to_numpy()
    in_family = np.isin(symbols, family)
    cacng8 = symbols == LEFTOVER_GENE
    above = order["rho"].to_numpy(float) > 0
    ax.plot(k, -np.log10(k * q / len(p)), color=RED, ls="--", lw=1, label=f"BH, q {q}")

    # the sign of each gene's rho as the marker: up above zero, down below
    kinds = (
        (~in_family, LIGHT_GREY, 12, "other genes"),
        (in_family & ~cacng8, DARK_GREY, 18, "the family"),
    )
    for mine, colour, size, label in kinds:
        for sign, marker in ((above, "^"), (~above, "v")):
            show = mine & sign
            ax.scatter(
                k[show],
                -np.log10(p[show]),
                s=size,
                marker=marker,
                color=colour,
                lw=0,
                label=f"{label}, rho {'above' if marker == '^' else 'below'} zero",
            )
    ax.scatter(k[cacng8], -np.log10(p[cacng8]), s=40, color=RED, lw=0, label="Cacng8")
    ax.set_xscale("log")
    ax.set_xlabel("genes, smallest p first")
    ax.set_ylabel("-log10 spatial p")
    ax.legend(loc="lower left", fontsize=7, frameon=False)
    tidy(ax)


# the reliability below which a gene's Allen experiments are said to disagree, as
# the gene table's figures say it (02)
def named_rows(tests: pd.DataFrame) -> dict[str, pd.Series]:
    """The rows of named_tests.csv by what they are.

    Tier 1, and tier 2's group tests and check rows on the leftover and on the map.
    """
    by = tests.set_index(["map", "test"])
    return dict(
        first=tests[tests["tier"] == "1"].iloc[0],
        spatial=by.loc[("leftover", SPATIAL_TEST)],
        matched=by.loc[("leftover", MATCHED_TEST)],
        without=by.loc[("leftover", WITHOUT_TEST)],
        outside=by.loc[("leftover", OUTSIDE_TEST)],
        on_map=by.loc[("nano", SPATIAL_TEST)],
        on_map_matched=by.loc[("nano", MATCHED_TEST)],
    )


def leftover_takeaway(
    t: dict[str, pd.Series], genes: pd.DataFrame, q: float, n_surrogates: int
) -> str:
    """The line under figure 11's title, each verdict following its p.

    `genes` is leftover_genes.csv: the genes past BH over all of them, and how many
    of those run below zero.
    """
    past = genes[genes["q_all"] < q]
    n_below = int((past["rho"] < 0).sum())
    follows = "follows" if t["first"]["p"] < ALPHA else "does not follow"
    if t["spatial"]["p"] < ALPHA and t["matched"]["p"] >= ALPHA:
        group = (
            "the family does so beyond the surrogates, not beyond other postsynaptic "
            "genes"
        )
    elif t["spatial"]["p"] < ALPHA:
        group = "the family does so beyond the surrogates and other postsynaptic genes"
    else:
        group = "the family as a group does not"
    return (
        f"Cacng8 {follows} what Gria1 and synapse density leave "
        f"({p_text(t['first']['p'], n_surrogates)}), a re-test: its p against two "
        f"earlier leftovers was seen before it was named; {group}; exploratory, "
        f"{len(past)} of {len(genes)} genes pass BH, {n_below} of them below zero"
    )


def tier_one_panel(
    ax: plt.Axes,
    profile: dict[str, float],
    leftover: pd.Series,
    set_table: pd.DataFrame,
    first: pd.Series,
    n_surrogates: int,
) -> None:
    """A of figure 11: Cacng8 against the leftover, the structures furthest off named."""
    x, y, shared = against_gene(
        ax,
        profile,
        leftover,
        group_of(set_table),
        ("Cacng8, rank among structures", "leftover (ranks above prediction)"),
        ranks=False,
    )
    order = np.argsort(y)
    acronyms = acronym_of(set_table)
    named = list(order[-3:]) + list(order[:2])
    spread_labels(ax, [(x[i], y[i], acronyms.get(shared[i], ""), "0.2") for i in named])
    panel_title(
        ax,
        "A",
        "Tier 1: Cacng8",
        f"rho {first['first']:+.2f}, {p_text(first['p'], n_surrogates)}",
    )


def plot_leftover(
    table: pd.DataFrame,
    leftover: pd.Series,
    profile: dict[str, float],
    tests: pd.DataFrame,
    nulls: dict[str, dict[str, np.ndarray]],
    genes: pd.DataFrame,
    set_table: pd.DataFrame,
    q: float,
    n_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 11: the three tiers named in advance against the leftover.

    Cacng8, the AMPA receptor complex family and every other gene, against what Gria1 and
    synapse density leave.

    `table` is top_genes.csv, `leftover` the main model's leftover on the structures
    of its fit, `profile` Cacng8's, `tests` named_tests.csv, `nulls` per map the
    surrogates' medians of the family, `genes` leftover_genes.csv; `n_surrogates`
    counts the leftover's surrogates.
    """
    t = named_rows(tests)
    family = table[table["family"]].reset_index(drop=True)
    past = genes[genes["q_all"] < q]
    n_pass, n_below = len(past), int((past["rho"] < 0).sum())
    fig = plt.figure(figsize=(17, 6.8))
    heading(fig, "leftover", leftover_takeaway(t, genes, q, n_surrogates))

    # A, tier 1: Cacng8 against the leftover
    ax = fig.add_axes([0.04, 0.21, 0.19, 0.55])
    tier_one_panel(ax, profile, leftover, set_table, t["first"], n_surrogates)

    # B and C, tier 2: the family as a group, against the surrogates and against
    # postsynaptic genes of the same expression
    ax = fig.add_axes([0.29, 0.21, 0.19, 0.55])
    group_null_panel(ax, t["spatial"], nulls["leftover"]["spatial"], RED)
    panel_title(
        ax,
        "B",
        f"Tier 2: the family ({len(family)} genes)",
        f"median {t['spatial']['first']:+.3f} (red), "
        f"{p_text(t['spatial']['p'], n_surrogates)}; without Cacng8\n"
        f"{p_text(t['without']['p'], n_surrogates)}",
    )
    ax = fig.add_axes([0.54, 0.21, 0.19, 0.55])
    paired_strips(
        ax, family, "leftover_rho", "control_leftover_rho", np.random.default_rng(0)
    )
    ax.set_ylabel("Spearman rho with the leftover")
    if t["outside"]["controls_changed"] == 0:
        outside = "no control is a model term, so the check row is the same test"
    else:
        outside = (
            f"controls outside the model's terms {t['outside']['difference']:+.3f}, "
            f"p = {t['outside']['p']:.2f}"
        )
    panel_title(
        ax,
        "C",
        "The family against matched genes",
        f"{t['matched']['difference']:+.3f} between medians, p = "
        f"{t['matched']['p']:.2f};\n{outside}",
    )

    # D, tier 3: every gene
    ax = fig.add_axes([0.79, 0.21, 0.19, 0.55])
    bh_panel(ax, genes, list(family["symbol"]), q)
    panel_title(
        ax,
        "D",
        "Tier 3: every gene, exploratory",
        f"{n_pass} of {len(genes)} cross the line, {n_below} of them below zero",
    )
    footer(
        fig,
        [
            "How to read: the leftover's surrogates each go through the same fit; tier 1 "
            "is one test, uncorrected; tier 2 is the family as a group, against the "
            "surrogates and against expression-matched postsynaptic genes, each one "
            "test.",
            "A: dots coloured by group of divisions (cortex orange, hippocampal "
            "formation red, thalamus blue, other grey matter).",
            f"D: with {n_surrogates:,} surrogates the smallest p is "
            f"{1 / (n_surrogates + 1):.4f}, just under the line for the first gene "
            f"({q} / {len(genes)}); the model's own genes (Gria1, the density genes) "
            "stay in the BH, as the rule set it.",
            f"Every gene and gene set: {figure_ref('leftover_genes')}; the family member "
            f"by member and on the map: {figure_ref('ampa_family')}.",
        ],
    )
    return saved(fig, save)


def family_takeaway(
    t: dict[str, pd.Series], past: pd.DataFrame, n_surrogates: int
) -> str:
    """The line under figure 11s2's title, each verdict following its p.

    `past` holds the members past BH within the family, whose sign the line says.
    """
    cacng8 = "follows" if t["first"]["p"] < ALPHA else "does not follow"
    beyond = "follows it beyond" if t["spatial"]["p"] < ALPHA else "does not pass"
    than = "more" if t["matched"]["p"] < ALPHA else "no more"
    within = "no member passes BH within it"
    if len(past):
        within = f"past BH within it: {', '.join(past['symbol'])}"
        if (past["leftover_rho"] < 0).all():
            within += ", each below zero"
        elif (past["leftover_rho"] > 0).all():
            within += ", each above zero"
    return (
        f"Cacng8 {cacng8} the leftover ({p_text(t['first']['p'], n_surrogates)}, a "
        f"re-test: its p against two earlier leftovers was seen before it was named); "
        f"the family {beyond} the surrogates "
        f"(p = {t['spatial']['p']:.3f}), {than} than postsynaptic genes of its "
        f"expression (p = {t['matched']['p']:.3f}); {within}"
    )


def on_map_panels(
    fig: plt.Figure,
    rows: pd.DataFrame,
    t: dict[str, pd.Series],
    nulls: dict[str, dict[str, np.ndarray]],
    rng: np.random.Generator,
    n_map_surrogates: int,
) -> None:
    """D and E of figure 11s2: the family's group tests on the map itself.

    A description of the family rather than a test of the leftover.
    """
    ax = fig.add_axes([0.5, 0.13, 0.2, 0.3])
    group_null_panel(ax, t["on_map"], nulls["nano"]["spatial"], NANO)
    panel_title(
        ax,
        "D",
        "On the map itself (a description)",
        f"median {t['on_map']['first']:+.3f} (orange); "
        f"{p_text(t['on_map']['p'], n_map_surrogates)}",
    )
    ax = fig.add_axes([0.77, 0.13, 0.2, 0.3])
    paired_strips(ax, rows, "rho", "control_rho", rng)
    ax.set_ylabel("Spearman rho with the nano map")
    panel_title(
        ax,
        "E",
        "On the map, against the same controls",
        f"{t['on_map_matched']['difference']:+.3f} between medians, p = "
        f"{t['on_map_matched']['p']:.3f}",
    )


def plot_ampa_family(
    rows: pd.DataFrame,
    tests: pd.DataFrame,
    nulls: dict[str, dict[str, np.ndarray]],
    members: pd.DataFrame,
    q: float,
    t_max: float,
    n_surrogates: int,
    n_map_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 11s2: the AMPA receptor complex against what Gria1 and density leave.

    `rows` holds the family's rows of top_genes.csv, `tests` named_tests.csv, `nulls`
    per map the surrogates' medians and the label null of the group tests, and
    `members` family_members.csv; `n_surrogates` and `n_map_surrogates` count the
    leftover's and the map's surrogates.
    """
    rows = rows.sort_values("leftover_rho", ascending=False).reset_index(drop=True)
    t = named_rows(tests)
    passing = rows[rows["leftover_q_family"] < q]
    past = passing["symbol"]
    untested = members.loc[~members["tested"], "symbol"]
    in_model = ", ".join(rows.loc[rows["in_model"] != "", "symbol"]) or "none"
    shaky = rows["reliability"] < LOW_RELIABILITY
    unreliable = ", ".join(rows.loc[shaky, "symbol"]) or "none"
    fig = plt.figure(figsize=(16, 11.5))
    heading(fig, "ampa_family", family_takeaway(t, passing, n_surrogates))

    # A: every member against the leftover
    ax = fig.add_axes([0.14, 0.12, 0.27, 0.74])
    family_bars(ax, rows, q, t_max, n_surrogates)
    panel_title(
        ax,
        "A",
        "Each member against the leftover",
        f"pale: 95% of {n_surrogates} surrogates, each through the same fit; bold: "
        f"past BH\nwithin the family ({len(past)} of {len(rows)}); not tested: "
        f"{', '.join(untested)}",
    )

    # B and C: the group against the surrogates, and against matched controls
    ax = fig.add_axes([0.5, 0.58, 0.2, 0.27])
    group_null_panel(ax, t["spatial"], nulls["leftover"]["spatial"], RED)
    panel_title(
        ax,
        "B",
        "As a group, against the surrogates",
        f"median {t['spatial']['first']:+.3f} (red); "
        f"{p_text(t['spatial']['p'], n_surrogates)}",
    )
    rng = np.random.default_rng(0)
    ax = fig.add_axes([0.77, 0.58, 0.2, 0.27])
    paired_strips(ax, rows, "leftover_rho", "control_leftover_rho", rng)
    ax.set_ylabel("Spearman rho with the leftover")
    panel_title(
        ax,
        "C",
        "Against genes of the same expression",
        f"{t['matched']['difference']:+.3f} between medians, p = "
        f"{t['matched']['p']:.3f}\n(labels permuted; ±{t['matched']['critical']:.2f} "
        "would give p < 0.05)",
    )

    # D and E: the same on the map itself, to describe the family
    on_map_panels(fig, rows, t, nulls, rng, n_map_surrogates)
    footer(
        fig,
        [
            "How to read: the family was named before this model was fitted, from "
            "sources outside this analysis: S, the native AMPA receptor complexes "
            "of Schwenk et al. 2012 (Figure 1D, Table S2); GO, GO:0032281 AMPA",
            "glutamate receptor complex; P, the partner subunits (dark blue); minus "
            "Gria1, the abundance term. Tier 1 is Cacng8 alone (red); tier 2 the family "
            "as a group (B, C), then gene by gene, BH within the family (A, bold).",
            f"In the model's own terms, so near zero against the leftover by "
            f"construction: {in_model}. Allen experiments that disagree (reliability "
            f"below {LOW_RELIABILITY}): {unreliable}.",
        ],
    )
    return saved(fig, save)


# ===== Top-gene sheets =====


def numbers_block(row: pd.Series) -> list[str]:
    """The lines of a gene's sheet that its panels do not show."""
    reliability = row["reliability"]
    lines = [
        f"tier against the leftover: {row['tier']}",
        f"Allen: {int(row['n_experiments'])} experiments, reliability "
        + ("-" if not np.isfinite(reliability) else f"{reliability:.2f}"),
        f"with the map: rank {int(row['rank_all'])} of all genes, q {row['q_all']:.3f}",
        f"over the robustness variants: rank {int(row['rank_variants_min'])} to "
        f"{int(row['rank_variants_max'])}, rho {row['rho_variants_min']:+.2f} to "
        f"{row['rho_variants_max']:+.2f}",
        f"inside divisions: {row['rho_within']:+.2f}, q {row['q_within']:.3f}",
        f"with the leftover: rank {int(row['leftover_rank'])}, q over every gene "
        f"{row['leftover_q_all']:.3f}",
        f"the main model on its {int(row['n_added'])} structures leaves "
        f"{row['left_main']:.1%}, with the gene {row['left_with_gene']:.1%}",
    ]
    if np.isfinite(row["rho_psd95"]):
        lines.append(
            f"PSD95 puncta, {int(row['n_psd95'])} structures: rho {row['rho_psd95']:+.2f}"
        )
    if row["in_model"]:
        lines.append(f"in the main model: {row['in_model']}")
    if row["family"]:
        lines += [
            f"q within the family: {row['leftover_q_family']:.3f}",
            f"family sources: {row['family_sources']}",
            f"matched control: {row['matched_control']}",
        ]
    return lines


def annotation_lines(row: pd.Series, width: int) -> list[str]:
    """What the gene table says the gene is, wrapped to `width` characters."""
    out = []
    for label, column in (
        ("role in the ontology panel", "ontology_role"),
        ("GO terms that put it there", "go_panel_terms"),
        ("GO terms at the synapse", "go_synapse_terms"),
        ("gene sets", "gene_sets"),
        ("P9's category (a label)", "p9_category"),
    ):
        text = row[column] if isinstance(row[column], str) and row[column] else "-"
        out += textwrap.wrap(f"{label}: {text}", width, subsequent_indent="    ")
    return out


def sheet_heading(
    fig: plt.Figure, row: pd.Series, n_surrogates: int, n_map_surrogates: int
) -> None:
    """A gene sheet's two top lines: who the gene is, then its numbers.

    Why it is here; then its numbers with the map, inside divisions, with the leftover and
    what it takes.
    """
    symbol = row["symbol"]
    why = []
    if row["top"]:
        why.append("past the map's null after BH")
    if symbol == LEFTOVER_GENE:
        why.append("the one gene named for the leftover")
    elif row["named"]:
        why.append("the gene of the stained protein")
    if row["family"]:
        why.append("a member of the AMPA receptor complex family")
    fig.text(
        0.5,
        0.975,
        f"{symbol}, {row['name']}: {'; '.join(why)}",
        ha="center",
        va="top",
        fontsize=12,
        color=named_colour(symbol),
    )
    if row["in_model"] == "abundance":
        taken = "a term of the main model"
    else:
        taken = (
            f"takes {points(row['taken'])} points of the reproducible map (p = "
            f"{row['p_taken']:.3f} plain, {row['p_taken_alike']:.3f} alike)"
        )
    fig.text(
        0.5,
        0.945,
        f"with the map {row['rho']:+.2f} "
        f"({p_text(row['p_spatial'], n_map_surrogates)}, rank "
        f"{int(row['rank_all'])}); inside divisions {row['rho_within']:+.2f}; with the "
        f"leftover {row['leftover_rho']:+.2f} "
        f"({p_text(row['leftover_p'], n_surrogates)}); {taken}",
        ha="center",
        va="top",
        fontsize=9,
        color=DARK_GREY,
    )


def taken_panel(ax: plt.Axes, row: pd.Series, taken_nulls: dict[str, np.ndarray]) -> None:
    """E of a gene sheet: what the gene takes of the leftover, against its nulls.

    Maps of its smoothness put in its place: plain surrogates, and maps alike to the
    model.
    """
    own = row["in_model"] == "abundance"
    plain, alike = 100 * taken_nulls["plain"], 100 * taken_nulls["alike"]
    bins = np.linspace(min(plain.min(), alike.min()), max(plain.max(), alike.max()), 41)
    ax.hist(plain, bins=bins, color=NULL_BAND, label="plain surrogates of the gene")

    # a term of the model has no remainder, so maps alike to the model are the gene
    if not own:
        ax.hist(
            alike,
            bins=bins,
            histtype="step",
            color=MID_GREY,
            lw=1.1,
            label="maps alike to the model",
        )
    ax.axvline(100 * row["taken"], color=RED, lw=1.8, label="the gene")
    ax.set_xlabel("points of the reproducible map taken from the leftover")
    ax.set_ylabel(f"maps of its smoothness ({len(plain)} each)")
    ax.legend(
        loc="upper right", fontsize=7, frameon=True, framealpha=0.9, edgecolor="none"
    )
    tidy(ax)
    if own:
        numbers = "a term of the main model: adding it again takes nothing"
    else:
        numbers = (
            f"{points(row['taken'])} points; p = {row['p_taken']:.3f} against plain "
            f"surrogates, {row['p_taken_alike']:.3f} against maps alike"
        )
    panel_title(ax, "E", "Added to the main model", numbers)


def likeness_panel(ax: plt.Axes, row: pd.Series) -> None:
    """F of a gene sheet: the gene's rho with the map and with what predicts it.

    With the map, inside divisions, with each term of the main model and its prediction,
    and with the leftover.
    """
    items = [("the map", row["rho"]), ("inside divisions", row["rho_within"])]
    items += [(label, row[column]) for column, label, _, _ in LIKENESS if column in row]
    items.append(("the leftover", row["leftover_rho"]))
    ax.barh(range(len(items)), [float(v) for _, v in items], color=MID_GREY, height=0.6)
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([label for label, _ in items], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_xlim(-1, 1)
    ax.set_xlabel("Spearman rho")
    tidy(ax)
    panel_title(ax, "F", "What it is like")


def plot_top_gene_sheet(
    row: pd.Series,
    map_values: pd.Series,
    profile: dict[str, float],
    leftover: pd.Series,
    taken_nulls: dict[str, np.ndarray],
    set_table: pd.DataFrame,
    lab: np.ndarray,
    names: dict[int, str],
    plane: int,
    n_surrogates: int,
    n_map_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """One gene characterised on one sheet.

    Its map, with the map and the leftover, what it takes, what it is like, and what the
    gene table says it is.

    `row` is the gene's row of top_genes.csv, `leftover` the main model's leftover
    on the structures of its fit, `taken_nulls` what maps of the gene's smoothness
    take in its place (plain, alike); `n_surrogates` and `n_map_surrogates` count
    the leftover's and the map's surrogates. PNG only, for review.
    """
    symbol = row["symbol"]
    groups = group_of(set_table)
    fig = plt.figure(figsize=(16, 10))
    sheet_heading(fig, row, n_surrogates, n_map_surrogates)

    # A: the gene's map as ranks on the plane
    shared = [s for s in map_values.index if np.isfinite(profile.get(s, np.nan))]
    ranks = ranks01(np.array([profile[s] for s in shared]))
    ax = fig.add_axes([0.01, 0.5, 0.27, 0.36])
    rank_plane(fig, ax, dict(zip(shared, ranks)), lab, names)
    panel_title(
        ax,
        "A",
        f"{symbol} as ranks, CCF plane {plane}",
        f"{len(shared)} declared structures; flat grey: others",
    )

    # B and C: with the map and with the leftover
    ax = fig.add_axes([0.37, 0.53, 0.25, 0.33])
    against_gene(ax, profile, map_values, groups, (f"{symbol}, rank", "nano map, rank"))
    panel_title(
        ax, "B", "With the map", f"rho {row['rho']:+.2f}, {len(shared)} structures"
    )
    ax = fig.add_axes([0.71, 0.53, 0.25, 0.33])
    against_gene(
        ax, profile, leftover, groups, (f"{symbol}, rank", "leftover (ranks)"), False
    )
    panel_title(
        ax,
        "C",
        "With what Gria1 and density leave",
        f"rho {row['leftover_rho']:+.2f}; {p_text(row['leftover_p'], n_surrogates)}",
    )

    # D: the numbers, and what the gene table says it is; E and F beside them
    ax = fig.add_axes([0.01, 0.03, 0.3, 0.4])
    ax.axis("off")
    lines = numbers_block(row) + [""] + annotation_lines(row, 64)
    ax.text(0, 1, "\n".join(lines), fontsize=7.5, va="top", ha="left", linespacing=1.4)
    panel_title(ax, "D", "Numbers, and what the gene table says")
    taken_panel(fig.add_axes([0.37, 0.08, 0.25, 0.3]), row, taken_nulls)
    likeness_panel(fig.add_axes([0.71, 0.08, 0.25, 0.3]), row)
    return saved(fig, save, eps=False)
