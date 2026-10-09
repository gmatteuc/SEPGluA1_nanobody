"""The overview of the ISH line, and April's headline then and now.

Figure 00 puts the line on one page: the question, the two parts of the argument
with their key numbers and what each stands on, the limit, and the map of the
figures. Figure 15 and 15s show what is left of April's headline: the gene order,
the category p under each choice it rests on, and the groups against the null.

Called by run_ish_overview.py.
"""

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle
from scipy.stats import spearmanr

from sepmap.ish.figure_index import FIGURE_NUMBERS, QUESTIONS, figure_ref
from sepmap.ish.plotting.shared import footer, heading, p_text, panel_title, saved
from sepmap.plotting import (
    DARK_BLUE,
    DARK_GREY,
    LIGHT_GREY,
    MID_GREY,
    NANO,
    NANO_DOT,
    NULL_BAND,
    PAIR_LINE,
    RED,
    tidy,
)

# ===== 00 The ISH line on one page =====

# the figure's width and margins, in inches; the height follows from the text
OVERVIEW_WIDTH = 16.0
OVERVIEW_MARGIN = 0.45

# the tint of each block, a palette colour and its alpha, and the colour of its key
# numbers
PART_TINTS = {
    "question": (LIGHT_GREY, 0.35),
    "part 1": (NANO, 0.22),
    "part 2": (NULL_BAND, 0.75),
    "limit": (LIGHT_GREY, 0.35),
}
KEY_NUMBER_COLOURS = {"part 1": NANO_DOT, "part 2": DARK_BLUE}

# the font sizes of the overview: a part's claim, its key numbers, the text beside
# them, where it stands, and a line of the figure map
CLAIM_SIZE = 10.5
KEY_SIZE = 15.0
TEXT_SIZE = 9.0
STANDS_SIZE = 8.8
MAP_SIZE = 8.3

# the groups of the figure map in each of its three columns, by their place in the walk
MAP_COLUMNS = ((0, 1), (2,), (3, 4, 5))


def wrapped(text: str, width_in: float, size: float, bold: bool = False) -> list[str]:
    """`text` wrapped to a column `width_in` inches wide at font `size` points.

    A character of the default sans font is on average a little over half the
    font size wide, bold a little more, which keeps a line inside its column with
    a little to spare.
    """
    share = 0.6 if bold else 0.54
    chars = max(int(width_in / (share * size / 72)), 10)
    lines = []
    for paragraph in text.split("\n"):
        lines += textwrap.wrap(paragraph, chars) or [""]
    return lines


def line_height(size: float) -> float:
    """The height of one line of text at font `size`, in inches."""
    return 1.45 * size / 72


def text_height(text: str, width_in: float, size: float, bold: bool = False) -> float:
    """The height `text` takes wrapped to `width_in` inches, in inches."""
    return len(wrapped(text, width_in, size, bold)) * line_height(size)


def text_block(
    ax: plt.Axes, x: float, y: float, width_in: float, text: str, size: float, **style
) -> float:
    """`text` wrapped at (x, y) inches from the top left; returns the height it takes."""
    bold = style.get("fontweight") == "bold"
    lines = wrapped(text, width_in, size, bold)
    ax.text(x, y, "\n".join(lines), fontsize=size, va="top", zorder=3, **style)
    return len(lines) * line_height(size)


def tinted_box(ax: plt.Axes, x: float, y: float, w: float, h: float, part: str) -> None:
    """A rounded box in the tint of its part, in the overview's inch coordinates."""
    colour, alpha = PART_TINTS[part]
    style = "round,pad=0,rounding_size=0.06"
    ax.add_patch(
        FancyBboxPatch((x, y), w, h, boxstyle=style, fc=colour, alpha=alpha, ec="none")
    )
    ax.add_patch(
        FancyBboxPatch((x, y), w, h, boxstyle=style, fc="none", ec=MID_GREY, lw=0.6)
    )


def part_height(part: dict, width: float) -> float:
    """The height of one part's box: its claim, its key numbers, where it stands."""
    h = 0.15 + text_height(part["claim"], width - 0.3, CLAIM_SIZE, bold=True) + 0.12
    for _, text in part["numbers"]:
        beside = text_height(text, width - 2.05, TEXT_SIZE)
        h += max(line_height(KEY_SIZE), beside) + 0.12
    return h + text_height(part["stands"], width - 0.3, STANDS_SIZE) + 0.2


def draw_part(
    ax: plt.Axes, x: float, y: float, width: float, height: float, part: dict
) -> None:
    """One part of the argument in its tinted box, the key numbers large on the left."""
    tinted_box(ax, x, y, width, height, part["part"])
    yy = y + 0.15
    yy += text_block(
        ax, x + 0.15, yy, width - 0.3, part["claim"], CLAIM_SIZE, fontweight="bold"
    )
    yy += 0.12
    for key, text in part["numbers"]:
        ax.text(
            x + 0.15,
            yy,
            key,
            fontsize=KEY_SIZE,
            fontweight="bold",
            color=KEY_NUMBER_COLOURS[part["part"]],
            va="top",
            zorder=3,
        )
        beside = text_block(ax, x + 1.9, yy + 0.03, width - 2.05, text, TEXT_SIZE)
        yy += max(line_height(KEY_SIZE), beside) + 0.12
    text_block(
        ax,
        x + 0.15,
        yy,
        width - 0.3,
        part["stands"],
        STANDS_SIZE,
        color=DARK_GREY,
        fontstyle="italic",
    )


def map_line(label: str, question: str, details: str) -> str:
    """One figure of the map: its number, its question, its detailed versions."""
    text = f"{label}  {question}"
    if details:
        text += f"  ({details})"
    return text


def map_column_height(groups: list, width: float) -> float:
    """The height of one column of the figure map."""
    h = 0.0
    for _, figures in groups:
        h += line_height(MAP_SIZE + 0.7) + 0.05
        for figure in figures:
            h += text_height(map_line(*figure), width, MAP_SIZE)
        h += 0.15
    return h


def draw_figure_map(ax: plt.Axes, top: float, figure_map: list) -> float:
    """The figures in three columns, a group per part of the walk.

    Returns the inch below it.

    `figure_map` holds (heading, [(label, question, detailed labels)]) per group.
    """
    ax.text(
        OVERVIEW_MARGIN,
        top,
        "The figures (a main figure, and in brackets its detailed versions)",
        fontsize=10,
        fontweight="bold",
        va="top",
    )
    top += 0.35
    width = (OVERVIEW_WIDTH - 2 * OVERVIEW_MARGIN) / len(MAP_COLUMNS)
    bottom = top
    for c, places in enumerate(MAP_COLUMNS):
        x = OVERVIEW_MARGIN + c * width
        y = top
        for heading_text, figures in (figure_map[i] for i in places):
            ax.text(
                x, y, heading_text, fontsize=MAP_SIZE + 0.7, fontweight="bold", va="top"
            )
            y += line_height(MAP_SIZE + 0.7) + 0.05
            for figure in figures:
                y += text_block(ax, x, y, width - 0.25, map_line(*figure), MAP_SIZE)
            y += 0.15
        bottom = max(bottom, y)
    return bottom


def labelled_box(
    ax: plt.Axes, y: float, h: float, part: str, label: str, text: str, size: float
) -> None:
    """A box the overview's width: its label in bold on the left, its text beside it."""
    width = OVERVIEW_WIDTH - 2 * OVERVIEW_MARGIN
    tinted_box(ax, OVERVIEW_MARGIN, y, width, h, part)
    ax.text(
        OVERVIEW_MARGIN + 0.15,
        y + 0.12,
        label,
        fontsize=10.5,
        fontweight="bold",
        va="top",
    )
    text_block(ax, OVERVIEW_MARGIN + 1.75, y + 0.12, width - 1.9, text, size)


def plot_overview(
    content: dict, figure_map: list, save: Path | None = None
) -> plt.Figure:
    """Figure 00: the question, the two parts of the argument, the limit, the figures.

    Each part with its key numbers; the map of the figures last.

    `content` holds the question, the two parts (part, claim, numbers as (key
    number, text), stands) and the limit (ish.overview.overview_content);
    `figure_map` the groups of figures (ish.overview.figure_map). The figure is as
    tall as its text: every size is worked out in inches first.
    """
    full = OVERVIEW_WIDTH - 2 * OVERVIEW_MARGIN
    half = (full - 0.2) / 2
    question_h = text_height(content["question"], full - 1.9, 10) + 0.25
    parts_h = max(part_height(part, half) for part in content["parts"])
    limit_h = text_height(content["limit"], full - 1.9, 9.5) + 0.25
    width = full / len(MAP_COLUMNS)
    columns = [[figure_map[i] for i in places] for places in MAP_COLUMNS]
    map_h = 0.35 + max(map_column_height(groups, width - 0.25) for groups in columns)
    height = 1.0 + question_h + 0.15 + parts_h + 0.15 + limit_h + 0.3 + map_h + 0.6

    fig = plt.figure(figsize=(OVERVIEW_WIDTH, height))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, OVERVIEW_WIDTH)
    ax.set_ylim(height, 0)
    ax.axis("off")
    ax.text(
        OVERVIEW_WIDTH / 2,
        0.3,
        f"{FIGURE_NUMBERS['overview']}.  {QUESTIONS['overview']}",
        ha="center",
        va="top",
        fontsize=13,
    )
    ax.text(
        OVERVIEW_WIDTH / 2,
        0.62,
        "every number is this run's (tables/numbers_for_the_text.csv); the story, with "
        "what each result means and does not, is docs/ISH_ANALYSIS.md",
        ha="center",
        va="top",
        fontsize=9,
        color=DARK_GREY,
    )

    # the question, the two parts side by side, the limit
    y = 1.0
    labelled_box(ax, y, question_h, "question", "The question", content["question"], 10)
    y += question_h + 0.15
    for k, part in enumerate(content["parts"]):
        draw_part(ax, OVERVIEW_MARGIN + k * (half + 0.2), y, half, parts_h, part)
    y += parts_h + 0.15
    labelled_box(ax, y, limit_h, "limit", "The limit", content["limit"], 9.5)
    y += limit_h + 0.3

    # the map of the figures
    draw_figure_map(ax, y, figure_map)
    return saved(fig, save)


# ===== 16 and 16s April's headline, then and now =====

# the genes named on the violins of the appendix
HEADLINE_NAMED = (
    "Cacng8",
    "Grm5",
    "Gria1",
    "Gria2",
    "Gria3",
    "Gria4",
    "Dlg2",
    "Nptx1",
    "Htr3a",
    "Chrm1",
    "Slc17a6",
    "Th",
    "Glra1",
    "Calb2",
    "Olig2",
    "Arc",
    "Cnih3",
    "Hcn1",
    "Nrgn",
)

# P9's kernel for its violins: Gaussian, bandwidth 0.12 in rho
VIOLIN_BANDWIDTH = 0.12


def half_violin(ax: plt.Axes, x0: float, values: np.ndarray, colour, side: int) -> None:
    """Half a violin at `x0`, to the left (side -1) or the right (+1), P9's kernel."""
    if len(values) < 2:
        return
    h = VIOLIN_BANDWIDTH
    grid = np.linspace(values.min() - 2 * h, values.max() + 2 * h, 200)
    density = np.exp(-0.5 * ((grid[:, None] - values[None, :]) / h) ** 2).sum(axis=1)
    density = density / density.max() * 0.17
    ax.fill_betweenx(grid, x0, x0 + side * density, color=colour, alpha=0.28, lw=0)


def headline_violins(ax: plt.Axes, headline: pd.DataFrame, order: tuple) -> None:
    """A: each group's April violin (grey, left) beside today's (orange, right)."""
    ax.axhline(0, color="0.8", lw=0.7, zorder=0)
    for i, (key, _) in enumerate(order):
        mine = headline[headline["group"] == key]
        april = mine.dropna(subset=["rho_april"])
        today = mine.dropna(subset=["rho_today"])
        xa, xt = i - 0.2, i + 0.2
        half_violin(ax, xa, april["rho_april"].to_numpy(), MID_GREY, -1)
        half_violin(ax, xt, today["rho_today"].to_numpy(), NANO, +1)
        both = mine.dropna(subset=["rho_april", "rho_today"])
        for _, r in both.iterrows():
            ax.plot([xa, xt], [r["rho_april"], r["rho_today"]], color=PAIR_LINE, lw=0.5)
        ax.scatter([xa] * len(april), april["rho_april"], s=16, color=DARK_GREY, lw=0)
        ax.scatter([xt] * len(today), today["rho_today"], s=16, color=NANO_DOT, lw=0)
        for x, values, colour in (
            (xa, april["rho_april"], DARK_GREY),
            (xt, today["rho_today"], NANO_DOT),
        ):
            if len(values):
                m = float(np.median(values))
                ax.plot([x - 0.13, x + 0.13], [m, m], color=colour, lw=2.4)
        placed = []
        for value, symbol in sorted(
            zip(today["rho_today"], today["symbol"]), key=lambda t: t[0]
        ):
            if symbol not in HEADLINE_NAMED:
                continue
            y = max(value, placed[-1] + 0.05) if placed else value
            placed.append(y)
            colour = DARK_BLUE if symbol.startswith("Gria") else "0.2"
            ax.text(xt + 0.05, y, symbol, fontsize=6.5, va="center", color=colour)
        ax.text(
            i,
            -0.95,
            f"n = {len(april)} / {len(today)}",
            ha="center",
            fontsize=7,
            color=DARK_GREY,
        )
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([label for _, label in order], fontsize=8.5)
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylim(-1.0, 1.05)
    ax.set_ylabel("Spearman rho of the gene with the adult map")
    ax.legend(
        handles=[
            plt.Line2D(
                [],
                [],
                ls="",
                marker="o",
                mfc=DARK_GREY,
                mec="none",
                label="April (22 April, eroded Spearman, P9's structures)",
            ),
            plt.Line2D(
                [],
                [],
                ls="",
                marker="o",
                mfc=NANO_DOT,
                mec="none",
                label="today (the declared structures, QC, merged profiles)",
            ),
            plt.Line2D([], [], color=PAIR_LINE, lw=1, label="the same gene"),
            plt.Line2D([], [], color="0.3", lw=2.4, label="the group's median"),
        ],
        loc="lower center",
        ncol=4,
        fontsize=7.5,
        bbox_to_anchor=(0.5, 1.0),
    )
    tidy(ax)


def headline_scatter(ax: plt.Axes, headline: pd.DataFrame, letter: str = "B") -> float:
    """B: each gene's rho in April against today's; returns their agreement."""
    both = headline.dropna(subset=["rho_april", "rho_today"])
    x = both["rho_april"].to_numpy()
    y = both["rho_today"].to_numpy()
    lim = (-0.75, 0.95)
    ax.plot(lim, lim, color=MID_GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.axhline(0, color="0.88", lw=0.6, zorder=0)
    ax.axvline(0, color="0.88", lw=0.6, zorder=0)
    sub = both["symbol"].str.startswith("Gria").to_numpy()
    ax.scatter(x[~sub], y[~sub], s=28, color=DARK_GREY, lw=0, alpha=0.85, zorder=2)
    ax.scatter(x[sub], y[sub], s=34, color=DARK_BLUE, lw=0, zorder=3)
    moved = set(both["symbol"].iloc[np.argsort(-np.abs(y - x))[:5]])
    for symbol in moved | {"Cacng8", "Gria1", "Gria4", "Dlg2"}:
        k = list(both["symbol"]).index(symbol)
        colour = DARK_BLUE if symbol.startswith("Gria") else "0.15"
        ax.annotate(
            symbol,
            (x[k], y[k]),
            xytext=(4, 2),
            textcoords="offset points",
            fontsize=7,
            color=colour,
        )
    agreement = spearmanr(x, y).statistic
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("April (eroded Spearman)")
    ax.set_ylabel("today")
    panel_title(
        ax,
        letter,
        f"Gene by gene ({len(both)} genes in both)",
        f"the order of the genes agrees at rho {agreement:.2f}; subunits dark blue",
    )
    tidy(ax)
    return agreement


def anova_panel(
    ax: plt.Axes, anova: pd.DataFrame, f_info: dict, letter: str = "C"
) -> None:
    """C: the ANOVA p under each choice, and the same F against the surrogates."""
    rows = list(anova.itertuples(index=False))
    labels = [f"{r.label} ({r.n_genes})" for r in rows]
    labels.append("today, F against the map's surrogates\n(co-expressed genes kept)")
    ps = [r.p for r in rows] + [f_info["p_spatial"]]
    ys = np.arange(len(ps))[::-1]
    for k, (y, p) in enumerate(zip(ys, ps)):
        ax.axhline(y, color="0.93", lw=0.8, zorder=0)
        if k == len(ps) - 1:
            ax.scatter([p], [y], s=60, marker="D", color="0.1", zorder=3, lw=0)
        else:
            colour = DARK_GREY if rows[k].kind == "april" else NANO_DOT
            ax.scatter([p], [y], s=44, color=colour, zorder=3, lw=0)
        ax.text(p * 1.3, y, f"{p:.3f}", va="center", fontsize=7.5, color="0.2")
    ax.axvline(0.05, color=RED, lw=1, zorder=0)
    ax.text(0.05, len(ps) - 0.4, " 0.05", color=RED, fontsize=7.5, va="bottom")
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1.5)
    ax.set_ylim(-0.7, len(ps) - 0.2)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("one-way ANOVA p across the groups (log scale)")
    panel_title(
        ax,
        letter,
        "The group p under each choice it rests on",
        "dots treat co-expressed genes as independent draws, as April did; the "
        "diamond does not",
    )
    tidy(ax)


def headline_null_panel(
    ax: plt.Axes, headline: pd.DataFrame, groups: pd.DataFrame, n_surrogates: int
) -> None:
    """D: today's groups, each median against the band of its surrogates' medians."""
    rng = np.random.default_rng(0)
    ax.axhline(0, color="0.8", lw=0.7, zorder=0)
    labels = []
    for i, g in enumerate(groups.itertuples(index=False)):
        mine = headline[(headline["group"] == g.group) & headline["rho_today"].notna()]
        if g.tested:
            ax.add_patch(
                Rectangle(
                    (i - 0.4, g.null_lo),
                    0.8,
                    g.null_hi - g.null_lo,
                    color=NULL_BAND,
                    lw=0,
                    zorder=0,
                )
            )
        jitter = rng.uniform(-0.15, 0.15, len(mine))
        ax.scatter(i + jitter, mine["rho_today"], s=16, color=NANO_DOT, lw=0, zorder=2)
        ax.plot([i - 0.3, i + 0.3], [g.median_rho] * 2, color="0.1", lw=2, zorder=3)
        after = "\nsplit written\nafter looking" if g.group == "aux_forebrain" else ""
        if g.tested:
            labels.append(
                f"{g.label}\nn = {g.n_genes}\n{p_text(g.p_spatial, n_surrogates)}, "
                f"q = {g.q:.3f}{after}"
            )
        else:
            labels.append(f"{g.label}\nn = {g.n_genes}\nnot tested, no band{after}")
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.set_ylabel("today's rho with the adult map")
    tidy(ax)


def plot_april_headline(
    headline: pd.DataFrame,
    anova: pd.DataFrame,
    groups: pd.DataFrame,
    f_info: dict,
    n_surrogates: int,
    q: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 16: what is left of April's headline.

    The gene order, the category p under each choice it rests on, and today's groups
    against the null.

    The arguments are those of plot_april_headline_detail but the groups' order.
    """
    first = anova.iloc[0]
    tested = groups[groups["tested"]]
    n_past = int((tested["q"] < q).sum())
    fig = plt.figure(figsize=(16, 11))
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1, 0.9],
        width_ratios=[1, 1.25],
        hspace=0.5,
        wspace=0.45,
        left=0.06,
        right=0.98,
        top=0.86,
        bottom=0.13,
    )
    agreement = headline_scatter(fig.add_subplot(grid[0, 0]), headline, letter="A")
    anova_panel(fig.add_subplot(grid[0, 1]), anova, f_info, letter="B")
    ax = fig.add_subplot(grid[1, :])
    headline_null_panel(ax, headline, groups, n_surrogates)
    panel_title(
        ax,
        "C",
        "Today's groups against the null",
        f"{n_past} of {len(tested)} tested groups past it after BH; pale: 95% of each "
        "group's median over the surrogates",
    )
    heading(
        fig,
        "april_headline",
        f"April's gene order reproduces (rho {agreement:.2f}); its category p of "
        f"{first['p']:.3f} rests on a split written after looking and on genes "
        f"treated as independent: against the surrogates "
        f"{p_text(f_info['p_spatial'], n_surrogates)}",
    )
    footer(
        fig,
        [
            "How to read: the groups are P9's categories of gene_targets.csv, with "
            "'auxiliary' split by hand into Aux forebrain and Aux other; the surrogates "
            "keep co-expressed genes together, since every gene meets the same "
            "surrogate of the map.",
            f"April's violins beside today's: {figure_ref('april_headline_detail')}.",
        ],
    )
    return saved(fig, save)


def plot_april_headline_detail(
    headline: pd.DataFrame,
    anova: pd.DataFrame,
    groups: pd.DataFrame,
    f_info: dict,
    order: tuple,
    n_surrogates: int,
    q: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 16s: April's category violins, recomputed today and against the null.

    `headline` holds P9's genes with their group and both rho (ish.overview's
    april_headline.csv), `anova` the ANOVA under each choice, `groups` each group's
    median and null band today, `f_info` the ANOVA's F against the surrogates,
    `order` the groups as April sorted them, (key, label), `q` the BH level.
    """
    fig = plt.figure(figsize=(16, 17))
    grid = fig.add_gridspec(
        3,
        2,
        height_ratios=[1.25, 1, 0.9],
        width_ratios=[1, 1.25],
        hspace=0.42,
        wspace=0.45,
        left=0.06,
        right=0.98,
        top=0.9,
        bottom=0.12,
    )
    ax = fig.add_subplot(grid[0, :])
    headline_violins(ax, headline, order)
    ax.set_title(
        "A.  April's ten groups: its headline (grey, left of each pair) beside the "
        "same genes today (orange, right); April's order",
        loc="left",
        fontsize=9,
        pad=30,
    )
    agreement = headline_scatter(fig.add_subplot(grid[1, 0]), headline)
    anova_panel(fig.add_subplot(grid[1, 1]), anova, f_info)
    ax = fig.add_subplot(grid[2, :])
    headline_null_panel(ax, headline, groups, n_surrogates)
    tested = groups[groups["tested"]]
    n_past = int((tested["q"] < q).sum())
    panel_title(
        ax,
        "D",
        "Today's groups against the null: each median against where the median of "
        "the same genes falls with 95% of the map's surrogates",
        f"{n_past} of {len(tested)} tested groups past the null at BH q < {q}; the "
        f"groups' F against the surrogates: {p_text(f_info['p_spatial'], n_surrogates)}",
    )
    first = anova.iloc[0]
    plain = anova.iloc[1]
    heading(
        fig,
        "april_headline_detail",
        f"April: {first['n_genes']} genes, group ANOVA p {first['p']:.3f} with the hand "
        f"split, {plain['p']:.2f} without; the gene order reproduces today at rho "
        f"{agreement:.2f}; against the surrogates "
        f"{p_text(f_info['p_spatial'], n_surrogates)}",
    )
    footer(
        fig,
        [
            "How to read: the groups are P9's categories of gene_targets.csv, written "
            "while the genes were chosen, with 'auxiliary' split by hand into Aux "
            "forebrain (Cacng8, Cacng3, Cnih2, Cnih3, Grm5) and Aux other.",
            "An ANOVA over genes treats co-expressed genes as independent draws; the "
            "surrogates keep them together, since every gene meets the same surrogate "
            "of the map.",
            "What would mean what: groups past their bands and an F past its null: "
            "the kinds of genes differ beyond the brain's smoothness, though a group "
            "written after looking (Aux forebrain) cannot count as a test.",
            "An F inside its null: what reproduces from April is the order of the "
            "genes, not a difference between the groups; the question is asked again "
            f"on sets fixed in advance ({figure_ref('gene_sets')}).",
        ],
    )
    return saved(fig, save)
