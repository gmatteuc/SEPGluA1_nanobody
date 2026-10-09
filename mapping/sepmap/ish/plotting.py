"""The figures of the ISH analysis: the numbered guide, the QC sheets, the gene sheets.

The guided figures are numbered in the order of the argument (03_beyond.png,
04_beyond_where.png, ...). A main figure answers one question in two to four
panels: its title is the number and the question, under it one line that says what
to take from it, computed from the tables, and at its foot a line or two in grey on
how to read it. Its detailed version, a supplementary figure, carries the same
number with an s (03s_beyond_budget.png; 10s1 and 10s2 where a main figure has
two), with every panel and number of the analysis. A takeaway that is a verdict
follows its numbers, so a rerun cannot leave a figure saying what its numbers
contradict.

The drawing follows docs/STYLE.md: maps of intensities and ranks in hot on black,
the atlas dark grey beneath and the area borders on top, structures left out flat
grey; scatters of structures with 35-point dots coloured by group of divisions;
bars whose grey says how reliable a value is; the subunits in dark blue. The
modules that compute draw nothing: every figure function takes tables and returns
the figure, saved as PNG and EPS at 150 dpi when `save` is given.

Called by the run scripts of the ISH line.
"""

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
from scipy.stats import rankdata, spearmanr

from sepmap.ish.gene_sets import AMPA_FAMILY, LEFTOVER_GENE, detectable
from sepmap.ish.section_qc import ISH_QC, SECTION_AXIS
from sepmap.plotting import (
    AUTO,
    AUTO_DOT,
    DARK_BLUE,
    DARK_GREY,
    DENSITY_BLUE,
    DIVISION_GROUP,
    DIVISION_GROUP_COLOURS,
    GROUP_COLOURS,
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

# the guided figures in the order of the argument: the question (00), the inputs (01,
# 02), part 1, how much of the map Gria1 and synapse density leave (03, 04), part 2,
# the genes against the map and against the leftover (05 to 11), the controls (12 to
# 14), the limit (15) and April's headline (16). A main figure's detailed version
# carries its number with an s. A figure's title, its file name and the figures'
# references to each other read this table; the docstrings and the run scripts'
# headers that name the files follow it by hand
FIGURES = {
    "overview": "00",
    "structures": "01",
    "genes": "02",
    "genes_detail": "02s",
    "beyond": "03",
    "beyond_budget": "03s",
    "beyond_where": "04",
    "one_comparison": "05",
    "one_comparison_detail": "05s",
    "spatial_null": "06",
    "spatial_null_detail": "06s",
    "top_genes": "07",
    "gene_ranking": "07s",
    "cacng8_gria1": "08",
    "between_within": "09",
    "between_within_detail": "09s",
    "gene_kinds": "10",
    "gene_sets": "10s1",
    "localisation": "10s2",
    "leftover": "11",
    "leftover_genes": "11s1",
    "ampa_family": "11s2",
    "autofluorescence": "12",
    "autofluorescence_detail": "12s",
    "robustness": "13",
    "robustness_detail": "13s",
    "synaptome": "14",
    "synaptome_detail": "14s",
    "green_channel": "15",
    "green_channel_detail": "15s",
    "april_headline": "16",
    "april_headline_detail": "16s",
}

# the question each guided figure answers: its title, and its heading in the index of
# the figures (ish.overview)
QUESTIONS = {
    "overview": "The ISH line on one page",
    "structures": "Which structures enter every comparison?",
    "genes": "Which genes, and how good is one Allen map?",
    "genes_detail": "The genes in detail: experiments, section QC, reliability, what "
    "was left out",
    "beyond": "Do Gria1 expression and synapse density explain the map?",
    "beyond_budget": "Part 1 in detail: each predictor, the budget, the calibration, "
    "the check rows",
    "beyond_where": "Where does the map sit above or below what they predict?",
    "one_comparison": "What is a gene's rho with the map?",
    "one_comparison_detail": "One comparison step by step: the maps as measured, as "
    "ranks, for three genes",
    "spatial_null": "How large a rho do unrelated smooth maps give?",
    "spatial_null_detail": "The spatial null in detail: surrogates on a plane, two "
    "genes against their null",
    "top_genes": "Which genes follow the map, and what kind of maps are they?",
    "gene_ranking": "P9's 100 genes against the map, the Cacng8 - Gria1 gap in "
    "detail, the two genes with the tissue",
    "cacng8_gria1": "Does the map follow Cacng8 more closely than Gria1?",
    "between_within": "Do genes follow the map inside divisions, or only between them?",
    "between_within_detail": "Between and within divisions in detail: five genes "
    "division by division, the choice of null",
    "gene_kinds": "Do the genes that set surface receptor follow the map better than "
    "other genes?",
    "gene_sets": "The gene sets in detail: every set, the two contrasts named in "
    "advance, where each set comes from",
    "localisation": "The localisation test in detail: its null, the positive control, "
    "the matching, every test of the design",
    "leftover": "Does any gene follow what Gria1 and synapse density leave?",
    "leftover_genes": "Every gene and every gene set against the leftover",
    "ampa_family": "The AMPA receptor complex family against the leftover, member by "
    "member, and on the map",
    "autofluorescence": "Is the gene ranking the label's or the tissue's?",
    "autofluorescence_detail": "The autofluorescence map in detail: the genes' rho with "
    "each map, the counts past each null",
    "robustness": "Does the ranking depend on the choices made?",
    "robustness_detail": "The choices made, with every gene under four of them",
    "synaptome": "Does a measured synapse density change part 1?",
    "synaptome_detail": "The measured synapse density in detail: coverage, the two "
    "hemispheres, agreement with every map",
    "green_channel": "Does the green channel report the tagged receptor, or the tissue?",
    "green_channel_detail": "The three channels in detail: one adult's raw planes, each "
    "channel against Gria1",
    "april_headline": "What is left of April's headline?",
    "april_headline_detail": "April's ten groups then and now, and against the null",
}


# ===== Shared pieces =====


def figure_file(key: str) -> str:
    """The file name of a guided figure: its number and its key, 07s_gene_ranking.png."""
    return f"{FIGURES[key]}_{key}.png"


def figure_ref(key: str) -> str:
    """How a figure names another: 'figure 06', 'figure 10s2'."""
    return f"figure {FIGURES[key]}"


def heading(fig: plt.Figure, key: str, line: str) -> None:
    """The figure's number and question, and under it one line computed from the
    tables: what to take from a main figure, the numbers of a detailed one.

    Placed in inches from the top, so a short main figure and a tall detailed one
    leave the same room between the two lines.
    """
    height = fig.get_figheight()
    title = f"{FIGURES[key]}.  {QUESTIONS[key]}"
    fig.text(0.5, 1 - 0.16 / height, title, ha="center", va="top", fontsize=12)
    fig.text(
        0.5,
        1 - 0.47 / height,
        line,
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


def saved(fig: plt.Figure, path: Path | None, eps: bool = True) -> plt.Figure:
    """Save the figure as PNG, and with `eps` as EPS, when a path is given; return it.

    The QC sheets are PNG only: they are for review, not for a paper, and 750 of
    them as EPS would take about 2 GB.
    """
    if path is not None:
        save_figure(fig, Path(path), DPI, eps=eps)
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
            (n_few, "white", "////"),
            (n_white, "#e6e6e6", None),
            (n_catch, "white", "...."),
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
        Patch(color=c, label=f"kept: {g}") for g, c in DIVISION_GROUP_COLOURS.items()
    ]
    handles += [
        Patch(
            facecolor="white", edgecolor=MID_GREY, hatch="////", label=f"left out: {few}"
        ),
        Patch(color="#e6e6e6", label="left out: not grey matter"),
        Patch(
            facecolor="white",
            edgecolor=MID_GREY,
            hatch="....",
            label="left out: a catch-all label",
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
        "structures",
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


# ===== QC sheets =====

# the colours of a section's status in the QC figures: not judged white, ok light
# grey, flagged red; a dim section kept as a true absence is hatched over light grey,
# and a section dim against one side only (a step in expression, kept) is marked
# with a dark grey dot
STATUS_CODES = {"ok": 1, "flagged": 2, "absence kept": 3, "step": 4}
STATUS_COLOURS = ["white", "#e6e6e6", RED]

# the reason section_qc gives a section kept at a step in expression
STEP_REASON = "dim against one side only"


def status_matrix(sections: pd.DataFrame, experiments: list[str]) -> np.ndarray:
    """Experiments by sections, each cell a status code (0 for not judged)."""
    n_sections = int(sections["section"].max()) + 1
    out = np.zeros((len(experiments), n_sections), dtype=int)
    row_of = {e: i for i, e in enumerate(experiments)}
    reasons = sections["reason"].fillna("").astype(str)
    for eid, k, status, reason in zip(
        sections["experiment_id"], sections["section"], sections["status"], reasons
    ):
        if eid in row_of:
            code = STATUS_CODES.get(status, 0)
            if status == "ok" and reason.startswith(STEP_REASON):
                code = STATUS_CODES["step"]
            out[row_of[eid], int(k)] = code
    return out


def draw_status(ax: plt.Axes, matrix: np.ndarray) -> None:
    """Draw a status matrix: white, light grey, red, hatching for kept absences, a dot
    on a section kept at a step."""
    shown = np.where(matrix == 3, 1, matrix)
    shown = np.where(matrix == 4, 1, shown)
    ax.imshow(
        shown,
        cmap=ListedColormap(STATUS_COLOURS),
        vmin=-0.5,
        vmax=2.5,
        aspect="auto",
        interpolation="nearest",
    )
    for i, k in zip(*np.nonzero(matrix == 3)):
        ax.add_patch(
            Rectangle(
                (k - 0.5, i - 0.5), 1, 1, fill=False, hatch="////", edgecolor=RED, lw=0
            )
        )
    rows, cols = np.nonzero(matrix == 4)
    ax.scatter(cols, rows, s=9, color=DARK_GREY, lw=0, zorder=3)


def status_handles() -> list:
    """The legend of a status matrix."""
    return [
        Patch(
            facecolor="white",
            edgecolor=MID_GREY,
            lw=0.5,
            label="not judged: too "
            "little brain, no data, or neighbours at the noise level",
        ),
        Patch(color="#e6e6e6", label="judged and kept"),
        Patch(color=RED, label="dim against both sides: set missing"),
        Patch(
            facecolor="#e6e6e6",
            edgecolor=RED,
            hatch="////",
            lw=0,
            label="dim, kept as a true absence (proposed, to review)",
        ),
        plt.Line2D(
            [],
            [],
            ls="",
            marker="o",
            ms=4,
            color=DARK_GREY,
            label="dim against one side only: a step in expression, kept",
        ),
    ]


def plot_flagged(
    sections: pd.DataFrame,
    summary: pd.DataFrame,
    p9_genes: set[str],
    save: Path | None = None,
) -> plt.Figure:
    """Every experiment with a flagged section or a kept absence, one row each.

    The sheet to review the exceptions list from: sections as columns along each
    experiment's own axis, flagged sections red, dim sections kept as true absence
    hatched; P9's genes in bold.
    """
    ok = summary[summary["grid"] == "ok"]
    reasons = sections["reason"].fillna("").astype(str)
    stepped = set(sections.loc[reasons.str.startswith(STEP_REASON), "experiment_id"])
    keep = (ok["n_flagged"] > 0) | (ok["n_absence_kept"] > 0)
    shown = ok[keep | ok["experiment_id"].isin(stepped)]
    shown = shown.sort_values(["symbol", "experiment_id"])
    experiments = list(shown["experiment_id"])
    matrix = status_matrix(sections, experiments)
    height = 2.9 + 0.15 * len(experiments)
    fig = plt.figure(figsize=(11, height))
    ax = fig.add_axes([0.2, 1.7 / height, 0.75, 1 - 2.6 / height])
    draw_status(ax, matrix)
    labels = [
        f"{s}  {e} ({p[0]})"
        for s, e, p in zip(shown["symbol"], shown["experiment_id"], shown["plane"])
    ]
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=6.5)
    for tick, symbol in zip(ax.get_yticklabels(), shown["symbol"]):
        if symbol in p9_genes:
            tick.set_fontweight("bold")
    ax.set_xlabel(
        "section along the experiment's own axis, 200 um (AP for coronal, c; "
        "ML for sagittal, s)"
    )
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    n_sections = int(shown["n_flagged"].sum())
    n_kept = int(shown["n_absence_kept"].sum())
    n_steps = int((matrix == STATUS_CODES["step"]).sum())
    fig.text(
        0.5,
        1 - 0.25 / height,
        f"Section QC: {len(shown)} of {len(ok)} experiments with a dim section; "
        f"{n_sections} sections set missing in {int((shown['n_flagged'] > 0).sum())} "
        f"experiments (red), {n_kept} kept as true absence (hatched), {n_steps} kept "
        "at a step (dot); P9's genes in bold",
        ha="center",
        va="top",
        fontsize=10,
    )
    fig.legend(
        handles=status_handles(),
        loc="lower center",
        ncol=2,
        fontsize=7.5,
        bbox_to_anchor=(0.5, 0.15 / height),
        frameon=False,
    )
    return saved(fig, save, eps=False)


def middle_section(table: pd.DataFrame) -> int:
    """The judged, unflagged section nearest the middle of the judged ones."""
    good = table.loc[table["status"] == "ok", "section"].to_numpy()
    if good.size == 0:
        return int(table["section"].iloc[len(table) // 2])
    middle = (good.min() + good.max()) / 2
    return int(good[np.argmin(np.abs(good - middle))])


def section_plane(volume: np.ndarray, k: int, axis: int) -> np.ndarray:
    """One section of an (AP, DV, ML) volume, dorsal up: (DV, ML) or (DV, AP)."""
    if axis == 0:
        return volume[k]
    return volume[:, :, k].T


def profile_panel(ax: plt.Axes, table: pd.DataFrame, shown: int, fraction: float):
    """A: each section's median energy, the flag line, flagged and kept sections."""
    k = table["section"].to_numpy()
    median = table["median_energy"].to_numpy()
    status = table["status"].to_numpy()
    colours = [RED if s == "flagged" else MID_GREY for s in status]
    ax.bar(k, np.nan_to_num(median), color=colours, width=0.85)
    top = np.nanmax(median) if np.isfinite(median).any() else 1.0
    if not top > 0:
        top = 1.0

    # a flagged or kept section is often near zero, so it is also marked over the
    # whole height: a red triangle on the axis, or a hatched band
    for kk, s in zip(k, status):
        if s == "flagged":
            ax.plot([kk], [0], marker="^", color=RED, ms=6, clip_on=False)
        elif s == "absence kept":
            ax.add_patch(
                Rectangle(
                    (kk - 0.45, 0),
                    0.9,
                    top * 1.12,
                    fill=False,
                    hatch="////",
                    edgecolor=RED,
                    lw=0,
                )
            )
    judged = np.isin(status, ["ok", "flagged", "absence kept"])
    line = np.where(judged, fraction * table["local_reference"].to_numpy(), np.nan)
    ax.step(k, line, where="mid", color=DARK_GREY, ls="--", lw=0.9)
    ax.plot([shown], [top * 1.05], marker="v", color=DARK_GREY, ms=6)
    ax.set_xlim(k.min() - 1, k.max() + 1)
    ax.set_ylim(0, top * 1.12)
    ax.set_xlabel(f"section along {table['axis'].iloc[0]}, 200 um")
    ax.set_ylabel("median energy in the brain")
    tidy(ax)


def plot_qc_sheet(
    table: pd.DataFrame,
    summary: dict,
    grid: np.ndarray,
    template: np.ndarray,
    labels: np.ndarray,
    save: Path | None = None,
) -> plt.Figure:
    """One experiment's QC sheet: its section profile, and the orientation check.

    `table` holds its rows of section_qc.csv, `summary` its row of
    experiment_qc.csv, `grid` the energy (AP, DV, ML) cut to the CCF on the 200 um
    grid; `template` and `labels` are the CCF template and its structures (codes
    of structures.name_volume) on the 20 um grid, ten times finer, so that the
    borders drawn are those of structures and not of single 200 um voxels.
    """
    axis = SECTION_AXIS[summary["plane"]]
    shown = middle_section(table)
    fig = plt.figure(figsize=(12, 7.4))
    flagged = summary["flagged_sections"] or "none"
    kept = summary["absence_sections"] or "none"
    fig.text(
        0.5,
        0.975,
        f"{summary['symbol']}, Allen experiment {summary['experiment_id']} "
        f"({summary['plane']}): sections set missing: {flagged}; kept as true "
        f"absence: {kept}",
        ha="center",
        va="top",
        fontsize=11,
    )

    # A: the profile along the section axis
    ax = fig.add_axes([0.07, 0.57, 0.9, 0.3])
    profile_panel(ax, table, shown, ISH_QC["local_fraction"])
    panel_title(
        ax,
        "A",
        "Median energy of each section; dashed: "
        f"{ISH_QC['local_fraction']} x the median of the {ISH_QC['neighbours']} "
        "judged sections on each side",
        f"{summary['n_judged']} sections judged; red: below the dashed line and below "
        f"{ISH_QC['local_fraction']} x the brightest of the {ISH_QC['neighbours']} on "
        "each side, set missing; hatched: dim but kept as a true absence; empty: no "
        "data; grey triangle: the section below",
    )

    # B to D: the orientation check at one section; a 200 um section k is 20 um
    # plane 10 k, and each of its voxels covers 10 x 10 pixels of the finer grid
    fine = 10 * shown
    atlas = section_plane(template, fine, axis).astype(float)
    lab = section_plane(labels, fine, axis)
    energy = section_plane(grid, shown, axis)
    energy_fine = np.repeat(np.repeat(energy, 10, axis=0), 10, axis=1)
    finite = grid[np.isfinite(grid)]
    vmax = float(np.percentile(finite, 99)) if finite.size else 1.0
    ax_b = fig.add_axes([0.03, 0.07, 0.28, 0.36])
    ax_c = fig.add_axes([0.35, 0.07, 0.28, 0.36])
    ax_d = fig.add_axes([0.67, 0.07, 0.28, 0.36])
    ax_b.imshow(atlas, cmap="gray", interpolation="nearest")
    ax_c.imshow(
        np.ma.masked_invalid(energy),
        cmap=hot_cut(),
        vmin=0,
        vmax=vmax,
        interpolation="nearest",
    )
    ax_c.set_facecolor("k")
    image = draw_plane(ax_d, energy_fine, lab, hot_cut(), 0, vmax)
    colour_bar(fig, ax_d, image, "expression energy")
    for a in (ax_b, ax_c):
        a.set_xticks([])
        a.set_yticks([])
    if axis == 0:
        what = "coronal, dorsal up"
    else:
        what = "sagittal, dorsal up, anterior left"
    ax_b.set_title(f"B.  CCF template at section {shown}", loc="left", fontsize=9)
    ax_b.set_xlabel(what, fontsize=8)
    ax_c.set_title("C.  The ISH section, as on its grid", loc="left", fontsize=9)
    ax_d.set_title("D.  The same, with the CCF's structures", loc="left", fontsize=9)
    return saved(fig, save, eps=False)


# ===== 02 The genes, and how good their maps are =====

# where a gene comes from, in the order panel A stacks them
SOURCES = ("P9's panel only", "both panels", "ontology panel only")

# a gene in a GO set and a marker set is drawn in the marker set
STACK_ORDER = (
    "subunits",
    "localisation",
    "GABAergic markers",
    "glia",
    "other postsynaptic",
    "presynaptic",
)

# the genes named on the reliability panels: the subunit, the top of P9's ranking,
# and the scaffold second to it
NAMED_GENES = ("Gria1", "Cacng8", "Dlg2")


def gene_source(genes: pd.DataFrame) -> pd.Series:
    """Which panel lists each gene: P9's only, both, or the ontology panel's only."""
    in_p9 = genes["p9_gene"].astype(bool)
    in_ontology = genes["ontology_role"].fillna("") != ""
    source = np.where(in_p9 & in_ontology, SOURCES[1], SOURCES[2])
    source = np.where(in_p9 & ~in_ontology, SOURCES[0], source)
    return pd.Series(source, index=genes.index)


def first_set(text: str) -> str:
    """The set a gene is drawn in: the first of STACK_ORDER it belongs to."""
    sets = [s.strip() for s in str(text).split(";") if s.strip()]
    for name in STACK_ORDER:
        if name in sets:
            return name
    return "in no set"


def union_panel(ax: plt.Axes, genes: pd.DataFrame) -> None:
    """A: the two panels and their union, each part stacked by gene set."""
    source = gene_source(genes)
    drawn = genes["gene_sets"].fillna("").map(first_set)
    colours = dict(SET_COLOURS)
    colours["in no set"] = "#ececec"
    for i, part in enumerate(SOURCES):
        left = 0
        for name in STACK_ORDER + ("in no set",):
            n = int(((source == part) & (drawn == name)).sum())
            if n:
                ax.barh(i, n, left=left, color=colours[name], height=0.6)
            left += n
        ax.text(left + 4, i, f"{left}", va="center", fontsize=8)
    ax.set_yticks(range(len(SOURCES)))
    ax.set_yticklabels(SOURCES)
    ax.invert_yaxis()
    ax.set_xlabel("genes with a usable experiment")
    handles = [Patch(color=colours[n], label=n) for n in STACK_ORDER + ("in no set",)]
    ax.legend(handles=handles, loc="upper right", fontsize=7)
    panel_title(
        ax,
        "A",
        "The genes, by panel and by gene set",
        f"{len(genes)} genes; one in a GO set and a marker set drawn in the marker set",
    )
    tidy(ax)


def experiments_panel(ax: plt.Axes, experiments: pd.DataFrame) -> None:
    """B: usable experiments per gene, by the planes of section they cover."""
    used = experiments[~experiments["excluded"]]
    per_gene = used.groupby("symbol")["plane"].agg(lambda p: " ".join(sorted(set(p))))
    count = used.groupby("symbol")["experiment_id"].count().clip(upper=4)
    kinds = {
        "coronal only": DARK_GREY,
        "sagittal only": "#c8c8c8",
        "both planes": MID_GREY,
    }
    bottom = np.zeros(4)
    for kind, colour in kinds.items():
        if kind == "both planes":
            mine = per_gene == "coronal sagittal"
        else:
            mine = per_gene == kind.split()[0]
        n = np.array([int(((count == k) & mine).sum()) for k in (1, 2, 3, 4)])
        ax.bar([1, 2, 3, 4], n, bottom=bottom, color=colour, width=0.7, label=kind)
        bottom += n
    for k, n in zip((1, 2, 3, 4), bottom):
        ax.text(k, n + 3, f"{int(n)}", ha="center", fontsize=8)
    ax.set_ylim(0, bottom.max() * 1.12)
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(["1", "2", "3", "4 or more"])
    ax.set_xlabel("usable Allen experiments of the gene")
    ax.set_ylabel("genes")
    ax.legend(loc="upper right", fontsize=7)
    n_two = int((count >= 2).sum())
    panel_title(
        ax,
        "B",
        "Experiments per gene",
        f"{len(used)} experiments used, {n_two} genes measured more than once",
    )
    tidy(ax)


def qc_panel(
    ax: plt.Axes,
    sections: pd.DataFrame,
    summary: pd.DataFrame,
    p9_experiments: set[str],
) -> None:
    """C: P9's own experiments with a dim section, one row each, named."""
    ok = summary[summary["grid"] == "ok"]
    dim = ok[(ok["n_flagged"] > 0) | (ok["n_absence_kept"] > 0)]
    shown = dim[dim["experiment_id"].isin(p9_experiments)].sort_values("symbol")
    matrix = status_matrix(sections, list(shown["experiment_id"]))
    draw_status(ax, matrix)
    ax.set_yticks(range(len(shown)))
    ax.set_yticklabels(shown["symbol"], fontsize=8)
    ax.set_xlabel("section along AP, 200 um (P9's experiments are coronal)")
    with_grid = len(p9_experiments & set(ok["experiment_id"]))
    panel_title(
        ax,
        "C",
        "P9's experiments with a dim section",
        f"{len(shown)} of the {with_grid} P9 experiments with a usable grid (F lists "
        f"the others);\n{len(dim) - len(shown)} more among the other experiments "
        "(qc/00_flagged.png)",
    )


def reliability_panel(ax: plt.Axes, rel: pd.DataFrame, letter: str = "D") -> None:
    """D: how reliable one Allen map is, the named genes marked."""
    values = rel["reliability"].dropna()
    ax.hist(values, bins=np.linspace(-0.4, 1.0, 36), color=MID_GREY)
    ax.axvline(0.3, color=DARK_GREY, ls="--", lw=0.9)
    top = ax.get_ylim()[1]
    by = rel.set_index("symbol")["reliability"]
    named = [(g, by.get(g, np.nan)) for g in NAMED_GENES]
    named = [(g, v) for g, v in named if np.isfinite(v)]
    groups: list[list[tuple[str, float]]] = []
    for gene, value in sorted(named, key=lambda t: t[1]):
        if groups and abs(value - groups[-1][-1][1]) < 0.02:
            groups[-1].append((gene, value))
        else:
            groups.append([(gene, value)])
    ax.set_ylim(0, top * 1.2)
    for i, group in enumerate(groups):
        value = float(np.mean([v for _, v in group]))
        ax.axvline(value, color="0.1", lw=1.0, zorder=3)
        ax.text(
            value - 0.01,
            top * (1.12 - 0.08 * i),
            f"{', '.join(g for g, _ in group)}: {value:.2f}",
            fontsize=7.5,
            ha="right",
            va="center",
        )
    q1, median, q3 = np.percentile(values, [25, 50, 75])
    ax.set_xlabel("reliability: median Spearman between two experiments of a gene")
    ax.set_ylabel("genes")
    panel_title(
        ax,
        letter,
        "How reliable one Allen map is",
        f"{len(values)} genes; median {median:.2f} (quartiles {q1:.2f} to {q3:.2f}); "
        f"{int((values < 0.3).sum())} below 0.3 (dashed)",
    )
    tidy(ax)


def expression_panel(ax: plt.Axes, rel: pd.DataFrame, subunits: set[str]) -> None:
    """E: reliability against how strongly the gene is expressed."""
    have = rel.dropna(subset=["reliability"])
    have = have[have["median_energy"] > 0]
    x = np.log10(have["median_energy"].to_numpy())
    y = have["reliability"].to_numpy()
    sub = have["symbol"].isin(subunits).to_numpy()
    ax.scatter(x[~sub], y[~sub], s=20, color=DARK_GREY, alpha=0.6, linewidths=0)
    ax.scatter(x[sub], y[sub], s=46, color=DARK_BLUE, linewidths=0, zorder=4)
    for xi, yi, symbol in zip(x[sub], y[sub], have["symbol"][sub]):
        if symbol not in NAMED_GENES:
            ax.annotate(
                symbol,
                (xi, yi),
                xytext=(5, -9),
                textcoords="offset points",
                fontsize=7.5,
                color=DARK_BLUE,
            )
    offsets = {"Gria1": (6, -3), "Cacng8": (-40, 8), "Dlg2": (6, 8)}
    for gene in NAMED_GENES:
        mine = (have["symbol"] == gene).to_numpy()
        if mine.any():
            ax.annotate(
                gene,
                (x[mine][0], y[mine][0]),
                textcoords="offset points",
                xytext=offsets[gene],
                fontsize=7.5,
                arrowprops=dict(arrowstyle="-", color=MID_GREY, lw=0.5),
            )
    rho = pd.Series(x).corr(pd.Series(y), method="spearman")
    ax.axhline(0.3, color=DARK_GREY, ls="--", lw=0.9)
    ax.set_xlabel("log10 median expression energy of the gene")
    ax.set_ylabel("reliability")
    panel_title(
        ax,
        "E",
        "Reliability against expression",
        f"Spearman rho {rho:+.2f} over {len(have)} genes; subunits dark blue",
    )
    tidy(ax)


def reason_kind(reason: str) -> str:
    """A short kind of reason for an experiment left out, to group them by."""
    if reason.startswith("grid not on disk"):
        return "no grid from Allen"
    if reason.startswith("grid is"):
        return "grid in a box of its own"
    if reason.startswith("data in"):
        return "data in too few structures"
    return reason


def repair_panel(ax: plt.Axes, experiments: pd.DataFrame) -> None:
    """F: the experiments left out, by reason, and what the repair added."""
    ax.axis("off")
    dropped = experiments[experiments["excluded"]].sort_values("symbol")
    lines = [
        f"Experiments left out ({len(dropped)}; * P9's own;",
        " (+) another experiment of the gene is used):",
    ]
    kinds = dropped["exclude_reason"].map(reason_kind)
    used_genes = set(experiments.loc[~experiments["excluded"], "symbol"])
    for kind in dict.fromkeys(kinds):
        mine = dropped[kinds == kind]
        genes = [
            f"{g}{'*' if p9 else ''}{' (+)' if g in used_genes else ''}"
            for g, p9 in zip(mine["symbol"], mine["p9_experiment"])
        ]
        text = f"{kind}: {', '.join(genes)}"
        lines.extend(
            textwrap.wrap(text, 58, initial_indent="  ", subsequent_indent="    ")
        )
    lines.append("")
    lines.append("Repaired with the gene's other Allen experiments:")
    repair = experiments[experiments["repair_experiment"].astype(bool)]
    for symbol, mine in repair.groupby("symbol"):
        used = int((~mine["excluded"]).sum())
        planes = mine["plane"].value_counts()
        text = ", ".join(f"{n} {plane}" for plane, n in planes.items())
        lines.append(f"  {symbol}: {used} of {len(mine)} used ({text})")
    lines.append("  Sst: none needed, the ontology panel has two")
    text = "\n".join(lines)
    ax.text(0, 1, text, va="top", ha="left", fontsize=7.5, family="monospace")
    ax.set_title("F.  What was left out, and the repair", loc="left", fontsize=9)


def plot_genes(
    genes: pd.DataFrame,
    experiments: pd.DataFrame,
    rel: pd.DataFrame,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 02: the genes, their Allen experiments, and how reliable one map is.

    `genes` holds a row per gene with a usable experiment (p9_gene, ontology_role,
    gene_sets), `experiments` the gene table (excluded, plane), `rel` the
    reliability per gene.
    """
    fig = plt.figure(figsize=(16, 6.2))
    n_used = int((~experiments["excluded"]).sum())
    median = float(rel["reliability"].median())
    heading(
        fig,
        "genes",
        f"{len(genes)} genes and {n_used} usable Allen experiments; two experiments of "
        f"one gene agree at a median rho of {median:.2f}, which caps any gene's rho "
        "with the nano map",
    )
    union_panel(fig.add_axes([0.11, 0.17, 0.25, 0.62]), genes)
    experiments_panel(fig.add_axes([0.43, 0.17, 0.2, 0.62]), experiments)
    reliability_panel(fig.add_axes([0.7, 0.17, 0.27, 0.62]), rel, letter="C")
    footer(
        fig,
        [
            "How to read: one table holds P9's 100 genes and the ontology panel's 390, "
            "a row per Allen experiment; reliability is how well two experiments of "
            "the same gene, two Allen mice, agree across structures.",
        ],
    )
    return saved(fig, save)


def plot_genes_detail(
    genes: pd.DataFrame,
    experiments: pd.DataFrame,
    sections: pd.DataFrame,
    summary: pd.DataFrame,
    rel: pd.DataFrame,
    subunits: set[str],
    save: Path | None = None,
) -> plt.Figure:
    """Figure 02s: the genes, their experiments, the section QC and the reliability.

    `genes` holds a row per gene with a usable experiment (p9_gene, ontology_role,
    gene_sets), `experiments` the gene table (excluded, plane), `sections` and
    `summary` the section QC, `rel` the reliability per gene.
    """
    fig = plt.figure(figsize=(16, 10.5))
    n_listed = experiments["symbol"].nunique()
    n_used = int((~experiments["excluded"]).sum())
    heading(
        fig,
        "genes_detail",
        f"{n_listed} genes and {len(experiments)} experiments listed; {len(genes)} "
        f"genes and {n_used} experiments usable; "
        f"{int(experiments['excluded'].sum())} experiments left out (F)",
    )
    ax_a = fig.add_axes([0.11, 0.62, 0.25, 0.27])
    ax_b = fig.add_axes([0.43, 0.62, 0.2, 0.27])
    ax_c = fig.add_axes([0.72, 0.55, 0.26, 0.31])
    ax_d = fig.add_axes([0.06, 0.17, 0.26, 0.3])
    ax_e = fig.add_axes([0.40, 0.17, 0.24, 0.3])
    ax_f = fig.add_axes([0.70, 0.17, 0.28, 0.3])
    union_panel(ax_a, genes)
    experiments_panel(ax_b, experiments)
    p9_experiments = set(
        experiments.loc[experiments["p9_experiment"].astype(bool), "experiment_id"]
    )
    qc_panel(ax_c, sections, summary, p9_experiments)
    reliability_panel(ax_d, rel)
    expression_panel(ax_e, rel, subunits)
    repair_panel(ax_f, experiments)
    fig.legend(
        handles=status_handles()[:4],
        loc="upper left",
        bbox_to_anchor=(0.7, 0.2),
        fontsize=7,
        title="C, the sections",
        title_fontsize=7.5,
        frameon=False,
    )
    footer(
        fig,
        [
            "How to read: one table holds P9's 100 genes and the ontology panel's 390, "
            "a row per Allen experiment (A9). A section five-fold dimmer than the "
            "sections on both sides of it is set missing, never filled in (A2); a "
            "gene's profile is the mean rank of its usable experiments.",
            "Reliability is how well two experiments of the same gene, two Allen mice, "
            "agree across structures; a gene measured once has none.",
            "What would mean what: a gene's correlation with the nano map is capped by "
            "its own reliability, so a low rho of an unreliable gene says little, and "
            "genes measured once are only as good as one Allen mouse.",
        ],
    )
    return saved(fig, save)


# ===== Shared pieces of the analyses =====

# the structures named on the scatters of figure 05, by acronym: the top of the map,
# a thalamic relay nucleus, the striatum, the barrel field
NAMED_STRUCTURES = ("CA1", "VPM", "CP", "SSp-bfd")


def ordinal(value: str | float) -> str:
    """A rank as words write it: '1st', '2nd', '11th', '17th'."""
    k = int(float(value))
    suffix = "th"
    if k % 100 not in (11, 12, 13):
        suffix = {1: "st", 2: "nd", 3: "rd"}.get(k % 10, "th")
    return f"{k}{suffix}"


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


def scatter_groups(ax: plt.Axes, x: np.ndarray, y: np.ndarray, groups: list[str]):
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
    dot, in order of height; one that would sit on an earlier one (closer than a
    line of text, and within a name's width) is pushed down, and joined to its dot
    by a thin line when moved.
    """
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    gap = abs(y1 - y0) * 0.032 * fontsize / 7
    width = abs(x1 - x0) * 0.14
    placed: list[tuple[float, float]] = []
    for x, y, text, colour in sorted(points, key=lambda t: -t[1]):
        ty = y
        while any(abs(px - x) < width and abs(py - ty) < gap for px, py in placed):
            ty -= gap * 0.6
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


# ===== 05 What one comparison is =====


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
        f"spatial {p_text(r['p_spatial'], n_surrogates)}; rank {int(r['rank_p9'])} "
        f"of P9's {row['n_p9']} genes",
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
    """Figure 05: a gene's rho with the map, on a gene that follows it and one that
    does not.

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
    `n_p9`, the genes of P9's panel ranked.
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
            rate = float((sub[column] < 0.05).mean())
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


def false_positive_rates(calibration: pd.DataFrame, column: str) -> dict[str, float]:
    """{design: share of tests with p < 0.05} of one kind of p in the calibration."""
    rates = {}
    for design in dict.fromkeys(calibration["design"]):
        sub = calibration[calibration["design"] == design]
        rates[design] = float((sub[column] < 0.05).mean())
    return rates


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
    ordinary = false_positive_rates(calibration, "p_ordinary")
    spatial = false_positive_rates(calibration, "p_spatial")
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
    spatial = false_positive_rates(calibration, "p_spatial")
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
    past = (table["p_spatial"] < 0.05).to_numpy()
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
    """B: the Cacng8 - Gria1 gap against its two nulls, with its adult interval and
    the pairings of experiments as rugs beneath.

    `null` holds the null gaps of maps related to both genes alike (first row) and
    of maps unrelated to both (second row).
    """
    merged = gap[gap["kind"] == "merged profiles"].iloc[0]
    pairs = gap[gap["kind"] == "experiment pairing"]
    first, second = merged["first"], merged["second"]
    equal, unrelated = null[0], null[1]
    bins = np.linspace(-0.8, 0.8, 81)
    ax.axvspan(
        merged["equal_lo"],
        merged["equal_hi"],
        ymin=0.19,
        color=NULL_BAND,
        alpha=0.5,
        lw=0,
    )
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
        f"maps related to both alike: {p_text(merged['p_equal'], n_surrogates)} (shaded: "
        f"their 95%, {merged['equal_lo']:+.2f} to {merged['equal_hi']:+.2f}); "
        f"unrelated maps: {p_text(merged['p_spatial'], n_surrogates)}",
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
        ("p < 0.05\n(past its band)", lambda t: t["p_spatial"] < 0.05),
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
    """Figure 07s: P9's genes ranked against the map, with the null; the gap; the two
    genes with the tissue.

    `ranking` is gene_ranking.csv, `gap` gap.csv and `gap_null` the gap of every
    surrogate.
    """
    nano = ranking[(ranking["map"] == "nano") & ranking["p9_gene"]]
    nano = nano.sort_values("rho", ascending=False)
    auto = ranking[ranking["map"] == "auto"].set_index("symbol")["rho"]
    n_pass = int((nano["q_p9"] < q).sum())
    n_band = int((nano["p_spatial"] < 0.05).sum())
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
NAMED_AUTO = ("Cacng8", "Gria1", "Aqp4")

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
) -> None:
    """A: each gene's rho with autofluorescence against its rho with nano."""
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
    named += [g for g in NAMED_AUTO if g in genes and g not in named]
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
    panel_title(
        ax,
        letter,
        "Each gene against both maps",
        f"{len(genes)} genes (P9's dark, subunits blue);\nthe two gene orders agree "
        f"at rho {agree:+.2f}",
    )
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
    n_surrogates: int,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 12: whether Gria1 and Cacng8 follow the label or the tissue, adult by
    adult, against each map's null, and gene by gene.

    `ranking` is gene_ranking.csv (both maps), `per_adult` the per-adult rho table.
    """
    genes = ("Gria1", "Cacng8")
    mine = per_adult[per_adult["symbol"].isin(genes)]

    # the adults in which nano's rho is above autofluorescence's for both genes
    nano_above = mine["rho_nano"] > mine["rho_auto"]
    above = int(nano_above.groupby(mine["mouse"]).all().sum())
    auto = ranking[(ranking["map"] == "auto") & ranking["symbol"].isin(genes)]
    inside = "inside" if (auto["p_spatial"] >= 0.05).all() else "not all inside"
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
    auto_scatter(fig.add_axes([0.7, 0.13, 0.28, 0.64]), ranking, subunits, letter="C")
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


def variant_rows(summary: pd.DataFrame) -> np.ndarray:
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
    panel_title(ax, "B", "Where the two genes sit", "")
    tidy(ax)


def gap_rows_panel(
    ax: plt.Axes, summary: pd.DataFrame, y: np.ndarray, band: tuple[float, float]
) -> None:
    """C: the Cacng8 - Gria1 gap under each variant, the primary's null shaded."""
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
        f"pale blue: 95% of the primary gap's null\n(maps related to both alike, "
        f"{band[0]:+.2f} to {band[1]:+.2f})",
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
    """Figure 13: the gene order, Cacng8's and Gria1's places and their gap under
    other statistics, borders, readings, inputs and structure sets.

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
        f" or more and Cacng8 stays {first}; Gria1's place moves "
        f"({ordinal(others['rank_p9_Gria1'].min())} to "
        f"{ordinal(others['rank_p9_Gria1'].max())}), so it is quoted with its null",
    )
    y = variant_rows(summary)
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
    y = variant_rows(summary)
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
            "a conclusion depends on; a gap that stays inside the pale band under "
            "every choice is a lead these maps cannot resolve.",
        ],
    )
    return saved(fig, save)


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
        rate = float((calibration[column] < 0.05).mean())
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
    """A: one gene against the real map and against a map that knows only each
    structure's division, in the two axes `rects` place."""
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
    """Legend entries of the within panel: filled past the null, hollow not; greys."""
    entries = (
        (dict(color=DARK_GREY), "past the within null (BH, all genes)"),
        (dict(mfc="white", mec=DARK_GREY), "not past it"),
        (dict(color=DARK_GREY), "dark: P9's genes; light: the others"),
        (dict(color=DARK_BLUE), "AMPA receptor subunits"),
    )
    return [
        plt.Line2D([], [], ls="", marker="o", ms=5, label=label, **style)
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
    """C: every gene's whole-brain rho against its mean rho inside divisions, filled
    where the within rho passes its null, the median gene's null band shaded."""
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
    """A of figure 10: each set fixed in advance, a dot per gene, its median against
    the band where the same genes' median falls with 95% of the surrogates."""
    rng = np.random.default_rng(0)
    nano = tests[tests["map"] == "nano"].set_index("gene_set")
    labels = []
    for k, name in enumerate(order):
        test = nano.loc[name]
        rows = members[members["gene_set"] == name]
        set_column(ax, k, rows, test, SET_COLOURS[name], rng)
        if test["tested"]:
            stats = p_text(test["p_spatial"], n_surrogates)
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
    """Figure 10: whether the genes that set surface receptor follow the map better
    than other genes: the gene sets against the null, the localisation genes against
    postsynaptic genes of the same expression, and what that test could find.

    `members` and `tests` are gene_sets.csv and set_tests.csv, `order` the sets in
    drawing order; `table`, `summary`, `pairs` and `power` those of
    plot_localisation.
    """
    nano = tests[(tests["map"] == "nano") & tests["tested"]]
    n_pass = int((nano["q"] < q).sum())
    main = summary.set_index("test").loc["matched controls"]
    found = detectable(power, "power_labels")
    if n_pass:
        sets = f"{n_pass} of {len(nano)} gene sets pass their null"
    else:
        sets = "No gene set passes its null"
    if main["p_labels"] < 0.05:
        local = "the localisation genes stand above postsynaptic genes of the same"
    else:
        local = "the localisation genes do no better than postsynaptic genes of the same"
    fig = plt.figure(figsize=(16, 6.4))
    heading(
        fig,
        "gene_kinds",
        f"{sets}, and {local} expression ({main['difference']:+.3f}, p = "
        f"{main['p_labels']:.2f}); a difference of {found:+.2f} would have been found",
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
    matched = sorted(set(pairs.values()))
    ax = fig.add_axes([0.56, 0.2, 0.17, 0.56])
    strips(
        ax,
        [
            ("localisation", by.loc[loc, "rho_partial"].to_numpy(), RED),
            ("matched\ncontrols", by.loc[matched, "rho_partial"].to_numpy(), DARK_GREY),
        ],
        np.random.default_rng(0),
        "partial rho with the map, subunits removed",
    )
    panel_title(
        ax,
        "B",
        "Localisation against matched controls",
        f"{main['difference']:+.3f} between medians, p = {main['p_labels']:.2f}",
    )

    # C: what the test could find
    ax = fig.add_axes([0.79, 0.2, 0.19, 0.56])
    power_panel(ax, power, float(main["difference"]))
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
    """E: the share of maps on which each test finds a planted localisation effect."""
    for column, colour, label in (
        ("power_labels", DARK_GREY, "labels permuted (the test)"),
        ("power_spatial", NULL_BAND, "against the remainder's surrogates"),
    ):
        edge = DARK_GREY if column == "power_labels" else "#7f9cc9"
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
    """Figure 10s2: localisation genes against expression-matched postsynaptic
    controls.

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
    matched = sorted(set(pairs.values()))
    test = summary.set_index("test")
    main = test.loc["matched controls"]
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
        "; ".join(status)
        + f".  The test: {len(loc)} localisation genes against {len(matched)} "
        f"expression-matched controls, {main['difference']:+.3f} (p = "
        f"{main['p_labels']:.3f}; spatial {p_text(p_spatial, n_surrogates)}){found}",
    )
    ax = fig.add_axes([0.06, 0.58, 0.22, 0.29])
    strips(
        ax,
        [
            ("localisation", by.loc[loc, "rho_partial"].to_numpy(), RED),
            ("matched\ncontrols", by.loc[matched, "rho_partial"].to_numpy(), DARK_GREY),
            ("all\ncontrols", by.loc[ctrl, "rho_partial"].to_numpy(), LIGHT_GREY),
        ],
        rng,
        "partial rho with the map, subunit composite removed",
    )
    panel_title(
        ax,
        "A",
        "The test",
        f"medians {main['median_first']:+.3f} against {main['median_second']:+.3f}",
    )
    ax = fig.add_axes([0.36, 0.58, 0.26, 0.29])
    label_null_panel(ax, nulls["matched controls"], main)
    panel_title(
        ax,
        "B",
        "Its null, labels permuted between the two sets",
        f"p = {main['p_labels']:.3f}; dashed: a difference beyond "
        f"±{main['critical']:.3f} reaches p < 0.05\nagainst maps whose remainder is a "
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
            ("matched\ncontrols", np.log10(level[matched].to_numpy() + 1e-3), DARK_GREY),
            ("all\ncontrols", np.log10(level[ctrl].to_numpy() + 1e-3), LIGHT_GREY),
        ],
        rng,
        "log10 median expression energy",
    )
    panel_title(ax, "D", "The matching", "a quiet gene correlates with nothing")
    if power is not None:
        ax = fig.add_axes([0.36, 0.18, 0.26, 0.26])
        power_panel(ax, power, float(main["difference"]))
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


# ===== 03 and 03s How much of the map Gria1 and synapse density predict =====

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

# the variants of figure 03s G by their key in variants.csv, short; the number of
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


def budget_panel(ax: plt.Axes, n: dict) -> None:
    """D: the variance budget of the main model, of its four-subunit variant and of
    control F.

    The main model's leftover carries a bracket with its share and jackknife
    interval, the calibration floor hatched at its start, and a dark blue mark where
    the leftover would begin if the map left as much as a map of one Allen Gria1
    experiment does.
    """
    steps = n["steps"]
    floor = n["floor"]
    end_model = steps[2]
    four = n["variants"].set_index("key").loc["four_subunits"]
    four_steps = [
        four["abundance"],
        four["abundance"] + four["plus_density"],
        1 - four["left"],
    ]
    end_four = four_steps[2]
    rows = [
        (
            2,
            f"the main model\n({n['abundance']}, synapse density,\nautofluorescence)",
            [
                ("abundance", 0, steps[0], f"{n['abundance']}\n{steps[0]:.0%}"),
                (
                    "density",
                    steps[0],
                    min(steps[1], end_model),
                    f"+ density\n{steps[1] - steps[0]:.0%}",
                ),
                ("autofluorescence", steps[1], end_model, ""),
                ("left", end_model, 1, ""),
            ],
        ),
        (
            1,
            "check row: Gria1 to Gria4\nin place of Gria1, four terms",
            [
                ("abundance", 0, four_steps[0], f"Gria1-4\n{four_steps[0]:.0%}"),
                (
                    "density",
                    four_steps[0],
                    min(four_steps[1], end_four),
                    f"+ density\n{four_steps[1] - four_steps[0]:.0%}",
                ),
                ("autofluorescence", four_steps[1], end_four, ""),
                ("left", end_four, 1, f"left\n{1 - end_four:.0%}"),
            ],
        ),
        (
            0,
            f"control F: components of the {n['f_genes']} genes\nmeasured in every "
            f"structure (most often {n['f_k']})",
            [
                ("components", 0, n["f_share"], f"{n['f_share']:.0%}"),
                ("left", n["f_share"], 1, f"left\n{1 - n['f_share']:.0%}"),
            ],
        ),
    ]
    for y, _, parts in rows:
        budget_bar(ax, y, parts)

    # the floors hatched at the start of the leftovers (the model's, control F's own)
    hatched(ax, 2, end_model, floor["left_median"])
    f_floor = float(np.mean(n["f_floor"]))
    hatched(ax, 0, n["f_share"], f_floor)

    # the model's leftover: a bracket with its share, the interval of where it starts
    left = 1 - end_model
    lo, hi = n["left_ci"]
    ax.plot([end_model, 1], [2.36, 2.36], color="0.1", lw=1.2, zorder=5)
    for x in (end_model, 1):
        ax.plot([x, x], [2.3, 2.36], color="0.1", lw=1.2, zorder=5)
    ax.text(
        1.0,
        2.4,
        f"left {left:.0%} (95% {lo:.0%} to {hi:.0%}, jackknife over structures),\n"
        f"of which the floor {floor['left_median']:.0%} (hatched)",
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="0.1",
        zorder=6,
    )
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

    # the benchmark of a map that is one Allen Gria1 experiment
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
    auto = end_model - steps[1]
    ax.annotate(
        f"+ autofluorescence {auto:+.1%} (alone "
        f"{n['partition']['autofluorescence']:+.0%})",
        xy=(min(steps[1], end_model), 2.28),
        xytext=(0.55, 2.5),
        fontsize=7.4,
        color="0.2",
        ha="right",
        va="bottom",
        arrowprops=dict(arrowstyle="-", color="0.4", lw=0.7),
        zorder=6,
    )
    ax.text(
        n["f_share"] + f_floor / 2,
        -0.42,
        f"hatched: control F's own floor, {n['f_floor'][0]:.0%} to {n['f_floor'][1]:.0%}",
        fontsize=7,
        ha="left",
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
        f"map's reliability and the fit (R2 {n['r2']:.2f})",
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
        ax.text(
            r["left"],
            i - 0.32,
            f"{r['left']:.0%}",
            fontsize=7,
            ha="center",
            va="bottom",
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
        f"left, 95% over structures: {lo:.0%} to {hi:.0%}",
        fontsize=7.5,
        va="center",
        ha="right",
    )
    ax.text(
        0,
        -0.48,
        f"autofluorescence adds {steps[2] - steps[1]:+.1%}",
        fontsize=7.5,
        va="center",
        color=DARK_GREY,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.7, 0.45)
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
    """C of figure 03: what the main model leaves of nano, beside maps whose answer is
    known, on the calibration's structures.

    A dot per draw of made-up adults, the median as a bar; nano is the mean over the
    two halves of the Allen experiments its predictors come from.
    """
    random = calibration[calibration["folds"] == "random"]
    rows = (
        ("abundance and density", "Gria1 + synapse density only\n(the floor)", DARK_GREY),
        ("nano", "the nano map", RED),
        ("Gria1 mRNA", "one Allen Gria1 experiment\n(the benchmark)", DARK_BLUE),
    )
    rng = np.random.default_rng(0)
    for y, (kind, _, colour) in enumerate(rows):
        if kind == "nano":
            ax.scatter(n["nano_cal_left"], y, s=80, color=colour, zorder=3)
            ax.text(
                n["nano_cal_left"],
                y - 0.3,
                f"{n['nano_cal_left']:.0%}",
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
    ax.set_xlim(0, max(0.5, float(random["left"].max()) + 0.04))
    ticks = np.arange(0, ax.get_xlim()[1] + 1e-9, 0.1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map left by the same model")
    tidy(ax)


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
    """Figure 03: what Gria1 expression and synapse density predict of the map, and
    whether what they leave is real.

    The arguments are those of plot_beyond_budget, its detailed version.
    """
    n = numbers
    left = 1 - n["steps"][2]
    floor = n["floor"]["left_median"]
    above = "above" if n["nano_cal_left"] > floor else "no larger than"
    fig = plt.figure(figsize=(15, 10.5))
    heading(
        fig,
        "beyond",
        f"Gria1 and synapse density predict {n['steps'][2]:.0%} of the map's "
        f"reproducible pattern; {left:.0%} is left, it replicates across mice "
        f"({n['rep_left']:.2f}), and it is {above} what Allen-to-Allen mismatch alone "
        f"leaves ({floor:.0%})",
    )

    # A: the map against the main model's held-out prediction
    ax = fig.add_axes([0.07, 0.53, 0.3, 0.36])
    map_scatter(
        ax,
        held_out["model"],
        y,
        groups,
        "what Gria1 and synapse density predict (held out), rank",
    )
    order = np.argsort(residual)
    xr, yr = ranks01(held_out["model"]), ranks01(y)
    extremes = list(order[:NAMED_LEFTOVER]) + list(order[-NAMED_LEFTOVER:])
    spread_labels(ax, [(xr[i], yr[i], acronyms[i], "0.15") for i in extremes])
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
        "bent (x, x², x³)",
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
        f"{len(replication)} ways to split the {n['n_adults']} adults into two fives",
    )
    footer(
        fig,
        [
            "How to read: the map and every predictor are ranks across structures; a "
            "share is the R² on structures the fit has not seen, over what two halves "
            "of the cohort reproduce.",
            "The floor and the benchmark are maps whose answer is known, made from half "
            "of each gene's Allen experiments with ten made-up adults as noisy as ours. "
            f"In detail, with the check rows: {figure_ref('beyond_budget')}.",
        ],
    )
    return saved(fig, save)


def plot_beyond_budget(
    structures: list[str],
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
    """Figure 03s: how much of the map Gria1 and synapse density predict, in detail.

    `y` is the map's ranks on `structures`, `gria1` Gria1's, `held_out` each
    structure's prediction by the density model and by the whole model from fits
    that never saw it, `residual` the model's leftover; `numbers` holds the numbers
    of adult.beyond_figures, `calibration` and `replication` the tables of
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

    # A to C: the map against Gria1, against density, against the whole model
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
    map_scatter(
        ax,
        held_out["model"],
        y,
        groups,
        "what the model predicts (held out), rank",
    )
    order = np.argsort(residual)
    xr, yr = ranks01(held_out["model"]), ranks01(y)
    spread_labels(
        ax,
        [
            (xr[i], yr[i], acronyms[i], "0.15")
            for i in list(order[:NAMED_LEFTOVER]) + list(order[-NAMED_LEFTOVER:])
        ],
    )
    panel_title(
        ax,
        "C",
        "The nano map against the whole model",
        f"{n['abundance']}, synapse density and autofluorescence, bent: "
        f"{n['steps'][2]:.0%} of\nthe reproducible map; the spread off the diagonal is "
        "the leftover",
    )

    # D: the budget; G: its variants
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

    # E: the calibration; F: the replication
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
    footer(
        fig,
        [
            "How to read: the map and every predictor are ranks across structures; each "
            "predictor enters as x, x^2 and x^3, and each share is the held-out R2 "
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
    tidy(ax)


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


# ===== 07 The genes that follow the map =====

# the band behind the rows of Gria1 and Cacng8, across every panel of figure 07
ROW_TINT = "0.93"

# what a gene is set beside in figure 07 C and on its sheet: column, label, colour
# and marker; a column the main model does not have is skipped
LIKENESS = (
    ("rho_Gria1", "Gria1", DARK_BLUE, "o"),
    ("rho_markers", "markers", DENSITY_BLUE, "o"),
    ("rho_psd_pc1", "psd_pc1", DENSITY_BLUE, "s"),
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
    """The rows of figure 12: the genes past the map's null and the named, best first."""
    rows = table[table["top"] | table["named"]]
    return rows.sort_values("rho", ascending=False).reset_index(drop=True)


def own_term(rows: pd.DataFrame) -> np.ndarray:
    """Whether each gene is a term of the main model itself (Gria1, the abundance).

    Its rho with the leftover is then near zero by construction, and its null so
    narrow that a p says nothing; a gene inside a composite term (a marker, a gene
    of psd_pc1) is not.
    """
    return (rows["in_model"] == "abundance").to_numpy()


def follows_leftover(rows: pd.DataFrame) -> np.ndarray:
    """Whether each gene that is not a term itself follows the leftover at p < 0.05."""
    return (rows["leftover_p"] < 0.05).to_numpy() & ~own_term(rows)


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
        ax.scatter(
            rows[column],
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
    alike to the model take (ish.top_genes); a bar is dark past the band, mid grey
    past the tick only.
    """
    y = np.arange(len(rows))
    tint_rows(ax, rows)
    taken = 100 * rows["taken"].to_numpy(float)
    plain = 100 * rows["taken_null_hi"].to_numpy(float)
    alike = 100 * rows["taken_null_alike_hi"].to_numpy(float)
    row_bands(ax, np.zeros(len(rows)), plain)
    free = ~own_term(rows)
    colours = []
    for p_plain, p_alike, outside in zip(rows["p_taken"], rows["p_taken_alike"], free):
        if outside and p_plain < 0.05:
            colours.append(DARK_GREY)
        elif outside and p_alike < 0.05:
            colours.append(MID_GREY)
        else:
            colours.append(LIGHT_GREY)
    ax.barh(y, taken, color=colours, height=0.62, zorder=2)
    ax.scatter(alike, y, marker="|", s=120, color="0.05", linewidths=1.6, zorder=3)
    for i, value in enumerate(taken):
        ax.text(
            max(value, plain[i], alike[i]) + 0.3,
            i,
            points(value / 100),
            fontsize=6.5,
            va="center",
            color=DARK_GREY,
        )
    ax.axvline(0, color="0.3", lw=0.6, zorder=1)
    ax.set_xlim(min(0.0, float(np.min(taken)) - 0.5), float(np.max(plain)) * 1.7)
    ax.set_xlabel("points of the reproducible map")
    gene_rows_axis(ax, rows, labels=False)
    tidy(ax)


def plot_top_genes(
    table: pd.DataFrame,
    q: float,
    t_max: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 07: the genes that follow the map, what kind of map each is, and what
    each takes of the leftover.

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
        "pale: 95% of maps with its smoothness;\nred diamond: also follows the leftover",
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
        "pale: 95% of plain surrogates;\ntick: 95% of maps alike to the model",
    )
    footer(
        fig,
        [
            "How to read: the leftover is what Gria1, synapse density and "
            "autofluorescence leave; D adds the gene to that model and counts the "
            "points of the reproducible map it takes, against maps of its smoothness "
            "in its place.",
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
    """Figure 08: the map against Gria1 and against Cacng8, and Cacng8's lead over
    Gria1 against maps that follow both alike.

    `table` is top_genes.csv indexed by symbol, `map_values` the map on the declared
    structures, `gap` the merged row of gap.csv and `gap_null` the gaps of maps that
    follow both genes alike; `n_map_surrogates` counts the map's surrogates.
    """
    cacng8, gria1 = table.loc["Cacng8"], table.loc["Gria1"]
    groups = group_of(set_table)
    if gap["p_equal"] < 0.05:
        verdict = "a lead past what maps that follow both alike give"
    else:
        verdict = (
            "but maps that follow both alike give such a lead too often to tell "
            "which it follows more"
        )
    fig = plt.figure(figsize=(16, 6.4))
    heading(
        fig,
        "cacng8_gria1",
        f"The map follows Cacng8 ({cacng8['rho']:+.2f}) more closely than Gria1 "
        f"({gria1['rho']:+.2f}), {verdict} (p = {gap['p_equal']:.3f})",
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
    ax.hist(gap_null, bins=np.linspace(-0.6, 0.6, 61), color=NULL_BAND)
    ax.axvspan(gap["equal_lo"], gap["equal_hi"], color=NULL_BAND, alpha=0.4, lw=0)
    ax.axvline(gap["gap"], color=RED, lw=1.8)
    ax.set_xlim(-0.6, 0.6)
    ax.set_xlabel("rho(Cacng8) - rho(Gria1)")
    ax.set_ylabel("maps that follow both alike")
    panel_title(
        ax,
        "C",
        "Cacng8's lead over Gria1",
        f"{gap['gap']:+.2f} on {int(gap['n_structures'])} structures (red), "
        f"{p_text(gap['p_equal'], n_map_surrogates)}",
    )
    tidy(ax)
    footer(
        fig,
        [
            "How to read: one dot per structure, coloured by group of divisions. C: maps "
            "made of both genes' patterns alike plus a surrogate of the map (shaded: "
            f"their 95%); {figure_ref('gene_ranking')} B adds the adults' interval and "
            "each pairing of Allen experiments.",
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
    """D of figure 11: every gene's spatial p against the leftover, smallest first,
    beside the line Benjamini-Hochberg asks a gene to cross; the family dark, Cacng8
    red."""
    order = genes.sort_values("p_spatial")
    p = order["p_spatial"].to_numpy(float)
    k = np.arange(1, len(p) + 1)
    symbols = order["symbol"].to_numpy()
    in_family = np.isin(symbols, family)
    cacng8 = symbols == LEFTOVER_GENE
    ax.plot(k, -np.log10(k * q / len(p)), color=RED, ls="--", lw=1, label=f"BH, q {q}")
    ax.scatter(k[~in_family], -np.log10(p[~in_family]), s=10, color=LIGHT_GREY, lw=0)
    ax.scatter(
        k[in_family & ~cacng8],
        -np.log10(p[in_family & ~cacng8]),
        s=16,
        color=DARK_GREY,
        lw=0,
        label="the family",
    )
    ax.scatter(k[cacng8], -np.log10(p[cacng8]), s=40, color=RED, lw=0, label="Cacng8")
    ax.set_xscale("log")
    ax.set_xlabel("genes, smallest p first")
    ax.set_ylabel("-log10 spatial p")
    ax.legend(loc="upper right", fontsize=7, frameon=False)
    tidy(ax)


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
    """Figure 11: the three tiers named on 8 October against what Gria1 and synapse
    density leave: Cacng8, the AMPA receptor complex family, every other gene.

    `table` is top_genes.csv, `leftover` the main model's leftover on the structures
    of its fit, `profile` Cacng8's, `tests` named_tests.csv, `nulls` per map the
    surrogates' medians of the family, `genes` leftover_genes.csv; `n_surrogates`
    counts the leftover's surrogates.
    """
    by = tests.set_index(["map", "test"])
    first = tests[tests["tier"] == "1"].iloc[0]
    spatial = by.loc[("leftover", "the family's median against the surrogates")]
    matched = by.loc[("leftover", "the family against matched postsynaptic controls")]
    family = table[table["family"]].reset_index(drop=True)
    n_pass = int((genes["q_all"] < q).sum())
    follows = "follows" if first["p"] < 0.05 else "does not follow"
    if spatial["p"] < 0.05 and matched["p"] >= 0.05:
        group = (
            "the family does so beyond the surrogates, not beyond other postsynaptic "
            "genes"
        )
    elif spatial["p"] < 0.05:
        group = "the family does so beyond the surrogates and other postsynaptic genes"
    else:
        group = "the family as a group does not"
    fig = plt.figure(figsize=(17, 6.2))
    heading(
        fig,
        "leftover",
        f"Cacng8, named in advance, {follows} what Gria1 and synapse density leave "
        f"({p_text(first['p'], n_surrogates)}); {group}; {n_pass} of {len(genes)} "
        "genes pass once every gene is corrected for",
    )

    # A, tier 1: Cacng8 against the leftover
    ax = fig.add_axes([0.04, 0.15, 0.19, 0.6])
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

    # B and C, tier 2: the family as a group, against the surrogates and against
    # postsynaptic genes of the same expression
    ax = fig.add_axes([0.29, 0.15, 0.19, 0.6])
    group_null_panel(ax, spatial, nulls["leftover"]["spatial"], RED)
    panel_title(
        ax,
        "B",
        f"Tier 2: the family ({len(family)} genes)",
        f"median {spatial['first']:+.3f} (red), {p_text(spatial['p'], n_surrogates)}",
    )
    ax = fig.add_axes([0.54, 0.15, 0.19, 0.6])
    paired_strips(
        ax, family, "leftover_rho", "control_leftover_rho", np.random.default_rng(0)
    )
    ax.set_ylabel("Spearman rho with the leftover")
    panel_title(
        ax,
        "C",
        "The family against matched genes",
        f"{matched['difference']:+.3f} between medians, p = {matched['p']:.2f}",
    )

    # D, tier 3: every gene
    ax = fig.add_axes([0.79, 0.15, 0.19, 0.6])
    bh_panel(ax, genes, list(family["symbol"]), q)
    panel_title(
        ax,
        "D",
        "Tier 3: every gene, exploratory",
        f"{n_pass} of {len(genes)} cross the line",
    )
    footer(
        fig,
        [
            "How to read: the leftover's surrogates each go through the same fit; tier 1 "
            "is one test, uncorrected; tier 2 is the family as a group, against the "
            "surrogates and against expression-matched postsynaptic genes, each one "
            "test.",
            f"Every gene and gene set: {figure_ref('leftover_genes')}; the family member "
            f"by member and on the map: {figure_ref('ampa_family')}.",
        ],
    )
    return saved(fig, save)


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
    by = tests.set_index(["map", "test"])
    first = tests[tests["tier"] == "1"].iloc[0]
    spatial = by.loc[("leftover", "the family's median against the surrogates")]
    matched = by.loc[("leftover", "the family against matched postsynaptic controls")]
    on_map = by.loc[("nano", "the family's median against the surrogates")]
    on_map_matched = by.loc[("nano", "the family against matched postsynaptic controls")]
    past = rows.loc[rows["leftover_q_family"] < q, "symbol"]
    beyond = "follows it beyond" if spatial["p"] < 0.05 else "does not pass"
    than = "more" if matched["p"] < 0.05 else "no more"
    cacng8 = "follows" if first["p"] < 0.05 else "does not follow"
    untested = members.loc[~members["tested"], "symbol"]
    fig = plt.figure(figsize=(16, 11.5))
    heading(
        fig,
        "ampa_family",
        f"Cacng8, named in advance, {cacng8} the leftover "
        f"({p_text(first['p'], n_surrogates)}); the family {beyond} the surrogates "
        f"(p = {spatial['p']:.3f}), {than} than postsynaptic genes of its expression "
        f"(p = {matched['p']:.3f})",
    )

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
    group_null_panel(ax, spatial, nulls["leftover"]["spatial"], RED)
    panel_title(
        ax,
        "B",
        "As a group, against the surrogates",
        f"median {spatial['first']:+.3f} (red); {p_text(spatial['p'], n_surrogates)}",
    )
    rng = np.random.default_rng(0)
    ax = fig.add_axes([0.77, 0.58, 0.2, 0.27])
    paired_strips(ax, rows, "leftover_rho", "control_leftover_rho", rng)
    ax.set_ylabel("Spearman rho with the leftover")
    panel_title(
        ax,
        "C",
        "Against genes of the same expression",
        f"{matched['difference']:+.3f} between medians, p = {matched['p']:.3f}\n"
        f"(labels permuted; ±{matched['critical']:.2f} would give p < 0.05)",
    )

    # D and E: the same on the map itself, to describe the family
    ax = fig.add_axes([0.5, 0.13, 0.2, 0.3])
    group_null_panel(ax, on_map, nulls["nano"]["spatial"], NANO)
    panel_title(
        ax,
        "D",
        "On the map itself (a description)",
        f"median {on_map['first']:+.3f} (orange); "
        f"{p_text(on_map['p'], n_map_surrogates)}",
    )
    ax = fig.add_axes([0.77, 0.13, 0.2, 0.3])
    paired_strips(ax, rows, "rho", "control_rho", rng)
    ax.set_ylabel("Spearman rho with the nano map")
    panel_title(
        ax,
        "E",
        "On the map, against the same controls",
        f"{on_map_matched['difference']:+.3f} between medians, p = "
        f"{on_map_matched['p']:.3f}",
    )
    footer(
        fig,
        [
            "How to read: the family was named on 8 October, before the main model ran, "
            "from sources outside this analysis: S, the native AMPA receptor complexes "
            "of Schwenk et al. 2012 (Figure 1D, Table S2); GO, GO:0032281 AMPA",
            "glutamate receptor complex; P, the partner subunits (dark blue); minus "
            "Gria1, the abundance term. Tier 1 is Cacng8 alone (red); tier 2 the family "
            "as a group (B, C), then gene by gene, BH within the family (A, bold).",
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
        f"{row['leftover_q_all']:.2f}",
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
    """One gene characterised: its map, with the map and the leftover, what it takes,
    what it is like, and what the gene table says it is.

    `row` is the gene's row of top_genes.csv, `leftover` the main model's leftover
    on the structures of its fit, `taken_nulls` what maps of the gene's smoothness
    take in its place (plain, alike); `n_surrogates` and `n_map_surrogates` count
    the leftover's and the map's surrogates. PNG only, for review.
    """
    symbol = row["symbol"]
    groups = group_of(set_table)
    own = row["in_model"] == "abundance"
    why = []
    if row["top"]:
        why.append("past the map's null after BH")
    if symbol == LEFTOVER_GENE:
        why.append("named for the leftover in advance")
    elif row["named"]:
        why.append("the gene of the stained protein")
    if row["family"]:
        why.append("a member of the AMPA receptor complex family")
    fig = plt.figure(figsize=(16, 10))
    fig.text(
        0.5,
        0.975,
        f"{symbol}, {row['name']}: {'; '.join(why)}",
        ha="center",
        va="top",
        fontsize=12,
        color=named_colour(symbol),
    )
    fig.text(
        0.5,
        0.945,
        f"with the map {row['rho']:+.2f} "
        f"({p_text(row['p_spatial'], n_map_surrogates)}, rank "
        f"{int(row['rank_all'])}); inside divisions {row['rho_within']:+.2f}; with the "
        f"leftover {row['leftover_rho']:+.2f} "
        f"({p_text(row['leftover_p'], n_surrogates)}); "
        + (
            "a term of the main model"
            if own
            else f"takes {points(row['taken'])} points of the reproducible map (p = "
            f"{row['p_taken']:.3f} plain, {row['p_taken_alike']:.3f} alike)"
        ),
        ha="center",
        va="top",
        fontsize=9,
        color=DARK_GREY,
    )

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

    # D: what it takes of the leftover, against maps of its smoothness
    ax = fig.add_axes([0.37, 0.08, 0.25, 0.3])
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
    ax.legend(loc="upper right", fontsize=7, frameon=False)
    tidy(ax)
    if own:
        numbers = "a term of the main model: adding it again takes nothing"
    else:
        numbers = (
            f"{points(row['taken'])} points; p = {row['p_taken']:.3f} against plain "
            f"surrogates, {row['p_taken_alike']:.3f} against maps alike"
        )
    panel_title(ax, "D", "Added to the main model", numbers)

    # E: what it is like
    ax = fig.add_axes([0.71, 0.08, 0.25, 0.3])
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
    panel_title(ax, "E", "What it is like")

    # F: the numbers, and what the gene table says it is
    ax = fig.add_axes([0.01, 0.03, 0.3, 0.4])
    ax.axis("off")
    lines = numbers_block(row) + [""] + annotation_lines(row, 64)
    ax.text(0, 1, "\n".join(lines), fontsize=7.5, va="top", ha="left", linespacing=1.4)
    panel_title(ax, "F", "Numbers, and what the gene table says")
    return saved(fig, save, eps=False)


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
    """Per adult, which channel follows which: SEP with autofluorescence, SEP with
    nano, nano with autofluorescence, joined by adult."""
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


def tinted_box(ax: plt.Axes, x: float, y: float, w: float, h: float, part: str):
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


def draw_part(ax: plt.Axes, x: float, y: float, width: float, height: float, part: dict):
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
    """The figures in three columns, a group per part of the walk; returns the inch
    below it.

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
    """Figure 00: the question, the two parts of the argument with their key numbers,
    the limit, and the map of the figures.

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
        f"{FIGURES['overview']}.  {QUESTIONS['overview']}",
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
    """Figure 16: what is left of April's headline: the gene order, the category p
    under each choice it rests on, and today's groups against the null.

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
    """C of figure 14: per density measure, what it predicts alone and what the main
    model leaves with it, on the structures PSD95 covers."""
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


def plot_synaptome(
    density: pd.DataFrame,
    coverage: pd.DataFrame,
    markers: pd.Series,
    variants: pd.DataFrame,
    min_coverage: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 14: whether the measured synapse density covers the fit, and whether it
    predicts the map better than the mRNA it would replace.

    `density` and `coverage` are the tables of run_synaptome, `markers` the marker
    composite per structure of the fit, `variants` beyond/variants.csv.
    """
    table = density.set_index("structure")
    n_fit, n_measured, needed = coverage_counts(table, min_coverage)
    measured = table[table["in_fit"] & table["measured"]]
    rows = variants.set_index("key")
    psd95, panel = rows.loc["psd95"], rows.loc["panel"]
    rule = "not below" if n_measured >= needed else "below"
    better = "better" if psd95["density_alone"] > panel["density_alone"] else "less well"
    fig = plt.figure(figsize=(16, 5.8))
    heading(
        fig,
        "synaptome",
        f"PSD95 puncta cover {n_measured / n_fit:.0%} of the fit, {rule} the "
        f"{min_coverage:.0%} the rule asks; where measured they predict the map {better} "
        f"than the mRNA panel ({psd95['density_alone']:.0%} against "
        f"{panel['density_alone']:.0%}) and leave {psd95['left']:.0%} against "
        f"{panel['left']:.0%}",
    )
    ax = fig.add_axes([0.07, 0.17, 0.24, 0.62])
    coverage_panel(ax, coverage)
    panel_title(
        ax,
        "A",
        "Structures of the fit with a measured density",
        f"{n_measured} of {n_fit}; the rule asks for {needed}",
    )
    ax = fig.add_axes([0.39, 0.17, 0.22, 0.62])
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
    ax = fig.add_axes([0.76, 0.17, 0.21, 0.62])
    density_rows_panel(ax, variants)
    panel_title(
        ax,
        "C",
        "Each density measure in the model",
        f"on the {int(psd95['n_structures'])} structures PSD95 covers",
    )
    footer(
        fig,
        [
            "How to read: PSD95 punctum density from one adult mouse (Zhu et al. 2018, "
            "as shared by Hansen et al.), placed in the CCF structures; the mRNA panel "
            "is the main model's density. In detail, with the two hemispheres and every "
            f"density's agreement: {figure_ref('synaptome_detail')}.",
        ],
    )
    return saved(fig, save)


def plot_synaptome_detail(
    density: pd.DataFrame,
    coverage: pd.DataFrame,
    agreement: pd.DataFrame,
    markers: pd.Series,
    min_coverage: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 14s: the measured synapse density, how much of the fit it covers, how
    it compares.

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
    if n_measured >= needed:
        verdict = "it replaces the mRNA density terms in the main model"
    else:
        verdict = (
            "the mRNA terms stay in the main model, and PSD95 is a variant on the "
            f"{n_measured} structures it covers"
        )
    heading(
        fig,
        "synaptome_detail",
        f"PSD95 puncta cover {n_measured / n_fit:.0%} of the fit, "
        f"{'not below' if n_measured >= needed else 'below'} the {min_coverage:.0%} "
        f"the rule asks: {verdict}",
    )
    footer(
        fig,
        [
            "Zhu et al. 2018, one adult male mouse, as Hansen et al. share it. PSD95 "
            "density: the mean of the 30 subtypes whose puncta hold PSD95, each scaled "
            "0..1 as shared;",
            "never a punctum's intensity or size. A structure is measured when it or one "
            "of its parts was sampled; a region above several structures is never spread "
            "onto",
            "them. Every rho is Spearman over the structures both maps have.",
        ],
    )
    return saved(fig, save)
