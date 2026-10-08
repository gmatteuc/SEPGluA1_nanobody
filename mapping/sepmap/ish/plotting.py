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

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
from scipy.stats import rankdata, spearmanr

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

# the guided figures in reading order, which is the order of the argument: the question
# (00), the inputs (01, 02), how much of the map receptor mRNA and synaptic density
# leave (03, 04), what the leftover looks like through the genes (05 to 10), the
# controls and limits (11 to 13), and the April headline as an appendix (14). A
# figure's title, its file name and the figures' references to each other read this
# table; the docstrings and the run scripts' headers that name the files follow it by
# hand
FIGURES = {
    "overview": 0,
    "structures": 1,
    "genes": 2,
    "beyond_budget": 3,
    "beyond_where": 4,
    "one_comparison": 5,
    "spatial_null": 6,
    "gene_ranking": 7,
    "gene_sets": 8,
    "localisation": 9,
    "between_within": 10,
    "autofluorescence": 11,
    "robustness": 12,
    "green_channel": 13,
    "april_headline": 14,
}

# the question each guided figure answers: its title, and its heading in the index of
# the figures (ish.overview)
QUESTIONS = {
    "overview": "The ISH line on one page: the question, the argument, and which figure "
    "answers what",
    "structures": "Which structures enter every comparison, and how much does the "
    "declared reference move zref?",
    "genes": "Which genes, which Allen experiments, and how trustworthy is each map?",
    "beyond_budget": "How much of the nano map do receptor mRNA and synaptic density "
    "predict, and is what they leave real?",
    "beyond_where": "Where does the leftover live, and does any gene's map follow it?",
    "one_comparison": "What does a gene's rho with the map mean, concretely?",
    "spatial_null": "How large a rho do unrelated smooth maps give, and do the "
    "surrogates have the nano map's smoothness?",
    "gene_ranking": "Which genes' maps order the structures as the nano map does, "
    "beyond what a map with its smoothness would?",
    "gene_sets": "Do kinds of genes defined before looking match the map better than "
    "others, beyond the null?",
    "localisation": "Once receptor abundance is removed, do the genes that put AMPA "
    "receptors at the membrane predict the map better than other postsynaptic genes?",
    "between_within": "Does a gene follow the map inside divisions, or only through the "
    "contrast between them?",
    "autofluorescence": "Would the tissue's own autofluorescence, in the same sections, "
    "give the same gene ranking?",
    "robustness": "Does the order of the genes, and where Cacng8 and Gria1 sit, change "
    "with the choices made?",
    "green_channel": "Does the green channel report the tagged receptor, or the tissue?",
    "april_headline": "What is left of April's headline, the category violins and their "
    "p = 0.032?",
}


# ===== Shared pieces =====


def figure_file(key: str) -> str:
    """The file name of a guided figure: its number and its key, 07_gene_ranking.png."""
    return f"{FIGURES[key]:02d}_{key}.png"


def figure_ref(key: str) -> str:
    """How a figure names another: 'figure 06'."""
    return f"figure {FIGURES[key]:02d}"


def heading(fig: plt.Figure, key: str, numbers: str) -> None:
    """The figure's number and question, and under it the line of computed numbers."""
    title = f"{FIGURES[key]}.  {QUESTIONS[key]}"
    fig.text(0.5, 0.985, title, ha="center", va="top", fontsize=12)
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
# grey, flagged red; a dim section kept as a true absence is hatched over light grey
STATUS_CODES = {"ok": 1, "flagged": 2, "absence kept": 3}
STATUS_COLOURS = ["white", "#e6e6e6", RED]


def status_matrix(sections: pd.DataFrame, experiments: list[str]) -> np.ndarray:
    """Experiments by sections, each cell a status code (0 for not judged)."""
    n_sections = int(sections["section"].max()) + 1
    out = np.zeros((len(experiments), n_sections), dtype=int)
    row_of = {e: i for i, e in enumerate(experiments)}
    for eid, k, status in zip(
        sections["experiment_id"], sections["section"], sections["status"]
    ):
        if eid in row_of:
            out[row_of[eid], int(k)] = STATUS_CODES.get(status, 0)
    return out


def draw_status(ax: plt.Axes, matrix: np.ndarray) -> None:
    """Draw a status matrix: white, light grey, red, and hatching for kept absences."""
    shown = np.where(matrix == 3, 1, matrix)
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
    shown = ok[(ok["n_flagged"] > 0) | (ok["n_absence_kept"] > 0)]
    shown = shown.sort_values(["symbol", "experiment_id"])
    experiments = list(shown["experiment_id"])
    matrix = status_matrix(sections, experiments)
    height = 2.2 + 0.15 * len(experiments)
    fig = plt.figure(figsize=(11, height))
    ax = fig.add_axes([0.2, 1.0 / height, 0.75, 1 - 1.9 / height])
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
    fig.text(
        0.5,
        1 - 0.25 / height,
        f"Section QC: {len(shown)} of {len(ok)} experiments with a dim section; "
        f"{n_sections} sections set missing (red), {n_kept} kept as true absence "
        "(hatched); P9's genes in bold",
        ha="center",
        va="top",
        fontsize=10,
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
        f"{summary['n_judged']} sections judged; red: set missing; hatched: dim "
        "but kept as a true absence; empty: no data; grey triangle: the section below",
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
        "both planes": DENSITY_BLUE,
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
    panel_title(
        ax,
        "C",
        "P9's experiments with a dim section",
        "red: set missing; hatched: kept as true absence (proposed)\n"
        f"{len(shown)} of P9's {len(p9_experiments & set(ok['experiment_id']))}; "
        f"{len(dim) - len(shown)} more in the other experiments "
        "(qc/00_flagged.png)",
    )


def reliability_panel(ax: plt.Axes, rel: pd.DataFrame) -> None:
    """D: how reliable one Allen map is, the named genes marked."""
    values = rel["reliability"].dropna()
    ax.hist(values, bins=np.linspace(-0.4, 1.0, 36), color=MID_GREY)
    ax.axvline(0.3, color=DARK_GREY, ls="--", lw=0.9)
    top = ax.get_ylim()[1]
    for i, gene in enumerate(NAMED_GENES):
        value = rel.set_index("symbol")["reliability"].get(gene, np.nan)
        if np.isfinite(value):
            colour = DARK_BLUE if gene == "Gria1" else RED
            ax.axvline(value, color=colour, lw=1.2)
            ax.text(
                value - 0.01,
                top * (0.92 - 0.09 * i),
                f"{gene} {value:.2f}",
                ha="right",
                fontsize=7.5,
            )
    q1, median, q3 = np.percentile(values, [25, 50, 75])
    ax.set_xlabel("reliability: median Spearman between two experiments of a gene")
    ax.set_ylabel("genes")
    panel_title(
        ax,
        "D",
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
    colours = [DARK_BLUE if s in subunits else DARK_GREY for s in have["symbol"]]
    ax.scatter(x, y, s=20, c=colours, alpha=0.7, linewidths=0)
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
    lines = [f"Experiments left out ({len(dropped)}; * P9's own):"]
    kinds = dropped["exclude_reason"].map(reason_kind)
    for kind in dict.fromkeys(kinds):
        mine = dropped[kinds == kind]
        genes = [
            f"{g}*" if p9 else g for g, p9 in zip(mine["symbol"], mine["p9_experiment"])
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
    sections: pd.DataFrame,
    summary: pd.DataFrame,
    rel: pd.DataFrame,
    subunits: set[str],
    save: Path | None = None,
) -> plt.Figure:
    """Figure 02: the genes, their experiments, the section QC and the reliability.

    `genes` holds a row per gene with a usable experiment (p9_gene, ontology_role,
    gene_sets), `experiments` the gene table (excluded, plane), `sections` and
    `summary` the section QC, `rel` the reliability per gene.
    """
    fig = plt.figure(figsize=(16, 10.5))
    n_listed = experiments["symbol"].nunique()
    n_used = int((~experiments["excluded"]).sum())
    heading(
        fig,
        "genes",
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
    footer(
        fig,
        [
            "How to read: one table holds P9's 100 genes and the ontology panel's 390, "
            "a row per Allen experiment (A9). A section five-fold dimmer than its "
            "neighbours is set missing, never filled in (A2); a gene's profile is the "
            "mean rank of its usable experiments.",
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
            ax.annotate(
                acronym,
                (x[i], y[i]),
                xytext=(4, -9),
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


def plot_one_comparison(
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
    """Figure 05: what one gene's rho with the map is, for nano and a few genes.

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
        "one_comparison",
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


def rates_panel(ax: plt.Axes, calibration: pd.DataFrame) -> None:
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
        ax, "D", "False positives on maps with no relation", f"{n_tests} tests per design"
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
    for ax, (title, values) in zip(axes, panels):
        image = draw_plane(ax, paint(lab, values, names), lab, hot_cut(), RANK_FLOOR, 1)
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


def plot_spatial_null(
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
    """Figure 06: why a spatial null, and the null itself.

    `map_ranks` and each of `surrogate_ranks` give a 0-1 rank per structure; `genes`,
    when given, holds (name, observed rho, null rhos, p) for panel E, which needs
    the gene ranking and is left out otherwise.
    """
    fig = plt.figure(figsize=(16, 10.5))
    spatial = {}
    for design in dict.fromkeys(calibration["design"]):
        sub = calibration[calibration["design"] == design]
        spatial[design] = float((sub["p_spatial"] < 0.05).mean())
    heading(
        fig,
        "spatial_null",
        f"{n_surrogates} surrogates of the nano map on {len(map_ranks)} structures; "
        f"spatial p below 0.05 in {spatial.get('pair', np.nan):.1%} of pairs of random "
        f"maps and {spatial.get('map', np.nan):.1%} of random maps against nano "
        "(5% expected)",
    )
    null_rho_panel(fig.add_axes([0.05, 0.6, 0.25, 0.28]), calibration)
    variogram_panel(fig.add_axes([0.37, 0.6, 0.27, 0.28]), variogram, "nano")
    rates_panel(fig.add_axes([0.72, 0.6, 0.25, 0.28]), calibration)
    plane_axes = [fig.add_axes([0.02 + 0.155 * i, 0.12, 0.15, 0.36]) for i in range(4)]
    surrogate_planes(fig, plane_axes, map_ranks, surrogate_ranks, lab, names)
    plane_axes[0].text(
        0,
        1.15,
        f"C.  The nano map and three of its surrogates, CCF plane {plane}",
        transform=plane_axes[0].transAxes,
        fontsize=9,
    )
    if genes:
        gene_null_panel(fig.add_axes([0.75, 0.17, 0.22, 0.28]), genes)
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


# ===== 07 Which expression maps look like the nano map =====


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
        ("dot", "yellow dot: rho with the autofluorescence map"),
        ("bold", f"bold name: past the null, BH within P9's genes, q < {q}"),
        ("blue", "dark blue name: an AMPA receptor subunit"),
    ]
    for k, (kind, text) in enumerate(rows):
        yy = 0.9 - k * 0.15
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
    """B: the Cacng8 - Gria1 gap, its spatial null, its adult interval, the pairings."""
    merged = gap[gap["kind"] == "merged profiles"].iloc[0]
    pairs = gap[gap["kind"] == "experiment pairing"]
    first, second = merged["first"], merged["second"]
    bins = np.linspace(-0.8, 0.8, 65)
    ax.hist(null, bins=bins, color=NULL_BAND, label="surrogates of the map")
    top = ax.get_ylim()[1]
    ax.axvline(merged["gap"], color=RED, lw=1.8, label="observed")
    ax.plot(
        [merged["boot_lo"], merged["boot_hi"]],
        [top * 1.08] * 2,
        color=DARK_GREY,
        lw=2.2,
        solid_capstyle="butt",
    )
    ax.text(
        merged["boot_hi"] + 0.02,
        top * 1.08,
        "95% over resampled adults",
        fontsize=7,
        va="center",
    )
    ax.scatter(
        pairs["gap"],
        [top * 1.2] * len(pairs),
        marker="|",
        s=90,
        color="0.1",
        linewidths=1.4,
    )
    ax.text(
        pairs["gap"].max() + 0.02,
        top * 1.2,
        f"each pairing of their Allen experiments ({len(pairs)})",
        fontsize=7,
        va="center",
    )
    ax.set_ylim(0, top * 1.3)
    ax.set_xlabel(f"rho({first}) - rho({second}), on the structures both have")
    ax.set_ylabel("surrogates")
    ax.legend(loc="upper left", fontsize=7)
    panel_title(
        ax,
        "B",
        f"The {first} - {second} gap",
        f"{merged['gap']:+.3f} on {int(merged['n_structures'])} structures "
        f"({merged['rho_first']:+.2f} against {merged['rho_second']:+.2f}); spatial "
        f"{p_text(merged['p_spatial'], n_surrogates)}",
    )
    tidy(ax)


def counts_panel(ax: plt.Axes, ranking: pd.DataFrame, q: float) -> None:
    """C: genes past the null at the BH level, nano against autofluorescence."""
    groups = (("p9", "P9's genes"), ("all", "every gene"))
    for i, (within, label) in enumerate(groups):
        for j, (name, colour) in enumerate((("nano", NANO), ("auto", AUTO))):
            mine = ranking[ranking["map"] == name]
            if within == "p9":
                mine = mine[mine["p9_gene"]]
            n_pass = int((mine[f"q_{within}"] < q).sum())
            x = i + (j - 0.5) * 0.38
            ax.bar(x, n_pass, width=0.36, color=colour)
            ax.text(
                x, n_pass + 0.6, f"{n_pass} of {len(mine)}", ha="center", fontsize=7.5
            )
    ax.set_xticks([0, 1])
    ax.set_xticklabels([g[1] for g in groups])
    ax.set_ylabel(f"genes past the null (BH q < {q})")
    ax.legend(
        handles=[
            Patch(color=NANO, label="nano map"),
            Patch(color=AUTO, label="autofluorescence map"),
        ],
        loc="upper left",
    )
    panel_title(
        ax,
        "C",
        "How many genes pass",
        "BH within P9's genes, and within all "
        f"(which ones: {figure_ref('autofluorescence')})",
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
    """Figure 07: P9's genes ranked against the map, with the null; the gap; counts.

    `ranking` is gene_ranking.csv, `gap` gap.csv and `gap_null` the gap of every
    surrogate.
    """
    nano = ranking[(ranking["map"] == "nano") & ranking["p9_gene"]]
    nano = nano.sort_values("rho", ascending=False)
    auto = ranking[ranking["map"] == "auto"].set_index("symbol")["rho"]
    n_pass = int((nano["q_p9"] < q).sum())
    fig = plt.figure(figsize=(16, 17))
    heading(
        fig,
        "gene_ranking",
        f"P9's {len(nano)} genes on the declared structures; {n_pass} past the "
        f"spatial null at BH q < {q} within them; {n_surrogates} surrogates per gene",
    )
    ax_a = fig.add_axes([0.08, 0.1, 0.36, 0.785])
    gene_bars(ax_a, nano, auto, subunits, t_max, q)
    ax_a.set_title(
        "A.  P9's genes, best first\nthe column at the right: the gene's Allen "
        "reliability ('-': measured once)\n\n",
        loc="left",
        fontsize=9,
    )
    bar_key(fig.add_axes([0.56, 0.78, 0.4, 0.14]), t_max, q)
    gap_panel(fig.add_axes([0.58, 0.47, 0.38, 0.24]), gap, gap_null, n_surrogates)
    counts_panel(fig.add_axes([0.58, 0.15, 0.26, 0.22]), ranking, q)
    footer(
        fig,
        [
            "How to read: a bar past its pale band is a rho that fewer than 5% of maps "
            "with the nano map's smoothness reach with that gene; BH then allows for "
            f"testing {len(nano)} genes.",
            "The gap asks whether the map follows Cacng8 (a TARP, which brings AMPA "
            "receptors to the surface and holds them at synapses) more closely than "
            "Gria1, the receptor's own mRNA.",
            "What would mean what: Gria1 past its band and autofluorescence not: the "
            "map follows receptor expression, not the tissue. The gap past its null: "
            "the map follows a regulator of surface receptor more closely than the "
            "receptor's mRNA,",
            "as a surface-fraction reading predicts (a correspondence, not a "
            "measurement: Cacng8 is rich in hippocampus, as the map is). The gap inside "
            "its null: these maps cannot tell the two apart. Nothing past the null: the "
            "ranking is descriptive only.",
        ],
    )
    return saved(fig, save)


# ===== 11 Is the ranking the nanobody's or the tissue's =====

# genes always named on the scatter of figure 11, beside those far from the diagonal
NAMED_AUTO = ("Cacng8", "Gria1", "Aqp4")

# the gene lists of figure 11 D: lines per column, and the width of a column in the
# panel's width
PASSING_LINES = 14
PASSING_WIDTH = 0.19


def auto_scatter(
    ax: plt.Axes, ranking: pd.DataFrame, subunits: set[str], n_named: int = 8
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
    for g in named:
        i = genes.index(g)
        ax.annotate(
            g,
            (x[i], y[i]),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=7,
            color=gene_colour(g, subunits),
        )
    agree = spearmanr(x, y).statistic
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("rho with the nano map")
    ax.set_ylabel("rho with the autofluorescence map")
    panel_title(
        ax,
        "A",
        "Each gene against both maps",
        f"{len(genes)} genes (P9's dark, subunits blue); the two gene orders agree "
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
        ax.hist(mine["rho"], bins=bins, color=colour, alpha=0.75, label=label)
        n_all = int((mine["q_all"] < q).sum())
        p9 = mine[mine["p9_gene"]]
        n_p9 = int((p9["q_p9"] < q).sum())
        lines.append(
            f"{label}: median rho {mine['rho'].median():+.2f}; past the null "
            f"{n_all} of {len(mine)} (P9's {n_p9} of {len(p9)})"
        )
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_xlabel("Spearman rho with the map")
    ax.set_ylabel("genes")
    ax.legend(loc="upper left")
    panel_title(ax, "B", "The genes' rho with each map", "\n".join(lines))
    tidy(ax)


def per_adult_panel(ax: plt.Axes, per_adult: pd.DataFrame, genes: tuple[str, ...]):
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
    panel_title(ax, "C", "Adult by adult", f"{n_mice} adults, a line joins one brain")
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
        f"D.  The genes past each map's null (BH within every gene, q < {q}); "
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
    """Figure 11: the same ranking on the autofluorescence map of the same sections.

    `ranking` is gene_ranking.csv (both maps), `per_adult` the per-adult rho table.
    """
    fig = plt.figure(figsize=(16, 11))
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
    heading(
        fig,
        "autofluorescence",
        "; ".join(parts),
    )
    auto_scatter(fig.add_axes([0.05, 0.4, 0.31, 0.47]), ranking, subunits)
    rho_distributions(fig.add_axes([0.43, 0.56, 0.24, 0.29]), ranking, q)
    per_adult_panel(
        fig.add_axes([0.75, 0.56, 0.22, 0.29]), per_adult, ("Gria1", "Cacng8")
    )
    passing_panel(fig.add_axes([0.43, 0.1, 0.55, 0.34]), ranking, q)
    footer(
        fig,
        [
            "How to read: the autofluorescence map is the unlabelled channel of the "
            "same sections, read as nano is (zref over the declared set) and tested "
            "with surrogates of its own (A8); a line in C joins one brain.",
            "What would mean what: a gene past the nano null and not past the "
            "autofluorescence one (in every adult, C) follows the label and not the "
            "tissue. Genes past the autofluorescence null share the tissue's own "
            "pattern;",
            "the nano ranking is the label's own as far as its gene order departs from "
            "the autofluorescence one (A).",
        ],
    )
    return saved(fig, save)


# ===== 12 Does the ranking depend on the choices made =====

# the colour of each kind of choice in figure 12
KIND_COLOURS = {
    "primary": "0.1",
    "statistic": RED,
    "borders": DENSITY_BLUE,
    "reading": NANO,
    "inputs": DARK_GREY,
    "structures": "#3a6db5",
    "5 October": MID_GREY,
}

# the variants drawn against the primary in panel C, one per kind of choice
SMALL_MULTIPLES = ("pearson_log2", "eroded_both", "ratio", "every_structure")


def agreement_panel(ax: plt.Axes, summary: pd.DataFrame) -> None:
    """A: each variant's gene order against the primary's, over P9's genes."""
    y = np.arange(len(summary))
    for i, (_, r) in enumerate(summary.iterrows()):
        colour = KIND_COLOURS.get(r["kind"], DARK_GREY)
        ax.plot([0.5, r["agreement_p9"]], [i, i], color=LIGHT_GREY, lw=1, zorder=1)
        ax.scatter(r["agreement_p9"], i, s=40, color=colour, zorder=2)
        ax.text(
            r["agreement_p9"] + 0.012,
            i,
            f"{r['agreement_p9']:.3f}",
            va="center",
            fontsize=7.5,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(summary["label"], fontsize=8)
    for tick, kind in zip(ax.get_yticklabels(), summary["kind"]):
        tick.set_color(KIND_COLOURS.get(kind, DARK_GREY))
    ax.set_ylim(len(summary) - 0.5, -0.5)
    ax.set_xlim(0.5, 1.08)
    ax.set_xlabel("agreement with the primary order\n(Spearman over P9's genes)")
    panel_title(ax, "A", "The order of the genes", "")
    tidy(ax)


def ranks_panel(ax: plt.Axes, summary: pd.DataFrame) -> None:
    """B: where Cacng8 and Gria1 sit among P9's genes under each variant."""
    for gene, marker, colour in (("Cacng8", "o", RED), ("Gria1", "s", DARK_BLUE)):
        ax.scatter(
            summary[f"rank_p9_{gene}"],
            np.arange(len(summary)),
            marker=marker,
            s=36,
            color=colour,
            label=gene,
            zorder=2,
        )
    for i, (_, r) in enumerate(summary.iterrows()):
        ax.text(
            r["rank_p9_Gria1"] + 1,
            i,
            f"{int(r['rank_p9_Gria1'])}",
            va="center",
            fontsize=7,
            color=DARK_BLUE,
        )
    ax.set_ylim(len(summary) - 0.5, -0.5)
    ax.set_yticks([])
    top = float(np.nanmax(summary["rank_p9_Gria1"]))
    ax.set_xlim(0, top + 6)
    ax.set_xlabel("rank among P9's genes (1 = highest rho)")
    ax.legend(loc="lower right", fontsize=7)
    panel_title(ax, "B", "Where the two genes sit", "")
    tidy(ax)


def gap_rows_panel(ax: plt.Axes, summary: pd.DataFrame) -> None:
    """B, right: the Cacng8 - Gria1 gap under each variant."""
    y = np.arange(len(summary))
    ax.barh(y, summary["gap"], color=MID_GREY, height=0.55)
    for i, g in enumerate(summary["gap"]):
        ax.text(max(g, 0) + 0.005, i, f"{g:+.2f}", va="center", fontsize=7)
    ax.axvline(0, color="0.3", lw=0.6)
    ax.set_ylim(len(summary) - 0.5, -0.5)
    ax.set_yticks([])
    ax.set_xlim(min(0, summary["gap"].min()) - 0.02, summary["gap"].max() + 0.08)
    ax.set_xlabel("rho(Cacng8) - rho(Gria1)")
    tidy(ax)


def multiples_panel(
    axes: list[plt.Axes],
    robustness: pd.DataFrame,
    summary: pd.DataFrame,
    subunits: set[str],
) -> None:
    """C: each gene's rho under four variants against its rho in the primary."""
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
    robustness: pd.DataFrame,
    subunits: set[str],
    save: Path | None = None,
) -> plt.Figure:
    """Figure 12: the ranking under other statistics, borders, readings, inputs, sets.

    `summary` is robustness_summary.csv in the order of its rows, `robustness`
    ranking_robustness.csv.
    """
    fig = plt.figure(figsize=(16, 11.5))
    others = summary[summary["variant"] != "primary"]
    lo, hi = int(others["rank_p9_Cacng8"].min()), int(others["rank_p9_Cacng8"].max())
    cacng8 = f"Cacng8 {lo} to {hi}"
    if hi == 1:
        cacng8 = "Cacng8 first in every variant"
    heading(
        fig,
        "robustness",
        f"{len(others)} variants of the primary ranking; the gene order agrees at "
        f"{others['agreement_p9'].min():.2f} to {others['agreement_p9'].max():.3f} "
        f"over P9's genes; Gria1 ranks {int(others['rank_p9_Gria1'].min())} to "
        f"{int(others['rank_p9_Gria1'].max())}, {cacng8}",
    )
    agreement_panel(fig.add_axes([0.27, 0.5, 0.24, 0.38]), summary)
    ranks_panel(fig.add_axes([0.56, 0.5, 0.2, 0.38]), summary)
    gap_rows_panel(fig.add_axes([0.8, 0.5, 0.16, 0.38]), summary)
    fig.axes[-1].set_title("the Cacng8 - Gria1 gap", loc="left", fontsize=9)
    axes = [fig.add_axes([0.06 + 0.235 * k, 0.1, 0.19, 0.27]) for k in range(4)]
    multiples_panel(axes, robustness, summary, subunits)
    axes[0].text(
        0,
        1.25,
        "C.  Gene by gene, four of the variants against the primary (P9's genes "
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
            "a conclusion depends on.",
        ],
    )
    return saved(fig, save)


# ===== 10 What a whole-brain rho is made of =====

# the colour of each detail gene in figure 10 D
DETAIL_COLOURS = {
    "Cacng8": RED,
    "Gria1": DARK_BLUE,
    "Grm5": NANO,
    "Dlg2": DARK_GREY,
    "Aqp4": MID_GREY,
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
    for g in named:
        if g in x.index:
            ax.annotate(
                g,
                (x[g], y[g]),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=7,
                color=gene_colour(g, subunits),
            )
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
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
        colour = DETAIL_COLOURS.get(gene, DARK_GREY)
        for i, division in enumerate(divisions):
            if division in mine.index:
                ax.scatter(
                    mine.loc[division, "rho"],
                    i + offsets[k],
                    s=6 + 1.2 * mine.loc[division, "n_structures"],
                    color=colour,
                    linewidths=0,
                    alpha=0.9,
                )
        ax.scatter(
            means.get(gene, np.nan),
            len(divisions) + offsets[k],
            s=40,
            marker="D",
            color=colour,
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


def top_within_panel(ax: plt.Axes, table: pd.DataFrame, n_shown: int = 15) -> None:
    """F: the genes highest inside divisions, with their whole-brain rho."""
    ax.axis("off")
    ranked = table[table["rho_within"].notna()].sort_values("rho_within", ascending=False)
    lines = ["gene       within  whole"]
    lines += [
        f"{s:9s}  {r['rho_within']:+.2f}   {r['rho']:+.2f}{' *' if r['p9_gene'] else ''}"
        for s, r in ranked.head(n_shown).iterrows()
    ]
    ax.text(0, 1, "\n".join(lines), va="top", fontsize=7.5, family="monospace")
    ax.set_title(
        "F.  The genes highest inside divisions; * P9's genes", loc="left", fontsize=9
    )


def plot_between_within(
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
    """Figure 10: what a whole-brain rho is made of, between and within divisions.

    `within` and `detail` are the tables of run_ish_divisions, `map_values` the map
    on the declared structures and `coarse` its division-only version, `example`
    the profile of the gene of panel A, `calibration` the within rho of random maps.
    """
    table = within.set_index("symbol")
    p9 = table["p9_gene"]
    have = table["rho_within"].notna()
    agree_div = spearmanr(table["rho"], table["rho_division_only"]).statistic
    agree_within = spearmanr(
        table.loc[have, "rho"], table.loc[have, "rho_within"]
    ).statistic
    n_pass = int((table["q_all_spatial"] < q).sum())
    n_pass_p9 = int((table.loc[p9, "q_p9_spatial"] < q).sum())
    fig = plt.figure(figsize=(16, 12.5))
    heading(
        fig,
        "between_within",
        f"{len(table)} genes; median rho {table['rho'].median():+.2f} over the whole "
        f"brain, {table['rho_within'].median():+.2f} inside divisions; {n_pass} genes "
        f"past the within null at BH q < {q} within all genes, {n_pass_p9} within "
        "P9's",
    )

    # A: the example gene against the real map and against the division-only map
    groups = group_of(set_table)
    shared = [s for s in map_values.index if s in example]
    y = np.array([example[s] for s in shared])
    g = [groups[s] for s in shared]
    ax0 = fig.add_axes([0.05, 0.58, 0.15, 0.27])
    ax1 = fig.add_axes([0.22, 0.58, 0.15, 0.27])
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

    # B and C: every gene
    ax_b = fig.add_axes([0.44, 0.55, 0.24, 0.32])
    named_b = ["Cacng8", "Gria1", "Aqp4"]
    genes_scatter(ax_b, table["rho"], table["rho_division_only"], p9, subunits, named_b)
    ax_b.set_ylabel("rho with the division-only map")
    panel_title(
        ax_b,
        "B",
        "Most of a whole-brain rho is the contrast between divisions",
        f"the two gene orders agree at rho {agree_div:.2f} ({len(table)} genes)",
    )
    ax_c = fig.add_axes([0.74, 0.55, 0.24, 0.32])
    band = (table["null_lo"].median(), table["null_hi"].median())
    ax_c.axhspan(*band, color=NULL_BAND, lw=0, zorder=0)
    top_within = list(table.loc[have, "rho_within"].nlargest(5).index)
    passed = table["q_all_spatial"] < q
    genes_scatter(
        ax_c,
        table.loc[have, "rho"],
        table.loc[have, "rho_within"],
        p9[have],
        subunits,
        list(dict.fromkeys(top_within + list(genes))),
        filled=passed[have],
    )
    ax_c.set_ylabel(f"mean rho inside divisions ({min_structures}+ structures each)")
    panel_title(
        ax_c,
        "C",
        "Inside divisions the order changes, and rho shrinks",
        f"gene orders agree at rho {agree_within:.2f}; filled: past the within\n"
        "null (BH, all genes); pale band: the median gene's null, 95%",
    )

    # D: the detail genes division by division; E: the two nulls on random maps;
    # F: the top of the within ranking
    ax_d = fig.add_axes([0.08, 0.1, 0.28, 0.35])
    detail_panel(ax_d, detail, within, genes)
    panel_title(
        ax_d,
        "D",
        "Division by division, five genes",
        "dot size: structures in the division; diamond: the weighted mean",
    )
    ax_e = fig.add_axes([0.47, 0.13, 0.14, 0.29])
    within_rates_panel(ax_e, calibration)
    panel_title(
        ax_e,
        "E",
        "Which null for the within rho",
        f"{len(calibration)} random smooth maps\nwith no relation to the map",
    )
    top_within_panel(fig.add_axes([0.72, 0.1, 0.26, 0.33]), table)
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


# ===== 08 Which kinds of genes match =====

# sets of at least this many genes get a violin outline; smaller ones only dots
VIOLIN_MIN = 15

# sets of at most this many genes have every gene named; larger ones their top three
NAME_ALL_MAX = 16


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
    for the nano map (None for the context group, which has no test).
    """
    rho = rows["rho_nano"].to_numpy(float)
    if test is not None:
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
    edges = ["0.05" if p else ("0.5" if hollow else colour) for p in p9]
    ax.scatter(
        x + jitter,
        rho,
        s=16,
        facecolors=face,
        edgecolors=edges,
        linewidths=0.6,
        zorder=3,
    )
    ax.plot([x - 0.32, x + 0.32], [np.median(rho)] * 2, color="0.1", lw=2, zorder=4)

    # the genes named: all of a small set, the top three and the lowest of a large one
    order = np.argsort(rho)[::-1]
    if len(rho) <= NAME_ALL_MAX:
        named = order
    else:
        named = np.concatenate([order[:3], order[-1:]])
    for i in named:
        ax.text(
            x + jitter[i] + 0.06,
            rho[i],
            rows["symbol"].iloc[i],
            fontsize=6,
            va="center",
            color=DARK_BLUE if colour == DARK_BLUE else "0.15",
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
            stats = "too few genes to test"
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
    """B: one contrast named in advance against its surrogate differences."""
    bins = np.linspace(-0.8, 0.8, 65)
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
    """Figure 08: the gene sets fixed in advance against the map and its null.

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
        ax_b = fig.add_axes([0.08 + 0.3 * k, 0.09, 0.24, 0.16])
        contrast_panel(ax_b, row, contrast_nulls[row["contrast"]], n_surrogates)
    fig.text(
        0.06,
        0.335,
        "B.  The contrasts named in advance, against the surrogates (red: observed)",
        fontsize=9,
    )
    ax_c = fig.add_axes([0.68, 0.06, 0.3, 0.24])
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
            "What would mean what: postsynaptic above presynaptic and glia beyond the "
            "null: the map is postsynaptic-like, as any glutamate receptor label should "
            "be (a sanity check). Localisation genes (transport, anchoring, auxiliary "
            "subunits) past their band:",
            "the genes that set surface receptor follow the map, as a surface-fraction "
            "reading predicts; a set inside its band: its genes look like the map no "
            "more than a random smooth map would.",
        ],
    )
    return saved(fig, save)


# ===== 09 Localisation genes against matched controls =====


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
    """The label-permutation null of a difference, the observed one, the detectable."""
    ax.hist(null, bins=60, color=LIGHT_GREY)
    ax.axvline(row["difference"], color=RED, lw=1.8)
    for side in (-1, 1):
        ax.axvline(side * row["detectable"], color=DARK_GREY, ls="--", lw=0.9)
    ax.set_xlabel("difference of the medians, labels permuted")
    ax.set_ylabel("permutations")
    tidy(ax)


def plot_localisation(
    table: pd.DataFrame,
    summary: pd.DataFrame,
    nulls: dict[str, np.ndarray],
    pairs: dict[str, str],
    p_spatial: float,
    n_surrogates: int,
    panel_table: pd.DataFrame | None = None,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 09: localisation genes against expression-matched postsynaptic controls.

    `table` is localisation_test.csv, `summary` localisation_summary.csv, `nulls`
    each test's label null, `pairs` the matching, `p_spatial` the matched
    difference's spatial p; `panel_table`, when given, the same genes with the
    control pool of 5 October, whose positive control panel C draws beside.
    """
    rng = np.random.default_rng(0)
    by = table.set_index("symbol")
    loc = list(table.loc[table["side"] == "localisation", "symbol"])
    ctrl = list(table.loc[table["side"] == "control", "symbol"])
    matched = sorted(set(pairs.values()))
    test = summary.set_index("test")
    main = test.loc["matched controls"]
    positive = test.loc["positive control"]
    fig = plt.figure(figsize=(16, 10))
    heading(
        fig,
        "localisation",
        f"{len(loc)} localisation genes against {len(matched)} expression-matched "
        f"controls: difference {main['difference']:+.3f}, p = {main['p_labels']:.3f} "
        f"(spatial {p_text(p_spatial, n_surrogates)}); positive control "
        f"{positive['difference']:+.3f}, p = {positive['p_labels']:.4f}",
    )
    ax = fig.add_axes([0.06, 0.5, 0.22, 0.36])
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
    ax = fig.add_axes([0.36, 0.5, 0.26, 0.36])
    label_null_panel(ax, nulls["matched controls"], main)
    panel_title(
        ax,
        "B",
        "Its null, labels permuted between the two sets",
        f"p = {main['p_labels']:.3f}; dashed: the difference detectable at p < 0.05 "
        f"(±{main['detectable']:.3f})\nthe same difference against the surrogates of "
        f"the map: spatial {p_text(p_spatial, n_surrogates)}",
    )
    ax = fig.add_axes([0.69, 0.5, 0.28, 0.36])
    groups = []
    pools = [("GO", table, positive)]
    if panel_table is not None:
        pools.append(
            ("5 Oct.", panel_table, test.loc["5 October's controls, positive control"])
        )
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
        f"{row['difference']:+.3f}, p = {row['p_labels']:.4f}"
        for label, _, row in pools
    ]
    panel_title(
        ax, "C", "The positive control: a difference that must exist", "\n".join(lines)
    )
    ax = fig.add_axes([0.06, 0.15, 0.22, 0.24])
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
    ax = fig.add_axes([0.36, 0.11, 0.6, 0.3])
    ax.axis("off")
    lines = [
        f"{r['test']:40s} {r['median_first']:+.3f} ({int(r['n_first']):3d})  "
        f"{r['median_second']:+.3f} ({int(r['n_second']):3d})  {r['difference']:+.3f}  "
        f"{r['p_labels']:.4f}"
        for _, r in summary.iterrows()
    ]
    top = table[table["side"] == "localisation"].nlargest(5, "rho_partial")
    lines += ["", "localisation genes with the highest partial rho:"]
    lines += [f"  {r['symbol']:10s} {r['rho_partial']:+.3f}" for _, r in top.iterrows()]
    header = f"{'test':40s} {'first':12s} {'second':12s} {'diff.':7s} p, labels"
    ax.text(
        0, 1, header + "\n" + "\n".join(lines), va="top", fontsize=7.5, family="monospace"
    )
    ax.set_title(
        "E.  Every test of the design (first: localisation genes, or the reproducible "
        "controls)",
        loc="left",
        fontsize=9,
    )
    footer(
        fig,
        [
            "How to read: partial rho is the rank correlation of a gene with the map "
            "once the mean rank of the four AMPA subunits has been regressed out of "
            "both; each localisation gene is paired with the control closest in",
            "expression (controls: the other postsynaptic genes of "
            f"{figure_ref('gene_sets')}; 5 Oct.: the ontology panel's "
            "postsynaptic-density genes, as on 5 October), and labels are permuted "
            "between the two sets.",
            "What would mean what: localisation above its matched controls: the genes "
            "that set surface receptor predict the part of the map receptor abundance "
            "does not, as a surface-fraction reading predicts.",
            "A positive control found and the localisation difference not: the "
            "negative is informative, any well-measured postsynaptic gene predicts the "
            "map about as well; a positive control missed: the test says nothing.",
        ],
    )
    return saved(fig, save)


# ===== 03 How much of the map receptor mRNA and synaptic density predict =====

# the steps of the variance budget: abundance in the subunits' dark blue, density in
# the pale blue of the plan, autofluorescence in its channel's yellow, what is left
# in red; control F's components in mid grey
BUDGET_COLOURS = {
    "abundance": DARK_BLUE,
    "density": DENSITY_BLUE,
    "autofluorescence": AUTO,
    "components": MID_GREY,
    "left": RED,
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
    """The nano map's ranks against a predictor's, dots by group, the identity dashed."""
    n = len(y)
    scatter_groups(ax, x, y, groups)
    ax.plot([1, n], [1, n], color=MID_GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.set_xlim(-0.05 * n, 1.05 * n)
    ax.set_ylim(-0.05 * n, 1.05 * n)
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("nano map, rank among structures")
    tidy(ax)


def budget_bar(
    ax: plt.Axes,
    y: float,
    parts: list[tuple[str, float, str]],
    height: float = 0.56,
) -> list[float]:
    """One stacked bar of shares from 0 to 1; returns where each part starts.

    `parts` holds (colour key, share, label); a label is written inside a part wide
    enough to hold it. A negative share is drawn as nothing.
    """
    starts, x = [], 0.0
    for key, share, label in parts:
        w = max(share, 0.0)
        starts.append(x)
        ax.barh(
            y,
            w,
            left=x,
            height=height,
            color=BUDGET_COLOURS[key],
            edgecolor="white",
            linewidth=1.5,
            zorder=2,
        )
        if w > 0.06 and label:
            ink = "white" if key in ("abundance", "left") else "0.1"
            ax.text(
                x + w / 2,
                y,
                label,
                ha="center",
                va="center",
                fontsize=8,
                color=ink,
                fontweight="bold",
                zorder=4,
            )
        x += w
    return starts


def budget_panel(ax: plt.Axes, n: dict) -> None:
    """D: the variance budget of the model, the April model and control F."""
    steps, april = n["steps"], n["april_steps"]
    floor = n["floor"]
    rows = [
        (
            2,
            "the model\n(Gria1-4 as four terms)",
            [
                ("abundance", steps[0], f"Gria1-4\n{steps[0]:.0%}"),
                ("density", steps[1] - steps[0], f"+ density\n{steps[1] - steps[0]:.0%}"),
                ("autofluorescence", steps[2] - steps[1], ""),
                ("left", 1 - steps[2], ""),
            ],
        ),
        (
            1,
            "April model\n(Gria1-4 averaged into one term)",
            [
                ("abundance", april[0], f"Gria1-4\n{april[0]:.0%}"),
                ("density", april[1] - april[0], f"+ density\n{april[1] - april[0]:.0%}"),
                ("autofluorescence", april[2] - april[1], ""),
                ("left", 1 - april[2], f"left\n{1 - april[2]:.0%}"),
            ],
        ),
        (
            0,
            f"control F: the whole gene table\n({n['f_k']} components of "
            f"{n['f_genes']} genes)",
            [
                ("components", n["f_share"], f"{n['f_share']:.0%}"),
                ("left", 1 - n["f_share"], f"left\n{1 - n['f_share']:.0%}"),
            ],
        ),
    ]
    for y, _, parts in rows:
        budget_bar(ax, y, parts)

    # the model's leftover: the floor hatched inside it, its label in the rest, and
    # above the bar the interval of its edge over resampled structures
    edge = steps[2]
    ax.barh(
        2,
        floor["left_median"],
        left=edge,
        height=0.56,
        fill=False,
        hatch="////",
        edgecolor="white",
        linewidth=0,
        zorder=3,
    )
    rest = edge + floor["left_median"]
    if 1 - rest > 0.05:
        ax.text(
            (rest + 1) / 2,
            2,
            f"left\n{1 - edge:.0%}",
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            fontweight="bold",
            zorder=4,
        )
    lo, hi = n["left_ci"]
    ax.plot([1 - hi, 1 - lo], [2.4, 2.4], color="0.1", lw=1.6, zorder=5)
    for x in (1 - hi, 1 - lo):
        ax.plot([x, x], [2.35, 2.45], color="0.1", lw=1.2, zorder=5)
    ax.text(
        1.0,
        2.47,
        f"bracket: where 'left' begins, 95% over\nresampled structures: {lo:.0%} to "
        f"{hi:.0%} left",
        fontsize=7.5,
        va="bottom",
        ha="right",
        color="0.1",
    )
    ax.annotate(
        f"hatched: the calibration floor, {floor['left_median']:.0%} "
        f"({floor['left_lo']:.0%} to {floor['left_hi']:.0%}): what a map made only of "
        "receptor mRNA and\ndensity leaves when predicted from other Allen "
        "experiments (panel E)",
        xy=(edge + floor["left_median"] / 2, 2.28),
        xytext=(0.3, 2.95),
        fontsize=7.5,
        color="0.2",
        arrowprops=dict(arrowstyle="-", color="0.4", lw=0.7),
        ha="left",
        va="bottom",
    )
    auto = steps[2] - steps[1]
    ax.annotate(
        f"+ autofluorescence {auto:+.1%} (alone: "
        f"{n['partition']['autofluorescence']:+.0%})",
        xy=(steps[1] + max(auto, 0) / 2, 2.28),
        xytext=(steps[1] - 0.04, 2.47),
        fontsize=7.4,
        color="0.2",
        ha="right",
        va="bottom",
        arrowprops=dict(arrowstyle="-", color="0.4", lw=0.7),
        zorder=6,
    )
    ax.set_yticks([r[0] for r in rows])
    ax.set_yticklabels([r[1] for r in rows], fontsize=8.5)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.45, 3.35)
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
    """E: the share left against the leftover's replication, known maps and nano."""
    kinds = (
        (
            "abundance and density",
            "s",
            DARK_GREY,
            "known answer: the map is receptor mRNA + density",
        ),
        ("Gria1 mRNA", "^", DARK_BLUE, "known answer: the map is one Gria1 experiment"),
    )
    for kind, marker, colour, label in kinds:
        mine = calibration[calibration["map"] == kind]
        ax.scatter(
            mine["left"],
            mine["replication"],
            s=18,
            marker=marker,
            color=colour,
            alpha=0.55,
            linewidths=0,
            zorder=2,
            label=f"{label} ({len(mine)} draws)",
        )
        ax.scatter(
            [mine["left"].median()],
            [mine["replication"].median()],
            s=90,
            marker=marker,
            facecolor="white",
            edgecolor=colour,
            linewidths=1.6,
            zorder=3,
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
        label=f"nano, the model ({n['n_structures']} structures; 95% over structures)",
    )
    ax.set_xlim(0, max(0.5, calibration["left"].max() + 0.04))
    ax.set_ylim(min(0.6, calibration["replication"].min() - 0.03), 1.0)
    ticks = np.arange(0, ax.get_xlim()[1] + 1e-9, 0.1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    ax.set_xlabel("share of the reproducible map left")
    ax.set_ylabel("the leftover's replication (half against half)")
    ax.legend(loc="lower right", fontsize=7)
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
    """Figure 03: how much of the map receptor mRNA and synaptic density predict.

    `y` is the map's ranks on `structures`, `gria1` Gria1's, `held_out` each
    structure's prediction by the density model and by the whole model from fits
    that never saw it, `residual` the model's leftover; `numbers` holds the numbers
    of adult.beyond_figures, `calibration` and `replication` the tables of
    beyond_calibration and beyond_density.
    """
    n = numbers
    left = 1 - n["steps"][2]
    fig = plt.figure(figsize=(16, 18))
    heading(
        fig,
        "beyond_budget",
        f"{n['n_structures']} grey-matter structures, {n['n_adults']} adults; the four "
        f"subunits, density and autofluorescence predict {n['steps'][2]:.0%} of the "
        f"reproducible map; {left:.0%} left ({n['left_ci'][0]:.0%} to "
        f"{n['left_ci'][1]:.0%} over structures); calibration floor "
        f"{n['floor']['left_median']:.0%}",
    )

    # A to C: the map against Gria1, against density, against the whole model
    top = 0.705
    size = 0.22
    rho_g = spearmanr(y, gria1).statistic
    ax = fig.add_axes([0.06, top, size, size * 16 / 18])
    map_scatter(ax, gria1, y, groups, "Gria1 mRNA, rank among structures")
    panel_title(
        ax,
        "A",
        "The nano map against Gria1 mRNA",
        f"rho {rho_g:+.2f}; Gria1 alone, bent, predicts "
        f"{n['partition']['gria1']:.0%} of the reproducible map",
    )
    ax.legend(
        handles=group_handles(),
        loc="upper left",
        fontsize=7,
        frameon=True,
        framealpha=0.85,
        edgecolor="none",
    )
    ax = fig.add_axes([0.385, top, size, size * 16 / 18])
    map_scatter(
        ax,
        held_out["density"],
        y,
        groups,
        "what synaptic density predicts (held out, rank units)",
    )
    panel_title(
        ax,
        "B",
        "The nano map against synaptic density",
        f"{n['n_markers']} marker genes and the first component of {n['psd_genes']} "
        f"postsynaptic-density\ngenes, bent: {n['partition']['density']:.0%} of the "
        "reproducible map",
    )
    ax = fig.add_axes([0.71, top, size, size * 16 / 18])
    map_scatter(
        ax,
        held_out["model"],
        y,
        groups,
        "what the model predicts (held out, rank units)",
    )
    order = np.argsort(residual)
    for i in list(order[:NAMED_LEFTOVER]) + list(order[-NAMED_LEFTOVER:]):
        ax.annotate(
            acronyms[i],
            (held_out["model"][i], y[i]),
            xytext=(4, 2),
            textcoords="offset points",
            fontsize=7,
            color="0.15",
        )
    panel_title(
        ax,
        "C",
        "The nano map against the whole model",
        f"Gria1-4, density and autofluorescence, bent: {n['steps'][2]:.0%}\nof the "
        "reproducible map; the spread off the diagonal is the\nleftover "
        f"({figure_ref('beyond_where')})",
    )

    # D: the budget
    ax = fig.add_axes([0.22, 0.45, 0.74, 0.17])
    budget_panel(ax, n)
    fig.text(
        0.06,
        0.64,
        "D.  The variance budget, on structures the fit has not seen: what the four "
        "subunits predict, what synaptic density adds, what is left",
        fontsize=9,
        va="bottom",
    )

    # E: the calibration; F: the replication
    ax = fig.add_axes([0.07, 0.1, 0.38, 0.26])
    calibration_panel(ax, calibration, n)
    floor, gria = n["floor"], n["gria1"]
    panel_title(
        ax,
        "E",
        "The same model on maps whose answer is known",
        f"made of receptor mRNA + density: {floor['left_median']:.0%} left "
        f"({floor['left_lo']:.0%} to {floor['left_hi']:.0%}); one Gria1 experiment: "
        f"{gria['left_median']:.0%} ({gria['left_lo']:.0%} to {gria['left_hi']:.0%});\n"
        f"nano predicted the same way: "
        f"{min(n['nano_cal'][h] for h in ('A', 'B')):.0%} to "
        f"{max(n['nano_cal'][h] for h in ('A', 'B')):.0%}",
    )
    ax = fig.add_axes([0.6, 0.1, 0.36, 0.26])
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
            "predictor enters as x, x^2 and x^3, and each share is the cross-validated "
            "R2 over the ceiling, the share of the map that two halves of the cohort "
            "reproduce.",
            "The calibration builds maps whose answer is known from half of each "
            "gene's Allen experiments, gives them ten made-up adults as noisy as ours, "
            "and predicts them from the other half of the experiments.",
            "What would mean what: a leftover well above the floor: part of the map is "
            "set by something receptor mRNA and synaptic density do not predict. A "
            "leftover no larger than a one-gene map's: its size is within what one "
            "Allen",
            "experiment disagreeing with another produces. Either way the leftover is "
            "what the model does not predict, not a measurement of the surface "
            "fraction; its replication follows from the map's reliability (dashed, F).",
        ],
    )
    return saved(fig, save)


# ===== 04 Where the leftover lives =====

# the structures shown in the bars, each way
N_LEFTOVER_BARS = 12

# the genes drawn against the leftover: the best, and these named in any case
N_LEFTOVER_GENES = 15
LEFTOVER_NAMED = ("Cacng8", "Gria1", "Dlg2")


def leftover_planes(fig: plt.Figure, grid, planes: list[dict], span: float) -> list:
    """A: map, prediction and leftover on each plane; returns the two images to key."""
    images = []
    titles = ("the nano map", "what receptor mRNA and density predict", "the leftover")
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
        [
            f"{r['acronym']}  {textwrap.shorten(r['structure'], 34, placeholder='...')}"
            for _, r in show.iterrows()
        ],
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
    """C: the genes closest to the leftover and the named ones, with their null bands."""
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
        text = f"{r['symbol']}  (rank {int(r['rank_all'])})"
        if r["in_model"]:
            text += f"  in the model: {r['in_model']}"
        labels.append(text)
    ax.set_yticklabels(labels, fontsize=7)
    for tick, (_, r) in zip(ax.get_yticklabels(), show.iterrows()):
        tick.set_color(gene_colour(r["symbol"], subunits))
        if r["q_all"] < q:
            tick.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlim(-0.6, 0.6)
    ax.set_xlabel("Spearman rho of the gene with the leftover")
    tidy(ax)


def leftover_sets(ax: plt.Axes, genes: pd.DataFrame, sets: pd.DataFrame) -> None:
    """D: the gene sets of analysis 3 against the leftover, each with its null band."""
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
            stats = f"p = {test['p_spatial']:.3f}, q = {test['q']:.2f}"
        else:
            stats = "too few to test"
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
    genes: pd.DataFrame,
    sets: pd.DataFrame,
    numbers: dict,
    n_surrogates: int,
    t_max: float,
    t_max_genes: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 04: where the leftover lives, and whether any gene's map follows it.

    `planes` holds per plane its labels and the map, prediction and leftover
    painted (adult.beyond_figures.plane_images); `residuals`, `genes` and `sets` are
    residual_by_structure.csv, leftover_genes.csv and leftover_sets.csv; `t_max`
    the t at which a structure's bar turns black, `t_max_genes` the same for a
    gene's rho over its SD across resampled adults.
    """
    n = numbers
    q = n["q"]
    fig = plt.figure(figsize=(16, 19))
    top_up = ", ".join(residuals["acronym"].head(4))
    top_down = ", ".join(residuals["acronym"].tail(4)[::-1])
    heading(
        fig,
        "beyond_where",
        f"most above prediction: {top_up}; most below: {top_down}; genes past the "
        f"leftover's null at BH q < {q}: {n['n_pass']} of {n['n_genes']}",
    )

    # A: the three planes
    span = float(np.nanmax([np.nanmax(np.abs(p["leftover"])) for p in planes]))
    grid = fig.add_gridspec(
        len(planes),
        3,
        left=0.04,
        right=0.6,
        top=0.905,
        bottom=0.47,
        hspace=0.12,
        wspace=0.04,
    )
    images = leftover_planes(fig, grid, planes, span)
    cax = fig.add_axes([0.06, 0.445, 0.3, 0.008])
    cb = fig.colorbar(images[0], cax=cax, orientation="horizontal")
    cb.ax.set_xlim(0, 1)
    cb.set_label("rank among the structures of the fit (0 low, 1 high)", fontsize=8)
    cax = fig.add_axes([0.43, 0.445, 0.15, 0.008])
    cb = fig.colorbar(images[1], cax=cax, orientation="horizontal")
    cb.set_label("leftover (ranks): red above prediction", fontsize=8)
    fig.text(
        0.04,
        0.925,
        f"A.  The map, the prediction and the leftover on {len(planes)} coronal "
        "planes\nflat grey: structures the fit does not use",
        fontsize=9,
        va="bottom",
    )

    # B: the structures with the largest leftovers
    ax = fig.add_axes([0.75, 0.47, 0.22, 0.435])
    leftover_bars(ax, residuals, t_max)
    panel_title(
        ax,
        "B",
        f"The {N_LEFTOVER_BARS} largest leftovers each way",
        f"grey: t of the ten adults' own leftovers (black at {t_max:g})",
    )

    # C: the genes against the leftover; D: the gene sets
    subunits = {
        s for s, t in zip(genes["symbol"], genes["gene_sets"]) if "subunits" in str(t)
    }
    ax = fig.add_axes([0.2, 0.08, 0.25, 0.29])
    leftover_gene_bars(ax, genes, subunits, t_max_genes, q)
    cacng8 = genes.set_index("symbol").loc["Cacng8"]
    panel_title(
        ax,
        "C",
        "The genes closest to the leftover",
        f"pale blue: 95% of rho with the leftover's {n_surrogates} surrogates;\n"
        f"grey: rho over its SD across resampled adults (black at {t_max_genes:g});\n"
        f"{n['n_pass']} of {n['n_genes']} genes past the null at BH q < {q}; Cacng8 "
        f"{cacng8['rho']:+.2f} ({p_text(cacng8['p_spatial'], n_surrogates)})",
    )
    ax = fig.add_axes([0.56, 0.08, 0.41, 0.29])
    leftover_sets(ax, genes, sets)
    panel_title(
        ax,
        "D",
        "The gene sets of analysis 3 against the leftover",
        "pale blue: 95% of the set's median rho with the surrogates; spatial p, and "
        "q by BH over the tested sets",
    )
    footer(
        fig,
        [
            "How to read: the leftover is the nano rank minus the rank receptor mRNA, "
            "synaptic density and autofluorescence predict (in-sample fit); a gene's rho "
            "with it is tested against surrogates with the leftover's own smoothness.",
            "Genes in the model (the subunits, the markers, the density component's "
            "genes) correlate with the leftover near zero by construction. None of the "
            "tests on the leftover was named before it was seen.",
            "What would mean what: a gene or a set past its band follows what receptor "
            "mRNA and density leave, a lead on what the leftover is; nothing past its "
            "band: no map in the gene table looks like the leftover.",
            "A claim about one structure needs its own null; the bars say where the "
            "leftover is largest and how steady it is across the adults, not that a "
            "structure stands out.",
        ],
    )
    return saved(fig, save)


# ===== 13 What the green channel reports =====

# the three channels as figure 13 names and draws them: name, what it records, the
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
        ax.text(x + 0.26, v.mean(), format(v.mean(), fmt), fontsize=7.5, va="center")
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


def plot_green_channel(
    rows: pd.DataFrame,
    images: dict[str, np.ndarray],
    lab: np.ndarray,
    mouse: str,
    plane: int,
    n_structures: int,
    gene: str,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 13: whether the green channel reports the tag or the tissue.

    `rows` is sep_channel_check.csv (one row per adult), `images` the raw channels
    (sig, sep, auto) of `mouse` on CCF plane `plane` with its labels `lab`;
    `gene` is the gene the channels are compared with (Gria1).
    """
    fig = plt.figure(figsize=(16, 16))
    sep_auto = rows["rho_sep_auto"]
    nano_auto = rows["rho_nano_auto"]
    heading(
        fig,
        "green_channel",
        f"{len(rows)} adults, {n_structures} declared structures: SEP follows "
        f"autofluorescence at {sep_auto.min():.2f} to {sep_auto.max():.2f} in every "
        f"adult, nano at {nano_auto.min():.2f} to {nano_auto.max():.2f}",
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
    paired_columns(
        ax,
        rows,
        [
            ("SEP with\nautofluorescence", "rho_sep_auto", SEP_DOT),
            ("SEP with\nnano", "rho_sep_nano", SEP_DOT),
            ("nano with\nautofluorescence", "rho_nano_auto", NANO_DOT),
        ],
    )
    ax.set_xticks(range(3))
    ax.set_xticklabels(
        ["SEP with\nautofluorescence", "SEP with\nnano", "nano with\nautofluorescence"]
    )
    ax.axhline(0, color="0.8", lw=0.7, zorder=0)
    ax.set_xlim(-0.6, 2.8)
    ax.set_ylim(
        min(
            -0.05,
            float(rows[["rho_sep_auto", "rho_sep_nano", "rho_nano_auto"]].min().min())
            - 0.05,
        ),
        1.0,
    )
    ax.set_ylabel("Spearman rho across structures, per adult")
    panel_title(
        ax,
        "C",
        "Whom each channel follows, adult by adult",
        "one dot per adult, lines join the same adult; bar: the mean",
    )
    tidy(ax)

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
    paired_columns(
        ax,
        rows,
        [
            ("nano", "range_nano", NANO_DOT),
            ("autofluorescence", "range_auto", AUTO_DOT),
            ("SEP", "range_sep", SEP_DOT),
        ],
        fmt=".2f",
    )
    ax.set_xticks(range(3))
    ax.set_xticklabels(["nano", "autofluorescence", "SEP"])
    ax.set_xlim(-0.6, 2.8)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("p90 - p10 across structures (log2)")
    panel_title(
        ax,
        "D",
        "How much each channel varies across the brain",
        "a channel reporting the receptor should vary about as much as nano",
    )
    tidy(ax)
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

# the columns of the overview's table: heading, width (inches), font size
OVERVIEW_COLUMNS = (
    ("step (figures)", 2.5, 9.0),
    ("question", 3.1, 8.5),
    ("what this run says", 5.9, 8.0),
    ("what stays open", 3.3, 8.0),
)

# the tint of each part of the walk: a palette colour and its alpha
PART_TINTS = {
    "inputs": (LIGHT_GREY, 0.35),
    "part 1": (NANO, 0.22),
    "part 2": (NULL_BAND, 0.75),
    "controls": (LIGHT_GREY, 0.35),
    "limit": (LIGHT_GREY, 0.35),
    "appendix": ("white", 1.0),
}

# the tint of each paragraph of the argument, by its heading
ARGUMENT_TINTS = {
    "The question": "inputs",
    "Part 1": "part 1",
    "Part 2": "part 2",
    "The limit": "limit",
}


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


def tinted_box(ax: plt.Axes, x: float, y: float, w: float, h: float, part: str):
    """A rounded box in the tint of its part, in the overview's inch coordinates."""
    colour, alpha = PART_TINTS[part]
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0,rounding_size=0.06",
            fc=colour,
            alpha=alpha,
            ec="none",
            zorder=1,
        )
    )
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0,rounding_size=0.06",
            fc="none",
            ec=MID_GREY,
            lw=0.6,
            zorder=2,
        )
    )


def overview_cells(row: dict) -> list[tuple[list[str], float, str, str]]:
    """The four cells of one row: wrapped lines, font size, weight, colour."""
    widths = [c[1] - 0.25 for c in OVERVIEW_COLUMNS]
    sizes = [c[2] for c in OVERVIEW_COLUMNS]
    title = wrapped(row["title"], widths[0], sizes[0], bold=True)
    word = "figures" if "," in row["figures"] else "figure"
    title += ["", f"{word} {row['figures']}"]
    numbers = []
    for line in row["lines"]:
        part = wrapped(line, widths[2] - 0.12, sizes[2])
        numbers += ["- " + part[0]] + ["  " + p for p in part[1:]]
    return [
        (title, sizes[0], "bold", "0.1"),
        (wrapped(row["question"], widths[1], sizes[1]), sizes[1], "normal", "0.1"),
        (numbers, sizes[2], "normal", "0.1"),
        (wrapped(row["open"], widths[3], sizes[3]), sizes[3], "normal", DARK_GREY),
    ]


def argument_block(ax: plt.Axes, top: float, argument: list[tuple[str, str]]) -> float:
    """The argument at the top of the overview; returns the inch below it."""
    width = OVERVIEW_WIDTH - 2 * OVERVIEW_MARGIN
    size = 9.5
    y = top
    for heading_text, body in argument:
        lines = wrapped(body, width - 1.9, size)
        h = len(lines) * line_height(size) + 0.18
        tinted_box(ax, OVERVIEW_MARGIN, y, width, h, ARGUMENT_TINTS[heading_text])
        ax.text(
            OVERVIEW_MARGIN + 0.15,
            y + 0.1,
            heading_text,
            fontsize=10,
            fontweight="bold",
            va="top",
            zorder=3,
        )
        ax.text(
            OVERVIEW_MARGIN + 1.7,
            y + 0.1,
            "\n".join(lines),
            fontsize=size,
            va="top",
            linespacing=1.3,
            zorder=3,
        )
        y += h + 0.07
    return y


def table_block(ax: plt.Axes, top: float, rows: list[dict]) -> float:
    """The overview's table, a row per step of the walk; returns the inch below it."""
    x = OVERVIEW_MARGIN
    for heading_text, width, _ in OVERVIEW_COLUMNS:
        ax.text(x + 0.1, top, heading_text, fontsize=9.5, fontweight="bold", va="top")
        x += width + 0.1
    y = top + 0.32
    for row in rows:
        cells = overview_cells(row)
        h = max(len(c[0]) * line_height(c[1]) for c in cells) + 0.22
        x = OVERVIEW_MARGIN
        for (lines, size, weight, colour), (_, width, _) in zip(cells, OVERVIEW_COLUMNS):
            tinted_box(ax, x, y, width, h, row["part"])
            ax.text(
                x + 0.12,
                y + 0.11,
                "\n".join(lines),
                fontsize=size,
                fontweight=weight,
                color=colour,
                va="top",
                linespacing=1.3,
                zorder=3,
            )
            x += width + 0.1
        y += h + 0.08
    return y


def items_block(ax: plt.Axes, top: float, items: tuple) -> float:
    """The A-items the ISH line touches, in two columns; returns the inch below them."""
    ax.text(
        OVERVIEW_MARGIN,
        top,
        "The open analysis items of docs/REFACTOR_COVERAGE.md that the ISH line "
        "touches, and where they stand",
        fontsize=9.5,
        fontweight="bold",
        va="top",
    )
    half = (len(items) + 1) // 2
    width = (OVERVIEW_WIDTH - 2 * OVERVIEW_MARGIN) / 2
    for k, (item, text) in enumerate(items):
        x = OVERVIEW_MARGIN + (k // half) * width
        y = top + 0.32 + (k % half) * line_height(8.5)
        ax.text(x, y, f"{item:<4} {text}", fontsize=8, family="monospace", va="top")
    return top + 0.32 + half * line_height(8.5)


def plot_overview(
    argument: list[tuple[str, str]],
    rows: list[dict],
    items: tuple,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 00: the question, the argument, and which figure answers what.

    `argument` holds (heading, text) of the box at the top, `rows` the table's rows
    (ish.overview.overview_rows), `items` the A-items and their state. The figure is
    as tall as its text: every size is worked out in inches first.
    """
    # the height, from a dry run of the layout on a throwaway axes
    probe = plt.figure(figsize=(OVERVIEW_WIDTH, 1))
    ax = probe.add_axes([0, 0, 1, 1])
    bottom = items_block(
        ax, table_block(ax, argument_block(ax, 1.0, argument), rows) + 0.25, items
    )
    plt.close(probe)
    height = bottom + 1.0

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
        fontsize=12,
    )
    ax.text(
        OVERVIEW_WIDTH / 2,
        0.62,
        "every number below is this run's (tables/numbers_for_the_text.csv); the "
        "verdict in words is in docs/ISH_ANALYSIS.md",
        ha="center",
        va="top",
        fontsize=9,
        color=DARK_GREY,
    )
    y = argument_block(ax, 1.0, argument)
    y = table_block(ax, y + 0.25, rows)
    y = items_block(ax, y + 0.25, items)
    ax.text(
        OVERVIEW_MARGIN,
        y + 0.3,
        "How to read: one row per step of the walk, in the figures' order; the tint "
        "says which part of the argument a row serves (orange: part 1, blue: part 2, "
        "grey: inputs, controls and the limit).\nWhat would mean what: part 1 stands "
        "when the leftover is well above the calibration floor and replicates; part 2 "
        "corroborates the surface-fraction reading when the surface-regulating genes "
        "pass their nulls and their matched controls.",
        fontsize=8,
        color=NOTE_GREY,
        va="top",
        linespacing=1.5,
    )
    return saved(fig, save)


# ===== 14 April's headline, then and now =====

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


def headline_scatter(ax: plt.Axes, headline: pd.DataFrame) -> float:
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
        "B",
        f"Gene by gene ({len(both)} genes in both)",
        f"the order of the genes agrees at rho {agreement:.2f}; subunits dark blue",
    )
    tidy(ax)
    return agreement


def anova_panel(ax: plt.Axes, anova: pd.DataFrame, f_info: dict) -> None:
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
        "C",
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
        if g.tested:
            labels.append(
                f"{g.label}\nn = {g.n_genes}\n{p_text(g.p_spatial, n_surrogates)}, "
                f"q = {g.q:.3f}"
            )
        else:
            labels.append(f"{g.label}\nn = {g.n_genes}\ntoo few genes to test")
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
    order: tuple,
    n_surrogates: int,
    q: float,
    save: Path | None = None,
) -> plt.Figure:
    """Figure 14: April's category violins, recomputed today and against the null.

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
        "april_headline",
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
