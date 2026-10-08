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
from matplotlib.patches import Patch, Rectangle
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

# the guided figures and their place in the walk; a figure's title, its file name and
# the figures' references to each other all read this table, so the walk is
# renumbered here alone
FIGURES = {
    "structures": 1,
    "genes": 2,
    "one_comparison": 3,
    "spatial_null": 4,
    "gene_ranking": 5,
    "autofluorescence": 6,
    "robustness": 7,
    "between_within": 8,
    "gene_sets": 9,
    "localisation": 10,
}


# ===== Shared pieces =====


def figure_file(key: str) -> str:
    """The file name of a guided figure: its number and its key, 05_gene_ranking.png."""
    return f"{FIGURES[key]:02d}_{key}.png"


def figure_ref(key: str) -> str:
    """How a figure names another: 'figure 04'."""
    return f"figure {FIGURES[key]:02d}"


def heading(fig: plt.Figure, key: str, question: str, numbers: str) -> None:
    """The figure's number and question, and under it the line of computed numbers."""
    title = f"{FIGURES[key]}.  {question}"
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
        "Which genes, which Allen experiments, and how trustworthy is each map?",
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

# the structures named on the scatters of figure 03, by acronym: the top of the map,
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


# ===== 03 What one comparison is =====


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
    """Figure 03: what one gene's rho with the map is, for nano and a few genes.

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
        "What does a gene's rho with the map mean, concretely?",
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


# ===== 04 The spatial null =====


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
    """Figure 04: why a spatial null, and the null itself.

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
        "How large a rho do unrelated smooth maps give, and do the surrogates have "
        "the nano map's smoothness?",
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


# ===== 05 Which expression maps look like the nano map =====


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
    """Figure 05: P9's genes ranked against the map, with the null; the gap; counts.

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
        "Which genes' maps order the structures as the nano map does, beyond what a "
        "map with its smoothness would?",
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


# ===== 06 Is the ranking the nanobody's or the tissue's =====

# genes always named on the scatter of figure 06, beside those far from the diagonal
NAMED_AUTO = ("Cacng8", "Gria1", "Aqp4")

# the gene lists of figure 06 D: lines per column, and the width of a column in the
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
    """Figure 06: the same ranking on the autofluorescence map of the same sections.

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
        "Would the tissue's own autofluorescence, in the same sections, give the same "
        "gene ranking?",
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
