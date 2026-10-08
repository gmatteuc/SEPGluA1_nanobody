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
from scipy.stats import rankdata

from sepmap.ish.section_qc import ISH_QC, SECTION_AXIS
from sepmap.plotting import (
    DARK_BLUE,
    DARK_GREY,
    DENSITY_BLUE,
    DIVISION_GROUP,
    DIVISION_GROUP_COLOURS,
    GROUP_COLOURS,
    MID_GREY,
    NO_DATA_GREY,
    PAIR_LINE,
    RED,
    SET_COLOURS,
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


# ===== Shared pieces =====


def heading(fig: plt.Figure, number: int, question: str, numbers: str) -> None:
    """The figure's number and question, and under it the line of computed numbers."""
    fig.text(0.5, 0.985, f"{number}.  {question}", ha="center", va="top", fontsize=12)
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
        1,
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
        2,
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
