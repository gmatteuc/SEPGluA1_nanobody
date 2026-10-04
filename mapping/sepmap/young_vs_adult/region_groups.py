"""Region statistics pooled by cortical system and by layer, each brain on its own atlas.

The per-mouse measurements of young_vs_adult.region_plot, pooled to ask coarser
questions than one area at a time:

    by system   primary and higher-order visual, somatosensory and auditory, then
                frontal, motor, retrosplenial, lateral/insular and the subcortical
                divisions. A shift spread thinly over many small areas shows,
                instead of being split into pieces that are each not significant,
                and the distinction between primary and higher-order areas that
                the critical-period argument rests on is explicit.
    by layer    within each cortical system, supragranular (1, 2/3), granular (4)
                and infragranular (5, 6a, 6b), read from the ontology's
                substructure names: the laminar pattern seen by eye in the videos,
                measured.

Nothing is warped: a group mean is the voxel-weighted mean over the group's labels
in that brain, a ratio's over the voxels where its smoothed reference is positive
(region_plot.label_sums). The readings and their references are those of
region_plot: ratio, nano per autofluorescence; sepratio, nano per unit SEP, which
was meant as surface receptor per unit receptor expressed and is not, the green
channel being mostly autofluorescence here (adult.sep_channel_check); cref and
subref, relative to the brain's isocortex and to its subcortex (TH, HY, PAL, MB, P
and MY: no isocortex, OLF, HPF, CTXsp, STR, CB, fiber tracts, ventricles or
unassigned labels); zref, range-matched to each brain's own spread. So is the test,
Mann-Whitney of the young group against the adults, uncorrected in the figure, with
BH q-values in the CSV.

Writes into comparisons_v2/young_vs_adult/:

    group_stats.csv   one row per reading, grouping and group: the medians, both
                      tests, the P20-only contrast and the naive-rws null
    group_plot.png    the systems, one dot per mouse
    laminar_plot.png  the layers within each cortical system

Run by run_region_groups.py.
"""

import csv
import math
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from sepmap.config import SETTINGS
from sepmap.plotting import GROUP_COLOURS, save_figure
from sepmap.volumes.cohort import (
    NAIVE,
    RWS,
    SIGNED_READINGS,
    YOUNG_P20,
)
from sepmap.volumes.per_mouse import CSV_MAP, DATA, MICE, annotation_20
from sepmap.young_vs_adult.region_plot import (
    ADULTS,
    GROUPS,
    LABEL,
    NOT_SUBCORTEX,
    READINGS,
    bh_fdr,
    label_sums,
    mannwhitney,
    welch,
)

# the smallest group (and structure) kept, and the brains a group needs to be tested
REGION_TABLES = SETTINGS["region_tables"]
REGION_GROUPS = SETTINGS["region_groups"]

OUT = DATA / "comparisons_v2" / "young_vs_adult"

# cortical systems, primary apart from higher order: a thalamorecipient primary area
# and its higher-order neighbours mature on different schedules, which the
# critical-period argument rests on. A sensory system names in brackets what it holds.
SYSTEMS = [
    ("primary visual (VISp)", lambda a: a == "VISp"),
    # RL and AL sit on the visual-somatosensory border and look modulated in the
    # maps: a group of their own rather than averaged into the belt around them
    ("associative VT (RL+AL)", lambda a: a in ("VISrl", "VISal")),
    # every other visual area; VISC, the visceral area, is not one and goes to
    # lateral/insular
    (
        "higher visual (all other)",
        lambda a: a.startswith("VIS") and a not in ("VISp", "VISC", "VISrl", "VISal"),
    ),
    ("primary somatosensory (SSp)", lambda a: a.startswith("SSp")),
    ("higher somatosensory (SSs)", lambda a: a == "SSs"),
    ("primary auditory (AUDp)", lambda a: a == "AUDp"),
    ("higher auditory (AUDd/po/v)", lambda a: a.startswith("AUD") and a != "AUDp"),
    # anterior cingulate, orbital, prelimbic, infralimbic, frontal pole and dorsal
    # peduncular areas
    (
        "frontal",
        lambda a: a.startswith(("ACA", "ORB")) or a in ("PL", "ILA", "FRP", "DP"),
    ),
    ("motor", lambda a: a in ("MOp", "MOs")),
    ("retrosplenial", lambda a: a.startswith("RSP")),
    # agranular insular, gustatory, visceral, temporal association, perirhinal and
    # ectorhinal areas
    (
        "lateral/insular",
        lambda a: a.startswith("AI") or a in ("GU", "VISC", "TEa", "PERI", "ECT"),
    ),
]

# the systems of the laminar figure, the ones the layer question is about: with
# every system on it, the figure would be unreadable
LAMINAR_SYSTEMS = (
    "primary visual (VISp)",
    "associative VT (RL+AL)",
    "higher visual (all other)",
    "primary somatosensory (SSp)",
    "higher somatosensory (SSs)",
    "primary auditory (AUDp)",
    "higher auditory (AUDd/po/v)",
    "frontal",
    "retrosplenial",
)

# the subcortical side, taken straight from the ontology's divisions
DIVISIONS = [
    ("thalamus", "TH"),
    ("striatum", "STR"),
    ("pallidum", "PAL"),
    ("hippocampus", "HPF"),
    ("hypothalamus", "HY"),
    ("midbrain", "MB"),
    ("olfactory", "OLF"),
    ("cortical subplate", "CTXsp"),
]

# the layers of each depth band, as layer_of reads them from the substructure names
LAYERS = [
    ("supragranular", ("1", "2", "3", "2/3")),
    ("granular", ("4",)),
    ("infragranular", ("5", "6", "6a", "6b")),
]


def layer_of(substructure_name: str) -> str | None:
    """'Primary visual area, layer 2/3' -> '2/3'; anything without a layer -> None."""
    # the token after 'layer': a number, an optional '/number' (2/3), an optional a or b
    m = re.search(r"layer\s*([0-9]+(?:/[0-9]+)?[ab]?)", substructure_name, re.I)
    return m.group(1) if m else None


def load_parcellation_terms() -> tuple[dict[int, str], dict[int, str], dict[int, str]]:
    """Structure and division acronyms and substructure names, by parcellation index.

    Read from the parcellation term membership table (CSV_MAP).
    """
    stru, divi, sub = {}, {}, {}
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        # a row ties an index to one term of one term set: keep the structure and
        # division acronyms, and the substructure name, which carries the layer
        for r in csv.DictReader(fh):
            i = int(r["parcellation_index"])
            if r["parcellation_term_set_name"] == "structure":
                stru[i] = r["parcellation_term_acronym"]
            elif r["parcellation_term_set_name"] == "division":
                divi[i] = r["parcellation_term_acronym"]
            elif r["parcellation_term_set_name"] == "substructure":
                sub[i] = r["parcellation_term_name"]
    return stru, divi, sub


def define_groups(
    stru: dict[int, str], divi: dict[int, str], layer: dict[int, str | None]
) -> dict[tuple[str, str], set[int]]:
    """The systems, the subcortical divisions and the layers within LAMINAR_SYSTEMS.

    Each group is a set of parcellation indices, keyed by (grouping, name); groups
    with no index are dropped.
    """
    groups = {}

    # each cortical system: the isocortical indices whose structure its rule accepts
    for name, pred in SYSTEMS:
        groups[("system", name)] = {
            i for i, a in stru.items() if divi.get(i) == "Isocortex" and pred(a)
        }

    # each subcortical division whole, under the same grouping as the systems
    for name, div in DIVISIONS:
        groups[("system", name)] = {i for i in stru if divi.get(i) == div}

    # the three depth bands within each system of the laminar figure
    for name, pred in SYSTEMS:
        if name not in LAMINAR_SYSTEMS:
            continue
        for lname, tokens in LAYERS:
            groups[("layer", f"{name} {lname}")] = {
                i
                for i, a in stru.items()
                if divi.get(i) == "Isocortex" and pred(a) and layer.get(i) in tokens
            }

    # drop the groups no index falls in
    groups = {k: v for k, v in groups.items() if v}
    return groups


def group_cells(
    groups: dict[tuple[str, str], set[int]], sums: dict[str, np.ndarray]
) -> dict[tuple[str, str], tuple | None]:
    """The voxel-weighted means of each group large enough, None for the others.

    `sums` is region_plot.label_sums' table per parcellation index; a ratio's mean
    is over the voxels that have one, NaN when none has.
    """
    nlab = len(sums["n"])
    cells = {}
    for key, ids in groups.items():
        # the group's indices this brain's atlas reaches: its sums stop at its
        # largest label
        ids = [i for i in ids if i < nlab]

        # tissue voxels, and the voxels that have each ratio
        c = sums["n"][ids].sum()
        c_rat = sums["ratio_n"][ids].sum()
        c_sep = sums["sepratio_n"][ids].sum()

        # the means, if the group is as large as the smallest structure the tables
        # keep (region_tables.min_vox20); None otherwise
        if c >= REGION_TABLES["min_vox20"]:
            cells[key] = (
                int(c),
                sums["sig"][ids].sum() / c,
                sums["ratio"][ids].sum() / c_rat if c_rat else np.nan,
                sums["sepratio"][ids].sum() / c_sep if c_sep else np.nan,
            )
        else:
            cells[key] = None
    return cells


def references(
    stru: dict[int, str],
    divi: dict[int, str],
    n: np.ndarray,
    s_sig: np.ndarray,
    nlab: int,
) -> dict[str, float]:
    """The two references: mean sig of the isocortex, and of the subcortex.

    The subcortex leaves out the divisions in NOT_SUBCORTEX.
    """
    # the indices of each, among those this brain's atlas reaches
    iso = [i for i in stru if divi.get(i) == "Isocortex" and i < nlab]
    sub_ids = [i for i in stru if divi.get(i, "") not in NOT_SUBCORTEX and i < nlab]

    # the voxel-weighted mean sig over each
    return {
        "cref": s_sig[iso].sum() / n[iso].sum(),
        "subref": s_sig[sub_ids].sum() / n[sub_ids].sum(),
    }


def structure_sig_means(
    stru: dict[int, str], n: np.ndarray, s_sig: np.ndarray
) -> dict[str, float]:
    """The brain's own distribution over structures, for the range match.

    Taken per structure so that it does not depend on the grouping; structures
    under region_tables.min_vox20 voxels are left out.
    """
    # tissue voxels and summed sig per structure, its layers pooled; label 0 (outside
    # the brain) and an index with no structure term are skipped
    by_struct = defaultdict(lambda: [0, 0.0])
    for i in np.nonzero(n)[0]:
        if i == 0 or i not in stru:
            continue
        by_struct[stru[i]][0] += int(n[i])
        by_struct[stru[i]][1] += s_sig[i]

    # the mean of each structure large enough
    return {
        k: v[1] / v[0] for k, v in by_struct.items() if v[0] >= REGION_TABLES["min_vox20"]
    }


def group_means(
    mice: list[str],
    groups: dict[tuple[str, str], set[int]],
    stru: dict[int, str],
    divi: dict[int, str],
) -> tuple[dict[str, dict], dict[str, dict[str, float]], dict[str, dict[str, float]]]:
    """Per mouse: the group means, the two references and the structure means of sig.

    Returns (per, refs, struct_mean). A group's cell is (voxels, mean sig, mean
    ratio, mean sepratio), None under region_tables.min_vox20 voxels; the structure
    means, of the structures with at least as many, feed the range match. Stops if
    a brain has no SEP channel while the sepratio reading is in force.
    """
    anns, per, refs, struct_mean = {}, {}, {}, {}
    for mouse in mice:
        # the brain's own atlas, loaded once for all the brains on it
        atlas_key = MICE[mouse][1]
        if atlas_key not in anns:
            anns[atlas_key] = annotation_20(atlas_key)

        # sums per parcellation index, then the groups, the two references and the
        # structure means
        sums = label_sums(mouse, anns[atlas_key])
        n, s_sig = sums["n"], sums["sig"]
        per[mouse] = group_cells(groups, sums)
        refs[mouse] = references(stru, divi, n, s_sig, len(n))
        struct_mean[mouse] = structure_sig_means(stru, n, s_sig)

        # one line per brain: how many groups are large enough to have a value
        print(
            f"{mouse:20s} {sum(v is not None for v in per[mouse].values())}"
            f"/{len(groups)} groups",
            flush=True,
        )
    return per, refs, struct_mean


def range_match(
    mice: list[str],
    struct_mean: dict[str, dict[str, float]],
    refs: dict[str, dict[str, float]],
) -> dict[str, tuple[float, float]]:
    """Per mouse the median and p90-p10 spread of log2 cortex-relative structure means.

    Taken over the structures every brain has, so the spread does not depend on
    which regions the sections happened to cover.
    """
    # the structures every brain has, so all spreads are taken over the same set
    common = set.intersection(*[set(struct_mean[m]) for m in mice])
    norm = {}
    for m in mice:
        # log2 of each structure's mean over the brain's isocortex mean; a mean at or
        # below the background has no log and is left out
        v = np.array(
            [
                math.log2(struct_mean[m][k] / refs[m]["cref"])
                for k in sorted(common)
                if struct_mean[m][k] > 0
            ]
        )

        # the median and the p90-p10 spread, floored so a flat brain cannot divide
        # by zero
        p10, med, p90 = np.percentile(v, [10, 50, 90])
        norm[m] = (med, max(p90 - p10, 1e-6))
    return norm


def value(
    reading: str,
    mouse: str,
    key: tuple[str, str],
    per: dict[str, dict],
    norm: dict[str, tuple[float, float]],
    refs: dict[str, dict[str, float]],
) -> float | None:
    """Value of `reading` for `mouse` in group `key` (log2 or range-matched), or None."""
    # None where the group is too small in this brain
    cell = per[mouse].get(key)
    if cell is None:
        return None
    _, m_sig, m_rat, m_sep = cell

    # a signed reading (zref): the group's log2 mean over the isocortex, minus the
    # brain's median, over its spread; None for a mean at or below the background
    if reading in SIGNED_READINGS:
        if m_sig <= 0:
            return None
        med, spread = norm[mouse]
        return (math.log2(m_sig / refs[mouse]["cref"]) - med) / spread

    # ratio and sepratio are already ratios against a channel; the rest divide sig
    # by a single number measured on this brain
    if reading in ("ratio", "sepratio"):
        v = m_rat if reading == "ratio" else m_sep
    else:
        v = m_sig / refs[mouse][reading]
    return math.log2(v) if v > 0 else None


def group_stats(
    groups: dict[tuple[str, str], set[int]],
    mice: list[str],
    per: dict[str, dict],
    norm: dict[str, tuple[float, float]],
    refs: dict[str, dict[str, float]],
) -> list[dict]:
    """Per group and reading the young-against-adult tests, as rows of group_stats.csv.

    A group is tested when at least region_groups.min_young young and min_adult
    adult brains have a value. The
    Mann-Whitney BH q is taken within each reading and grouping.
    """
    rows = []
    for key in groups:
        for reading, _ in READINGS:
            # the values of each group and subgroup, a brain without one (the group
            # too small there) left out
            v = {m: value(reading, m, key, per, norm, refs) for m in mice}
            v = {m: x for m, x in v.items() if x is not None}
            yo = [v[m] for m in GROUPS["young"] if m in v]
            y20 = [v[m] for m in YOUNG_P20 if m in v]
            ad = [v[m] for m in ADULTS if m in v]
            nv = [v[m] for m in NAIVE if m in v]
            rw = [v[m] for m in RWS if m in v]

            # tested only with enough young and adult brains that have a value
            if (
                len(yo) < REGION_GROUPS["min_young"]
                or len(ad) < REGION_GROUPS["min_adult"]
            ):
                continue

            # the P20 brains alone, only with at least three, so the median does not
            # rest on one or two
            if len(y20) >= 3:
                diff_p20only = np.median(y20) - np.median(ad)
            else:
                diff_p20only = float("nan")

            # naive minus rws: the size of an adult-adult difference, which carries no
            # developmental meaning
            if nv and rw:
                naive_minus_rws = np.median(nv) - np.median(rw)
            else:
                naive_minus_rws = float("nan")

            # the medians, the young-adult differences and both tests, then the two
            # contrasts to read them against
            rows.append(
                dict(
                    reading=reading,
                    grouping=key[0],
                    group=key[1],
                    n_young=len(yo),
                    n_adult=len(ad),
                    young_median=np.median(yo),
                    adult_median=np.median(ad),
                    diff_median=np.median(yo) - np.median(ad),
                    diff_mean=np.mean(yo) - np.mean(ad),
                    mannwhitney_p=mannwhitney(yo, ad),
                    welch_p=welch(yo, ad),
                    diff_P20only=diff_p20only,
                    naive_minus_rws=naive_minus_rws,
                )
            )

    # the Mann-Whitney q-values, within each reading and grouping
    for reading, _ in READINGS:
        for grouping in ("system", "layer"):
            idx = [
                i
                for i, r in enumerate(rows)
                if r["reading"] == reading and r["grouping"] == grouping
            ]
            for i, q in zip(idx, bh_fdr([rows[i]["mannwhitney_p"] for i in idx])):
                rows[i]["mannwhitney_q_BH"] = q
    return rows


def write_group_stats(rows: list[dict]) -> None:
    """Write group_stats.csv into OUT, floats to four decimals."""
    with open(OUT / "group_stats.csv", "w", newline="", encoding="utf-8") as fh:
        # the columns in the order of the first row, which every row shares
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(
                {k: (f"{x:.4f}" if isinstance(x, float) else x) for k, x in r.items()}
            )


def print_group_table(
    groups: dict[tuple[str, str], set[int]],
    rows: list[dict],
    star: dict[tuple[str, str, str], str],
) -> None:
    """Print the median log2 young - adult of every group and reading, with its stars."""
    # the header, a column per reading
    print(
        f"\n{'group':22s} "
        + " ".join(f"{r:>12s}" for r, _ in READINGS)
        + "   (log2 young - adult, medians)"
    )

    # the systems, then the layers, one line per group
    for grouping in ("system", "layer"):
        print(f"--- by {grouping}")
        for key in [k for k in groups if k[0] == grouping]:
            # one cell per reading
            cells = []
            for reading, _ in READINGS:
                # the row of this group and reading, None if it was not tested
                r = next(
                    (
                        r
                        for r in rows
                        if r["reading"] == reading
                        and r["group"] == key[1]
                        and r["grouping"] == grouping
                    ),
                    None,
                )

                # its median difference and stars, or -- for a group not tested
                if r:
                    cells.append(
                        f"{r['diff_median']:+8.2f}{star[(reading, grouping, key[1])]:<3s}"
                    )
                else:
                    cells.append(f"{'--':>11s}")
            print(f"  {key[1]:22s} " + " ".join(cells))


def draw_group_dots(
    ax: plt.Axes,
    reading: str,
    g: str,
    keys: list[tuple[str, str]],
    per: dict[str, dict],
    norm: dict[str, tuple[float, float]],
    refs: dict[str, dict[str, float]],
) -> None:
    """Draw group `g`'s mice side by side for each of `keys`, its median as a bar."""
    # the mice spread across a column: the young right of each tick, the two adult
    # groups overlapping left of it, told apart by colour
    ms = GROUPS[g]
    jit = np.linspace(-0.2, 0.2, len(ms))
    xo = 0.26 if g == "young" else -0.1
    mids = []
    for i, key in enumerate(keys):
        # a dot per mouse with a value in this group, and the group's median
        ys = [value(reading, m, key, per, norm, refs) for m in ms]
        ys = [y for y in ys if y is not None]
        for j, y in enumerate(ys):
            ax.plot(
                i + jit[j] + xo,
                y,
                "o",
                ms=4.5,
                color=GROUP_COLOURS[g],
                alpha=0.9,
                mec="none",
            )
        mids.append(np.median(ys) if ys else np.nan)

    # the medians as horizontal bars, in one call so the group has one legend entry
    ax.plot(
        np.arange(len(keys)) + xo,
        mids,
        "_",
        ms=16,
        mew=2.4,
        color=GROUP_COLOURS[g],
        label=LABEL[g],
    )


def draw_stars(
    ax: plt.Axes,
    reading: str,
    keys: list[tuple[str, str]],
    star: dict[tuple[str, str, str], str],
) -> None:
    """Draw the young-vs-adult rank-sum stars, just under the top of the panel."""
    # 14% more room above the data, so the stars do not sit on the dots
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.14 * (hi - lo))

    # each group's stars, in the young colour, centred on its tick
    lo, hi = ax.get_ylim()
    for i, key in enumerate(keys):
        st = star.get((reading, key[0], key[1]), "")
        if st:
            ax.text(
                i,
                hi - 0.04 * (hi - lo),
                st,
                ha="center",
                va="top",
                fontsize=12,
                color=GROUP_COLOURS["young"],
            )


def dotplot(
    keys: list[tuple[str, str]],
    labels: list[str],
    fname: str,
    title: str,
    star: dict[tuple[str, str, str], str],
    per: dict[str, dict],
    norm: dict[str, tuple[float, float]],
    refs: dict[str, dict[str, float]],
) -> None:
    """Dot plot of the groups `keys` into OUT/`fname`, a panel per reading."""
    # a panel per reading on a shared x axis, wider with more groups
    fig, axes = plt.subplots(
        len(READINGS),
        1,
        figsize=(max(9, 0.62 * len(keys) + 4), 3.6 * len(READINGS)),
        sharex=True,
    )
    for ax, (reading, rtitle) in zip(axes, READINGS):
        # each group's mice side by side, and its median as a bar
        for g in ("naive", "rws", "young"):
            draw_group_dots(ax, reading, g, keys, per, norm, refs)

        # stars for the young-vs-adult rank-sum test, just under the top of the panel
        draw_stars(ax, reading, keys, star)

        # zero: equal to the reading's reference (for zref, the brain's median)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(rtitle, fontsize=10.5, loc="left")

        # the y label: what each reading is relative to
        ax.set_ylabel(
            {
                "ratio": "log2  nano / auto",
                "sepratio": "log2  nano / SEP",
                "cref": "log2  vs own isocortex",
                "subref": "log2  vs subcortex",
                "zref": "range-matched",
            }[reading],
            fontsize=10,
        )
        ax.grid(axis="y", lw=0.3, alpha=0.6)
        ax.set_xlim(-0.7, len(keys) - 0.3)

    # the legend once, on the top panel; the group names under the bottom one
    axes[0].legend(
        loc="lower left",
        fontsize=9,
        frameon=True,
        framealpha=0.9,
        edgecolor="none",
        ncol=3,
    )
    axes[-1].set_xticks(range(len(keys)))
    axes[-1].set_xticklabels(labels, rotation=55, ha="right", fontsize=9)
    fig.suptitle(title, fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    save_figure(fig, OUT / fname, dpi=105)
    plt.close(fig)


def plot_groups(
    groups: dict[tuple[str, str], set[int]],
    star: dict[tuple[str, str, str], str],
    per: dict[str, dict],
    norm: dict[str, tuple[float, float]],
    refs: dict[str, dict[str, float]],
) -> None:
    """Draw group_plot.png (the systems) and laminar_plot.png (the layers)."""
    # the systems, named on the axis as in the table
    sys_keys = [k for k in groups if k[0] == "system"]
    dotplot(
        sys_keys,
        [k[1] for k in sys_keys],
        "group_plot.png",
        "Young vs adult by system, one dot per mouse, "
        "each brain on the atlas of its own age\n"
        "Bars are group medians.  * p<0.05, ** p<0.01 "
        f"(Mann-Whitney {len(GROUPS['young'])} vs {len(ADULTS)}, uncorrected)",
        star,
        per,
        norm,
        refs,
    )

    # the layers, their bands shortened to layer numbers on the axis
    lay_keys = [k for k in groups if k[0] == "layer"]
    dotplot(
        lay_keys,
        [
            k[1]
            .replace(" supragranular", " L1-3")
            .replace(" granular", " L4")
            .replace(" infragranular", " L5-6")
            for k in lay_keys
        ],
        "laminar_plot.png",
        "Young vs adult by cortical layer within each system "
        "(layers from the ontology)\n"
        "Bars are group medians.  * p<0.05, ** p<0.01 "
        f"(Mann-Whitney {len(GROUPS['young'])} vs {len(ADULTS)}, uncorrected)",
        star,
        per,
        norm,
        refs,
    )


def main() -> None:
    """Measure the group means per mouse, test young against adult, write and draw."""
    # index -> (structure acronym, division acronym, layer or None)
    stru, divi, sub = load_parcellation_terms()
    layer = {i: layer_of(nm) for i, nm in sub.items()}

    # every group is a set of parcellation indices
    groups = define_groups(stru, divi, layer)

    # per mouse: the voxel-weighted mean of sig and of the ratios in every group,
    # plus the two references, all on that brain's own atlas
    mice = [m for g in GROUPS.values() for m in g]
    per, refs, struct_mean = group_means(mice, groups, stru, divi)

    # range match of each brain
    norm = range_match(mice, struct_mean, refs)

    # young-against-adult tests, the table, and the stars of the figures
    rows = group_stats(groups, mice, per, norm, refs)
    write_group_stats(rows)

    # two stars under p < 0.01, one under 0.05, Mann-Whitney uncorrected as the
    # figure titles say
    star = {}
    for r in rows:
        if r["mannwhitney_p"] < 0.01:
            mark = "**"
        elif r["mannwhitney_p"] < 0.05:
            mark = "*"
        else:
            mark = ""
        star[(r["reading"], r["grouping"], r["group"])] = mark

    # summary table, printed
    print_group_table(groups, rows, star)

    # figures
    plot_groups(groups, star, per, norm, refs)
    print("\nwrote group_plot.png, laminar_plot.png, group_stats.csv")
