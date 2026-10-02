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
in that brain. The readings and their references are those of region_plot: ratio,
nano per autofluorescence; sepratio, nano per unit SEP, which was meant as surface
receptor per unit receptor expressed and is not, the green channel being mostly
autofluorescence here (adult.sep_channel_check); cref and subref, relative to the
brain's isocortex and to its subcortex without HPF and STR; zref, range-matched to
each brain's own spread. So is the test, Mann-Whitney of the young group against
the adults, uncorrected in the figure, with BH q-values in the CSV.

Writes into comparisons_v2/young_vs_adult/:

    group_stats.csv   one row per reading, grouping and group: the medians, both
                      tests, the P20-only contrast and the naive-rws null
    group_plot.png    the systems, one dot per mouse
    laminar_plot.png  the layers within each cortical system

Run by run_region_groups.py.
"""

import csv
import math
import os
import re
from collections import defaultdict

import matplotlib
import numpy as np
from scipy.ndimage import gaussian_filter

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sepmap.volumes.cohort import (
    NAIVE,
    RATIO_CLIP,
    RWS,
    SIGNED_READINGS,
    YOUNG_P16,
    YOUNG_P20,
    YOUNG_P22,
)
from sepmap.volumes.per_mouse import CSV_MAP, DATA, MICE, annotation_20
from sepmap.volumes.per_mouse import OUT as PER_MOUSE
from sepmap.young_vs_adult.region_plot import (
    COL,
    NOT_SUBCORTEX,
    READINGS,
    bh_fdr,
    mannwhitney,
    welch,
)

OUT = os.path.join(DATA, "comparisons_v2", "young_vs_adult")

# the groups, as in region_plot: every young brain (P16, P20, P22) against the
# adults, naive and rws
GROUPS = {"young": YOUNG_P20 + YOUNG_P16 + YOUNG_P22, "naive": NAIVE, "rws": RWS}
ADULTS = NAIVE + RWS

# legend labels, with the counts computed from the lists
LABEL = {
    "young": f"young P16-P22 (n = {len(YOUNG_P20) + len(YOUNG_P16) + len(YOUNG_P22)})",
    "naive": f"adult naive (n = {len(NAIVE)})",
    "rws": f"adult rws (n = {len(RWS)})",
}

# cortical systems, primary apart from higher order: a thalamorecipient primary area
# and its higher-order neighbours mature on different schedules, which the
# critical-period argument rests on. A sensory system names in brackets what it holds.
SYSTEMS = [
    ("primary visual (VISp)", lambda a: a == "VISp"),
    # RL and AL sit on the visual-somatosensory border and look modulated in the
    # maps: a group of their own rather than averaged into the belt around them
    ("associative VT (RL+AL)", lambda a: a in ("VISrl", "VISal")),
    (
        "higher visual (all other)",
        lambda a: a.startswith("VIS") and a not in ("VISp", "VISC", "VISrl", "VISal"),
    ),
    ("primary somatosensory (SSp)", lambda a: a.startswith("SSp")),
    ("higher somatosensory (SSs)", lambda a: a == "SSs"),
    ("primary auditory (AUDp)", lambda a: a == "AUDp"),
    ("higher auditory (AUDd/po/v)", lambda a: a.startswith("AUD") and a != "AUDp"),
    (
        "frontal",
        lambda a: a.startswith(("ACA", "ORB")) or a in ("PL", "ILA", "FRP", "DP"),
    ),
    ("motor", lambda a: a in ("MOp", "MOs")),
    ("retrosplenial", lambda a: a.startswith("RSP")),
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


def save_figure(fig, path):
    """Save `fig` as a PNG at `path` and as an EPS beside it.

    Windows refuses to overwrite a PNG that an image viewer holds open. The
    figure then goes to <name>_new.png, with a note, so a run that writes
    several figures does not lose the rest because one of them was being
    looked at.

    The EPS is what goes into a figure for a paper. PostScript has no
    transparency, so the image layers are rasterised and composited by Agg
    first; otherwise a no-data region, transparent here, would come out opaque
    black instead of showing the ground beneath it. Text, lines and axes stay
    vector, the part that has to be editable.
    """
    try:
        fig.savefig(path, dpi=105)
    except OSError:
        alt = path.replace(".png", "_new.png")
        fig.savefig(alt, dpi=105)
        print(
            f"  NOTE: {os.path.basename(path)} is open elsewhere; "
            f"wrote {os.path.basename(alt)} instead",
            flush=True,
        )

    # the EPS, with the image layers rasterised
    eps = os.path.splitext(path)[0] + ".eps"
    for ax in fig.axes:
        for im in ax.images:
            im.set_rasterized(True)
    try:
        fig.savefig(eps, dpi=105, facecolor=fig.get_facecolor(), format="eps")
    except OSError:
        print(
            f"  NOTE: {os.path.basename(eps)} is open elsewhere; "
            "the PNG was still written",
            flush=True,
        )


def layer_of(substructure_name):
    """'Primary visual area, layer 2/3' -> '2/3'; anything without a layer -> None."""
    m = re.search(r"layer\s*([0-9]+(?:/[0-9]+)?[ab]?)", substructure_name, re.I)
    return m.group(1) if m else None


def load_parcellation_terms():
    """Structure and division acronyms and substructure names, by parcellation index.

    Read from the parcellation term membership table (CSV_MAP).
    """
    stru, divi, sub = {}, {}, {}
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            i = int(r["parcellation_index"])
            if r["parcellation_term_set_name"] == "structure":
                stru[i] = r["parcellation_term_acronym"]
            elif r["parcellation_term_set_name"] == "division":
                divi[i] = r["parcellation_term_acronym"]
            elif r["parcellation_term_set_name"] == "substructure":
                sub[i] = r["parcellation_term_name"]
    return stru, divi, sub


def define_groups(stru, divi, layer):
    """The systems, the subcortical divisions and the layers within LAMINAR_SYSTEMS.

    Each group is a set of parcellation indices, keyed by (grouping, name); groups
    with no index are dropped.
    """
    groups = {}
    for name, pred in SYSTEMS:
        groups[("system", name)] = {
            i for i, a in stru.items() if divi.get(i) == "Isocortex" and pred(a)
        }
    for name, div in DIVISIONS:
        groups[("system", name)] = {i for i in stru if divi.get(i) == div}
    for name, pred in SYSTEMS:
        if name not in LAMINAR_SYSTEMS:
            continue
        for lname, tokens in LAYERS:
            groups[("layer", f"{name} {lname}")] = {
                i
                for i, a in stru.items()
                if divi.get(i) == "Isocortex" and pred(a) and layer.get(i) in tokens
            }
    groups = {k: v for k, v in groups.items() if v}
    return groups


def per_unit(ref, sig, tissue):
    """`sig` per unit of the channel `ref` within `tissue`, clipped to RATIO_CLIP.

    Voxel by voxel; the denominator is smoothed by one 20 um voxel and
    mask-normalised, as in volumes.cohort and young_vs_adult.region_plot.
    """
    ref_s = gaussian_filter(np.where(tissue, ref, 0), 1.0) / np.maximum(
        gaussian_filter(tissue.astype(np.float32), 1.0), 1e-3
    )
    return np.clip(
        np.where(ref_s > 0, sig / np.maximum(ref_s, 1e-3), 0), -RATIO_CLIP, RATIO_CLIP
    )


def group_means(mice, groups, stru, divi):
    """Per mouse: the group means, the two references and the structure means of sig.

    Returns (per, refs, struct_mean). A group's cell is (voxels, mean sig, mean
    ratio, mean sepratio), None under 250 voxels; the structure means, of the
    structures with at least 250 voxels, feed the range match. Stops if a brain has
    no SEP channel while the sepratio reading is in force.
    """
    anns, per, refs, struct_mean = {}, {}, {}, {}
    for mouse in mice:
        atlas_key = MICE[mouse][1]
        if atlas_key not in anns:
            anns[atlas_key] = annotation_20(atlas_key)
        ann = anns[atlas_key]
        z = np.load(os.path.join(PER_MOUSE, mouse + ".npz"))
        sig = z["sig"].astype(np.float32)
        auto = z["auto"].astype(np.float32)
        tissue = z["tissue"]
        if any(r == "sepratio" for r, _ in READINGS) and "sep" not in z.files:
            raise SystemExit(
                f"{mouse}: no SEP channel in its per-mouse file. Run\n"
                f"  run_add_sep_channel.m for this brain, then run_per_mouse.py,\n"
                f"  or drop the reading with V2_READINGS."
            )

        # sums per parcellation index of the voxels, sig and the two ratios
        ratio = per_unit(auto, sig, tissue)
        sepratio = (
            per_unit(z["sep"].astype(np.float32), sig, tissue)
            if "sep" in z.files
            else np.zeros_like(sig)
        )
        lab = ann[tissue]
        nlab = int(ann.max()) + 1
        n = np.bincount(lab, minlength=nlab)
        s_sig = np.bincount(lab, weights=sig[tissue], minlength=nlab)
        s_rat = np.bincount(lab, weights=ratio[tissue], minlength=nlab)
        s_sep = np.bincount(lab, weights=sepratio[tissue], minlength=nlab)

        # the voxel-weighted means of each group with at least 250 voxels
        per[mouse] = {}
        for key, ids in groups.items():
            ids = [i for i in ids if i < nlab]
            c = n[ids].sum()
            per[mouse][key] = (
                (int(c), s_sig[ids].sum() / c, s_rat[ids].sum() / c, s_sep[ids].sum() / c)
                if c >= 250
                else None
            )

        # the two references: mean sig of the isocortex, and of the subcortex
        # without the divisions in NOT_SUBCORTEX
        iso = [i for i in stru if divi.get(i) == "Isocortex" and i < nlab]
        sub_ids = [i for i in stru if divi.get(i, "") not in NOT_SUBCORTEX and i < nlab]
        refs[mouse] = {
            "cref": s_sig[iso].sum() / n[iso].sum(),
            "subref": s_sig[sub_ids].sum() / n[sub_ids].sum(),
        }

        # the brain's own distribution over structures, for the range match, taken
        # per structure so that it does not depend on the grouping
        by_struct = defaultdict(lambda: [0, 0.0])
        for i in np.nonzero(n)[0]:
            if i == 0 or i not in stru:
                continue
            by_struct[stru[i]][0] += int(n[i])
            by_struct[stru[i]][1] += s_sig[i]
        struct_mean[mouse] = {k: v[1] / v[0] for k, v in by_struct.items() if v[0] >= 250}
        print(
            f"{mouse:20s} {sum(v is not None for v in per[mouse].values())}"
            f"/{len(groups)} groups",
            flush=True,
        )
    return per, refs, struct_mean


def range_match(mice, struct_mean, refs):
    """Per mouse the median and p90-p10 spread of log2 cortex-relative structure means.

    Taken over the structures every brain has, so the spread does not depend on
    which regions the sections happened to cover.
    """
    common = set.intersection(*[set(struct_mean[m]) for m in mice])
    norm = {}
    for m in mice:
        v = np.array(
            [
                math.log2(struct_mean[m][k] / refs[m]["cref"])
                for k in sorted(common)
                if struct_mean[m][k] > 0
            ]
        )
        p10, med, p90 = np.percentile(v, [10, 50, 90])
        norm[m] = (med, max(p90 - p10, 1e-6))
    return norm


def value(reading, mouse, key, per, norm, refs):
    """Value of `reading` for `mouse` in group `key` (log2 or range-matched), or None."""
    cell = per[mouse].get(key)
    if cell is None:
        return None
    _, m_sig, m_rat, m_sep = cell
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


def group_stats(groups, mice, per, norm, refs):
    """Per group and reading the young-against-adult tests, as rows of group_stats.csv.

    A group is tested when at least 3 young and 5 adult brains have a value. The
    Mann-Whitney BH q is taken within each reading and grouping.
    """
    rows = []
    for key in groups:
        for reading, _ in READINGS:
            # the values of each group and subgroup
            v = {m: value(reading, m, key, per, norm, refs) for m in mice}
            v = {m: x for m, x in v.items() if x is not None}
            yo = [v[m] for m in GROUPS["young"] if m in v]
            y20 = [v[m] for m in YOUNG_P20 if m in v]
            ad = [v[m] for m in ADULTS if m in v]
            nv = [v[m] for m in NAIVE if m in v]
            rw = [v[m] for m in RWS if m in v]
            if len(yo) < 3 or len(ad) < 5:
                continue
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
                    diff_P20only=(np.median(y20) - np.median(ad))
                    if len(y20) >= 3
                    else float("nan"),
                    naive_minus_rws=(np.median(nv) - np.median(rw))
                    if nv and rw
                    else float("nan"),
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


def write_group_stats(rows):
    """Write group_stats.csv into OUT, floats to four decimals."""
    with open(
        os.path.join(OUT, "group_stats.csv"), "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(
                {k: (f"{x:.4f}" if isinstance(x, float) else x) for k, x in r.items()}
            )


def print_group_table(groups, rows, star):
    """Print the median log2 young - adult of every group and reading, with its stars."""
    print(
        f"\n{'group':22s} "
        + " ".join(f"{r:>12s}" for r, _ in READINGS)
        + "   (log2 young - adult, medians)"
    )
    for grouping in ("system", "layer"):
        print(f"--- by {grouping}")
        for key in [k for k in groups if k[0] == grouping]:
            cells = []
            for reading, _ in READINGS:
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
                cells.append(
                    f"{r['diff_median']:+8.2f}{star[(reading, grouping, key[1])]:<3s}"
                    if r
                    else f"{'--':>11s}"
                )
            print(f"  {key[1]:22s} " + " ".join(cells))


def dotplot(keys, labels, fname, title, star, per, norm, refs):
    """Dot plot of the groups `keys` into OUT/`fname`, a panel per reading."""
    fig, axes = plt.subplots(
        len(READINGS),
        1,
        figsize=(max(9, 0.62 * len(keys) + 4), 3.6 * len(READINGS)),
        sharex=True,
    )
    for ax, (reading, rtitle) in zip(axes, READINGS):
        # each group's mice side by side, and its median as a bar
        for g in ("naive", "rws", "young"):
            ms = GROUPS[g]
            jit = np.linspace(-0.2, 0.2, len(ms))
            xo = 0.26 if g == "young" else -0.1
            mids = []
            for i, key in enumerate(keys):
                ys = [value(reading, m, key, per, norm, refs) for m in ms]
                ys = [y for y in ys if y is not None]
                for j, y in enumerate(ys):
                    ax.plot(
                        i + jit[j] + xo,
                        y,
                        "o",
                        ms=4.5,
                        color=COL[g],
                        alpha=0.9,
                        mec="none",
                    )
                mids.append(np.median(ys) if ys else np.nan)
            ax.plot(
                np.arange(len(keys)) + xo,
                mids,
                "_",
                ms=16,
                mew=2.4,
                color=COL[g],
                label=LABEL[g],
            )

        # stars for the young-vs-adult rank-sum test, just under the top of the panel
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.14 * (hi - lo))
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
                    color=COL["young"],
                )
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(rtitle, fontsize=10.5, loc="left")
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
    save_figure(fig, os.path.join(OUT, fname))
    plt.close(fig)


def plot_groups(groups, star, per, norm, refs):
    """Draw group_plot.png (the systems) and laminar_plot.png (the layers)."""
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


def main():
    """Group means per mouse, young-against-adult tests, group_stats.csv, both plots."""
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
    star = {
        (r["reading"], r["grouping"], r["group"]): (
            "**"
            if r["mannwhitney_p"] < 0.01
            else ("*" if r["mannwhitney_p"] < 0.05 else "")
        )
        for r in rows
    }

    # summary table, printed
    print_group_table(groups, rows, star)

    # figures
    plot_groups(groups, star, per, norm, refs)
    print("\nwrote group_plot.png, laminar_plot.png, group_stats.csv")
