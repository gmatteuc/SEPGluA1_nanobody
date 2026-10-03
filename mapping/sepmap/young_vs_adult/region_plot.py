"""Per-mouse region statistics and dot plot, each brain measured on its own atlas.

This is where the numbers to quote come from. A region-wise comparison needs only
every voxel's label, and the DeMBA annotations were remapped to the Allen
parcellation_index when they were built, so VISp in a P16 brain is measured on
the P16 annotation, VISp in a P20 brain on the P20 one, and only the per-mouse
means are compared. Pooling ages is therefore clean here, as it is not for a
voxelwise map, where the young brains must first be carried into the CCF
(volumes.to_ccf). Nothing is warped.

Per mouse and structure, the mean of sig, of sig/auto and of sig/SEP over the
tissue voxels, and from them per mouse:

    ratio      nano per unit autofluorescence, no reference. Not an absolute
               measure: the young cortex is 2.0 log2 below the adult in nano and
               1.0 log2 below it in autofluorescence, so the denominator carries
               an age effect of its own (see volumes.cohort).
    sepratio   nano per unit SEP, meant as membrane receptor per unit receptor
               expressed, SEP being the tag on GluA1 itself. It is not that:
               adult.sep_channel_check finds the green channel dominated by
               autofluorescence in this tissue, so the reading behaves as a
               second nano over autofluorescence.
    cref       sig over the brain's isocortex mean: a share of the cortex.
    subref     sig over the mean of the subcortex without HPF and STR, the two
               structures that would dominate the scale.
    zref       range-matched: the cortex-relative values minus the brain's
               median over structures, divided by its own p90-p10 spread.

zref gives every brain the same level and the same dynamic range, so the question
becomes where a region sits within its own brain's range. The young brain is
flatter, p90-p10 = 0.89 +- 0.25 log2 against 1.79 +- 0.24 in adults, and no
one-number reference can touch that: dividing by cortex, subcortex or hippocampus
shifts every point equally and only moves the zero. The price is that the
compression is defined away, so zref can show a reordering but says nothing about
amplitude.

Per structure, the young group is tested against the adults (Welch), with the
naive-minus-rws difference beside it as the size of a difference that carries no
developmental meaning. The young group pools every registered young brain
whatever its age, each one mouse like any other in the median and the tests; the
P20-only contrast stays in the CSV, so what the other ages do to the answer can be
checked. The figure marks each region with the Mann-Whitney (rank-sum) test,
uncorrected: ranks rather than means, because the samples are small and log
ratios need not be normal. region_stats.csv carries the Welch p, the Mann-Whitney
p and the Benjamini-Hochberg q of each, so a corrected reading is one column away.

Writes region_means_per_mouse.csv, region_stats.csv and region_plot.png into
comparisons_v2/young_vs_adult/.

Run by run_region_plot.py.
"""

import csv
import math
from collections import defaultdict
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np

from sepmap.config import SETTINGS
from sepmap.plotting import GROUP_COLOURS, save_figure
from sepmap.volumes.cohort import (
    MODES,
    NAIVE,
    NOT_SUBCORTEX,
    RWS,
    SIGNED_READINGS,
    YOUNG_P16,
    YOUNG_P20,
    YOUNG_P22,
    per_unit,
)
from sepmap.volumes.per_mouse import DATA, MICE, annotation_20, structure_terms
from sepmap.volumes.per_mouse import OUT as PER_MOUSE

OUT = DATA / "comparisons_v2" / "young_vs_adult"

# the smallest structure kept, and the brains a structure needs to be tested
REGION_TABLES = SETTINGS["region_tables"]
REGION_PLOT = SETTINGS["region_plot"]

# areas of the printed table and the figure: cortex, then after "|" the subcortex
AREAS = [
    "VISp",
    "VISl",
    "VISal",
    "VISrl",
    "VISpm",
    "VISam",
    "SSp-bfd",
    "SSp-ul",
    "SSp-ll",
    "SSp-m",
    "SSp-n",
    "SSs",
    "AUDp",
    "AUDd",
    "MOp",
    "MOs",
    "RSPd",
    "RSPv",
    "ACAd",
    "ACAv",
    "PL",
    "ILA",
    "ORBl",
    "|",
    "VPM",
    "VPL",
    "LGd",
    "LP",
    "CP",
    "ACB",
    "CA1",
    "CA3",
    "DG",
    "GPe",
    "PVH",
    "ZI",
]

# the readings and their panel titles
READINGS = [
    ("ratio", "nanobody / autofluorescence, both background-subtracted  (log2)"),
    (
        "sepratio",
        "nanobody / SEP  -  NOT a surface fraction: SEP is mostly autofluorescence "
        "here  (log2)",
    ),
    (
        "cref",
        "background-subtracted nanobody, relative to the mouse's own isocortex  (log2)",
    ),
    (
        "subref",
        "background-subtracted nanobody, relative to subcortex excluding HPF and STR  "
        "(log2)",
    ),
    (
        "zref",
        "range-matched: cortex-relative, then centred and scaled by each brain's own "
        "spread",
    ),
]

# only the readings in force (volumes.cohort: settings.toml, or V2_READINGS for one
# run), so one setting covers the chain; young_vs_adult.region_groups follows this list
READINGS = [r for r in READINGS if r[0] in MODES]

# the groups: every young brain (P16, P20, P22) against the adults, naive and rws
GROUPS = {"young": YOUNG_P20 + YOUNG_P16 + YOUNG_P22, "naive": NAIVE, "rws": RWS}
ADULTS = NAIVE + RWS

# legend labels, with the counts computed from the lists
LABEL = {
    "young": f"young P16-P22 (n = {len(YOUNG_P20) + len(YOUNG_P16) + len(YOUNG_P22)})",
    "naive": f"adult naive (n = {len(NAIVE)})",
    "rws": f"adult rws (n = {len(RWS)})",
}


# the columns of region_means_per_mouse.csv and of region_stats.csv
PER_MOUSE_COLUMNS = [
    "reading",
    "group",
    "cohort",
    "mouse",
    "structure",
    "acronym",
    "division",
    "n_vox20",
    "log2_value",
]
STATS_COLUMNS = [
    "reading",
    "structure",
    "acronym",
    "division",
    "n_young",
    "n_adult",
    "young_mean_log2",
    "adult_mean_log2",
    "diff_log2",
    "diff_median_log2",
    "welch_p",
    "mannwhitney_p",
    "welch_q_BH",
    "mannwhitney_q_BH",
    "diff_log2_P20only",
    "welch_p_P20only",
    "diff_log2_P16_single",
    "naive_minus_rws_log2",
]


def stars(p: float) -> str:
    """The stars of a p-value: ** under 0.01, * under 0.05, otherwise none."""
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def bh_fdr(p: list[float] | np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values for one family of tests; NaN where p is NaN.

    A few hundred structures are tested per reading, so a handful of p < 0.05
    is expected from noise alone. The q-value is what should be quoted for
    anything other than the regions named in advance (SS, VIS, prefrontal).
    """
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if ok.sum() == 0:
        return q

    # p times n over its rank, made monotone from the largest p down, at most 1
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    n = len(ranked)
    adj = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def mannwhitney(a: list[float], b: list[float]) -> float:
    """Two-sided rank-sum p for two independent samples, or NaN if too small."""
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    from scipy.stats import mannwhitneyu

    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def welch(a: list[float], b: list[float]) -> float:
    """Two-sided Welch t-test p for two independent samples.

    NaN when a sample has fewer than two values or neither varies.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    if va + vb == 0:
        return float("nan")

    # t with unequal variances, and the Welch-Satterthwaite degrees of freedom
    t = (a.mean() - b.mean()) / math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va**2 / (len(a) - 1) + vb**2 / (len(b) - 1))
    from scipy.stats import t as tdist

    return float(2 * tdist.sf(abs(t), df))


def structure_means(mice: list[str], names: dict[int, str]) -> dict[str, dict]:
    """Per mouse and structure: tissue voxels, and the mean of sig, ratio and sepratio.

    Each brain on its own atlas; structures under region_tables.min_vox20 voxels
    are left out. The
    layers of an area are separate parcellation indices with one structure name,
    and are pooled under that name. Stops if a brain has no SEP channel while the
    sepratio reading is in force.
    """
    anns, per = {}, {}
    for mouse in mice:
        cohort, atlas_key = MICE[mouse][:2]
        if atlas_key not in anns:
            anns[atlas_key] = annotation_20(atlas_key)
        ann = anns[atlas_key]
        z = np.load(PER_MOUSE / (mouse + ".npz"))
        sig = z["sig"].astype(np.float32)
        auto = z["auto"].astype(np.float32)
        tissue = z["tissue"]
        if "sepratio" in MODES and "sep" not in z.files:
            raise ValueError(
                f"{mouse}: no SEP channel in its per-mouse file. Run\n"
                "  run_add_sep_channel.m for this brain, then run_per_mouse.py,\n"
                "  or drop the reading with V2_READINGS."
            )

        # sums per parcellation index of the voxels, sig and the two ratios
        ratio = per_unit(sig, auto, tissue)
        if "sep" in z.files:
            sepratio = per_unit(sig, z["sep"].astype(np.float32), tissue)
        else:
            sepratio = np.zeros_like(sig)
        lab = ann[tissue]
        nlab = int(ann.max()) + 1
        n = np.bincount(lab, minlength=nlab)
        s_sig = np.bincount(lab, weights=sig[tissue], minlength=nlab)
        s_rat = np.bincount(lab, weights=ratio[tissue], minlength=nlab)
        s_sep = np.bincount(lab, weights=sepratio[tissue], minlength=nlab)

        # pooled per structure name, then the means of the structures large enough
        d = defaultdict(lambda: [0, 0.0, 0.0, 0.0])
        for idx in np.nonzero(n)[0]:
            if idx == 0:
                continue
            key = names.get(int(idx), f"id{idx}")
            d[key][0] += int(n[idx])
            d[key][1] += s_sig[idx]
            d[key][2] += s_rat[idx]
            d[key][3] += s_sep[idx]
        per[mouse] = {
            k: (v[0], v[1] / v[0], v[2] / v[0], v[3] / v[0])
            for k, v in d.items()
            if v[0] >= REGION_TABLES["min_vox20"]
        }
        print(f"{mouse:20s} {atlas_key:10s} {len(per[mouse])} structures", flush=True)
    return per


def ref(
    m: str,
    pred: Callable[[str], bool],
    per: dict[str, dict],
    meta: dict[str, tuple[str, str]],
) -> float:
    """Voxel-weighted mean sig of mouse `m` over structures whose division passes `pred`.

    NaN when no structure does.
    """
    s = c = 0.0
    for k, (n, ms, _, _) in per[m].items():
        if pred(meta.get(k, ("", ""))[1]):
            s += ms * n
            c += n
    return s / c if c else float("nan")


def brain_references(
    mice: list[str], per: dict[str, dict], meta: dict[str, tuple[str, str]]
) -> dict[str, dict[str, float]]:
    """Per mouse the two one-number references: mean sig of isocortex and subcortex."""
    # one number per brain, so every reading is a pure scale and region ratios within
    # a brain survive it exactly; only the question asked changes (module docstring)
    refs = {
        m: {
            "cref": ref(m, lambda d: d == "Isocortex", per, meta),
            "subref": ref(m, lambda d: d not in NOT_SUBCORTEX, per, meta),
        }
        for m in mice
    }
    return refs


def range_match(
    mice: list[str], per: dict[str, dict], refs: dict[str, dict[str, float]]
) -> dict[str, tuple[float, float]]:
    """Per mouse the median and the p90-p10 spread of log2 cortex-relative sig.

    Taken over the structures every brain has, and printed.
    """
    # the same structures in every brain, or the spread would depend on which
    # regions the sections happened to cover
    common = set.intersection(*[set(per[m]) for m in mice])
    norm = {}
    for m in mice:
        v = np.array(
            [
                math.log2(per[m][k][1] / refs[m]["cref"])
                for k in sorted(common)
                if per[m][k][1] > 0
            ]
        )
        p10, med, p90 = np.percentile(v, [10, 50, 90])
        norm[m] = (med, max(p90 - p10, 1e-6))
    print(
        "dynamic range per brain (p90-p10 of log2 over %d shared structures):"
        % len(common)
    )
    for m in mice:
        print(f"  {m:20s} median {norm[m][0]:+.2f}   spread {norm[m][1]:.2f}")
    return norm


def value(
    reading: str,
    m: str,
    k: str | None,
    per: dict[str, dict],
    norm: dict[str, tuple[float, float]],
    refs: dict[str, dict[str, float]],
) -> float | None:
    """Value of `reading` for mouse `m` in structure `k`, or None.

    log2 for every reading but the signed ones, which are range-matched. None when
    the structure is missing in that brain or its value is not positive.
    """
    if k is None or k not in per[m]:
        return None
    n, ms, mr, msep = per[m][k]
    if reading in SIGNED_READINGS:
        if ms <= 0:
            return None
        med, spread = norm[m]
        return (math.log2(ms / refs[m]["cref"]) - med) / spread

    # ratio and sepratio are already ratios, taken voxel by voxel against a
    # channel; the rest divide sig by a single number measured on this brain
    if reading in ("ratio", "sepratio"):
        v = mr if reading == "ratio" else msep
    else:
        v = ms / refs[m][reading]
    return math.log2(v) if v > 0 else None


def structure_row(
    reading: str,
    k: str,
    meta: dict[str, tuple[str, str]],
    yo: list[float],
    y20: list[float],
    p16: list[float],
    ad: list[float],
    nv: list[float],
    rw: list[float],
) -> dict:
    """One structure's test row, in the columns of region_stats.csv but the q-values.

    `yo`, `y20`, `p16`, `ad`, `nv` and `rw` are the values of the young, the P20
    and P16 young, the adults, the naive and the RWS adults.
    """
    return dict(
        reading=reading,
        structure=k,
        acronym=meta[k][0],
        division=meta[k][1],
        n_young=len(yo),
        n_adult=len(ad),
        young_mean_log2=np.mean(yo),
        adult_mean_log2=np.mean(ad),
        diff_log2=np.mean(yo) - np.mean(ad),
        diff_median_log2=np.median(yo) - np.median(ad),
        welch_p=welch(yo, ad),
        mannwhitney_p=mannwhitney(yo, ad),
        diff_log2_P20only=(np.mean(y20) - np.mean(ad)) if len(y20) >= 2 else float("nan"),
        welch_p_P20only=welch(y20, ad) if len(y20) >= 2 else float("nan"),
        diff_log2_P16_single=(p16[0] - np.mean(ad)) if p16 else float("nan"),
        naive_minus_rws_log2=(np.mean(nv) - np.mean(rw)) if nv and rw else float("nan"),
    )


def with_q_values(rows_st: list[dict]) -> list[dict]:
    """The test rows with the BH q of the Welch and of the Mann-Whitney p added.

    The q-values are taken within each reading, so the brain-wide lists can be
    read honestly; in the CSV they go after the Mann-Whitney p.
    """
    q_welch, q_mw = {}, {}
    for reading, _ in READINGS:
        idx = [i for i, r in enumerate(rows_st) if r["reading"] == reading]
        for store, col in ((q_welch, "welch_p"), (q_mw, "mannwhitney_p")):
            for i, qi in zip(idx, bh_fdr([rows_st[i][col] for i in idx])):
                store[i] = qi
    return [
        {**r, "welch_q_BH": q_welch[i], "mannwhitney_q_BH": q_mw[i]}
        for i, r in enumerate(rows_st)
    ]


def region_rows(
    mice: list[str],
    meta: dict[str, tuple[str, str]],
    group_of: dict[str, str],
    per: dict[str, dict],
    norm: dict[str, tuple[float, float]],
    refs: dict[str, dict[str, float]],
) -> tuple[list[dict], list[dict]]:
    """The per-mouse rows, and per structure and reading the young-against-adult tests.

    A structure is tested when at least region_plot.min_young young and min_adult
    adult brains have a value.
    The test rows carry the BH q of the Welch and of the Mann-Whitney p, taken
    within each reading.
    """
    rows_pm, rows_st = [], []
    for k in sorted(meta, key=lambda k: (meta[k][1], meta[k][0])):
        for reading, _ in READINGS:
            v = {m: value(reading, m, k, per, norm, refs) for m in mice}
            v = {m: x for m, x in v.items() if x is not None}

            # one row per mouse, in the columns of region_means_per_mouse.csv
            for m, x in v.items():
                rows_pm.append(
                    dict(
                        reading=reading,
                        group=group_of[m],
                        cohort=MICE[m][0],
                        mouse=m,
                        structure=k,
                        acronym=meta[k][0],
                        division=meta[k][1],
                        n_vox20=per[m][k][0],
                        log2_value=x,
                    )
                )

            # the values of each group and subgroup
            yo = [v[m] for m in GROUPS["young"] if m in v]
            y20 = [v[m] for m in YOUNG_P20 if m in v]
            p16 = [v[m] for m in YOUNG_P16 if m in v]
            ad = [v[m] for m in ADULTS if m in v]
            nv = [v[m] for m in NAIVE if m in v]
            rw = [v[m] for m in RWS if m in v]
            if len(yo) < REGION_PLOT["min_young"] or len(ad) < REGION_PLOT["min_adult"]:
                continue

            # one row per structure, in the columns of region_stats.csv but the
            # q-values, which go in after the Mann-Whitney p below
            rows_st.append(structure_row(reading, k, meta, yo, y20, p16, ad, nv, rw))

    # q-values within each reading, so the brain-wide lists can be read honestly
    rows_st = with_q_values(rows_st)
    return rows_pm, rows_st


def write_tables(rows_pm: list[dict], rows_st: list[dict]) -> None:
    """Write region_means_per_mouse.csv and region_stats.csv into OUT.

    Values to four decimals, but the text columns and the two counts.
    """
    with open(
        OUT / "region_means_per_mouse.csv", "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.DictWriter(fh, fieldnames=PER_MOUSE_COLUMNS)
        w.writeheader()
        w.writerows([{**r, "log2_value": f"{r['log2_value']:.4f}"} for r in rows_pm])
    as_is = STATS_COLUMNS[:6]
    with open(OUT / "region_stats.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=STATS_COLUMNS)
        w.writeheader()
        w.writerows(
            [
                {c: r[c] if c in as_is else f"{r[c]:.4f}" for c in STATS_COLUMNS}
                for r in rows_st
            ]
        )


def print_cortex_table(rows_st: list[dict]) -> None:
    """Print log2(young / adult) of each area in AREAS and reading, with its stars."""
    st = {(r["reading"], r["acronym"]): r for r in rows_st}
    print(
        f"\nCORTEX  log2(young / adult), young = {len(GROUPS['young'])} mice (P20 + P16) "
        f"vs {len(ADULTS)} adults "
        "(* p<0.05, ** p<0.01, Mann-Whitney, uncorrected; q in the CSV). "
        "P20only = without the P16 brain; P16 = that brain alone; "
        "naive-rws = the null scale."
    )
    print(
        f"  {'area':9s} "
        + " ".join(f"{r:>10s}" for r, _ in READINGS)
        + f" {'P20only':>9s} {'P16':>7s} {'naive-rws':>10s}"
    )
    for a in AREAS:
        if a == "|":
            print("  " + "-" * 60)
            continue
        cells = []
        for reading, _ in READINGS:
            r = st.get((reading, a))
            if r is None:
                cells.append(f"{'--':>10s}")
                continue

            # stars from the rank-sum p
            star = stars(r["mannwhitney_p"])
            cells.append(f"{r['diff_log2']:+7.2f}{star:3s}")
        r = st.get(("cref", a))
        tail = ""
        if r:
            tail = (
                f"{r['diff_log2_P20only']:+9.2f} {r['diff_log2_P16_single']:+7.2f} "
                f"{r['naive_minus_rws_log2']:+10.2f}"
            )
        print(f"  {a:9s} " + " ".join(cells) + " " + tail)


def draw_group_dots(
    ax: plt.Axes,
    reading: str,
    g: str,
    xs: list[int],
    by_acro: dict[str, str],
    per: dict[str, dict],
    norm: dict[str, tuple[float, float]],
    refs: dict[str, dict[str, float]],
) -> None:
    """Draw group `g`'s mice side by side in each of AREAS, and its median as a bar."""
    ms = GROUPS[g]
    jit = np.linspace(-0.22, 0.22, len(ms))
    xo = 0.28 if g == "young" else -0.1
    mids = []
    for i, a in enumerate(AREAS):
        if a == "|":
            mids.append(np.nan)
            continue
        k = by_acro.get(a)
        ys = []
        for j, m in enumerate(ms):
            y = value(reading, m, k, per, norm, refs)
            if y is None:
                continue
            ys.append(y)
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
    ax.plot(
        np.array(xs) + xo,
        [mids[i] for i in xs],
        "_",
        ms=14,
        mew=2.2,
        color=GROUP_COLOURS[g],
        label=LABEL[g],
    )


def draw_stars(
    ax: plt.Axes, reading: str, star_of: dict[tuple[str, str], str]
) -> tuple[float, float]:
    """Draw the rank-sum stars just under the top of the panel; returns its y limits."""
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.12 * (hi - lo))
    lo, hi = ax.get_ylim()
    for i, a in enumerate(AREAS):
        st = star_of.get((reading, a), "")
        if st:
            ax.text(
                i,
                hi - 0.04 * (hi - lo),
                st,
                ha="center",
                va="top",
                fontsize=11,
                color=GROUP_COLOURS["young"],
            )
    return lo, hi


def draw_divide(ax: plt.Axes, lo: float, hi: float) -> None:
    """Draw the zero line and the cortex-subcortex divide, labelled."""
    ax.axhline(0, color="k", lw=0.6)
    sep = AREAS.index("|")
    ax.axvline(sep, color="k", lw=0.6, ls=":")
    box = dict(facecolor="w", edgecolor="none", alpha=0.85, pad=1.5)
    ax.text(
        sep - 0.5,
        lo + 0.02 * (hi - lo),
        "cortex",
        ha="right",
        va="bottom",
        fontsize=9,
        color="#333",
        bbox=box,
    )
    ax.text(
        sep + 0.5,
        lo + 0.02 * (hi - lo),
        "subcortex",
        ha="left",
        va="bottom",
        fontsize=9,
        color="#333",
        bbox=box,
    )


def plot_regions(
    rows_st: list[dict],
    by_acro: dict[str, str],
    per: dict[str, dict],
    norm: dict[str, tuple[float, float]],
    refs: dict[str, dict[str, float]],
) -> None:
    """Draw region_plot.png: a panel per reading, a dot per mouse in each of AREAS."""
    # every mouse a dot of the same size, the P16 brain included; the bar is the
    # group median, to match the rank-sum test that sets the stars
    star_of = {(r["reading"], r["acronym"]): stars(r["mannwhitney_p"]) for r in rows_st}
    ylab = {
        "ratio": "log2  nano / auto",
        "sepratio": "log2  nano / SEP",
        "cref": "log2  relative to own isocortex",
        "subref": "log2  relative to subcortex",
        "zref": "range-matched (median 0, spread 1)",
    }
    fig, axes = plt.subplots(
        len(READINGS), 1, figsize=(15, 3.8 * len(READINGS)), sharex=True
    )
    xs = [i for i, a in enumerate(AREAS) if a != "|"]
    for ax, (reading, title) in zip(axes, READINGS):
        # each group's mice, side by side, and its median as a bar
        for g in ("naive", "rws", "young"):
            draw_group_dots(ax, reading, g, xs, by_acro, per, norm, refs)

        # stars for the young-vs-adult rank-sum test, just under the top of the
        # panel; the zero line, and the cortex-subcortex divide
        lo, hi = draw_stars(ax, reading, star_of)
        draw_divide(ax, lo, hi)
        ax.set_title(title, fontsize=10.5, loc="left")
        ax.set_ylabel(ylab[reading], fontsize=10)
        ax.grid(axis="y", lw=0.3, alpha=0.6)
        ax.set_xlim(-0.8, len(AREAS) - 0.2)
    axes[0].legend(
        loc="lower left",
        fontsize=9,
        frameon=True,
        framealpha=0.9,
        edgecolor="none",
        ncol=3,
    )
    axes[-1].set_xticks(range(len(AREAS)))
    axes[-1].set_xticklabels(
        ["" if a == "|" else a for a in AREAS], rotation=60, ha="right", fontsize=9
    )
    fig.suptitle(
        "Young vs adult, nano channel: one dot per mouse, "
        "each brain measured on the atlas of its own age\n"
        f"Bars are group medians.  * p<0.05, ** p<0.01, Mann-Whitney "
        f"{len(GROUPS['young'])} vs {len(ADULTS)}, uncorrected",
        fontsize=11.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    save_figure(fig, OUT / "region_plot.png", dpi=105)


def main() -> None:
    """Measure the region means per mouse, test young against adult, write and draw."""
    # structure names, acronyms and divisions of the ontology
    names, acro, divi = structure_terms()

    # per mouse: structure means, each brain on its own atlas
    mice = [m for g in GROUPS.values() for m in g]
    per = structure_means(mice, names)
    meta = {nm: (acro[idx], divi.get(idx, "")) for idx, nm in names.items()}
    by_acro = {meta[k][0]: k for k in meta}
    group_of = {m: g for g, ms in GROUPS.items() for m in ms}

    # the references of each brain, and its range match
    refs = brain_references(mice, per, meta)
    norm = range_match(mice, per, refs)

    # per structure: the values of each mouse and the tests, as two tables
    rows_pm, rows_st = region_rows(mice, meta, group_of, per, norm, refs)
    write_tables(rows_pm, rows_st)

    # cortex summary, printed
    print_cortex_table(rows_st)

    # the dot plot
    plot_regions(rows_st, by_acro, per, norm, refs)
    print("\nwrote", OUT / "region_plot.png")
