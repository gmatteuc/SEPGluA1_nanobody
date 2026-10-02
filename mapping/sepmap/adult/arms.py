"""The three channel arms of the measurement argument, one region table per adult.

The gene ranking alone cannot say whether the nanobody reports receptor on the
membrane or receptor anywhere: Gria1, Cacng8, Dlg2 and Grip1 are all postsynaptic,
so no grouping of genes (ours, GO's or SynGO's) separates the two. A contrast
between channels could, and the SEP channel that run_add_sep_channel carries into
registered space gives three arms:

    sepauto    SEP / autofluorescence     total receptor, however it is localised
    ratio      nano / autofluorescence    receptor at the surface
    sepratio   nano / SEP                 the surface fraction

In log space

    sepratio = ratio - sepauto

up to the difference between a mean of ratios and a ratio of means, so
autofluorescence, and anything else common to both channels, cancels out of the
third. The predictions, stated before the numbers existed: sepauto tracks Gria1
expression best of the three, ratio sits in between, and sepratio tracks the
trafficking and anchoring machinery best and Gria1 least.

The premise does not hold in this tissue. SEP is superecliptic pHluorin on GluA1,
so in a living cell it reports the surface pool; fixed, cleared and mounted tissue
has lost its pH gradients, so the green channel was expected to report the
receptor wherever it sits. It does not: adult.sep_channel_check finds that it
tracks the autofluorescence channel at rho 0.79 +- 0.04 across the ten adults,
with a dynamic range of 0.95 log2 against nano's 1.93; whatever tag survives the
protocol, autofluorescence dominates what is left. The table is still correct
arithmetic, but sepauto is not total receptor and sepratio is not a surface
fraction. They are kept because they are how that was established, and because
ratio and the consistency check are needed either way; read
adult.sep_channel_check before using either.

The arithmetic is the same as young_vs_adult.region_plot's: the same 20 um
annotation, the same smallest structure, the same mask-normalised smoothing of the
denominator, structures keyed by name over their layer indices. The module checks
itself against that table for the two arms both compute, and prints the
difference, since a silent divergence would invalidate every comparison
downstream.

Writes, in adult_v2/arms/ under the data root:

    region_means_arms.csv    arm x mouse x structure, log2
    arms_consistency.png     the self-check: against young_vs_adult.region_plot,
                             and the log-space identity between the three arms

Run by run_adult_arms.py.
"""

import csv
import math
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from sepmap.config import SETTINGS
from sepmap.plotting import RED, tidy
from sepmap.volumes.cohort import NAIVE, RWS, per_unit
from sepmap.volumes.per_mouse import DATA, MICE, annotation_20, structure_terms
from sepmap.volumes.per_mouse import OUT as PER_MOUSE

# the smallest structure kept, as in young_vs_adult.region_plot
REGION_TABLES = SETTINGS["region_tables"]

OUT = os.path.join(DATA, "adult_v2", "arms")

# young_vs_adult.region_plot's per-mouse table, which the shared arms must match
EXISTING = os.path.join(
    DATA, "comparisons_v2", "young_vs_adult", "region_means_per_mouse.csv"
)

# the ten adults, naive and rws pooled
ADULTS = NAIVE + RWS

# the three arms, and labels for them that no figure of this module uses
ARMS = ("sepauto", "ratio", "sepratio")
LABEL = {
    "sepauto": "SEP / autofluorescence   (total receptor)",
    "ratio": "nano / autofluorescence   (surface receptor)",
    "sepratio": "nano / SEP   (surface fraction)",
}

# the two arms young_vs_adult.region_plot also computes, and must agree on
SHARED = {"ratio": "ratio", "sepratio": "sepratio"}


def mouse_table(
    mouse: str, names: dict[int, str]
) -> dict[str, tuple[int, dict[str, float]]]:
    """The mean of each arm per structure for one adult: {name: (n voxels, {arm: mean})}.

    The layer indices of a structure are pooled by name; a structure under
    region_tables.min_vox20 voxels is left out.
    """
    z = np.load(os.path.join(PER_MOUSE, mouse + ".npz"))
    if "sep" not in z.files:
        raise ValueError(
            f"{mouse}: no SEP channel. Run run_add_sep_channel.m, "
            f"then run_per_mouse.py, for this brain."
        )
    sig = z["sig"].astype(np.float32)
    auto = z["auto"].astype(np.float32)
    sep = z["sep"].astype(np.float32)
    tissue = z["tissue"]
    annotation = annotation_20(MICE[mouse][1])

    vol = {
        "sepauto": per_unit(sep, auto, tissue),
        "ratio": per_unit(sig, auto, tissue),
        "sepratio": per_unit(sig, sep, tissue),
    }

    # voxel count and sum of each arm per annotation index
    labels = annotation[tissue]
    n_labels = int(annotation.max()) + 1
    n = np.bincount(labels, minlength=n_labels)
    sums = {
        arm: np.bincount(labels, weights=vol[arm][tissue], minlength=n_labels)
        for arm in ARMS
    }

    # pool the indices by structure name (index 0 is outside the brain)
    acc = defaultdict(lambda: [0] + [0.0] * len(ARMS))
    for idx in np.nonzero(n)[0]:
        if idx == 0:
            continue
        a = acc[names.get(int(idx), f"id{idx}")]
        a[0] += int(n[idx])
        for j, arm in enumerate(ARMS, 1):
            a[j] += sums[arm][idx]
    return {
        k: (v[0], {arm: v[j] / v[0] for j, arm in enumerate(ARMS, 1)})
        for k, v in acc.items()
        if v[0] >= REGION_TABLES["min_vox20"]
    }


def check_against_existing(rows: list[dict]) -> dict[str, list[float]]:
    """Check that the two arms young_vs_adult.region_plot also computes are identical.

    They are computed here from the same per-mouse files with the same
    arithmetic, so anything above rounding means the two modules have drifted
    apart and nothing downstream can be trusted; a RuntimeError stops the run
    then. Returns the absolute differences per arm, or an empty dict when that
    table does not exist.
    """
    if not os.path.exists(EXISTING):
        print("no existing table to check against -- skipped")
        return {}
    ours = {(r["arm"], r["mouse"], r["structure"]): float(r["log2_value"]) for r in rows}
    diffs = defaultdict(list)
    with open(EXISTING, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            key = (SHARED.get(r["reading"], ""), r["mouse"], r["structure"])
            if key in ours:
                diffs[key[0]].append(abs(ours[key] - float(r["log2_value"])))

    # that table is written to four decimals, so half of the last digit is the
    # most the two can differ by if they are computing the same thing
    tol = 0.5e-4
    print(
        f"\nagainst run_region_plot (bound is {tol:.1e}, half the last digit it stores):"
    )
    ok = True
    for arm, d in sorted(diffs.items()):
        worst = max(d)
        ok &= worst <= tol
        print(
            f"  {arm:9s} n = {len(d):5d}   max |difference| = {worst:.3e}   "
            f"{'agrees' if worst <= tol else 'DRIFTED'}"
        )
    if not ok:
        raise RuntimeError(
            "the two scripts no longer compute the same thing -- stop here"
        )
    return diffs


def figure(
    per: dict[str, dict[str, tuple[int, dict[str, float]]]], diffs: dict[str, list[float]]
) -> None:
    """Draw the self-check: the shared arms' agreement and the log-space identity."""
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9))

    # the differences of the shared arms, jittered, on a log axis (an exact zero
    # is drawn at 1e-17); each arm's jitter restarts from seed 0
    ax = axes[0]
    if diffs:
        for i, (arm, d) in enumerate(sorted(diffs.items())):
            ax.scatter(
                np.full(len(d), i)
                + np.random.default_rng(0).uniform(-0.12, 0.12, len(d)),
                np.maximum(d, 1e-17),
                s=5,
                facecolor="0.4",
                edgecolor="none",
            )
        ax.set_yscale("log")
        ax.set_xticks(range(len(diffs)))
        ax.set_xticklabels(sorted(diffs), fontsize=8)
        ax.axhline(1e-9, color=RED, lw=0.8, ls="--")
        ax.set_ylabel("|this script - run_region_plot|  (log2 units)", fontsize=8)
    ax.set_title("the two shared arms agree", fontsize=9)

    # log(nano/sep) should equal log(nano/auto) - log(sep/auto), except that
    # each arm is a mean of voxelwise ratios rather than a ratio of means
    gap = []
    for mouse, table in per.items():
        for k, (_, m) in table.items():
            if min(m.values()) > 0:
                gap.append(
                    math.log2(m["sepratio"])
                    - (math.log2(m["ratio"]) - math.log2(m["sepauto"]))
                )
    axes[1].hist(gap, bins=60, color="0.6", edgecolor="0.3", linewidth=0.4)
    axes[1].axvline(0, color=RED, lw=0.9)
    axes[1].set_xlabel(
        "log2(nano/SEP)  -  [log2(nano/auto) - log2(SEP/auto)]", fontsize=8
    )
    axes[1].set_title(
        f"Jensen gap: median {np.median(gap):+.3f}, "
        f"p5-p95 {np.percentile(gap, 5):+.2f} "
        f"to {np.percentile(gap, 95):+.2f}",
        fontsize=9,
    )
    axes[1].set_ylabel("structures x mice", fontsize=8)

    # the two channels against each other in the first mouse, one dot per structure
    ax = axes[2]
    mouse = sorted(per)[0]
    t = per[mouse]
    ks = [k for k in t if min(t[k][1].values()) > 0]
    ax.scatter(
        [math.log2(t[k][1]["sepauto"]) for k in ks],
        [math.log2(t[k][1]["ratio"]) for k in ks],
        s=9,
        facecolor="0.55",
        edgecolor="0.25",
        linewidth=0.3,
    )
    ax.set_xlabel("log2 SEP / auto", fontsize=8)
    ax.set_ylabel("log2 nano / auto", fontsize=8)
    ax.set_title(f"the two channels, {mouse}", fontsize=9)

    for ax in axes:
        tidy(ax)
    fig.tight_layout()
    path = os.path.join(OUT, "arms_consistency.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"{path}")


def main() -> None:
    """Write the arms table of the ten adults, check it, and draw the check."""
    os.makedirs(OUT, exist_ok=True)
    names, acro, divi = structure_terms()
    meta = {nm: (acro[idx], divi.get(idx, "")) for idx, nm in names.items()}
    group = {m: ("naive" if m in NAIVE else "rws") for m in ADULTS}

    # one row per arm, mouse and structure with a positive mean
    per, rows = {}, []
    for mouse in ADULTS:
        per[mouse] = mouse_table(mouse, names)
        print(f"{mouse:20s} {len(per[mouse])} structures", flush=True)
        for k, (n, m) in per[mouse].items():
            for arm in ARMS:
                if m[arm] > 0:
                    rows.append(
                        dict(
                            arm=arm,
                            group=group[mouse],
                            mouse=mouse,
                            structure=k,
                            acronym=meta.get(k, ("", ""))[0],
                            division=meta.get(k, ("", ""))[1],
                            n_vox20=n,
                            log2_value=f"{math.log2(m[arm]):.6f}",
                        )
                    )

    path = os.path.join(OUT, "region_means_arms.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n{len(rows):,} rows -> {path}")

    # the check stops the run before the figure when the shared arms have drifted
    diffs = check_against_existing(rows)
    figure(per, diffs)
