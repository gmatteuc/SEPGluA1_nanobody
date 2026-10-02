"""Whether the green channel reports tagged receptor or the tissue.

The three-arm argument of ish.arms rests on one premise, that ex vivo the SEP
channel reports total GluA1, and the arm built on it behaved strangely: SEP / auto
correlates with Gria1 expression at -0.11, while plain nano reaches +0.62. Either
the premise or the arm is wrong, and that is settled on the channels themselves
before any ratio is taken: one mean per structure per adult, in log2, straight from
the per-mouse files, with no denominator anywhere. Three things decide it.

    dynamic range   how much a channel varies across the brain. Total receptor is
                    not flat (hippocampus against thalamus is a large, well known
                    difference), so a channel reporting it should have a range
                    comparable to the nanobody's.
    who it tracks   if the green channel is tag, its spatial profile should look
                    like the nanobody's and like Gria1's; if it is tissue, like
                    the autofluorescence channel's.
    what is left    log2(SEP) regressed on log2(auto), and the residual asked
                    whether it still carries a receptor map: if the tag survives
                    under the autofluorescence, an unmixed SEP could still serve
                    as the total-receptor arm, and this says how much is there.

The answer is a property of this tissue, fixed, cleared and mounted, not of
SEP-GluA1 mice. Superecliptic pHluorin is quenched in acidic compartments in a
living cell, which makes it a surface reporter in vivo; after fixation and
clearing the pH gradients are gone, and whatever green emission survives the
protocol is what these volumes contain.

Writes, in adult_v2/arms/ under the data root:

    sep_channel_check.csv    the per-mouse numbers behind the figure
    sep_channel_check.png    the figure

Run by run_sep_channel_check.py.
"""

import csv
import math
import os
from collections import defaultdict

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from sepmap.volumes.cohort import NAIVE, RWS
from sepmap.volumes.per_mouse import CSV_MAP, DATA, MICE, annotation_20
from sepmap.volumes.per_mouse import OUT as PER_MOUSE

OUT = os.path.join(DATA, "adult_v2", "arms")
GENES = os.path.join(DATA, "adult_v2", "ish", "gene_region_table.csv")

# the ten adults, naive and rws pooled
ADULTS = NAIVE + RWS

# the three raw channels of the per-mouse files, and names for them that no
# figure of this module uses
CHANNELS = ("sig", "auto", "sep")
NICE = {"sig": "nano", "auto": "autofluo", "sep": "SEP (green)"}

# smallest structure kept, in 20 um voxels, as in young_vs_adult.region_plot
MIN_VOX = 250

# ISH voxels a structure needs for its Gria1 value to be used
MIN_ISH_VOXELS = 10

# the gene the channels are compared with
CONTROL = "Gria1"


def structure_names() -> dict[int, str]:
    """Structure name of each parcellation index."""
    names = {}
    with open(CSV_MAP, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["parcellation_term_set_name"] == "structure":
                names[int(row["parcellation_index"])] = row["parcellation_term_name"]
    return names


def mouse_channels(mouse: str, names: dict[int, str]) -> dict[str, dict[str, float]]:
    """{channel: {structure: log2 mean}} for one adult, raw, no denominators.

    The layer indices of a structure are pooled by name; a structure under MIN_VOX
    voxels, or with no signal, is left out.
    """
    z = np.load(os.path.join(PER_MOUSE, mouse + ".npz"))
    tissue = z["tissue"]
    annotation = annotation_20(MICE[mouse][1])
    labels = annotation[tissue]
    n = np.bincount(labels, minlength=int(annotation.max()) + 1)

    out = {}
    for key in CHANNELS:
        s = np.bincount(
            labels, weights=z[key].astype(np.float32)[tissue], minlength=len(n)
        )

        # pool the indices by structure name (index 0 is outside the brain)
        acc = defaultdict(lambda: [0, 0.0])
        for i in np.nonzero(n)[0]:
            if i == 0:
                continue
            a = acc[names.get(int(i), f"id{i}")]
            a[0] += int(n[i])
            a[1] += s[i]
        out[key] = {
            k: math.log2(t / c) for k, (c, t) in acc.items() if c >= MIN_VOX and t > 0
        }
    return out


def gria1_profile() -> dict[str, float]:
    """Gria1 ISH mean per structure, where the structure has MIN_ISH_VOXELS voxels."""
    prof = {}
    with open(GENES, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["symbol"] == CONTROL and int(r["n_voxels"]) >= MIN_ISH_VOXELS:
                prof[r["structure"]] = float(r["ish_mean"])
    return prof


def residual(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """log2(SEP) with the part predicted by log2(auto) taken out, linearly."""
    a = np.column_stack([x, np.ones_like(x)])
    return y - a @ np.linalg.lstsq(a, y, rcond=None)[0]


def main() -> None:
    """Measure the three channels per adult, write the table, print and draw it."""
    os.makedirs(OUT, exist_ok=True)

    # structure means of the three channels per adult
    names = structure_names()
    per = {}
    for mouse in ADULTS:
        per[mouse] = mouse_channels(mouse, names)
        print(f"{mouse:20s} {len(per[mouse]['sep'])} structures", flush=True)

    # per adult: range, correlations between channels and with Gria1, over the
    # structures all three channels have
    profile = gria1_profile()
    rows = []
    for mouse in ADULTS:
        ch = per[mouse]
        common = sorted(set.intersection(*[set(ch[k]) for k in CHANNELS]))
        v = {k: np.array([ch[k][s] for s in common]) for k in CHANNELS}
        with_gria1 = [s for s in common if s in profile]
        gria1 = np.array([profile[s] for s in with_gria1])
        idx = [common.index(s) for s in with_gria1]
        res = residual(v["sep"], v["auto"])
        row = dict(
            mouse=mouse,
            n_structures=len(common),
            range_nano=float(np.diff(np.percentile(v["sig"], [10, 90]))[0]),
            range_auto=float(np.diff(np.percentile(v["auto"], [10, 90]))[0]),
            range_sep=float(np.diff(np.percentile(v["sep"], [10, 90]))[0]),
            rho_sep_auto=float(spearmanr(v["sep"], v["auto"]).statistic),
            rho_sep_nano=float(spearmanr(v["sep"], v["sig"]).statistic),
            rho_nano_auto=float(spearmanr(v["sig"], v["auto"]).statistic),
            rho_nano_gria=float(spearmanr(v["sig"][idx], gria1).statistic),
            rho_auto_gria=float(spearmanr(v["auto"][idx], gria1).statistic),
            rho_sep_gria=float(spearmanr(v["sep"][idx], gria1).statistic),
            rho_sepresid_gria=float(spearmanr(res[idx], gria1).statistic),
            rho_sepresid_nano=float(spearmanr(res, v["sig"]).statistic),
        )
        rows.append(row)

    path = os.path.join(OUT, "sep_channel_check.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(
            {k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()}
            for r in rows
        )
    print(f"\n{len(rows)} mice -> {path}")

    def col(k):
        """One column of the table as an array, a value per adult."""
        return np.array([r[k] for r in rows])

    def say(k):
        """A column's mean +- SD over the adults, as printed."""
        return f"{np.mean(col(k)):+.3f} +- {np.std(col(k)):.3f}"

    # the summary over adults, mean +- SD
    print(
        f"\ndynamic range across structures, p90-p10 of log2, "
        f"mean over {len(rows)} adults"
    )
    for k, label in (
        ("range_nano", "nano"),
        ("range_auto", "autofluo"),
        ("range_sep", "SEP"),
    ):
        print(f"  {label:10s} {np.mean(col(k)):.2f} +- {np.std(col(k)):.2f}")
    print("\nwhat the green channel tracks (per mouse, over structures)")
    print(f"  SEP  ~ autofluo   {say('rho_sep_auto')}")
    print(f"  SEP  ~ nano       {say('rho_sep_nano')}")
    print(
        f"  nano ~ autofluo   {say('rho_nano_auto')}   "
        "<- the nano channel is its own thing"
    )
    print(f"\nagainst {CONTROL} expression")
    for k, label in (
        ("rho_nano_gria", "nano"),
        ("rho_auto_gria", "autofluo"),
        ("rho_sep_gria", "SEP"),
        ("rho_sepresid_gria", "SEP minus autofluo"),
    ):
        print(f"  {label:20s} {say(k)}")
    print(f"\n  SEP minus autofluo, against nano: {say('rho_sepresid_nano')}")

    # figure: the ranges, the green channel against the other two in the first
    # adult, and the correlations with Gria1; one jitter generator for the figure
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.0))
    rng = np.random.default_rng(0)

    ax = axes[0]
    for i, k in enumerate(("range_nano", "range_auto", "range_sep")):
        v = col(k)
        ax.scatter(
            np.full(len(v), i) + rng.uniform(-0.1, 0.1, len(v)),
            v,
            s=18,
            facecolor="#c0392b" if k == "range_sep" else "0.6",
            edgecolor="0.25",
            linewidth=0.4,
            zorder=2,
        )
        ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color="0.15", lw=1.7, zorder=3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["nano", "autofluo", "SEP"], fontsize=8)
    ax.set_ylabel("p90 - p10 across structures (log2)", fontsize=8)
    ax.set_ylim(bottom=0)
    ax.set_title("how much each channel varies\nacross the brain", fontsize=9)

    ax = axes[1]
    mouse = ADULTS[0]
    ch = per[mouse]
    common = sorted(set.intersection(*[set(ch[k]) for k in CHANNELS]))
    x = np.array([ch["auto"][s] for s in common])
    y = np.array([ch["sep"][s] for s in common])
    ax.scatter(x, y, s=10, facecolor="0.55", edgecolor="0.25", linewidth=0.3)
    ax.set_xlabel("log2 autofluorescence", fontsize=8)
    ax.set_ylabel("log2 SEP", fontsize=8)
    ax.set_title(
        f"{mouse}\nSEP against autofluo, rho = {spearmanr(x, y).statistic:+.2f}",
        fontsize=9,
    )

    ax = axes[2]
    nano = np.array([ch["sig"][s] for s in common])
    ax.scatter(nano, y, s=10, facecolor="0.55", edgecolor="0.25", linewidth=0.3)
    ax.set_xlabel("log2 nano", fontsize=8)
    ax.set_ylabel("log2 SEP", fontsize=8)
    ax.set_title(
        f"SEP against nano, rho = {spearmanr(nano, y).statistic:+.2f}", fontsize=9
    )

    ax = axes[3]
    keys = ("rho_nano_gria", "rho_auto_gria", "rho_sep_gria", "rho_sepresid_gria")
    labels = ["nano", "autofluo", "SEP", "SEP minus\nautofluo"]
    for i, k in enumerate(keys):
        v = col(k)
        ax.scatter(
            np.full(len(v), i) + rng.uniform(-0.1, 0.1, len(v)),
            v,
            s=18,
            facecolor="#c0392b" if "sep" in k else "0.6",
            edgecolor="0.25",
            linewidth=0.4,
            zorder=2,
        )
        ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color="0.15", lw=1.7, zorder=3)
    ax.axhline(0, color="0.8", lw=0.7, zorder=0)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(f"Spearman with {CONTROL} expression", fontsize=8)
    ax.set_title(
        "if the green channel were total receptor,\nit would beat nano here",
        fontsize=9,
    )

    for ax in axes:
        ax.tick_params(labelsize=7)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.suptitle(
        "The green channel in fixed, cleared tissue: one dot per adult, "
        "structure means with no denominator anywhere",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path = os.path.join(OUT, "sep_channel_check.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"\n{path}")
