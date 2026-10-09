"""What the green channel reports in this tissue: the tag, or the tissue (analysis 5).

Each adult brain was imaged in three channels: nano (the nanobody against the SEP
tag, on sections that were not permeabilised), SEP (the tag's own green
fluorescence) and auto (a channel with no label, the tissue's autofluorescence).
The three-channel design of September read receptor at the membrane against
receptor anywhere from ratios of these channels, and assumed that ex vivo the green
channel reports the tagged receptor wherever it sits. That premise is checked here
on the raw channels, before any ratio is taken: one mean per structure and adult, in
log2, from the adult profiles (adult.profiles), with no denominator anywhere. Three
things decide it.

    dynamic range   how much a channel varies across the brain (p90 - p10 of log2
                    over the declared structures). Receptor is not flat
                    (hippocampus against thalamus is a large, well known difference),
                    so a channel reporting it should have a range comparable to the
                    nanobody's.
    who it tracks   if the green channel is the tag, its profile across structures
                    should look like the nanobody's and like Gria1 mRNA's; if it is
                    the tissue, like the autofluorescence channel's.
    what is left    log2(SEP) regressed on log2(auto), and the residual asked
                    whether it still follows the nanobody and Gria1: how much tag
                    survives under the autofluorescence.

The answer is a property of this tissue, fixed and mounted, neither cleared nor
permeabilised, not of SEP-GluA1 mice. Superecliptic pHluorin is quenched in acidic
compartments in a living cell, which makes it a surface reporter in vivo; after
fixation the pH gradients are gone, and whatever green emission survives the
protocol is what these volumes contain. Whether surface and all receptor differ
across the brain needs a total-GluA1 stain on the same brains; these data cannot say.

The inputs are those of the rest of the ISH line: the declared structures
(sepmap.structures), the channel means of adult.profiles, and Gria1's merged
profile from the gene table (ish.gene_table, section QC applied). The September
version took every structure each adult measured and Gria1's single experiment of
the 100-gene panel; the arithmetic is unchanged.

Writes, in adult_v2/ish_analysis/green_channel/ under the data root:

    sep_channel_check.csv    per adult: the ranges and the correlations
    sep_channel_check.png    the working figure (guided figure 17 is the one to show)

Run by run_sep_channel_check.py, which also draws guided figure 17.
"""

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap.adult import profiles
from sepmap.adult.beyond_density import residual
from sepmap.config import DATA, SETTINGS
from sepmap.ish import gene_table
from sepmap.plotting import RED, tidy
from sepmap.structures import declared_structures

# the gene the channels are compared with (Gria1)
ISH = SETTINGS["ish"]

OUT = DATA / "adult_v2" / "ish_analysis" / "green_channel"
TABLE = OUT / "sep_channel_check.csv"

# the ten adults, naive and rws pooled
ADULTS = profiles.ADULTS

# the three raw channels, by the name the per-mouse files give them, and the column
# of the adult profiles that holds each
CHANNELS = ("sig", "auto", "sep")
COLUMNS = {"sig": "nano_mean", "auto": "auto_mean", "sep": "sep_mean"}


def channel_means(
    table: pd.DataFrame, structures: list[str] | None = None
) -> dict[str, dict[str, dict[str, float]]]:
    """{adult: {channel: {structure: log2 mean}}}, raw, no denominators.

    `table` is the per-adult table of adult.profiles (structures under
    region_tables.min_vox20 voxels already left out); a mean not above background
    is left out. With `structures`, only those.
    """
    if structures is not None:
        table = table[table["structure"].isin(set(structures))]
    out = {}
    for mouse in ADULTS:
        mine = table[table["mouse"] == mouse]
        out[mouse] = {
            key: {
                s: math.log2(v)
                for s, v in zip(mine["structure"], mine[COLUMNS[key]])
                if v > 0
            }
            for key in CHANNELS
        }
    return out


def gria1_profile() -> dict[str, float]:
    """Gria1's merged profile from the gene table: {structure: mean rank}."""
    return gene_table.load_profiles()[ISH["control_gene"]]


def channel_rows(
    per: dict[str, dict[str, dict[str, float]]], profile: dict[str, float]
) -> list[dict]:
    """Per adult: the channels' ranges, their correlations, and with Gria1.

    Over the structures all three channels have in that adult; one row of
    sep_channel_check.csv per adult.
    """
    rows = []
    for mouse in ADULTS:
        ch = per[mouse]
        common = sorted(set.intersection(*[set(ch[k]) for k in CHANNELS]))
        v = {k: np.array([ch[k][s] for s in common]) for k in CHANNELS}
        with_gria1 = [s for s in common if s in profile]
        gria1 = np.array([profile[s] for s in with_gria1])
        idx = [common.index(s) for s in with_gria1]
        res = residual(v["sep"], [v["auto"]])
        rows.append(
            dict(
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
        )
    return rows


def write_table(rows: list[dict]) -> None:
    """Write sep_channel_check.csv, floats to four decimals, and say where."""
    table = pd.DataFrame(rows)
    table.to_csv(TABLE, index=False, float_format="%.4f")
    print(f"\n{len(rows)} adults -> {TABLE}")


def load_table() -> pd.DataFrame:
    """The table run_sep_channel_check wrote."""
    if not TABLE.exists():
        raise FileNotFoundError(f"{TABLE} not found: run run_sep_channel_check.py first")
    return pd.read_csv(TABLE)


def col(rows: list[dict], k: str) -> np.ndarray:
    """One column of the table as an array, a value per adult."""
    return np.array([r[k] for r in rows])


def say(rows: list[dict], k: str) -> str:
    """A column's mean +- SD over the adults, and its range, as printed."""
    v = col(rows, k)
    return f"{np.mean(v):+.3f} +- {np.std(v):.3f}  ({v.min():+.2f} to {v.max():+.2f})"


def print_summary(rows: list[dict]) -> None:
    """Print the summary over the adults."""
    print(f"\ndynamic range across structures, p90-p10 of log2, over {len(rows)} adults")
    for k, label in (
        ("range_nano", "nano"),
        ("range_auto", "autofluo"),
        ("range_sep", "SEP"),
    ):
        print(f"  {label:10s} {np.mean(col(rows, k)):.2f} +- {np.std(col(rows, k)):.2f}")
    print("\nwhat the green channel tracks (per adult, over structures)")
    print(f"  SEP  ~ autofluo   {say(rows, 'rho_sep_auto')}")
    print(f"  SEP  ~ nano       {say(rows, 'rho_sep_nano')}")
    print(f"  nano ~ autofluo   {say(rows, 'rho_nano_auto')}")
    print(f"\nagainst {ISH['control_gene']} mRNA")
    for k, label in (
        ("rho_nano_gria", "nano"),
        ("rho_auto_gria", "autofluo"),
        ("rho_sep_gria", "SEP"),
        ("rho_sepresid_gria", "SEP minus autofluo"),
    ):
        print(f"  {label:20s} {say(rows, k)}")
    print(f"\n  SEP minus autofluo, against nano: {say(rows, 'rho_sepresid_nano')}")


def panel_ranges(ax: plt.Axes, rows: list[dict], rng: np.random.Generator) -> None:
    """Draw how much each channel varies across the brain, a dot per adult."""
    for i, k in enumerate(("range_nano", "range_auto", "range_sep")):
        v = col(rows, k)
        ax.scatter(
            np.full(len(v), i) + rng.uniform(-0.1, 0.1, len(v)),
            v,
            s=18,
            facecolor=RED if k == "range_sep" else "0.6",
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


def panels_first_adult(
    axes: np.ndarray, per: dict[str, dict[str, dict[str, float]]]
) -> None:
    """Draw the green channel against autofluorescence and against nano, first adult."""
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


def panel_gria1(ax: plt.Axes, rows: list[dict], rng: np.random.Generator) -> None:
    """Draw each channel's correlation with Gria1 mRNA, a dot per adult."""
    keys = ("rho_nano_gria", "rho_auto_gria", "rho_sep_gria", "rho_sepresid_gria")
    labels = ["nano", "autofluo", "SEP", "SEP minus\nautofluo"]
    for i, k in enumerate(keys):
        v = col(rows, k)
        ax.scatter(
            np.full(len(v), i) + rng.uniform(-0.1, 0.1, len(v)),
            v,
            s=18,
            facecolor=RED if "sep" in k else "0.6",
            edgecolor="0.25",
            linewidth=0.4,
            zorder=2,
        )
        ax.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color="0.15", lw=1.7, zorder=3)
    ax.axhline(0, color="0.8", lw=0.7, zorder=0)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(f"Spearman with {ISH['control_gene']} mRNA", fontsize=8)
    ax.set_title(
        "a channel that reports the tagged receptor\nshould follow Gria1 at least as "
        "nano does",
        fontsize=9,
    )


def figure(per: dict[str, dict[str, dict[str, float]]], rows: list[dict]) -> None:
    """Draw the working figure: the ranges, the first adult, the Gria1 correlations.

    One jitter generator for the figure, drawn from in panel order.
    """
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.0))
    rng = np.random.default_rng(0)
    panel_ranges(axes[0], rows, rng)
    panels_first_adult(axes, per)
    panel_gria1(axes[3], rows, rng)
    for ax in axes:
        tidy(ax)
    fig.suptitle(
        "The green channel in fixed, mounted tissue: one dot per adult, "
        "structure means with no denominator anywhere",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path = OUT / "sep_channel_check.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"\n{path}")


def main_inputs() -> tuple[dict, dict[str, float], list[str]]:
    """The channel means on the declared structures, Gria1's profile, the structures."""
    structures = declared_structures()
    per = channel_means(profiles.load_per_mouse(), structures)
    return per, gria1_profile(), structures


def main() -> tuple[dict, list[dict]]:
    """Measure the three channels per adult, write the table, print and draw it.

    Returns the channel means and the rows, for guided figure 17.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    per, profile, structures = main_inputs()
    print(
        f"{len(structures)} declared structures, {len(ADULTS)} adults; "
        f"{ISH['control_gene']} has a profile in "
        f"{sum(s in profile for s in structures)} of them"
    )
    rows = channel_rows(per, profile)
    write_table(rows)
    print_summary(rows)
    figure(per, rows)
    return per, rows
