"""What the green channel reports in this tissue: the tag, or the tissue (analysis 5).

Each adult brain was imaged in three channels: nano (the nanobody against the SEP
tag, on sections that were not permeabilised), SEP (the tag's own green
fluorescence) and auto (a channel with no label, the tissue's autofluorescence).
Reading receptor at the membrane against receptor anywhere from ratios of these
channels assumes that ex vivo the green channel reports the tagged receptor
wherever it sits. That premise is checked here on the raw channels, before any
ratio is taken: one mean per structure and adult, in log2, from the adult profiles
(adult.profiles), with no denominator anywhere. Three things decide it.

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
profile from the gene table (ish.gene_table, section QC applied).

Writes, in adult_v2/ish_analysis/ under the data root:

    green_channel/sep_channel_check.csv    per adult: the ranges and the correlations
    tables/numbers_green_channel.csv       the numbers of analysis 5, for the text

Run by run_sep_channel_check.py, which draws the working figure
(green_channel/sep_channel_check.png) and guided figures 15 and 15s.
"""

import math

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap.adult.beyond_density import residual
from sepmap.adult.profiles import ADULTS, CHANNELS, load_per_mouse
from sepmap.config import SETTINGS
from sepmap.ish import gene_table
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.structures import ISH_OUT, declared_structures

# the gene the channels are compared with (Gria1)
ISH = SETTINGS["ish"]

OUT = ISH_OUT / "green_channel"
TABLE = OUT / "sep_channel_check.csv"
NUMBERS = numbers_path("green_channel")


def channel_means(
    table: pd.DataFrame, subset: list[str]
) -> dict[str, dict[str, dict[str, float]]]:
    """{adult: {channel: {structure: log2 mean}}}, raw, no denominators.

    `table` is the per-adult table of adult.profiles (structures under
    region_tables.min_vox20 voxels already left out), cut to the structures of
    `subset`; a mean not above background is left out. The channels are named as
    the per-mouse files name them (sig for nano, auto, sep).
    """
    table = table[table["structure"].isin(set(subset))]
    out = {}
    for mouse in ADULTS:
        mine = table[table["mouse"] == mouse]
        out[mouse] = {
            key: {
                s: math.log2(v)
                for s, v in zip(mine["structure"], mine[f"{name}_mean"])
                if v > 0
            }
            for name, key in CHANNELS.items()
        }
    return out


def gria1_profile() -> dict[str, float]:
    """Gria1's merged profile from the gene table: {structure: mean rank}."""
    return gene_table.load_profiles()[ISH["control_gene"]]


def log2_range(values: np.ndarray) -> float:
    """p90 - p10 of log2 values: how much a channel varies across structures."""
    return float(np.diff(np.percentile(values, [10, 90]))[0])


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
        common = sorted(set.intersection(*[set(values) for values in ch.values()]))
        v = {k: np.array([ch[k][s] for s in common]) for k in ch}
        with_gria1 = [s for s in common if s in profile]
        gria1 = np.array([profile[s] for s in with_gria1])
        idx = [common.index(s) for s in with_gria1]
        res = residual(v["sep"], [v["auto"]])
        rows.append(
            dict(
                mouse=mouse,
                n_structures=len(common),
                range_nano=log2_range(v["sig"]),
                range_auto=log2_range(v["auto"]),
                range_sep=log2_range(v["sep"]),
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


def column(rows: list[dict], key: str) -> np.ndarray:
    """One column of the table as an array, a value per adult."""
    return np.array([r[key] for r in rows])


def mean_sd_range(rows: list[dict], key: str) -> str:
    """A column's mean +- SD over the adults, and its range, as printed."""
    v = column(rows, key)
    return f"{np.mean(v):+.3f} +- {np.std(v):.3f}  ({v.min():+.2f} to {v.max():+.2f})"


def print_summary(rows: list[dict]) -> None:
    """Print the summary over the adults."""
    print(f"\ndynamic range across structures, p90-p10 of log2, over {len(rows)} adults")
    for k, label in (
        ("range_nano", "nano"),
        ("range_auto", "autofluo"),
        ("range_sep", "SEP"),
    ):
        v = column(rows, k)
        print(f"  {label:10s} {np.mean(v):.2f} +- {np.std(v):.2f}")
    print("\nwhat the green channel tracks (per adult, over structures)")
    print(f"  SEP  ~ autofluo   {mean_sd_range(rows, 'rho_sep_auto')}")
    print(f"  SEP  ~ nano       {mean_sd_range(rows, 'rho_sep_nano')}")
    print(f"  nano ~ autofluo   {mean_sd_range(rows, 'rho_nano_auto')}")
    print(f"\nagainst {ISH['control_gene']} mRNA")
    for k, label in (
        ("rho_nano_gria", "nano"),
        ("rho_auto_gria", "autofluo"),
        ("rho_sep_gria", "SEP"),
        ("rho_sepresid_gria", "SEP minus autofluo"),
    ):
        print(f"  {label:20s} {mean_sd_range(rows, k)}")
    against_nano = mean_sd_range(rows, "rho_sepresid_nano")
    print(f"\n  SEP minus autofluo, against nano: {against_nano}")


def numbers_table(rows: pd.DataFrame, n_structures: int) -> pd.DataFrame:
    """numbers_green_channel.csv: the numbers of analysis 5 that the text quotes.

    `rows` is sep_channel_check.csv's table: per column, its mean, lowest and highest
    over the adults.
    """
    out = [
        ("adults", len(rows), "adults measured"),
        ("structures", n_structures, "declared structures"),
    ]
    for name in rows.columns:
        if name in ("mouse", "n_structures"):
            continue
        v = rows[name]
        out += [
            (f"{name}_mean", round(float(v.mean()), 3), f"{name}, mean of the adults"),
            (f"{name}_min", round(float(v.min()), 3), f"{name}, lowest adult"),
            (f"{name}_max", round(float(v.max()), 3), f"{name}, highest adult"),
        ]
    return numbers_frame(out)


def main_inputs() -> tuple[dict, dict[str, float], list[str]]:
    """The channel means on the declared structures, Gria1's profile, the structures."""
    structures = declared_structures()
    per = channel_means(load_per_mouse(), structures)
    return per, gria1_profile(), structures


def main() -> tuple[dict, list[dict]]:
    """Measure the three channels per adult, write the table and print it.

    Returns the channel means and the rows, for the figures.
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
    return per, rows
