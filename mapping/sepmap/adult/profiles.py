"""Per adult and structure: the three channels, plain and eroded, and the adult profile.

Every analysis of the ISH line compares one adult profile, a value per structure,
with the Allen maps. It is built here once, from the per-mouse files of
volumes.per_mouse, so the gene ranking, the autofluorescence control, the green
channel check and the beyond analysis all read the same numbers:

    channel means   per adult and structure, the mean of nano (sig),
                    autofluorescence and SEP over the tissue voxels, in counts
                    above each channel's off-tissue background; the layers of an
                    area pooled under its name; structures under
                    region_tables.min_vox20 voxels left out, as in region_plot
    eroded means    the same over the structure's voxels that have no
                    six-neighbour of another structure (or outside the brain): one
                    20 um voxel peeled off every border. A structure with fewer than
                    min_vox20 voxels left keeps its plain mean, and says so
                    (eroded_is_plain). The registration is good to a few voxels, so
                    the border of a small structure is partly its neighbour's
                    signal; the eroded mean shows whether that matters
    zref            per channel, with the declared set as reference
                    (structures.zref): nano over the brain's own isocortex mean of
                    nano, autofluorescence over its own isocortex mean of
                    autofluorescence. The eroded zref uses the plain reference, so
                    only the structure means change between the two
    profile         per structure, the mean over the adults that measure it

The channel means take about ten minutes for the ten adults, so they are cached
in adult_per_mouse.csv beside the zref columns; run_structure_set.py --recompute
redoes them.

The stored zref of young_vs_adult/region_means_per_mouse.csv used the 17 brains'
shared structures as reference. stored_reference recovers each adult's median and
spread from that table (cref = median + spread x zref, exact to the four decimals
written), so figure 01 can show how far the declared reference moves them.

Run by run_structure_set.py.
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap.config import SETTINGS
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.structures import TABLES, zref
from sepmap.volumes.cohort import NAIVE, RWS
from sepmap.volumes.per_mouse import OUT as PER_MOUSE
from sepmap.young_vs_adult.region_plot import REGION_MEANS

# the smallest structure a brain's table keeps; the adults a declared structure needs
REGION_TABLES = SETTINGS["region_tables"]
STRUCTURES = SETTINGS["structures"]

PER_MOUSE_TABLE = TABLES / "adult_per_mouse.csv"
PROFILE = TABLES / "adult_profile.csv"
REFERENCE = TABLES / "zref_reference.csv"
NUMBERS = numbers_path("structure_set")

# the ten adults, naive then rws, with their group
ADULTS = NAIVE + RWS
GROUP = {**{m: "naive" for m in NAIVE}, **{m: "rws" for m in RWS}}

# the channels of the per-mouse files, by the name the tables give them
CHANNELS = {"nano": "sig", "auto": "auto", "sep": "sep"}

# the structures left out are counted in two groups: those seen in fewer than this
# many adults, the few-mice structures of the hindbrain, and those seen in more but
# not in all
FEW_ADULTS = 5

# the columns of the cache, in order
MEAN_COLUMNS = [
    "mouse",
    "group",
    "structure",
    "acronym",
    "division",
    "n_vox20",
    "n_vox20_eroded",
    "eroded_is_plain",
] + [f"{c}_mean{e}" for e in ("", "_eroded") for c in CHANNELS]


def structure_meta(
    names: dict[int, str], acronyms: dict[int, str], divisions: dict[int, str]
) -> dict[str, tuple[str, str]]:
    """{structure name: (acronym, division)}, as region_plot builds it."""
    return {name: (acronyms[idx], divisions.get(idx, "")) for idx, name in names.items()}


def mouse_rows(
    mouse: str,
    labels: np.ndarray,
    interior: np.ndarray,
    label_names: list[str],
    meta: dict[str, tuple[str, str]],
) -> list[dict]:
    """The plain and eroded channel means of one adult, one row per structure.

    `labels` holds the structure codes of the CCF on the 20 um grid
    (structures.name_volume) and `interior` its voxels with no border.
    """
    z = np.load(PER_MOUSE / (mouse + ".npz"))
    tissue = z["tissue"]
    codes = labels[tissue]
    inner = interior[tissue]
    n_labels = len(label_names)

    # voxels and channel sums per structure, plain and eroded
    n = np.bincount(codes, minlength=n_labels)
    n_eroded = np.bincount(codes[inner], minlength=n_labels)
    sums, sums_eroded = {}, {}
    for name, key in CHANNELS.items():
        values = z[key].astype(np.float32)[tissue]
        sums[name] = np.bincount(codes, weights=values, minlength=n_labels)
        sums_eroded[name] = np.bincount(
            codes[inner], weights=values[inner], minlength=n_labels
        )

    rows = []
    for code in np.nonzero(n >= REGION_TABLES["min_vox20"])[0]:
        if code == 0:
            continue
        structure = label_names[code]
        acronym, division = meta.get(structure, ("", ""))

        # a structure eroded below the smallest size kept keeps its plain mean
        plain = n_eroded[code] < REGION_TABLES["min_vox20"]
        row = dict(
            mouse=mouse,
            group=GROUP[mouse],
            structure=structure,
            acronym=acronym,
            division=division,
            n_vox20=int(n[code]),
            n_vox20_eroded=int(n_eroded[code]),
            eroded_is_plain=bool(plain),
        )
        for name in CHANNELS:
            row[f"{name}_mean"] = sums[name][code] / n[code]
        for name in CHANNELS:
            if plain:
                row[f"{name}_mean_eroded"] = row[f"{name}_mean"]
            else:
                row[f"{name}_mean_eroded"] = sums_eroded[name][code] / n_eroded[code]
        rows.append(row)
    return rows


def channel_table(
    labels: np.ndarray, interior: np.ndarray, label_names: list[str], meta: dict
) -> pd.DataFrame:
    """Plain and eroded channel means of the ten adults, a row per adult and structure."""
    rows = []
    for mouse in ADULTS:
        mine = mouse_rows(mouse, labels, interior, label_names, meta)
        print(f"  {mouse:16s} {len(mine)} structures", flush=True)
        rows.extend(mine)
    return pd.DataFrame(rows, columns=MEAN_COLUMNS)


def read_table(path: Path) -> pd.DataFrame:
    """A per-adult table as written, eroded_is_plain back to a boolean."""
    table = pd.read_csv(
        path, keep_default_na=False, na_values=[""], float_precision="round_trip"
    )
    table["eroded_is_plain"] = table["eroded_is_plain"].astype(str) == "True"
    return table


def check_channel_table(table: pd.DataFrame) -> None:
    """Stop when the table's adults are not the ten of the cohort table, in order."""
    mice = list(dict.fromkeys(table["mouse"]))
    if mice != ADULTS:
        raise ValueError(
            f"{PER_MOUSE_TABLE} holds the adults {mice}, but the cohort table has "
            f"{ADULTS}. Run run_structure_set.py --recompute."
        )


def load_channel_table() -> pd.DataFrame:
    """The channel means cached in adult_per_mouse.csv, checked against the adults."""
    return load_per_mouse()[MEAN_COLUMNS]


def isocortex_mean(rows: pd.DataFrame, channel: str) -> float:
    """Voxel-weighted mean of a channel over a brain's isocortex structures.

    `rows` is one brain's part of the channel table; as region_plot.ref, over the
    structures the table keeps (at least region_tables.min_vox20 voxels).
    """
    iso = rows[rows["division"] == "Isocortex"]
    return float((iso[f"{channel}_mean"] * iso["n_vox20"]).sum() / iso["n_vox20"].sum())


def log2_over(values: pd.Series, reference: float) -> pd.Series:
    """log2 of `values` over `reference`, NaN where a value is not above background."""
    return values.where(values > 0).div(reference).map(math.log2, na_action="ignore")


def with_zref(table: pd.DataFrame, structures: list[str]) -> tuple[pd.DataFrame, list]:
    """The channel table with cref and zref columns, and each brain's reference.

    Adds, per row, cref_nano and cref_auto (log2 over the brain's own isocortex
    mean of that channel), zref_nano, zref_auto and their eroded versions (the
    eroded means over the same isocortex mean, scaled by the plain reference).
    Returns the new table and one dict per adult and channel: median, spread and
    the structures they were taken over.
    """
    out = table.copy()
    references = []
    for mouse in ADULTS:
        mine = out["mouse"] == mouse
        rows = out[mine].set_index("structure")
        for channel in ("nano", "auto"):
            iso = isocortex_mean(rows, channel)
            v = log2_over(rows[f"{channel}_mean"], iso)
            v_eroded = log2_over(rows[f"{channel}_mean_eroded"], iso)
            z, median, spread = zref(v, structures)
            out.loc[mine, f"cref_{channel}"] = v.to_numpy()
            out.loc[mine, f"zref_{channel}"] = z.to_numpy()
            out.loc[mine, f"zref_{channel}_eroded"] = (
                (v_eroded - median) / spread
            ).to_numpy()
            references.append(
                dict(
                    mouse=mouse,
                    group=GROUP[mouse],
                    channel=channel,
                    isocortex_mean=iso,
                    median=median,
                    spread=spread,
                    n_reference=int(v.reindex(structures).notna().sum()),
                )
            )
    return out, references


def stored_reference() -> dict[str, tuple[float, float]]:
    """{adult: (median, spread)} of the stored zref of region_plot.

    The stored table has cref and zref per structure, and zref = (cref - median) /
    spread, so a least-squares line of cref on zref gives the median (intercept)
    and the spread (slope) the young-against-adult tables used.
    """
    table = pd.read_csv(REGION_MEANS)
    table = table[table["mouse"].isin(ADULTS)]
    out = {}
    for mouse in ADULTS:
        mine = table[table["mouse"] == mouse]
        cref = mine[mine["reading"] == "cref"].set_index("structure")["log2_value"]
        z = mine[mine["reading"] == "zref"].set_index("structure")["log2_value"]
        both = pd.concat([cref.rename("cref"), z.rename("zref")], axis=1).dropna()
        design = np.column_stack([np.ones(len(both)), both["zref"]])
        coef, *_ = np.linalg.lstsq(design, both["cref"].to_numpy(), rcond=None)
        out[mouse] = (float(coef[0]), float(coef[1]))
    return out


def reference_table(references: list[dict]) -> pd.DataFrame:
    """Each adult's nano and autofluorescence reference, the stored one beside nano."""
    stored = stored_reference()
    rows = []
    for r in references:
        row = dict(r)
        if r["channel"] == "nano":
            row["median_stored"], row["spread_stored"] = stored[r["mouse"]]
        else:
            row["median_stored"] = row["spread_stored"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def adult_profile(per_mouse: pd.DataFrame, set_table: pd.DataFrame) -> pd.DataFrame:
    """Per structure of the adult table, the mean over the adults that measure it.

    Columns: structure, acronym, division, in_set, n_adults (with a nano zref),
    then the mean of zref_nano, zref_auto, their eroded versions, cref_nano and
    cref_auto, and the SD across adults of zref_nano. In the declared set every
    structure has all ten adults.
    """
    value_columns = [
        "zref_nano",
        "zref_auto",
        "zref_nano_eroded",
        "zref_auto_eroded",
        "cref_nano",
        "cref_auto",
    ]
    grouped = per_mouse.groupby("structure")
    means = grouped[value_columns].mean()
    out = set_table[["structure", "acronym", "division", "in_set"]].set_index("structure")
    out["n_adults"] = (
        grouped["zref_nano"].count().reindex(out.index).fillna(0).astype(int)
    )
    out = out.join(means)
    out["zref_nano_sd"] = grouped["zref_nano"].std()
    return out.reset_index()


def load_profile() -> pd.DataFrame:
    """The adult profile run_structure_set wrote, indexed by structure."""
    if not PROFILE.exists():
        raise FileNotFoundError(f"{PROFILE} not found: run run_structure_set.py first")
    table = pd.read_csv(PROFILE, keep_default_na=False, na_values=[""])
    table["in_set"] = table["in_set"].astype(str) == "True"
    return table.set_index("structure")


def load_per_mouse() -> pd.DataFrame:
    """The per-adult table run_structure_set wrote, with its zref columns."""
    if not PER_MOUSE_TABLE.exists():
        raise FileNotFoundError(
            f"{PER_MOUSE_TABLE} not found: run run_structure_set.py first"
        )
    table = read_table(PER_MOUSE_TABLE)
    check_channel_table(table)
    return table


def stored_agreement(profile: pd.DataFrame) -> float:
    """Spearman of the adults' mean stored zref with the new one, over the set.

    The stored zref is region_plot's, with the 17 brains' shared structures as
    reference; `profile` is the adult profile, indexed by structure.
    """
    table = pd.read_csv(REGION_MEANS)
    stored = table[(table["reading"] == "zref") & table["mouse"].isin(ADULTS)]
    stored = stored.groupby("structure")["log2_value"].mean()
    in_set = profile[profile["in_set"]]
    return float(spearmanr(stored.reindex(in_set.index), in_set["zref_nano"]).statistic)


def numbers_table(
    set_table: pd.DataFrame,
    reference: pd.DataFrame,
    centroids: pd.DataFrame,
    agreement: float,
) -> pd.DataFrame:
    """numbers_structure_set.csv: the numbers of step 13 that the text quotes.

    `set_table` is structure_set.csv, `reference` each adult's zref reference
    (reference_table), `agreement` the cohort map's order before and after the
    declared reference (stored_agreement).
    """
    nano = reference[reference["channel"] == "nano"]
    shift = (nano["median"] - nano["median_stored"]).abs()
    n_adults = set_table["n_adults"]
    min_adults = STRUCTURES["min_adults"]
    rows = [
        ("structures_in_table", len(set_table), "structures in the adult table"),
        (
            "structures_all_adults",
            int((n_adults >= min_adults).sum()),
            "measured in all the adults the rule asks for",
        ),
        (
            f"structures_below_{FEW_ADULTS}_adults",
            int((n_adults < FEW_ADULTS).sum()),
            "left out: measured in fewer than five adults",
        ),
        (
            f"structures_{FEW_ADULTS}_to_9_adults",
            int(((n_adults >= FEW_ADULTS) & (n_adults < min_adults)).sum()),
            "left out: measured in five or more adults, but not in all",
        ),
        ("structures_declared", int(set_table["in_set"].sum()), "the declared set"),
        (
            "zref_zero_shift_min",
            round(shift.min(), 4),
            "smallest shift of an adult's zero",
        ),
        (
            "zref_zero_shift_max",
            round(shift.max(), 4),
            "largest shift of an adult's zero",
        ),
        (
            "cohort_map_agreement",
            round(agreement, 5),
            "Spearman of the cohort map before and after A1, declared set",
        ),
        (
            "centroids_ml_max_mm",
            round(float(centroids["ml_mm"].max()), 3),
            "largest ML centroid, mm (one hemisphere: below 5.7)",
        ),
    ]
    in_set = set_table[set_table["in_set"]]
    for division, n in in_set["division"].value_counts().items():
        rows.append(
            (f"declared_{division}", int(n), f"declared structures in {division}")
        )
    return numbers_frame(rows)
