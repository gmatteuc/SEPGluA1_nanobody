"""The gene ranking under other reasonable choices, against the primary (A3).

The primary ranking (ish.gene_ranking) is one set of choices: Spearman, the nano
map as zref with the declared reference, each gene's merged profile after section
QC, full structure means, the declared structures. Each row here
changes one of them, or several where a row reproduces an earlier route, and asks
whether the order of the genes, and where Cacng8 and Gria1 sit, moves:

    statistic   Pearson on log2: the map's zref against the log2 of the gene's
                energy (ish_analysis.log2_floor added first), each experiment
                centred on its median structure before a gene's experiments are
                averaged, since an experiment's scale is its own
    borders     eroded means: the ISH side by one 200 um voxel (ish.regions), the
                nano side by one 20 um voxel (adult.profiles), and both
    reading     the stored ratio reading (nano over autofluorescence) and the
                stored zref of the young-against-adult tables (the 17 brains'
                shared structures as reference), each the mean over the ten adults
    inputs      no section QC (the flagged sections measured as they are); P9's
                single experiment per gene instead of the merged profile
    structures  P9's nine divisions (the structures of those divisions in the adult
                table, seen in any number of adults, no catch-all labels), and
                every structure of the adult table
    5 October   the route as it ran on 5 October: stored zref, every structure, P9's
                one experiment per gene, no QC (the region table of
                adult_v2/ish/, frozen), as ish.compare correlated them

Each variant keeps the primary's rule for a gene's structures (ish.min_voxels
voxels of data, ish.min_structures structures). Its agreement with the primary is
the Spearman of the genes' rho, over P9's genes (the genes every row has) and over
all genes the row has. No row carries a p value: the null is drawn on the primary
map's structures, and the rows ask about the order, not about significance.

Run by run_ish_robustness.py.
"""

from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from sepmap.config import DATA, SETTINGS
from sepmap.ish import gene_ranking, gene_table, section_qc
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.ish.reliability import merge
from sepmap.structures import TABLES
from sepmap.young_vs_adult.region_plot import REGION_MEANS

# the voxels a gene value needs, the structures a correlation needs; P9's divisions;
# the floor of the log2 row
ISH = SETTINGS["ish"]
ISH_ANALYSIS = SETTINGS["ish_analysis"]
STRUCTURES = SETTINGS["structures"]

ROBUSTNESS = TABLES / "ranking_robustness.csv"
SUMMARY = TABLES / "robustness_summary.csv"
NUMBERS = numbers_path("robustness")

# the region table of 5 October, frozen: P9's panel, one experiment per gene, no QC
FROZEN_REGIONS = DATA / "adult_v2" / "ish" / "gene_region_table.csv"

# the adult groups of the stored tables
ADULT_GROUPS = ("naive", "rws")

# the rows of the table, in the order the figure draws them: name, kind of choice,
# label
VARIANTS = (
    ("primary", "primary", "primary: Spearman, zref (A1), QC, merged, declared set"),
    ("pearson_log2", "statistic", "Pearson on log2"),
    ("eroded_ish", "borders", "ISH means eroded"),
    ("eroded_nano", "borders", "nano means eroded"),
    ("eroded_both", "borders", "both eroded"),
    ("ratio", "reading", "ratio reading (nano / autofluorescence)"),
    ("stored_zref", "reading", "zref before A1 (17 brains' reference)"),
    ("no_qc", "inputs", "no section QC"),
    ("p9_experiment", "inputs", "P9's single experiment per gene"),
    ("p9_divisions", "structures", "P9's nine divisions, any number of adults"),
    ("every_structure", "structures", "every structure of the adult table"),
    ("before_build", "5 October", "the route of 5 October (stored zref, all, no QC)"),
)


# ===== Gene profiles =====


def log2_profiles(
    per_experiment: dict[str, dict[str, dict[str, float]]],
) -> dict[str, dict[str, float]]:
    """Each gene's mean over its experiments of log2 energy, centred per experiment.

    log2(energy + ish_analysis.log2_floor) minus the experiment's median over its
    structures, so experiments of one gene with different scales can be averaged.
    """
    floor = ISH_ANALYSIS["log2_floor"]
    out = {}
    for gene, experiments in per_experiment.items():
        acc = defaultdict(list)
        for values in experiments.values():
            keys = sorted(values)
            logs = np.log2(np.array([values[k] for k in keys]) + floor)
            logs = logs - np.median(logs)
            for k, v in zip(keys, logs):
                acc[k].append(v)
        out[gene] = {k: float(np.mean(v)) for k, v in acc.items()}
    return out


def merge_experiments(
    per_experiment: dict[str, dict[str, dict[str, float]]],
) -> dict[str, dict[str, float]]:
    """Each gene's merged profile, the mean of its experiments' 0-1 ranks."""
    return {gene: merge(experiments)[0] for gene, experiments in per_experiment.items()}


def single_experiment(
    per_experiment: dict[str, dict[str, dict[str, float]]], chosen: dict[str, str]
) -> dict[str, dict[str, float]]:
    """{gene: profile of the experiment `chosen` for it}, for the genes that have it."""
    out = {}
    for gene, eid in chosen.items():
        if eid and eid in per_experiment.get(gene, {}):
            out[gene] = per_experiment[gene][eid]
    return out


def unflagged_region(region: pd.DataFrame, recomputed: pd.DataFrame) -> pd.DataFrame:
    """The region table with the flagged experiments' rows measured without QC.

    `recomputed` holds the region rows of the experiments that had a section set
    missing, measured again with nothing set missing; every other experiment is
    the same with or without QC, so its rows are kept.
    """
    redone = set(recomputed["experiment_id"])
    kept = region[~region["experiment_id"].isin(redone)]
    return pd.concat([kept, recomputed], ignore_index=True)


def region_without_qc(table: pd.DataFrame, region: pd.DataFrame) -> pd.DataFrame:
    """The region table with the flagged experiments measured again, nothing missing.

    `table` is the gene table, `region` the region table of the experiments it uses.
    """
    flags = section_qc.flagged_sections(section_qc.load_experiment_qc())
    flagged = table[table["experiment_id"].isin(flags) & ~table["excluded"]]
    redone = gene_table.region_table(flagged, {})
    print(f"no QC: {len(flagged)} flagged experiments measured again", flush=True)
    return unflagged_region(region, redone)


def frozen_profiles() -> dict[str, dict[str, float]]:
    """{gene: {structure: energy}} of the frozen 5 October table, read as compare did."""
    table = pd.read_csv(FROZEN_REGIONS)
    usable = table[table["n_voxels"] >= ISH["min_voxels"]]
    out = defaultdict(dict)
    for symbol, structure, value in zip(
        usable["symbol"], usable["structure"], usable["ish_mean"]
    ):
        out[symbol][structure] = float(value)
    return dict(out)


# ===== Maps and structure sets =====


def stored_map(reading: str) -> pd.Series:
    """A stored reading of region_plot, the mean over the adults that have it."""
    table = pd.read_csv(REGION_MEANS)
    mine = table[table["group"].isin(ADULT_GROUPS) & (table["reading"] == reading)]
    return mine.groupby("structure")["log2_value"].mean()


def p9_division_structures(set_table: pd.DataFrame) -> list[str]:
    """The adult table's structures in P9's nine divisions, catch-all labels left out."""
    mine = set_table[
        set_table["division"].isin(STRUCTURES["p9_divisions"])
        & ~set_table["structure"].str.lower().str.contains("unassigned")
    ]
    return sorted(mine["structure"])


def variant_inputs(
    table: pd.DataFrame,
    genes: pd.DataFrame,
    profile: pd.DataFrame,
    set_table: pd.DataFrame,
    declared: list[str],
) -> dict[str, tuple[pd.Series, dict, list[str], str]]:
    """{variant: (map, gene profiles, structures, statistic)} of every row of VARIANTS.

    `table` is the gene table, `genes` its one row per gene, `profile` the adult
    profile and `set_table` the structure set; the row without section QC measures
    the flagged experiments again (about a minute).
    """
    region = gene_table.load_region_table()
    used = set(table.loc[~table["excluded"], "experiment_id"])
    region = region[region["experiment_id"].isin(used)]
    per_experiment = gene_table.experiment_profiles(region)
    eroded = merge_experiments(gene_table.experiment_profiles(region, "ish_mean_eroded"))
    no_qc = merge_experiments(
        gene_table.experiment_profiles(region_without_qc(table, region))
    )
    p9_own = dict(zip(genes["symbol"], genes["p9_experiment_id"]))
    merged = gene_table.load_profiles()
    zref = profile["zref_nano"]
    zref_eroded = profile["zref_nano_eroded"]
    stored_zref = stored_map("zref")
    every = sorted(set_table["structure"])
    inputs = {
        "primary": (zref, merged, declared, "spearman"),
        "pearson_log2": (zref, log2_profiles(per_experiment), declared, "pearson"),
        "eroded_ish": (zref, eroded, declared, "spearman"),
        "eroded_nano": (zref_eroded, merged, declared, "spearman"),
        "eroded_both": (zref_eroded, eroded, declared, "spearman"),
        "ratio": (stored_map("ratio"), merged, declared, "spearman"),
        "stored_zref": (stored_zref, merged, declared, "spearman"),
        "no_qc": (zref, no_qc, declared, "spearman"),
        "p9_experiment": (
            zref,
            single_experiment(per_experiment, p9_own),
            declared,
            "spearman",
        ),
        "p9_divisions": (zref, merged, p9_division_structures(set_table), "spearman"),
        "every_structure": (zref, merged, every, "spearman"),
        "before_build": (
            stored_zref,
            frozen_profiles(),
            sorted(stored_zref.index),
            "spearman",
        ),
    }
    missing = [name for name, _, _ in VARIANTS if name not in inputs]
    if missing:
        raise ValueError(f"VARIANTS lists {missing}, which variant_inputs does not build")
    return inputs


# ===== Correlations =====


def correlate(
    map_values: pd.Series,
    profiles: dict[str, dict[str, float]],
    structures: list[str],
    method: str = "spearman",
) -> pd.DataFrame:
    """rho of every gene with the map over `structures`, one row per gene.

    Only the structures where both have a value enter; a gene with fewer than
    ish.min_structures of them gets no row. `method` is spearman or pearson.
    """
    statistic = {"spearman": spearmanr, "pearson": pearsonr}[method]
    have = map_values.reindex(structures).dropna()
    rows = []
    for gene, values in sorted(profiles.items()):
        shared = [s for s in have.index if s in values and np.isfinite(values[s])]
        if len(shared) < ISH["min_structures"]:
            continue
        x = have[shared].to_numpy(float)
        y = np.array([values[s] for s in shared])
        rows.append(
            dict(
                symbol=gene,
                rho=float(statistic(x, y).statistic),
                n_structures=len(shared),
            )
        )
    return pd.DataFrame(rows)


def gap_on_shared(
    map_values: pd.Series,
    first: dict[str, float],
    second: dict[str, float],
    structures: list[str],
    method: str = "spearman",
) -> float:
    """rho(first) - rho(second) on the structures both genes and the map have."""
    statistic = {"spearman": spearmanr, "pearson": pearsonr}[method]
    have = map_values.reindex(structures).dropna()
    shared = [s for s in have.index if s in first and s in second]
    x = have[shared].to_numpy(float)
    a = np.array([first[s] for s in shared])
    b = np.array([second[s] for s in shared])
    return float(statistic(x, a).statistic - statistic(x, b).statistic)


def variant_rows(
    name: str,
    map_values: pd.Series,
    profiles: dict[str, dict[str, float]],
    structures: list[str],
    p9_genes: set[str],
    method: str = "spearman",
) -> pd.DataFrame:
    """One variant's rho per gene, with its rank among P9's genes and among all."""
    table = correlate(map_values, profiles, structures, method)
    table.insert(0, "variant", name)
    p9 = table["symbol"].isin(p9_genes)
    table["p9_gene"] = p9
    table["rank_all"] = table["rho"].rank(ascending=False, method="min").astype(int)
    table["rank_p9"] = table.loc[p9, "rho"].rank(ascending=False, method="min")
    return table


def agreement(
    table: pd.DataFrame, primary: pd.DataFrame, genes: set[str] | None
) -> tuple[float, int]:
    """Spearman of the genes' rho in a variant and in the primary, and the genes."""
    a = table.set_index("symbol")["rho"]
    b = primary.set_index("symbol")["rho"]
    both = sorted(set(a.index) & set(b.index))
    if genes is not None:
        both = [g for g in both if g in genes]
    if len(both) < 3:
        return np.nan, len(both)
    return float(spearmanr(a[both], b[both]).statistic), len(both)


def summary_row(
    variant: tuple[str, str, str],
    table: pd.DataFrame,
    primary: pd.DataFrame,
    p9_genes: set[str],
    gap: float,
) -> dict:
    """One row of robustness_summary.csv: agreement, where Cacng8 and Gria1 sit, gap."""
    name, kind, label = variant
    agree_p9, n_p9 = agreement(table, primary, p9_genes)
    agree_all, n_all = agreement(table, primary, None)
    by = table.set_index("symbol")
    row = dict(
        variant=name,
        kind=kind,
        label=label,
        n_genes=len(table),
        n_p9_genes=int(table["p9_gene"].sum()),
        agreement_p9=agree_p9,
        n_agreement_p9=n_p9,
        agreement_all=agree_all,
        n_agreement_all=n_all,
        gap=gap,
    )
    for gene in gene_ranking.GAP_GENES:
        row[f"rho_{gene}"] = by["rho"].get(gene, np.nan)
        row[f"rank_p9_{gene}"] = by["rank_p9"].get(gene, np.nan)
        row[f"rank_all_{gene}"] = by["rank_all"].get(gene, np.nan)
    return row


def variant_tables(
    inputs: dict[str, tuple[pd.Series, dict, list[str], str]], p9_genes: set[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ranking_robustness.csv and robustness_summary.csv, from variant_inputs.

    The primary row must be analysis 1's ranking, and the run stops if it is not.
    The gap is taken where a variant has both genes, NaN where it lacks one.
    """
    parts, gaps = [], {}
    first, second = gene_ranking.GAP_GENES
    for name, _, _ in VARIANTS:
        map_values, gene_profiles, structure_list, method = inputs[name]
        parts.append(
            variant_rows(
                name, map_values, gene_profiles, structure_list, p9_genes, method
            )
        )
        gaps[name] = np.nan
        if first in gene_profiles and second in gene_profiles:
            gaps[name] = gap_on_shared(
                map_values,
                gene_profiles[first],
                gene_profiles[second],
                structure_list,
                method,
            )
    per_gene = pd.concat(parts, ignore_index=True)
    primary = per_gene[per_gene["variant"] == "primary"]
    gene_ranking.check_matches_ranking(
        primary.set_index("symbol")["rho"], "the primary row"
    )
    summary = pd.DataFrame(
        [
            summary_row(
                v, per_gene[per_gene["variant"] == v[0]], primary, p9_genes, gaps[v[0]]
            )
            for v in VARIANTS
        ]
    )
    return per_gene, summary


def numbers_table(summary: pd.DataFrame) -> pd.DataFrame:
    """numbers_robustness.csv: the numbers of this step that the text quotes."""
    rows = []
    for _, r in summary.iterrows():
        name = r["variant"]
        rows += [
            (f"agreement_p9_{name}", round(r["agreement_p9"], 4), r["label"]),
            (f"agreement_all_{name}", round(r["agreement_all"], 4), r["label"]),
            (f"rank_p9_Gria1_{name}", r["rank_p9_Gria1"], r["label"]),
            (f"rank_p9_Cacng8_{name}", r["rank_p9_Cacng8"], r["label"]),
            (f"rho_Gria1_{name}", round(r["rho_Gria1"], 3), r["label"]),
            (f"rho_Cacng8_{name}", round(r["rho_Cacng8"], 3), r["label"]),
            (f"gap_{name}", round(r["gap"], 3), r["label"]),
        ]
    return numbers_frame(rows)


def load_summary() -> pd.DataFrame:
    """The summary run_ish_robustness wrote."""
    if not SUMMARY.exists():
        raise FileNotFoundError(f"{SUMMARY} not found: run run_ish_robustness.py first")
    return pd.read_csv(SUMMARY)


def load_robustness() -> pd.DataFrame:
    """The per-gene table run_ish_robustness wrote."""
    if not ROBUSTNESS.exists():
        raise FileNotFoundError(
            f"{ROBUSTNESS} not found: run run_ish_robustness.py first"
        )
    table = pd.read_csv(ROBUSTNESS)
    table["p9_gene"] = table["p9_gene"].astype(str) == "True"
    return table
