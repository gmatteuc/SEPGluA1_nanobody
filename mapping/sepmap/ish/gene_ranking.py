"""Analysis 1: which expression maps order the structures as the nano map does.

Each gene of the gene table is correlated with the adult map on the declared
structures it has (its merged profile, ish.gene_table: structures with at least
ish.min_voxels voxels of data), and the correlation is tested against the spatial
null of ish.spatial_null (A7):

    rho          Spearman over the declared structures the gene has; a gene with
                 fewer than ish.min_structures of them gets no row
    spatial p    two-sided, against the map's surrogates cut to the same
                 structures; with 10,000 surrogates the smallest p is 1e-4
    q            Benjamini-Hochberg, within P9's 100 genes (the genes the figures
                 name, as in April) and within every gene of the table
    null band    the 2.5th and 97.5th percentiles of the gene's null rho: the band
                 a rho must leave to pass at 0.05
    spread       the SD of rho over ish_analysis.n_boot_mice resamples of the ten
                 adults with replacement, each resample's map the mean of its
                 adults' zref; rho over that SD (t) is the grey of a bar, how far the
                 rho depends on which adults were measured

The autofluorescence map of the same brains (each brain over its own isocortex
mean of autofluorescence, scaled over the declared set as zref is) goes through the
same steps with its own surrogates (A8). Were the ranking the tissue's rather than
the label's, autofluorescence would order the genes as nano does.

The Cacng8 - Gria1 gap (named in advance, S5) is rho(Cacng8) - rho(Gria1) on the
structures both genes have. The question is whether the map is more closely related
to one gene than to the other, so its null is a map related to both alike: each
null map is

    c z(Cacng8 + Gria1) + sqrt(1 - c^2) z(surrogate)

the two genes' standardised ranks added, plus a surrogate of the nano map for the
rest, with c set so that the null maps' mean rho with the two genes equals the
observed mean of the two rhos. The gap of every null map gives the distribution a
gap takes when the map follows both genes equally, with the brain's smoothness;
correlations near 0.7 vary less than correlations near 0, so this null is narrower
than that of a map unrelated to both. That second null (the surrogates themselves,
the two genes' own correlation kept) is kept beside it as a conservative bound.
Beside them, the gap's interval over the adult bootstrap, and the gap for each
pairing of the two genes' Allen experiments, each one Allen mouse.

Every gene's rho with every surrogate is kept (null_rho.npz), so the set tests of
analysis 3 read the same null as the genes.

Run by run_ish_gene_ranking.py.
"""

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, rankdata, spearmanr

from sepmap.config import SETTINGS
from sepmap.ish.spatial_null import null_rho, spatial_p
from sepmap.structures import TABLES

# the voxels a gene value needs and the structures a correlation needs; the BH
# level and the adult resamples
ISH = SETTINGS["ish"]
ISH_ANALYSIS = SETTINGS["ish_analysis"]

RANKING = TABLES / "gene_ranking.csv"
GAP = TABLES / "gap.csv"
PER_ADULT = TABLES / "gene_ranking_per_adult.csv"
NULL_RHO = TABLES / "null_rho.npz"

# the maps, by name: the column of the adult tables that holds each
MAPS = {"nano": "zref_nano", "auto": "zref_auto"}

# the two genes of the gap, named in advance (S5)
GAP_GENES = ("Cacng8", "Gria1")

# the percentiles of a null distribution that a rho must leave to pass at 0.05
BAND = (2.5, 97.5)


# ===== Inputs =====


def gene_vectors(
    profiles: dict[str, dict[str, float]], structures: list[str]
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """{gene: (indices into `structures`, the gene's values there)}.

    A gene with fewer than ish.min_structures of the structures is left out.
    """
    position = {s: i for i, s in enumerate(structures)}
    out = {}
    for gene, values in sorted(profiles.items()):
        have = [s for s in structures if s in values and np.isfinite(values[s])]
        if len(have) < ISH["min_structures"]:
            continue
        columns = np.array([position[s] for s in have])
        out[gene] = (columns, np.array([values[s] for s in have], dtype=float))
    return out


def adult_matrix(
    per_mouse: pd.DataFrame, column: str, structures: list[str], adults: list[str]
) -> np.ndarray:
    """Adults x structures of one column of the per-adult table, in the given order.

    Every declared structure is measured in every adult, so a missing value means
    the table and the set disagree, and the run stops.
    """
    table = per_mouse.pivot(index="mouse", columns="structure", values=column)
    table = table.reindex(index=adults, columns=structures)
    if table.isna().any().any():
        missing = int(table.isna().sum().sum())
        raise ValueError(
            f"{missing} values of {column} missing for the declared structures; the "
            "per-adult table and the set disagree. Run run_structure_set.py again."
        )
    return table.to_numpy(float)


def bootstrap_maps(adults: np.ndarray, n: int | None = None, seed: int = 0) -> np.ndarray:
    """Maps of resampled cohorts: each the mean of the adults drawn with replacement.

    `adults` is adults x structures; returns n x structures.
    """
    if n is None:
        n = ISH_ANALYSIS["n_boot_mice"]
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, adults.shape[0], size=(n, adults.shape[0]))
    return adults[draws].mean(axis=1)


# ===== Gene by gene =====


def rank_genes(
    map_values: np.ndarray,
    surr: np.ndarray,
    boot: np.ndarray,
    vectors: dict[str, tuple[np.ndarray, np.ndarray]],
) -> tuple[pd.DataFrame, np.ndarray]:
    """rho, spatial p, null band and adult spread of every gene against one map.

    `map_values` holds the map on the declared structures, `surr` its surrogates
    and `boot` its adult resamples, both with one structure per column in the same
    order. Returns one row per gene and the null rho, genes x surrogates.
    """
    rows, nulls = [], []
    for gene, (columns, values) in vectors.items():
        test = spearmanr(map_values[columns], values)
        null = null_rho(surr[:, columns], values)
        spread = null_rho(boot[:, columns], values)
        sd = float(spread.std(ddof=1))
        lo, hi = np.percentile(null, BAND)
        rows.append(
            dict(
                symbol=gene,
                n_structures=len(columns),
                rho=float(test.statistic),
                p_ordinary=float(test.pvalue),
                p_spatial=spatial_p(test.statistic, null),
                null_lo=float(lo),
                null_hi=float(hi),
                boot_sd=sd,
                boot_lo=float(np.percentile(spread, BAND[0])),
                boot_hi=float(np.percentile(spread, BAND[1])),
                t_boot=float(abs(test.statistic) / sd) if sd > 0 else np.nan,
            )
        )
        nulls.append(null.astype(np.float32))
    return pd.DataFrame(rows), np.vstack(nulls)


def with_q_and_ranks(table: pd.DataFrame, p9_genes: set[str]) -> pd.DataFrame:
    """Add BH q values and ranks, within P9's genes and within all genes.

    Ranks count from the highest rho (1); a gene outside P9's panel has no q_p9 and
    no rank_p9.
    """
    out = table.copy()
    p9 = out["symbol"].isin(p9_genes)
    out["q_all"] = false_discovery_control(out["p_spatial"], method="bh")
    out["q_p9"] = np.nan
    out.loc[p9, "q_p9"] = false_discovery_control(out.loc[p9, "p_spatial"], method="bh")
    out["rank_all"] = out["rho"].rank(ascending=False, method="min").astype(int)
    out["rank_p9"] = out.loc[p9, "rho"].rank(ascending=False, method="min")
    return out


def ranking_table(tables: dict[str, pd.DataFrame], genes: pd.DataFrame) -> pd.DataFrame:
    """gene_ranking.csv: one row per map and gene, with q, ranks and the gene's labels.

    `tables` holds rank_genes' table of each map, `genes` one row per gene
    (ish.gene_table.per_gene).
    """
    p9_genes = set(genes.loc[genes["p9_gene"], "symbol"])
    labels = genes.set_index("symbol")[
        ["p9_gene", "p9_category", "gene_sets", "reliability", "n_experiments_used"]
    ]
    parts = []
    for name, table in tables.items():
        part = with_q_and_ranks(table, p9_genes)
        part.insert(0, "map", name)
        parts.append(part.join(labels, on="symbol"))
    return pd.concat(parts, ignore_index=True)


def passing(ranking: pd.DataFrame, map_name: str, within: str) -> pd.DataFrame:
    """The genes past the null at ish_analysis.q for one map, BH within p9 or all."""
    mine = ranking[ranking["map"] == map_name]
    return mine[mine[f"q_{within}"] < ISH_ANALYSIS["q"]]


def per_adult_table(
    maps: dict[str, np.ndarray],
    vectors: dict[str, tuple[np.ndarray, np.ndarray]],
    adults: list[str],
    group: dict[str, str],
) -> pd.DataFrame:
    """rho of every gene with each adult's own map, one row per adult and gene.

    `maps` holds adults x structures for each map name; columns rho_<map>.
    """
    rows = []
    for i, mouse in enumerate(adults):
        for gene, (columns, values) in vectors.items():
            row = dict(mouse=mouse, group=group[mouse], symbol=gene)
            for name, matrix in maps.items():
                row[f"rho_{name}"] = float(
                    spearmanr(matrix[i, columns], values).statistic
                )
            rows.append(row)
    return pd.DataFrame(rows)


# ===== The Cacng8 - Gria1 gap =====


def standardise_rows(x: np.ndarray) -> np.ndarray:
    """Ranks along the last axis, centred and scaled to unit SD."""
    r = rankdata(x, axis=-1).astype(float)
    r = r - r.mean(axis=-1, keepdims=True)
    return r / r.std(axis=-1, keepdims=True)


# surrogates used to set the weight of the shared part of the equal null, and the
# bisection steps that set it
EQUAL_FIT = 1000
EQUAL_STEPS = 40


def equal_null(
    first: np.ndarray, second: np.ndarray, surr: np.ndarray, target: float
) -> tuple[np.ndarray, float]:
    """Null gaps of maps equally related to both genes, and the weight c used.

    Each null map is c z(z(first) + z(second)) + sqrt(1 - c^2) z(surrogate), on
    ranks; c is found by bisection so that the mean, over the first EQUAL_FIT
    surrogates, of the null maps' average rho with the two genes equals `target`.
    """
    shared = standardise_rows(standardise_rows(first) + standardise_rows(second))
    noise = standardise_rows(surr)

    def mean_rho(c: float, rows: np.ndarray) -> np.ndarray:
        """The null maps' mean rho with the two genes, at weight c."""
        maps = c * shared[None, :] + np.sqrt(1 - c**2) * rows
        return (null_rho(maps, first) + null_rho(maps, second)) / 2

    lo, hi = 0.0, 1.0
    for _ in range(EQUAL_STEPS):
        c = (lo + hi) / 2
        if mean_rho(c, noise[:EQUAL_FIT]).mean() < target:
            lo = c
        else:
            hi = c
    c = (lo + hi) / 2
    maps = c * shared[None, :] + np.sqrt(1 - c**2) * noise
    return null_rho(maps, first) - null_rho(maps, second), c


def gap_row(
    map_values: np.ndarray,
    surr: np.ndarray,
    boot: np.ndarray,
    first: dict[str, float],
    second: dict[str, float],
    structures: list[str],
) -> tuple[dict, np.ndarray, np.ndarray]:
    """The gap rho(first) - rho(second) on the structures both have, and its nulls.

    Returns a row (n_structures, rho of each, gap, its p against maps equally
    related to both genes and against maps unrelated to both, both bands, the adult
    bootstrap interval), the null gaps of the equal null and those of the unrelated
    one, one per surrogate. The p counts gaps as large either way, the test fixed
    in advance; the equal null is skewed, so beside it the share of its maps where
    the first gene leads by at least the gap, and where the second does.
    """
    shared = [
        i
        for i, s in enumerate(structures)
        if np.isfinite(first.get(s, np.nan)) and np.isfinite(second.get(s, np.nan))
    ]
    columns = np.array(shared)
    a = np.array([first[structures[i]] for i in shared])
    b = np.array([second[structures[i]] for i in shared])
    rho_a = float(spearmanr(map_values[columns], a).statistic)
    rho_b = float(spearmanr(map_values[columns], b).statistic)
    unrelated = null_rho(surr[:, columns], a) - null_rho(surr[:, columns], b)
    equal, weight = equal_null(a, b, surr[:, columns], (rho_a + rho_b) / 2)
    boot_gap = null_rho(boot[:, columns], a) - null_rho(boot[:, columns], b)
    gap = rho_a - rho_b
    row = dict(
        n_structures=len(shared),
        rho_first=rho_a,
        rho_second=rho_b,
        gap=gap,
        p_equal=spatial_p(gap, equal),
        equal_lo=float(np.percentile(equal, BAND[0])),
        equal_hi=float(np.percentile(equal, BAND[1])),
        equal_weight=weight,
        equal_first_as_large=float(np.mean(equal >= abs(gap))),
        equal_second_as_large=float(np.mean(equal <= -abs(gap))),
        p_spatial=spatial_p(gap, unrelated),
        null_lo=float(np.percentile(unrelated, BAND[0])),
        null_hi=float(np.percentile(unrelated, BAND[1])),
        boot_lo=float(np.percentile(boot_gap, BAND[0])),
        boot_hi=float(np.percentile(boot_gap, BAND[1])),
    )
    return row, equal, unrelated


def gap_table(
    map_values: np.ndarray,
    surr: np.ndarray,
    boot: np.ndarray,
    merged: dict[str, dict[str, float]],
    per_experiment: dict[str, dict[str, dict[str, float]]],
    structures: list[str],
    genes: tuple[str, str] = GAP_GENES,
) -> tuple[pd.DataFrame, np.ndarray]:
    """gap.csv: the gap on the merged profiles, then for each pairing of experiments.

    `per_experiment` is {gene: {experiment: {structure: energy}}}
    (ish.gene_table.experiment_profiles). Returns the table and the null gaps of
    the merged profiles, equal null then unrelated null, 2 x surrogates.
    """
    first, second = genes
    row, equal, unrelated = gap_row(
        map_values, surr, boot, merged[first], merged[second], structures
    )
    rows = [
        dict(kind="merged profiles", first_experiment="", second_experiment="", **row)
    ]
    for ea, pa in sorted(per_experiment[first].items()):
        for eb, pb in sorted(per_experiment[second].items()):
            pair, _, _ = gap_row(map_values, surr, boot, pa, pb, structures)
            rows.append(
                dict(
                    kind="experiment pairing",
                    first_experiment=ea,
                    second_experiment=eb,
                    **pair,
                )
            )
    out = pd.DataFrame(rows)
    out.insert(1, "first", first)
    out.insert(2, "second", second)
    return out, np.vstack([equal, unrelated])


# ===== Reading back =====


def load_ranking() -> pd.DataFrame:
    """The table run_ish_gene_ranking wrote, p9_gene as a boolean."""
    if not RANKING.exists():
        raise FileNotFoundError(f"{RANKING} not found: run run_ish_gene_ranking.py first")
    table = pd.read_csv(RANKING, keep_default_na=False, na_values=[""])
    table["p9_gene"] = table["p9_gene"].astype(str) == "True"
    for column in ("p9_category", "gene_sets"):
        table[column] = table[column].fillna("")
    return table


def load_null_rho(map_name: str) -> tuple[np.ndarray, list[str]]:
    """Every gene's rho with every surrogate of one map, and the genes of its rows.

    `map_name` gap gives the Cacng8 - Gria1 gap of every null map, two rows: maps
    equally related to both genes, then the surrogates of the nano map.
    """
    if not NULL_RHO.exists():
        raise FileNotFoundError(
            f"{NULL_RHO} not found: run run_ish_gene_ranking.py first"
        )
    with np.load(NULL_RHO) as z:
        return z[map_name], [str(g) for g in z["genes"]]
