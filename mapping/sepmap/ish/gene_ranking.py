"""Analysis 1: which expression maps order the structures as the nano map does.

Each gene of the gene table is correlated with the adult map on the declared
structures it has (its merged profile, ish.gene_table: structures with at least
ish.min_voxels voxels of data), and the correlation is tested against the spatial
null of ish.spatial_null:

    rho          Spearman over the declared structures the gene has; a gene with
                 fewer than ish.min_structures of them gets no row
    spatial p    two-sided, against the map's surrogates cut to the same
                 structures; the smallest p is 1 / (surrogates + 1)
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
same steps with its own surrogates. Were the ranking the tissue's rather than the
label's, autofluorescence would order the genes as nano does.

The Cacng8 - Gria1 gap (named in advance) is rho(Cacng8) - rho(Gria1) on the
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

from sepmap.adult import profiles
from sepmap.config import SETTINGS
from sepmap.ish import gene_table
from sepmap.ish.gene_sets import LEFTOVER_GENE
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.ish.spatial_null import ALPHA, BAND, null_rho, spatial_p
from sepmap.structures import TABLES

# the voxels a gene value needs and the structures a correlation needs; the BH
# level and the adult resamples
ISH = SETTINGS["ish"]
ISH_ANALYSIS = SETTINGS["ish_analysis"]

RANKING = TABLES / "gene_ranking.csv"
GAP = TABLES / "gap.csv"
PER_ADULT = TABLES / "gene_ranking_per_adult.csv"
NULL_RHO = TABLES / "null_rho.npz"
NUMBERS = numbers_path("gene_ranking")

# the maps, by name: the column of the adult tables that holds each
MAPS = {"nano": "zref_nano", "auto": "zref_auto"}

# the two genes of the gap, named in advance
GAP_GENES = ("Cacng8", "Gria1")

# the row of gap.csv on the genes' merged profiles; the others pair two experiments
MERGED = "merged profiles"

# the genes figure 05 and the text quote: the top of P9's ranking, the gene of the
# stained protein, and an astrocyte gene as the control
QUOTED_GENES = ("Cacng8", "Gria1", "Aqp4")

# surrogates used to set the weight of the shared part of the equal null, and the
# bisection steps that set it
EQUAL_FIT = 1000
EQUAL_STEPS = 40


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


def map_rows(ranking: pd.DataFrame, map_name: str) -> pd.DataFrame:
    """One map's rows of gene_ranking.csv, indexed by gene."""
    return ranking[ranking["map"] == map_name].set_index("symbol")


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


def rank_against_maps(
    per_mouse: pd.DataFrame,
    profile: pd.DataFrame,
    declared: list[str],
    surr: dict[str, np.ndarray],
    vectors: dict[str, tuple[np.ndarray, np.ndarray]],
    genes: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Every gene against each map: rho, spatial p, null band, spread over adults.

    `per_mouse` and `profile` are the adult tables (adult.profiles), `surr` each
    map's surrogates on the `declared` structures, `genes` the gene table's one row
    per gene. Returns gene_ranking.csv's table, each map's null rho (genes x
    surrogates) and each map's adults x structures.
    """
    results, nulls, adults = {}, {}, {}
    for k, (name, column) in enumerate(MAPS.items()):
        adults[name] = adult_matrix(per_mouse, column, declared, profiles.ADULTS)
        boot = bootstrap_maps(adults[name], seed=k)
        map_values = profile.loc[declared, column].to_numpy(float)
        results[name], nulls[name] = rank_genes(map_values, surr[name], boot, vectors)
        print(f"  {name}: {len(results[name])} genes tested", flush=True)
    return ranking_table(results, genes), nulls, adults


def check_matches_ranking(rho: pd.Series, what: str) -> None:
    """Stop unless `rho` (by gene) is the nano map's rho of gene_ranking.csv.

    `what` names the rho in the message: a step that recomputes the whole-brain rho
    must find what analysis 1 wrote, or the two read different inputs.
    """
    nano = map_rows(load_ranking(), "nano")["rho"]
    largest = float((rho - nano).abs().max())
    if largest > 1e-9:
        raise ValueError(
            f"{what} differs from gene_ranking.csv by up to {largest:.2e}; "
            "run run_ish_gene_ranking.py first"
        )


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
    in advance, so the gap passes only beyond equal_two_sided, the 95th percentile of
    the null's gaps taken either way; the equal null is skewed, so beside it the
    share of its maps where the first gene leads by at least the gap, and where the
    second does, and its central 95% (equal_lo, equal_hi), which describe it and are
    not the test.
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
        equal_two_sided=float(np.percentile(np.abs(equal), 100 * (1 - ALPHA))),
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
    rows = [dict(kind=MERGED, first_experiment="", second_experiment="", **row)]
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


def gap_with_nulls(
    table: pd.DataFrame,
    merged: dict[str, dict[str, float]],
    profile: pd.DataFrame,
    declared: list[str],
    surr: dict[str, np.ndarray],
    nano_adults: np.ndarray,
) -> tuple[pd.DataFrame, np.ndarray]:
    """gap_table on the nano map: merged, per pairing of experiments, and its nulls.

    `table` is the gene table, `merged` every gene's merged profile, `nano_adults`
    the nano map's adults x structures, resampled with the seed of the ranking.
    """
    boot = bootstrap_maps(nano_adults, seed=0)
    return gap_table(
        profile.loc[declared, "zref_nano"].to_numpy(float),
        surr["nano"],
        boot,
        merged,
        gene_table.used_experiment_profiles(table),
        declared,
    )


def merged_gap(gap: pd.DataFrame) -> pd.Series:
    """The row of gap.csv on the two genes' merged profiles."""
    return gap[gap["kind"] == MERGED].iloc[0]


# ===== The numbers for the text =====


def numbers_table(
    ranking: pd.DataFrame, gap: pd.DataFrame, per_adult: pd.DataFrame
) -> pd.DataFrame:
    """numbers_gene_ranking.csv: the numbers of this step that the text quotes."""
    q = ISH_ANALYSIS["q"]
    rows = []
    for name in MAPS:
        mine = ranking[ranking["map"] == name]
        p9 = mine[mine["p9_gene"]]
        rows += [
            (f"{name}_genes", len(mine), f"genes correlated with the {name} map"),
            (f"{name}_p9_genes", len(p9), f"P9's genes correlated with the {name} map"),
            (
                f"{name}_pass_p9",
                int((p9["q_p9"] < q).sum()),
                f"P9's genes past the null, BH within P9's genes, q < {q}",
            ),
            (
                f"{name}_pass_all",
                int((mine["q_all"] < q).sum()),
                f"genes past the null, BH within all genes, q < {q}",
            ),
            (
                f"{name}_median_rho",
                round(float(mine["rho"].median()), 3),
                "median rho over all genes",
            ),
        ]
        for gene in QUOTED_GENES:
            r = mine[mine["symbol"] == gene].iloc[0]
            rows += [
                (f"{name}_rho_{gene}", round(r["rho"], 3), f"{gene}'s rho"),
                (f"{name}_p_{gene}", round(r["p_spatial"], 5), f"{gene}'s spatial p"),
                (f"{name}_q_p9_{gene}", round(r["q_p9"], 5), f"{gene}'s q within P9"),
                (f"{name}_rank_p9_{gene}", r["rank_p9"], f"{gene}'s rank among P9's"),
                (f"{name}_rank_all_{gene}", r["rank_all"], f"{gene}'s rank among all"),
            ]
    rows += gap_numbers(gap)
    for gene in (ISH["control_gene"], LEFTOVER_GENE):
        mine = per_adult[per_adult["symbol"] == gene]
        for name in MAPS:
            rows += [
                (
                    f"per_adult_{name}_{gene}_min",
                    round(mine[f"rho_{name}"].min(), 3),
                    f"lowest rho of one adult's {name} map with {gene}",
                ),
                (
                    f"per_adult_{name}_{gene}_max",
                    round(mine[f"rho_{name}"].max(), 3),
                    f"highest rho of one adult's {name} map with {gene}",
                ),
            ]
    return numbers_frame(rows)


def gap_numbers(gap: pd.DataFrame) -> list[tuple]:
    """The numbers of the Cacng8 - Gria1 gap: merged, its two nulls, the pairings."""
    merged = merged_gap(gap)
    pairs = gap[gap["kind"] == "experiment pairing"]
    return [
        ("gap", round(merged["gap"], 3), "rho(Cacng8) - rho(Gria1), merged profiles"),
        ("gap_structures", int(merged["n_structures"]), "structures both genes have"),
        ("gap_p", round(merged["p_equal"], 5), "its p, maps related to both alike"),
        ("gap_equal_lo", round(merged["equal_lo"], 3), "2.5% of that null"),
        ("gap_equal_hi", round(merged["equal_hi"], 3), "97.5% of that null"),
        (
            "gap_two_sided",
            round(merged["equal_two_sided"], 3),
            "the lead either way the test fixed in advance needs: 95% of |null|",
        ),
        (
            "gap_equal_first_as_large",
            round(merged["equal_first_as_large"], 4),
            "share of that null where Cacng8 leads by at least the gap",
        ),
        (
            "gap_equal_second_as_large",
            round(merged["equal_second_as_large"], 4),
            "share of that null where Gria1 leads by at least the gap",
        ),
        ("gap_p_unrelated", round(merged["p_spatial"], 5), "its p, unrelated maps"),
        ("gap_unrelated_lo", round(merged["null_lo"], 3), "2.5% of that null"),
        ("gap_unrelated_hi", round(merged["null_hi"], 3), "97.5% of that null"),
        ("gap_boot_lo", round(merged["boot_lo"], 3), "2.5% over resampled adults"),
        ("gap_boot_hi", round(merged["boot_hi"], 3), "97.5% over resampled adults"),
        ("gap_pairs_min", round(pairs["gap"].min(), 3), "smallest over pairings"),
        ("gap_pairs_max", round(pairs["gap"].max(), 3), "largest over pairings"),
        ("gap_pairs_p_min", round(pairs["p_equal"].min(), 5), "smallest p over pairings"),
        ("gap_pairs_p_max", round(pairs["p_equal"].max(), 5), "largest p over pairings"),
    ]


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


def load_gap() -> pd.DataFrame:
    """The table of the Cacng8 - Gria1 gap that run_ish_gene_ranking wrote."""
    if not GAP.exists():
        raise FileNotFoundError(f"{GAP} not found: run run_ish_gene_ranking.py first")
    return pd.read_csv(GAP)


def save_null_rho(
    genes: list[str], nulls: dict[str, np.ndarray], gap_null: np.ndarray
) -> None:
    """Write null_rho.npz: each map's null rho (genes x surrogates) and the gap's."""
    np.savez(NULL_RHO, genes=np.array(genes), gap=gap_null, **nulls)


def check_null_rho() -> None:
    """Stop unless analysis 1 has written null_rho.npz."""
    if not NULL_RHO.exists():
        raise FileNotFoundError(
            f"{NULL_RHO} not found: run run_ish_gene_ranking.py first"
        )


def load_null_rho(map_name: str) -> tuple[np.ndarray, list[str]]:
    """Every gene's rho with every surrogate of one map, and the genes of its rows."""
    check_null_rho()
    with np.load(NULL_RHO) as z:
        return z[map_name], [str(g) for g in z["genes"]]


def load_gap_null() -> np.ndarray:
    """The Cacng8 - Gria1 gap of every null map, from null_rho.npz.

    Two rows: maps equally related to both genes, then the surrogates of the nano map.
    """
    check_null_rho()
    with np.load(NULL_RHO) as z:
        return z["gap"]
