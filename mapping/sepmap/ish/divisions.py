"""Analysis 2: does a gene follow the map inside divisions, or only between them?

A whole-brain rho mixes two things: the contrast between divisions (cortex and
hippocampus high, thalamus and hypothalamus low) and the order of the structures
inside a division. Brain maps share the first almost whatever they measure, so a
gene can rank high through it alone. Two numbers per gene separate them:

    division-only rho   rho of the gene with a map that knows only each
                        structure's division: every declared structure takes the
                        median of the map over its division's declared structures.
                        A gene whose whole-brain rho is close to this one gets it
                        from the contrast between divisions
    within rho          Spearman inside each division where the gene has at least
                        ish_analysis.min_division_structures declared structures,
                        averaged with the structure counts as weights. Every
                        division then has the same scale, and no contrast between
                        divisions can enter

The within rho is tested two ways:

    shuffle null    the map's values shuffled inside each division,
                    ish_analysis.n_perm_within times: keeps which division a value
                    sits in and destroys the order inside it. It ignores that near
                    structures inside a division are alike, so it is too narrow;
                    the cheap first null of the ISH discussion
                    (docs/history/ISH_DISCUSSION.md), reported beside
    spatial null    the same weighted within rho for every surrogate of the map
                    (ish.spatial_null), whose smoothness inside divisions is the
                    map's; the p the figures use

Both are two-sided and corrected by Benjamini-Hochberg within P9's genes and within
all genes, as in analysis 1. Both nulls are checked on random smooth maps with the
map's smoothness (the spatial null's calibration fields), which have no relation to
the map: a null that keeps its 5% there can be read.

Run by run_ish_divisions.py.
"""

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, rankdata, spearmanr

from sepmap.config import SETTINGS
from sepmap.ish import spatial_null
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.ish.spatial_null import ALPHA, BAND, spatial_p
from sepmap.structures import TABLES

# the structures a division needs to enter, the shuffles, the BH level
ISH_ANALYSIS = SETTINGS["ish_analysis"]

WITHIN = TABLES / "within_division.csv"
WITHIN_DETAIL = TABLES / "within_division_detail.csv"
CALIBRATION = TABLES / "within_calibration.csv"
NUMBERS = numbers_path("divisions")

# the genes of the gene sheets and of figure 09s D: the top of P9's ranking, the
# subunit, a metabotropic receptor and a scaffold near the top, and an astrocyte
# gene as the control
DETAIL_GENES = ("Cacng8", "Gria1", "Grm5", "Dlg2", "Aqp4")

# the gene of panel A of figures 09 and 09s, the top of P9's ranking
EXAMPLE_GENE = "Cacng8"


def division_only(map_values: np.ndarray, divisions: np.ndarray) -> np.ndarray:
    """The map with every structure set to the median of its division."""
    out = np.empty_like(map_values, dtype=float)
    for division in np.unique(divisions):
        mine = divisions == division
        out[mine] = np.median(map_values[mine])
    return out


def division_parts(
    columns: np.ndarray, divisions: np.ndarray, min_structures: int | None = None
) -> list[tuple[str, np.ndarray]]:
    """The divisions a gene enters, and the positions of its values in each.

    `columns` are the gene's structures as indices into the declared list,
    `divisions` the division of every declared structure. Returns (division,
    positions into `columns`) for the divisions with at least
    ish_analysis.min_division_structures of them.
    """
    if min_structures is None:
        min_structures = ISH_ANALYSIS["min_division_structures"]
    mine = divisions[columns]
    out = []
    for division in sorted(set(mine)):
        positions = np.nonzero(mine == division)[0]
        if len(positions) >= min_structures:
            out.append((division, positions))
    return out


def _centred_ranks(values: np.ndarray) -> np.ndarray:
    """Ranks along the last axis, centred and scaled to unit length."""
    r = rankdata(values, axis=-1)
    r = r - r.mean(axis=-1, keepdims=True)
    return r / np.sqrt(np.square(r).sum(axis=-1, keepdims=True))


def within_rho(
    maps: np.ndarray, gene: np.ndarray, parts: list[tuple[str, np.ndarray]]
) -> np.ndarray:
    """The weighted mean of Spearman inside each division, for each row of `maps`.

    `maps` is maps x the gene's structures (one map may be a 1-D array), `gene` the
    gene's values there, `parts` the divisions it enters (division_parts).
    """
    maps = np.atleast_2d(maps)
    total = np.zeros(maps.shape[0])
    weights = 0
    for _, positions in parts:
        g = _centred_ranks(gene[positions])
        m = _centred_ranks(maps[:, positions])
        total += len(positions) * (m @ g)
        weights += len(positions)
    return total / weights


def shuffled_maps(
    map_values: np.ndarray,
    parts: list[tuple[str, np.ndarray]],
    n: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """`n` copies of the map, each with its values shuffled inside every division."""
    maps = np.tile(map_values, (n, 1))
    for _, positions in parts:
        maps[:, positions] = rng.permuted(maps[:, positions], axis=1)
    return maps


def shuffled_within_rho(
    map_values: np.ndarray,
    gene: np.ndarray,
    parts: list[tuple[str, np.ndarray]],
    n: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """The within rho of `n` maps shuffled inside each division."""
    return within_rho(shuffled_maps(map_values, parts, n, rng), gene, parts)


def gene_rows(
    map_values: np.ndarray,
    divisions: np.ndarray,
    surr: np.ndarray,
    vectors: dict[str, tuple[np.ndarray, np.ndarray]],
    seed: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per gene: division-only rho, within rho and its two nulls; per division rho.

    `map_values`, `divisions` and the columns of `surr` follow the declared list;
    `vectors` is ish.gene_ranking.gene_vectors. A gene in no division with enough
    structures gets NaN for the within columns. Returns the per-gene table and the
    per-division one.
    """
    rng = np.random.default_rng(seed)
    n_perm = ISH_ANALYSIS["n_perm_within"]
    coarse = division_only(map_values, divisions)
    rows, detail = [], []
    for gene, (columns, values) in vectors.items():
        parts = division_parts(columns, divisions)
        row = dict(
            symbol=gene,
            n_structures=len(columns),
            rho=float(spearmanr(map_values[columns], values).statistic),
            rho_division_only=float(spearmanr(coarse[columns], values).statistic),
            n_divisions=len(parts),
            n_structures_within=int(sum(len(p) for _, p in parts)),
        )
        if parts:
            observed = float(within_rho(map_values[columns], values, parts)[0])
            shuffled = shuffled_within_rho(
                map_values[columns], values, parts, n_perm, rng
            )
            spatial = within_rho(surr[:, columns], values, parts)
            row.update(
                rho_within=observed,
                p_within_shuffle=spatial_p(observed, shuffled),
                p_within_spatial=spatial_p(observed, spatial),
                shuffle_lo=float(np.percentile(shuffled, BAND[0])),
                shuffle_hi=float(np.percentile(shuffled, BAND[1])),
                null_lo=float(np.percentile(spatial, BAND[0])),
                null_hi=float(np.percentile(spatial, BAND[1])),
            )
            for division, positions in parts:
                x = map_values[columns][positions]
                detail.append(
                    dict(
                        symbol=gene,
                        division=division,
                        n_structures=len(positions),
                        rho=float(spearmanr(x, values[positions]).statistic),
                    )
                )
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(detail)


def with_q(table: pd.DataFrame, p9_genes: set[str]) -> pd.DataFrame:
    """Add BH q of both within nulls, within P9's genes and within all genes."""
    out = table.copy()
    p9 = out["symbol"].isin(p9_genes)
    have = out["p_within_spatial"].notna()
    out.insert(1, "p9_gene", p9)
    for null in ("spatial", "shuffle"):
        column = f"p_within_{null}"
        out[f"q_all_{null}"] = np.nan
        out.loc[have, f"q_all_{null}"] = false_discovery_control(
            out.loc[have, column], method="bh"
        )
        out[f"q_p9_{null}"] = np.nan
        out.loc[have & p9, f"q_p9_{null}"] = false_discovery_control(
            out.loc[have & p9, column], method="bh"
        )
    return out


def agreement(table: pd.DataFrame, column: str, genes: pd.Series | None = None) -> float:
    """Spearman over genes of the whole-brain rho and `column`.

    `genes` marks the genes taken (P9's, say), every gene with a value by default.
    """
    have = table[column].notna()
    if genes is not None:
        have = have & genes
    mine = table[have]
    return float(spearmanr(mine["rho"], mine[column]).statistic)


def calibration(
    map_values: np.ndarray,
    divisions: np.ndarray,
    surr: np.ndarray,
    fields: np.ndarray,
    seed: int = 0,
) -> pd.DataFrame:
    """Random smooth maps tested as genes inside divisions, one row per map.

    `fields` holds maps with no relation to the map (ish.spatial_null.random_fields)
    on every declared structure. Per map: its within rho with the map, and the p
    of the shuffle null (ish_analysis.n_perm_within_calibration shuffles) and of
    the spatial null.
    """
    rng = np.random.default_rng(seed)
    n_perm = ISH_ANALYSIS["n_perm_within_calibration"]
    columns = np.arange(len(map_values))
    parts = division_parts(columns, divisions)
    rows = []
    for k, field in enumerate(fields):
        observed = float(within_rho(map_values, field, parts)[0])
        shuffled = shuffled_within_rho(map_values, field, parts, n_perm, rng)
        spatial = within_rho(surr, field, parts)
        rows.append(
            dict(
                test=k,
                rho_within=observed,
                p_shuffle=spatial_p(observed, shuffled),
                p_spatial=spatial_p(observed, spatial),
            )
        )
    return pd.DataFrame(rows)


def within_calibration(
    map_values: np.ndarray,
    divisions: np.ndarray,
    surr: np.ndarray,
    centroids: pd.DataFrame,
    seed: int = 1,
) -> pd.DataFrame:
    """The calibration on fields with the map's smoothness, as the spatial null's are.

    `centroids` are the declared structures' (sepmap.structures), in the order of
    the map; ish_analysis.n_calibration_within fields of an exponential covariance
    fitted to the map's variogram (ish.spatial_null.fit_exponential), drawn with
    `seed`.
    """
    d = spatial_null.distance_matrix(centroids)
    params = spatial_null.fit_exponential(map_values, d)
    fields = spatial_null.random_fields(
        d, params, ISH_ANALYSIS["n_calibration_within"], np.random.default_rng(seed)
    )
    return calibration(map_values, divisions, surr, fields)


def numbers_table(
    within: pd.DataFrame, detail: pd.DataFrame, calibration_table: pd.DataFrame
) -> pd.DataFrame:
    """numbers_divisions.csv: the numbers of this step that the text quotes."""
    q = ISH_ANALYSIS["q"]
    have = within[within["rho_within"].notna()]
    p9 = within["p9_gene"]
    rows = [
        ("genes", len(within), "genes split into between and within"),
        ("genes_within", len(have), "genes in a division with enough structures"),
        (
            "agreement_division_only_all",
            round(agreement(within, "rho_division_only"), 3),
            "Spearman over genes, whole-brain against division-only rho",
        ),
        (
            "agreement_division_only_p9",
            round(agreement(within, "rho_division_only", p9), 3),
            "the same over P9's genes",
        ),
        (
            "agreement_within_all",
            round(agreement(within, "rho_within"), 3),
            "Spearman over genes, whole-brain against within rho",
        ),
        (
            "agreement_within_p9",
            round(agreement(within, "rho_within", p9), 3),
            "the same over P9's genes",
        ),
        ("median_rho", round(within["rho"].median(), 3), "median whole-brain rho"),
        (
            "median_rho_within",
            round(have["rho_within"].median(), 3),
            "median within rho",
        ),
        (
            "median_rho_p9",
            round(within.loc[p9, "rho"].median(), 3),
            "median whole-brain rho, P9's genes",
        ),
        (
            "median_rho_within_p9",
            round(within.loc[p9, "rho_within"].median(), 3),
            "median within rho, P9's genes",
        ),
    ]
    n_maps = len(calibration_table)
    for null in ("spatial", "shuffle"):
        rows += [
            (
                f"pass_all_{null}",
                int((within[f"q_all_{null}"] < q).sum()),
                f"genes past the within {null} null, BH within all, q < {q}",
            ),
            (
                f"pass_p9_{null}",
                int((within[f"q_p9_{null}"] < q).sum()),
                f"P9's genes past the within {null} null, BH within P9's",
            ),
            (
                f"calibration_{null}",
                round(float((calibration_table[f"p_{null}"] < ALPHA).mean()), 4),
                f"random maps with p < {ALPHA} by the {null} null ({n_maps})",
            ),
        ]
    for gene in DETAIL_GENES:
        r = within[within["symbol"] == gene].iloc[0]
        rows += [
            (f"rho_division_only_{gene}", round(r["rho_division_only"], 3), gene),
            (f"rho_within_{gene}", round(r["rho_within"], 3), gene),
            (f"p_within_spatial_{gene}", round(r["p_within_spatial"], 5), gene),
            (f"p_within_shuffle_{gene}", round(r["p_within_shuffle"], 5), gene),
        ]
    divisions_used = sorted(set(detail["division"]))
    rows.append(("divisions_used", " ".join(divisions_used), "divisions entered"))
    return numbers_frame(rows)


def load_within() -> tuple[pd.DataFrame, pd.DataFrame]:
    """The two tables run_ish_divisions wrote: per gene, and per gene and division."""
    if not WITHIN.exists():
        raise FileNotFoundError(f"{WITHIN} not found: run run_ish_divisions.py first")
    table = pd.read_csv(WITHIN)
    table["p9_gene"] = table["p9_gene"].astype(str) == "True"
    return table, pd.read_csv(WITHIN_DETAIL)
