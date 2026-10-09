"""The genes that follow the map most, characterised; the tests named for the leftover.

Part 2 of the argument (docs/ISH_ANALYSIS.md): the map is more than Gria1 expression
and synapse density, and the reading the data support is the surface fraction of the
receptor. The genes that follow the map most are its corroboration, so each is
described on every axis the ISH line has, and a reader sees what kind of map it is:

    the map            rho with the nano map, its spatial p, q over every gene and
                       rank (analysis 1, ish.gene_ranking)
    inside divisions   the mean rho inside divisions and its spatial p (analysis 2,
                       ish.divisions)
    the choices        the range of its rho over the robustness variants (A3,
                       ish.robustness), and of its rank over the variants that hold
                       every gene (a variant of P9's genes alone ranks among fewer)
    its Allen map      how many experiments, and how well they agree (reliability)
    Gria1 or density   rho with Gria1, each density term, autofluorescence and the
                       main model's prediction on the structures of the fit, and
                       with the PSD95 punctum density where it is measured: whether
                       the gene is just a Gria1-like or density-like map
    the leftover       rho with the main model's leftover and its spatial p, each
                       surrogate put through the same fit (leftover_genes.csv of
                       adult.beyond_density)
    what it takes      the main model with the gene's ranks added (x, x^2, x^3),
                       scored on held-out structures as the main model is: the share
                       of the reproducible map it takes from the leftover, against
                       surrogate maps with the gene's own smoothness added instead
    annotation         its name, the GO terms that put it in the ontology panel, its
                       GO cellular-component terms at the synapse, its gene sets and
                       P9's category (a label), all from the gene table's records

The genes: those past the spatial null after BH over every gene in analysis 1,
Gria1 and Cacng8 always, and every member of the AMPA receptor complex family that
the leftover table holds, whether it passes or not.

The tests against the leftover were named on 8 October 2026, before the main model
ran (ish.gene_sets holds the family, its sources and what was seen before):

    1  Cacng8      its rho with the leftover; the uncorrected spatial p is the test
    2  the family  its median rho with the leftover against the same genes' median
                   over the leftover's surrogates (ish.gene_sets.set_test); and its
                   rho against as many genes of the other postsynaptic set outside
                   the family, each the unused one closest in log median energy
                   (ish.panel_test.greedy_match), labels permuted
                   (ish.panel_test.two_sample), one test each; then gene by gene,
                   BH within the family
    3  the rest    every other gene, exploratory, BH over every gene (q_all of
                   leftover_genes.csv)

The same three are read on the nano map itself, to describe the family, not as tests
of the leftover.

Why the share taken has a null of its own: three more columns always win a little
held-out variance by chance, and a smooth map wins more than a rough one, since
its neighbours carry the same pattern into the test folds. So a gene's gain is set
beside what maps of the gene's own smoothness gain in its place; p is the share of
them that take at least as much (one-sided: only taking more is evidence), the
observed gain counted once. Two kinds of map, as the Cacng8 - Gria1 gap has two
nulls (ish.gene_ranking):

    plain   the gene's values in a new order with its variogram (ish.spatial_null).
            Such a map is unrelated to the main model, so more of it is new to the
            model than of a gene that shares the model's pattern, as the genes
            that follow the map do: a wide null, the conservative bound
    alike   the gene's fit on the main model's columns kept, its remainder replaced
            by a surrogate of the remainder as large: maps that relate to the model
            as the gene does, whose part beyond it is a smooth map unrelated to the
            leftover. Added on 9 October, after the plain null's numbers were seen

The main model is refitted on the structures of the fit the gene has, every
predictor and the ceiling rebuilt there as the jackknife of adult.beyond_figures
does, and the gene, the main model and every map meet the same folds.

Run by run_ish_top_genes.py.
"""

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, rankdata, spearmanr

from sepmap.adult import beyond_density as bd
from sepmap.adult import synaptome
from sepmap.config import SETTINGS
from sepmap.ish import gene_table, spatial_null
from sepmap.ish.gene_sets import (
    AMPA_FAMILY,
    GO_AMPA_COMPLEX,
    GO_AMPA_COMPLEX_RELEASE,
    LEFTOVER_GENE,
    PARTNER_SUBUNITS,
    SCHWENK_2012,
    set_test,
)
from sepmap.ish.panel_test import greedy_match, two_sample
from sepmap.structures import TABLES

# the structures a gene needs, the BH level, the surrogates added in place of a gene
ISH = SETTINGS["ish"]
ISH_ANALYSIS = SETTINGS["ish_analysis"]
TOP_GENES = SETTINGS["top_genes"]

TOP_TABLE = TABLES / "top_genes.csv"
NAMED_TESTS = TABLES / "named_tests.csv"
FAMILY_TABLE = TABLES / "family_members.csv"
NUMBERS = TABLES / "numbers_top_genes.csv"

# characterised whatever their rank: the gene of the stained protein, and the gene
# named for the leftover
NAMED_GENES = ("Gria1", LEFTOVER_GENE)

# where each member of the family comes from, as top_genes.csv writes it
FAMILY_SOURCES = (
    ("Schwenk 2012", tuple(SCHWENK_2012)),
    (f"GO:0032281 ({GO_AMPA_COMPLEX_RELEASE})", GO_AMPA_COMPLEX),
    ("partner subunit", PARTNER_SUBUNITS),
)

# the cellular component under which a gene's GO terms count as synaptic
SYNAPSE = "GO:0045202"

# the percentiles of a null distribution that a value must leave to pass at 0.05
BAND = (2.5, 97.5)


# ===== The genes =====


def chosen_genes(nano: pd.DataFrame, leftover: pd.DataFrame, q: float) -> pd.DataFrame:
    """The genes characterised, best with the map first, and why each is here.

    `nano` is the nano map's rows of gene_ranking.csv, `leftover` leftover_genes.csv.
    top: past the null after BH over every gene; named: NAMED_GENES; family: a
    member of the AMPA receptor complex family that the leftover table holds. tier
    is the tier of the tests against the leftover (ish.gene_sets) the gene is read
    in.
    """
    top = set(nano.loc[nano["q_all"] < q, "symbol"])
    family = set(AMPA_FAMILY) & set(leftover["symbol"])
    rows = []
    for symbol in sorted(top | set(NAMED_GENES) | family):
        if symbol == LEFTOVER_GENE:
            tier = "1 named in advance"
        elif symbol in family:
            tier = "2 the family"
        else:
            tier = "3 exploratory"
        rows.append(
            dict(
                symbol=symbol,
                tier=tier,
                top=symbol in top,
                named=symbol in NAMED_GENES,
                family=symbol in family,
            )
        )
    out = pd.DataFrame(rows)
    rho = nano.set_index("symbol")["rho"]
    order = out["symbol"].map(rho).sort_values(ascending=False).index
    return out.loc[order].reset_index(drop=True)


def family_members(genes: pd.DataFrame, leftover: pd.DataFrame) -> pd.DataFrame:
    """family_members.csv: every member of the family, its sources, tested or why not.

    `genes` is the gene table, one row per gene (ish.gene_table.per_gene); a member
    is tested when the leftover table holds it.
    """
    listed = set(genes["symbol"])
    usable = set(genes.loc[genes["n_experiments_used"] > 0, "symbol"])
    tested = set(leftover["symbol"])
    rows = []
    for symbol in AMPA_FAMILY:
        sources = [name for name, members in FAMILY_SOURCES if symbol in members]
        if symbol in tested:
            reason = ""
        elif symbol not in listed:
            reason = "not in the gene table: in neither panel"
        elif symbol not in usable:
            reason = "no usable Allen experiment"
        else:
            reason = f"fewer than {ISH['min_structures']} structures of the fit"
        rows.append(
            dict(
                symbol=symbol,
                protein=SCHWENK_2012.get(symbol, ""),
                sources="; ".join(sources),
                tested=reason == "",
                reason=reason,
            )
        )
    return pd.DataFrame(rows)


# ===== What kind of map each gene is =====


def named_terms(terms: list[str], names: dict[str, str]) -> str:
    """GO ids with their names, joined: 'GO:0014069 postsynaptic density; ...'."""
    return "; ".join(f"{t} {names.get(t, '')}".strip() for t in terms)


def annotation_table(
    symbols: list[str],
    documentation: pd.DataFrame,
    records: dict[str, dict],
    parents: dict[str, set[str]],
    names: dict[str, str],
) -> pd.DataFrame:
    """Per gene, what the gene table says it is.

    `documentation` is gene_documentation.csv, `records` the genes' mygene records,
    `parents` and `names` GO's relations and term names (ish.gene_table.load_obo).
    go_panel_terms are the GO terms that put the gene in the ontology panel;
    go_synapse_terms its direct cellular-component annotations (NOT left out) at
    or below GO:0045202 synapse, which carry SynGO's curation where it reached GO.
    """
    doc = documentation.set_index("symbol")
    rows = []
    for symbol in symbols:
        panel = doc["ontology_go_terms"].get(symbol, "")
        panel = panel.split() if isinstance(panel, str) else []
        direct = gene_table.go_annotations(records.get(symbol, {}), "CC")
        synaptic = sorted(
            t for t in direct if SYNAPSE in gene_table.with_ancestors({t}, parents)
        )
        rows.append(
            dict(
                symbol=symbol,
                name=doc["name"].get(symbol, ""),
                ontology_role=doc["ontology_role"].get(symbol, ""),
                go_panel_terms=named_terms(panel, names),
                go_synapse_terms=named_terms(synaptic, names),
                gene_sets=doc["gene_sets"].get(symbol, ""),
                p9_category=doc["p9_category"].get(symbol, ""),
            )
        )
    return pd.DataFrame(rows).fillna("")


def robustness_ranges(robust: pd.DataFrame, symbols: list[str]) -> pd.DataFrame:
    """Per gene, its rho and rank over the robustness variants, the primary left out.

    `robust` is ranking_robustness.csv. rho over every variant that has the gene;
    rank over the variants that hold every gene of the primary.
    """
    others = robust[robust["variant"] != "primary"]
    n_primary = int((robust["variant"] == "primary").sum())
    size = others.groupby("variant")["symbol"].size()
    full = set(size.index[size == n_primary])
    rows = []
    for symbol in symbols:
        mine = others[others["symbol"] == symbol]
        ranked = mine[mine["variant"].isin(full)]
        rows.append(
            dict(
                symbol=symbol,
                n_variants=len(mine),
                rho_variants_min=float(mine["rho"].min()),
                rho_variants_max=float(mine["rho"].max()),
                rank_variants_min=int(ranked["rank_all"].min()),
                rank_variants_max=int(ranked["rank_all"].max()),
            )
        )
    return pd.DataFrame(rows)


def gene_on(inputs: bd.Inputs, symbol: str) -> tuple[list[int], np.ndarray]:
    """The positions of the structures of the fit the gene has, and its values there."""
    profile = inputs.expr.get(symbol, {})
    columns = [
        i for i, s in enumerate(inputs.structures) if np.isfinite(profile.get(s, np.nan))
    ]
    values = np.array([profile[inputs.structures[i]] for i in columns], dtype=float)
    return columns, values


def term_names(terms: dict[str, tuple[str, ...]]) -> list[str]:
    """The main model's predictors by name, in the order of its budget."""
    return [name for group in bd.ORDER for name in terms[group]]


def predictor_rhos(inputs: bd.Inputs, symbols: list[str]) -> pd.DataFrame:
    """Per gene, rho with each term of the main model, its prediction, and PSD95.

    On the structures of the fit the gene has (n_fit): rho with Gria1, each density
    term, autofluorescence and the main model's in-sample prediction; and on those
    of them where the PSD95 punctum density is measured (n_psd95), rho with it.
    """
    covariates, _, _ = bd.covariates_for(inputs)
    y = bd.full_map(inputs.nano)
    prediction = y - bd.residual(y, bd.model(covariates, inputs.terms))
    psd95 = inputs.synapses.reindex(columns=[synaptome.MEASURE])[synaptome.MEASURE]
    psd95 = psd95.to_numpy(float)
    rows = []
    for symbol in symbols:
        columns, values = gene_on(inputs, symbol)
        row = dict(symbol=symbol, n_fit=len(columns))
        for name in term_names(inputs.terms):
            row[f"rho_{name}"] = float(
                spearmanr(values, covariates[name][columns]).statistic
            )
        row["rho_prediction"] = float(spearmanr(values, prediction[columns]).statistic)
        measured = np.isfinite(psd95[columns])
        row["n_psd95"] = int(measured.sum())
        row["rho_psd95"] = float(
            spearmanr(values[measured], psd95[columns][measured]).statistic
        )
        rows.append(row)
    return pd.DataFrame(rows)


def remainder_surrogates(
    gene: np.ndarray, design: np.ndarray, d: np.ndarray, n: int, seed: int
) -> np.ndarray:
    """Maps that relate to the model as the gene does, with a surrogate remainder.

    The gene's least-squares fit on `design` (the model's columns and a constant)
    is kept, and its remainder replaced by a surrogate of it (ish.spatial_null),
    made orthogonal to the design and scaled to the remainder's SD; maps x
    structures.
    """
    hat = design @ np.linalg.pinv(design)
    fitted = hat @ gene
    rest = gene - fitted
    surr = spatial_null.surrogates(rest, d, n, seed=seed).astype(float)
    surr = surr - surr @ hat.T
    surr = surr / surr.std(axis=1, keepdims=True) * rest.std()
    return fitted[None, :] + surr


def added_share(
    inputs: bd.Inputs,
    symbol: str,
    centroids: pd.DataFrame,
    n_surrogates: int | None = None,
    seed: int = 0,
) -> tuple[dict, dict[str, np.ndarray]]:
    """What the gene takes of the leftover when added to the main model, and its nulls.

    On the structures of the fit the gene has, the main model and the main model
    with the gene's ranks bent (x, x^2, x^3) are scored on the same held-out folds
    as a share of the ceiling there; taken is the main model's share left minus the
    share left with the gene. Two nulls, top_genes.n_added_surrogates maps each,
    put in the gene's place and ranked as it is:

        plain   surrogates of the gene's map. They are unrelated to the main model,
                so more of each is new to it than of a gene that shares the model's
                pattern: a wide null, the conservative bound
        alike   maps that relate to the model as the gene does, its remainder a
                surrogate (remainder_surrogates): the null of a gene whose part
                beyond the model is a smooth map unrelated to the leftover

    `centroids` holds each structure's one-hemisphere centroid (ap_mm, dv_mm,
    ml_mm; sepmap.structures.load_centroids). Returns the row (n_added, left_main,
    left_with_gene, taken, and per null its median, 95th percentile and one-sided
    p) and the nulls' values by name.
    """
    if n_surrogates is None:
        n_surrogates = TOP_GENES["n_added_surrogates"]
    columns, values = gene_on(inputs, symbol)
    sub = inputs
    if len(columns) < len(inputs.structures):
        sub = bd.restrict(inputs, [inputs.structures[i] for i in columns])
    _, explainable = bd.ceiling(sub.nano, bd.half_splits())
    y = bd.full_map(sub.nano)
    covariates, _, _ = bd.covariates_for(sub)
    xs = bd.model(covariates, sub.terms)
    labels = bd.fold_labels(len(y))

    def left_with(extra: np.ndarray | None) -> float:
        """The share of the reproducible map left, with `extra` bent and added."""
        used = xs if extra is None else xs + bd.flexible([rankdata(extra)])
        return 1 - bd.cv_r2(y, used, labels) / explainable

    left_main = left_with(None)
    left_gene = left_with(values)
    taken = left_main - left_gene
    gene = rankdata(values)
    d = spatial_null.distance_matrix(centroids.loc[sub.structures])
    design = np.column_stack(xs + [np.ones(len(y))])
    maps = dict(
        plain=spatial_null.surrogates(gene, d, n_surrogates, seed=seed),
        alike=remainder_surrogates(gene, design, d, n_surrogates, seed + 1),
    )
    row = dict(
        n_added=len(y),
        left_main=left_main,
        left_with_gene=left_gene,
        taken=taken,
    )
    nulls = {}
    for name, surr in maps.items():
        null = np.array([left_main - left_with(s) for s in surr])
        suffix = "" if name == "plain" else f"_{name}"
        row[f"taken_null{suffix}_median"] = float(np.median(null))
        row[f"taken_null{suffix}_hi"] = float(np.percentile(null, 95))
        row[f"p_taken{suffix}"] = float((np.sum(null >= taken) + 1) / (len(null) + 1))
        nulls[name] = null
    return row, nulls


# ===== The tests named for the leftover =====


def family_controls(
    family: list[str], pool: list[str], level: dict[str, float]
) -> dict[str, str]:
    """{family gene: its control}: as many genes of the pool, matched on expression.

    The localisation test's matching (ish.panel_test.greedy_match): the family's
    genes in order of median energy, highest first, each taking the unused pool
    gene closest in log median energy. The pool is the other postsynaptic set
    outside the family.
    """
    pool = [g for g in pool if g not in family]
    return greedy_match(family, pool, level)


def tier_one(leftover: pd.DataFrame) -> dict:
    """Tier 1: Cacng8's rho with the leftover, its null band and its spatial p."""
    r = leftover.set_index("symbol").loc[LEFTOVER_GENE]
    return dict(
        tier="1",
        map="leftover",
        test=f"{LEFTOVER_GENE} against the leftover's surrogates",
        role="the test",
        n_first=1,
        n_second=np.nan,
        first=float(r["rho"]),
        second=np.nan,
        difference=np.nan,
        null_lo=float(r["null_lo"]),
        null_hi=float(r["null_hi"]),
        critical=np.nan,
        p=float(r["p_spatial"]),
    )


def group_tests(
    family: list[str],
    pairs: dict[str, str],
    maps: dict[str, tuple[pd.Series, np.ndarray, list[str]]],
    seed: int = 0,
) -> tuple[pd.DataFrame, dict[str, dict[str, np.ndarray]]]:
    """Tier 2 as a group, on each map: against the surrogates and against controls.

    `maps` holds per map name (leftover, nano) every gene's rho (by symbol), its
    null rho (genes x surrogates) and the genes of the null's rows. Per map two
    rows: the family's median rho against the median of the same genes over the
    surrogates (second, the null's median), and the median difference from the
    matched controls (`pairs`) with its label-permutation p, `critical` the
    difference beyond which that p falls below 0.05. Each label test has a
    generator of its own. On the leftover the rows are the test; on the nano map a
    description. Returns the rows and, per map, the surrogates' medians (spatial)
    and the label null (labels).
    """
    controls = [pairs[g] for g in family if g in pairs]
    seeds = np.random.SeedSequence(seed).spawn(len(maps))
    rows, nulls = [], {}
    for (name, (rho, null, null_genes)), child in zip(maps.items(), seeds):
        role = "the test" if name == "leftover" else "description"
        row_of = {g: i for i, g in enumerate(null_genes)}
        spatial, null_median = set_test(family, rho, null, row_of)
        rows.append(
            dict(
                tier="2",
                map=name,
                test="the family's median against the surrogates",
                role=role,
                n_first=len(family),
                n_second=np.nan,
                first=spatial["median_rho"],
                second=float(np.median(null_median)),
                difference=spatial["median_rho"] - float(np.median(null_median)),
                null_lo=spatial["null_lo"],
                null_hi=spatial["null_hi"],
                critical=np.nan,
                p=spatial["p_spatial"],
            )
        )
        a = rho[family].to_numpy(float)
        b = rho[controls].to_numpy(float)
        difference, p, label_null = two_sample(a, b, np.random.default_rng(child))
        rows.append(
            dict(
                tier="2",
                map=name,
                test="the family against matched postsynaptic controls",
                role=role,
                n_first=len(a),
                n_second=len(b),
                first=float(np.median(a)),
                second=float(np.median(b)),
                difference=difference,
                null_lo=float(np.percentile(label_null, BAND[0])),
                null_hi=float(np.percentile(label_null, BAND[1])),
                critical=float(np.percentile(np.abs(label_null), 95)),
                p=p,
            )
        )
        nulls[name] = dict(spatial=null_median, labels=label_null)
    return pd.DataFrame(rows), nulls


def within_family_q(p: pd.Series) -> pd.Series:
    """Benjamini-Hochberg over the family's genes only, in the order given."""
    return pd.Series(false_discovery_control(p.to_numpy(float), method="bh"), p.index)


# ===== The table =====


def characterise(parts: list[pd.DataFrame]) -> pd.DataFrame:
    """top_genes.csv: the parts, each one row per gene with a symbol, side by side.

    The first part (chosen_genes) sets the genes and their order.
    """
    out = parts[0]
    for part in parts[1:]:
        out = out.merge(part, on="symbol", how="left")
    return out


def ranking_columns(nano: pd.DataFrame) -> pd.DataFrame:
    """The columns of analysis 1 that top_genes.csv keeps, from gene_ranking.csv."""
    columns = [
        "symbol",
        "n_experiments_used",
        "reliability",
        "n_structures",
        "rho",
        "p_spatial",
        "q_all",
        "rank_all",
        "null_lo",
        "null_hi",
        "t_boot",
    ]
    return nano[columns].rename(columns={"n_experiments_used": "n_experiments"})


def within_columns(within: pd.DataFrame) -> pd.DataFrame:
    """The columns of analysis 2 that top_genes.csv keeps, from within_division.csv."""
    out = within[
        [
            "symbol",
            "rho_division_only",
            "rho_within",
            "p_within_spatial",
            "q_all_spatial",
            "null_lo",
            "null_hi",
        ]
    ]
    return out.rename(
        columns={
            "p_within_spatial": "p_within",
            "q_all_spatial": "q_within",
            "null_lo": "within_null_lo",
            "null_hi": "within_null_hi",
        }
    )


def leftover_columns(leftover: pd.DataFrame) -> pd.DataFrame:
    """The columns of the leftover that top_genes.csv keeps, from leftover_genes.csv."""
    out = leftover[
        [
            "symbol",
            "in_model",
            "rho",
            "p_spatial",
            "q_all",
            "rank_all",
            "null_lo",
            "null_hi",
            "t_boot",
        ]
    ]
    return out.rename(
        columns={
            "rho": "leftover_rho",
            "p_spatial": "leftover_p",
            "q_all": "leftover_q_all",
            "rank_all": "leftover_rank",
            "null_lo": "leftover_null_lo",
            "null_hi": "leftover_null_hi",
            "t_boot": "leftover_t_boot",
        }
    )


def load_top_genes() -> pd.DataFrame:
    """The table run_ish_top_genes wrote, its flags back to booleans."""
    if not TOP_TABLE.exists():
        raise FileNotFoundError(f"{TOP_TABLE} not found: run run_ish_top_genes.py first")
    table = pd.read_csv(TOP_TABLE, keep_default_na=False, na_values=[""])
    for column in ("top", "named", "family"):
        table[column] = table[column].astype(str) == "True"
    for column in ("in_model", "family_sources", "gene_sets", "p9_category"):
        table[column] = table[column].fillna("")
    return table
