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
    Gria1 or density   rho with Gria1, the density term and the main model's
                       prediction on the structures of the fit, with autofluorescence
                       there, and with the PSD95 punctum density where it is
                       measured: whether the gene is just a Gria1-like or
                       density-like map
    the leftover       rho with the main model's leftover and its spatial p, each
                       surrogate put through the same fit (leftover_genes.csv of
                       adult.beyond_density)
    what it takes      the main model with the gene's ranks added as one more
                       straight term, as each of the model's is, scored on held-out
                       structures as the main model is: the share
                       of the reproducible map it takes from the leftover, against
                       surrogate maps with the gene's own smoothness added instead
    annotation         its name, the GO terms that put it in the ontology panel, its
                       GO cellular-component terms at the synapse, its gene sets and
                       P9's category (a label), all from the gene table's records

The genes: those past the spatial null after BH over every gene in analysis 1,
Gria1 and Cacng8 always, and every member of the AMPA receptor complex family that
the leftover table holds, whether it passes or not.

The tests against the leftover were named in advance, before the main model ran
(ish.gene_sets holds the family, its sources and what was seen before):

    1  Cacng8      its rho with the leftover; the uncorrected spatial p is the test
    2  the family  its median rho with the leftover against the same genes' median
                   over the leftover's surrogates (ish.gene_sets.set_test); and its
                   rho against as many genes of the other postsynaptic set outside
                   the family, each the unused one closest in log median energy
                   (ish.panel_test.greedy_match), labels permuted
                   (ish.panel_test.two_sample), one test each; then gene by gene,
                   BH within the family
    3  the rest    every other gene, exploratory, BH over every gene (q_all of
                   leftover_genes.csv), the model's own genes (Gria1, the density genes)
                   among them as the rules set it

The same three are read on the nano map itself, to describe the family, not as tests
of the leftover.

What the tests rest on is set out in docs/ISH_ANALYSIS.md (section 5.7): Cacng8's
p against the leftover of the four-subunit model was seen before it was named, so
tier 1 re-tests a result already seen on a near-identical leftover; Dlg4 and
Camk2a, which the first version's density markers held, are no terms of the main
model (the rule of the density genes leaves them out as AMPA-receptor-linked); and
each surrogate of
the leftover, the model projected out of it, is rougher at short range than the
leftover, so the tests are read again against smoother Gaussian fields projected
the same way (smooth_null_check), and the family again against controls from
outside the model's terms (check_tests).

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
            leftover. It was added after the plain null's numbers were seen

The main model is refitted on the structures of the fit the gene has, every
predictor and the ceiling rebuilt there as the jackknife of adult.beyond_figures
does, and the gene, the main model and every map meet the same folds.

Run by run_ish_top_genes.py.
"""

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, rankdata, spearmanr

from sepmap.adult import beyond_density
from sepmap.adult import synaptome
from sepmap.config import SETTINGS
from sepmap.ish import gene_ranking, gene_sets, gene_table, spatial_null
from sepmap.ish.gene_sets import (
    AMPA_FAMILY,
    GO_AMPA_COMPLEX,
    GO_AMPA_COMPLEX_RELEASE,
    LEFTOVER_GENE,
    PARTNER_SUBUNITS,
    SCHWENK_2012,
    set_test,
)
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.ish.panel_test import greedy_match, two_sample
from sepmap.ish.spatial_null import ALPHA, BAND, null_rho, spatial_p
from sepmap.structures import TABLES

# the structures a gene needs, the BH level, the surrogates added in place of a gene
ISH = SETTINGS["ish"]
ISH_ANALYSIS = SETTINGS["ish_analysis"]
TOP_GENES = SETTINGS["top_genes"]

TOP_TABLE = TABLES / "top_genes.csv"
NAMED_TESTS = TABLES / "named_tests.csv"
FAMILY_TABLE = TABLES / "family_members.csv"
NUMBERS = numbers_path("top_genes")
NULL_CHECK = TABLES / "leftover_null_check.csv"

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

# the Allen experiments of the family members the gene table does not hold, counted
# when the family was named (ish.gene_sets): none, or some outside both panels, where
# adding a gene would change the rules
UNTABLED_EXPERIMENTS = {
    "Olfm3": 0,
    "Prrt1": 0,
    "Prrt2": 0,
    "Shisa8": 0,
    "Gsg1l": 2,
    "Rap2b": 2,
}

# the columns of a gene's row that the numbers file carries for Cacng8 and Gria1:
# column, decimals, what
GENE_NUMBERS = (
    ("rho", 3, "rho with the map"),
    ("p_spatial", 6, "its spatial p"),
    ("rank_all", 0, "its rank among all genes"),
    ("rho_within", 3, "mean rho inside divisions"),
    ("p_within", 6, "its spatial p"),
    ("rho_variants_min", 3, "lowest rho over the robustness variants"),
    ("rho_variants_max", 3, "highest rho over the robustness variants"),
    ("rank_variants_min", 0, "best rank over the variants holding every gene"),
    ("rank_variants_max", 0, "worst rank over the variants holding every gene"),
    ("reliability", 3, "agreement of its Allen experiments"),
    ("rho_Gria1", 3, "rho with Gria1, structures of the fit"),
    ("rho_density", 3, "rho with the density term"),
    ("rho_autofluo", 3, "rho with autofluorescence, structures of the fit"),
    ("rho_prediction", 3, "rho with the main model's prediction"),
    ("rho_psd95", 3, "rho with PSD95 punctum density, where measured"),
    ("leftover_rho", 3, "rho with the main model's leftover"),
    ("leftover_p", 6, "its spatial p, each surrogate through the same fit"),
    ("leftover_rank", 0, "its rank against the leftover"),
    ("left_main", 4, "share of the reproducible map the main model leaves there"),
    ("left_with_gene", 4, "the same with the gene added"),
    ("taken", 4, "share of the reproducible map the gene takes from the leftover"),
    ("taken_null_hi", 4, "95% of plain surrogates of the gene take less than this"),
    ("p_taken", 6, "share of those maps that take as much, one-sided"),
    ("taken_null_alike_hi", 4, "95% of maps alike to the model take less than this"),
    ("p_taken_alike", 6, "share of those maps that take as much, one-sided"),
)

# the named tests of tier 2 as named_tests.csv writes them
SPATIAL_TEST = "the family's median against the surrogates"
MATCHED_TEST = "the family against matched postsynaptic controls"
WITHOUT_TEST = f"the family without {LEFTOVER_GENE} against the surrogates"
OUTSIDE_TEST = "the family against matched controls outside the model"


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
            tier = "1 named for the leftover"
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


def choose_genes(nano: pd.DataFrame, leftover: pd.DataFrame, q: float) -> tuple:
    """The genes characterised, the family's members, and each member's control.

    `nano` is the nano map's rows of gene_ranking.csv, `leftover` leftover_genes.csv.
    Returns the chosen genes (chosen_genes), family_members.csv's table, the members
    tested, the pool of other postsynaptic genes against the leftover, the pairs of
    member and control (family_controls) and every gene's median energy.
    """
    table = gene_table.load_gene_table()
    genes = gene_table.per_gene(table)
    chosen = chosen_genes(nano, leftover, q)
    members = family_members(genes, leftover)
    family = list(members.loc[members["tested"], "symbol"])
    _, level = gene_table.gene_levels(table)
    tested = set(leftover["symbol"])
    postsynaptic = gene_sets.members_from_table(genes)["other postsynaptic"]
    pool = [g for g in postsynaptic if g in tested]
    pairs = family_controls(family, pool, level)
    return chosen, members, family, pool, pairs, level


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
        elif symbol not in listed and UNTABLED_EXPERIMENTS.get(symbol, 0) == 0:
            reason = "no Allen experiment"
        elif symbol not in listed:
            n = UNTABLED_EXPERIMENTS[symbol]
            reason = (
                f"{n} Allen experiments, but in neither panel of the gene table; adding "
                "it would change the rules"
            )
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
                go_panel_terms=gene_table.go_term_list(panel, names),
                go_synapse_terms=gene_table.go_term_list(synaptic, names),
                gene_sets=doc["gene_sets"].get(symbol, ""),
                p9_category=doc["p9_category"].get(symbol, ""),
            )
        )
    return pd.DataFrame(rows).fillna("")


def annotation(symbols: list[str]) -> tuple[pd.DataFrame, str]:
    """annotation_table of `symbols` from the gene table's cached records and GO.

    Never asks a server: the records and go-basic.obo are those step 15 cached.
    Returns the table and the release of go-basic.obo.
    """
    documentation = pd.read_csv(
        gene_table.DOCUMENTATION, encoding="utf-8-sig", keep_default_na=False
    )
    records = gene_table.mygene_records(symbols, offline=True)
    parents, names, release = gene_table.load_obo(offline=True)
    table = annotation_table(symbols, documentation, records, parents, names)
    return table, release


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


def gene_on_fit(
    inputs: beyond_density.Inputs, symbol: str
) -> tuple[list[int], np.ndarray]:
    """The positions of the structures of the fit the gene has, and its values there."""
    profile = inputs.expr.get(symbol, {})
    columns = [
        i for i, s in enumerate(inputs.structures) if np.isfinite(profile.get(s, np.nan))
    ]
    values = np.array([profile[inputs.structures[i]] for i in columns], dtype=float)
    return columns, values


def predictor_names(terms: dict[str, tuple[str, ...]]) -> list[str]:
    """The main model's predictors by name, in the order the model lists them.

    Autofluorescence follows, a description: it is no term of the main model.
    """
    names = [name for group in beyond_density.ORDER for name in terms.get(group, ())]
    return names + [n for n in ("autofluo",) if n not in names]


def predictor_rhos(inputs: beyond_density.Inputs, symbols: list[str]) -> pd.DataFrame:
    """Per gene, rho with each term of the main model, its prediction, and PSD95.

    On the structures of the fit the gene has (n_fit): rho with Gria1, the density
    term, autofluorescence and the main model's in-sample prediction; and on those
    of them where the PSD95 punctum density is measured (n_psd95), rho with it.
    """
    covariates, _, _ = beyond_density.covariates_for(inputs)
    y = beyond_density.full_map(inputs.nano)
    prediction = y - beyond_density.residual(
        y, beyond_density.model_columns(covariates, inputs.terms)
    )
    psd95 = inputs.synapses.reindex(columns=[synaptome.MEASURE])[synaptome.MEASURE]
    psd95 = psd95.to_numpy(float)
    rows = []
    for symbol in symbols:
        columns, values = gene_on_fit(inputs, symbol)
        row = dict(symbol=symbol, n_fit=len(columns))
        for name in predictor_names(inputs.terms):
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
    hat = beyond_density.projection(design)
    fitted = hat @ gene
    rest = gene - fitted
    surr = spatial_null.surrogates(rest, d, n, seed=seed).astype(float)
    surr = surr - surr @ hat.T
    surr = surr / surr.std(axis=1, keepdims=True) * rest.std()
    return fitted[None, :] + surr


def added_share(
    inputs: beyond_density.Inputs,
    symbol: str,
    centroids: pd.DataFrame,
    n_surrogates: int | None = None,
    seed: int = 0,
) -> tuple[dict, dict[str, np.ndarray]]:
    """What the gene takes of the leftover when added to the main model, and its nulls.

    On the structures of the fit the gene has, the main model and the main model
    with the gene's ranks as one more straight term are scored on the same folds
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
    columns, values = gene_on_fit(inputs, symbol)
    sub = inputs
    if len(columns) < len(inputs.structures):
        sub = beyond_density.restrict(inputs, [inputs.structures[i] for i in columns])
    _, explainable = beyond_density.ceiling(sub.nano, beyond_density.half_splits())
    y = beyond_density.full_map(sub.nano)
    covariates, _, _ = beyond_density.covariates_for(sub)
    xs = beyond_density.model_columns(covariates, sub.terms)
    labels = beyond_density.fold_labels(len(y))

    def left_with(extra: np.ndarray | None) -> float:
        """The share of the reproducible map left, with `extra` ranked and added."""
        used = xs if extra is None else xs + [rankdata(extra)]
        return 1 - beyond_density.cv_r2(y, used, labels) / explainable

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


def spatial_row(
    map_name: str,
    test: str,
    role: str,
    genes: list[str],
    spatial: dict,
    null_median: np.ndarray,
) -> dict:
    """A row of named_tests.csv for a group's median rho against the surrogates.

    `spatial` is ish.gene_sets.set_test's row of `genes`, `null_median` their median
    over each surrogate; second is the null's median.
    """
    second = float(np.median(null_median))
    return dict(
        tier="2",
        map=map_name,
        test=test,
        role=role,
        n_first=len(genes),
        n_second=np.nan,
        first=spatial["median_rho"],
        second=second,
        difference=spatial["median_rho"] - second,
        null_lo=spatial["null_lo"],
        null_hi=spatial["null_hi"],
        critical=np.nan,
        p=spatial["p_spatial"],
    )


def matched_row(
    map_name: str,
    test: str,
    role: str,
    first: np.ndarray,
    second: np.ndarray,
    rng: np.random.Generator,
) -> tuple[dict, np.ndarray]:
    """A row of named_tests.csv for a group against its controls, labels permuted.

    `first` and `second` are the two sides' rho; critical is the difference beyond
    which the label-permutation p falls below 0.05. Returns the row and the null.
    """
    difference, p, label_null = two_sample(first, second, rng)
    row = dict(
        tier="2",
        map=map_name,
        test=test,
        role=role,
        n_first=len(first),
        n_second=len(second),
        first=float(np.median(first)),
        second=float(np.median(second)),
        difference=difference,
        null_lo=float(np.percentile(label_null, BAND[0])),
        null_hi=float(np.percentile(label_null, BAND[1])),
        critical=float(np.percentile(np.abs(label_null), 95)),
        p=p,
    )
    return row, label_null


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
        rows.append(spatial_row(name, SPATIAL_TEST, role, family, spatial, null_median))
        a = rho[family].to_numpy(float)
        b = rho[controls].to_numpy(float)
        rng = np.random.default_rng(child)
        row, label_null = matched_row(name, MATCHED_TEST, role, a, b, rng)
        rows.append(row)
        nulls[name] = dict(spatial=null_median, labels=label_null)
    return pd.DataFrame(rows), nulls


def check_tests(
    family: list[str],
    leftover: pd.DataFrame,
    null: np.ndarray,
    null_genes: list[str],
    pool: list[str],
    level: dict[str, float],
    seed: int = 1,
) -> pd.DataFrame:
    """Two check rows of tier 2 on the leftover, added after the tests ran.

    The family without Cacng8 against the same surrogates, so a reader sees whether
    the group test says more than tier 1; and the family against controls matched
    afresh from the pool's genes that are not in the main model (in_model empty),
    since a control inside the density term is projected out of the leftover with
    the model. `leftover` is leftover_genes.csv, `null` its genes x surrogates.
    """
    rho = leftover.set_index("symbol")["rho"]
    row_of = {g: i for i, g in enumerate(null_genes)}
    without = [g for g in family if g != LEFTOVER_GENE]
    spatial, null_median = set_test(without, rho, null, row_of)
    in_model = leftover.set_index("symbol")["in_model"]
    outside = [g for g in pool if not in_model.get(g, "")]
    pairs = family_controls(family, outside, level)
    a = rho[family].to_numpy(float)
    b = rho[[pairs[g] for g in family if g in pairs]].to_numpy(float)
    rng = np.random.default_rng(seed)
    rows = [
        spatial_row("leftover", WITHOUT_TEST, "check", without, spatial, null_median),
        matched_row("leftover", OUTSIDE_TEST, "check", a, b, rng)[0],
    ]
    return pd.DataFrame(rows)


def named_tests(
    family: list[str],
    pairs: dict[str, str],
    leftover: pd.DataFrame,
    null: dict,
    nano: pd.DataFrame,
    nano_null: dict,
    pool: list[str],
    level: dict[str, float],
) -> tuple[pd.DataFrame, dict[str, dict[str, np.ndarray]]]:
    """named_tests.csv: tier 1 and tier 2 on the leftover; the group tests on the map.

    Tier 2 as a group with its check rows; on the map as a description.

    `null` and `nano_null` hold each map's null rho (rho) and the genes of its rows
    (genes). Returns the table and the group tests' nulls (group_tests).
    """
    maps = dict(
        leftover=(leftover.set_index("symbol")["rho"], null["rho"], null["genes"]),
        nano=(nano.set_index("symbol")["rho"], nano_null["rho"], nano_null["genes"]),
    )
    group, group_nulls = group_tests(family, pairs, maps)
    checks = check_tests(family, leftover, null["rho"], null["genes"], pool, level)
    tier_one_row = pd.DataFrame([tier_one(leftover)])
    tests = pd.concat([tier_one_row, group, checks], ignore_index=True)
    return tests, group_nulls


def neighbour_agreement(maps: np.ndarray, d: np.ndarray, k: int) -> np.ndarray:
    """Per map (a row), Spearman of each structure with the mean of its k nearest.

    How smooth a map is at short range: near 1 when neighbours carry the same value,
    near 0 for a shuffled map. `d` is the structures' distance matrix.
    """
    nearest = np.argsort(d, axis=1)[:, 1 : k + 1]
    neighbours = maps[:, nearest].mean(axis=2)
    a = rankdata(maps, axis=1)
    b = rankdata(neighbours, axis=1)
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean(axis=1, keepdims=True)
    return (a * b).sum(axis=1) / np.sqrt((a**2).sum(axis=1) * (b**2).sum(axis=1))


def null_check_row(
    name: str,
    range_mm: float,
    maps: np.ndarray,
    leftover: pd.DataFrame,
    vectors: dict[str, tuple[np.ndarray, np.ndarray]],
    family: list[str],
    d: np.ndarray,
) -> dict:
    """The tests against one null of the leftover.

    Cacng8's p, the family's p, and how many genes fall below p 0.05, with how smooth the
    null's maps are.
    """
    rho = leftover.set_index("symbol")["rho"]
    nulls = {
        g: null_rho(maps[:, columns], values) for g, (columns, values) in vectors.items()
    }
    p = {g: spatial_p(rho[g], nulls[g]) for g in vectors}
    family = [g for g in family if g in nulls]
    family_null = np.median([nulls[g] for g in family], axis=0)
    return dict(
        null=name,
        range_mm=range_mm,
        n_maps=len(maps),
        neighbour_rho=float(
            np.median(neighbour_agreement(maps, d, TOP_GENES["neighbours"]))
        ),
        p_cacng8=p[LEFTOVER_GENE],
        p_family=spatial_p(float(np.median(rho[family])), family_null),
        genes_p05=int(sum(v < ALPHA for v in p.values())),
        genes=len(p),
    )


def smooth_null_check(
    inputs: beyond_density.Inputs,
    residual: np.ndarray,
    surrogates: np.ndarray,
    leftover: pd.DataFrame,
    family: list[str],
    centroids: pd.DataFrame,
    seed: int = 0,
) -> pd.DataFrame:
    """The leftover's tests against nulls smoother than its own surrogates.

    The step's surrogates are the leftover's, each with the main model projected
    out, which leaves them rougher at short range than the leftover. So the tests
    are run again against Gaussian random fields with an exponential covariance of
    each range in top_genes.null_check_ranges_mm (ish.spatial_null.random_fields),
    top_genes.null_check_maps of them, each projected through the same fit. Rows:
    the leftover itself (its smoothness only), the leftover's surrogates and each range.
    """
    s = inputs.structures
    d = spatial_null.distance_matrix(centroids.loc[s])
    covariates, _, _ = beyond_density.covariates_for(inputs)
    columns = beyond_density.model_columns(covariates, inputs.terms)
    hat = beyond_density.projection(np.column_stack(columns + [np.ones(len(s))]))
    vectors = gene_ranking.gene_vectors(inputs.expr, s)
    smooth = neighbour_agreement(residual[None, :], d, TOP_GENES["neighbours"])[0]
    rows = [
        dict(null="the leftover itself", range_mm=np.nan, n_maps=1, neighbour_rho=smooth),
        null_check_row(
            "the leftover's surrogates", np.nan, surrogates, leftover, vectors, family, d
        ),
    ]
    rng = np.random.default_rng(seed)
    for range_mm in TOP_GENES["null_check_ranges_mm"]:
        fields = spatial_null.random_fields(
            d, (0.0, 1.0, float(range_mm)), TOP_GENES["null_check_maps"], rng
        )
        fields = fields - fields @ hat.T
        rows.append(
            null_check_row(
                "Gaussian fields", float(range_mm), fields, leftover, vectors, family, d
            )
        )
    return pd.DataFrame(rows)


def within_family_q(p: pd.Series) -> pd.Series:
    """Benjamini-Hochberg over the family's genes only, in the order given."""
    return pd.Series(false_discovery_control(p.to_numpy(float), method="bh"), p.index)


# ===== The table =====


def characterise(parts: list[pd.DataFrame]) -> pd.DataFrame:
    """top_genes.csv: the parts, each one row per gene with a symbol, side by side.

    The first part (chosen_genes) sets the genes and their order. A text column of
    a part that does not hold every gene (a family member's sources and control,
    a gene's term of the model) is empty, not NaN, for the others.
    """
    out = parts[0]
    for part in parts[1:]:
        out = out.merge(part, on="symbol", how="left")
    for column in ("family_sources", "in_model", "matched_control"):
        out[column] = out[column].fillna("")
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


def control_columns(
    pairs: dict[str, str],
    leftover: pd.DataFrame,
    nano: pd.DataFrame,
    level: dict[str, float],
) -> pd.DataFrame:
    """Per family member, its matched control, both levels and the control's rho.

    `leftover` is leftover_genes.csv, `nano` the map's rows of gene_ranking.csv,
    `level` each gene's median energy.
    """
    left_rho = leftover.set_index("symbol")["rho"]
    map_rho = nano.set_index("symbol")["rho"]
    rows = []
    for gene, control in pairs.items():
        rows.append(
            dict(
                symbol=gene,
                median_energy=level.get(gene, np.nan),
                matched_control=control,
                control_median_energy=level.get(control, np.nan),
                control_rho=map_rho[control],
                control_leftover_rho=left_rho[control],
            )
        )
    return pd.DataFrame(rows)


def family_q_table(
    family: list[str], leftover: pd.DataFrame, nano: pd.DataFrame
) -> pd.DataFrame:
    """Tier 2 gene by gene: BH within the family, against the leftover and the map."""
    p_leftover = leftover.set_index("symbol")["p_spatial"].reindex(family)
    p_map = nano.set_index("symbol")["p_spatial"].reindex(family)
    return pd.DataFrame(
        dict(
            symbol=family,
            leftover_q_family=within_family_q(p_leftover).to_numpy(),
            q_family=within_family_q(p_map).to_numpy(),
        )
    )


# ===== The numbers for the text =====


def top_numbers(table: pd.DataFrame, q: float) -> list[tuple]:
    """The numbers of the genes past the map's null, and of Cacng8 and Gria1."""
    top = table[table["top"]]
    free = top[top["in_model"] != "abundance"]
    follow = free[free["leftover_p"] < ALPHA]
    take = free[free["p_taken"] < ALPHA]
    take_alike = free[free["p_taken_alike"] < ALPHA]
    rows = [
        ("characterised", len(table), "genes characterised"),
        ("top_genes", len(top), f"genes past the map's null, BH over all, q < {q}"),
        ("top_list", " ".join(top["symbol"]), "those genes, best first"),
        (
            "top_within",
            int((top["q_within"] < q).sum()),
            f"of them, past the null inside divisions, q < {q}",
        ),
        (
            "top_follow_leftover",
            len(follow),
            "of them, not a term of the model, rho with the leftover at p < 0.05",
        ),
        ("top_follow_leftover_list", " ".join(follow["symbol"]), "those genes"),
        (
            "top_take",
            len(take),
            "of them, not a term of the model, taking more than 95% of its "
            "plain surrogates",
        ),
        ("top_take_list", " ".join(take["symbol"]), "those genes"),
        (
            "top_take_alike",
            len(take_alike),
            "of them, not a term of the model, taking more than 95% of maps alike to it",
        ),
        ("top_take_alike_list", " ".join(take_alike["symbol"]), "those genes"),
        (
            "added_surrogates",
            TOP_GENES["n_added_surrogates"],
            "maps of a gene's smoothness added in its place",
        ),
    ]
    by_gene = table.set_index("symbol")
    for gene in NAMED_GENES:
        for column, digits, what in GENE_NUMBERS:
            value = by_gene.loc[gene, column]

            # adding 0.0 turns a rounded -0.0 into 0.0
            value = int(value) if digits == 0 else round(float(value), digits) + 0.0
            rows.append((f"{column}_{gene}", value, f"{gene}: {what}"))
    return rows


def group_numbers(
    table: pd.DataFrame, tier2: pd.DataFrame, map_name: str, q: float
) -> list[tuple]:
    """The numbers of tier 2's group tests on one map, and its members past BH there.

    BH runs within the family. `tier2` is named_tests.csv's tier 2, by map and test.
    """
    spatial = tier2.loc[(map_name, SPATIAL_TEST)]
    matched = tier2.loc[(map_name, MATCHED_TEST)]
    key = f"family_{map_name}"
    column = "leftover_q_family" if map_name == "leftover" else "q_family"
    rho = "leftover_rho" if map_name == "leftover" else "rho"
    past = table[table["family"] & (table[column] < q)]
    return [
        (f"{key}_median", round(spatial["first"], 3), "the family's median rho"),
        (f"{key}_null_median", round(spatial["second"], 3), "over the surrogates"),
        (f"{key}_null_lo", round(spatial["null_lo"], 3), "2.5% of that null"),
        (f"{key}_null_hi", round(spatial["null_hi"], 3), "97.5% of that null"),
        (f"{key}_p", round(spatial["p"], 6), "the family's spatial p"),
        (
            f"{key}_controls_median",
            round(matched["second"], 3),
            "the matched controls' median rho",
        ),
        (
            f"{key}_controls_difference",
            round(matched["difference"], 3),
            "family minus controls, medians",
        ),
        (
            f"{key}_controls_critical",
            round(matched["critical"], 3),
            "the difference that would give p < 0.05",
        ),
        (f"{key}_controls_p", round(matched["p"], 6), "label-permutation p"),
        (f"{key}_pass_within", len(past), f"members past BH within the family, q < {q}"),
        (f"{key}_pass_within_list", " ".join(past["symbol"]), "those members"),
        (
            f"{key}_pass_within_rho",
            " ".join(f"{v:+.3f}" for v in past[rho]),
            "their rho, in the same order",
        ),
    ]


def check_numbers(tier2: pd.DataFrame) -> list[tuple]:
    """The numbers of tier 2's check rows on the leftover."""
    without = tier2.loc[("leftover", WITHOUT_TEST)]
    outside = tier2.loc[("leftover", OUTSIDE_TEST)]
    return [
        (
            "family_leftover_without_cacng8_median",
            round(without["first"], 3),
            "check: the family's median without Cacng8",
        ),
        ("family_leftover_without_cacng8_p", round(without["p"], 6), "check: its p"),
        (
            "family_leftover_outside_controls",
            int(outside["n_second"]),
            "check: controls matched from genes outside the model's terms",
        ),
        (
            "family_leftover_outside_controls_median",
            round(outside["second"], 3),
            "check: their median rho with the leftover",
        ),
        (
            "family_leftover_outside_difference",
            round(outside["difference"], 3),
            "check: family minus those controls, medians",
        ),
        (
            "family_leftover_outside_p",
            round(outside["p"], 6),
            "check: label-permutation p",
        ),
    ]


def family_numbers(
    table: pd.DataFrame, members: pd.DataFrame, tests: pd.DataFrame, q: float
) -> list[tuple]:
    """The numbers of the family: its members, its group tests and their check rows."""
    untested = members[~members["tested"]]
    tier2 = tests[tests["tier"] == "2"].set_index(["map", "test"])
    rows = [
        ("family_listed", len(members), "members of the AMPA receptor complex family"),
        ("family_tested", int(members["tested"].sum()), "of them tested"),
        ("family_not_tested", " ".join(untested["symbol"]), "not tested"),
        (
            "family_controls",
            int(tier2.loc[("leftover", MATCHED_TEST), "n_second"]),
            "matched postsynaptic controls",
        ),
    ]
    for map_name in ("leftover", "nano"):
        rows += group_numbers(table, tier2, map_name, q)
    return rows + check_numbers(tier2)


def null_check_numbers(check: pd.DataFrame) -> list[tuple]:
    """The numbers of the leftover's tests against smoother nulls."""
    rows = []
    for r in check.itertuples():
        if r.null == "the leftover itself":
            key = "null_check_leftover"
        elif np.isfinite(r.range_mm):
            key = f"null_check_{r.range_mm:g}mm"
        else:
            key = "null_check_surrogates"
        rows.append(
            (f"{key}_neighbour_rho", round(r.neighbour_rho, 3), f"{r.null}: smoothness")
        )
        if r.null == "the leftover itself":
            continue
        rows += [
            (f"{key}_p_cacng8", round(r.p_cacng8, 6), f"{r.null}: Cacng8's p"),
            (f"{key}_p_family", round(r.p_family, 6), f"{r.null}: the family's p"),
            (f"{key}_genes_p05", int(r.genes_p05), f"{r.null}: genes below p 0.05"),
        ]
    return rows


def numbers_table(
    table: pd.DataFrame,
    members: pd.DataFrame,
    tests: pd.DataFrame,
    check: pd.DataFrame,
    q: float,
) -> pd.DataFrame:
    """numbers_top_genes.csv: the numbers of this step that the text quotes."""
    rows = top_numbers(table, q)
    rows += family_numbers(table, members, tests, q)
    rows += null_check_numbers(check)
    return numbers_frame(rows)


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
