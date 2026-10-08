"""The kinds of genes of analysis 3, fixed before any correlation with the map is read.

The April result grouped genes by categories written while the genes were being
chosen (gene_targets.csv), and split one of them by hand after looking; its
p = 0.032 rested on that split. Here the groups come from annotation nobody here
wrote, or from a cited marker list, and this file was committed before any gene
was correlated with the map on the new inputs (git log shows the order):

    subunits             GO:0004971 AMPA glutamate receptor activity, intersected
                         with GO:0032281 AMPA glutamate receptor complex, minus
                         Grid1 and Grid2 (ish.panel_build's rule, role subunit)
    localisation         GO:0099072, GO:0099645, GO:0097113, GO:0098970, and
                         GO:0032281 minus the subunits (panel_build's rule, role
                         localisation): the genes that bring AMPA receptors to the
                         postsynaptic membrane and hold them there
    other postsynaptic   annotated (cellular component) to GO:0098794 postsynapse
                         or any term below it, by is_a or part_of, and not to the
                         presynapse; not a subunit or localisation gene
    presynaptic          annotated to GO:0098793 presynapse or below, not to the
                         postsynapse; not a subunit or localisation gene, so the
                         four GO sets share no gene
    GABAergic markers    GABA synthesis and vesicular transport (Gad1, Gad2,
                         Slc32a1), and the GABAergic subclass markers of Tasic et
                         al. 2018 that have a grid here (of Lamp5, Sncg, Serpinf1,
                         Vip, Sst, Pvalb and Meis2: Sst and Pvalb)
    glia                 P9's control_glia, unchanged (Aldh1l1, Aqp4, Gfap, Olig2,
                         Sox10, Mbp, Cx3cr1)

GO is read with its ancestry, from a dated go-basic.obo (ish.gene_table), so a gene
annotated to "postsynaptic density" counts as postsynaptic. Annotations with a NOT
qualifier are left out. SynGO's curated annotations enter GO, so the GO sets carry
those that reached the records of mygene.info; SynGO's own release is not read (it
ships as spreadsheets, and venv_atlas has no reader for them). GO has no cell-type
terms, so the two control sets come from a cited list and from P9's own glia list,
neither chosen by correlation. A gene may sit in one GO set and one marker set
(Slc32a1 is presynaptic by GO and a GABAergic marker); a contrast between two sets
leaves out the genes in both, and says so.

Two contrasts are named in advance (the plan of 8 October, before any correlation on
the new inputs): the three postsynaptic sets pooled against the presynaptic set, and
against glia. A postsynaptic-like map, as any glutamate receptor label should give,
puts the first side above the second in both.

One more group is drawn beside the six for context and never tested: genes GO
annotates to both the presynapse and the postsynapse. The rule "not both" leaves the
classic vesicle genes (Syn1, Vamp2, Stx1a, Bsn, Syt1, Cplx1), Camk2a and Slc17a7 out
of the presynaptic set, which keeps 16 genes, mostly inhibitory and neuromodulatory;
a reader looks for them there. It was added after the sets' membership was read, and
before any set was correlated with the map.

The tests of analysis 3 read each gene's rho with the map and with every surrogate
of it (ish.gene_ranking), so every gene of a set meets the same surrogate and the
co-expression of a set's genes is kept in its null:

    set median     a set's median rho, against the median of the same genes' rho
                   with every surrogate; two-sided, BH over the six sets; sets with
                   fewer than ish_analysis.min_set_genes genes are drawn, not tested
    contrast       the difference of the two sides' medians against the same
                   difference for every surrogate; beside it, a permutation of the
                   labels between the sides, which treats co-expressed genes as
                   independent draws (the kind of null of the word test of April)
    localisation   ish.panel_test's design on the new inputs: each gene's partial
                   rho with the map once the subunit composite (the mean rank of
                   the subunits) is removed from both; the localisation genes
                   against the same number of other postsynaptic genes matched on
                   expression (greedily, on log median energy), labels permuted
                   between the two; a positive control (control genes with a
                   reproducible map against those without, on |rho|) that must
                   come out, or the test detects nothing; and the matched
                   difference against the surrogates, each surrogate in place of
                   the map

Membership is computed by run_ish_gene_table.py, with the gene table; the context
group by run_ish_gene_sets.py, from the same cached GO records. The tests are run by
run_ish_gene_sets.py.
"""

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, rankdata, spearmanr

from sepmap.config import SETTINGS
from sepmap.ish.panel_test import greedy_match, partial, two_sample
from sepmap.ish.spatial_null import spatial_p
from sepmap.structures import TABLES

# the genes a set needs to be tested, the BH level; the structures, permutations and
# reliability of the localisation test
ISH_ANALYSIS = SETTINGS["ish_analysis"]
ISH_PANEL_TEST = SETTINGS["ish_panel_test"]

MEMBERS = TABLES / "gene_sets.csv"
SET_TESTS = TABLES / "set_tests.csv"
CONTRAST_TESTS = TABLES / "contrasts.csv"
LOCALISATION = TABLES / "localisation_test.csv"
LOCALISATION_SUMMARY = TABLES / "localisation_summary.csv"

# the percentiles of a null distribution that a value must leave to pass at 0.05
BAND = (2.5, 97.5)

# the sets in the order the figures draw them
SET_ORDER = (
    "subunits",
    "localisation",
    "other postsynaptic",
    "presynaptic",
    "GABAergic markers",
    "glia",
)

# the roles of ish.panel_build that define the first two sets
PANEL_ROLE_SETS = {"subunits": "subunit", "localisation": "localisation"}

# the two cellular components of the next two sets
POSTSYNAPSE = "GO:0098794"
PRESYNAPSE = "GO:0098793"

# GABA synthesis and vesicular transport, and the GABAergic subclass markers named by
# Tasic et al. 2018 (Nature 563:72); only those with a grid enter the set
GABA_MACHINERY = ("Gad1", "Gad2", "Slc32a1")
TASIC_GABA_MARKERS = ("Lamp5", "Sncg", "Serpinf1", "Vip", "Sst", "Pvalb", "Meis2")

# P9's control_glia, unchanged
GLIA = ("Aldh1l1", "Aqp4", "Gfap", "Olig2", "Sox10", "Mbp", "Cx3cr1")

# where each set comes from, as the gene table and the figures write it
GENE_SETS = {
    "subunits": "GO:0004971 and GO:0032281, minus Grid1 and Grid2 (panel role subunit)",
    "localisation": (
        "GO:0099072, GO:0099645, GO:0097113, GO:0098970, GO:0032281 minus the "
        "subunits (panel role localisation)"
    ),
    "other postsynaptic": (
        "GO:0098794 postsynapse or below, not presynapse, not in the two sets above"
    ),
    "presynaptic": (
        "GO:0098793 presynapse or below, not postsynapse, not in the first two sets"
    ),
    "GABAergic markers": "Gad1, Gad2, Slc32a1, and Tasic et al. 2018 subclass markers",
    "glia": "P9's control_glia: Aldh1l1, Aqp4, Gfap, Olig2, Sox10, Mbp, Cx3cr1",
}

# the contrasts named in advance: the postsynaptic sets pooled, against the presynaptic
# set and against glia; a gene on both sides of a contrast is left out of it
POSTSYNAPTIC_SETS = ("subunits", "localisation", "other postsynaptic")
CONTRASTS = {
    "postsynaptic against presynaptic": (POSTSYNAPTIC_SETS, ("presynaptic",)),
    "postsynaptic against glia": (POSTSYNAPTIC_SETS, ("glia",)),
}

# the group drawn for context and never tested: annotated to both synapse sides
CONTEXT_SET = "pre- and postsynaptic"
CONTEXT_RULE = (
    "GO:0098793 presynapse and GO:0098794 postsynapse (or below) both, not in the "
    "first two sets; context, not tested"
)


def gene_sets(
    genes: list[str], role: dict[str, str], components: dict[str, set[str]]
) -> dict[str, list[str]]:
    """The members of each set among `genes`, {set: sorted symbols}.

    `role` is each gene's role in the ontology panel ('' for a gene outside it),
    `components` each gene's cellular-component terms with their ancestors (empty
    for a gene without a GO record).
    """
    out = {}
    for name, panel_role in PANEL_ROLE_SETS.items():
        out[name] = sorted(g for g in genes if role.get(g, "") == panel_role)
    by_role = set(out["subunits"]) | set(out["localisation"])
    post = {g for g in genes if POSTSYNAPSE in components.get(g, set())}
    pre = {g for g in genes if PRESYNAPSE in components.get(g, set())}
    out["other postsynaptic"] = sorted(post - pre - by_role)
    out["presynaptic"] = sorted(pre - post - by_role)
    markers = GABA_MACHINERY + TASIC_GABA_MARKERS
    out["GABAergic markers"] = sorted(g for g in genes if g in markers)
    out["glia"] = sorted(g for g in genes if g in GLIA)
    return {name: out[name] for name in SET_ORDER}


def context_set(
    genes: list[str], role: dict[str, str], components: dict[str, set[str]]
) -> list[str]:
    """The genes annotated to both the presynapse and the postsynapse, sorted.

    Subunit and localisation genes are left out, as from the two GO sets; the
    arguments are those of gene_sets.
    """
    by_role = {g for g in genes if role.get(g, "") in PANEL_ROLE_SETS.values()}
    both = {
        g
        for g in genes
        if {POSTSYNAPSE, PRESYNAPSE} <= components.get(g, set()) and g not in by_role
    }
    return sorted(both)


def contrast_sides(
    members: dict[str, list[str]], contrast: str
) -> tuple[list[str], list[str], list[str]]:
    """The genes of each side of a contrast named in advance, and those left out.

    A gene in a set of both sides (a marker set beside a GO set) is in neither, and
    comes back in the third list.
    """
    first, second = CONTRASTS[contrast]
    a = {g for name in first for g in members[name]}
    b = {g for name in second for g in members[name]}
    both = a & b
    return sorted(a - both), sorted(b - both), sorted(both)


# ===== Sets against the map =====


def members_from_table(genes: pd.DataFrame) -> dict[str, list[str]]:
    """{set: sorted symbols} from the gene_sets column of the per-gene table."""
    out = {name: [] for name in SET_ORDER}
    for symbol, text in zip(genes["symbol"], genes["gene_sets"]):
        for name in str(text).split(";"):
            name = name.strip()
            if name in out:
                out[name].append(symbol)
    return {name: sorted(symbols) for name, symbols in out.items()}


def member_table(members: dict[str, list[str]], ranking: pd.DataFrame) -> pd.DataFrame:
    """gene_sets.csv: one row per set and gene, with the gene's rho with both maps."""
    by_map = {
        name: ranking[ranking["map"] == name].set_index("symbol")
        for name in ("nano", "auto")
    }
    rows = []
    for name, symbols in members.items():
        for symbol in symbols:
            if symbol not in by_map["nano"].index:
                continue
            nano = by_map["nano"].loc[symbol]
            rows.append(
                dict(
                    gene_set=name,
                    symbol=symbol,
                    p9_gene=bool(nano["p9_gene"]),
                    rho_nano=nano["rho"],
                    rho_auto=by_map["auto"].loc[symbol, "rho"],
                    q_all_nano=nano["q_all"],
                )
            )
    return pd.DataFrame(rows)


def set_test(
    genes: list[str], rho: pd.Series, null: np.ndarray, row_of: dict[str, int]
) -> tuple[dict, np.ndarray]:
    """One set's median rho against its median over the surrogates.

    `rho` is every gene's rho with the map, `null` genes x surrogates of rho with
    every surrogate, `row_of` each gene's row there. Returns the row of the table
    and the null medians.
    """
    rows = [row_of[g] for g in genes]
    observed = float(np.median(rho[genes]))
    null_median = np.median(null[rows], axis=0)
    row = dict(
        n_genes=len(genes),
        median_rho=observed,
        null_lo=float(np.percentile(null_median, BAND[0])),
        null_hi=float(np.percentile(null_median, BAND[1])),
        p_spatial=spatial_p(observed, null_median),
    )
    return row, null_median


def set_tests(
    members: dict[str, list[str]],
    rankings: dict[str, pd.Series],
    nulls: dict[str, np.ndarray],
    null_genes: list[str],
) -> pd.DataFrame:
    """set_tests.csv: per map and set, its median rho, null band, p and q.

    `rankings` and `nulls` hold, per map, every gene's rho and its null rho
    (genes x surrogates, rows in the order of `null_genes`). A set with fewer than
    ish_analysis.min_set_genes genes is listed with tested False and no p; q is
    Benjamini-Hochberg over the tested sets of a map.
    """
    row_of = {g: i for i, g in enumerate(null_genes)}
    parts = []
    for map_name, rho in rankings.items():
        part = []
        for name, symbols in members.items():
            genes = [g for g in symbols if g in row_of]
            row, _ = set_test(genes, rho, nulls[map_name], row_of)
            row = dict(map=map_name, gene_set=name, **row)
            row["tested"] = len(genes) >= ISH_ANALYSIS["min_set_genes"]
            if not row["tested"]:
                row["p_spatial"] = np.nan
            part.append(row)
        part = pd.DataFrame(part)
        tested = part["tested"]
        part["q"] = np.nan
        part.loc[tested, "q"] = false_discovery_control(
            part.loc[tested, "p_spatial"], method="bh"
        )
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def contrast_test(
    first: list[str],
    second: list[str],
    rho: pd.Series,
    null: np.ndarray,
    row_of: dict[str, int],
    rng: np.random.Generator,
) -> tuple[dict, np.ndarray]:
    """The difference of two sets' median rho, against the surrogates and labels.

    Returns the row (both medians, the difference, its spatial p and band, its
    label-permutation p) and the surrogate differences.
    """
    a = rho[first].to_numpy(float)
    b = rho[second].to_numpy(float)
    observed = float(np.median(a) - np.median(b))
    null_a = np.median(null[[row_of[g] for g in first]], axis=0)
    null_b = np.median(null[[row_of[g] for g in second]], axis=0)
    spatial = null_a - null_b
    _, p_labels, _ = two_sample(a, b, rng)
    row = dict(
        n_first=len(first),
        n_second=len(second),
        median_first=float(np.median(a)),
        median_second=float(np.median(b)),
        difference=observed,
        p_spatial=spatial_p(observed, spatial),
        null_lo=float(np.percentile(spatial, BAND[0])),
        null_hi=float(np.percentile(spatial, BAND[1])),
        p_labels=p_labels,
    )
    return row, spatial


def contrast_tests(
    members: dict[str, list[str]],
    rho: pd.Series,
    null: np.ndarray,
    null_genes: list[str],
    seed: int = 0,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """contrasts.csv: each contrast named in advance, and its surrogate differences.

    Each contrast's label permutation has a generator of its own, so its p does not
    depend on which other contrasts ran.
    """
    row_of = {g: i for i, g in enumerate(null_genes)}
    seeds = np.random.SeedSequence(seed).spawn(len(CONTRASTS))
    rows, nulls = [], {}
    for name, child in zip(CONTRASTS, seeds):
        first, second, both = contrast_sides(members, name)
        first = [g for g in first if g in row_of]
        second = [g for g in second if g in row_of]
        rng = np.random.default_rng(child)
        row, nulls[name] = contrast_test(first, second, rho, null, row_of, rng)
        rows.append(dict(contrast=name, n_left_out=len(both), **row))
    return pd.DataFrame(rows), nulls


# ===== Localisation against matched controls =====


def subunit_composite(
    subunits: list[str], profiles: dict[str, dict[str, float]], structures: list[str]
) -> tuple[list[str], np.ndarray]:
    """The structures every subunit has, and the mean of the subunits' ranks there."""
    common = [s for s in structures if all(s in profiles[g] for g in subunits)]
    ranks = [rankdata([profiles[g][s] for s in common]) for g in subunits]
    return common, np.mean(ranks, axis=0)


def partial_rows(
    genes: dict[str, str],
    common: list[str],
    map_values: np.ndarray,
    composite: np.ndarray,
    profiles: dict[str, dict[str, float]],
    reliability: dict[str, float],
    level: dict[str, float],
) -> pd.DataFrame:
    """Plain and partial rho of each gene with the map, the composite removed.

    `genes` maps each gene to its side (localisation or control) and `map_values`
    holds the map on the `common` structures; a gene with fewer than
    ish_panel_test.min_structures of them has no row.
    """
    rows = []
    for gene, side in genes.items():
        positions = [i for i, s in enumerate(common) if s in profiles[gene]]
        if len(positions) < ISH_PANEL_TEST["min_structures"]:
            continue
        values = np.array([profiles[gene][common[i]] for i in positions])
        x = map_values[positions]
        rows.append(
            dict(
                symbol=gene,
                side=side,
                n_structures=len(positions),
                rho=float(spearmanr(x, values).statistic),
                rho_partial=partial(x, values, [composite[positions]]),
                reliability=reliability.get(gene, np.nan),
                median_energy=level.get(gene, np.nan),
            )
        )
    return pd.DataFrame(rows)


def matched_controls(table: pd.DataFrame, level: dict[str, float]) -> dict[str, str]:
    """{localisation gene: the unused control closest in log median energy}.

    ish.panel_test.greedy_match: the localisation genes in order of expression,
    highest first, each taking the nearest control left.
    """
    loc = list(table.loc[table["side"] == "localisation", "symbol"])
    ctrl = list(table.loc[table["side"] == "control", "symbol"])
    return greedy_match(loc, ctrl, level)


def residual_ranks(values: np.ndarray, design: np.ndarray) -> np.ndarray:
    """Ranks along the last axis minus their least-squares fit on `design`."""
    r = rankdata(values, axis=-1).astype(float)
    hat = design @ np.linalg.pinv(design)
    return r - r @ hat.T


def null_partial(surr: np.ndarray, gene: np.ndarray, composite: np.ndarray) -> np.ndarray:
    """The partial rho of the gene with every surrogate (a row), the composite removed.

    As ish.panel_test.partial: ranks, the composite's ranks and a constant
    regressed out of both, Pearson on the residuals.
    """
    design = np.column_stack([rankdata(composite), np.ones(len(composite))])
    a = residual_ranks(surr, design)
    b = residual_ranks(gene, design)
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean()
    return (a @ b) / np.sqrt(np.square(a).sum(axis=1) * np.square(b).sum())


def localisation_tests(
    table: pd.DataFrame, pairs: dict[str, str], seed: int = 0
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """localisation_summary.csv: the tests of ish.panel_test on the new inputs.

    `table` is partial_rows' table, `pairs` each localisation gene's matched
    control. Rows: against the matched controls (the test), against all controls,
    before partialling, on the genes with a reproducible map only, and the positive
    control. Each test has a generator of its own. Returns the table and each
    test's label-permutation null.
    """
    seeds = np.random.SeedSequence(seed).spawn(5)
    by = table.set_index("symbol")
    loc = list(table.loc[table["side"] == "localisation", "symbol"])
    ctrl = list(table.loc[table["side"] == "control", "symbol"])
    matched = sorted(set(pairs.values()))
    min_rel = ISH_PANEL_TEST["min_reliability"]
    good_loc = [g for g in loc if by.loc[g, "reliability"] >= min_rel]
    good_ctrl = [g for g in matched if by.loc[g, "reliability"] >= min_rel]
    have = [g for g in ctrl if np.isfinite(by.loc[g, "reliability"])]
    cut = float(np.median(by.loc[have, "reliability"]))
    high = [g for g in have if by.loc[g, "reliability"] >= cut]
    low = [g for g in have if by.loc[g, "reliability"] < cut]
    tests = (
        ("matched controls", loc, matched, "rho_partial", False),
        ("all controls", loc, ctrl, "rho_partial", False),
        ("before partialling", loc, matched, "rho", False),
        (f"reliability {min_rel} or more", good_loc, good_ctrl, "rho_partial", False),
        ("positive control", high, low, "rho", True),
    )
    rows, nulls = [], {}
    for (name, first, second, column, absolute), child in zip(tests, seeds):
        a = by.loc[first, column].to_numpy(float)
        b = by.loc[second, column].to_numpy(float)
        if absolute:
            a, b = np.abs(a), np.abs(b)
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        difference, p, null = two_sample(a, b, np.random.default_rng(child))
        nulls[name] = null
        rows.append(
            dict(
                test=name,
                statistic=f"|{column}|" if absolute else column,
                n_first=len(a),
                n_second=len(b),
                median_first=float(np.median(a)),
                median_second=float(np.median(b)),
                difference=difference,
                p_labels=p,
                detectable=float(np.percentile(np.abs(null), 95)),
                reliability_cut=cut if absolute else np.nan,
            )
        )
    return pd.DataFrame(rows), nulls


def matched_spatial(
    table: pd.DataFrame,
    pairs: dict[str, str],
    surr: np.ndarray,
    common: list[str],
    composite: np.ndarray,
    profiles: dict[str, dict[str, float]],
) -> tuple[float, np.ndarray]:
    """The matched difference against the surrogates, each one in place of the map.

    `surr` holds the surrogates on the `common` structures, one per row. Returns
    the two-sided p and the surrogate differences.
    """
    loc = list(table.loc[table["side"] == "localisation", "symbol"])
    matched = sorted(set(pairs.values()))
    partials = {}
    for gene in loc + matched:
        positions = [i for i, s in enumerate(common) if s in profiles[gene]]
        values = np.array([profiles[gene][common[i]] for i in positions])
        partials[gene] = null_partial(surr[:, positions], values, composite[positions])
    first = np.median([partials[g] for g in loc], axis=0)
    second = np.median([partials[g] for g in matched], axis=0)
    by = table.set_index("symbol")
    observed = float(
        np.median(by.loc[loc, "rho_partial"]) - np.median(by.loc[matched, "rho_partial"])
    )
    null = first - second
    return spatial_p(observed, null), null


def load_tables() -> dict[str, pd.DataFrame]:
    """The tables run_ish_gene_sets wrote, by name."""
    paths = dict(
        members=MEMBERS,
        sets=SET_TESTS,
        contrasts=CONTRAST_TESTS,
        localisation=LOCALISATION,
        summary=LOCALISATION_SUMMARY,
    )
    for path in paths.values():
        if not path.exists():
            raise FileNotFoundError(f"{path} not found: run run_ish_gene_sets.py first")
    return {name: pd.read_csv(path) for name, path in paths.items()}
