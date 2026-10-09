"""What is left of April's headline: P9's gene categories, then and now.

P9 grouped its genes by the categories of gene_targets.csv, written while the genes
were being chosen, split "auxiliary" by hand into Aux forebrain (Cacng8, Cacng3,
Cnih2, Cnih3 and Grm5, as written into the P9 script) and Aux other, and tested the
ten groups with a one-way ANOVA, p = 0.032. Measured here on today's inputs and
against the map's surrogates (the appendix figure, 16):

    April's rho     the eroded Spearman of the frozen gene_panel_summary.csv of
                    adult_matlab/run_compare_with_allen_ish.m (22 April), the
                    statistic of April's headline figure
    today's rho     ish.gene_ranking's, on the declared structures
    choices         the ANOVA with and without the hand split, without the two
                    stretched grids (Olig2, Calb2), and today on two other structure
                    sets and with P9's single experiment per gene (rows of
                    ish.robustness); each treats co-expressed genes as independent
                    draws, as April did
    the null        each group's median rho against the median of the same genes'
                    rho with every surrogate of the map (ish.gene_sets.set_test), and
                    the ANOVA's F against the F of every surrogate; co-expressed genes
                    meet the same surrogate, so they stay together in the null

Run by run_ish_overview.py.
"""

import csv

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, false_discovery_control, spearmanr

from sepmap.config import DATA, SETTINGS
from sepmap.ish.gene_sets import set_test
from sepmap.ish.numbers import numbers_frame
from sepmap.structures import TABLES

# the genes a group needs to be tested, the BH level
ISH_ANALYSIS = SETTINGS["ish_analysis"]

HEADLINE = TABLES / "april_headline.csv"
HEADLINE_GROUPS = TABLES / "april_groups.csv"
HEADLINE_ANOVA = TABLES / "april_anova.csv"

# the frozen table of P9 (22 April), the source of April's rho
P9_SUMMARY = (
    DATA
    / "comparisons"
    / "merged_naive_rws_vs_ish_summary_nosmooth"
    / "gene_panel_summary.csv"
)

# the hand split of P9's "auxiliary", as written into the P9 script
AUX_FOREBRAIN = ("Cacng8", "Cacng3", "Cnih2", "Cnih3", "Grm5")

# P9's ten groups in the order of its headline (sorted by their median there)
HEADLINE_ORDER = (
    ("aux_forebrain", "Aux forebrain"),
    ("plasticity", "Plasticity"),
    ("AMPAR_core", "AMPAR core"),
    ("scaffold", "Scaffold"),
    ("excitatory", "Excitatory"),
    ("trafficking", "Trafficking"),
    ("control_struct", "Ctrl struct"),
    ("control_inhib", "Ctrl inhib"),
    ("aux_other", "Aux other"),
    ("control_glia", "Ctrl glia"),
)

# the two grids P9 stretched onto the atlas, which sat in a box of their own
STRETCHED = ("Olig2", "Calb2")

# the robustness rows the ANOVA is run on: the structure sets and P9's experiments
ANOVA_ROWS = (
    ("every_structure", "today, every structure of the adult table"),
    ("p9_divisions", "today, P9's nine divisions, any number of adults"),
    ("p9_experiment", "today, P9's single experiment per gene"),
)


def load_p9_headline() -> pd.DataFrame:
    """P9's genes as its frozen table holds them: symbol, category, rho_april.

    rho_april is the eroded Spearman, the statistic of the headline; genes without
    one (no usable grid in April) are left out. Raises FileNotFoundError when the
    table is missing, rather than draw the appendix without April.
    """
    if not P9_SUMMARY.exists():
        raise FileNotFoundError(
            f"P9's gene table, the source of April's headline, is missing: {P9_SUMMARY}"
        )
    rows = []
    with open(P9_SUMMARY, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["r_spearman_ero"] in ("", "NaN"):
                continue
            rows.append(
                dict(
                    symbol=r["symbol"],
                    category=r["category"],
                    rho_april=float(r["r_spearman_ero"]),
                )
            )
    return pd.DataFrame(rows)


def headline_group(symbol: str, category: str) -> str:
    """The group of P9's headline: its category, with "auxiliary" split by hand."""
    if category != "auxiliary":
        return category
    if symbol in AUX_FOREBRAIN:
        return "aux_forebrain"
    return "aux_other"


def headline_table(p9: pd.DataFrame, ranking: pd.DataFrame) -> pd.DataFrame:
    """april_headline.csv: P9's genes with their group, April's rho and today's.

    One row per gene of P9's panel that either side measured; rho_today is NaN for
    a gene today's ranking lacks, rho_april for one April lacked.
    """
    nano = ranking[(ranking["map"] == "nano") & ranking["p9_gene"]].set_index("symbol")
    april = p9.set_index("symbol")
    rows = []
    for symbol in sorted(set(april.index) | set(nano.index)):
        if symbol in nano.index:
            category = nano.loc[symbol, "p9_category"]
        else:
            category = april.loc[symbol, "category"]
        rows.append(
            dict(
                symbol=symbol,
                category=category,
                group=headline_group(symbol, category),
                rho_april=april["rho_april"].get(symbol, np.nan),
                rho_today=nano["rho"].get(symbol, np.nan),
            )
        )
    return pd.DataFrame(rows)


def anova_p(rho: pd.Series, groups: pd.Series) -> tuple[float, int]:
    """The one-way ANOVA p across the groups, and the genes it used (rho not NaN)."""
    ok = rho.notna() & groups.notna()
    samples = [rho[ok & (groups == g)].to_numpy() for g in groups[ok].unique()]
    return float(f_oneway(*samples).pvalue), int(ok.sum())


def anova_table(headline: pd.DataFrame, robustness_rows: pd.DataFrame) -> pd.DataFrame:
    """april_anova.csv: the headline's ANOVA under each choice it depends on.

    `robustness_rows` is ranking_robustness.csv. Columns: label, kind (april or
    today), n_genes, p.
    """
    h = headline.set_index("symbol")
    by_hand = h["group"]
    plain = h["category"]
    no_stretched = h["rho_april"].drop(list(STRETCHED), errors="ignore")
    variants = [
        ("April's headline (hand split)", "april", h["rho_april"], by_hand),
        ("April without the hand split", "april", h["rho_april"], plain),
        ("April without Olig2 and Calb2 (stretched)", "april", no_stretched, by_hand),
        ("today, declared structures (hand split)", "today", h["rho_today"], by_hand),
        ("today, without the hand split", "today", h["rho_today"], plain),
    ]
    for variant, label in ANOVA_ROWS:
        mine = robustness_rows[robustness_rows["variant"] == variant]
        rho = mine.set_index("symbol")["rho"].reindex(h.index)
        variants.append((label, "today", rho, by_hand))
    rows = []
    for label, kind, rho, groups in variants:
        p, n = anova_p(rho, groups.reindex(rho.index))
        rows.append(dict(label=label, kind=kind, n_genes=n, p=p))
    return pd.DataFrame(rows)


def f_statistic(values: np.ndarray, group_of_row: np.ndarray) -> np.ndarray:
    """The one-way ANOVA F of every column of `values` (genes x columns) across groups.

    `group_of_row` gives each row's group; the same F as scipy's f_oneway, for
    thousands of columns at once. A column with no spread inside the groups has F
    infinite (NaN when it has none between them either).
    """
    groups = np.unique(group_of_row)
    n, k = values.shape[0], len(groups)
    grand = values.mean(axis=0)
    between = np.zeros(values.shape[1])
    within = np.zeros(values.shape[1])
    for g in groups:
        part = values[group_of_row == g]
        mean = part.mean(axis=0)
        between += len(part) * (mean - grand) ** 2
        within += ((part - mean) ** 2).sum(axis=0)
    between, within = between / (k - 1), within / (n - k)
    out = np.full(values.shape[1], np.nan)
    np.divide(between, within, out=out, where=within > 0)
    out[(within == 0) & (between > 0)] = np.inf
    return out


def headline_null(
    headline: pd.DataFrame, null: np.ndarray, null_genes: list[str]
) -> tuple[pd.DataFrame, dict]:
    """April's groups on today's rho, against the map's surrogates.

    `null` is genes x surrogates of rho with the nano map's surrogates, rows in the
    order of `null_genes` (ish.gene_ranking). Returns april_groups.csv (per group:
    genes, median rho, the band of its median over the surrogates, spatial p and
    BH q over the groups tested; as for the gene sets, a group with fewer than
    ish_analysis.min_set_genes genes is drawn, not tested) and the ANOVA's F with
    its spatial p: the share of surrogates whose F across the same groups is at
    least the observed one, counted once (one-sided, since any difference between
    groups raises F).
    """
    row_of = {g: i for i, g in enumerate(null_genes)}
    today = headline[headline["rho_today"].notna() & headline["symbol"].isin(row_of)]
    rho = today.set_index("symbol")["rho_today"]
    rows = []
    for key, label in HEADLINE_ORDER:
        genes = sorted(today.loc[today["group"] == key, "symbol"])
        if genes:
            row, _ = set_test(genes, rho, null, row_of)
        else:
            row = dict(n_genes=0, median_rho=np.nan, null_lo=np.nan, null_hi=np.nan)
            row["p_spatial"] = np.nan
        rows.append(dict(group=key, label=label, **row))
    groups = pd.DataFrame(rows)
    groups["tested"] = groups["n_genes"] >= ISH_ANALYSIS["min_set_genes"]
    groups.loc[~groups["tested"], "p_spatial"] = np.nan
    groups["q"] = np.nan
    tested = groups["tested"]
    groups.loc[tested, "q"] = false_discovery_control(
        groups.loc[tested, "p_spatial"], method="bh"
    )

    order = [row_of[g] for g in today["symbol"]]
    group_of_row = today["group"].to_numpy()
    f_observed = float(f_statistic(rho.to_numpy()[:, None], group_of_row)[0])
    f_null = f_statistic(null[order], group_of_row)
    p = (int(np.sum(f_null >= f_observed)) + 1) / (len(f_null) + 1)
    return groups, dict(f=f_observed, p_spatial=p, n_genes=len(today), f_null=f_null)


def headline_numbers(
    headline: pd.DataFrame, anova: pd.DataFrame, groups: pd.DataFrame, f: dict
) -> pd.DataFrame:
    """The numbers of the appendix for the text, one row each (step 'overview')."""
    both = headline.dropna(subset=["rho_april", "rho_today"])
    agreement = spearmanr(both["rho_april"], both["rho_today"]).statistic
    rows = [
        ("april_genes", int(headline["rho_april"].notna().sum()), "P9's genes in April"),
        ("today_genes", int(headline["rho_today"].notna().sum()), "P9's genes today"),
        ("april_today_genes", len(both), "genes with both"),
        ("april_today_agreement", round(agreement, 3), "Spearman over those genes"),
    ]
    for i, r in anova.iterrows():
        rows.append((f"anova_p_{i + 1}", round(r["p"], 4), r["label"]))
    rows += [
        ("anova_f_today", round(f["f"], 3), "today's F across April's groups"),
        ("anova_p_spatial", round(f["p_spatial"], 4), "that F against the surrogates"),
    ]
    tested = groups[groups["tested"]]
    q = ISH_ANALYSIS["q"]
    rows += [
        ("groups_tested", len(tested), "April's groups with enough genes to test"),
        ("groups_past_q", int((tested["q"] < q).sum()), f"of them past BH q < {q}"),
    ]
    return numbers_frame(rows)
