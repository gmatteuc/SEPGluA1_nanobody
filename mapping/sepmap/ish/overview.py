"""The ISH line gathered: the numbers for the text, April's headline, the figure index.

Every run script of the ISH line writes the numbers it computed to
tables/numbers_<step>.csv. They are gathered here into one table with the step each
comes from, so docs/ISH_ANALYSIS.md quotes every number from one file:

    numbers_for_the_text.csv    step, name, value, what
    numbers_for_the_text.txt    the same as aligned lines, for reading

What is left of April's headline is measured here too (the appendix figure). P9
grouped its genes by the categories of gene_targets.csv, written while the genes were
being chosen, split "auxiliary" by hand into Aux forebrain (Cacng8, Cacng3, Cnih2,
Cnih3 and Grm5, as written into the P9 script) and Aux other, and tested the ten
groups with a one-way ANOVA, p = 0.032:

    April's rho     the eroded Spearman of the frozen gene_panel_summary.csv of
                    adult_matlab/run_compare_with_allen_ish.m (22 April), the
                    statistic of the figure Sami saw
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

Last, the index of the guided figures (figures/README.md) and the rows of the
overview figure: for each figure the question, what to look at and what to take
from it, with the numbers of this run. The figures' order is the argument's: the
question, the inputs, how much of the map receptor mRNA and synaptic density leave,
what the leftover looks like through the genes, the controls and the limits.

Run by run_ish_overview.py.
"""

import csv
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, false_discovery_control, spearmanr

from sepmap.config import DATA, SETTINGS
from sepmap.ish.gene_sets import set_test
from sepmap.ish.plotting import FIGURES, QUESTIONS, figure_file, figure_ref
from sepmap.structures import TABLES

# the BH level the figures and the text count genes and sets at
ISH_ANALYSIS = SETTINGS["ish_analysis"]

NUMBERS = TABLES / "numbers_for_the_text.csv"
NUMBERS_TXT = TABLES / "numbers_for_the_text.txt"
HEADLINE = TABLES / "april_headline.csv"
HEADLINE_GROUPS = TABLES / "april_groups.csv"
HEADLINE_ANOVA = TABLES / "april_anova.csv"
INDEX = DATA / "adult_v2" / "ish_analysis" / "figures" / "README.md"
BEYOND = DATA / "adult_v2" / "ish_analysis" / "beyond"

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


# ===== The numbers for the text =====


def gather_numbers(tables: Path | None = None) -> pd.DataFrame:
    """Every numbers_<step>.csv of the run in one table: step, name, value, what.

    The values stay as the steps wrote them (text), so nothing is rounded twice.
    """
    if tables is None:
        tables = TABLES
    parts = []
    for path in sorted(tables.glob("numbers_*.csv")):
        if path == NUMBERS:
            continue
        part = pd.read_csv(path, dtype=str, keep_default_na=False)
        part.insert(0, "step", path.stem.removeprefix("numbers_"))
        parts.append(part)
    if not parts:
        raise FileNotFoundError(f"no numbers_*.csv in {tables}: run steps 13 to 26 first")
    return pd.concat(parts, ignore_index=True)


def lookup(numbers: pd.DataFrame) -> dict[str, str]:
    """{'step.name': value}; a name written twice by one step keeps its first value."""
    out = {}
    for step, name, value in zip(numbers["step"], numbers["name"], numbers["value"]):
        out.setdefault(f"{step}.{name}", value)
    return out


def text_lines(numbers: pd.DataFrame) -> list[str]:
    """numbers_for_the_text.txt: one aligned line per number, step by step."""
    width = max(len(f"{s}.{n}") for s, n in zip(numbers["step"], numbers["name"]))
    lines = []
    for step, part in numbers.groupby("step", sort=False):
        lines.append(f"[{step}]")
        for _, r in part.iterrows():
            key = f"{step}.{r['name']}"
            lines.append(f"  {key:<{width}}  {r['value']:<12}  {r['what']}")
        lines.append("")
    return lines


# ===== April's headline =====


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
        ("April as Sami saw it (hand split)", "april", h["rho_april"], by_hand),
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
    thousands of columns at once.
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
    return (between / (k - 1)) / (within / (n - k))


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
        row, _ = set_test(genes, rho, null, row_of)
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
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


def leftover_extremes(n_each: int = 4) -> pd.DataFrame:
    """The structures furthest above and below prediction (analysis 4), for the text."""
    table = pd.read_csv(BEYOND / "regression_table.csv")
    table = table.sort_values("residual", ascending=False)
    rows = [
        ("leftover_above", " ".join(table["acronym"].head(n_each)), "most above"),
        ("leftover_below", " ".join(table["acronym"].tail(n_each)[::-1]), "most below"),
        ("leftover_top_ranks", round(table["residual"].iloc[0], 1), "largest, ranks"),
    ]
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


def ranking_extras(
    ranking: pd.DataFrame, per_adult: pd.DataFrame, summary: pd.DataFrame
) -> pd.DataFrame:
    """Numbers of the gene by gene figures that their steps do not write.

    `ranking`, `per_adult` and `summary` are gene_ranking.csv,
    gene_ranking_per_adult.csv and robustness_summary.csv.
    """
    nano = ranking[ranking["map"] == "nano"].set_index("symbol")
    auto = ranking[ranking["map"] == "auto"].set_index("symbol")
    q = ISH_ANALYSIS["q"]
    agreement = spearmanr(nano["rho"], auto["rho"].reindex(nano.index)).statistic
    passing = nano[nano["q_all"] < q].sort_values("rho", ascending=False)
    auto_top = auto[auto["q_all"] < q].sort_values("rho", ascending=False)
    rows = [
        ("nano_auto_agreement", round(agreement, 3), "nano and auto gene orders"),
        ("nano_pass_genes", " ".join(passing.index), "genes past the nano null"),
        ("auto_top_genes", " ".join(auto_top.index[:5]), "highest past the auto null"),
    ]
    for gene in ("Cacng8", "Gria1", "Aqp4", "Dlg2"):
        mine = per_adult[per_adult["symbol"] == gene]
        above = int((mine["rho_nano"] > mine["rho_auto"]).sum())
        rows += [
            (f"n_structures_{gene}", int(nano.loc[gene, "n_structures"]), gene),
            (f"rho_{gene}", round(nano.loc[gene, "rho"], 3), gene),
            (f"rank_p9_{gene}", int(nano.loc[gene, "rank_p9"]), gene),
            (f"adults_nano_above_auto_{gene}", above, f"of {len(mine)} adults"),
        ]
    pair = per_adult[per_adult["symbol"].isin(["Cacng8", "Gria1"])]
    pair = pair.pivot(index="mouse", columns="symbol", values="rho_nano")
    leads = int((pair["Cacng8"] > pair["Gria1"]).sum())
    rows.append(("adults_cacng8_above_gria1", leads, f"of {len(pair)} adults"))
    others = summary[summary["variant"] != "primary"]
    rows += [
        ("robustness_variants", len(others), "variants of the primary ranking"),
        ("robustness_agreement_min", round(others["agreement_p9"].min(), 3), "lowest"),
        ("robustness_gria1_rank_min", int(summary["rank_p9_Gria1"].min()), "Gria1"),
        ("robustness_gria1_rank_max", int(summary["rank_p9_Gria1"].max()), "Gria1"),
        ("robustness_cacng8_rank_max", int(summary["rank_p9_Cacng8"].max()), "Cacng8"),
        ("robustness_gap_min", round(summary["gap"].min(), 3), "the gap"),
        ("robustness_gap_max", round(summary["gap"].max(), 3), "the gap"),
    ]
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


def gene_extras(tables: Path | None = None) -> pd.DataFrame:
    """Lists of genes the text names that their steps write only as tables.

    The genes past the within-division null (within_division.csv), the genes past
    the nano null counted by gene set (gene_ranking.csv), and the localisation
    genes with the highest partial rho against the GO control pool
    (localisation_test.csv).
    """
    if tables is None:
        tables = TABLES
    q = ISH_ANALYSIS["q"]
    within = pd.read_csv(tables / "within_division.csv")
    past = within[within["q_all_spatial"] < q].sort_values("rho_within", ascending=False)
    ranking = pd.read_csv(tables / "gene_ranking.csv", keep_default_na=False)
    nano = ranking[(ranking["map"] == "nano") & (ranking["q_all"].astype(float) < q)]
    first_set = nano["gene_sets"].str.split(";").str[0].replace("", "none")
    counts = first_set.value_counts()
    local = pd.read_csv(tables / "localisation_test.csv")
    local = local[(local["side"] == "localisation") & local["pool"].str.contains("GO")]
    top = local.sort_values("rho_partial", ascending=False)["symbol"].head(5)
    rows = [
        ("within_pass_genes", " ".join(past["symbol"]), "past the within null"),
        ("localisation_top_partial", " ".join(top), "highest partial rho"),
    ]
    for name, count in counts.items():
        key = name.replace(" ", "_")
        rows.append((f"nano_pass_set_{key}", int(count), f"genes past the null: {name}"))
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


# ===== The index of the figures =====


def num(n: dict[str, str], key: str) -> float:
    """One number of the run, by 'step.name'."""
    return float(n[key])


def pct(n: dict[str, str], key: str) -> str:
    """A share of the run as a whole percentage: '27%'."""
    return f"{num(n, key):.0%}"


def ordinal(value: str | float) -> str:
    """A rank as words write it: '1st', '2nd', '11th', '17th'."""
    k = int(float(value))
    if k % 100 in (11, 12, 13):
        return f"{k}th"
    return f"{k}{ {1: 'st', 2: 'nd', 3: 'rd'}.get(k % 10, 'th') }"


def inside_or_past(p: float, q: float = 0.05) -> str:
    """'past its null' or 'inside its null', by the p and the level."""
    if p < q:
        return "past its null"
    return "inside its null"


def input_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the figures of the inputs."""
    return {
        "overview": (
            "The argument in the box at the top; then one row per step of the walk: "
            "the question, the figures that answer it, the numbers of this run, and "
            "what stays open.",
            "Which figure answers which question, and where the argument stands.",
        ),
        "structures": (
            "A, the rule as a funnel; C, the adult map on the declared structures; D, "
            "how far each adult's zero and spread of zref move with the declared "
            "reference.",
            f"{int(num(n, 'structure_set.structures_declared'))} structures enter every "
            "comparison: grey matter measured in all ten adults, fixed by rule before "
            "any correlation. The declared reference moves each adult's zero by "
            f"{num(n, 'structure_set.zref_zero_shift_min'):.3f} to "
            f"{num(n, 'structure_set.zref_zero_shift_max'):.3f} and keeps the cohort "
            f"map's order (rho {num(n, 'structure_set.cohort_map_agreement'):.5f}).",
        ),
        "genes": (
            "B, how many Allen experiments each gene has; C, the sections set missing "
            "in P9's experiments (red), and those kept as true absence (hatched); D, "
            "how well two Allen experiments of one gene agree.",
            f"{n['gene_table.genes_listed']} genes and "
            f"{n['gene_table.experiments_used']} usable Allen experiments. Section QC "
            f"set {n['section_qc.sections_flagged']} failed sections missing in "
            f"{n['section_qc.experiments_flagged']} experiments, none filled in. Two "
            "experiments of the same gene agree at a median rho of "
            f"{num(n, 'gene_table.reliability_median'):.2f} "
            f"({n['gene_table.genes_reliability']} genes measured more than once), "
            "which caps any gene's rho with the nano map.",
        ),
    }


def beyond_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the figures of the leftover."""
    above = ", ".join(n["overview.leftover_above"].split())
    below = ", ".join(n["overview.leftover_below"].split())
    return {
        "beyond_budget": (
            "A to C, the nano map against Gria1 mRNA, against synaptic density and "
            "against the whole model; D, the variance budget, where red is what is "
            "left and the hatching the calibration floor; E, the leftover beside maps "
            "whose answer is known; F, one half of the cohort's leftover against the "
            "other's.",
            f"Gria1 mRNA alone predicts {pct(n, 'beyond.share_gria1')} of the map's "
            f"reproducible pattern, the four subunits {pct(n, 'beyond.share_abundance')}"
            f", synaptic density alone {pct(n, 'beyond.share_density')}; together, "
            f"with autofluorescence, {pct(n, 'beyond.share_model')}. "
            f"{pct(n, 'beyond.left')} is left ({pct(n, 'beyond.left_lo')} to "
            f"{pct(n, 'beyond.left_hi')} over resampled structures), and one half of "
            "the cohort's leftover agrees with the other's at "
            f"{num(n, 'beyond.replication_leftover'):.2f}. A map made only of receptor "
            f"mRNA and density leaves {pct(n, 'beyond.floor')} "
            f"({pct(n, 'beyond.floor_lo')} to {pct(n, 'beyond.floor_hi')}); a map "
            f"that is one Allen Gria1 experiment leaves {pct(n, 'beyond.gria1_map_left')}"
            f" ({pct(n, 'beyond.gria1_map_left_lo')} to "
            f"{pct(n, 'beyond.gria1_map_left_hi')}). The whole gene table, with no "
            f"hypothesis, leaves {1 - num(n, 'beyond.control_f_share'):.0%}.",
        ),
        "beyond_where": (
            "A, the map, the prediction and the leftover on three planes (red above "
            "prediction, blue below); B, the structures furthest from prediction; C "
            "and D, every gene and gene set against the leftover, each against the "
            "leftover's own surrogates.",
            f"Most above prediction: {above}; most below: {below}. "
            f"{n['beyond.leftover_genes_pass']} of {n['beyond.leftover_genes']} genes "
            "follow the leftover beyond its null; the closest is "
            f"{n['beyond.leftover_top_gene']} "
            f"({num(n, 'beyond.leftover_top_rho'):+.2f}, p "
            f"{num(n, 'beyond.leftover_top_p'):.2f} before correction), Cacng8 "
            f"{num(n, 'beyond.leftover_rho_Cacng8'):+.2f} (p "
            f"{num(n, 'beyond.leftover_p_Cacng8'):.2f}).",
        ),
    }


def gene_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the gene by gene figures."""
    gap_p = num(n, "gene_ranking.gap_p")
    return {
        "one_comparison": (
            "Each row reduces a map to one rank per declared structure (middle) and "
            "puts the gene's ranks against the map's (right); the dots are coloured "
            "by group of divisions.",
            f"Cacng8 {num(n, 'overview.rho_Cacng8'):+.2f} on "
            f"{n['overview.n_structures_Cacng8']} structures, Gria1 "
            f"{num(n, 'overview.rho_Gria1'):+.2f} on "
            f"{n['overview.n_structures_Gria1']}, the astrocyte gene Aqp4 "
            f"{num(n, 'overview.rho_Aqp4'):+.2f}. Much of a whole-brain rho is cortex "
            "and hippocampus against thalamus: the null and the test inside "
            "divisions come next.",
        ),
        "spatial_null": (
            "A, how widely unrelated smooth maps correlate; B, the surrogates' "
            "variogram on the map's; D, how often the ordinary and the spatial p fall "
            "below 0.05 on maps with no relation.",
            "Unrelated smooth maps correlate with an SD of rho "
            f"{num(n, 'spatial_null.sd_rho_pair'):.2f}, against "
            f"{num(n, 'spatial_null.sd_rho_shuffled_pair'):.2f} when one is "
            "shuffled. Random maps tested against nano give an ordinary p below 0.05 "
            f"in {pct(n, 'spatial_null.fpr_ordinary_map')} of tests, a spatial p in "
            f"{num(n, 'spatial_null.fpr_spatial_map'):.1%} (5% expected). The "
            "surrogates match the map's variogram within a median "
            f"{pct(n, 'spatial_null.variogram_misfit_nano')} over the matched range.",
        ),
        "gene_ranking": (
            "A, each gene's bar against its pale null band (a bar past its band "
            "passes at 0.05), the yellow dot its rho with autofluorescence; B, the "
            "Cacng8 - Gria1 gap against its null; C, how many genes pass.",
            f"{n['gene_ranking.nano_pass_all']} of {n['gene_ranking.nano_genes']} "
            "genes pass the spatial null at BH q < 0.05 "
            f"({n['gene_ranking.nano_pass_p9']} of P9's "
            f"{n['gene_ranking.nano_p9_genes']}). Cacng8 is first "
            f"({num(n, 'gene_ranking.nano_rho_Cacng8'):+.2f}) and Gria1 "
            f"{ordinal(n['overview.rank_p9_Gria1'])} of P9's genes "
            f"({num(n, 'gene_ranking.nano_rho_Gria1'):+.2f}), both past the null. "
            f"Cacng8 leads Gria1 by {num(n, 'gene_ranking.gap'):+.2f} "
            f"({num(n, 'gene_ranking.gap_boot_lo'):+.2f} to "
            f"{num(n, 'gene_ranking.gap_boot_hi'):+.2f} over resampled adults), and "
            f"the lead is {inside_or_past(gap_p)} (spatial p {gap_p:.2f}).",
        ),
        "gene_sets": (
            "A, one column per set fixed in advance: its median (black bar) against "
            "the band where the median of the same genes falls with 95% of the "
            "surrogates; B, the two contrasts named in advance.",
            "Localisation genes (transport, anchoring, auxiliary subunits; "
            f"{n['gene_sets.set_size_localisation']}) have a median rho of "
            f"{num(n, 'gene_sets.set_median_nano_localisation'):+.2f} (spatial p "
            f"{num(n, 'gene_sets.set_p_nano_localisation'):.3f}), other postsynaptic "
            f"genes {num(n, 'gene_sets.set_median_nano_other postsynaptic'):+.2f} (p "
            f"{num(n, 'gene_sets.set_p_nano_other postsynaptic'):.3f}), presynaptic "
            f"{num(n, 'gene_sets.set_median_nano_presynaptic'):+.2f}, glia "
            f"{num(n, 'gene_sets.set_median_nano_glia'):+.2f}. Postsynaptic genes "
            "sit above presynaptic ones by "
            f"{num(n, 'gene_sets.postsynaptic_against_presynaptic_difference'):+.2f} "
            "(spatial p "
            f"{num(n, 'gene_sets.postsynaptic_against_presynaptic_p_spatial'):.3f}).",
        ),
        "localisation": (
            "A, each gene's rho with the map once the subunits are removed, "
            "localisation genes against controls matched on expression; B, its null; "
            "C, the positive control, which must come out for the test to say "
            "anything.",
            "Localisation against matched controls: "
            f"{num(n, 'gene_sets.localisation_matched_controls_difference'):+.3f}, p "
            f"{num(n, 'gene_sets.localisation_matched_controls_p'):.2f} (spatial p "
            f"{num(n, 'gene_sets.localisation_p_spatial'):.2f}). The positive control "
            "with the GO control pool gives "
            f"{num(n, 'gene_sets.localisation_positive_control_difference'):+.3f} (p "
            f"{num(n, 'gene_sets.localisation_positive_control_p'):.2f}); with the "
            "control pool of 5 October it gives "
            f"{num(n, 'gene_sets.localisation_5_Octobers_controls_positive_control_difference'):+.3f}"  # noqa: E501
            " (p "
            f"{num(n, 'gene_sets.localisation_5_Octobers_controls_positive_control_p'):.3f}"  # noqa: E501
            "), and localisation against those controls "
            f"{num(n, 'gene_sets.localisation_5_Octobers_controls_matched_controls_difference'):+.3f}"  # noqa: E501
            " (p "
            f"{num(n, 'gene_sets.localisation_5_Octobers_controls_matched_controls_p'):.2f}"  # noqa: E501
            ").",
        ),
        "between_within": (
            "B, whole-brain rho against rho with a map that knows only each "
            "structure's division; C, whole-brain rho against the mean rho inside "
            "divisions (filled: past the within null); D, five genes division by "
            "division; E, why the surrogates and not a shuffle inside divisions.",
            "A map that knows only the division orders the genes as the real map "
            f"does ({num(n, 'divisions.agreement_division_only_all'):.2f}). Inside "
            "divisions the agreement falls to "
            f"{num(n, 'divisions.agreement_within_all'):.2f} and the median rho from "
            f"{num(n, 'divisions.median_rho'):+.2f} to "
            f"{num(n, 'divisions.median_rho_within'):+.2f}; "
            f"{n['divisions.pass_all_spatial']} genes pass the within null, Cacng8 "
            f"({num(n, 'divisions.rho_within_Cacng8'):+.2f}), Dlg2 "
            f"({num(n, 'divisions.rho_within_Dlg2'):+.2f}) and Gria1 "
            f"({num(n, 'divisions.rho_within_Gria1'):+.2f}) among them.",
        ),
    }


def control_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the controls, the limits and April."""
    auto_top = ", ".join(n["overview.auto_top_genes"].split())
    return {
        "autofluorescence": (
            "A, each gene's rho with the autofluorescence map against its rho with "
            "nano; C, adult by adult, nano and autofluorescence against Gria1 and "
            "Cacng8; D, the genes past each map's null.",
            "Autofluorescence of the same sections correlates with Gria1 at "
            f"{num(n, 'gene_ranking.auto_rho_Gria1'):+.2f} (p "
            f"{num(n, 'gene_ranking.auto_p_Gria1'):.2f}) and with Cacng8 at "
            f"{num(n, 'gene_ranking.auto_rho_Cacng8'):+.2f} (p "
            f"{num(n, 'gene_ranking.auto_p_Cacng8'):.2f}); adult by adult, nano's rho "
            f"with Gria1 is {num(n, 'gene_ranking.per_adult_nano_Gria1_min'):.2f} to "
            f"{num(n, 'gene_ranking.per_adult_nano_Gria1_max'):.2f}, "
            "autofluorescence's "
            f"{num(n, 'gene_ranking.per_adult_auto_Gria1_min'):+.2f} to "
            f"{num(n, 'gene_ranking.per_adult_auto_Gria1_max'):+.2f}. The tissue has "
            f"a pattern of its own: {n['gene_ranking.auto_pass_all']} genes pass its "
            f"null (highest {auto_top}), and the two gene orders agree at "
            f"{num(n, 'overview.nano_auto_agreement'):.2f}.",
        ),
        "robustness": (
            "A, how well each variant's gene order agrees with the primary's; B, "
            "where Cacng8 and Gria1 sit and the gap between them, variant by variant.",
            f"Over {n['overview.robustness_variants']} variants (statistic, borders, "
            "reading, inputs, structure set, and the route of 5 October) the gene "
            "order agrees with the primary at "
            f"{num(n, 'overview.robustness_agreement_min'):.2f} to 1.00; Cacng8's "
            f"lowest rank is {n['overview.robustness_cacng8_rank_max']}, Gria1 ranks "
            f"{n['overview.robustness_gria1_rank_min']} to "
            f"{n['overview.robustness_gria1_rank_max']}, and the gap stays between "
            f"{num(n, 'overview.robustness_gap_min'):+.2f} and "
            f"{num(n, 'overview.robustness_gap_max'):+.2f}.",
        ),
        "green_channel": (
            "A and C, which channel follows which, adult by adult; B, one adult's raw "
            "channels; D, how much each channel varies across structures.",
            "In every adult the green (SEP) channel follows autofluorescence (rho "
            f"{num(n, 'green_channel.rho_sep_auto_min'):.2f} to "
            f"{num(n, 'green_channel.rho_sep_auto_max'):.2f}) and varies across "
            f"structures about as little ({num(n, 'green_channel.range_sep_mean'):.2f}"
            f" against {num(n, 'green_channel.range_auto_mean'):.2f}, p90 - p10 of "
            f"log2; nano {num(n, 'green_channel.range_nano_mean'):.2f}). It does not "
            "report total receptor, so the surface fraction cannot be measured in "
            "these brains.",
        ),
        "april_headline": (
            "A, April's groups (grey) beside today's (orange), a line per gene; B, "
            "gene by gene; C, the ANOVA p under each choice it rests on; D, today's "
            "groups against the band of their median over the surrogates.",
            "April's gene order reproduces (rho "
            f"{num(n, 'overview.april_today_agreement'):.2f} over "
            f"{n['overview.april_today_genes']} genes). The group p was "
            f"{num(n, 'overview.anova_p_1'):.3f} with the hand split and "
            f"{num(n, 'overview.anova_p_2'):.2f} without; today "
            f"{num(n, 'overview.anova_p_4'):.3f} and {num(n, 'overview.anova_p_5'):.3f}"
            ". Against the map's surrogates, which keep co-expressed genes together, "
            f"the groups differ with p {num(n, 'overview.anova_p_spatial'):.2f}.",
        ),
    }


def meanings(n: dict[str, str]) -> dict[str, str]:
    """What each figure means, in one or two sentences; a verdict follows its numbers."""
    q = ISH_ANALYSIS["q"]
    left = num(n, "beyond.left")
    floor = num(n, "beyond.floor")
    one_gene = num(n, "beyond.gria1_map_left")
    than_one_gene = "no more than" if left <= one_gene else "more than"
    calibrated = 0.03 <= num(n, "spatial_null.fpr_spatial_map") <= 0.08
    gap_p = num(n, "gene_ranking.gap_p")
    pre_p = num(n, "gene_sets.postsynaptic_against_presynaptic_p_spatial")
    go_control = num(n, "gene_sets.localisation_positive_control_p")
    october = "gene_sets.localisation_5_Octobers_controls_"
    october_control = num(n, october + "positive_control_p")
    october_matched = num(n, october + "matched_controls_p")
    within = [
        g
        for g in ("Cacng8", "Dlg2", "Gria1")
        if num(n, f"divisions.p_within_spatial_{g}") < q
    ]
    sets_past = sum(
        1
        for key, value in n.items()
        if key.startswith("gene_sets.set_q_nano_") and value and float(value) < q
    )
    out = {
        "overview": "The argument and where it stands, on one page.",
        "structures": "Every comparison runs on structures every adult measures; the "
        "declared reference shifts each brain's zref without reordering the map.",
        "genes": "A gene measured once is only as good as one Allen mouse; a low rho "
        "of an unreliable gene says little.",
        "beyond_budget": f"{left:.0%} of the reproducible map is not predicted by "
        f"receptor mRNA or synaptic density, reproducibly across mice: "
        f"{left / floor:.1f} times what such a map leaves through Allen-to-Allen "
        f"mismatch, and {than_one_gene} what one Allen Gria1 experiment leaves.",
        "one_comparison": "A whole-brain rho mixes fine agreement with the contrast "
        "between divisions; the next figures separate them.",
    }
    if int(float(n["beyond.leftover_genes_pass"])) == 0:
        out["beyond_where"] = (
            "No gene's map and no gene set looks like the leftover beyond chance: "
            "where it sits is known, what it is these maps do not say."
        )
    else:
        out["beyond_where"] = "Some genes follow the leftover beyond chance."
    if calibrated:
        out["spatial_null"] = (
            "The spatial p keeps about 5% false positives where the ordinary one "
            "does not; every p of the gene analyses is a spatial p."
        )
    else:
        out["spatial_null"] = "The spatial p is not calibrated; no claim rests on it."
    if gap_p < q:
        lead = "and follows Cacng8 more closely than Gria1 beyond chance"
    else:
        lead = "but these maps cannot tell which of the two it follows more closely"
    out["gene_ranking"] = (
        f"The map follows a TARP's pattern and its own subunit's mRNA beyond "
        f"chance, {lead}: consistent with the surface-fraction reading, not "
        "evidence for it."
    )
    kind = "postsynaptic-like" if pre_p < q else "not clearly postsynaptic-like"
    passing = f"{sets_past} gene sets pass their null"
    if sets_past == 0:
        passing = "No gene set passes its null"
    localisation_q = n.get("gene_sets.set_q_nano_localisation", "")
    if localisation_q and float(localisation_q) < q:
        surface = "the genes that set surface receptors stand out beyond the null"
    else:
        surface = (
            "the genes that set surface receptors do not stand out from other "
            "postsynaptic genes"
        )
    out["gene_sets"] = (
        f"{passing}; the map is {kind}, as a glutamate receptor label should be, and "
        f"{surface}."
    )
    if go_control >= q and october_control < q and october_matched >= q:
        out["localisation"] = (
            "With the GO control pool the positive control fails, so that test alone "
            "says nothing; with the pool of 5 October it works, and the localisation "
            "genes are still no better than their matched controls."
        )
    else:
        out["localisation"] = (
            "Read the positive control first: the localisation difference counts only "
            "where it comes out."
        )
    if within:
        names = within[-1]
        if len(within) > 1:
            names = ", ".join(within[:-1]) + " and " + within[-1]
        out["between_within"] = (
            f"{names} follow the map inside divisions too, not only through the "
            "contrast between cortex and thalamus: the stronger claim."
        )
    else:
        out["between_within"] = "No detail gene follows the map inside divisions."
    auto_p = num(n, "gene_ranking.auto_p_Gria1")
    above = int(float(n["overview.adults_nano_above_auto_Gria1"]))
    if auto_p >= q and above == 10:
        out["autofluorescence"] = (
            "The ranking is the label's, not the tissue's: autofluorescence does not "
            "follow Gria1 beyond its null, falls below nano for Gria1 in every adult, "
            "and orders the genes its own way."
        )
    else:
        out["autofluorescence"] = "The tissue shares part of the label's ranking."
    out["robustness"] = (
        "Cacng8 first holds under every choice; Gria1's rank moves with them, so it "
        "is quoted with its null, not as a place."
    )
    if num(n, "green_channel.rho_sep_auto_min") > 0.5:
        out["green_channel"] = (
            "The green channel reports mostly the tissue in every adult: these "
            "brains cannot measure the surface fraction; a total-GluA1 stain would."
        )
    else:
        out["green_channel"] = "The green channel does not simply follow the tissue."
    if num(n, "overview.anova_p_spatial") >= q:
        out["april_headline"] = (
            "What reproduces from April is the gene order, Cacng8 first; the "
            "difference between its categories does not survive the null."
        )
    else:
        out["april_headline"] = "April's categories differ beyond the null."
    return out


# what each part of the walk is, for the index; the figures of each, in order
PARTS = (
    ("The question", ("overview",)),
    ("The inputs", ("structures", "genes")),
    (
        "Part 1: the map is not explained by receptor mRNA or synaptic density",
        ("beyond_budget", "beyond_where"),
    ),
    (
        "Part 2: what the leftover looks like, through the genes",
        (
            "one_comparison",
            "spatial_null",
            "gene_ranking",
            "gene_sets",
            "localisation",
            "between_within",
        ),
    ),
    ("Controls and limits", ("autofluorescence", "robustness", "green_channel")),
    ("Appendix", ("april_headline",)),
)

# the run script that draws each guided figure
DRAWN_BY = {
    "overview": "run_ish_overview.py",
    "structures": "run_structure_set.py",
    "genes": "run_ish_gene_table.py",
    "beyond_budget": "run_beyond_figures.py",
    "beyond_where": "run_beyond_figures.py",
    "one_comparison": "run_ish_gene_ranking.py",
    "spatial_null": "run_ish_gene_ranking.py",
    "gene_ranking": "run_ish_gene_ranking.py",
    "gene_sets": "run_ish_gene_sets.py",
    "localisation": "run_ish_gene_sets.py",
    "between_within": "run_ish_divisions.py",
    "autofluorescence": "run_ish_gene_ranking.py",
    "robustness": "run_ish_robustness.py",
    "green_channel": "run_sep_channel_check.py",
    "april_headline": "run_ish_overview.py",
}

INDEX_HEAD = """# The ISH analysis, figure by figure

Written by `run_ish_overview.py` from the tables of this run, so every number below
is this run's. The story, with what each result means and does not mean, is
`docs/ISH_ANALYSIS.md` in the code repository.

The figures follow one argument in two parts. Part 1: the adult nano map across
structures is not satisfactorily explained by Gria1 expression or by synaptic density
(receptor mRNA, synaptic markers, postsynaptic-density genes); a sizeable,
reproducible part is left over. Part 2: the genes that regulate surface AMPA
receptors (Cacng8, a TARP; trafficking and scaffolding genes) are tested as
corroboration of the reading the data support, the surface fraction of the receptor,
against abundance genes, unrelated genes and the spatial null. The surface fraction
is an interpretation, not a measurement: a total-GluA1 stain on the same brains would
measure it, and the green channel cannot ({limit}).
"""

INDEX_TAIL = """## Sheets

- `qc/00_flagged.png`: every experiment with a section set missing, one row each,
  sections as columns (red set missing, hatched kept as true absence). The sheet to
  review `mapping/ish_section_exceptions.csv` from.
- `qc/<gene>_<experiment>.png`: one sheet per Allen experiment: its median energy
  section by section along its own axis, with the local reference and the flags, and
  the orientation check at its middle section (`run_ish_section_qc.py --sheets`).
- `genes/<gene>.png`: Cacng8, Gria1, Grm5, Dlg2 and Aqp4: the nano and gene rank
  maps, the scatter of ranks with one fitted line per division, the whole-brain and
  within-division rho with their p (`run_ish_divisions.py --sheets`).
"""


def figure_index(n: dict[str, str]) -> str:
    """figures/README.md: the guided walk, each figure with its question and caption."""
    walk = input_walk(n) | beyond_walk(n) | gene_walk(n) | control_walk(n)
    meaning = meanings(n)
    lines = [INDEX_HEAD.format(limit=figure_ref("green_channel"))]
    for part, keys in PARTS:
        lines.append(f"## {part}\n")
        for key in keys:
            look, take = walk[key]
            name = figure_file(key)
            lines += [
                f"### {FIGURES[key]:02d}. {QUESTIONS[key]}\n",
                f"![{name}]({name})\n",
                f"**Look at.** {look}\n",
                f"**Take from it.** {take}\n",
                f"**What it means.** {meaning[key]}\n",
                f"Drawn by `{DRAWN_BY[key]}`; also as `{name.replace('.png', '.eps')}`."
                "\n",
            ]
    lines.append(INDEX_TAIL)
    return "\n".join(lines)


# ===== The overview figure =====


def argument_text(n: dict[str, str]) -> list[tuple[str, str]]:
    """The box at the top of the overview: the question, the two parts, the limit."""
    return [
        (
            "The question",
            "The adult nano map orders the brain's structures the same way in every "
            "half of the cohort (half against half, rho "
            f"{num(n, 'beyond.half_agreement'):.2f}). Is that order how much GluA1 "
            "mRNA a structure makes, or how many synapses it has?",
        ),
        (
            "Part 1",
            "Not satisfactorily. Gria1 mRNA alone predicts "
            f"{pct(n, 'beyond.share_gria1')} of the reproducible map; the four "
            "subunits, synaptic density and autofluorescence together "
            f"{pct(n, 'beyond.share_model')}. {pct(n, 'beyond.left')} is left "
            f"({pct(n, 'beyond.left_lo')} to {pct(n, 'beyond.left_hi')}), it "
            "replicates across halves of the cohort at "
            f"{num(n, 'beyond.replication_leftover'):.2f}, and a map made only of "
            f"receptor mRNA and density would leave {pct(n, 'beyond.floor')} "
            f"(figures {FIGURES['beyond_budget']:02d} and "
            f"{FIGURES['beyond_where']:02d}).",
        ),
        (
            "Part 2",
            "The map is therefore something else; the reading the data support is "
            "the surface fraction of the receptor, shaped by trafficking and "
            "scaffolding. If so, the genes that regulate surface AMPA receptors "
            "(Cacng8, trafficking and scaffolding genes) should track the map better "
            "than abundance genes or unrelated genes, beyond the spatial null "
            f"(figures {FIGURES['one_comparison']:02d} to "
            f"{FIGURES['between_within']:02d}).",
        ),
        (
            "The limit",
            "The surface fraction is an interpretation, not a measurement: a "
            "total-GluA1 stain on the same brains would measure it; the green channel "
            f"cannot (figure {FIGURES['green_channel']:02d}).",
        ),
    ]


def figs(*keys: str) -> str:
    """The numbers of some guided figures, as the overview lists them: '03, 04'."""
    return ", ".join(f"{FIGURES[k]:02d}" for k in keys)


def first_rows(n: dict[str, str]) -> list[dict]:
    """The overview's rows of the inputs and of part 1."""
    return [
        dict(
            part="inputs",
            title="The inputs",
            figures=figs("structures", "genes"),
            question="Which structures, which genes, and how good is one Allen map?",
            lines=[
                f"{int(num(n, 'structure_set.structures_declared'))} structures: grey "
                "matter measured in all 10 adults",
                f"{n['gene_table.genes_listed']} genes, "
                f"{n['gene_table.experiments_used']} Allen experiments; "
                f"{n['section_qc.sections_flagged']} failed sections set missing in "
                f"{n['section_qc.experiments_flagged']} experiments, none filled in",
                "two experiments of a gene agree at a median rho of "
                f"{num(n, 'gene_table.reliability_median'):.2f} "
                f"({n['gene_table.genes_reliability']} genes)",
            ],
            open="The true absences of the exceptions list are proposed, to review on "
            "qc/00_flagged.png. The adults' ages are not recorded (Allen: P56).",
        ),
        dict(
            part="part 1",
            title="Part 1: not explained by receptor mRNA or synaptic density",
            figures=figs("beyond_budget", "beyond_where"),
            question="How much of the map do receptor mRNA and synaptic density "
            "predict, and is what they leave real?",
            lines=[
                f"Gria1 alone {pct(n, 'beyond.share_gria1')}, Gria1-4 "
                f"{pct(n, 'beyond.share_abundance')}, density alone "
                f"{pct(n, 'beyond.share_density')} of the reproducible map",
                "all together, with autofluorescence, "
                f"{pct(n, 'beyond.share_model')}: {pct(n, 'beyond.left')} left "
                f"({pct(n, 'beyond.left_lo')} to {pct(n, 'beyond.left_hi')})",
                "the leftover replicates across halves of the cohort at "
                f"{num(n, 'beyond.replication_leftover'):.2f} (the map's reliability "
                f"alone implies {num(n, 'beyond.replication_implied'):.2f})",
                f"calibration floor {pct(n, 'beyond.floor')} "
                f"({pct(n, 'beyond.floor_lo')} to {pct(n, 'beyond.floor_hi')}); one "
                f"Allen Gria1 experiment as the map: {pct(n, 'beyond.gria1_map_left')}",
                "the whole gene table, no hypothesis: "
                f"{1 - num(n, 'beyond.control_f_share'):.0%} left; "
                f"{n['beyond.controls_passed']} of 7 controls pass",
            ],
            open="What the leftover is: a total-GluA1 stain separates the surface "
            "fraction from translation, turnover, subunit composition and nanobody "
            "access. A claim about one structure needs its own null.",
        ),
    ]


def gene_rows(n: dict[str, str]) -> list[dict]:
    """The overview's rows of part 2: gene by gene, kinds of genes, inside divisions."""
    gap_p = num(n, "gene_ranking.gap_p")
    october = "gene_sets.localisation_5_Octobers_controls_positive_control_p"
    return [
        dict(
            part="part 2",
            title="Part 2: gene by gene",
            figures=figs("one_comparison", "spatial_null", "gene_ranking"),
            question="Which genes' maps order the structures as the nano map does, "
            "beyond a map with its smoothness? Does Cacng8 lead Gria1?",
            lines=[
                f"{n['gene_ranking.nano_pass_all']} of {n['gene_ranking.nano_genes']} "
                "genes past the spatial null (BH q < 0.05); "
                f"{n['gene_ranking.nano_pass_p9']} of P9's "
                f"{n['gene_ranking.nano_p9_genes']}",
                f"Cacng8 first, {num(n, 'gene_ranking.nano_rho_Cacng8'):+.2f}; Gria1 "
                f"{num(n, 'gene_ranking.nano_rho_Gria1'):+.2f}, "
                f"{ordinal(n['overview.rank_p9_Gria1'])} of P9's genes",
                f"Cacng8 - Gria1 gap {num(n, 'gene_ranking.gap'):+.2f} "
                f"({num(n, 'gene_ranking.gap_boot_lo'):+.2f} to "
                f"{num(n, 'gene_ranking.gap_boot_hi'):+.2f} over adults): "
                f"{inside_or_past(gap_p)}, spatial p {gap_p:.2f}",
                "the spatial p is calibrated: "
                f"{num(n, 'spatial_null.fpr_spatial_map'):.1%} false positives at "
                f"0.05 (ordinary p: {pct(n, 'spatial_null.fpr_ordinary_map')})",
            ],
            open="A background panel of a few thousand Allen genes, to place the "
            "panel among genes nobody chose.",
        ),
        dict(
            part="part 2",
            title="Part 2: kinds of genes",
            figures=figs("gene_sets", "localisation"),
            question="Do the genes that set surface AMPA receptors (trafficking, "
            "anchoring, auxiliary subunits) match the map better than others?",
            lines=[
                "localisation genes "
                f"({n['gene_sets.set_size_localisation']}): median "
                f"{num(n, 'gene_sets.set_median_nano_localisation'):+.2f}, p "
                f"{num(n, 'gene_sets.set_p_nano_localisation'):.3f}; other "
                "postsynaptic "
                f"{num(n, 'gene_sets.set_median_nano_other postsynaptic'):+.2f}, p "
                f"{num(n, 'gene_sets.set_p_nano_other postsynaptic'):.3f}",
                "postsynaptic above presynaptic by "
                f"{num(n, 'gene_sets.postsynaptic_against_presynaptic_difference'):+.2f}"
                ", spatial p "
                f"{num(n, 'gene_sets.postsynaptic_against_presynaptic_p_spatial'):.3f}",
                "localisation against matched controls, abundance removed: "
                f"{num(n, 'gene_sets.localisation_matched_controls_difference'):+.3f}"
                f", p {num(n, 'gene_sets.localisation_matched_controls_p'):.2f}",
                "positive control: GO pool p "
                f"{num(n, 'gene_sets.localisation_positive_control_p'):.2f}; "
                f"5 October's pool p {num(n, october):.4f}",
            ],
            open="A positive control that comes out with the GO control pool.",
        ),
        dict(
            part="part 2",
            title="Part 2: inside divisions, and the leftover itself",
            figures=figs("between_within", "beyond_where"),
            question="Does a gene follow the map inside divisions, and does any gene "
            "follow the leftover?",
            lines=[
                "a division-only map orders the genes as the map does: "
                f"{num(n, 'divisions.agreement_division_only_all'):.2f}; inside "
                f"divisions {num(n, 'divisions.agreement_within_all'):.2f}",
                f"{n['divisions.pass_all_spatial']} genes past the within null: "
                f"Cacng8 {num(n, 'divisions.rho_within_Cacng8'):+.2f}, Dlg2 "
                f"{num(n, 'divisions.rho_within_Dlg2'):+.2f}, Gria1 "
                f"{num(n, 'divisions.rho_within_Gria1'):+.2f}",
                f"against the leftover: {n['beyond.leftover_genes_pass']} of "
                f"{n['beyond.leftover_genes']} genes past its null; Cacng8 "
                f"{num(n, 'beyond.leftover_rho_Cacng8'):+.2f}, localisation median "
                f"{num(n, 'beyond.leftover_set_localisation'):+.2f}",
            ],
            open="",
        ),
    ]


def last_rows(n: dict[str, str]) -> list[dict]:
    """The overview's rows of the controls, the limit and the appendix."""
    return [
        dict(
            part="controls",
            title="Controls: the tissue, and the choices made",
            figures=figs("autofluorescence", "robustness"),
            question="Is the ranking the label's or the tissue's? Does it hang on a "
            "choice of statistic, border, reading, input or structure set?",
            lines=[
                "autofluorescence with Gria1 "
                f"{num(n, 'gene_ranking.auto_rho_Gria1'):+.2f} (p "
                f"{num(n, 'gene_ranking.auto_p_Gria1'):.2f}); nano above it in "
                f"{n['overview.adults_nano_above_auto_Gria1']} of 10 adults",
                f"{n['gene_ranking.auto_pass_all']} genes past the autofluorescence "
                "null; the two gene orders agree at "
                f"{num(n, 'overview.nano_auto_agreement'):.2f}",
                f"{n['overview.robustness_variants']} variants: gene order agrees at "
                f"{num(n, 'overview.robustness_agreement_min'):.2f} to 1.00; Gria1 "
                f"{ordinal(n['overview.robustness_gria1_rank_min'])} to "
                f"{ordinal(n['overview.robustness_gria1_rank_max'])}",
            ],
            open="A no-primary or knockout control; DAPI through the same chain.",
        ),
        dict(
            part="limit",
            title="The limit: the green channel",
            figures=figs("green_channel"),
            question="Does the green channel report total receptor, so the surface "
            "fraction can be measured?",
            lines=[
                "SEP follows autofluorescence at "
                f"{num(n, 'green_channel.rho_sep_auto_min'):.2f} to "
                f"{num(n, 'green_channel.rho_sep_auto_max'):.2f} in every adult; nano "
                f"at {num(n, 'green_channel.rho_nano_auto_min'):.2f} to "
                f"{num(n, 'green_channel.rho_nano_auto_max'):.2f}",
                "spread across structures (log2 p90 - p10): nano "
                f"{num(n, 'green_channel.range_nano_mean'):.2f}, autofluorescence "
                f"{num(n, 'green_channel.range_auto_mean'):.2f}, SEP "
                f"{num(n, 'green_channel.range_sep_mean'):.2f}",
            ],
            open="A total-GluA1 stain or autoradiography on some of the same brains; "
            "the green channel's filter sets.",
        ),
        dict(
            part="appendix",
            title="Appendix: April's headline",
            figures=figs("april_headline"),
            question="What is left of the category violins and their p = 0.032?",
            lines=[
                "the gene order reproduces: rho "
                f"{num(n, 'overview.april_today_agreement'):.2f} over "
                f"{n['overview.april_today_genes']} genes",
                f"group ANOVA p {num(n, 'overview.anova_p_1'):.3f} with the hand "
                f"split, {num(n, 'overview.anova_p_2'):.2f} without; today "
                f"{num(n, 'overview.anova_p_4'):.3f}",
                "against the map's surrogates (co-expressed genes kept together): p "
                f"{num(n, 'overview.anova_p_spatial'):.2f}",
            ],
            open="",
        ),
    ]


def overview_rows(n: dict[str, str]) -> list[dict]:
    """The rows of the overview figure, in the order of the walk.

    Each row: part, title, figures, question, lines (the numbers of this run, one
    short line each) and open (fixed text).
    """
    return first_rows(n) + gene_rows(n) + last_rows(n)


# the A-items of docs/REFACTOR_COVERAGE.md that the ISH line touches, and their state
ITEMS = (
    ("A1", "one declared structure set and zref reference: built for the ISH line"),
    ("A2", "section QC, failed sections set missing: built; true absences proposed"),
    ("A3", "one adult profile; Spearman, Pearson, borders, structure sets: built"),
    ("A4", "which structures sit reliably above the median: open"),
    ("A5", "nano against autofluorescence, structure by structure: open"),
    ("A6", "agreement within divisions, per-gene sheets: built"),
    ("A7", "a spatial null for every rho, the gap and the sets: built"),
    ("A8", "every gene against the autofluorescence map: built"),
    ("A9", "the gene table and its documentation: built"),
    ("A10", "per-gene videos: optional, not built"),
)
