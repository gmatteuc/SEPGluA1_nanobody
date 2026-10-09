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
question, the inputs, how much of the map Gria1 and synapse density leave, what the
leftover looks like through the genes, the controls and the limits.

Run by run_ish_overview.py.
"""

import csv
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, false_discovery_control, spearmanr

from sepmap.config import DATA, SETTINGS
from sepmap.ish.gene_sets import set_test
from sepmap.ish.plotting import FIGURES, QUESTIONS, figure_file, figure_ref, ordinal
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
        raise FileNotFoundError(f"no numbers_*.csv in {tables}: run steps 13 to 29 first")
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
    ranking: pd.DataFrame,
    per_adult: pd.DataFrame,
    summary: pd.DataFrame,
    robustness: pd.DataFrame,
) -> pd.DataFrame:
    """Numbers of the gene by gene figures that their steps do not write.

    `ranking`, `per_adult`, `summary` and `robustness` are gene_ranking.csv,
    gene_ranking_per_adult.csv, robustness_summary.csv and ranking_robustness.csv.
    """
    nano = ranking[ranking["map"] == "nano"].set_index("symbol")
    auto = ranking[ranking["map"] == "auto"].set_index("symbol")
    q = ISH_ANALYSIS["q"]
    agreement = spearmanr(nano["rho"], auto["rho"].reindex(nano.index)).statistic
    passing = nano[nano["q_all"] < q].sort_values("rho", ascending=False)
    auto_top = auto[auto["q_all"] < q].sort_values("rho", ascending=False)
    pearson = robustness[robustness["variant"] == "pearson_log2"].set_index("symbol")
    pearson_rank = pearson["rho"].rank(ascending=False, method="min").get("Cacng8", 0)
    unreliable = passing[passing["reliability"] < 0.3]
    rows = [
        ("nano_auto_agreement", round(agreement, 3), "nano and auto gene orders"),
        ("nano_pass_genes", " ".join(passing.index), "genes past the nano null"),
        (
            "nano_pass_unreliable",
            " ".join(f"{g}:{r:.2f}" for g, r in unreliable["reliability"].items()),
            "of them, genes whose Allen experiments disagree (reliability below 0.3)",
        ),
        ("auto_top_genes", " ".join(auto_top.index[:5]), "highest past the auto null"),
        (
            "nano_band_p9",
            int((nano.loc[nano["p9_gene"], "p_spatial"] < 0.05).sum()),
            "P9's genes past their band, p < 0.05 before BH",
        ),
        (
            "nano_band_all",
            int((nano["p_spatial"] < 0.05).sum()),
            "genes past their band, p < 0.05 before BH",
        ),
        (
            "auto_band_all",
            int((auto["p_spatial"] < 0.05).sum()),
            "genes past the autofluorescence band, p < 0.05 before BH",
        ),
        (
            "nano_rank_all_Cacng8_pearson",
            int(pearson_rank),
            "Cacng8's rank among all genes under Pearson on log2",
        ),
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


def inside_or_past(p: float, q: float = 0.05) -> str:
    """'past its null' or 'inside its null', by the p and the level."""
    if p < q:
        return "past its null"
    return "inside its null"


def upper_first(text: str) -> str:
    """`text` with its first letter capital, the rest as it is."""
    return text[:1].upper() + text[1:]


def sets_past(n: dict[str, str]) -> tuple[int, int]:
    """The gene sets past their null on the map at the BH level, and those tested."""
    q = ISH_ANALYSIS["q"]
    values = [
        float(value)
        for key, value in n.items()
        if key.startswith("gene_sets.set_q_nano_") and value
    ]
    return sum(1 for v in values if v < q), len(values)


def input_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the overview and the inputs."""
    return {
        "overview": (
            "The question; the two parts of the argument, each with its key numbers "
            "and where it stands; the limit; and the map of the figures.",
            "Which figure answers which question, and where the argument stands.",
        ),
        "structures": (
            "A, the rule as a funnel; B, kept and left out by division; C, the adult "
            "map on the declared structures; D, how far each adult's zero and spread "
            "of zref move with the declared reference.",
            f"{int(num(n, 'structure_set.structures_declared'))} structures enter every "
            "comparison: grey matter measured in all ten adults, fixed by rule before "
            "any correlation. The declared reference moves each adult's zero by "
            f"{num(n, 'structure_set.zref_zero_shift_min'):.3f} to "
            f"{num(n, 'structure_set.zref_zero_shift_max'):.3f} and keeps the cohort "
            f"map's order (rho {num(n, 'structure_set.cohort_map_agreement'):.5f}).",
        ),
        "genes": (
            "A, the genes by panel and gene set; B, how many Allen experiments each "
            "gene has; C, how well two Allen experiments of one gene agree.",
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


def variants_text(n: dict[str, str]) -> str:
    """Part 1's check rows that change the model, in two sentences."""
    return (
        "With Gria1 to Gria4 in place of Gria1, "
        f"{pct(n, 'beyond.variant_left_four_subunits')} is left. On the "
        f"{n['beyond.variant_structures_psd95']} structures where PSD95 puncta are "
        "measured, the measured density leaves "
        f"{pct(n, 'beyond.variant_left_psd95')}, the mRNA panel "
        f"{pct(n, 'beyond.variant_left_panel')}, both together "
        f"{pct(n, 'beyond.variant_left_psd95_and_panel')}."
    )


def beyond_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the figures of part 1."""
    above = ", ".join(n["overview.leftover_above"].split())
    below = ", ".join(n["overview.leftover_below"].split())
    return {
        "beyond": (
            "A, the map against what Gria1, synapse density and autofluorescence "
            "predict for structures the fit has not seen, the structures furthest off "
            "named; B, what Gria1 predicts, what synapse density adds and what is "
            "left, with its interval; C, what the same model leaves of nano, of a map "
            "made only of Gria1 and synapse density (the floor) and of a map that is "
            "one Allen Gria1 experiment (the benchmark); D, one half of the cohort's "
            "leftover against the other's, over every split.",
            f"Gria1 mRNA alone predicts {pct(n, 'beyond.share_abundance')} of the map's "
            f"reproducible pattern, synapse density alone "
            f"{pct(n, 'beyond.share_density')}; the main model, Gria1, synapse density "
            f"and autofluorescence together, {pct(n, 'beyond.share_model')}. "
            f"{pct(n, 'beyond.left')} is left ({pct(n, 'beyond.left_lo')} to "
            f"{pct(n, 'beyond.left_hi')} over resampled structures), and one half of "
            "the cohort's leftover agrees with the other's at "
            f"{num(n, 'beyond.replication_leftover'):.2f}. On the "
            f"{n['beyond.calibration_structures']} structures of the calibration, nano "
            f"leaves {pct(n, 'beyond.nano_calibration_left')}, the floor "
            f"{pct(n, 'beyond.floor')} (nano minus it "
            f"{num(n, 'beyond.nano_minus_floor'):+.0%}, "
            f"{num(n, 'beyond.nano_minus_floor_lo'):+.0%} to "
            f"{num(n, 'beyond.nano_minus_floor_hi'):+.0%}), the benchmark "
            f"{pct(n, 'beyond.gria1_map_left')} (nano minus it "
            f"{num(n, 'beyond.nano_minus_gria1'):+.0%}, "
            f"{num(n, 'beyond.nano_minus_gria1_lo'):+.0%} to "
            f"{num(n, 'beyond.nano_minus_gria1_hi'):+.0%}). " + variants_text(n),
        ),
        "beyond_where": (
            "A, the map, the prediction and the leftover on three planes (red above "
            "prediction, blue below); B, the structures furthest from prediction, "
            "grey by how steady they are across the adults.",
            f"Most above prediction: {above}; most below: {below}.",
        ),
    }


def map_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the figures of the genes and the map."""
    gap_p = num(n, "gene_ranking.gap_p")
    return {
        "one_comparison": (
            "A and B, the nano map and Cacng8's Allen map as one rank per structure "
            "on the same plane; C and D, a gene's ranks against the map's, for Cacng8, "
            "which follows it, and Aqp4, an astrocyte gene, which does not.",
            f"Cacng8 {num(n, 'overview.rho_Cacng8'):+.2f} on "
            f"{n['overview.n_structures_Cacng8']} structures, Gria1 "
            f"{num(n, 'overview.rho_Gria1'):+.2f} on "
            f"{n['overview.n_structures_Gria1']}, Aqp4 "
            f"{num(n, 'overview.rho_Aqp4'):+.2f}. Much of a whole-brain rho is cortex "
            "and hippocampus against thalamus: the null and the test inside "
            "divisions come next.",
        ),
        "spatial_null": (
            "A, how widely unrelated smooth maps correlate, against shuffled ones; B, "
            "the surrogates' variogram on the map's; C, how often the ordinary and the "
            "spatial p fall below 0.05 on maps with no relation.",
            "Unrelated smooth maps correlate with an SD of rho "
            f"{num(n, 'spatial_null.sd_rho_pair'):.2f}, against "
            f"{num(n, 'spatial_null.sd_rho_shuffled_pair'):.2f} when one is "
            "shuffled. Random maps tested against nano give an ordinary p below 0.05 "
            f"in {pct(n, 'spatial_null.fpr_ordinary_map')} of tests, a spatial p in "
            f"{num(n, 'spatial_null.fpr_spatial_map'):.1%} (5% expected). The "
            "surrogates match the map's variogram within a median "
            f"{pct(n, 'spatial_null.variogram_misfit_nano')} over the matched range.",
        ),
        "top_genes": (
            "A, each gene's rho with the map against its pale null band, a red "
            "diamond where the gene also follows the leftover; B, its mean rho inside "
            "divisions; C, its rho with Gria1, the density terms and the main model's "
            "prediction; D, what it takes of the leftover when added to the main "
            "model, against maps of its smoothness.",
            top_genes_take(n),
        ),
        "cacng8_gria1": (
            "A and B, the map against Gria1 and against Cacng8, one dot per structure; "
            "C, Cacng8's lead over Gria1 against maps that follow both alike.",
            f"The map follows Cacng8 at {num(n, 'top_genes.rho_Cacng8'):+.2f} and "
            f"Gria1 at {num(n, 'top_genes.rho_Gria1'):+.2f} (inside divisions "
            f"{num(n, 'top_genes.rho_within_Cacng8'):+.2f} and "
            f"{num(n, 'top_genes.rho_within_Gria1'):+.2f}); the two genes agree at "
            f"{num(n, 'top_genes.rho_Gria1_Cacng8'):+.2f}. Cacng8 leads in "
            f"{n['overview.adults_cacng8_above_gria1']} of 10 adults, by "
            f"{num(n, 'gene_ranking.gap'):+.2f} "
            f"({num(n, 'gene_ranking.gap_boot_lo'):+.2f} to "
            f"{num(n, 'gene_ranking.gap_boot_hi'):+.2f} over resampled adults), a lead "
            f"{inside_or_past(gap_p)} of maps that follow both alike (p {gap_p:.2f}; "
            f"against unrelated maps p {num(n, 'gene_ranking.gap_p_unrelated'):.2f}).",
        ),
        "between_within": (
            "A, Cacng8 against the real map and against a map that knows only each "
            "structure's division; B, every gene's whole-brain rho against its rho "
            "with that division-only map; C, its whole-brain rho against its mean rho "
            "inside divisions (filled: past the within null).",
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


def top_genes_take(n: dict[str, str]) -> str:
    """The numbers of figure 07: the genes past the map's null and what they take."""
    follow = n["top_genes.top_follow_leftover_list"].split()
    alike = n["top_genes.top_take_alike_list"].split()
    plain = n["top_genes.top_take_list"].split()
    return (
        f"{n['top_genes.top_genes']} of {n['gene_ranking.nano_genes']} genes pass the "
        f"map's null after BH, {n['top_genes.top_within']} of them inside divisions "
        f"too; Cacng8 is first ({num(n, 'top_genes.rho_Cacng8'):+.2f}) and Gria1 "
        f"{ordinal(n['top_genes.rank_all_Gria1'])} "
        f"({num(n, 'top_genes.rho_Gria1'):+.2f}). {len(follow)} also follow the "
        f"leftover before correction ({', '.join(follow) or 'none'}). Added to the "
        f"main model, {len(alike)} take more of it than 95% of maps that relate to the "
        f"model as they do ({', '.join(alike) or 'none'}), {len(plain)} more than 95% "
        f"of their plain surrogates ({', '.join(plain) or 'none'}). Cacng8 takes "
        f"{num(n, 'top_genes.taken_Cacng8'):.1%} of the reproducible map (p "
        f"{num(n, 'top_genes.p_taken_alike_Cacng8'):.3f} against maps alike, "
        f"{num(n, 'top_genes.p_taken_Cacng8'):.3f} against plain surrogates)."
    )


def kinds_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from figures 10 and 11."""
    october = "gene_sets.localisation_5_Octobers_controls_positive_control_p"
    past, tested = sets_past(n)
    within_family = n["top_genes.family_leftover_pass_within_list"].split()
    glia_q = n.get("beyond.leftover_set_q_glia", "")
    glia = f"{float(glia_q):.2f}" if glia_q else "not tested"
    return {
        "gene_kinds": (
            "A, the gene sets fixed in advance, a dot per gene, each median against "
            "the band where it falls with 95% of the surrogates; B, the localisation "
            "genes' partial rho against their expression-matched postsynaptic "
            "controls; C, how often that test finds a planted difference of each size.",
            f"{past} of {tested} gene sets pass their null. Localisation genes "
            f"({n['gene_sets.set_size_localisation']}) have a median rho of "
            f"{num(n, 'gene_sets.set_median_nano_localisation'):+.2f} (spatial p "
            f"{num(n, 'gene_sets.set_p_nano_localisation'):.3f}), other postsynaptic "
            f"genes {num(n, 'gene_sets.set_median_nano_other postsynaptic'):+.2f}. "
            "Against matched controls, the subunit composite removed: "
            f"{num(n, 'gene_sets.localisation_matched_controls_difference'):+.3f} (p "
            f"{num(n, 'gene_sets.localisation_matched_controls_p'):.2f}); the test "
            "finds a planted difference of "
            f"{num(n, 'gene_sets.localisation_detectable_labels'):+.3f} in 80% of "
            "maps. Its positive control: GO pool p "
            f"{num(n, 'gene_sets.localisation_positive_control_p'):.2f}, pool of 5 "
            f"October p {num(n, october):.3f}.",
        ),
        "leftover": (
            "A, tier 1: Cacng8 against what Gria1 and synapse density leave; B and C, "
            "tier 2: the AMPA receptor complex family's median against the leftover's "
            "surrogates, and against expression-matched postsynaptic genes; D, tier "
            "3: every gene's p against the line BH asks it to cross.",
            f"Cacng8 {num(n, 'top_genes.leftover_rho_Cacng8'):+.3f}, spatial p "
            f"{num(n, 'top_genes.leftover_p_Cacng8'):.4f} "
            f"({ordinal(n['top_genes.leftover_rank_Cacng8'])} of "
            f"{n['beyond.leftover_genes']} by rho). The family "
            f"({n['top_genes.family_tested']} of {n['top_genes.family_listed']} members "
            "with a usable map): median "
            f"{num(n, 'top_genes.family_leftover_median'):+.3f}, p "
            f"{num(n, 'top_genes.family_leftover_p'):.3f} against the surrogates; "
            f"{num(n, 'top_genes.family_leftover_controls_difference'):+.3f} above "
            f"{n['top_genes.family_controls']} matched postsynaptic genes, p "
            f"{num(n, 'top_genes.family_leftover_controls_p'):.2f}; past BH within the "
            f"family: {', '.join(within_family) or 'none'}. Every gene: "
            f"{n['beyond.leftover_genes_pass']} past BH, the closest "
            f"{n['beyond.leftover_top_gene']} "
            f"({num(n, 'beyond.leftover_top_rho'):+.2f}, p "
            f"{num(n, 'beyond.leftover_top_p'):.2f}). The glia set, the lead of the "
            f"first, four-subunit run, has q {glia} here.",
        ),
    }


def control_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the controls, the limit and April."""
    return {
        "autofluorescence": (
            "A, adult by adult, nano's and autofluorescence's rho with Gria1 and "
            "Cacng8; B, the two genes against each map's own null; C, every gene's rho "
            "with autofluorescence against its rho with nano.",
            "Autofluorescence of the same sections correlates with Gria1 at "
            f"{num(n, 'gene_ranking.auto_rho_Gria1'):+.2f} (p "
            f"{num(n, 'gene_ranking.auto_p_Gria1'):.2f}) and with Cacng8 at "
            f"{num(n, 'gene_ranking.auto_rho_Cacng8'):+.2f} (p "
            f"{num(n, 'gene_ranking.auto_p_Cacng8'):.2f}); adult by adult, nano's rho "
            f"with Gria1 is {num(n, 'gene_ranking.per_adult_nano_Gria1_min'):.2f} to "
            f"{num(n, 'gene_ranking.per_adult_nano_Gria1_max'):.2f}, "
            "autofluorescence's "
            f"{num(n, 'gene_ranking.per_adult_auto_Gria1_min'):+.2f} to "
            f"{num(n, 'gene_ranking.per_adult_auto_Gria1_max'):+.2f}. The two gene "
            f"orders agree at {num(n, 'overview.nano_auto_agreement'):.2f}.",
        ),
        "robustness": (
            "A, how well each variant's gene order agrees with the primary's; B, "
            "where Cacng8 and Gria1 sit; C, the gap between them against the band of "
            "its null.",
            f"Over {n['overview.robustness_variants']} variants (statistic, borders, "
            "reading, inputs, structure set, and the route of 5 October) the gene "
            "order agrees with the primary at "
            f"{num(n, 'overview.robustness_agreement_min'):.2f} to 1.00; Cacng8 is "
            "first of P9's genes in every variant, Gria1 ranks "
            f"{n['overview.robustness_gria1_rank_min']} to "
            f"{n['overview.robustness_gria1_rank_max']}, and the gap stays between "
            f"{num(n, 'overview.robustness_gap_min'):+.2f} and "
            f"{num(n, 'overview.robustness_gap_max'):+.2f}.",
        ),
        "synaptome": (
            "A, the structures of the fit with a measured PSD95 density, by division; "
            "B, the PSD95 density against the marker mRNA composite it would replace; "
            "C, per density measure, what it predicts alone and what the main model "
            "leaves with it, on the structures PSD95 covers.",
            f"PSD95 puncta (Zhu et al. 2018) are measured in "
            f"{n['synaptome.fit_measured']} of the {n['synaptome.fit']} structures of "
            f"the fit ({pct(n, 'synaptome.fit_share')}); the rule fixed beforehand "
            f"asks for {n['synaptome.fit_needed']}. Where measured, PSD95 orders the "
            "structures as the marker composite does at "
            f"{num(n, 'synaptome.rho_fit_markers'):.2f}, and the two hemispheres of the "
            f"one mouse agree at {num(n, 'synaptome.hemispheres_rho'):.2f}. On those "
            f"{n['beyond.variant_structures_psd95']} structures PSD95 alone predicts "
            f"{pct(n, 'beyond.variant_density_alone_psd95')} of the map, the mRNA panel "
            f"{pct(n, 'beyond.variant_density_alone_panel')}; the main model leaves "
            f"{pct(n, 'beyond.variant_left_psd95')} with PSD95, "
            f"{pct(n, 'beyond.variant_left_panel')} with the panel, "
            f"{pct(n, 'beyond.variant_left_psd95_and_panel')} with both.",
        ),
        "green_channel": (
            "A, the three channels and how each pair agrees across structures; B, "
            "adult by adult; C, how much each channel varies across structures.",
            "In every adult the green (SEP) channel follows autofluorescence (rho "
            f"{num(n, 'green_channel.rho_sep_auto_min'):.2f} to "
            f"{num(n, 'green_channel.rho_sep_auto_max'):.2f}; nano follows it at only "
            f"{num(n, 'green_channel.rho_nano_auto_min'):.2f} to "
            f"{num(n, 'green_channel.rho_nano_auto_max'):.2f}) and varies across "
            f"structures about as little ({num(n, 'green_channel.range_sep_mean'):.2f}"
            f" against {num(n, 'green_channel.range_auto_mean'):.2f}, p90 - p10 of "
            f"log2; nano {num(n, 'green_channel.range_nano_mean'):.2f}).",
        ),
        "april_headline": (
            "A, each gene's rho in April against today's; B, the category ANOVA p "
            "under each choice it rests on; C, today's groups against the band of "
            "their median over the surrogates.",
            "April's gene order reproduces (rho "
            f"{num(n, 'overview.april_today_agreement'):.2f} over "
            f"{n['overview.april_today_genes']} genes). The group p was "
            f"{num(n, 'overview.anova_p_1'):.3f} with the hand split and "
            f"{num(n, 'overview.anova_p_2'):.2f} without; on today's rho "
            f"{num(n, 'overview.anova_p_4'):.3f} and {num(n, 'overview.anova_p_5'):.3f}"
            ". Today's F against the map's surrogates, which keep co-expressed genes "
            f"together, gives p {num(n, 'overview.anova_p_spatial'):.2f}.",
        ),
    }


def part1_verdict(n: dict[str, str]) -> str:
    """Where part 1 stands, in words that follow the numbers."""
    lo = num(n, "beyond.nano_minus_floor_lo")
    floor = num(n, "beyond.floor")
    gria1 = num(n, "beyond.gria1_map_left")
    over_gria1 = num(n, "beyond.nano_minus_gria1")
    over_gria1_lo = num(n, "beyond.nano_minus_gria1_lo")
    if lo >= 0.05:
        above = f"well above the floor ({floor:.0%}) on the same structures"
    elif lo > 0:
        above = (
            f"above the floor ({floor:.0%}) on the same structures, though the "
            "interval of the difference starts near zero"
        )
    else:
        above = (
            f"above the floor ({floor:.0%}) at its point value, with an interval "
            "that reaches it"
        )
    if over_gria1_lo > 0:
        benchmark = "larger than what a map of one Allen Gria1 experiment leaves"
    elif over_gria1 > 0:
        benchmark = (
            "larger than what a map of one Allen Gria1 experiment leaves at its point "
            "value, with an interval that reaches it"
        )
    else:
        benchmark = "no larger than what a map of one Allen Gria1 experiment leaves"
    return (
        f"Where it stands: the leftover is real and reproducible, {above}, and "
        f"{benchmark} ({gria1:.0%})."
    )


def family_text(n: dict[str, str]) -> str:
    """What the AMPA receptor complex family does with the leftover, in words that
    follow its two group tests."""
    q = ISH_ANALYSIS["q"]
    spatial = num(n, "top_genes.family_leftover_p")
    matched = num(n, "top_genes.family_leftover_controls_p")
    if spatial < q and matched >= q:
        return (
            "the AMPA receptor complex family follows it beyond the surrogates but "
            "no more than postsynaptic genes of the same expression"
        )
    if spatial < q:
        return (
            "the AMPA receptor complex family follows it beyond the surrogates and "
            "beyond matched postsynaptic genes"
        )
    if matched < q:
        return (
            "the AMPA receptor complex family follows it beyond matched postsynaptic "
            "genes, not beyond the surrogates"
        )
    return "the AMPA receptor complex family as a group does not"


def named_tests_text(n: dict[str, str]) -> str:
    """The tests named in advance for the leftover, in words that follow the numbers."""
    cacng8 = num(n, "top_genes.leftover_p_Cacng8")
    first = "follows it" if cacng8 < ISH_ANALYSIS["q"] else "does not follow it"
    return (
        f"Against what Gria1 and synapse density leave, Cacng8, named in advance, "
        f"{first} (p {cacng8:.4f}); {family_text(n)} (p "
        f"{num(n, 'top_genes.family_leftover_p'):.3f} and "
        f"{num(n, 'top_genes.family_leftover_controls_p'):.3f})."
    )


def part2_verdict(n: dict[str, str]) -> str:
    """Where part 2 stands, in words that follow the numbers."""
    q = ISH_ANALYSIS["q"]
    gap_p = num(n, "gene_ranking.gap_p")
    past, _ = sets_past(n)
    local_p = num(n, "gene_sets.localisation_matched_controls_p")
    if gap_p < q or past or local_p < q:
        return (
            "Where it stands: some of the genes that set surface receptor single the "
            f"map out beyond the null; read {figure_ref('cacng8_gria1')} and "
            f"{figure_ref('gene_kinds')} for which."
        )
    return (
        "Where it stands: consistent with the surface-fraction reading, not singling "
        "it out: the genes that follow the map most are maps much like Gria1 and "
        "synapse density, Cacng8's lead over Gria1 is inside its null, and the "
        "localisation genes do no better than matched controls."
    )


def gene_meanings(n: dict[str, str]) -> dict[str, str]:
    """What the figures of the genes mean, each verdict following its numbers."""
    q = ISH_ANALYSIS["q"]
    gap_p = num(n, "gene_ranking.gap_p")
    pre_p = num(n, "gene_sets.postsynaptic_against_presynaptic_p_spatial")
    glia_p = num(n, "gene_sets.postsynaptic_against_glia_p_spatial")
    local_p = num(n, "gene_sets.localisation_matched_controls_p")
    found = num(n, "gene_sets.localisation_detectable_labels")
    follow = n["top_genes.top_follow_leftover_list"].split()
    within = [
        g
        for g in ("Cacng8", "Dlg2", "Gria1")
        if num(n, f"divisions.p_within_spatial_{g}") < q
    ]
    out = {
        "one_comparison": "A whole-brain rho mixes fine agreement with the contrast "
        "between divisions; the next figures separate them.",
    }
    if 0.03 <= num(n, "spatial_null.fpr_spatial_map") <= 0.08:
        out["spatial_null"] = (
            "The spatial p keeps about 5% false positives where the ordinary one "
            "does not; every p of the gene analyses is a spatial p."
        )
    else:
        out["spatial_null"] = "The spatial p is not calibrated; no claim rests on it."
    out["top_genes"] = (
        "The genes that follow the map most are maps much like Gria1 and synapse "
        "density, which is why the main model predicts most of the map; they describe "
        "it, they are not separate evidence."
    )
    if follow:
        out["top_genes"] += (
            f" {len(follow)} of them also follow part of what the model leaves; only "
            "Cacng8's test was named in advance."
        )
    if gap_p < q:
        out["cacng8_gria1"] = "The map follows Cacng8 more closely than Gria1."
    else:
        out["cacng8_gria1"] = (
            "These maps cannot tell whether the nano map follows Cacng8 more closely "
            "than Gria1; Cacng8's own test is against the leftover "
            f"({figure_ref('leftover')})."
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
    if pre_p < q and glia_p < q:
        kind = "the map is postsynaptic-like by the criterion named in advance"
    elif pre_p < q or glia_p < q:
        kind = "one of the two contrasts named in advance passes, not both"
    else:
        kind = "neither contrast named in advance passes"
    if local_p >= q:
        out["gene_kinds"] = (
            f"The genes that set surface receptor do not stand out: {kind}, and no "
            f"localisation advantage larger than {found:+.2f} exists beyond the null."
        )
    else:
        out["gene_kinds"] = (
            f"The localisation genes stand above their matched controls; {kind}."
        )
    out["leftover"] = named_tests_text(n)
    if num(n, "top_genes.family_leftover_controls_p") >= q:
        out["leftover"] += (
            " What the family shares with the leftover is postsynaptic, not "
            "particular to the AMPA receptor complex."
        )
    if int(num(n, "beyond.leftover_genes_pass")) == 0:
        out["leftover"] += " No gene passes once every gene is corrected for."
    return out


def other_meanings(n: dict[str, str]) -> dict[str, str]:
    """What the inputs, part 1, the controls, the limit and April's headline mean."""
    q = ISH_ANALYSIS["q"]
    out = {
        "overview": "Every number on it is this run's; the story in words is "
        "docs/ISH_ANALYSIS.md.",
        "structures": "Every comparison runs on structures every adult measures; the "
        "declared reference shifts each brain's zref without reordering the map.",
        "genes": "A gene measured once is only as good as one Allen mouse; a low rho "
        "of an unreliable gene says little.",
        "beyond": upper_first(part1_verdict(n).removeprefix("Where it stands: "))
        + " It is what the model does not predict, not a measurement of the surface "
        "fraction: the words are 'not predicted by Gria1 expression or synapse "
        "density', not 'beyond gene expression'.",
        "beyond_where": "The departure from prediction sits in particular structures, "
        "steadily across the adults; a claim about one structure needs its own null.",
        "robustness": "Cacng8 first of P9's genes holds under every choice; Gria1's "
        "rank moves with them, so it is quoted with its null, not as a place.",
    }
    auto_p = num(n, "gene_ranking.auto_p_Gria1")
    above = int(num(n, "overview.adults_nano_above_auto_Gria1"))
    if auto_p >= q and above == 10:
        out["autofluorescence"] = (
            "Gria1 and Cacng8 follow the label, not the tissue, in every adult; the "
            "tissue has a gene pattern of its own, a different one."
        )
    else:
        out["autofluorescence"] = "The tissue shares part of the label's ranking."
    psd95 = num(n, "beyond.variant_density_alone_psd95")
    panel = num(n, "beyond.variant_density_alone_panel")
    if n["synaptome.in_main_model"] == "True":
        out["synaptome"] = "The measured density replaces the mRNA panel in the model."
    elif psd95 < panel:
        out["synaptome"] = (
            "The measured density covers too little of the fit to replace the mRNA "
            "panel, and where it exists it predicts the map less well; by keeping the "
            "panel, the main model is the harder test of part 1."
        )
    else:
        out["synaptome"] = (
            "The measured density covers too little of the fit to replace the mRNA "
            "panel, though where it exists it predicts the map at least as well."
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


# the main figures of each part of the walk, in order
PARTS = (
    ("The question", ("overview",)),
    ("The inputs", ("structures", "genes")),
    (
        "Part 1: the map is not explained by Gria1 expression or synapse density",
        ("beyond", "beyond_where"),
    ),
    (
        "Part 2: what else it is: the genes that follow the map and its leftover",
        (
            "one_comparison",
            "spatial_null",
            "top_genes",
            "cacng8_gria1",
            "between_within",
            "gene_kinds",
            "leftover",
        ),
    ),
    ("Controls", ("autofluorescence", "robustness", "synaptome")),
    ("The limit", ("green_channel",)),
    ("April's headline", ("april_headline",)),
)

# the detailed versions of each main figure, and what each adds
SUPPLEMENTS = {
    "genes": (
        (
            "genes_detail",
            "the section QC of P9's experiments, reliability against expression, and "
            "what was left out and repaired",
        ),
    ),
    "beyond": (
        (
            "beyond_budget",
            "the map against each predictor, the budget beside its four-subunit check "
            "row and control F, every draw of the calibration, and the leftover under "
            "every check row, fold and structure set",
        ),
    ),
    "one_comparison": (
        (
            "one_comparison_detail",
            "the maps as measured and as ranks for nano, Cacng8, Gria1 and Aqp4, and "
            "the steps of one comparison",
        ),
    ),
    "spatial_null": (
        (
            "spatial_null_detail",
            "the nano map and three surrogates on a plane, and Cacng8's and Gria1's rho "
            "against their nulls",
        ),
    ),
    "top_genes": (
        (
            "gene_ranking",
            "P9's 100 genes one by one with their null bands and autofluorescence's "
            "rho, and the Cacng8 - Gria1 gap with the adults' interval and each pairing "
            "of Allen experiments",
        ),
    ),
    "between_within": (
        (
            "between_within_detail",
            "five genes division by division, the choice of null for the within rho, "
            "and the genes highest inside divisions",
        ),
    ),
    "gene_kinds": (
        (
            "gene_sets",
            "every set with its genes named, the two contrasts named in advance, and "
            "where each set comes from",
        ),
        (
            "localisation",
            "the label null, the positive control with both control pools, the "
            "matching on expression, and every test of the design",
        ),
    ),
    "leftover": (
        (
            "leftover_genes",
            "the genes closest to the leftover with their null bands, and every gene "
            "set against it",
        ),
        (
            "ampa_family",
            "each member of the family against the leftover, and the family on the "
            "map itself",
        ),
    ),
    "autofluorescence": (
        (
            "autofluorescence_detail",
            "the genes' rho with each map, how many pass each null at three "
            "thresholds, and the genes that pass",
        ),
    ),
    "robustness": (("robustness_detail", "every gene under four of the choices"),),
    "synaptome": (
        (
            "synaptome_detail",
            "the two hemispheres of the one mouse, and each density's agreement with "
            "the mRNA terms, Gria1 and the maps",
        ),
    ),
    "green_channel": (
        (
            "green_channel_detail",
            "one adult's raw channels on a plane, and each channel against Gria1, with "
            "what is left of SEP once autofluorescence is out",
        ),
    ),
    "april_headline": (
        ("april_headline_detail", "April's ten violins beside today's, gene by gene"),
    ),
}

# the run script that draws each guided figure
DRAWN_BY = {
    "overview": "run_ish_overview.py",
    "structures": "run_structure_set.py",
    "genes": "run_ish_gene_table.py",
    "genes_detail": "run_ish_gene_table.py",
    "beyond": "run_beyond_figures.py",
    "beyond_budget": "run_beyond_figures.py",
    "beyond_where": "run_beyond_figures.py",
    "one_comparison": "run_ish_gene_ranking.py",
    "one_comparison_detail": "run_ish_gene_ranking.py",
    "spatial_null": "run_ish_gene_ranking.py",
    "spatial_null_detail": "run_ish_gene_ranking.py",
    "top_genes": "run_ish_top_genes.py",
    "gene_ranking": "run_ish_gene_ranking.py",
    "cacng8_gria1": "run_ish_top_genes.py",
    "between_within": "run_ish_divisions.py",
    "between_within_detail": "run_ish_divisions.py",
    "gene_kinds": "run_ish_gene_sets.py",
    "gene_sets": "run_ish_gene_sets.py",
    "localisation": "run_ish_gene_sets.py",
    "leftover": "run_ish_top_genes.py",
    "leftover_genes": "run_beyond_figures.py",
    "ampa_family": "run_ish_top_genes.py",
    "autofluorescence": "run_ish_gene_ranking.py",
    "autofluorescence_detail": "run_ish_gene_ranking.py",
    "robustness": "run_ish_robustness.py",
    "robustness_detail": "run_ish_robustness.py",
    "synaptome": "run_beyond_figures.py",
    "synaptome_detail": "run_synaptome.py",
    "green_channel": "run_sep_channel_check.py",
    "green_channel_detail": "run_sep_channel_check.py",
    "april_headline": "run_ish_overview.py",
    "april_headline_detail": "run_ish_overview.py",
}


def index_head() -> str:
    """The head of figures/README.md: what it is, the argument, how the figures go."""
    return f"""# The ISH analysis, figure by figure

Written by `run_ish_overview.py` from the tables of this run, so every number below
is this run's. The story, with what each result means and does not mean, is
`docs/ISH_ANALYSIS.md` in the code repository.

The figures follow one argument in two parts. Part 1: the adult nano map across
structures is not explained by Gria1 expression and synapse density; a
reproducible part is left over, above what Allen-to-Allen mismatch alone leaves
(figures {FIGURES["beyond"]} and {FIGURES["beyond_where"]}). Part 2: it is
therefore something else, and the reading the data support is the surface
fraction of the receptor (trafficking, scaffolding); the genes that follow the
map and its leftover are the corroboration, Cacng8 (TARP gamma-8) first (figures
{FIGURES["one_comparison"]} to {FIGURES["leftover"]}). The surface fraction is
an interpretation, not a measurement: a total-GluA1 stain on the same brains
would measure it, and the green channel cannot ({figure_ref("green_channel")}).

Each main figure answers one question in two to four panels; its detailed version,
the same number with an s, holds every panel and number of the analysis.
"""


def index_tail() -> str:
    """The tail of figures/README.md: the sheets beside the guided figures."""
    return f"""## Sheets

- `qc/00_flagged.png`: every experiment with a section set missing, kept as a true
  absence, or kept at a step in expression, one row each, sections as columns (red
  set missing, hatched kept as true absence, a dot kept at a step). The sheet to review
  `mapping/ish_section_exceptions.csv` from.
- `qc/<gene>_<experiment>.png`: one sheet per Allen experiment: its median energy
  section by section along its own axis, with the local reference and the flags, and
  the orientation check at its middle section (`run_ish_section_qc.py --sheets`).
- `genes/<gene>.png`: Cacng8, Gria1, Grm5, Dlg2 and Aqp4: the nano and gene rank
  maps, the scatter of ranks with one fitted line per division, the whole-brain and
  within-division rho with their p (`run_ish_divisions.py --sheets`).
- `top_genes/<gene>.png`: every gene of figures {FIGURES["top_genes"]} and
  {FIGURES["ampa_family"]} (the genes past
  the map's null, Gria1, Cacng8 and the AMPA receptor complex family): its map,
  against the map and against the leftover, what it takes of the leftover against
  maps of its smoothness, what it is like, its numbers and GO terms
  (`run_ish_top_genes.py --sheets`; PNG only).
"""


def figure_index(n: dict[str, str]) -> str:
    """figures/README.md: the guided walk, each main figure with its question, what
    to look at and what to take from it, and its detailed versions."""
    walk = input_walk(n) | beyond_walk(n) | map_walk(n) | kinds_walk(n)
    walk |= control_walk(n)
    meaning = gene_meanings(n) | other_meanings(n)
    lines = [index_head()]
    for part, keys in PARTS:
        lines.append(f"## {part}\n")
        for key in keys:
            look, take = walk[key]
            name = figure_file(key)
            lines += [
                f"### {FIGURES[key]}. {QUESTIONS[key]}\n",
                f"![{name}]({name})\n",
                f"**Look at.** {look}\n",
                f"**Takeaway.** {take} {meaning[key]}\n",
            ]
            details = SUPPLEMENTS.get(key, ())
            for detail, adds in details:
                other = figure_file(detail)
                lines.append(f"- [{FIGURES[detail]}]({other}) in detail: {adds}.")
            if details:
                lines.append("")
            lines.append(
                f"Drawn by `{DRAWN_BY[key]}`; every figure also as an EPS of the same "
                "name.\n"
            )
    lines.append(index_tail())
    return "\n".join(lines)


# ===== The overview figure =====


def part1_content(n: dict[str, str]) -> dict:
    """The box of part 1: its claim, three key numbers, where it stands."""
    return dict(
        part="part 1",
        claim="Part 1. The map is not explained by Gria1 expression or synapse "
        "density: a reproducible part is left over",
        numbers=[
            (
                pct(n, "beyond.share_model"),
                "of the map's reproducible pattern is predicted by Gria1 and synapse "
                "density, on structures the fit has not seen (Gria1 alone "
                f"{pct(n, 'beyond.share_abundance')}) (figure {FIGURES['beyond']})",
            ),
            (
                f"{pct(n, 'beyond.left')} left",
                f"{pct(n, 'beyond.left_lo')} to {pct(n, 'beyond.left_hi')} over "
                "structures; the two halves of the cohort leave the same leftover "
                f"(rho {num(n, 'beyond.replication_leftover'):.2f})",
            ),
            (
                f"{pct(n, 'beyond.nano_calibration_left')} vs {pct(n, 'beyond.floor')}",
                "what the model leaves of nano, and of a map made only of Gria1 and "
                f"density (the floor): {num(n, 'beyond.nano_minus_floor'):+.0%} "
                f"({num(n, 'beyond.nano_minus_floor_lo'):+.0%} to "
                f"{num(n, 'beyond.nano_minus_floor_hi'):+.0%}); one Allen Gria1 "
                f"experiment leaves {pct(n, 'beyond.gria1_map_left')}",
            ),
        ],
        stands=part1_verdict(n),
    )


def part2_content(n: dict[str, str]) -> dict:
    """The box of part 2: its claim, three key numbers, where it stands."""
    within = ""
    if num(n, "top_genes.p_within_Cacng8") < ISH_ANALYSIS["q"]:
        within = ", inside divisions too"
    return dict(
        part="part 2",
        claim="Part 2. It is therefore something else: the reading the data support "
        "is the surface fraction of the receptor",
        numbers=[
            (
                f"{num(n, 'top_genes.rho_Cacng8'):+.2f}",
                "Cacng8 (TARP gamma-8, an AMPA receptor auxiliary subunit) follows the "
                f"map, {ordinal(n['top_genes.rank_all_Cacng8'])} of "
                f"{n['gene_ranking.nano_genes']} genes{within}; Gria1 "
                f"{num(n, 'top_genes.rho_Gria1'):+.2f} (figures "
                f"{FIGURES['top_genes']}, {FIGURES['cacng8_gria1']})",
            ),
            (
                f"p {num(n, 'top_genes.leftover_p_Cacng8'):.4f}",
                "Cacng8, named in advance, follows what Gria1 and synapse density "
                f"leave, and takes {num(n, 'top_genes.taken_Cacng8'):.1%} of the map "
                f"from it (figure {FIGURES['leftover']})",
            ),
            (
                f"p {num(n, 'top_genes.family_leftover_p'):.3f}",
                f"{upper_first(family_text(n))} (p "
                f"{num(n, 'top_genes.family_leftover_controls_p'):.2f})",
            ),
        ],
        stands=part2_verdict(n),
    )


def overview_content(n: dict[str, str]) -> dict:
    """The text of the overview figure: the question, the two parts, the limit."""
    return dict(
        question="The adult nano map orders the brain's structures the same way in "
        "both halves of the cohort (rho "
        f"{num(n, 'beyond.half_agreement'):.2f}). Is that order how much GluA1 mRNA a "
        "structure makes, or how many synapses it has? Allen's ISH maps measure both.",
        parts=[part1_content(n), part2_content(n)],
        limit="The surface fraction is an interpretation, not a measurement. The "
        "green channel follows autofluorescence in every adult (rho "
        f"{num(n, 'green_channel.rho_sep_auto_min'):.2f} to "
        f"{num(n, 'green_channel.rho_sep_auto_max'):.2f}), so it cannot give total "
        "receptor; a total-GluA1 stain on some of the same brains would measure it "
        f"(figure {FIGURES['green_channel']}).",
    )


def figure_map() -> list[tuple[str, list[tuple[str, str, str]]]]:
    """The groups of the overview's figure map: each main figure's number, question
    and detailed versions, by part of the walk, the overview itself left out."""
    groups = []
    for part, keys in PARTS[1:]:
        figures = []
        for key in keys:
            details = ", ".join(FIGURES[d] for d, _ in SUPPLEMENTS.get(key, ()))
            figures.append((FIGURES[key], QUESTIONS[key], details))
        groups.append((part.split(":")[0], figures))
    return groups
