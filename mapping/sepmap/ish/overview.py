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


def input_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the figures of the inputs."""
    return {
        "overview": (
            "The argument in the boxes at the top, each part with this run's numbers; "
            "then one row per step of the walk: the question, the figures that answer "
            "it, the numbers of this run, and what stays open.",
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
        "beyond_budget": (
            "A to C, the nano map against Gria1 mRNA, against synapse density and "
            "against the whole main model; D, the variance budget of the main model and "
            "of its four-subunit check row: red is what is left, hatched the "
            "calibration floor, the dashed mark where the leftover of a map that is one "
            "Allen Gria1 experiment would begin; E, the leftover beside maps whose "
            "answer is known (the two Gria1 clusters are the two halves of its "
            "experiments); F, one half of the cohort's leftover against the other's; "
            "G, the leftover under other folds, models and structures.",
            f"Gria1 mRNA alone predicts {pct(n, 'beyond.share_abundance')} of the map's "
            f"reproducible pattern, synapse density alone "
            f"{pct(n, 'beyond.share_density')}; the main model, Gria1, synapse density "
            f"and autofluorescence together, {pct(n, 'beyond.share_model')}. "
            f"{pct(n, 'beyond.left')} is left ({pct(n, 'beyond.left_lo')} to "
            f"{pct(n, 'beyond.left_hi')} over resampled structures), and one half of "
            "the cohort's leftover agrees with the other's at "
            f"{num(n, 'beyond.replication_leftover'):.2f}. " + variants_text(n) + " On "
            f"the {n['beyond.calibration_structures']} structures of the calibration, "
            f"nano leaves {pct(n, 'beyond.nano_calibration_left')}, a map made only of "
            f"Gria1 and synapse density {pct(n, 'beyond.floor')} (the floor; nano minus "
            f"it {num(n, 'beyond.nano_minus_floor'):+.0%}, "
            f"{num(n, 'beyond.nano_minus_floor_lo'):+.0%} to "
            f"{num(n, 'beyond.nano_minus_floor_hi'):+.0%}), and a map that is one Allen "
            f"Gria1 experiment {pct(n, 'beyond.gria1_map_left')} "
            f"({pct(n, 'beyond.gria1_map_left_A')} and "
            f"{pct(n, 'beyond.gria1_map_left_B')} from the two halves; nano minus it "
            f"{num(n, 'beyond.nano_minus_gria1'):+.0%}, "
            f"{num(n, 'beyond.nano_minus_gria1_lo'):+.0%} to "
            f"{num(n, 'beyond.nano_minus_gria1_hi'):+.0%}). Control F: the "
            f"components of the {n['beyond.control_f_genes']} genes measured in every "
            f"structure predict {pct(n, 'beyond.control_f_share')}, and on its own "
            f"calibration nano leaves {pct(n, 'beyond.control_f_nano_lo')} to "
            f"{pct(n, 'beyond.control_f_nano_hi')} against "
            f"{pct(n, 'beyond.control_f_floor_lo')} to "
            f"{pct(n, 'beyond.control_f_floor_hi')} for a map made of those genes.",
        ),
        "beyond_where": (
            "A, the map, the prediction and the leftover on three planes (red above "
            "prediction, blue below); B, the structures furthest from prediction, "
            "grey by how steady they are across the adults.",
            f"Most above prediction: {above}; most below: {below}.",
        ),
    }


def gene_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the figures of part 2."""
    gap_p = num(n, "gene_ranking.gap_p")
    october = "gene_sets.localisation_5_Octobers_controls_"
    four = "gene_sets.localisation_matched_controls_four_subunits_removed_"
    glia_q = n.get("beyond.leftover_set_q_glia", "")
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
            "A, each gene's bar against its pale null band: an open circle marks a "
            "gene past its band (p < 0.05 before correction), a bold name one past "
            "it after BH; the yellow dot is its rho with autofluorescence. B, the "
            "Cacng8 - Gria1 gap against maps related to both genes alike (filled) "
            "and against unrelated maps (outline), with the adults' interval and the "
            "pairings of experiments beneath. C, the two genes with the label and "
            "with the tissue.",
            f"{n['overview.nano_band_p9']} of P9's {n['gene_ranking.nano_p9_genes']} "
            f"genes are past their band, {n['gene_ranking.nano_pass_p9']} after BH; "
            f"of all {n['gene_ranking.nano_genes']}, {n['overview.nano_band_all']} and "
            f"{n['gene_ranking.nano_pass_all']}. Cacng8 is first "
            f"({num(n, 'gene_ranking.nano_rho_Cacng8'):+.2f}) and Gria1 "
            f"{ordinal(n['overview.rank_p9_Gria1'])} of P9's genes "
            f"({num(n, 'gene_ranking.nano_rho_Gria1'):+.2f}), both past the null. "
            f"Cacng8 leads Gria1 by {num(n, 'gene_ranking.gap'):+.2f} "
            f"({num(n, 'gene_ranking.gap_boot_lo'):+.2f} to "
            f"{num(n, 'gene_ranking.gap_boot_hi'):+.2f} over resampled adults); "
            f"against maps related to both alike the lead is {inside_or_past(gap_p)} "
            f"(p {gap_p:.2f}; band {num(n, 'gene_ranking.gap_equal_lo'):+.2f} to "
            f"{num(n, 'gene_ranking.gap_equal_hi'):+.2f}), against unrelated maps p "
            f"{num(n, 'gene_ranking.gap_p_unrelated'):.2f}.",
        ),
        "gene_sets": (
            "A, one column per set fixed in advance: its median (black bar) against "
            "the band where the median of the same genes falls with 95% of the "
            "surrogates (sets too small to test have none); B, the two contrasts "
            "named in advance, and under them the genes of the presynaptic set.",
            "Localisation genes (transport, anchoring, auxiliary subunits; "
            f"{n['gene_sets.set_size_localisation']}) have a median rho of "
            f"{num(n, 'gene_sets.set_median_nano_localisation'):+.2f} (spatial p "
            f"{num(n, 'gene_sets.set_p_nano_localisation'):.3f}), other postsynaptic "
            f"genes {num(n, 'gene_sets.set_median_nano_other postsynaptic'):+.2f} (p "
            f"{num(n, 'gene_sets.set_p_nano_other postsynaptic'):.3f}), presynaptic "
            f"{num(n, 'gene_sets.set_median_nano_presynaptic'):+.2f}, glia "
            f"{num(n, 'gene_sets.set_median_nano_glia'):+.2f}. Postsynaptic genes "
            "sit above the presynaptic set by "
            f"{num(n, 'gene_sets.postsynaptic_against_presynaptic_difference'):+.2f} "
            "(spatial p "
            f"{num(n, 'gene_sets.postsynaptic_against_presynaptic_p_spatial'):.3f}) "
            "and above glia by "
            f"{num(n, 'gene_sets.postsynaptic_against_glia_difference'):+.2f} (p "
            f"{num(n, 'gene_sets.postsynaptic_against_glia_p_spatial'):.3f}).",
        ),
        "localisation": (
            "C first: the positive control, with both control pools; then A and B, "
            "the test; E, how large a difference the test finds when one is planted; "
            "F, every test of the design.",
            "Positive control: GO pool "
            f"{num(n, 'gene_sets.localisation_positive_control_difference'):+.3f} (p "
            f"{num(n, 'gene_sets.localisation_positive_control_p'):.3f}), pool of 5 "
            f"October {num(n, october + 'positive_control_difference'):+.3f} (p "
            f"{num(n, october + 'positive_control_p'):.3f}). Localisation against "
            "matched controls, the subunit composite removed: "
            f"{num(n, 'gene_sets.localisation_matched_controls_difference'):+.3f} (p "
            f"{num(n, 'gene_sets.localisation_matched_controls_p'):.3f}; spatial p "
            f"{num(n, 'gene_sets.localisation_p_spatial'):.3f}); with the four subunits "
            f"removed as separate terms {num(n, four + 'difference'):+.3f} (p "
            f"{num(n, four + 'p'):.3f}). The test finds a planted difference of "
            f"{num(n, 'gene_sets.localisation_detectable_labels'):+.3f} in 80% of "
            "maps.",
        ),
        "between_within": (
            "B, whole-brain rho against rho with a map that knows only each "
            "structure's division; C, whole-brain rho against the mean rho inside "
            "divisions (filled: past the within null); D, five genes division by "
            "division; E, why the surrogates and not a shuffle inside divisions; F, "
            "the genes highest inside divisions with their q.",
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
        "leftover_genes": (
            "A, the genes closest to the leftover against its null, each surrogate "
            "put through the same fit; B, the gene sets against the leftover.",
            f"{n['beyond.leftover_genes_pass']} of {n['beyond.leftover_genes']} genes "
            "follow the leftover past BH ("
            f"{n['beyond.leftover_genes_p05']} below p 0.05 before correction); the "
            f"closest is {n['beyond.leftover_top_gene']} "
            f"({num(n, 'beyond.leftover_top_rho'):+.2f}, p "
            f"{num(n, 'beyond.leftover_top_p'):.4f}, q "
            f"{num(n, 'beyond.leftover_top_q'):.2f}); Cacng8 "
            f"{num(n, 'beyond.leftover_rho_Cacng8'):+.2f} (p "
            f"{num(n, 'beyond.leftover_p_Cacng8'):.4f}, "
            f"{ordinal(n['beyond.leftover_rank_Cacng8'])} of all). Glia "
            f"{num(n, 'beyond.leftover_set_glia'):+.2f} (q "
            f"{float(glia_q) if glia_q else float('nan'):.3f}), localisation "
            f"{num(n, 'beyond.leftover_set_localisation'):+.2f} (q "
            f"{num(n, 'beyond.leftover_set_q_localisation'):.2f}). Cacng8 was named "
            "for this leftover in advance, its uncorrected p the test, and the AMPA "
            f"receptor complex family with it ({figure_ref('ampa_family')}); the sets "
            "and the other genes describe it.",
        ),
    }


def top_gene_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from figures 12 to 14."""
    follow = n["top_genes.top_follow_leftover_list"].split()
    alike = n["top_genes.top_take_alike_list"].split()
    plain = n["top_genes.top_take_list"].split()
    untested = n["top_genes.family_not_tested"].split()
    past = n["top_genes.family_leftover_pass_within_list"].split()
    gap_p = num(n, "gene_ranking.gap_p")
    return {
        "top_genes": (
            "A, each gene's rho with the map against its pale null band, a red "
            "diamond where the gene also follows the leftover; B, its mean rho inside "
            "divisions; C, its rho with Gria1, the density terms and the main model's "
            "prediction; D, what it takes of the leftover when added to the main "
            "model, against maps of its smoothness.",
            f"{n['top_genes.top_genes']} genes pass the map's null after BH, "
            f"{n['top_genes.top_within']} of them inside divisions too. "
            f"{len(follow)} also follow the leftover before correction "
            f"({', '.join(follow) or 'none'}). Added to the main model, "
            f"{len(alike)} take more of it than 95% of maps that relate to the model "
            f"as they do ({', '.join(alike) or 'none'}), {len(plain)} more than 95% of "
            f"their plain surrogates ({', '.join(plain) or 'none'}). Cacng8 takes "
            f"{num(n, 'top_genes.taken_Cacng8'):.1%} of the reproducible map (p "
            f"{num(n, 'top_genes.p_taken_alike_Cacng8'):.3f} against maps alike, "
            f"{num(n, 'top_genes.p_taken_Cacng8'):.3f} against plain surrogates); its "
            "rho with the main model's prediction is "
            f"{num(n, 'top_genes.rho_prediction_Cacng8'):.2f}.",
        ),
        "cacng8_gria1": (
            "A and B, the map against Gria1 and against Cacng8, one dot per "
            "structure; C, what Gria1 and synapse density leave against Cacng8; D, "
            "Cacng8's lead over Gria1 against maps that follow both alike.",
            f"The map follows Cacng8 at {num(n, 'top_genes.rho_Cacng8'):+.2f} and "
            f"Gria1 at {num(n, 'top_genes.rho_Gria1'):+.2f} (inside divisions "
            f"{num(n, 'top_genes.rho_within_Cacng8'):+.2f} and "
            f"{num(n, 'top_genes.rho_within_Gria1'):+.2f}); the two genes agree at "
            f"{num(n, 'top_genes.rho_Gria1_Cacng8'):+.2f}. Cacng8's lead, "
            f"{num(n, 'gene_ranking.gap'):+.2f}, is {inside_or_past(gap_p)} of maps "
            f"that follow both alike (p {gap_p:.2f}). Against the leftover Cacng8 "
            f"gives {num(n, 'top_genes.leftover_rho_Cacng8'):+.2f}, p "
            f"{num(n, 'top_genes.leftover_p_Cacng8'):.4f}, the test named in advance.",
        ),
        "ampa_family": (
            "A, each member of the family against the leftover, with its null band; "
            "B, the family's median against the leftover's surrogates; C, against "
            "postsynaptic genes of the same expression; D and E, the same on the map "
            "itself, to describe the family.",
            f"{n['top_genes.family_tested']} of the family's "
            f"{n['top_genes.family_listed']} members have a usable map (not "
            f"{', '.join(untested)}). Against the leftover the family's median is "
            f"{num(n, 'top_genes.family_leftover_median'):+.3f}, spatial p "
            f"{num(n, 'top_genes.family_leftover_p'):.3f}; against "
            f"{n['top_genes.family_controls']} matched postsynaptic genes "
            f"{num(n, 'top_genes.family_leftover_controls_difference'):+.3f}, p "
            f"{num(n, 'top_genes.family_leftover_controls_p'):.3f}; "
            f"{len(past)} member past BH within the family "
            f"({', '.join(past) or 'none'}). On the map the median is "
            f"{num(n, 'top_genes.family_nano_median'):+.3f} (p "
            f"{num(n, 'top_genes.family_nano_p'):.3f}), against the same controls "
            f"{num(n, 'top_genes.family_nano_controls_difference'):+.3f} (p "
            f"{num(n, 'top_genes.family_nano_controls_p'):.3f}).",
        ),
    }


def control_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the controls, the limits and April."""
    auto_top = ", ".join(n["overview.auto_top_genes"].split())
    return {
        "autofluorescence": (
            "A, each gene's rho with the autofluorescence map against its rho with "
            "nano; C, adult by adult, nano and autofluorescence against Gria1 and "
            "Cacng8; D, how many genes each map passes at three thresholds.",
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
            f"a gene pattern of its own: past BH, {n['gene_ranking.auto_pass_all']} "
            f"genes for autofluorescence against {n['gene_ranking.nano_pass_all']} for "
            f"nano (highest {auto_top}); past their band before correction, "
            f"{n['overview.auto_band_all']} against {n['overview.nano_band_all']}; "
            "the two gene orders agree at "
            f"{num(n, 'overview.nano_auto_agreement'):.2f}.",
        ),
        "robustness": (
            "A, how well each variant's gene order agrees with the primary's; B, "
            "where Cacng8 and Gria1 sit; C, the gap between them against the band of "
            "its null; D, gene by gene under four variants.",
            f"Over {n['overview.robustness_variants']} variants (statistic, borders, "
            "reading, inputs, structure set, and the route of 5 October) the gene "
            "order agrees with the primary at "
            f"{num(n, 'overview.robustness_agreement_min'):.2f} to 1.00; Cacng8 is "
            "first of P9's genes in every variant (of all genes it is "
            f"{ordinal(n['overview.nano_rank_all_Cacng8_pearson'])} under Pearson), "
            f"Gria1 ranks {n['overview.robustness_gria1_rank_min']} to "
            f"{n['overview.robustness_gria1_rank_max']}, and the gap stays between "
            f"{num(n, 'overview.robustness_gap_min'):+.2f} and "
            f"{num(n, 'overview.robustness_gap_max'):+.2f}.",
        ),
        "green_channel": (
            "A and C, which channel follows which, adult by adult; B, one adult's raw "
            "channels; D, how much each channel varies across structures; E, each "
            "channel against Gria1, and what is left of SEP once its "
            "autofluorescence part is taken out.",
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
            "A, April's groups (grey) beside today's (orange), a line per gene; B, "
            "gene by gene; C, the ANOVA p under each choice it rests on; D, today's "
            "groups against the band of their median over the surrogates.",
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


def upper_first(text: str) -> str:
    """`text` with its first letter capital, the rest as it is."""
    return text[:1].upper() + text[1:]


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


def named_tests_text(n: dict[str, str]) -> str:
    """The tests named in advance for the leftover, in words that follow the numbers."""
    q = ISH_ANALYSIS["q"]
    cacng8 = num(n, "top_genes.leftover_p_Cacng8")
    spatial = num(n, "top_genes.family_leftover_p")
    matched = num(n, "top_genes.family_leftover_controls_p")
    first = "follows it" if cacng8 < q else "does not follow it"
    if spatial < q and matched >= q:
        family = (
            "the AMPA receptor complex family follows it beyond the surrogates but "
            "no more than postsynaptic genes of the same expression"
        )
    elif spatial < q:
        family = (
            "the AMPA receptor complex family follows it beyond the surrogates and "
            "beyond matched postsynaptic genes"
        )
    elif matched < q:
        family = (
            "the AMPA receptor complex family follows it beyond matched postsynaptic "
            "genes, not beyond the surrogates"
        )
    else:
        family = "the AMPA receptor complex family as a group does not"
    return (
        f"Against what Gria1 and synapse density leave, Cacng8, named in advance, "
        f"{first} (p {cacng8:.4f}); {family} (p {spatial:.3f} and {matched:.3f})."
    )


def part2_verdict(n: dict[str, str]) -> str:
    """Where part 2 stands, in words that follow the numbers."""
    q = ISH_ANALYSIS["q"]
    gap_p = num(n, "gene_ranking.gap_p")
    sets_past = sum(
        1
        for key, value in n.items()
        if key.startswith("gene_sets.set_q_nano_") and value and float(value) < q
    )
    local_p = num(n, "gene_sets.localisation_matched_controls_p")
    if gap_p < q or sets_past or local_p < q:
        return (
            "Where it stands: some of the genes that set surface receptor single the "
            "map out beyond the null; read the figures for which."
        )
    return (
        "Where it stands: consistent with the surface-fraction reading, and not "
        "singling it out: the "
        "map follows Cacng8 and Gria1 beyond the null, inside divisions too, but "
        "Cacng8's lead over Gria1 is inside the null of maps related to both alike, "
        "no gene set passes, and the localisation genes do no better than matched "
        "controls. " + named_tests_text(n)
    )


def top_gene_meanings(n: dict[str, str]) -> dict[str, str]:
    """What figures 12 to 14 mean, each verdict following its numbers."""
    q = ISH_ANALYSIS["q"]
    follow = n["top_genes.top_follow_leftover_list"].split()
    gap_p = num(n, "gene_ranking.gap_p")
    cacng8 = num(n, "top_genes.leftover_p_Cacng8")
    out = {}
    if follow:
        out["top_genes"] = (
            "The genes that follow the map most are maps much like Gria1 and synapse "
            "density, which is why the main model predicts most of the map; "
            f"{len(follow)} of them also follow part of what it leaves "
            f"({', '.join(follow)}). They describe the map; Cacng8's test is named "
            "in advance, the others' are not."
        )
    else:
        out["top_genes"] = (
            "The genes that follow the map most are maps much like Gria1 and synapse "
            "density; none of them follows what the main model leaves."
        )
    if gap_p >= q and cacng8 < q:
        out["cacng8_gria1"] = (
            "These maps cannot tell whether the nano map follows Cacng8 more "
            "closely than Gria1, but Cacng8, named in advance, follows what Gria1 "
            "and synapse density leave: the result that points the way the "
            "surface-fraction reading does, a correspondence and not a measurement."
        )
    elif cacng8 < q:
        out["cacng8_gria1"] = (
            "The map follows Cacng8 more closely than Gria1, and Cacng8 follows what "
            "Gria1 and synapse density leave."
        )
    else:
        out["cacng8_gria1"] = (
            "Cacng8 does not follow what Gria1 and synapse density leave beyond its null."
        )
    out["ampa_family"] = named_tests_text(n) + (
        " The family was defined by the proteomics of native complexes and by GO, "
        "before the main model ran."
    )
    return out


def meanings(n: dict[str, str]) -> dict[str, str]:
    """What each figure means, in one or two sentences; a verdict follows its numbers."""
    q = ISH_ANALYSIS["q"]
    calibrated = 0.03 <= num(n, "spatial_null.fpr_spatial_map") <= 0.08
    gap_p = num(n, "gene_ranking.gap_p")
    pre_p = num(n, "gene_sets.postsynaptic_against_presynaptic_p_spatial")
    glia_p = num(n, "gene_sets.postsynaptic_against_glia_p_spatial")
    go_control = num(n, "gene_sets.localisation_positive_control_p")
    october = "gene_sets.localisation_5_Octobers_controls_"
    october_control = num(n, october + "positive_control_p")
    local_p = num(n, "gene_sets.localisation_matched_controls_p")
    found = num(n, "gene_sets.localisation_detectable_labels")
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
        "beyond_budget": f"{num(n, 'beyond.left'):.0%} of the reproducible map is "
        "not predicted by Gria1 expression or synapse density. "
        + upper_first(part1_verdict(n).removeprefix("Where it stands: "))
        + " What the leftover is, these data do not say: it is what the model does "
        "not predict; the components of many panel genes (control F) predict most of "
        "the map, though no single gene follows the leftover, so the words are 'not "
        "predicted by Gria1 expression or synapse density', not 'beyond gene "
        "expression'.",
        "beyond_where": "The departure from prediction sits in particular structures, "
        "steadily across the adults; a claim about one structure needs its own null.",
        "one_comparison": "A whole-brain rho mixes fine agreement with the contrast "
        "between divisions; the next figures separate them.",
    }
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
        lead = (
            "but against maps related to both genes alike these maps cannot tell "
            "which of the two it follows more closely"
        )
    out["gene_ranking"] = (
        f"The map follows a TARP's pattern and its own subunit's mRNA beyond "
        f"chance, {lead}: consistent with the surface-fraction reading, not "
        "evidence for it."
    )
    if pre_p < q and glia_p < q:
        kind = "the map is postsynaptic-like by the criterion named in advance"
    elif pre_p < q or glia_p < q:
        kind = (
            "one of the two contrasts named in advance passes, so the criterion "
            "(both) is not met; the presynaptic set here is mostly cell-type markers"
        )
    else:
        kind = "neither contrast named in advance passes"
    passing = (
        f"{sets_past} gene sets pass their null"
        if sets_past
        else ("No gene set passes its null")
    )
    out["gene_sets"] = (
        f"{passing}; {kind}. Localisation genes sit no higher than other "
        f"postsynaptic genes ({figure_ref('localisation')} asks it at equal "
        "expression)."
    )
    if go_control >= q and october_control < q and local_p >= q:
        out["localisation"] = (
            "With the GO control pool the positive control fails; with the pool of "
            "5 October it works. The localisation genes do no better than their "
            f"matched controls, and the test would have found a difference of "
            f"{found:+.2f} in 80% of maps: no localisation advantage larger than that."
        )
    elif local_p >= q:
        out["localisation"] = (
            "The localisation genes do no better than their matched controls; no "
            f"advantage larger than {found:+.2f} (found in 80% of maps)."
        )
    else:
        out["localisation"] = "The localisation genes stand above their controls."
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
    glia_q = n.get("beyond.leftover_set_q_glia", "")
    glia = bool(glia_q) and float(glia_q) < q
    if int(float(n["beyond.leftover_genes_pass"])) == 0:
        out["leftover_genes"] = (
            "No single gene follows the leftover past BH over every gene"
            + (
                "; the glia set does, a test not named in advance, so an unplanned "
                "lead, not a finding"
                if glia
                else ", and no gene set does"
            )
            + ". The genes do not say what the leftover is."
        )
    else:
        out["leftover_genes"] = "Some genes follow the leftover past BH."
    cacng8_p = num(n, "beyond.leftover_p_Cacng8")
    if cacng8_p < q:
        out["leftover_genes"] += (
            " Cacng8, the one gene named in advance, follows it "
            f"(p {cacng8_p:.4f}, uncorrected as named)."
        )
    out |= top_gene_meanings(n)
    auto_p = num(n, "gene_ranking.auto_p_Gria1")
    above = int(float(n["overview.adults_nano_above_auto_Gria1"]))
    if auto_p >= q and above == 10:
        out["autofluorescence"] = (
            "Gria1 and Cacng8 follow the label, not the tissue, in every adult; the "
            "tissue has a gene pattern of its own, and the two gene orders differ "
            f"(rho {num(n, 'overview.nano_auto_agreement'):.2f})."
        )
    else:
        out["autofluorescence"] = "The tissue shares part of the label's ranking."
    out["robustness"] = (
        "Cacng8 first of P9's genes holds under every choice; Gria1's rank moves "
        "with them, so it is quoted with its null, not as a place."
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
        "Part 1: the map is not explained by Gria1 expression or synapse density",
        ("beyond_budget", "beyond_where"),
    ),
    (
        "Part 2: do the genes that set surface receptor follow the map better than "
        "abundance genes?",
        (
            "one_comparison",
            "spatial_null",
            "gene_ranking",
            "gene_sets",
            "localisation",
            "between_within",
            "leftover_genes",
            "top_genes",
            "cacng8_gria1",
            "ampa_family",
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
    "leftover_genes": "run_beyond_figures.py",
    "top_genes": "run_ish_top_genes.py",
    "cacng8_gria1": "run_ish_top_genes.py",
    "ampa_family": "run_ish_top_genes.py",
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
structures is not satisfactorily explained by Gria1 expression or by synapse density
(synaptic markers and postsynaptic-density genes, or PSD95 puncta where measured); a
sizeable, reproducible part is left over. Part 2: the reading the data support is
the surface fraction of the receptor, and the genes that regulate surface AMPA
receptors (Cacng8, a TARP; trafficking and scaffolding genes) are tested as
corroboration of it, against abundance genes, unrelated genes and the spatial null.
Three tests carry that corroboration: the Cacng8 - Gria1 gap ({gap}, panel B),
localisation genes against matched controls once the subunit composite is removed
({local}), and the genes against the leftover itself ({leftover}). The genes that
follow the map most are then described one by one ({top}), and the tests named in
advance for the leftover, Cacng8 and the AMPA receptor complex family, are
{family}. The surface fraction is an interpretation, not a measurement: a
total-GluA1 stain on the same brains would measure it, and the green channel
cannot ({limit}).
"""

INDEX_TAIL = """## Sheets

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
- `top_genes/<gene>.png`: every gene of figures 12 and 14 (the genes past the
  map's null, Gria1, Cacng8 and the AMPA receptor complex family): its map, against
  the map and against the leftover, what it takes of the leftover against maps of
  its smoothness, what it is like, its numbers and GO terms
  (`run_ish_top_genes.py --sheets`; PNG only).
"""


def figure_index(n: dict[str, str]) -> str:
    """figures/README.md: the guided walk, each figure with its question and caption."""
    walk = input_walk(n) | beyond_walk(n) | gene_walk(n) | top_gene_walk(n)
    walk |= control_walk(n)
    meaning = meanings(n)
    lines = [
        INDEX_HEAD.format(
            gap=figure_ref("gene_ranking"),
            local=figure_ref("localisation"),
            leftover=figure_ref("leftover_genes"),
            top=figure_ref("top_genes"),
            family=f"{figure_ref('cacng8_gria1')} and {figure_ref('ampa_family')}",
            limit=figure_ref("green_channel"),
        )
    ]
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
    """The boxes at the top of the overview: the question, the two parts, the limit.

    Each part gives its claim, this run's numbers, and where it stands in words
    worked out from those numbers (part1_verdict, part2_verdict).
    """
    gap_p = num(n, "gene_ranking.gap_p")
    sets = sum(
        1 for key, value in n.items() if key.startswith("gene_sets.set_q_nano_") and value
    )
    sets_past = sum(
        1
        for key, value in n.items()
        if key.startswith("gene_sets.set_q_nano_")
        and value
        and float(value) < ISH_ANALYSIS["q"]
    )
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
            "The map is not satisfactorily explained by Gria1 expression or synapse "
            f"density. Gria1 mRNA alone predicts {pct(n, 'beyond.share_abundance')} "
            "of the reproducible map; Gria1, synapse density and autofluorescence "
            f"together {pct(n, 'beyond.share_model')}. "
            f"{pct(n, 'beyond.left')} is left ({pct(n, 'beyond.left_lo')} to "
            f"{pct(n, 'beyond.left_hi')}) and it replicates across halves of the "
            f"cohort at {num(n, 'beyond.replication_leftover'):.2f}. On the same "
            f"structures nano leaves {pct(n, 'beyond.nano_calibration_left')}, a map "
            f"made only of Gria1 and synapse density {pct(n, 'beyond.floor')} (the "
            f"floor; difference {num(n, 'beyond.nano_minus_floor'):+.0%}, "
            f"{num(n, 'beyond.nano_minus_floor_lo'):+.0%} to "
            f"{num(n, 'beyond.nano_minus_floor_hi'):+.0%}), a map that is one Allen "
            f"Gria1 experiment {pct(n, 'beyond.gria1_map_left')} "
            f"(figures {FIGURES['beyond_budget']:02d} and "
            f"{FIGURES['beyond_where']:02d}).\n" + part1_verdict(n),
        ),
        (
            "Part 2",
            "The map is therefore something else; the reading the data support is "
            "the surface fraction of the receptor, shaped by trafficking and "
            "scaffolding. If so, the genes that regulate surface AMPA receptors "
            "(Cacng8, trafficking and scaffolding genes) should track the map better "
            "than abundance genes or unrelated genes, beyond the spatial null "
            f"(figures {FIGURES['one_comparison']:02d} to "
            f"{FIGURES['ampa_family']:02d}). This run: Cacng8 leads Gria1 by "
            f"{num(n, 'gene_ranking.gap'):+.2f} (p {gap_p:.2f} against maps related "
            f"to both alike); {sets_past} of {sets} gene sets pass; localisation "
            "against matched controls "
            f"{num(n, 'gene_sets.localisation_matched_controls_difference'):+.3f} "
            "(a difference of "
            f"{num(n, 'gene_sets.localisation_detectable_labels'):+.2f} would be "
            f"found); {n['beyond.leftover_genes_pass']} of "
            f"{n['beyond.leftover_genes']} genes follow the leftover past BH, and "
            "Cacng8, named for it in advance, has spatial p "
            f"{num(n, 'beyond.leftover_p_Cacng8'):.4f}; the AMPA receptor complex "
            f"family, p {num(n, 'top_genes.family_leftover_p'):.3f} against the "
            "surrogates and "
            f"{num(n, 'top_genes.family_leftover_controls_p'):.3f} against matched "
            "controls.\n" + part2_verdict(n),
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
            open="The one true absence of the exceptions list (Glra1, section 61) is "
            "proposed, to review on qc/00_flagged.png. The adults' ages are not "
            "recorded (Allen: P56).",
        ),
        dict(
            part="part 1",
            title="Part 1: not explained by Gria1 expression or synapse density",
            figures=figs("beyond_budget", "beyond_where"),
            question="How much of the map do Gria1 expression and synapse density "
            "predict, and is what they leave real?",
            lines=[
                f"Gria1 alone {pct(n, 'beyond.share_abundance')}, synapse density "
                f"alone {pct(n, 'beyond.share_density')} of the reproducible map",
                "the main model, Gria1, density and autofluorescence: "
                f"{pct(n, 'beyond.share_model')}, so {pct(n, 'beyond.left')} left "
                f"({pct(n, 'beyond.left_lo')} to {pct(n, 'beyond.left_hi')})",
                "check rows: Gria1 to Gria4 "
                f"{pct(n, 'beyond.variant_left_four_subunits')} left; PSD95 puncta "
                f"{pct(n, 'beyond.variant_left_psd95')} against the mRNA panel "
                f"{pct(n, 'beyond.variant_left_panel')} on the "
                f"{n['beyond.variant_structures_psd95']} structures PSD95 measures",
                "the leftover replicates across halves of the cohort at "
                f"{num(n, 'beyond.replication_leftover'):.2f} (the map's reliability "
                f"alone implies {num(n, 'beyond.replication_implied'):.2f})",
                f"same structures: nano {pct(n, 'beyond.nano_calibration_left')}, the "
                f"floor {pct(n, 'beyond.floor')} (difference "
                f"{num(n, 'beyond.nano_minus_floor'):+.0%}, 95% "
                f"{num(n, 'beyond.nano_minus_floor_lo'):+.0%} to "
                f"{num(n, 'beyond.nano_minus_floor_hi'):+.0%}), one Allen Gria1 "
                f"experiment as the map {pct(n, 'beyond.gria1_map_left')}",
                "control F, the components of "
                f"{n['beyond.control_f_genes']} genes: "
                f"{1 - num(n, 'beyond.control_f_share'):.0%} left, above its own floor; "
                f"{n['beyond.controls_passed']} of 7 controls pass",
            ],
            open="Which benchmark a claim uses: the floor or the one-experiment Gria1 "
            "map. What the leftover is: a total-GluA1 stain separates the surface "
            "fraction from translation, turnover, subunit composition and nanobody "
            "access. A claim about one structure needs its own null.",
        ),
    ]


def gene_rows(n: dict[str, str]) -> list[dict]:
    """The overview's rows of part 2."""
    gap_p = num(n, "gene_ranking.gap_p")
    october = "gene_sets.localisation_5_Octobers_controls_positive_control_p"
    unreliable = n.get("overview.nano_pass_unreliable", "")
    caveat = ""
    if unreliable:
        caveat = (
            " ("
            + ", ".join(
                f"{g.split(':')[0]}'s Allen experiments disagree, {g.split(':')[1]}"
                for g in unreliable.split()
            )
            + ")"
        )
    glia_q = n.get("beyond.leftover_set_q_glia", "")
    return [
        dict(
            part="part 2",
            title="Part 2: gene by gene",
            figures=figs("one_comparison", "spatial_null", "gene_ranking", "top_genes"),
            question="Which genes' maps order the structures as the nano map does, "
            "beyond a map with its smoothness? Does Cacng8 lead Gria1?",
            lines=[
                f"{n['gene_ranking.nano_pass_all']} of {n['gene_ranking.nano_genes']} "
                "genes past the spatial null after BH; "
                f"{n['gene_ranking.nano_pass_p9']} of P9's "
                f"{n['gene_ranking.nano_p9_genes']}{caveat}",
                f"Cacng8 first, {num(n, 'gene_ranking.nano_rho_Cacng8'):+.2f}; Gria1 "
                f"{num(n, 'gene_ranking.nano_rho_Gria1'):+.2f}, "
                f"{ordinal(n['overview.rank_p9_Gria1'])} of P9's genes",
                f"Cacng8 - Gria1 gap {num(n, 'gene_ranking.gap'):+.2f} "
                f"({num(n, 'gene_ranking.gap_boot_lo'):+.2f} to "
                f"{num(n, 'gene_ranking.gap_boot_hi'):+.2f} over adults): "
                f"{inside_or_past(gap_p)} of maps related to both alike, p "
                f"{gap_p:.2f}",
                "the spatial p is calibrated: "
                f"{num(n, 'spatial_null.fpr_spatial_map'):.1%} false positives at "
                f"0.05 (ordinary p: {pct(n, 'spatial_null.fpr_ordinary_map')})",
                f"of the {n['top_genes.top_genes']} past the null, "
                f"{n['top_genes.top_follow_leftover']} also follow the leftover; "
                f"Cacng8 takes {num(n, 'top_genes.taken_Cacng8'):.1%} of the map from "
                f"it (p {num(n, 'top_genes.p_taken_alike_Cacng8'):.3f} against maps "
                "alike to the model)",
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
                ", p "
                f"{num(n, 'gene_sets.postsynaptic_against_presynaptic_p_spatial'):.3f};"
                " above glia, p "
                f"{num(n, 'gene_sets.postsynaptic_against_glia_p_spatial'):.3f}",
                "localisation against matched controls, the subunit composite "
                "removed: "
                f"{num(n, 'gene_sets.localisation_matched_controls_difference'):+.3f}"
                f", p {num(n, 'gene_sets.localisation_matched_controls_p'):.2f}; "
                "found from "
                f"{num(n, 'gene_sets.localisation_detectable_labels'):+.3f} in 80%",
                "positive control: GO pool p "
                f"{num(n, 'gene_sets.localisation_positive_control_p'):.2f}; "
                f"5 October's pool p {num(n, october):.3f}",
            ],
            open="A positive control that comes out with the GO control pool.",
        ),
        dict(
            part="part 2",
            title="Part 2: inside divisions",
            figures=figs("between_within"),
            question="Does a gene follow the map inside divisions, or only through "
            "the contrast between them?",
            lines=[
                "a division-only map orders the genes as the map does: "
                f"{num(n, 'divisions.agreement_division_only_all'):.2f}; inside "
                f"divisions {num(n, 'divisions.agreement_within_all'):.2f}",
                f"{n['divisions.pass_all_spatial']} genes past the within null: "
                f"Cacng8 {num(n, 'divisions.rho_within_Cacng8'):+.2f}, Dlg2 "
                f"{num(n, 'divisions.rho_within_Dlg2'):+.2f}, Gria1 "
                f"{num(n, 'divisions.rho_within_Gria1'):+.2f}",
            ],
            open="",
        ),
        dict(
            part="part 2",
            title="Part 2: the leftover itself",
            figures=figs("leftover_genes", "cacng8_gria1", "ampa_family"),
            question="Does any gene, or any kind of gene, follow what Gria1 and "
            "synapse density leave?",
            lines=[
                f"{n['beyond.leftover_genes_pass']} of {n['beyond.leftover_genes']} "
                f"genes past BH ({n['beyond.leftover_genes_p05']} below p 0.05); "
                f"closest {n['beyond.leftover_top_gene']} "
                f"{num(n, 'beyond.leftover_top_rho'):+.2f}; Cacng8 "
                f"{num(n, 'beyond.leftover_rho_Cacng8'):+.2f} (p "
                f"{num(n, 'beyond.leftover_p_Cacng8'):.4f})",
                f"glia {num(n, 'beyond.leftover_set_glia'):+.2f} (q "
                f"{float(glia_q) if glia_q else float('nan'):.3f}), localisation "
                f"{num(n, 'beyond.leftover_set_localisation'):+.2f}; of these, only "
                "Cacng8 named in advance",
                "the AMPA receptor complex family, named in advance: median "
                f"{num(n, 'top_genes.family_leftover_median'):+.3f}, p "
                f"{num(n, 'top_genes.family_leftover_p'):.3f} against the "
                "surrogates, "
                f"{num(n, 'top_genes.family_leftover_controls_p'):.3f} against "
                "matched controls",
            ],
            open="The glia set's lead from the first, four-subunit model of 8 "
            "October: unplanned, not pursued.",
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
                "null after BH; the two gene orders agree at "
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
                "follows autofluorescence at only "
                f"{num(n, 'green_channel.rho_nano_auto_min'):.2f} to "
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
                f"split, {num(n, 'overview.anova_p_2'):.2f} without; on today's rho "
                f"{num(n, 'overview.anova_p_4'):.3f}",
                "today's F against the map's surrogates (co-expressed genes kept "
                f"together): p {num(n, 'overview.anova_p_spatial'):.2f}",
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
# on the branch that builds them (post-ish, not merged yet)
ITEMS = (
    ("A1", "one declared structure set and zref reference: built for the ISH line"),
    ("A2", "section QC, failed sections set missing: built; one absence proposed"),
    ("A3", "one adult profile; Spearman, Pearson, borders, structure sets: built"),
    ("A4", "which structures sit reliably above the median: open"),
    ("A5", "nano against autofluorescence, structure by structure: open"),
    ("A6", "agreement within divisions, per-gene sheets: built"),
    ("A7", "a spatial null for every rho, the gap and the sets: built"),
    ("A8", "every gene against the autofluorescence map: built"),
    ("A9", "the gene table and its documentation: built"),
    ("A10", "per-gene videos: optional, not built"),
)
