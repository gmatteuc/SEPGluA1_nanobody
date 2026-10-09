"""The ISH line on one page: the overview figure and the index of the figures.

The index of the guided figures (figures/README.md) and the text of the overview
figure: for each figure the question, what to look at and what to take from it,
with the numbers of this run (ish.numbers gathers them). A verdict follows its
numbers, so a rerun cannot leave the index saying what its numbers contradict. The
figures' order is the argument's: the question, the inputs, how much of the map
Gria1 and synapse density leave, what the leftover looks like through the genes,
the controls and the limits. A few numbers the index quotes are not written by
their own steps (the structures furthest from prediction, gene lists, counts over
the adults); they are computed here, for numbers_overview.csv.

Run by run_ish_overview.py.
"""

import pandas as pd
from scipy.stats import spearmanr

from sepmap.adult import beyond_regression
from sepmap.adult.profiles import ADULTS
from sepmap.config import SETTINGS
from sepmap.ish import divisions, gene_ranking, gene_sets
from sepmap.ish.figure_index import (
    DRAWN_BY,
    FIGURE_NUMBERS,
    PARTS,
    QUESTIONS,
    SUPPLEMENTS,
    figure_file,
    figure_ref,
)
from sepmap.ish.numbers import numbers_frame, ordinal
from sepmap.ish.spatial_null import ALPHA
from sepmap.structures import FIGURES

# the BH level the figures and the text count genes and sets at; the reliability
# below which a gene's Allen map does not reproduce
ISH_ANALYSIS = SETTINGS["ish_analysis"]
ISH_PANEL_TEST = SETTINGS["ish_panel_test"]

INDEX = FIGURES / "README.md"

# the false-positive rates of the spatial p on random maps (5% expected) within which
# the index reads it as calibrated
CALIBRATED = (0.03, 0.08)


# ===== Numbers the steps do not write =====


def leftover_extremes(n_each: int = 4) -> pd.DataFrame:
    """The structures furthest above and below prediction (analysis 4), for the text."""
    table = pd.read_csv(beyond_regression.REGRESSION)
    table = table.sort_values("residual", ascending=False)
    rows = [
        ("leftover_above", " ".join(table["acronym"].head(n_each)), "most above"),
        ("leftover_below", " ".join(table["acronym"].tail(n_each)[::-1]), "most below"),
        ("leftover_top_ranks", round(table["residual"].iloc[0], 1), "largest, ranks"),
    ]
    return numbers_frame(rows)


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
    nano = gene_ranking.map_rows(ranking, "nano")
    auto = gene_ranking.map_rows(ranking, "auto")
    q = ISH_ANALYSIS["q"]
    agreement = spearmanr(nano["rho"], auto["rho"].reindex(nano.index)).statistic
    passing = nano[nano["q_all"] < q].sort_values("rho", ascending=False)
    auto_top = auto[auto["q_all"] < q].sort_values("rho", ascending=False)
    pearson = robustness[robustness["variant"] == "pearson_log2"].set_index("symbol")
    pearson_rank = pearson["rho"].rank(ascending=False, method="min").get("Cacng8", 0)
    low = ISH_PANEL_TEST["min_reliability"]
    unreliable = passing[passing["reliability"] < low]
    rows = [
        ("nano_auto_agreement", round(agreement, 3), "nano and auto gene orders"),
        ("nano_pass_genes", " ".join(passing.index), "genes past the nano null"),
        (
            "nano_pass_unreliable",
            " ".join(f"{g}:{r:.2f}" for g, r in unreliable["reliability"].items()),
            f"of them, genes whose Allen experiments disagree (reliability below {low})",
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
    return numbers_frame(rows)


def gene_extras() -> pd.DataFrame:
    """Lists of genes the text names that their steps write only as tables.

    The genes past the within-division null (within_division.csv), the genes past
    the nano null counted by gene set (gene_ranking.csv), and the localisation
    genes with the highest partial rho against the GO control pool
    (localisation_test.csv).
    """
    q = ISH_ANALYSIS["q"]
    within, _ = divisions.load_within()
    past = within[within["q_all_spatial"] < q].sort_values("rho_within", ascending=False)
    nano = gene_ranking.map_rows(gene_ranking.load_ranking(), "nano")
    nano = nano[nano["q_all"] < q]
    first_set = nano["gene_sets"].str.split(";").str[0].replace("", "none")
    counts = first_set.value_counts()
    local = pd.read_csv(gene_sets.LOCALISATION)
    local = local[(local["side"] == "localisation") & local["pool"].str.contains("GO")]
    top = local.sort_values("rho_partial", ascending=False)["symbol"].head(5)
    rows = [
        ("within_pass_genes", " ".join(past["symbol"]), "past the within null"),
        ("localisation_top_partial", " ".join(top), "highest partial rho"),
    ]
    for name, count in counts.items():
        key = name.replace(" ", "_")
        rows.append((f"nano_pass_set_{key}", int(count), f"genes past the null: {name}"))
    return numbers_frame(rows)


# ===== The index of the figures =====


def num(n: dict[str, str], key: str) -> float:
    """One number of the run, by 'step.name'."""
    return float(n[key])


def pct(n: dict[str, str], key: str) -> str:
    """A share of the run as a whole percentage: '27%'."""
    return f"{num(n, key):.0%}"


def points_with_interval(n: dict[str, str], key: str) -> str:
    """A difference of shares in points with its 95% interval.

    '+7 points, 95% -12 to +27'.
    """
    lo, hi = num(n, f"{key}_lo"), num(n, f"{key}_hi")
    return f"{100 * num(n, key):+.0f} points, 95% {100 * lo:+.0f} to {100 * hi:+.0f}"


def gap_words(n: dict[str, str]) -> str:
    """The Cacng8 - Gria1 gap against maps that follow both alike, in words.

    The words follow its p (the test fixed in advance, a lead as large either way)
    and its band.
    """
    p = num(n, "gene_ranking.gap_p")
    first = num(n, "gene_ranking.gap_equal_first_as_large")
    second = num(n, "gene_ranking.gap_equal_second_as_large")
    shares = (
        f"{first:.1%} of those maps give a Cacng8 lead as large, {second:.1%} a Gria1 "
        f"lead; p {p:.3f} counts both"
    )
    if p < ALPHA:
        return f"a lead past what maps that follow both alike give ({shares})"
    if num(n, "gene_ranking.gap") > num(n, "gene_ranking.gap_equal_hi"):
        return (
            "a lead just past the 95% band of maps that follow both alike, inside the "
            f"test fixed in advance ({shares})"
        )
    return f"a lead inside what maps that follow both alike give ({shares})"


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
            f"Part 1: Gria1 and synapse density predict {pct(n, 'beyond.share_model')} "
            f"of the map's reproducible pattern and leave {pct(n, 'beyond.left')} "
            f"({pct(n, 'beyond.left_lo')} to {pct(n, 'beyond.left_hi')}). Part 2: "
            f"Cacng8 follows the map at {num(n, 'top_genes.rho_Cacng8'):+.2f} and what "
            f"the model leaves at p {num(n, 'top_genes.leftover_p_Cacng8'):.4f}.",
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


def check_margin(n: dict[str, str], key: str) -> str:
    """A check row's nano minus its floor in points, with its interval."""
    return points_with_interval(n, f"beyond.check_nano_minus_floor_{key}")


def checks_text(n: dict[str, str]) -> str:
    """Part 1's check rows that change the model, in two sentences."""
    return (
        "Curved, the bound of the straight model, it leaves "
        f"{pct(n, 'beyond.check_left_curved')} (nano minus its floor "
        f"{check_margin(n, 'curved')}); with autofluorescence "
        f"{pct(n, 'beyond.check_left_autofluorescence')}, with the first proposal "
        f"for density {pct(n, 'beyond.check_left_first_proposal')}, with the 11 marker "
        f"genes {pct(n, 'beyond.check_left_marker_panel')}, with psd_pc1 "
        f"{pct(n, 'beyond.check_left_psd_pc1')}, with Gria1 to Gria4 "
        f"{pct(n, 'beyond.check_left_four_subunits')}. On the "
        f"{n['beyond.check_structures_psd95']} structures where PSD95 puncta are "
        f"measured, PSD95 as the density term leaves {pct(n, 'beyond.check_left_psd95')}"
        f" (nano minus its floor {check_margin(n, 'psd95')}), the density genes "
        f"{pct(n, 'beyond.check_left_psd95_main')} ({check_margin(n, 'psd95_main')}); "
        f"on the {n['beyond.check_structures_large']} structures of 0.4 mm3 or more "
        f"{pct(n, 'beyond.check_left_large')} ({check_margin(n, 'large')}); with nano on "
        f"the Allen grid {pct(n, 'beyond.check_left_allen_grid')} "
        f"({check_margin(n, 'allen_grid')})."
    )


def lost_text(n: dict[str, str]) -> str:
    """The declared structures the rule of the fit leaves out, in words."""
    parts = []
    for key, value in n.items():
        if not key.startswith("beyond.lost_") or key.count("_") != 1:
            continue
        gene = key.removeprefix("beyond.lost_")
        divisions = [
            f"{k.removeprefix(f'beyond.lost_{gene}_')} {v}"
            for k, v in n.items()
            if k.startswith(f"beyond.lost_{gene}_") and not k.endswith("_structures")
        ]
        parts.append(f"{value} without {gene}: {', '.join(divisions)}")
    return "; ".join(parts)


def beyond_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from the figures of part 1."""
    above = ", ".join(n["overview.leftover_above"].split())
    below = ", ".join(n["overview.leftover_below"].split())
    genes = ", ".join(n["beyond.main_density_genes"].split())
    return {
        "beyond": (
            "A, the map against what Gria1 and synapse density predict for structures "
            "the fit has not seen, the structures furthest off named; B, the "
            "reproducible map in four parts, what only Gria1 predicts, what the two "
            "share, what only synapse density predicts and what is left, with its "
            "interval; C, what the same model leaves of nano and of a map made only of "
            "Gria1 and synapse density (the floor); D, one half of the cohort's "
            "leftover against the other's, over every split. Synapse density is the "
            f"mean rank of {genes}, postsynaptic genes chosen by their agreement with "
            f"PSD95 punctum density without the map ({figure_ref('density_markers')}). "
            f"The seven controls are in {figure_ref('beyond_controls')}, the check rows "
            f"in {figure_ref('beyond_budget')}.",
            f"On {n['beyond.structures']} structures (left out, {lost_text(n)}), Gria1 "
            f"mRNA alone predicts {pct(n, 'beyond.share_abundance')} of the map's "
            "reproducible pattern, synapse density alone "
            f"{pct(n, 'beyond.share_density')}, both {pct(n, 'beyond.share_model')}: "
            f"Gria1 only {num(n, 'beyond.part_gria1_only'):+.0%}, shared "
            f"{num(n, 'beyond.part_shared'):+.0%}, density only "
            f"{num(n, 'beyond.part_density_only'):+.0%}. {pct(n, 'beyond.left')} is "
            f"left ({pct(n, 'beyond.left_lo')} to {pct(n, 'beyond.left_hi')} over "
            "resampled structures), and one half of the cohort's leftover agrees with "
            f"the other's at {num(n, 'beyond.replication_leftover'):.2f}. The weights, "
            f"map and terms z-scored: Gria1 {num(n, 'beyond.weight_Gria1'):+.2f}, "
            f"density {num(n, 'beyond.weight_density'):+.2f}. On the "
            f"{n['beyond.calibration_structures']} structures of the calibration, nano "
            f"leaves {pct(n, 'beyond.nano_calibration_left')}, the floor "
            f"{pct(n, 'beyond.floor')}: nano minus the floor "
            f"{points_with_interval(n, 'beyond.nano_minus_floor')}. " + checks_text(n),
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
            f"{n['overview.adults_cacng8_above_gria1']} of {len(ADULTS)} adults, by "
            f"{num(n, 'gene_ranking.gap'):+.2f} "
            f"({num(n, 'gene_ranking.gap_boot_lo'):+.2f} to "
            f"{num(n, 'gene_ranking.gap_boot_hi'):+.2f} over resampled adults): "
            f"{gap_words(n)}; against unrelated maps p "
            f"{num(n, 'gene_ranking.gap_p_unrelated'):.2f}.",
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
        f"main model, {len(plain)} take more of it than 95% of their plain surrogates, "
        f"the null fixed first ({', '.join(plain) or 'none'}); {len(alike)} more than "
        "95% of maps that relate to the model as they do, a null added after seeing "
        f"({', '.join(alike) or 'none'}). Cacng8 takes "
        f"{100 * num(n, 'top_genes.taken_Cacng8'):.1f} points of the reproducible map (p "
        f"{num(n, 'top_genes.p_taken_Cacng8'):.3f} against plain surrogates; "
        f"{num(n, 'top_genes.p_taken_alike_Cacng8'):.3f} against maps alike)."
    )


def kinds_walk(n: dict[str, str]) -> dict[str, tuple[str, str]]:
    """What to look at and what to take from figures 10 and 11."""
    october = "gene_sets.localisation_5_Octobers_controls_positive_control_p"
    past, tested = sets_past(n)
    names = n["top_genes.family_leftover_pass_within_list"].split()
    signs = n["top_genes.family_leftover_pass_within_rho"].split()
    within_family = ", ".join(
        f"{g} ({float(r):+.2f})" for g, r in zip(names, signs, strict=True)
    )
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
            f"{n['beyond.leftover_genes']} by rho), a re-test: its p against the "
            "leftover of the four-subunit model was seen before it was named, and its "
            "p against the first version's before the second was decided. The family "
            f"({n['top_genes.family_tested']} of {n['top_genes.family_listed']} members "
            "with a usable map): median "
            f"{num(n, 'top_genes.family_leftover_median'):+.3f}, p "
            f"{num(n, 'top_genes.family_leftover_p'):.3f} against the surrogates "
            "(without Cacng8 p "
            f"{num(n, 'top_genes.family_leftover_without_cacng8_p'):.3f}); "
            f"{num(n, 'top_genes.family_leftover_controls_difference'):+.3f} above "
            f"{n['top_genes.family_controls']} matched postsynaptic genes, p "
            f"{num(n, 'top_genes.family_leftover_controls_p'):.2f} (against controls "
            "outside the model's terms p "
            f"{num(n, 'top_genes.family_leftover_outside_p'):.2f}); past BH within the "
            f"family: {within_family or 'none'}. Every gene: "
            f"{n['beyond.leftover_genes_pass']} past BH "
            f"({n['beyond.leftover_genes_pass_negative']} of them below zero, the lowest "
            f"{', '.join(n['beyond.leftover_genes_pass_lowest'].split()[:5])}); the "
            "highest rho, "
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
            f"{ordinal(n['overview.robustness_gria1_rank_min'])} to "
            f"{ordinal(n['overview.robustness_gria1_rank_max'])} of P9's "
            f"{n['gene_ranking.nano_p9_genes']} "
            f"({ordinal(n['top_genes.rank_variants_min_Gria1'])} to "
            f"{ordinal(n['top_genes.rank_variants_max_Gria1'])} of all genes over the "
            "variants that hold every gene), and the gap stays between "
            f"{num(n, 'overview.robustness_gap_min'):+.2f} and "
            f"{num(n, 'overview.robustness_gap_max'):+.2f}, around the edge of the "
            "primary's band.",
        ),
        "synaptome": (
            "A, the structures of the fit with a measured PSD95 density, by division; "
            "B, the PSD95 density against the density term; C, what the model leaves "
            "with the density genes and with PSD95 as its density term on the "
            "structures PSD95 covers, each with nano minus its own floor.",
            f"PSD95 puncta (Zhu et al. 2018) are measured in "
            f"{n['synaptome.fit_measured']} of the {n['synaptome.fit']} structures of "
            f"the fit ({pct(n, 'synaptome.fit_share')}). There the density genes "
            "order the structures as PSD95 does at "
            f"{num(n, 'synaptome.rho_fit_density'):.2f} (optimistic: they were chosen "
            "by it; held out "
            f"{num(n, 'density_markers.chosen_held_out'):.2f}), and the two hemispheres "
            f"of the one mouse agree at {num(n, 'synaptome.hemispheres_rho'):.2f}. On "
            f"the {n['beyond.check_structures_psd95']} structures with PSD95 and Gria1, "
            f"PSD95 alone predicts {pct(n, 'beyond.check_density_alone_psd95')} of the "
            "map, the density genes "
            f"{pct(n, 'beyond.check_density_alone_psd95_main')}; the model leaves "
            f"{pct(n, 'beyond.check_left_psd95')} with PSD95 (nano minus its floor "
            f"{check_margin(n, 'psd95')}), {pct(n, 'beyond.check_left_psd95_main')} "
            f"with the genes ({check_margin(n, 'psd95_main')}).",
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


def part1_standing(n: dict[str, str]) -> str:
    """Where part 1 stands, in words that follow the numbers: a clause, no capital."""
    lo = num(n, "beyond.nano_minus_floor_lo")
    floor = num(n, "beyond.floor")
    curved_lo = num(n, "beyond.check_nano_minus_floor_curved_lo")
    if lo >= 0.05:
        above = f"well above the floor ({floor:.0%}) on the same structures"
    elif lo > 0:
        above = (
            f"above the floor ({floor:.0%}) on the same structures, though the "
            "interval of the difference starts near zero"
        )
    elif num(n, "beyond.nano_minus_floor") > 0:
        above = (
            f"above the floor ({floor:.0%}) on the same structures at its point "
            f"value only, the interval of the difference reaching {lo:+.0%}"
        )
    else:
        above = f"no more than the floor ({floor:.0%}) on the same structures"
    if curved_lo > 0:
        bound = "and so is the curved model's against its own floor"
    else:
        bound = (
            "while the curved model's, the bound, has an interval against its own "
            f"floor that reaches {curved_lo:+.0%}"
        )
    return f"the leftover is reproducible, {above}, {bound}"


def part1_verdict(n: dict[str, str]) -> str:
    """Where part 1 stands, as the box of the overview figure says it."""
    return f"Where it stands: {part1_standing(n)}."


def family_text(n: dict[str, str]) -> str:
    """What the AMPA receptor complex family does with the leftover, in words.

    The words follow its two group tests.
    """
    spatial = num(n, "top_genes.family_leftover_p")
    matched = num(n, "top_genes.family_leftover_controls_p")
    if spatial < ALPHA and matched >= ALPHA:
        return (
            "the AMPA receptor complex family follows it beyond the surrogates but "
            "no more than postsynaptic genes of the same expression"
        )
    if spatial < ALPHA:
        return (
            "the AMPA receptor complex family follows it beyond the surrogates and "
            "beyond matched postsynaptic genes"
        )
    if matched < ALPHA:
        return (
            "the AMPA receptor complex family follows it beyond matched postsynaptic "
            "genes, not beyond the surrogates"
        )
    return "the AMPA receptor complex family as a group does not"


def named_tests_text(n: dict[str, str]) -> str:
    """The tests named in advance for the leftover, in words that follow the numbers."""
    cacng8 = num(n, "top_genes.leftover_p_Cacng8")
    first = "follows it" if cacng8 < ALPHA else "does not follow it"
    return (
        f"Against what Gria1 and synapse density leave, Cacng8 {first} (p "
        f"{cacng8:.4f}), a re-test of what was seen on 8 October; {family_text(n)} (p "
        f"{num(n, 'top_genes.family_leftover_p'):.3f} and "
        f"{num(n, 'top_genes.family_leftover_controls_p'):.3f})"
    )


def part2_verdict(n: dict[str, str]) -> str:
    """Where part 2 stands, in words that follow the numbers."""
    gap_p = num(n, "gene_ranking.gap_p")
    past, _ = sets_past(n)
    local_p = num(n, "gene_sets.localisation_matched_controls_p")
    if gap_p < ALPHA or past or local_p < ALPHA:
        return (
            "Where it stands: some of the genes that set surface receptor single the "
            f"map out beyond the null; read {figure_ref('cacng8_gria1')} and "
            f"{figure_ref('gene_kinds')} for which."
        )
    return (
        "Where it stands: consistent with the surface-fraction reading, not singling "
        "it out: the genes that follow the map most are maps much like Gria1 and "
        "synapse density, Cacng8's lead over Gria1 sits at the edge of its null, and "
        "the localisation genes do no better than matched controls."
    )


def gene_meanings(n: dict[str, str]) -> dict[str, str]:
    """What the figures of the genes mean, each verdict following its numbers."""
    gap_p = num(n, "gene_ranking.gap_p")
    pre_p = num(n, "gene_sets.postsynaptic_against_presynaptic_p_spatial")
    glia_p = num(n, "gene_sets.postsynaptic_against_glia_p_spatial")
    local_p = num(n, "gene_sets.localisation_matched_controls_p")
    found = num(n, "gene_sets.localisation_detectable_labels")
    follow = n["top_genes.top_follow_leftover_list"].split()
    within = [
        g
        for g in ("Cacng8", "Dlg2", "Gria1")
        if num(n, f"divisions.p_within_spatial_{g}") < ALPHA
    ]
    out = {
        "one_comparison": "A whole-brain rho mixes fine agreement with the contrast "
        "between divisions; the next figures separate them.",
    }
    if CALIBRATED[0] <= num(n, "spatial_null.fpr_spatial_map") <= CALIBRATED[1]:
        out["spatial_null"] = (
            "The spatial p keeps about 5% false positives where the ordinary one "
            "does not; every p of the gene analyses is a spatial p."
        )
    else:
        out["spatial_null"] = "The spatial p is not calibrated; no claim rests on it."
    out["top_genes"] = (
        "The genes that follow the map most are maps much like Gria1 and synapse "
        "density, so they describe the map rather than add evidence"
    )
    if follow:
        out["top_genes"] += (
            f", and of the {len(follow)} that also follow part of what the model leaves "
            "only Cacng8 was named for it"
        )
    out["top_genes"] += "."
    if gap_p < ALPHA:
        out["cacng8_gria1"] = "The map follows Cacng8 more closely than Gria1."
    else:
        out["cacng8_gria1"] = (
            "These maps do not settle whether the nano map follows Cacng8 more closely "
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
    if pre_p < ALPHA and glia_p < ALPHA:
        kind = "the map is postsynaptic-like by the criterion named in advance"
    elif pre_p < ALPHA or glia_p < ALPHA:
        kind = "one of the two contrasts named in advance passes, not both"
    else:
        kind = "neither contrast named in advance passes"
    if local_p >= ALPHA:
        out["gene_kinds"] = (
            f"The genes that set surface receptor do not stand out: {kind}, and no "
            f"localisation advantage larger than {found:+.2f} exists beyond the null."
        )
    else:
        out["gene_kinds"] = (
            f"The localisation genes stand above their matched controls; {kind}."
        )
    out["leftover"] = named_tests_text(n)
    spatial = num(n, "top_genes.family_leftover_p")
    matched = num(n, "top_genes.family_leftover_controls_p")
    if spatial < ALPHA <= matched:
        out["leftover"] += (
            ", so what the family shares with the leftover is postsynaptic, not "
            "particular to the AMPA receptor complex"
        )
    n_pass = int(num(n, "beyond.leftover_genes_pass"))
    if n_pass == 0:
        out["leftover"] += ", and no gene passes once every gene is corrected for"
    else:
        out["leftover"] += (
            f"; {n_pass} genes pass once every gene is corrected for, "
            f"{n['beyond.leftover_genes_pass_negative']} of them with a negative rho "
            "(nano below prediction where they are high), an exploratory finding"
        )
    out["leftover"] += "."
    return out


def other_meanings(n: dict[str, str]) -> dict[str, str]:
    """What the inputs, part 1, the controls, the limit and April's headline mean."""
    out = {
        "overview": "Every number on it is this run's; the story in words is "
        "docs/ISH_ANALYSIS.md.",
        "structures": "Every comparison runs on structures every adult measures; the "
        "declared reference shifts each brain's zref without reordering the map.",
        "genes": "A gene measured once is only as good as one Allen mouse; a low rho "
        "of an unreliable gene says little.",
        "beyond": upper_first(part1_standing(n))
        + "; it is what the model does not predict, not a measurement of the surface "
        "fraction.",
        "beyond_where": "The departure from prediction sits in particular structures, "
        "steadily across the adults; a claim about one structure needs its own null.",
        "robustness": "Cacng8 first of P9's genes holds under every choice; Gria1's "
        "rank moves with them, so it is quoted with its null, not as a place.",
    }
    auto_p = num(n, "gene_ranking.auto_p_Gria1")
    above = int(num(n, "overview.adults_nano_above_auto_Gria1"))
    if auto_p >= ALPHA and above == len(ADULTS):
        out["autofluorescence"] = (
            "Gria1 and Cacng8 follow the label, not the tissue, in every adult; the "
            "tissue has a gene pattern of its own, a different one."
        )
    else:
        out["autofluorescence"] = "The tissue shares part of the label's ranking."
    psd95 = num(n, "beyond.check_left_psd95")
    genes = num(n, "beyond.check_left_psd95_main")
    if psd95 > genes:
        out["synaptome"] = (
            "The measured density covers about half the fit; where it exists, the "
            "model leaves more with it than with the density genes, so the genes make "
            "part 1 no easier to pass."
        )
    else:
        out["synaptome"] = (
            "The measured density covers about half the fit; where it exists, the "
            "model leaves no more with it than with the density genes."
        )
    if num(n, "green_channel.rho_sep_auto_min") > 0.5:
        out["green_channel"] = (
            "The green channel reports mostly the tissue in every adult: these "
            "brains cannot measure the surface fraction; a total-GluA1 stain would."
        )
    else:
        out["green_channel"] = "The green channel does not simply follow the tissue."
    if num(n, "overview.anova_p_spatial") >= ALPHA:
        out["april_headline"] = (
            "What reproduces from April is the gene order, Cacng8 first; the "
            "difference between its categories does not survive the null."
        )
    else:
        out["april_headline"] = "April's categories differ beyond the null."
    return out


def index_head() -> str:
    """The head of figures/README.md: what it is, the argument, how the figures go."""
    first, where = FIGURE_NUMBERS["beyond"], FIGURE_NUMBERS["beyond_where"]
    genes, leftover = FIGURE_NUMBERS["one_comparison"], FIGURE_NUMBERS["leftover"]
    return f"""# The ISH analysis, figure by figure

Written by `run_ish_overview.py` from the tables of this run, so every number below
is this run's. The story, with what each result means and does not mean, is
`docs/ISH_ANALYSIS.md` in the code repository.

The figures follow one argument in two parts. Part 1 (figures {first} and
{where}): the adult nano map across structures is not fully
explained by Gria1 expression and synapse density; a reproducible part is left
over, and set against what Allen-to-Allen mismatch alone leaves (the calibration
floor). Part 2 (figures {genes} to
{leftover}): it is therefore something else, and the reading the data
support is the surface fraction
of the receptor (trafficking, scaffolding); the genes that follow the map and its
leftover are the corroboration, Cacng8 (TARP gamma-8) first. The surface fraction
is an interpretation, not a measurement: a total-GluA1 stain on the same brains
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
- `top_genes/<gene>.png`: every gene of figures {FIGURE_NUMBERS["top_genes"]} and
  {FIGURE_NUMBERS["ampa_family"]} (the genes past
  the map's null, Gria1, Cacng8 and the AMPA receptor complex family): its map,
  against the map and against the leftover, what it takes of the leftover against
  maps of its smoothness, what it is like, its numbers and GO terms
  (`run_ish_top_genes.py --sheets`; PNG only).
"""


def figure_index(n: dict[str, str]) -> str:
    """figures/README.md: the guided walk, a section per main figure.

    Each with its question, what to look at and what to take from it, and its
    detailed versions.
    """
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
                f"### {FIGURE_NUMBERS[key]}. {QUESTIONS[key]}\n",
                f"![{name}]({name})\n",
                f"**Look at.** {look}\n",
                f"**Numbers.** {take}\n",
                f"**Takeaway.** {meaning[key]}\n",
            ]
            details = SUPPLEMENTS.get(key, ())
            for detail, adds in details:
                other = figure_file(detail)
                lines.append(f"- [{FIGURE_NUMBERS[detail]}]({other}) in detail: {adds}.")
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
        claim="Part 1. The map is not fully explained by Gria1 expression and "
        "synapse density: a reproducible part is left over",
        numbers=[
            (
                pct(n, "beyond.share_model"),
                "of the map's reproducible pattern is predicted by Gria1 and synapse "
                "density, on structures the fit has not seen (Gria1 alone "
                f"{pct(n, 'beyond.share_abundance')}) "
                f"({figure_ref('beyond')})",
            ),
            (
                f"{pct(n, 'beyond.left')} left",
                f"{pct(n, 'beyond.left_lo')} to {pct(n, 'beyond.left_hi')} over "
                "structures; the two halves of the cohort leave the same leftover "
                f"(rho {num(n, 'beyond.replication_leftover'):.2f})",
            ),
            (
                f"{pct(n, 'beyond.nano_calibration_left')} vs {pct(n, 'beyond.floor')}",
                f"on the {n['beyond.calibration_structures']} structures of the "
                "calibration, what the model leaves of nano, and of a map made only of "
                "Gria1 and density (the floor): nano minus it "
                f"{points_with_interval(n, 'beyond.nano_minus_floor')}; curved, against "
                "its own floor, "
                f"{points_with_interval(n, 'beyond.check_nano_minus_floor_curved')}",
            ),
        ],
        stands=part1_verdict(n),
    )


def cacng8_follows(n: dict[str, str]) -> str:
    """Whether Cacng8 follows the leftover, as tier 1 says, and that it is a re-test."""
    if num(n, "top_genes.leftover_p_Cacng8") < ALPHA:
        verb = "follows"
    else:
        verb = "does not follow"
    return (
        f"Cacng8 {verb} what Gria1 and synapse density leave, a re-test of what was "
        "seen on 8 October"
    )


def taken_words(n: dict[str, str]) -> str:
    """What Cacng8 takes of the leftover against its plain surrogates, then maps alike.

    The plain surrogates are the null fixed first, quoted first; maps alike to the
    model, added after seeing, are the secondary line.
    """
    plain = num(n, "top_genes.p_taken_Cacng8")
    alike = num(n, "top_genes.p_taken_alike_Cacng8")
    first = "more than" if plain < ALPHA else "not more than"
    second = "more than" if alike < ALPHA else "not more than"
    return (
        f"{first} its plain surrogates take (p {plain:.3f}, the null fixed first); "
        f"{second} maps alike to the model (p {alike:.3f}, a null added after seeing)"
    )


def part2_content(n: dict[str, str]) -> dict:
    """The box of part 2: its claim, three key numbers, where it stands."""
    within = ""
    if num(n, "top_genes.p_within_Cacng8") < ALPHA:
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
                f"{FIGURE_NUMBERS['top_genes']}, {FIGURE_NUMBERS['cacng8_gria1']})",
            ),
            (
                f"p {num(n, 'top_genes.leftover_p_Cacng8'):.4f}",
                f"{cacng8_follows(n)} (figure {FIGURE_NUMBERS['leftover']}); added to "
                f"the model it takes {100 * num(n, 'top_genes.taken_Cacng8'):.1f} points "
                f"of the map, {taken_words(n)} (figure {FIGURE_NUMBERS['top_genes']})",
            ),
            (
                f"p {num(n, 'top_genes.family_leftover_p'):.3f}",
                f"{upper_first(family_text(n))} (p "
                f"{num(n, 'top_genes.family_leftover_controls_p'):.2f}); without "
                "Cacng8 the group's p is "
                f"{num(n, 'top_genes.family_leftover_without_cacng8_p'):.3f}",
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
        "structure makes, or how many synapses it has? Allen's ISH maps measure Gria1's "
        "mRNA, and synapse density through the mRNA of postsynaptic genes chosen by a "
        "measured synapse density.",
        parts=[part1_content(n), part2_content(n)],
        limit="The surface fraction is an interpretation, not a measurement. The "
        "green channel follows autofluorescence in every adult (rho "
        f"{num(n, 'green_channel.rho_sep_auto_min'):.2f} to "
        f"{num(n, 'green_channel.rho_sep_auto_max'):.2f}), so it cannot give total "
        "receptor; a total-GluA1 stain on some of the same brains would measure it "
        f"(figure {FIGURE_NUMBERS['green_channel']}).",
    )


def figure_map() -> list[tuple[str, list[tuple[str, str, str]]]]:
    """The groups of the overview's figure map, by part of the walk.

    Each main figure's number, question and detailed versions, the overview itself
    left out.
    """
    groups = []
    for part, keys in PARTS[1:]:
        figures = []
        for key in keys:
            details = ", ".join(FIGURE_NUMBERS[d] for d, _ in SUPPLEMENTS.get(key, ()))
            figures.append((FIGURE_NUMBERS[key], QUESTIONS[key], details))
        groups.append((part.split(":")[0], figures))
    return groups
