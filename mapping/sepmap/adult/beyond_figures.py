"""The numbers and the guided figures of part 1: 03, 03s1, 03s2, 04 and 11s1.

adult.beyond_density, adult.beyond_controls, adult.beyond_calibration,
adult.beyond_checks and adult.beyond_regression write the tables. This module adds
the intervals of the partition and the weights, over structures and over spatial
blocks (the jackknives of adult.beyond_calibration), gathers the numbers the text
quotes, and has adult.plotting draw the figures of part 1 (step 22 draws 03s3, the
choice of the density genes, and step 21 03s4, the measured synapse density):

    03   what Gria1 and synapse density predict of the map: the map against the
         main model's prediction, the reproducible map in four parts (Gria1 only,
         shared, density only, left) with the leftover's range, nano and the
         calibration floor on the same structures with their difference, and the
         density term against the measured PSD95 density that chose its genes
         (adult.density_markers)
    03s1 the check rows and the main model under other folds, each with what it
         leaves and nano minus its own floor; then the map against each term,
         every model's share, the four parts with their intervals and the
         weights, the calibration and the replication
    03s2 the seven controls of adult.beyond_controls, one panel each
    04   where the leftover lives: map, prediction and leftover on three coronal
         planes, the structures with the largest leftovers, and the leftover by
         division
    11s1 every gene of the gene table and every gene set against the leftover, with
         the leftover's spatial null, and each gene's rho with the leftover against
         its rho with the autofluorescence map

The leftover and the floor are compared on the same structures, the calibration's,
and their difference is resampled with them (beyond_calibration.paired_jackknife);
every check row has its own (adult.beyond_checks). The noise null of the
replication shuffles which structure is which in one half's leftover: what two
unrelated leftovers would give. It is quoted for scale only, since the leftover of a
reproducible map was always going to replicate (beyond_density.implied_replication).

Writes, under adult_v2/ish_analysis/ in the data root:

    beyond/jackknife.csv                per subsample: the ceiling, each model's
                                        share, the four parts and the weights
    beyond/jackknife_blocks.csv         the same with one spatial block left out at
                                        a time
    beyond/numbers_for_the_caption.txt  the figures' numbers as sentences
    tables/numbers_beyond.csv           the numbers of analysis 4, for the text
    figures/03_beyond.png               (and .eps, as every figure)
    figures/03s1_beyond_budget.png
    figures/03s2_beyond_controls.png
    figures/04_beyond_where.png
    figures/11s1_leftover_genes.png

Run by run_beyond_figures.py.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap import structures
from sepmap.adult import (
    beyond_calibration,
    beyond_checks,
    beyond_density,
    beyond_regression,
    density_markers,
)
from sepmap.adult import plotting as adult_plotting
from sepmap.adult.beyond_controls import (
    CONTROLS,
    FOLDS,
    GENE_SPACE,
    GENE_SPACE_CALIBRATION,
    GENE_SPACE_SUMMARY,
    READINGS,
    adult_pair_agreements,
    curvature_scores,
    group_leftovers,
    mean_sizes,
)
from sepmap.adult.profiles import ADULTS
from sepmap.config import SETTINGS
from sepmap.ish import gene_ranking
from sepmap.ish.figure_index import figure_file, figure_path
from sepmap.ish.gene_sets import LEFTOVER_GENE
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.ish.numbers import points as points_text
from sepmap.ish.plotting.shared import group_of
from sepmap.ish.spatial_null import ALPHA
from sepmap.plotting import paint

# the permutations of the noise null and the grey of the bars; the BH level; the grey
# of the genes' bars; the thresholds of the controls' verdicts, which 03s2 draws
BEYOND_FIGURES = SETTINGS["beyond_figures"]
ISH_ANALYSIS = SETTINGS["ish_analysis"]
ISH_FIGURES = SETTINGS["ish_figures"]
BEYOND_CONTROLS = SETTINGS["beyond_controls"]

JACKKNIFE = beyond_density.OUT / "jackknife.csv"
JACKKNIFE_BLOCKS = beyond_density.OUT / "jackknife_blocks.csv"
CAPTION = beyond_density.OUT / "numbers_for_the_caption.txt"
NUMBERS = numbers_path("beyond")

# the percentiles of a 95% interval
INTERVAL = (2.5, 97.5)

# the models whose share has an interval from the jackknife, by their key in
# variance_partition.csv
INTERVAL_SHARES = ("abundance", "density", "model")


def interval(values: np.ndarray | list[float]) -> tuple[float, float]:
    """The 95% interval of `values`."""
    lo, hi = np.percentile(np.asarray(values, dtype=float), INTERVAL)
    return float(lo), float(hi)


def noise_band(
    nano: np.ndarray,
    columns: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> tuple[float, float]:
    """95% of the agreement of two unrelated leftovers (first split, one shuffled)."""
    a, b = splits[0]
    ra = beyond_density.residual(beyond_density.half_map(nano, a), columns)
    rb = beyond_density.residual(beyond_density.half_map(nano, b), columns)
    null = [
        spearmanr(ra, rng.permutation(rb)).statistic
        for _ in range(BEYOND_FIGURES["n_perm"])
    ]
    return interval(null)


# ===== The numbers table and the caption =====


def structure_number_rows(n: dict) -> list[tuple]:
    """The structures of the fit: how many, and what the rule of the structures loses."""
    rows = [
        ("structures", n["n_structures"], "structures of the main model's fit"),
        ("structures_declared", n["n_declared"], "declared structures"),
        ("adults", n["n_adults"], "adults, naive and RWS pooled"),
        ("main_abundance", n["abundance"], "the main model's abundance term"),
        (
            "main_density_genes",
            " ".join(n["density_genes"]),
            "the genes whose mean rank is the density term",
        ),
        (
            "psd95_measured",
            n["psd95_measured"],
            "structures of the fit with a measured PSD95 density",
        ),
    ]
    for gene, divisions in n["lost"].items():
        rows.append(
            (
                f"lost_{gene}",
                sum(divisions.values()),
                f"declared structures left out: {gene} not measured",
            )
        )
        for division, count in divisions.items():
            rows.append((f"lost_{gene}_{division}", count, f"of them, in {division}"))
        rows.append(
            (
                f"lost_{gene}_structures",
                " ".join(n["lost_acronyms"][gene]),
                "those structures, by acronym",
            )
        )
    return rows


def interval_rows(
    key: str, over_structures: tuple[float, float], over_blocks: tuple[float, float]
) -> list[tuple]:
    """The rows of a value's two 95% intervals: <key>_lo, _hi, _blocks_lo, _blocks_hi."""
    return [
        (f"{key}_lo", round(over_structures[0], 6), "2.5%, jackknife over structures"),
        (f"{key}_hi", round(over_structures[1], 6), "97.5%, jackknife over structures"),
        (f"{key}_blocks_lo", round(over_blocks[0], 6), "2.5%, over spatial blocks"),
        (f"{key}_blocks_hi", round(over_blocks[1], 6), "97.5%, over spatial blocks"),
    ]


def model_number_rows(n: dict) -> list[tuple]:
    """The main model's numbers: ceiling, shares, parts, weights and replication."""
    rows = [
        (
            "half_agreement",
            round(n["half"], 4),
            f"half-cohort maps agree ({n['n_splits']} splits)",
        ),
        ("ceiling", round(n["ceiling"], 4), "Spearman-Brown: reproducible share"),
        ("ceiling_lo", round(n["ceiling_ci"][0], 6), "2.5%, jackknife over structures"),
        ("ceiling_hi", round(n["ceiling_ci"][1], 6), "97.5%, jackknife over structures"),
        (
            "jackknife_blocks",
            n["n_blocks"],
            "spatial blocks of the jackknife over blocks",
        ),
    ]
    for key, r in n["partition_table"].iterrows():
        rows.append(
            (
                f"share_{key}",
                round(r["share_of_ceiling"], 4),
                f"{r['model']}: CV R2 / ceiling",
            )
        )
        if key in n["shares_ci"]:
            rows += interval_rows(
                f"share_{key}", n["shares_ci"][key], n["shares_ci_blocks"][key]
            )
    for part, value in n["parts"].items():
        rows.append((f"part_{part}", round(value, 4), f"the reproducible map: {part}"))
        rows += interval_rows(
            f"part_{part}", n["parts_ci"][part], n["parts_ci_blocks"][part]
        )
    for term, value in n["weights"].items():
        rows.append(
            (f"weight_{term}", round(value, 4), f"{term}: weight, map and terms z-scored")
        )
        rows += interval_rows(
            f"weight_{term}", n["weights_ci"][term], n["weights_ci_blocks"][term]
        )
        rows.append(
            (
                f"rho_map_{term}",
                round(n["rho_terms"][term], 4),
                f"{term}: Spearman with the map",
            )
        )
    rows += [
        ("rho_terms", round(n["rho_between"], 4), "Gria1 against the density term"),
        ("left", round(n["left"], 4), "share of the reproducible map left"),
    ]
    rows += interval_rows("left", n["left_ci"], n["left_ci_blocks"])
    rows += [
        ("cv_r2", round(n["cv_r2"], 4), "the main model's CV R2"),
        ("in_sample_r2", round(n["r2"], 4), "the main model's in-sample R2"),
        ("terms", n["n_terms"], "the main model's columns, the intercept included"),
        ("replication_map", round(n["rep_map"], 4), "half-cohort maps, mean over splits"),
        ("replication_leftover", round(n["rep_left"], 4), "their leftovers, mean"),
        ("replication_leftover_min", round(n["rep_left_min"], 4), "lowest split"),
        ("replication_implied", round(n["implied"], 4), "implied by ceiling and fit"),
        ("noise_lo", round(n["noise"][0], 4), "2.5% of unrelated leftovers"),
        ("noise_hi", round(n["noise"][1], 4), "97.5% of unrelated leftovers"),
    ]
    return rows


def calibration_number_rows(n: dict) -> list[tuple]:
    """The calibration's numbers: the floor, and nano beside it."""
    rows = [
        ("calibration_structures", n["cal_n"], "structures of the calibration"),
        ("floor", round(n["floor"]["left_median"], 4), "Gria1-and-density map: left"),
        ("floor_lo", round(n["floor"]["left_lo"], 4), "lowest draw"),
        ("floor_hi", round(n["floor"]["left_hi"], 4), "highest draw"),
        (
            "floor_replication",
            round(n["floor"]["replication_median"], 4),
            "Gria1-and-density map: leftover replicates",
        ),
        (
            "floor_own_half",
            round(n["own_half"], 4),
            "that map read with its own half: left, the animals' noise alone",
        ),
        (
            "floor_two_halves",
            round(n["floor"]["left_median"] - n["own_half"], 4),
            "the floor less that: two Allen halves disagreeing",
        ),
        (
            "calibration_once",
            " ".join(n["once"]),
            "the main model's genes the same in both halves (measured once)",
        ),
    ]
    for h, value in n["nano_cal"].items():
        rows.append(
            (
                f"nano_calibration_left_{h}",
                round(value, 4),
                f"nano on the calibration's structures, predictors from {h}",
            )
        )
    rows.append(
        ("nano_calibration_left", round(n["nano_cal_left"], 4), "nano, mean of A and B")
    )
    rows += interval_rows("nano_calibration", n["nano_cal_ci"], n["nano_cal_ci_blocks"])
    rows += [
        (
            "nano_floor_ratio",
            round(n["nano_cal_left"] / n["floor"]["left_median"], 2),
            "nano's share left over the floor's, same structures",
        ),
        ("nano_minus_floor", round(n["minus_floor"], 4), "nano minus the floor"),
    ]
    rows += interval_rows(
        "nano_minus_floor", n["minus_floor_ci"], n["minus_floor_ci_blocks"]
    )
    rows += [
        ("blocks_nano_left", round(n["blocks"]["nano"], 4), "spatial blocks: nano"),
        ("blocks_floor", round(n["blocks"]["floor"], 4), "spatial blocks: the floor"),
    ]
    return rows


def check_number_rows(n: dict) -> list[tuple]:
    """The other folds and the check rows, each with its floor where it has one."""
    rows = []
    for r in n["folds"].itertuples():
        if r.kind == "main":
            continue
        rows.append((f"fold_left_{r.key}", round(r.left, 4), r.folds))
        if np.isfinite(r.lo):
            rows.append((f"fold_left_lo_{r.key}", round(r.lo, 4), "2.5%"))
            rows.append((f"fold_left_hi_{r.key}", round(r.hi, 4), "97.5%"))
    for r in n["checks"].itertuples():
        rows += [
            (f"check_structures_{r.key}", int(r.n_structures), r.check),
            (f"check_left_{r.key}", round(r.left, 4), f"{r.check}: share left"),
            (
                f"check_density_alone_{r.key}",
                round(r.density_alone, 4),
                f"{r.check}: density alone",
            ),
            (f"check_cal_structures_{r.key}", int(r.cal_structures), "its calibration"),
            (f"check_nano_cal_{r.key}", round(r.nano_cal_left, 4), "nano there"),
            (f"check_floor_{r.key}", round(r.floor, 4), "its floor"),
            (
                f"check_nano_minus_floor_{r.key}",
                round(r.nano_minus_floor, 4),
                "nano minus its floor",
            ),
            (
                f"check_nano_minus_floor_{r.key}_lo",
                round(r.nano_minus_floor_lo, 6),
                "2.5%, jackknife",
            ),
            (
                f"check_nano_minus_floor_{r.key}_hi",
                round(r.nano_minus_floor_hi, 6),
                "97.5%, jackknife",
            ),
            (
                f"check_nano_minus_floor_{r.key}_blocks_lo",
                round(r.nano_minus_floor_blocks_lo, 6),
                "2.5%, over spatial blocks",
            ),
            (
                f"check_nano_minus_floor_{r.key}_blocks_hi",
                round(r.nano_minus_floor_blocks_hi, 6),
                "97.5%, over spatial blocks",
            ),
            (f"check_once_{r.key}", r.once, "measured once, the same in both halves"),
        ]
        if np.isfinite(r.psd_genes):
            rows += [
                (f"check_psd_genes_{r.key}", int(r.psd_genes), "genes behind psd_pc1"),
                (
                    f"check_psd_genes_calibration_{r.key}",
                    int(r.psd_genes_calibration),
                    "those of its calibration: two halves each, one list for both",
                ),
            ]
    checks = n["checks"]
    rows += [
        ("check_rows", len(checks), "the check rows, each with its own floor"),
        (
            "checks_above_point",
            int((checks["nano_minus_floor"] > 0).sum()),
            f"of the {len(checks)} check rows, above their floor at the point value",
        ),
        (
            "checks_above_interval",
            int((checks["nano_minus_floor_lo"] > 0).sum()),
            "above it over the whole 95% over structures",
        ),
        (
            "checks_above_blocks",
            int((checks["nano_minus_floor_blocks_lo"] > 0).sum()),
            "above it over the whole 95% over spatial blocks",
        ),
        (
            "checks_reaching_zero_blocks",
            " ".join(checks.loc[checks["nano_minus_floor_blocks_lo"] <= 0, "key"]),
            "the rows whose 95% over spatial blocks reaches zero",
        ),
    ]
    return rows


def control_number_rows(n: dict) -> list[tuple]:
    """The controls' numbers, and the leftover by division."""
    rows = [
        ("controls_passed", n["controls_passed"], "of the seven controls"),
        ("controls_failed", " ".join(n["controls_failed"]), "the controls that do not"),
        ("control_e_straight", round(n["curvature"][0], 4), "control E: CV R2 straight"),
        ("control_e_curved", round(n["curvature"][1], 4), "control E: CV R2 curved"),
        ("control_e_quintic", round(n["curvature"][2], 4), "control E: to fifth powers"),
        ("control_f_genes", n["f_genes"], "control F: genes measured everywhere"),
        ("control_f_components", n["f_k"], "control F: components picked most often"),
        ("control_f_share", round(n["f_share"], 4), "control F: nested CV R2 / ceiling"),
        ("control_f_replication", round(n["f_rep"], 4), "control F: leftover replicates"),
        ("control_f_cal_structures", n["f_cal_n"], "control F's floor: structures"),
        ("control_f_cal_genes", n["f_cal_genes"], "control F's floor: genes"),
        ("control_f_nano_lo", round(n["f_nano"][0], 4), "control F: nano left, lowest"),
        ("control_f_nano_hi", round(n["f_nano"][1], 4), "control F: nano left, highest"),
        ("control_f_floor_lo", round(n["f_floor"][0], 4), "control F's floor, lowest"),
        ("control_f_floor_hi", round(n["f_floor"][1], 4), "control F's floor, highest"),
        (
            "control_f_curve_peak",
            n["f_curve_peak"],
            "control F: components of the highest held-out R2 (max_pcs the most tried)",
        ),
        (
            "control_f_curve_plateau",
            n["f_plateau"],
            f"control F: fewest components within {n['f_plateau_r2']} of that R2",
        ),
        (
            "control_e_quintic_left",
            round(1 - n["curvature"][2] / n["ceiling"], 4),
            "control E: share left with fifth powers",
        ),
        (
            "leftover_between_divisions",
            round(n["between"], 4),
            "share of the leftover's variance between divisions",
        ),
        (
            "leftover_divisions_above",
            " ".join(n["divisions_above"]),
            "divisions above prediction as a whole (95% over resampled adults above 0)",
        ),
        (
            "leftover_divisions_below",
            " ".join(n["divisions_below"]),
            "divisions below prediction as a whole",
        ),
    ]
    for r in n["divisions"].itertuples():
        rows.append(
            (
                f"leftover_division_{r.division}",
                round(r.mean, 1),
                f"mean leftover in ranks, {r.n_structures} structures",
            )
        )
    return rows


def leftover_number_rows(n: dict) -> list[tuple]:
    """The leftover against the tissue, and every gene and gene set against it."""
    rows = [
        (
            "leftover_rho_autofluorescence",
            round(n["rho_auto"], 3),
            "the leftover's Spearman with the autofluorescence map",
        ),
        (
            "leftover_genes_against_auto",
            round(n["genes_against_auto"], 3),
            "over the genes, rho with the leftover against rho with autofluorescence",
        ),
        (
            "leftover_pass_auto_positive",
            n["pass_auto_positive"],
            "genes past BH on the leftover with a positive rho with autofluorescence",
        ),
        (
            "leftover_pass_auto_bh",
            n["pass_auto_bh"],
            "of them, past autofluorescence's own null after BH",
        ),
        ("leftover_genes", n["n_genes"], "genes against the leftover"),
        ("leftover_genes_pass", n["n_pass"], f"past its null, BH q < {n['q']}"),
        (
            "leftover_genes_pass_negative",
            n["n_pass_negative"],
            "of them with a negative rho: nano below prediction where they are high",
        ),
        (
            "leftover_genes_pass_positive",
            " ".join(n["pass_positive"]),
            "those with a positive rho",
        ),
        (
            "leftover_genes_pass_lowest",
            " ".join(n["pass_lowest"]),
            "the ten of the lowest rho among them",
        ),
        ("leftover_top_gene", n["top_gene"], "the gene of the highest rho"),
        ("leftover_top_rho", round(n["top_rho"], 3), "its rho"),
        ("leftover_top_p", round(n["top_p"], 4), "its spatial p"),
        ("leftover_top_q", round(n["top_q"], 4), "its q over all genes"),
        ("leftover_genes_p05", n["n_p05"], "genes with a spatial p below 0.05"),
        ("leftover_rho_Cacng8", round(n["cacng8"]["rho"], 3), "Cacng8 with the leftover"),
        ("leftover_p_Cacng8", round(n["cacng8"]["p_spatial"], 4), "its spatial p"),
        ("leftover_rank_Cacng8", int(n["cacng8"]["rank_all"]), "its rank"),
    ]
    for r in n["sets"].itertuples():
        key = r.gene_set.replace(" ", "_")
        rows.append(
            (
                f"leftover_set_{key}",
                round(r.median_rho, 3),
                f"median rho with the leftover, {r.n_genes} genes",
            )
        )
        if r.tested:
            rows += [
                (f"leftover_set_p_{key}", round(r.p_spatial, 4), "its spatial p"),
                (f"leftover_set_q_{key}", round(r.q, 4), "its q over the tested sets"),
            ]
    return rows


def numbers_table(n: dict) -> pd.DataFrame:
    """numbers_beyond.csv: the numbers of analysis 4 that the text quotes.

    The jackknife intervals keep six decimals: the text quotes them as whole
    percentages, and a value such as 0.154996 kept as 0.155 would round to 16%,
    not 15%.
    """
    rows = structure_number_rows(n) + model_number_rows(n)
    rows += calibration_number_rows(n) + check_number_rows(n)
    rows += control_number_rows(n) + leftover_number_rows(n)
    return numbers_frame(rows)


def points_with_interval(value: float, lo: float, hi: float) -> str:
    """A difference of shares in points, with its 95% interval."""
    return f"{points_text(value)} points (95% {points_text(lo)} to {points_text(hi)})"


def span(interval: tuple[float, float]) -> str:
    """A 95% interval of a difference in points: '+7 to +48'."""
    return f"{points_text(interval[0])} to {points_text(interval[1])}"


def lost_text(n: dict) -> str:
    """What the rule of the structures loses, in words: '40 without Gria1: HY 17, ...'."""
    parts = []
    for gene, divisions in n["lost"].items():
        by = ", ".join(f"{d} {k}" for d, k in divisions.items())
        parts.append(f"{sum(divisions.values())} without {gene}: {by}")
    return "; ".join(parts)


def caption_lines(n: dict) -> list[str]:
    """numbers_for_the_caption.txt: the figures' numbers as sentences."""
    floor = n["floor"]
    parts, ci = n["parts"], n["parts_ci"]
    margins = []
    for r in n["checks"].itertuples():
        blocks = span((r.nano_minus_floor_blocks_lo, r.nano_minus_floor_blocks_hi))
        margin = points_with_interval(
            r.nano_minus_floor, r.nano_minus_floor_lo, r.nano_minus_floor_hi
        )
        margins.append(f"{margin}, over spatial blocks {blocks}")
    above = ", ".join(n["divisions_above"])
    below = ", ".join(n["divisions_below"])
    listed = "; ".join(
        f"{r.check}, {r.n_structures} structures: {r.left:.0%} left, nano minus its "
        f"floor {margin}"
        for r, margin in zip(n["checks"].itertuples(), margins)
    )
    weights = ", ".join(
        f"{k} {v:+.2f} ({n['weights_ci'][k][0]:+.2f} to {n['weights_ci'][k][1]:+.2f})"
        for k, v in n["weights"].items()
    )
    return [
        f"Analysis 4, on {n['n_structures']} of the {n['n_declared']} declared "
        f"grey-matter structures, those where {n['abundance']} and every density gene "
        f"are measured (left out, {lost_text(n)}), and {n['n_adults']} adults. The main "
        f"model: nano rank ~ {n['abundance']} rank + synapse-density rank, two straight "
        f"terms; synapse density is the mean rank of {', '.join(n['density_genes'])}, "
        "chosen by their agreement with PSD95 punctum density without the map.",
        "",
        f"Ceiling. Two halves of the cohort agree at {n['half']:.3f} "
        f"({n['n_splits']} splits), so "
        f"{n['ceiling']:.1%} of the map is reproducible (Spearman-Brown; 95% over "
        f"structures, jackknife, {n['ceiling_ci'][0]:.1%} to {n['ceiling_ci'][1]:.1%}).",
        "",
        f"Partition. On structures the fit has not seen, {n['abundance']} alone "
        f"predicts {n['shares']['abundance']:.0%} of the reproducible map, synapse "
        f"density alone {n['shares']['density']:.0%}, both {n['shares']['model']:.0%}. "
        f"{n['abundance']} only {parts['gria1_only']:+.0%} ({ci['gria1_only'][0]:+.0%} "
        f"to {ci['gria1_only'][1]:+.0%}), shared {parts['shared']:+.0%} "
        f"({ci['shared'][0]:+.0%} to {ci['shared'][1]:+.0%}), density only "
        f"{parts['density_only']:+.0%} ({ci['density_only'][0]:+.0%} to "
        f"{ci['density_only'][1]:+.0%}); {n['left']:.0%} is not predicted by Gria1 "
        f"expression or synapse density (95% over structures, jackknife, "
        f"{n['left_ci'][0]:.0%} to {n['left_ci'][1]:.0%}; over spatial blocks "
        f"{n['left_ci_blocks'][0]:.0%} to {n['left_ci_blocks'][1]:.0%}). The weights, "
        f"map and terms z-scored: {weights}.",
        "",
        f"Calibration ({n['cal_n']} structures where both halves of the Allen "
        f"experiments measure every gene of the model). A map made only of Gria1 and "
        f"synapse density, predicted from the other half of the experiments, leaves "
        f"{floor['left_median']:.0%} ({floor['left_lo']:.0%} to {floor['left_hi']:.0%} "
        f"over {n['n_draws']} draws), its leftover replicating at "
        f"{floor['replication_median']:.2f}; read with its own half it leaves "
        f"{n['own_half']:.0%}, the animals' noise. So the floor holds the disagreement "
        "of two halves of the Allen experiments, where nano, read with one half, meets "
        "one half's error and the mismatch of Allen's P56 mice with these brains (age, "
        "strain, the grid, registration), which the floor lacks: it errs low only if "
        "that mismatch costs more than one half's error, about half of the "
        f"{100 * (floor['left_median'] - n['own_half']):.0f} points two halves cost. "
        f"The nano map, predicted the same way on the same structures, leaves "
        f"{n['nano_cal_left']:.0%} ({n['nano_cal_ci'][0]:.0%} to "
        f"{n['nano_cal_ci'][1]:.0%}); nano minus the floor "
        f"{points_with_interval(n['minus_floor'], *n['minus_floor_ci'])}, jackknife "
        "over the structures with both recomputed; leaving out one spatial block at a "
        f"time, {span(n['minus_floor_ci_blocks'])}.",
        "",
        f"Check rows, each one change to the main model on its own structures with its "
        f"own floor: {listed}.",
        "",
        f"Replication. The leftover of one half-cohort agrees with the other's at "
        f"{n['rep_left']:.3f} (lowest split {n['rep_left_min']:.3f}); the map's "
        f"reliability and the fit alone imply {n['implied']:.3f}; two unrelated "
        f"leftovers fall between {n['noise'][0]:+.2f} and {n['noise'][1]:+.2f}.",
        "",
        f"Controls. {n['controls_passed']} of {n['n_controls']} pass"
        + (f" (not: {', '.join(n['controls_failed'])})" if n["controls_failed"] else "")
        + f". Control E: held out, straight {n['curvature'][0]:.3f}, curved "
        f"{n['curvature'][1]:.3f}, to fifth powers {n['curvature'][2]:.3f} (left "
        f"{1 - n['curvature'][2] / n['ceiling']:.0%}): the curved row is a check, not "
        f"a bound. Control F: the components of the {n['f_genes']} genes measured in "
        f"every structure (most often {n['f_k']}, picked inside each fold) predict "
        f"{n['f_share']:.0%} of the reproducible map; on its own calibration the nano "
        f"map leaves {n['f_nano'][0]:.0%} to {n['f_nano'][1]:.0%}, a map made of those "
        f"genes {n['f_floor'][0]:.0%} to {n['f_floor'][1]:.0%}; its held-out R2 is "
        f"within {n['f_plateau_r2']} of its highest from {n['f_plateau']} components.",
        "",
        f"By division: {n['between']:.0%} of the leftover's variance lies between "
        f"divisions; above prediction as a whole {above or 'none'}, below "
        f"{below or 'none'} (95% over resampled adults). The leftover follows the "
        f"autofluorescence map at {n['rho_auto']:+.2f}.",
        "",
        f"Genes against the leftover: {n['n_pass']} of {n['n_genes']} past its "
        f"spatial null at BH q < {n['q']}, {n['n_pass_negative']} of them with a "
        f"negative rho ({n['n_p05']} below p 0.05 before correction); the highest rho "
        f"is {n['top_gene']}'s ({n['top_rho']:+.2f}, spatial p {n['top_p']:.4f}, q "
        f"{n['top_q']:.2f}); Cacng8 {n['cacng8']['rho']:+.2f} (p "
        f"{n['cacng8']['p_spatial']:.4f}). Over every gene, the rho with the leftover "
        f"and the rho with autofluorescence agree at {n['genes_against_auto']:+.2f}; "
        f"of the genes past BH, {n['pass_auto_positive']} follow autofluorescence "
        f"positively and {n['pass_auto_bh']} pass its own null.",
        "",
        "Not shown by any of this: what the leftover is. It is the part of the map "
        "that Gria1 expression and synapse density do not predict.",
    ]


# ===== Gathering the numbers =====


def gene_space_row(summary: pd.DataFrame, calibration: pd.DataFrame) -> dict:
    """Control F: its genes, components, nested share, replication and own floor."""
    row = summary.iloc[0]
    nano = calibration.loc[calibration["map"] == beyond_calibration.NANO, "left"]
    known = calibration.loc[calibration["map"] == "gene space", "left"]
    return dict(
        f_genes=int(row["genes"]),
        f_k=int(row["k_most"]),
        f_k_range=(int(row["k_min"]), int(row["k_max"])),
        f_share=float(row["share"]),
        f_rep=float(row["replication"]),
        f_cal_n=int(calibration["n_structures"].iloc[0]),
        f_cal_genes=int(calibration["n_genes"].iloc[0]),
        f_nano=(float(nano.min()), float(nano.max())),
        f_floor=(float(known.min()), float(known.max())),
        f_curve_peak=int(row["curve_peak"]),
    )


def plateau(curve: pd.DataFrame, within: float) -> int:
    """The fewest components whose held-out R2 is within `within` of the curve's highest.

    `curve` is gene_space.csv: control F's held-out R2 per number of components.
    """
    reached = curve[curve["cv_r2"] >= curve["cv_r2"].max() - within]
    return int(reached["n_components"].min())


def load_tables() -> dict:
    """The tables of steps 23 to 26 that the numbers and the figures read."""
    sets = pd.read_csv(
        beyond_density.LEFTOVER_SETS, keep_default_na=False, na_values=[""]
    )
    sets["tested"] = sets["tested"].astype(str) == "True"
    residuals = pd.read_csv(
        beyond_density.RESIDUALS, keep_default_na=False, na_values=[""]
    )
    partition, parts, weights = beyond_density.load_partition()
    return dict(
        partition=partition,
        parts=parts,
        weights=weights,
        replication=pd.read_csv(beyond_density.REPLICATION),
        residuals=residuals,
        divisions=pd.read_csv(beyond_density.RESIDUAL_DIVISIONS),
        genes=beyond_density.load_leftover_genes(),
        sets=sets,
        ranking=gene_ranking.load_ranking(),
        controls=pd.read_csv(CONTROLS),
        f_curve=pd.read_csv(GENE_SPACE),
        f_calibration=pd.read_csv(GENE_SPACE_CALIBRATION, keep_default_na=False),
        f_summary=pd.read_csv(GENE_SPACE_SUMMARY),
        folds=pd.read_csv(FOLDS),
        checks=beyond_checks.load_check_rows(),
        calibration=beyond_calibration.load_calibration(),
        paired=beyond_calibration.load_jackknife(),
        paired_blocks=beyond_calibration.load_jackknife(blocks=True),
        regression=beyond_regression.load_regression(),
    )


def declared_rows(used: pd.DataFrame) -> pd.DataFrame:
    """The declared structures of structures_used.csv, in the fit or not."""
    return used[~used["reason"].str.startswith("not in the declared set")]


def lost_structures(inputs: beyond_density.Inputs) -> tuple[dict, dict]:
    """What the rule of the structures loses, by the first gene a structure lacks.

    Returns {gene: {division: count}} and {gene: acronyms}, the genes in the order
    the model lists them, a structure counted under the first it lacks.
    """
    declared = declared_rows(inputs.structures_used)
    lost = declared[~declared["used"]]
    counts, acronyms = {}, {}
    for gene in beyond_density.model_genes(inputs.terms):
        first = lost["missing_genes"].str.split().str[0] == gene
        mine = lost[first]
        if mine.empty:
            continue
        counts[gene] = mine["division"].value_counts().to_dict()
        acronyms[gene] = sorted(mine["acronym"])
    return counts, acronyms


def jackknife_intervals(jack: pd.DataFrame, n: dict) -> dict:
    """The 95% intervals over structures of the shares, the parts and the weights."""

    def around(column: str, point: float) -> tuple[float, float]:
        """The interval of `point` from its column of the jackknife."""
        return beyond_calibration.jackknife_interval(jack, jack[column], point)

    return dict(
        ceiling_ci=around("ceiling", n["ceiling"]),
        shares_ci={k: around(k, n["shares"][k]) for k in INTERVAL_SHARES},
        parts_ci={k: around(k, v) for k, v in n["parts"].items()},
        weights_ci={k: around(f"weight_{k}", v) for k, v in n["weights"].items()},
    )


def failed_controls(controls: pd.DataFrame) -> list[str]:
    """The letters of the controls whose verdict is not a pass."""
    failed = controls.loc[controls["verdict"] != "pass", "control"]
    return [c[0] for c in failed]


def division_numbers(
    inputs: beyond_density.Inputs,
    residual: np.ndarray,
    covariates: dict[str, np.ndarray],
    divisions: pd.DataFrame,
) -> dict:
    """The leftover by division, and against the autofluorescence map.

    The share of its variance between divisions, the table of residual_by_division.csv
    and the divisions whose 95% over resampled adults lies above or below zero; its
    Spearman with the autofluorescence map.
    """
    division = [inputs.division.get(s, "") for s in inputs.structures]
    above = divisions[divisions["lo"] > 0]
    below = divisions[divisions["hi"] < 0].iloc[::-1]
    return dict(
        between=beyond_density.between_share(residual, division),
        divisions=divisions,
        divisions_above=list(above["division"]),
        divisions_below=list(below["division"]),
        rho_auto=float(spearmanr(residual, covariates["autofluo"]).statistic),
    )


def model_numbers(
    inputs: beyond_density.Inputs,
    covariates: dict[str, np.ndarray],
    tables: dict,
    rng: np.random.Generator,
) -> dict:
    """The main model's numbers, with their intervals over structures and blocks.

    The structures, the ceiling, the shares, the four parts, the weights and the
    replication; writes the two jackknife tables. One generator feeds the jackknife
    of the partition, the noise null of the replication and then the jackknife over
    spatial blocks.
    """
    splits = beyond_density.half_splits()
    y = beyond_density.full_map(inputs.nano)
    columns = beyond_density.model_columns(covariates, inputs.terms)
    map_agreement, explainable = beyond_density.ceiling(inputs.nano, splits)
    partition = tables["partition"]
    replication = tables["replication"]
    jack = beyond_calibration.jackknife_partition(inputs, splits, rng)
    jack.to_csv(JACKKNIFE, index=False)
    noise = noise_band(inputs.nano, columns, splits, rng)
    blocks = beyond_calibration.jackknife_partition(inputs, splits, rng, blocks=True)
    blocks.to_csv(JACKKNIFE_BLOCKS, index=False)
    used = inputs.structures_used
    lost, acronyms = lost_structures(inputs)
    weights = tables["weights"]["weight"].to_dict()
    n = dict(
        n_structures=len(inputs.structures),
        n_declared=len(declared_rows(used)),
        n_adults=len(ADULTS),
        n_splits=len(splits),
        abundance=", ".join(inputs.terms["abundance"]),
        density_genes=list(beyond_density.DENSITY_GENES),
        lost=lost,
        lost_acronyms=acronyms,
        psd95_measured=int(used["psd95_measured"].sum()),
        half=float(np.mean(map_agreement)),
        ceiling=explainable,
        partition_table=partition,
        shares=partition["share_of_ceiling"].to_dict(),
        parts=tables["parts"]["share"].to_dict(),
        weights=weights,
        rho_terms=tables["weights"]["rho_with_map"].to_dict(),
        rho_between=float(spearmanr(*(covariates[t] for t in weights)).statistic),
        left=float(partition.loc["model", "left"]),
        cv_r2=float(partition.loc["model", "cv_r2"]),
        r2=float(partition.loc["model", "in_sample_r2"]),
        n_terms=int(partition.loc["model", "terms"]),
        rep_map=float(np.mean(map_agreement)),
        rep_left=float(replication["leftover_agreement"].mean()),
        rep_left_min=float(replication["leftover_agreement"].min()),
        implied=beyond_density.implied_replication(
            float(np.mean(map_agreement)),
            beyond_density.r_squared(y, columns),
            len(columns) + 1,
            len(y),
        ),
        noise=noise,
        folds=tables["folds"],
        checks=tables["checks"],
        controls_passed=int((tables["controls"]["verdict"] == "pass").sum()),
        controls_failed=failed_controls(tables["controls"]),
        n_controls=len(tables["controls"]),
        curvature=curvature_scores(y, covariates, inputs.terms),
        f_plateau=plateau(tables["f_curve"], BEYOND_CONTROLS["plateau_r2"]),
        f_plateau_r2=BEYOND_CONTROLS["plateau_r2"],
        n_blocks=len(blocks),
        **gene_space_row(tables["f_summary"], tables["f_calibration"]),
    )
    residual = beyond_density.residual(y, columns)
    n |= division_numbers(inputs, residual, covariates, tables["divisions"])
    n |= jackknife_intervals(jack, n)
    n |= {f"{k}_blocks": v for k, v in jackknife_intervals(blocks, n).items()}
    n["left_ci"] = n["parts_ci"]["left"]
    n["left_ci_blocks"] = n["parts_ci_blocks"]["left"]
    return n


def calibration_numbers(inputs: beyond_density.Inputs, tables: dict) -> dict:
    """The calibration's numbers: the floor and nano beside it, with the interval.

    Nano is read on the same structures; once lists the model's genes measured by one
    Allen experiment, the same in both halves (none, by the rule of the density
    genes, and Gria1 has two).
    """
    calibration = tables["calibration"]
    paired, paired_blocks = tables["paired"], tables["paired_blocks"]
    random_folds = calibration[calibration["folds"] == beyond_calibration.RANDOM]
    nano_rows = random_folds[random_folds["map"] == beyond_calibration.NANO]
    floor = beyond_calibration.left_summary(calibration)
    nano_cal_left = beyond_calibration.nano_left(calibration)
    minus_floor = nano_cal_left - floor["left_median"]
    _, split = beyond_calibration.half_profiles()
    blocks = beyond_calibration.BLOCKS

    def around(table: pd.DataFrame, column: str, point: float) -> tuple[float, float]:
        """The interval of `point` from its column of a paired jackknife."""
        return beyond_calibration.jackknife_interval(table, table[column], point)

    return dict(
        cal_n=int(calibration["n_structures"].iloc[0]),
        n_draws=int((random_folds["map"] == beyond_calibration.FLOOR_MAP).sum()),
        floor=floor,
        own_half=beyond_calibration.own_half_left(calibration),
        nano_cal=nano_rows.set_index("predictors_from")["left"].to_dict(),
        nano_cal_left=nano_cal_left,
        nano_cal_ci=around(paired, "nano", nano_cal_left),
        nano_cal_ci_blocks=around(paired_blocks, "nano", nano_cal_left),
        minus_floor=minus_floor,
        minus_floor_ci=around(paired, "nano_minus_floor", minus_floor),
        minus_floor_ci_blocks=around(paired_blocks, "nano_minus_floor", minus_floor),
        blocks=dict(
            nano=beyond_calibration.nano_left(calibration, blocks),
            floor=beyond_calibration.left_summary(
                calibration, beyond_calibration.FLOOR_MAP, blocks
            )["left_median"],
        ),
        once=beyond_calibration.measured_once(split, inputs.terms),
    )


def leftover_numbers(
    genes: pd.DataFrame, sets: pd.DataFrame, ranking: pd.DataFrame
) -> dict:
    """The numbers of every gene and gene set against the leftover.

    Beside them, how each gene's rho with the leftover stands to its rho with the
    autofluorescence map (gene_ranking.csv): over every gene, and for the genes past
    BH on the leftover, how many follow autofluorescence and pass its own null.
    """
    q = ISH_ANALYSIS["q"]
    past = genes[genes["q_all"] < q].sort_values("rho")
    auto = gene_ranking.map_rows(ranking, "auto")
    both = genes.set_index("symbol")[["rho"]].join(auto[["rho"]], rsuffix="_auto")
    past_auto = auto.reindex(past["symbol"])
    return dict(
        genes_against_auto=float(
            spearmanr(both["rho"], both["rho_auto"], nan_policy="omit").statistic
        ),
        pass_auto_positive=int((past_auto["rho"] > 0).sum()),
        pass_auto_bh=int((past_auto["q_all"] < q).sum()),
        n_genes=len(genes),
        n_pass=len(past),
        n_pass_negative=int((past["rho"] < 0).sum()),
        pass_positive=list(past.loc[past["rho"] > 0, "symbol"][::-1]),
        pass_lowest=list(past.loc[past["rho"] < 0, "symbol"][:10]),
        q=q,
        top_gene=genes.iloc[0]["symbol"],
        top_rho=float(genes.iloc[0]["rho"]),
        top_p=float(genes.iloc[0]["p_spatial"]),
        top_q=float(genes.iloc[0]["q_all"]),
        n_p05=int((genes["p_spatial"] < ALPHA).sum()),
        cacng8=genes.set_index("symbol").loc[LEFTOVER_GENE],
        sets=sets,
    )


# ===== The figures =====


def leftover_span(residual: pd.Series) -> float:
    """Figure 04's colour scale for the leftover, in ranks each way.

    The beyond_figures.leftover_span_pct percentile of |leftover|, rounded up to ten
    ranks, so about the share of the structures beyond it saturate and the rest keep
    their shades.
    """
    value = np.percentile(
        np.abs(residual.to_numpy(float)), BEYOND_FIGURES["leftover_span_pct"]
    )
    return float(10 * np.ceil(value / 10))


def plane_images(regression: pd.DataFrame) -> list[dict]:
    """The planes of figure 04: labels, and map, prediction and leftover painted."""
    n = len(regression)
    values = {
        "map": dict(
            zip(regression["structure"], (regression["observed_rank"] - 1) / (n - 1))
        ),
        "prediction": dict(
            zip(regression["structure"], (regression["predicted_rank"] - 1) / (n - 1))
        ),
        "leftover": dict(zip(regression["structure"], regression["residual"])),
    }
    out = []
    for plane in beyond_regression.coronal_planes():
        lab = plane["lab"]
        out.append(
            dict(
                ccf_plane=plane["ccf_plane"],
                lab=lab,
                **{k: paint(lab, v, plane["names"]) for k, v in values.items()},
            )
        )
    return out


def control_data(
    inputs: beyond_density.Inputs, covariates: dict[str, np.ndarray], tables: dict
) -> dict:
    """What figure 03s2 draws of the seven controls of adult.beyond_controls.

    Recomputed from the main model as that module computes it: the leftover with
    each structure's anterior-posterior centroid and mean volume, every pair of
    adults' leftovers, the naive and the RWS group's, and the held-out R2 straight,
    curved and to fifth powers. Read from its tables: control F's curve and own
    calibration, control G's readings and every verdict.
    """
    s = inputs.structures
    y = beyond_density.full_map(inputs.nano)
    columns = beyond_density.model_columns(covariates, inputs.terms)
    per_adult = beyond_density.per_adult_leftovers(inputs.nano, columns)
    naive, rws = group_leftovers(inputs.nano, columns)
    return dict(
        residual=beyond_density.residual(y, columns),
        ap=structures.centroid_xyz(s)[:, 0],
        sizes=mean_sizes(s),
        pairs=adult_pair_agreements(per_adult),
        groups=dict(naive=naive, rws=rws),
        curvature=curvature_scores(y, covariates, inputs.terms),
        curve=pd.read_csv(GENE_SPACE),
        f_calibration=tables["f_calibration"],
        readings=pd.read_csv(READINGS),
        controls=tables["controls"],
        thresholds=BEYOND_CONTROLS,
    )


def draw_beyond(
    inputs: beyond_density.Inputs,
    covariates: dict[str, np.ndarray],
    numbers: dict,
    tables: dict,
) -> None:
    """Draw figure 03 and its detailed versions 03s1 and 03s2.

    The map against its predictors, the partition, the floor, the density term
    against the measured synapse density, the check rows, the replication, and the
    seven controls (03s3, the choice of the density genes, is step 22's).
    """
    s = inputs.structures
    y = beyond_density.full_map(inputs.nano)
    columns = beyond_density.model_columns(covariates, inputs.terms)
    groups = group_of(structures.load_structure_set())
    shared = dict(
        y=y,
        held_out=beyond_density.cv_predict(y, columns),
        residual=beyond_density.residual(y, columns),
        groups=[groups[x] for x in s],
        acronyms=[inputs.acronym.get(x, "") for x in s],
        numbers=numbers,
        calibration=tables["calibration"],
    )
    fig = adult_plotting.plot_beyond(
        **shared,
        markers=density_markers.load_validation(),
        structure_groups=groups,
        save=figure_path("beyond"),
    )
    plt.close(fig)
    terms = {
        "Gria1": covariates[beyond_density.ABUNDANCE[0]],
        "density": covariates["density"],
    }
    fig = adult_plotting.plot_beyond_budget(
        terms=terms,
        replication=tables["replication"],
        **shared,
        save=figure_path("beyond_budget"),
    )
    plt.close(fig)
    fig = adult_plotting.plot_beyond_controls(
        control_data(inputs, covariates, tables),
        numbers,
        save=figure_path("beyond_controls"),
    )
    plt.close(fig)


def draw_figures(
    inputs: beyond_density.Inputs,
    covariates: dict[str, np.ndarray],
    numbers: dict,
    tables: dict,
) -> None:
    """Draw figures 03, 03s1, 03s2, 04 and 11s1 from the main model and tables."""
    draw_beyond(inputs, covariates, numbers, tables)

    # figure 04: where the leftover lives; figure 11s1: every gene against it
    fig = adult_plotting.plot_beyond_where(
        planes=plane_images(tables["regression"]),
        residuals=tables["residuals"],
        divisions=tables["divisions"],
        numbers=numbers,
        t_max=BEYOND_FIGURES["t_max_structures"],
        span=leftover_span(tables["residuals"]["residual"]),
        save=figure_path("beyond_where"),
    )
    plt.close(fig)
    with np.load(beyond_density.LEFTOVER_NULL) as null:
        n_surrogates = int(null["surrogates"].shape[0])
    auto = gene_ranking.map_rows(tables["ranking"], "auto")["rho"]
    fig = adult_plotting.plot_leftover_genes(
        genes=tables["genes"],
        sets=tables["sets"],
        auto=auto,
        numbers=numbers,
        n_surrogates=n_surrogates,
        t_max=ISH_FIGURES["t_max"],
        save=figure_path("leftover_genes"),
    )
    plt.close(fig)
    keys = (
        "beyond",
        "beyond_budget",
        "beyond_controls",
        "beyond_where",
        "leftover_genes",
    )
    drawn = [figure_file(k) for k in keys]
    print(f"figures: {', '.join(drawn)} in {structures.FIGURES}")


def main() -> None:
    """Compute the intervals, write the numbers and draw the figures of part 1.

    One seeded generator feeds every resampling, in the order of model_numbers.
    """
    rng = np.random.default_rng(0)
    inputs = beyond_density.load_inputs()
    covariates, _, _ = beyond_density.covariates_for(inputs)
    tables = load_tables()
    print(f"{len(inputs.structures)} structures, {len(ADULTS)} adults")

    # the numbers, with their intervals over structures, and the caption
    numbers = (
        model_numbers(inputs, covariates, tables, rng)
        | calibration_numbers(inputs, tables)
        | leftover_numbers(tables["genes"], tables["sets"], tables["ranking"])
    )
    numbers_table(numbers).to_csv(NUMBERS, index=False)
    lines = caption_lines(numbers)
    CAPTION.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))

    # the figures of part 1
    draw_figures(inputs, covariates, numbers, tables)
