"""The numbers and the guided figures of part 1: 03, 03s1, 03s2, 04, 11s1 and 14.

adult.beyond_density, adult.beyond_controls, adult.beyond_calibration and
adult.beyond_regression write the tables. This module adds the intervals over
structures (the jackknives of adult.beyond_calibration), gathers the numbers the
text quotes, and has adult.plotting draw the figures of part 1:

    03   what Gria1 and synapse density predict of the map: the map against the
         main model's prediction, the budget with the leftover's range, the
         leftover beside the calibration floor and the benchmark, its replication
    03s1 the same in detail: the map against each predictor, the budget beside its
         four-subunit check row and control F, the calibration, the leftover under
         other folds, structures and models
    03s2 the seven controls of adult.beyond_controls, one panel each
    04   where the leftover lives: map, prediction and leftover on three coronal
         planes, and the structures with the largest leftovers
    11s1 every gene of the gene table and every gene set against the leftover, with
         the leftover's spatial null
    14   the measured synapse density: how much of the fit it covers, against the
         mRNA it would replace, and what each density measure leaves (step 21
         draws its detailed version)

The leftover and the floor are compared on the same structures, the calibration's,
and their difference is resampled with them (beyond_calibration.paired_jackknife);
so are the check rows that change only the density measure, on the structures PSD95
covers (beyond_calibration.jackknife_density_rows). The noise null of the
replication shuffles which structure is which in one half's leftover: what two
unrelated leftovers would give. It is quoted for scale only, since the leftover of a
reproducible map was always going to replicate (beyond_density.implied_replication).

Writes, under adult_v2/ish_analysis/ in the data root:

    beyond/jackknife.csv                per subsample: the ceiling and each model's
                                        share of it
    beyond/jackknife_density_rows.csv   per subsample of the structures PSD95
                                        covers: each density row's share left and
                                        its density alone, and their differences
    beyond/numbers_for_the_caption.txt  the figures' numbers as sentences
    tables/numbers_beyond.csv           the numbers of analysis 4, for the text
    figures/03_beyond.png               (and .eps, as every figure)
    figures/03s1_beyond_budget.png
    figures/03s2_beyond_controls.png
    figures/04_beyond_where.png
    figures/11s1_leftover_genes.png
    figures/14_synaptome.png

Run by run_beyond_figures.py.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap import structures
from sepmap.adult import beyond_calibration, beyond_density, beyond_regression, synaptome
from sepmap.adult import plotting as adult_plotting
from sepmap.adult.beyond_controls import (
    CONTROLS,
    GENE_SPACE,
    GENE_SPACE_CALIBRATION,
    GENE_SPACE_SUMMARY,
    READINGS,
    VARIANTS,
    adult_pair_agreements,
    curvature_scores,
    group_leftovers,
    mean_sizes,
)
from sepmap.adult.profiles import ADULTS
from sepmap.config import SETTINGS
from sepmap.ish.figure_index import figure_file, figure_path
from sepmap.ish.gene_sets import LEFTOVER_GENE
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.ish.plotting import group_of
from sepmap.ish.spatial_null import ALPHA
from sepmap.plotting import paint

# the permutations of the noise null and the grey of the bars; the BH level; the grey
# of the genes' bars; the thresholds of the controls' verdicts, which 03s2 draws
BEYOND_FIGURES = SETTINGS["beyond_figures"]
ISH_ANALYSIS = SETTINGS["ish_analysis"]
ISH_FIGURES = SETTINGS["ish_figures"]
BEYOND_CONTROLS = SETTINGS["beyond_controls"]

JACKKNIFE = beyond_density.OUT / "jackknife.csv"
DENSITY_JACKKNIFE = beyond_density.OUT / "jackknife_density_rows.csv"
CAPTION = beyond_density.OUT / "numbers_for_the_caption.txt"
NUMBERS = numbers_path("beyond")

# the percentiles of a 95% interval
INTERVAL = (2.5, 97.5)


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


def model_number_rows(n: dict) -> list[tuple]:
    """The main model's numbers: structures, ceiling, budget and replication."""
    rows = [
        ("structures", n["n_structures"], "structures of the main model's fit"),
        ("adults", n["n_adults"], "adults, naive and RWS pooled"),
        ("main_abundance", n["abundance"], "the main model's abundance term"),
        ("main_density", n["density"], "the main model's synapse density terms"),
        (
            "psd95_measured",
            n["psd95_measured"],
            "structures of the fit with a measured PSD95 density",
        ),
        ("psd95_needed", n["psd95_needed"], "structures the rule asks PSD95 to cover"),
        (
            "half_agreement",
            round(n["half"], 4),
            f"half-cohort maps agree ({n['n_splits']} splits)",
        ),
        ("ceiling", round(n["ceiling"], 4), "Spearman-Brown: reproducible share"),
        ("ceiling_lo", round(n["ceiling_ci"][0], 6), "2.5%, jackknife over structures"),
        ("ceiling_hi", round(n["ceiling_ci"][1], 6), "97.5%, jackknife over structures"),
        ("psd_genes", n["psd_genes"], "postsynaptic-density genes behind psd_pc1"),
    ]
    for key, r in n["partition_table"].iterrows():
        rows.append(
            (
                f"share_{key}",
                round(r["share_of_ceiling"], 4),
                f"{r['model']}: CV R2 / ceiling",
            )
        )
    rows += [
        ("budget_abundance", round(n["steps"][0], 4), f"budget: {n['abundance']}"),
        ("budget_density", round(n["steps"][1] - n["steps"][0], 4), "budget: + density"),
        ("budget_auto", round(n["steps"][2] - n["steps"][1], 4), "budget: + autofluo"),
        ("left", round(1 - n["steps"][2], 4), "share of the reproducible map left"),
        ("left_lo", round(n["left_ci"][0], 6), "2.5%, jackknife over structures"),
        ("left_hi", round(n["left_ci"][1], 6), "97.5%, jackknife over structures"),
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
    """The calibration's numbers: the floor, the benchmark and nano beside them."""
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
        ("gria1_map_left", round(n["gria1"]["left_median"], 4), "Gria1 map: left"),
        ("gria1_map_left_lo", round(n["gria1"]["left_lo"], 4), "lowest draw"),
        ("gria1_map_left_hi", round(n["gria1"]["left_hi"], 4), "highest draw"),
        (
            "gria1_map_replication",
            round(n["gria1"]["replication_median"], 4),
            "Gria1 map: leftover replicates",
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
    rows += [
        (
            "calibration_markers_once",
            " ".join(n["markers_once"]),
            "the main model's genes and measures that are the same in both halves",
        ),
        (
            "calibration_psd_once",
            len(n["psd_once"]),
            f"of the {n['psd_genes']} genes behind psd_pc1, those measured once",
        ),
        ("calibration_psd_once_genes", " ".join(n["psd_once"]), "those genes"),
    ]
    for h, count in n["psd_halves"].items():
        rows.append(
            (
                f"calibration_psd_genes_{h}",
                count,
                f"genes behind psd_pc1 on the calibration's structures, from {h}",
            )
        )
    rows += [
        ("nano_calibration_left", round(n["nano_cal_left"], 4), "nano, mean of A and B"),
        ("nano_calibration_lo", round(n["nano_cal_ci"][0], 6), "2.5%, jackknife"),
        ("nano_calibration_hi", round(n["nano_cal_ci"][1], 6), "97.5%, jackknife"),
        (
            "nano_floor_ratio",
            round(n["nano_cal_left"] / n["floor"]["left_median"], 2),
            "nano's share left over the floor's, same structures",
        ),
        ("nano_minus_floor", round(n["minus_floor"], 4), "nano minus the floor"),
        ("nano_minus_floor_lo", round(n["minus_floor_ci"][0], 6), "2.5%, jackknife"),
        ("nano_minus_floor_hi", round(n["minus_floor_ci"][1], 6), "97.5%, jackknife"),
        ("nano_minus_gria1", round(n["minus_gria1"], 4), "nano minus the Gria1 map"),
        ("nano_minus_gria1_lo", round(n["minus_gria1_ci"][0], 6), "2.5%, jackknife"),
        ("nano_minus_gria1_hi", round(n["minus_gria1_ci"][1], 6), "97.5%, jackknife"),
        ("gria1_map_left_A", round(n["gria1_halves"]["A"], 4), "Gria1 map from half A"),
        ("gria1_map_left_B", round(n["gria1_halves"]["B"], 4), "Gria1 map from half B"),
        ("blocks_nano_left", round(n["blocks"]["nano"], 4), "spatial blocks: nano"),
        ("blocks_floor", round(n["blocks"]["floor"], 4), "spatial blocks: the floor"),
        ("blocks_gria1_map", round(n["blocks"]["gria1"], 4), "spatial blocks: Gria1 map"),
    ]
    return rows


def variant_number_rows(n: dict) -> list[tuple]:
    """The check rows and other folds, and the density rows' paired differences."""
    rows = []
    for r in n["variants"].itertuples():
        if r.kind == "main":
            continue
        rows.append((f"variant_left_{r.key}", round(r.left, 4), r.variant))
        rows.append((f"variant_structures_{r.key}", int(r.n_structures), r.variant))
        rows.append(
            (
                f"variant_density_alone_{r.key}",
                round(r.density_alone, 4),
                f"{r.variant}: density alone",
            )
        )
        if np.isfinite(r.lo):
            rows.append((f"variant_left_lo_{r.key}", round(r.lo, 4), "2.5%"))
            rows.append((f"variant_left_hi_{r.key}", round(r.hi, 4), "97.5%"))
    for name, (point, (lo, hi)) in n["density_contrasts"].items():
        rows += [
            (f"check_{name}", round(point, 4), "difference of two check rows, paired"),
            (f"check_{name}_lo", round(lo, 6), "2.5%, jackknife over their structures"),
            (f"check_{name}_hi", round(hi, 6), "97.5%, jackknife over their structures"),
        ]
    return rows


def control_and_gene_number_rows(n: dict) -> list[tuple]:
    """The controls' numbers, and every gene and gene set against the leftover."""
    rows = [
        ("controls_passed", n["controls_passed"], "of the seven controls"),
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
        ("leftover_genes", n["n_genes"], "genes against the leftover"),
        ("leftover_genes_pass", n["n_pass"], f"past its null, BH q < {n['q']}"),
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
    rows = model_number_rows(n) + calibration_number_rows(n) + variant_number_rows(n)
    rows += control_and_gene_number_rows(n)
    return numbers_frame(rows)


def contrast_text(n: dict, name: str) -> str:
    """A difference between check rows in points, with its interval: '+7 (-11 to +25)'."""
    point, (lo, hi) = n["density_contrasts"][name]
    return f"{100 * point:+.0f} points (95% {100 * lo:+.0f} to {100 * hi:+.0f})"


def caption_lines(n: dict) -> list[str]:
    """numbers_for_the_caption.txt: the figures' numbers as sentences."""
    floor = n["floor"]
    gria1 = n["gria1"]
    if n["psd95_measured"] >= n["psd95_needed"]:
        rule = "at least what the rule asks, so PSD95 puncta are the density"
    else:
        rule = "below what the rule asks, so the mRNA panel stays the density"
    checks = n["variants"][n["variants"]["kind"] == "model"]
    listed = "; ".join(
        f"{r.variant}, {r.n_structures} structures: {r.left:.0%} left"
        for r in checks.itertuples()
    )
    return [
        f"Analysis 4, on {n['n_structures']} grey-matter structures of the declared "
        f"set and {n['n_adults']} adults. The main model, fixed on 8 October 2026 "
        f"before it was run: {n['abundance']} + synapse density ({n['density']}) + "
        "autofluorescence, each as x, x^2 and x^3. PSD95 density is measured in "
        f"{n['psd95_measured']} of the {n['n_fit']} structures of the fit "
        f"({n['psd95_needed']} needed), {rule}.",
        "",
        f"Ceiling. Two halves of the cohort agree at {n['half']:.3f} "
        f"({n['n_splits']} splits), so "
        f"{n['ceiling']:.1%} of the map is reproducible (Spearman-Brown; 95% over "
        f"structures, jackknife, {n['ceiling_ci'][0]:.1%} to {n['ceiling_ci'][1]:.1%}).",
        "",
        f"Budget. On structures the fit has not seen, {n['abundance']} predicts "
        f"{n['steps'][0]:.0%} of the reproducible map, synapse density adds "
        f"{n['steps'][1] - n['steps'][0]:.0%} and autofluorescence "
        f"{n['steps'][2] - n['steps'][1]:+.0%}; {1 - n['steps'][2]:.0%} is not "
        f"predicted by Gria1 expression or synapse density (95% over structures, "
        f"jackknife, {n['left_ci'][0]:.0%} to {n['left_ci'][1]:.0%}). Density alone "
        f"predicts {n['partition']['density']:.0%}; Gria1 to Gria4 as four terms "
        f"alone {n['partition']['subunits']:.0%}.",
        "",
        f"Check rows, each one change to the main model on its own structures: {listed}. "
        "On the structures PSD95 covers the rows meet the same subsamples, so their "
        "differences have intervals of their own: the share left with PSD95 minus "
        f"that with the mRNA panel, {contrast_text(n, 'psd95_minus_panel_left')}; with "
        f"both, minus the panel, {contrast_text(n, 'both_minus_panel_left')}; the "
        "share density alone predicts, PSD95 minus the panel, "
        f"{contrast_text(n, 'psd95_minus_panel_alone')}.",
        "",
        f"Calibration ({n['cal_n']} structures where both halves of the Allen "
        f"experiments measure every subunit and marker). A map made only of Gria1 "
        f"and synapse density, predicted from the other half of the experiments, "
        f"leaves {floor['left_median']:.0%} ({floor['left_lo']:.0%} to "
        f"{floor['left_hi']:.0%} over {n['n_draws']} draws), its leftover replicating "
        f"at {floor['replication_median']:.2f}. The floor errs low: "
        f"{', '.join(n['markers_once'])} and {len(n['psd_once'])} of the "
        f"{n['psd_genes']} genes behind psd_pc1 are measured once and so the same in "
        "both halves, and it holds only the mismatch of one Allen map with another, "
        "not that of Allen's P56 mice with these brains. The nano map, predicted the "
        f"same way on the same structures, leaves {n['nano_cal_left']:.0%} "
        f"({n['nano_cal_ci'][0]:.0%} to {n['nano_cal_ci'][1]:.0%}), "
        f"{n['minus_floor']:+.0%} above the floor at its point value (95% "
        f"{n['minus_floor_ci'][0]:+.0%} to {n['minus_floor_ci'][1]:+.0%}, jackknife "
        "over the structures with both recomputed). A map of one Allen experiment of "
        f"Gria1 leaves {gria1['left_median']:.0%} ({n['gria1_halves']['A']:.0%} from "
        f"one half, {n['gria1_halves']['B']:.0%} from the other), replicating at "
        f"{gria1['replication_median']:.2f}; nano minus that, "
        f"{n['minus_gria1']:+.0%} ({n['minus_gria1_ci'][0]:+.0%} to "
        f"{n['minus_gria1_ci'][1]:+.0%}).",
        "",
        f"Replication. The leftover of one half-cohort agrees with the other's at "
        f"{n['rep_left']:.3f} (lowest split {n['rep_left_min']:.3f}); the map's "
        f"reliability and the fit alone imply {n['implied']:.3f}; two unrelated "
        f"leftovers fall between {n['noise'][0]:+.2f} and {n['noise'][1]:+.2f}.",
        "",
        f"Controls. {n['controls_passed']} of {n['n_controls']} pass. Control F: the "
        "components of "
        f"the {n['f_genes']} genes measured in every structure (most often "
        f"{n['f_k']}, picked inside each fold) predict {n['f_share']:.0%} of the "
        "reproducible map; on its own calibration the nano map leaves "
        f"{n['f_nano'][0]:.0%} to {n['f_nano'][1]:.0%}, a map made of those genes "
        f"{n['f_floor'][0]:.0%} to {n['f_floor'][1]:.0%}.",
        "",
        f"Genes against the leftover: {n['n_pass']} of {n['n_genes']} past its "
        f"spatial null at BH q < {n['q']} ({n['n_p05']} below p 0.05 before "
        f"correction); the highest rho is {n['top_gene']}'s ({n['top_rho']:+.2f}, "
        f"spatial p {n['top_p']:.4f}, q {n['top_q']:.2f}); Cacng8 "
        f"{n['cacng8']['rho']:+.2f} (p {n['cacng8']['p_spatial']:.4f}).",
        "",
        "Not shown by any of this: what the leftover is. It is the part of the map "
        "that Gria1 expression and synapse density do not predict.",
    ]


# ===== Gathering the numbers =====


def density_label(terms: dict[str, tuple[str, ...]], n_psd: int) -> str:
    """The main model's density terms in words, for figures 03 and 03s1."""
    parts = []
    if "markers" in terms["density"]:
        parts.append(f"{len(beyond_density.MARKERS)} marker genes")
    if "psd_pc1" in terms["density"]:
        parts.append(f"the first component of {n_psd} postsynaptic-density genes")
    parts += [f"{m} puncta" for m in terms["density"] if m in synaptome.MEASURES]
    return " and ".join(parts)


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


def load_tables() -> dict:
    """The tables of steps 23 to 26 that the numbers and the figures read."""
    sets = pd.read_csv(
        beyond_density.LEFTOVER_SETS, keep_default_na=False, na_values=[""]
    )
    sets["tested"] = sets["tested"].astype(str) == "True"
    residuals = pd.read_csv(
        beyond_density.RESIDUALS, keep_default_na=False, na_values=[""]
    )
    return dict(
        partition=pd.read_csv(beyond_density.PARTITION).set_index("key"),
        replication=pd.read_csv(beyond_density.REPLICATION),
        residuals=residuals,
        genes=beyond_density.load_leftover_genes(),
        sets=sets,
        controls=pd.read_csv(CONTROLS),
        f_calibration=pd.read_csv(GENE_SPACE_CALIBRATION, keep_default_na=False),
        f_summary=pd.read_csv(GENE_SPACE_SUMMARY),
        variants=pd.read_csv(VARIANTS),
        calibration=beyond_calibration.load_calibration(),
        paired=beyond_calibration.load_jackknife(),
        regression=beyond_regression.load_regression(),
    )


def model_numbers(
    inputs: beyond_density.Inputs,
    covariates: dict[str, np.ndarray],
    psd: list[str],
    tables: dict,
    rng: np.random.Generator,
) -> dict:
    """The main model's numbers, with their intervals over structures.

    The ceiling, the budget, the replication and the check rows. Writes the two
    jackknife tables. One generator feeds the jackknife of the budget, the noise
    null of the replication and the jackknife of the check rows, in that order.
    """
    splits = beyond_density.half_splits()
    y = beyond_density.full_map(inputs.nano)
    columns = beyond_density.model_columns(covariates, inputs.terms)
    map_agreement, explainable = beyond_density.ceiling(inputs.nano, splits)
    partition = tables["partition"]
    replication = tables["replication"]
    jack = beyond_calibration.jackknife_budget(inputs, splits, rng)
    jack.to_csv(JACKKNIFE, index=False)
    noise = noise_band(inputs.nano, columns, splits, rng)
    density_jack = beyond_calibration.jackknife_density_rows(inputs, splits, rng)
    density_jack.to_csv(DENSITY_JACKKNIFE, index=False)
    contrasts = beyond_calibration.density_contrasts(tables["variants"], density_jack)
    steps = beyond_density.budget(y, covariates, inputs.terms, explainable)
    used = inputs.structures_used
    n_fit = int(used["used"].sum())
    return dict(
        n_structures=len(inputs.structures),
        n_adults=len(ADULTS),
        n_splits=len(splits),
        abundance=", ".join(inputs.terms["abundance"]),
        density=", ".join(inputs.terms["density"]),
        density_label=density_label(inputs.terms, len(psd)),
        psd95_measured=int(used["psd95_measured"].sum()),
        n_fit=n_fit,
        psd95_needed=synaptome.n_needed(n_fit),
        half=float(np.mean(map_agreement)),
        ceiling=explainable,
        ceiling_ci=beyond_calibration.jackknife_interval(
            jack, jack["ceiling"], explainable
        ),
        psd_genes=len(psd),
        partition=partition["share_of_ceiling"].to_dict(),
        partition_table=partition,
        steps=steps,
        left_ci=beyond_calibration.jackknife_interval(jack, jack["left"], 1 - steps[2]),
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
        variants=tables["variants"],
        density_contrasts=contrasts,
        controls_passed=int((tables["controls"]["verdict"] == "pass").sum()),
        n_controls=len(tables["controls"]),
        **gene_space_row(tables["f_summary"], tables["f_calibration"]),
    )


def calibration_numbers(
    inputs: beyond_density.Inputs, psd: list[str], tables: dict
) -> dict:
    """The calibration's numbers: the floor, the benchmark and nano beside them.

    Nano is read on the same structures, with the intervals of the differences, and
    beside it what makes the floor err low. psd_halves counts the genes behind
    psd_pc1 on the calibration's structures as each half of the experiments builds
    it, and as the merged profiles do.
    """
    calibration = tables["calibration"]
    paired = tables["paired"]
    random_folds = calibration[calibration["folds"] == beyond_calibration.RANDOM]
    nano_rows = random_folds[random_folds["map"] == beyond_calibration.NANO]
    nano_cal = nano_rows.set_index("predictors_from")
    floor = beyond_calibration.left_summary(calibration, beyond_calibration.FLOOR_MAP)
    gria1 = beyond_calibration.left_summary(calibration, beyond_calibration.GRIA1_MAP)
    nano_cal_left = float(nano_cal.loc[list(beyond_calibration.HALVES), "left"].mean())
    minus_floor = nano_cal_left - floor["left_median"]
    minus_gria1 = nano_cal_left - gria1["left_median"]
    gria1_rows = random_folds[random_folds["map"] == beyond_calibration.GRIA1_MAP]
    halves, split = beyond_calibration.half_profiles()
    cal = beyond_density.restrict(
        inputs, beyond_calibration.calibration_structures(inputs.structures, halves)
    )
    psd_halves = {
        h: len(
            beyond_density.build_covariates(
                halves[h], cal.role, cal.auto, cal.structures
            )[1]
        )
        for h in beyond_calibration.HALVES
    }
    psd_halves[beyond_calibration.MERGED] = len(beyond_density.covariates_for(cal)[1])
    blocks = beyond_calibration.BLOCKS
    return dict(
        cal_n=int(calibration["n_structures"].iloc[0]),
        n_draws=int((random_folds["map"] == beyond_calibration.FLOOR_MAP).sum()),
        floor=floor,
        gria1=gria1,
        nano_cal=nano_cal["left"].to_dict(),
        nano_cal_left=nano_cal_left,
        nano_cal_ci=beyond_calibration.jackknife_interval(
            paired, paired["nano"], nano_cal_left
        ),
        minus_floor=minus_floor,
        minus_floor_ci=beyond_calibration.jackknife_interval(
            paired, paired["nano_minus_floor"], minus_floor
        ),
        minus_gria1=minus_gria1,
        minus_gria1_ci=beyond_calibration.jackknife_interval(
            paired, paired["nano_minus_gria1"], minus_gria1
        ),
        gria1_halves=gria1_rows.groupby("truth_from")["left"].median().to_dict(),
        blocks=dict(
            nano=beyond_calibration.left_summary(
                calibration, beyond_calibration.NANO, blocks
            )["left_median"],
            floor=beyond_calibration.left_summary(
                calibration, beyond_calibration.FLOOR_MAP, blocks
            )["left_median"],
            gria1=beyond_calibration.left_summary(
                calibration, beyond_calibration.GRIA1_MAP, blocks
            )["left_median"],
        ),
        markers_once=beyond_calibration.measured_once(split, inputs.terms),
        psd_once=[g for g in psd if g not in split],
        psd_halves=psd_halves,
    )


def leftover_numbers(genes: pd.DataFrame, sets: pd.DataFrame) -> dict:
    """The numbers of every gene and gene set against the leftover."""
    q = ISH_ANALYSIS["q"]
    return dict(
        n_genes=len(genes),
        n_pass=int((genes["q_all"] < q).sum()),
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
    adults' leftovers, the naive and the RWS group's, and the held-out R2 bent to
    each power. Read from its tables: control F's curve and own calibration,
    control G's readings and every verdict.
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

    The map against its predictors, the budget, the floor and the benchmark, the
    replication, and the seven controls.
    """
    s = inputs.structures
    y = beyond_density.full_map(inputs.nano)
    columns = beyond_density.model_columns(covariates, inputs.terms)
    density_columns = beyond_density.model_columns(covariates, inputs.terms, ("density",))
    groups = group_of(structures.load_structure_set())
    held_out = dict(
        density=beyond_density.cv_predict(y, density_columns),
        model=beyond_density.cv_predict(y, columns),
    )
    shared = dict(
        y=y,
        held_out=held_out,
        residual=beyond_density.residual(y, columns),
        groups=[groups[x] for x in s],
        acronyms=[inputs.acronym.get(x, "") for x in s],
        numbers=numbers,
        calibration=tables["calibration"],
        replication=tables["replication"],
    )
    fig = adult_plotting.plot_beyond(**shared, save=figure_path("beyond"))
    plt.close(fig)
    fig = adult_plotting.plot_beyond_budget(
        gria1=covariates[beyond_density.ABUNDANCE[0]],
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
    """Draw figures 03, 03s1, 03s2, 04, 11s1 and 14 from the main model and tables."""
    draw_beyond(inputs, covariates, numbers, tables)

    # figure 04: where the leftover lives; figure 11s1: every gene against it
    fig = adult_plotting.plot_beyond_where(
        planes=plane_images(tables["regression"]),
        residuals=tables["residuals"],
        numbers=numbers,
        t_max=BEYOND_FIGURES["t_max_structures"],
        save=figure_path("beyond_where"),
    )
    plt.close(fig)
    with np.load(beyond_density.LEFTOVER_NULL) as null:
        n_surrogates = int(null["surrogates"].shape[0])
    fig = adult_plotting.plot_leftover_genes(
        genes=tables["genes"],
        sets=tables["sets"],
        numbers=numbers,
        n_surrogates=n_surrogates,
        t_max=ISH_FIGURES["t_max"],
        save=figure_path("leftover_genes"),
    )
    plt.close(fig)

    # figure 14: the measured synapse density against the mRNA it would replace, with
    # the check rows that put it in the model
    fig = adult_plotting.plot_synaptome(
        density=synaptome.load_density().reset_index(),
        coverage=pd.read_csv(synaptome.COVERAGE),
        markers=pd.Series(covariates["markers"], index=inputs.structures),
        numbers=numbers,
        min_coverage=beyond_density.BEYOND["min_psd95_coverage"],
        save=figure_path("synaptome"),
    )
    plt.close(fig)
    keys = (
        "beyond",
        "beyond_budget",
        "beyond_controls",
        "beyond_where",
        "leftover_genes",
        "synaptome",
    )
    drawn = [figure_file(k) for k in keys]
    print(f"figures: {', '.join(drawn)} in {structures.FIGURES}")


def main() -> None:
    """Compute the intervals, write the numbers and draw the figures of part 1.

    One seeded generator feeds every resampling, in the order of model_numbers.
    """
    rng = np.random.default_rng(0)
    inputs = beyond_density.load_inputs()
    covariates, psd, _ = beyond_density.covariates_for(inputs)
    tables = load_tables()
    print(f"{len(inputs.structures)} structures, {len(ADULTS)} adults")

    # the numbers, with their intervals over structures, and the caption
    numbers = (
        model_numbers(inputs, covariates, psd, tables, rng)
        | calibration_numbers(inputs, psd, tables)
        | leftover_numbers(tables["genes"], tables["sets"])
    )
    numbers_table(numbers).to_csv(NUMBERS, index=False)
    lines = caption_lines(numbers)
    CAPTION.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))

    # the figures of part 1
    draw_figures(inputs, covariates, numbers, tables)
