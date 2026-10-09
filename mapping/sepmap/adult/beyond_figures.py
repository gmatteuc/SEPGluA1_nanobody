"""The guided figures of analysis 4, with their intervals: figures 03, 04, 11s1 and 14.

adult.beyond_density, adult.beyond_controls, adult.beyond_calibration and
adult.beyond_regression write the tables. This module adds the intervals, gathers
the numbers the text quotes, and has ish.plotting draw the figures of part 1:

    03   what Gria1 and synapse density predict of the map: the map against the
         main model's prediction, the budget with the leftover's range, the
         leftover beside the calibration floor and the benchmark, its replication
    03s  the same in detail: the map against each predictor, the budget beside its
         four-subunit check row and control F, the calibration, the leftover under
         other folds, structures and models
    04   where the leftover lives: map, prediction and leftover on three coronal
         planes, and the structures with the largest leftovers
    11s1 every gene of the gene table and every gene set against the leftover, with
         the leftover's spatial null
    14   the measured synapse density: how much of the fit it covers, against the
         mRNA it would replace, and what each density measure leaves (step 21
         draws its detailed version)

Intervals come from resampling structures, not animals: the claim is about where in
the brain the map departs from prediction, so what would differ in a repeat of the
analysis is which structures were measurable. A delete-d jackknife: each of
beyond_calibration.n_jackknife subsamples leaves out a share
beyond_calibration.jackknife_share of the structures, rebuilds every predictor on the
structures kept as the production run builds them (ranks, the marker composite and
the postsynaptic-density component over the subsample), takes the ceiling from
beyond_calibration.jackknife_splits splits of the adults, and scores each step of
the main model by the same cross-validation; the variance of the full-sample value
is (n - d) / (d N) times the spread of the subsamples' values about their mean, and
the interval is the value plus or minus 1.96 of its SD. A subsample draws no
structure twice, so no structure sits in a training and a test fold at once (a
bootstrap would let it, and flatter the fit). The uncertainty over animals is
carried by the ceiling and the half-cohort splits.

The leftover and the floor are compared on the same structures, the calibration's,
and their difference is resampled with them (beyond_calibration.paired_jackknife).

The noise null of the replication shuffles which structure is which in one half's
leftover: what two unrelated leftovers would give. It is quoted for scale only,
since the leftover of a reproducible map was always going to replicate
(beyond_density.implied_replication).

Writes, under adult_v2/ish_analysis/ in the data root:

    beyond/jackknife.csv                per subsample: the ceiling and each model's
                                        share of it
    beyond/numbers_for_the_caption.txt  the figures' numbers as sentences
    tables/numbers_beyond.csv           the numbers of analysis 4, for the text
    figures/03_beyond.png               (and .eps, as every figure)
    figures/03s_beyond_budget.png
    figures/04_beyond_where.png
    figures/11s1_leftover_genes.png
    figures/14_synaptome.png

Run by run_beyond_figures.py.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap.adult import beyond_calibration as bc
from sepmap.adult import beyond_density as bd
from sepmap.adult import synaptome
from sepmap.adult.beyond_controls import (
    CONTROLS,
    GENE_SPACE_CALIBRATION,
    GENE_SPACE_SUMMARY,
    VARIANTS,
)
from sepmap.adult.beyond_regression import PLANES, load_regression
from sepmap.config import DATA, SETTINGS
from sepmap.ish import plotting as ish_plotting
from sepmap.structures import TABLES, load_structure_set
from sepmap.volumes.per_mouse import annotation_20, structure_terms
from sepmap.volumes.to_ccf import CCF_CROP

# the permutations of the noise null and the grey of the bars; the jackknife's size
BEYOND_FIGURES = SETTINGS["beyond_figures"]
BEYOND_CALIBRATION = SETTINGS["beyond_calibration"]
ISH_ANALYSIS = SETTINGS["ish_analysis"]
ISH_FIGURES = SETTINGS["ish_figures"]

JACKKNIFE = bd.OUT / "jackknife.csv"
CAPTION = bd.OUT / "numbers_for_the_caption.txt"
NUMBERS = TABLES / "numbers_beyond.csv"
FIGURES = DATA / "adult_v2" / "ish_analysis" / "figures"

# the percentiles of a 95% interval
INTERVAL = (2.5, 97.5)


def interval(values) -> tuple[float, float]:
    """The 95% interval of `values`."""
    lo, hi = np.percentile(np.asarray(values, dtype=float), INTERVAL)
    return float(lo), float(hi)


def around(point: float, sd: float) -> tuple[float, float]:
    """The 95% interval of a value from its jackknife SD."""
    return float(point - 1.96 * sd), float(point + 1.96 * sd)


# ===== Intervals =====


def jackknife_budget(
    inputs: bd.Inputs,
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> pd.DataFrame:
    """The budget on beyond_calibration.n_jackknife subsamples of the structures.

    Per subsample: the ceiling on it and the cumulative share of each step of the
    main model (abundance, + density, + autofluorescence), each over that
    subsample's ceiling, with every predictor rebuilt on the structures kept.
    """
    n = len(inputs.structures)
    d = int(round(BEYOND_CALIBRATION["jackknife_share"] * n))
    rows = []
    for _ in range(BEYOND_CALIBRATION["n_jackknife"]):
        keep = np.sort(rng.choice(n, n - d, replace=False))
        sub = bd.restrict(inputs, [inputs.structures[i] for i in keep])
        few = [
            splits[i]
            for i in rng.choice(
                len(splits), BEYOND_CALIBRATION["jackknife_splits"], replace=False
            )
        ]
        _, explainable = bd.ceiling(sub.nano, few)
        y = bd.full_map(sub.nano)
        cov, _, _ = bd.covariates_for(sub)
        steps = bd.budget(y, cov, inputs.terms, explainable)
        rows.append(
            dict(
                ceiling=explainable,
                abundance=steps[0],
                abundance_density=steps[1],
                model=steps[2],
                left=1 - steps[2],
            )
        )
    out = pd.DataFrame(rows)
    out.insert(0, "n_left_out", d)
    out.insert(0, "n_structures", n)
    return out


def jackknife_interval(table: pd.DataFrame, column: str, point: float) -> tuple:
    """The 95% interval of a value from the jackknife table's column."""
    n, d = int(table["n_structures"].iloc[0]), int(table["n_left_out"].iloc[0])
    return around(point, bc.jackknife_sd(table[column], n, d))


def noise_band(
    nano: np.ndarray,
    xs: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> tuple[float, float]:
    """95% of the agreement of two unrelated leftovers (first split, one shuffled)."""
    a, b = splits[0]
    ra = bd.residual(bd.half_map(nano, a), xs)
    rb = bd.residual(bd.half_map(nano, b), xs)
    null = [
        spearmanr(ra, rng.permutation(rb)).statistic
        for _ in range(BEYOND_FIGURES["n_perm"])
    ]
    return interval(null)


# ===== Numbers =====


def numbers_table(n: dict) -> pd.DataFrame:
    """numbers_beyond.csv: the numbers of analysis 4 that the text quotes.

    The jackknife intervals keep six decimals: the text quotes them as whole
    percentages, and a value such as 0.154996 kept as 0.155 would round to 16%,
    not 15%.
    """
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
        ("half_agreement", round(n["half"], 4), "half-cohort maps agree (126 splits)"),
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
        ("nano_calibration_left", round(n["nano_cal_left"], 4), "nano, mean of A and B"),
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
    rows += [
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
        ("leftover_top_gene", n["top_gene"], "the gene closest to the leftover"),
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
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


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
        f"set and {n['n_adults']} adults. The main model, fixed on 8 October 2026: "
        f"{n['abundance']} + synapse density ({n['density']}) + autofluorescence, "
        f"each as x, x^2 and x^3. PSD95 density is measured in {n['psd95_measured']} "
        f"of the {n['n_fit']} structures of the fit ({n['psd95_needed']} needed), "
        f"{rule}.",
        "",
        f"Ceiling. Two halves of the cohort agree at {n['half']:.3f} (126 splits), so "
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
        f"Check rows, each one change to the main model on its own structures: {listed}.",
        "",
        f"Calibration ({n['cal_n']} structures where both halves of the Allen "
        f"experiments measure every subunit and marker). A map made only of Gria1 "
        f"and synapse density, predicted from the other half of the experiments, "
        f"leaves {floor['left_median']:.0%} ({floor['left_lo']:.0%} to "
        f"{floor['left_hi']:.0%} over {n['n_draws']} draws), its leftover replicating "
        f"at {floor['replication_median']:.2f}; {', '.join(n['markers_once'])} are "
        "measured once and so the same in both halves, so the floor errs low. The "
        f"nano map, predicted the same way on the same structures, leaves "
        f"{n['nano_cal_left']:.0%}, {n['minus_floor']:+.0%} above the floor (95% "
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
        f"Controls. {n['controls_passed']} of 7 pass. Control F: the components of "
        f"the {n['f_genes']} genes measured in every structure (most often "
        f"{n['f_k']}, picked inside each fold) predict {n['f_share']:.0%} of the "
        "reproducible map; on its own calibration the nano map leaves "
        f"{n['f_nano'][0]:.0%} to {n['f_nano'][1]:.0%}, a map made of those genes "
        f"{n['f_floor'][0]:.0%} to {n['f_floor'][1]:.0%}.",
        "",
        f"Genes against the leftover: {n['n_pass']} of {n['n_genes']} past its "
        f"spatial null at BH q < {n['q']} ({n['n_p05']} below p 0.05 before "
        f"correction); the closest is {n['top_gene']} ({n['top_rho']:+.2f}, spatial p "
        f"{n['top_p']:.4f}, q {n['top_q']:.2f}); Cacng8 {n['cacng8']['rho']:+.2f} "
        f"(p {n['cacng8']['p_spatial']:.4f}).",
        "",
        "Not shown by any of this: what the leftover is. It is the part of the map "
        "that Gria1 expression and synapse density do not predict.",
    ]


# ===== Main =====


def density_label(terms: dict[str, tuple[str, ...]], n_psd: int) -> str:
    """The main model's density terms in words, for figure 03s."""
    parts = []
    if "markers" in terms["density"]:
        parts.append(f"{len(bd.MARKERS)} marker genes")
    if "psd_pc1" in terms["density"]:
        parts.append(f"the first component of {n_psd} postsynaptic-density genes")
    parts += [f"{m} puncta" for m in terms["density"] if m in synaptome.MEASURES]
    return " and ".join(parts)


def gene_space_row(summary: pd.DataFrame, calibration: pd.DataFrame) -> dict:
    """Control F: its genes, components, nested share, replication and own floor."""
    row = summary.iloc[0]
    nano = calibration.loc[calibration["map"] == bc.NANO, "left"]
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


def plane_images(regression: pd.DataFrame) -> list[dict]:
    """The planes of figure 04: labels, and map, prediction and leftover painted."""
    names, _, _ = structure_terms()
    atlas = annotation_20("ccf")
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
    for p in PLANES:
        lab = np.asarray(atlas[p])
        out.append(
            dict(
                ccf_plane=2 * (p + CCF_CROP[0]),
                lab=lab,
                **{k: ish_plotting.paint(lab, v, names) for k, v in values.items()},
            )
        )
    return out


def draw_figures(
    inputs: bd.Inputs,
    covariates: dict[str, np.ndarray],
    y: np.ndarray,
    xs: list[np.ndarray],
    numbers: dict,
    calibration: pd.DataFrame,
    replication: pd.DataFrame,
    regression: pd.DataFrame,
    residuals: pd.DataFrame,
    genes: pd.DataFrame,
    sets: pd.DataFrame,
) -> None:
    """Draw figures 03 and 03s, 04, 11s1 and 14 from the main model and the tables."""
    s = inputs.structures

    # figure 03 and its detailed version: the map against its predictors, the budget,
    # the floor and the benchmark, the replication
    groups = ish_plotting.group_of(load_structure_set())
    held_out = dict(
        density=bd.cv_predict(y, bd.model(covariates, inputs.terms, ("density",))),
        model=bd.cv_predict(y, xs),
    )
    shared = dict(
        y=y,
        held_out=held_out,
        residual=bd.residual(y, xs),
        groups=[groups[x] for x in s],
        acronyms=[inputs.acronym.get(x, "") for x in s],
        numbers=numbers,
        calibration=calibration,
        replication=replication,
    )
    fig = ish_plotting.plot_beyond(
        **shared, save=FIGURES / ish_plotting.figure_file("beyond")
    )
    plt.close(fig)
    fig = ish_plotting.plot_beyond_budget(
        structures=s,
        gria1=covariates["Gria1"],
        **shared,
        save=FIGURES / ish_plotting.figure_file("beyond_budget"),
    )
    plt.close(fig)

    # figure 04: where the leftover lives; figure 11s1: every gene against it
    fig = ish_plotting.plot_beyond_where(
        planes=plane_images(regression),
        residuals=residuals,
        numbers=numbers,
        t_max=BEYOND_FIGURES["t_max_structures"],
        save=FIGURES / ish_plotting.figure_file("beyond_where"),
    )
    plt.close(fig)
    null = np.load(bd.LEFTOVER_NULL)
    fig = ish_plotting.plot_leftover_genes(
        genes=genes,
        sets=sets,
        numbers=numbers,
        n_surrogates=int(null["surrogates"].shape[0]),
        t_max=ISH_FIGURES["t_max"],
        save=FIGURES / ish_plotting.figure_file("leftover_genes"),
    )
    plt.close(fig)

    # figure 14: the measured synapse density against the mRNA it would replace, with
    # the check rows that put it in the model
    fig = ish_plotting.plot_synaptome(
        density=synaptome.load_density().reset_index(),
        coverage=pd.read_csv(synaptome.COVERAGE),
        markers=pd.Series(covariates["markers"], index=s),
        variants=numbers["variants"],
        min_coverage=bd.BEYOND["min_psd95_coverage"],
        save=FIGURES / ish_plotting.figure_file("synaptome"),
    )
    plt.close(fig)
    keys = ("beyond", "beyond_budget", "beyond_where", "leftover_genes", "synaptome")
    drawn = [ish_plotting.figure_file(k) for k in keys]
    print(f"figures: {', '.join(drawn)} in {FIGURES}")


def main() -> None:
    """Compute the intervals, write the numbers and draw the figures of part 1.

    One seeded generator feeds the jackknife and then the noise null.
    """
    rng = np.random.default_rng(0)
    inputs = bd.load_inputs()
    s = inputs.structures
    covariates, psd, _ = bd.covariates_for(inputs)
    xs = bd.model(covariates, inputs.terms)
    splits = bd.half_splits()
    y = bd.full_map(inputs.nano)
    map_agreement, explainable = bd.ceiling(inputs.nano, splits)
    print(f"{len(s)} structures, {len(bd.ADULTS)} adults, ceiling {explainable:.4f}")

    # the tables of steps 22 to 25
    partition = pd.read_csv(bd.PARTITION).set_index("key")
    replication = pd.read_csv(bd.REPLICATION)
    residuals = pd.read_csv(bd.RESIDUALS, keep_default_na=False, na_values=[""])
    genes = pd.read_csv(bd.LEFTOVER_GENES, keep_default_na=False, na_values=[""])
    for column in ("in_model", "gene_sets"):
        genes[column] = genes[column].fillna("")
    sets = pd.read_csv(bd.LEFTOVER_SETS, keep_default_na=False, na_values=[""])
    sets["tested"] = sets["tested"].astype(str) == "True"
    controls = pd.read_csv(CONTROLS)
    f_calibration = pd.read_csv(GENE_SPACE_CALIBRATION, keep_default_na=False)
    f_summary = pd.read_csv(GENE_SPACE_SUMMARY)
    variants = pd.read_csv(VARIANTS)
    calibration = bc.load_calibration()
    paired = bc.load_jackknife()
    regression = load_regression()

    # the intervals over structures, and the noise band of the replication
    jack = jackknife_budget(inputs, splits, rng)
    jack.to_csv(JACKKNIFE, index=False)
    noise = noise_band(inputs.nano, xs, splits, rng)
    steps = bd.budget(y, covariates, inputs.terms, explainable)
    implied = bd.implied_replication(
        float(np.mean(map_agreement)), bd.r_squared(y, xs), len(xs) + 1, len(y)
    )
    random_folds = calibration[calibration["folds"] == bc.RANDOM]
    nano_cal = random_folds[random_folds["map"] == bc.NANO].set_index("predictors_from")
    floor = bc.floor(calibration, bc.ABUNDANCE_DENSITY)
    gria1 = bc.floor(calibration, bc.GRIA1)
    nano_cal_left = float(nano_cal.loc[list(bc.HALVES), "left"].mean())
    minus_floor = nano_cal_left - floor["left_median"]
    minus_gria1 = nano_cal_left - gria1["left_median"]
    gria1_rows = random_folds[random_folds["map"] == bc.GRIA1]
    gria1_halves = gria1_rows.groupby("truth_from")["left"].median().to_dict()
    blocks = dict(
        nano=bc.floor(calibration, bc.NANO, bc.BLOCKS)["left_median"],
        floor=bc.floor(calibration, bc.ABUNDANCE_DENSITY, bc.BLOCKS)["left_median"],
        gria1=bc.floor(calibration, bc.GRIA1, bc.BLOCKS)["left_median"],
    )
    _, _, split = bc.experiment_halves(bc.per_experiment_profiles())
    markers_once = bc.measured_once(split, inputs.terms)
    cacng8 = genes.set_index("symbol").loc["Cacng8"]
    n_fit = int(inputs.rows["used"].sum())
    numbers = dict(
        n_structures=len(s),
        n_adults=len(bd.ADULTS),
        abundance=", ".join(inputs.terms["abundance"]),
        density=", ".join(inputs.terms["density"]),
        density_label=density_label(inputs.terms, len(psd)),
        psd95_measured=int(inputs.rows["psd95_measured"].sum()),
        n_fit=n_fit,
        psd95_needed=int(np.ceil(bd.BEYOND["min_psd95_coverage"] * n_fit)),
        half=float(np.mean(map_agreement)),
        ceiling=explainable,
        ceiling_ci=jackknife_interval(jack, "ceiling", explainable),
        psd_genes=len(psd),
        partition=partition["share_of_ceiling"].to_dict(),
        partition_table=partition,
        steps=steps,
        left_ci=jackknife_interval(jack, "left", 1 - steps[2]),
        cv_r2=float(partition.loc["model", "cv_r2"]),
        r2=float(partition.loc["model", "in_sample_r2"]),
        n_terms=int(partition.loc["model", "terms"]),
        rep_map=float(np.mean(map_agreement)),
        rep_left=float(replication["leftover_agreement"].mean()),
        rep_left_min=float(replication["leftover_agreement"].min()),
        implied=implied,
        noise=noise,
        cal_n=int(calibration["n_structures"].iloc[0]),
        n_draws=int((random_folds["map"] == bc.ABUNDANCE_DENSITY).sum()),
        floor=floor,
        gria1=gria1,
        nano_cal=nano_cal["left"].to_dict(),
        nano_cal_left=nano_cal_left,
        minus_floor=minus_floor,
        minus_floor_ci=jackknife_interval(paired, "nano_minus_floor", minus_floor),
        minus_gria1=minus_gria1,
        minus_gria1_ci=jackknife_interval(paired, "nano_minus_gria1", minus_gria1),
        gria1_halves=gria1_halves,
        blocks=blocks,
        markers_once=markers_once,
        variants=variants,
        controls_passed=int((controls["verdict"] == "pass").sum()),
        n_genes=len(genes),
        n_pass=int((genes["q_all"] < ISH_ANALYSIS["q"]).sum()),
        q=ISH_ANALYSIS["q"],
        top_gene=genes.iloc[0]["symbol"],
        top_rho=float(genes.iloc[0]["rho"]),
        top_p=float(genes.iloc[0]["p_spatial"]),
        top_q=float(genes.iloc[0]["q_all"]),
        n_p05=int((genes["p_spatial"] < 0.05).sum()),
        cacng8=cacng8,
        sets=sets,
        **gene_space_row(f_summary, f_calibration),
    )
    numbers_table(numbers).to_csv(NUMBERS, index=False)
    lines = caption_lines(numbers)
    with open(CAPTION, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    draw_figures(
        inputs=inputs,
        covariates=covariates,
        y=y,
        xs=xs,
        numbers=numbers,
        calibration=calibration,
        replication=replication,
        regression=regression,
        residuals=residuals,
        genes=genes,
        sets=sets,
    )
