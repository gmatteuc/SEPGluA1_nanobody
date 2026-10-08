"""The guided figures of analysis 4, with their intervals: figures 03 and 04.

adult.beyond_density, adult.beyond_controls, adult.beyond_calibration and
adult.beyond_regression write the tables. This module adds the intervals, gathers
the numbers the text quotes, and has ish.plotting draw the two figures of the guided
walk:

    11  how much of the nano map receptor mRNA and synaptic density predict, and
        whether what they leave is real: the map against Gria1, against synaptic
        density and against the whole model; the variance budget with the
        leftover's range and the calibration floor; the calibration; the leftover's
        replication across mice
    12  where the leftover lives: map, prediction and leftover on three coronal
        planes, the structures with the largest leftovers, and every gene of the
        gene table against the leftover, with the leftover's spatial null

Intervals come from resampling structures, not animals: the claim is about where in
the brain the map departs from prediction, so what would differ in a repeat of the
analysis is which structures were measurable. Each of beyond_figures.n_boot
replicates draws the structures with replacement, ranks the map and every predictor
again over the draw, scores each model by cross-validation as the point value is
scored, and takes the ceiling from beyond_figures.boot_splits splits of the adults on
the same draw. The uncertainty over animals is carried by the ceiling and the
half-cohort splits.

The bootstrap is April's. A structure drawn twice can sit in a training and in a
test fold, which flatters the fit; a replicate holds only about two thirds of the
structures as distinct ones, for a model of 22 terms, which handicaps it. On the
data of 8 October the two about cancel (the replicates' median sits by the point
value). Folds cut over distinct structures remove the first effect and keep the
second, which moves every replicate towards a larger leftover for a reason that is
the smaller sample's, not the map's, so that version is not used.

The noise null of the replication shuffles which structure is which in one half's
leftover: what two unrelated leftovers would give. It is quoted for scale only,
since the leftover of a reproducible map was always going to replicate
(beyond_density.implied_replication).

Writes, under adult_v2/ish_analysis/ in the data root:

    beyond/bootstrap.csv                per replicate: the ceiling and each model's
                                        share of it
    beyond/numbers_for_the_caption.txt  the figures' numbers as sentences
    tables/numbers_beyond.csv           the numbers of analysis 4, for the text
    figures/03_beyond_budget.png        (and .eps)
    figures/04_beyond_where.png         (and .eps)

Run by run_beyond_figures.py.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from sepmap.adult import beyond_calibration as bc
from sepmap.adult import beyond_density as bd
from sepmap.adult.beyond_controls import CONTROLS, GENE_SPACE
from sepmap.adult.beyond_regression import PLANES, load_regression
from sepmap.config import DATA, SETTINGS
from sepmap.ish import plotting as ish_plotting
from sepmap.structures import TABLES, load_structure_set
from sepmap.volumes.per_mouse import annotation_20, structure_terms
from sepmap.volumes.to_ccf import CCF_CROP

# the bootstrap replicates, the permutations of the noise null, the half-cohort splits
# used inside each replicate, and the grey of the structures' bars
BEYOND_FIGURES = SETTINGS["beyond_figures"]
ISH_ANALYSIS = SETTINGS["ish_analysis"]

BOOTSTRAP = bd.OUT / "bootstrap.csv"
CAPTION = bd.OUT / "numbers_for_the_caption.txt"
NUMBERS = TABLES / "numbers_beyond.csv"
FIGURES = DATA / "adult_v2" / "ish_analysis" / "figures"

# the percentiles of a 95% interval
INTERVAL = (2.5, 97.5)


def interval(values) -> tuple[float, float]:
    """The 95% interval of `values`."""
    lo, hi = np.percentile(np.asarray(values, dtype=float), INTERVAL)
    return float(lo), float(hi)


# ===== Intervals =====


def bootstrap_budget(
    inputs: bd.Inputs,
    covariates: dict[str, np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    rng: np.random.Generator,
) -> pd.DataFrame:
    """The budget over beyond_figures.n_boot resamples of the structures.

    Per replicate: the ceiling on the draw, the cumulative share of each step of
    the model (abundance, + density, + autofluorescence) and of the April model,
    each over that replicate's ceiling.
    """
    n = len(inputs.structures)
    rows = []
    for _ in range(BEYOND_FIGURES["n_boot"]):
        take = rng.integers(0, n, n)
        few = [
            splits[i]
            for i in rng.choice(len(splits), BEYOND_FIGURES["boot_splits"], replace=False)
        ]
        nano = inputs.nano[:, take]
        _, explainable = bd.ceiling(nano, few)
        y = bd.full_map(nano)
        cov = {k: rankdata(v[take]) for k, v in covariates.items()}
        steps = bd.budget(y, cov, explainable)
        april = bd.budget(y, cov, explainable, composite_abundance=True)
        rows.append(
            dict(
                ceiling=explainable,
                abundance=steps[0],
                abundance_density=steps[1],
                model=steps[2],
                left=1 - steps[2],
                april=april[2],
                april_left=1 - april[2],
            )
        )
    return pd.DataFrame(rows)


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
    """numbers_beyond.csv: the numbers of analysis 4 that the text quotes."""
    rows = [
        ("structures", n["n_structures"], "grey-matter structures of the fit"),
        ("adults", n["n_adults"], "adults, naive and RWS pooled"),
        ("half_agreement", round(n["half"], 4), "half-cohort maps agree (126 splits)"),
        ("ceiling", round(n["ceiling"], 4), "Spearman-Brown: reproducible share"),
        ("ceiling_lo", round(n["ceiling_ci"][0], 4), "2.5% over structures"),
        ("ceiling_hi", round(n["ceiling_ci"][1], 4), "97.5% over structures"),
        ("psd_genes", n["psd_genes"], "postsynaptic-density genes behind psd_pc1"),
    ]
    for key, what in (
        ("gria1", "Gria1 alone, bent"),
        ("abundance", "abundance alone (Gria1-4, four terms), bent"),
        ("density", "density alone (markers, psd_pc1), bent"),
        ("autofluorescence", "autofluorescence alone, bent"),
        ("straight", "the model, straight"),
        ("model", "the model: abundance + density + autofluorescence, bent"),
        ("april", "the April model (Gria1-4 composite), bent"),
    ):
        rows.append(
            (f"share_{key}", round(n["partition"][key], 4), f"{what}: CV R2 / ceiling")
        )
    rows += [
        ("budget_abundance", round(n["steps"][0], 4), "budget: abundance"),
        ("budget_density", round(n["steps"][1] - n["steps"][0], 4), "budget: + density"),
        ("budget_auto", round(n["steps"][2] - n["steps"][1], 4), "budget: + autofluo"),
        ("left", round(1 - n["steps"][2], 4), "share of the reproducible map left"),
        ("left_lo", round(n["left_ci"][0], 4), "2.5% over structures"),
        ("left_hi", round(n["left_ci"][1], 4), "97.5% over structures"),
        ("april_left", round(1 - n["partition"]["april"], 4), "left, April model"),
        ("april_left_lo", round(n["april_ci"][0], 4), "2.5% over structures"),
        ("april_left_hi", round(n["april_ci"][1], 4), "97.5% over structures"),
        ("cv_r2", round(n["cv_r2"], 4), "the model's CV R2"),
        ("in_sample_r2", round(n["r2"], 4), "the model's in-sample R2"),
        ("replication_map", round(n["rep_map"], 4), "half-cohort maps, mean over splits"),
        ("replication_leftover", round(n["rep_left"], 4), "their leftovers, mean"),
        ("replication_leftover_min", round(n["rep_left_min"], 4), "lowest split"),
        ("replication_implied", round(n["implied"], 4), "implied by ceiling and fit"),
        ("noise_lo", round(n["noise"][0], 4), "2.5% of unrelated leftovers"),
        ("noise_hi", round(n["noise"][1], 4), "97.5% of unrelated leftovers"),
        ("calibration_structures", n["cal_n"], "structures of the calibration"),
        ("floor", round(n["floor"]["left_median"], 4), "abundance-and-density map: left"),
        ("floor_lo", round(n["floor"]["left_lo"], 4), "lowest draw"),
        ("floor_hi", round(n["floor"]["left_hi"], 4), "highest draw"),
        (
            "floor_replication",
            round(n["floor"]["replication_median"], 4),
            "abundance-and-density map: leftover replicates",
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
        ("controls_passed", n["controls_passed"], "of the seven controls"),
        ("control_f_genes", n["f_genes"], "control F: genes measured everywhere"),
        ("control_f_components", n["f_k"], "control F: components of the best model"),
        ("control_f_share", round(n["f_share"], 4), "control F: CV R2 / ceiling"),
        ("control_f_replication", round(n["f_rep"], 4), "control F: leftover replicates"),
        ("leftover_genes", n["n_genes"], "genes against the leftover"),
        ("leftover_genes_pass", n["n_pass"], f"past its null, BH q < {n['q']}"),
        ("leftover_top_gene", n["top_gene"], "the gene closest to the leftover"),
        ("leftover_top_rho", round(n["top_rho"], 3), "its rho"),
        ("leftover_top_p", round(n["top_p"], 4), "its spatial p"),
        ("leftover_rho_Cacng8", round(n["cacng8"]["rho"], 3), "Cacng8 with the leftover"),
        ("leftover_p_Cacng8", round(n["cacng8"]["p_spatial"], 4), "its spatial p"),
        ("leftover_rank_Cacng8", int(n["cacng8"]["rank_all"]), "its rank"),
    ]
    for r in n["sets"].itertuples():
        rows.append(
            (
                f"leftover_set_{r.gene_set.replace(' ', '_')}",
                round(r.median_rho, 3),
                f"median rho with the leftover, {r.n_genes} genes; spatial p "
                + (f"{r.p_spatial:.4f}" if r.tested else "not tested"),
            )
        )
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


def caption_lines(n: dict) -> list[str]:
    """numbers_for_the_caption.txt: the figures' numbers as sentences."""
    floor = n["floor"]
    gria1 = n["gria1"]
    return [
        f"Analysis 4, on {n['n_structures']} grey-matter structures of the declared "
        f"set and {n['n_adults']} adults.",
        "",
        f"Ceiling. Two halves of the cohort agree at {n['half']:.3f} (126 splits), so "
        f"{n['ceiling']:.1%} of the map is reproducible (Spearman-Brown; 95% over "
        f"structures {n['ceiling_ci'][0]:.1%} to {n['ceiling_ci'][1]:.1%}).",
        "",
        f"Budget. On structures the fit has not seen, the four subunits predict "
        f"{n['steps'][0]:.0%} of the reproducible map, synaptic density adds "
        f"{n['steps'][1] - n['steps'][0]:.0%} and autofluorescence "
        f"{n['steps'][2] - n['steps'][1]:.0%}; {1 - n['steps'][2]:.0%} is not "
        f"predicted by receptor mRNA or synaptic density (95% over structures "
        f"{n['left_ci'][0]:.0%} to {n['left_ci'][1]:.0%}). Gria1 alone predicts "
        f"{n['partition']['gria1']:.0%}, density alone {n['partition']['density']:.0%}. "
        f"With the April composite of the subunits, {1 - n['partition']['april']:.0%} "
        "is left.",
        "",
        f"Calibration ({n['cal_n']} structures where both halves of the Allen "
        f"experiments measure every subunit and marker). A map made only of receptor "
        f"mRNA and synaptic density, predicted from the other half of the "
        f"experiments, leaves {floor['left_median']:.0%} ({floor['left_lo']:.0%} to "
        f"{floor['left_hi']:.0%} over {n['n_draws']} draws), its leftover replicating "
        f"at {floor['replication_median']:.2f}. The nano map, predicted the same way, "
        f"leaves {min(n['nano_cal'][h] for h in bc.HALVES):.0%} to "
        f"{max(n['nano_cal'][h] for h in bc.HALVES):.0%}. A map of one Allen experiment "
        f"of Gria1 leaves {gria1['left_median']:.0%} ({gria1['left_lo']:.0%} to "
        f"{gria1['left_hi']:.0%}), replicating at {gria1['replication_median']:.2f}.",
        "",
        f"Replication. The leftover of one half-cohort agrees with the other's at "
        f"{n['rep_left']:.3f} (lowest split {n['rep_left_min']:.3f}); the map's "
        f"reliability and the fit alone imply {n['implied']:.3f}; two unrelated "
        f"leftovers fall between {n['noise'][0]:+.2f} and {n['noise'][1]:+.2f}.",
        "",
        f"Controls. {n['controls_passed']} of 7 pass. Control F: {n['f_k']} components "
        f"of {n['f_genes']} genes predict {n['f_share']:.0%} of the reproducible map, "
        f"and their leftover replicates at {n['f_rep']:.2f}.",
        "",
        f"Genes against the leftover: {n['n_pass']} of {n['n_genes']} past its "
        f"spatial null at BH q < {n['q']}; the closest is {n['top_gene']} "
        f"({n['top_rho']:+.2f}, spatial p {n['top_p']:.3f}); Cacng8 "
        f"{n['cacng8']['rho']:+.2f} (p {n['cacng8']['p_spatial']:.2f}).",
        "",
        "Not shown by any of this: what the leftover is. It is the part of the map "
        "that receptor mRNA and synaptic density do not predict.",
    ]


# ===== Main =====


def gene_space_row(controls: pd.DataFrame, curve: pd.DataFrame) -> dict:
    """Control F's best model: its components, genes, share and replication."""
    row = controls.set_index("control").loc["F whole gene space", "number"]
    best = curve.loc[curve["cv_r2"].idxmax()]
    genes = int(row.split(" genes")[0])
    rep = float(row.rsplit("replication ", 1)[1])
    return dict(
        f_genes=genes,
        f_k=int(best["n_components"]),
        f_share=float(best["share_of_ceiling"]),
        f_rep=rep,
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


def main() -> None:
    """Compute the intervals, write the numbers and draw figures 03 and 04.

    One seeded generator feeds the bootstrap and then the noise null.
    """
    rng = np.random.default_rng(0)
    inputs = bd.load_inputs()
    s = inputs.structures
    covariates, psd, _ = bd.build_covariates(inputs.expr, inputs.role, inputs.auto, s)
    xs = bd.model(covariates)
    splits = bd.half_splits()
    y = bd.full_map(inputs.nano)
    map_agreement, explainable = bd.ceiling(inputs.nano, splits)
    print(f"{len(s)} structures, {len(bd.ADULTS)} adults, ceiling {explainable:.4f}")

    # the tables of steps 21 to 24
    partition = pd.read_csv(bd.PARTITION).set_index("key")
    replication = pd.read_csv(bd.REPLICATION)
    residuals = pd.read_csv(bd.RESIDUALS, keep_default_na=False, na_values=[""])
    genes = pd.read_csv(bd.LEFTOVER_GENES, keep_default_na=False, na_values=[""])
    for column in ("in_model", "gene_sets"):
        genes[column] = genes[column].fillna("")
    sets = pd.read_csv(bd.LEFTOVER_SETS, keep_default_na=False, na_values=[""])
    sets["tested"] = sets["tested"].astype(str) == "True"
    controls = pd.read_csv(CONTROLS)
    curve = pd.read_csv(GENE_SPACE)
    calibration = bc.load_calibration()
    regression = load_regression()

    # the intervals over structures, and the noise band of the replication
    boot = bootstrap_budget(inputs, covariates, splits, rng)
    boot.to_csv(BOOTSTRAP, index=False)
    noise = noise_band(inputs.nano, xs, splits, rng)
    steps = bd.budget(y, covariates, explainable)
    implied = bd.implied_replication(
        float(np.mean(map_agreement)), bd.r_squared(y, xs), len(xs) + 1, len(y)
    )
    nano_cal = calibration[calibration["map"] == bc.NANO].set_index("predictors_from")
    cacng8 = genes.set_index("symbol").loc["Cacng8"]
    numbers = dict(
        n_structures=len(s),
        n_adults=len(bd.ADULTS),
        n_markers=len(bd.MARKERS),
        half=float(np.mean(map_agreement)),
        ceiling=explainable,
        ceiling_ci=interval(boot["ceiling"]),
        psd_genes=len(psd),
        partition=partition["share_of_ceiling"].to_dict(),
        steps=steps,
        april_steps=bd.budget(y, covariates, explainable, composite_abundance=True),
        left_ci=interval(boot["left"]),
        april_ci=interval(boot["april_left"]),
        cv_r2=float(partition.loc["model", "cv_r2"]),
        r2=float(partition.loc["model", "in_sample_r2"]),
        rep_map=float(np.mean(map_agreement)),
        rep_left=float(replication["leftover_agreement"].mean()),
        rep_left_min=float(replication["leftover_agreement"].min()),
        implied=implied,
        noise=noise,
        cal_n=int(calibration["n_structures"].iloc[0]),
        n_draws=int((calibration["map"] == bc.ABUNDANCE_DENSITY).sum()),
        floor=bc.floor(calibration, bc.ABUNDANCE_DENSITY),
        gria1=bc.floor(calibration, bc.GRIA1),
        nano_cal=nano_cal["left"].to_dict(),
        controls_passed=int((controls["verdict"] == "pass").sum()),
        n_genes=len(genes),
        n_pass=int((genes["q_all"] < ISH_ANALYSIS["q"]).sum()),
        q=ISH_ANALYSIS["q"],
        top_gene=genes.iloc[0]["symbol"],
        top_rho=float(genes.iloc[0]["rho"]),
        top_p=float(genes.iloc[0]["p_spatial"]),
        cacng8=cacng8,
        sets=sets,
        **gene_space_row(controls, curve),
    )
    numbers_table(numbers).to_csv(NUMBERS, index=False)
    lines = caption_lines(numbers)
    with open(CAPTION, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # figure 03: the map against its predictors, the budget, the calibration, the
    # replication
    set_table = load_structure_set()
    groups = ish_plotting.group_of(set_table)
    held_out = dict(
        density=bd.cv_predict(y, bd.model(covariates, ("density",))),
        model=bd.cv_predict(y, xs),
    )
    fig = ish_plotting.plot_beyond_budget(
        structures=s,
        y=y,
        gria1=covariates["Gria1"],
        held_out=held_out,
        residual=bd.residual(y, xs),
        groups=[groups[x] for x in s],
        acronyms=[inputs.acronym.get(x, "") for x in s],
        numbers=numbers,
        calibration=calibration,
        replication=replication,
        save=FIGURES / ish_plotting.figure_file("beyond_budget"),
    )
    plt.close(fig)

    # figure 04: where the leftover lives, and the genes against it
    null = np.load(bd.LEFTOVER_NULL)
    fig = ish_plotting.plot_beyond_where(
        planes=plane_images(regression),
        residuals=residuals,
        genes=genes,
        sets=sets,
        numbers=numbers,
        n_surrogates=int(null["surrogates"].shape[0]),
        t_max=BEYOND_FIGURES["t_max_structures"],
        t_max_genes=BEYOND_FIGURES["t_max_genes"],
        save=FIGURES / ish_plotting.figure_file("beyond_where"),
    )
    plt.close(fig)
    print(f"figures 03 and 04 in {FIGURES}")
