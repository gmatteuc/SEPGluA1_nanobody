"""What receptor mRNA and synaptic density leave of the adult map (analysis 4).

The question: is the adult nano map just how much receptor mRNA a structure has, or
just how many synapses it has? Ranking genes against the map cannot tell, because a
map of pure synaptic density also puts the postsynaptic genes on top (analysis 3).
So the question is turned round: predict the map as well as possible from receptor
mRNA and synaptic density, and ask how much is left, whether the leftover is more
than one Allen map disagreeing with another, and what it looks like.

    nano map  ~  receptor mRNA + synaptic density + autofluorescence  ->  leftover

The steps, each printing its numbers (the working figures fig0 to fig3 draw them):

    0  the structures   the declared set (sepmap.structures), less the structures
                        where a subunit or marker gene has no ISH value
    1  the ceiling      how reproducible the map is: nothing can predict variance
                        that is not there, so every R2 is read against it
    2  the budget       what receptor mRNA, then synaptic density, then
                        autofluorescence predict, on structures the fit has not seen
    3  the leftover     whether it replicates across mice, and the value the ceiling
                        and the fit alone imply it would replicate at
    4  where it lives   which structures carry it and how consistently across
                        half-cohorts, and which genes' maps follow it, against the
                        leftover's own spatial null

Everything is one number per structure, on ranks:

    y            the ten adults' mean zref per structure (the declared reference,
                 A1), ranked
    abundance    the four AMPA receptor subunits Gria1, Gria2, Gria3 and Gria4, each
                 its own term. The April model averaged their ranks into one
                 composite; Gria4 runs against the map, so the average diluted
                 Gria1. The composite stays as a check row (the April model)
    density      the mean rank of the synaptic marker genes ([beyond] markers), and
                 the first principal component of the postsynaptic-density genes of
                 the ontology panel (role control_psd) that have a value in every
                 structure: whatever those genes have most in common
    autofluo     the ten adults' mean autofluorescence zref, the only predictor
                 measured in the same sections as the map
    bent         every predictor enters as x, x^2 and x^3. The rank relationships
                 are curved, and a straight line would leave the curvature in the
                 leftover; control E (adult.beyond_controls) checks that bending
                 further buys nothing
    CV R2        each structure predicted from a fit that never saw it: the
                 structures are shuffled once with a fixed seed, cut into five
                 folds, and each fold predicted from the other four. In-sample R2
                 always grows with terms, so the held-out one is quoted

The ceiling: the ten adults are split into two fives every possible way (126
splits), the two half-maps are correlated, and Spearman-Brown turns their mean
agreement r into the reliability of the ten-adult map,

    ceiling = 2 r / (1 + r)

A reliability is already a share of variance (true over observed), so the ceiling is
not squared; the April code squared it once more. Then

    explained    CV R2 / ceiling        the share of the reproducible map predicted
    left         1 - CV R2 / ceiling    the share not predicted: the leftover

The leftover is quoted as a range (adult.beyond_figures resamples the structures)
beside the calibration floor of adult.beyond_calibration: what the same model leaves
of a map that is exactly receptor mRNA and synaptic density, measured with other
Allen experiments. Part of any leftover is one Allen map disagreeing with another,
and the floor says how much.

The leftover replicating across mice is not separate evidence. If the map
replicates, what is left of it once a smooth fit is removed must replicate too;
implied_replication gives the value that the ceiling and the fit alone predict, and
the figures show it beside the observed one.

The words: the leftover is "not predicted by receptor mRNA or synaptic density",
never "beyond gene expression" (control F gives the model the whole gene panel).
Nothing here measures what the leftover is. The surface fraction of the receptor is
the reading the data support; translation, turnover, subunit composition or
nanobody access would land in the same place, and a total-GluA1 stain on the same
brains is what would tell them apart.

The choices, and why:

    structures, not voxels  the gene data exist only at 200 um and in another mouse,
                            so at the voxel level the comparison would be mostly
                            registration error
    ranks everywhere        Allen expression energy has an arbitrary scale per
                            experiment; ranks assume only that more mRNA gives more
                            signal
    the declared set        grey matter seen in all ten adults (A1): fibre tracts
                            have no synapses, and structures seen in a few mice do
                            not replicate
    naive and RWS pooled    the groups differ in whisker experience, not in what the
                            stain is; control D checks that the leftover does not
                            depend on the group

Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    structures_used.csv        every structure of the adult table, used or not, why
    variance_partition.csv     what each model predicts, against the ceiling
    replication.csv            per split of the adults, how well the two half-maps
                               and their two leftovers agree
    residual_by_structure.csv  per structure: map, prediction, leftover, and the
                               share of half-cohort leftovers with the same sign
    leftover_genes.csv         every gene's rho with the leftover, its spatial p
    leftover_sets.csv          the gene sets of analysis 3 against the leftover
    leftover_null.npz          the leftover's surrogates and every gene's rho with them
    fig0_structures.png, fig1_ceiling.png, fig2_covariates.png, fig3_residual.png

Run by run_beyond_density.py.
"""

import itertools
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from sepmap.adult import profiles
from sepmap.config import DATA, SETTINGS
from sepmap.ish import gene_ranking, gene_sets, gene_table, spatial_null
from sepmap.plotting import DARK_BLUE, RED, tidy
from sepmap.structures import load_centroids, load_structure_set
from sepmap.volumes.cohort import NAIVE, RWS

# the half-cohort size and the gene sets standing in for receptor mRNA and synaptic
# density; the agreement a leftover needs to count as replicating
BEYOND = SETTINGS["beyond"]
BEYOND_CONTROLS = SETTINGS["beyond_controls"]

OUT = DATA / "adult_v2" / "ish_analysis" / "beyond"
STRUCTURES_USED = OUT / "structures_used.csv"
PARTITION = OUT / "variance_partition.csv"
REPLICATION = OUT / "replication.csv"
RESIDUALS = OUT / "residual_by_structure.csv"
LEFTOVER_GENES = OUT / "leftover_genes.csv"
LEFTOVER_SETS = OUT / "leftover_sets.csv"
LEFTOVER_NULL = OUT / "leftover_null.npz"

# the ten adults, naive then rws, and the rows of each group in the adult matrices
ADULTS = profiles.ADULTS
NAIVE_ROWS = [ADULTS.index(m) for m in NAIVE]
RWS_ROWS = [ADULTS.index(m) for m in RWS]

# the genes standing in for the two explanations, as tuples
SUBUNITS = tuple(BEYOND["subunits"])
MARKERS = tuple(BEYOND["markers"])

# the groups of predictors in the order the budget adds them; abundance is the four
# subunits as separate terms
GROUPS = {
    "abundance": SUBUNITS,
    "density": ("markers", "psd_pc1"),
    "autofluorescence": ("autofluo",),
}
ORDER = tuple(GROUPS)

# the ontology panel's role whose genes make the density component
PSD_ROLE = "control_psd"


# ===== Loading =====


@dataclass
class Inputs:
    """What every step of analysis 4 reads, loaded once.

    `nano` and `auto` are adults x structures (zref of the nano and of the
    autofluorescence channel), in the order of `structures`; `expr` holds every
    gene's merged profile (mean ranks per structure), `role` its role in the
    ontology panel ('' outside it); `rows` lists every structure of the adult table,
    used by the fit or not, with the reason.
    """

    structures: list[str]
    nano: np.ndarray
    auto: np.ndarray
    expr: dict[str, dict[str, float]]
    role: dict[str, str]
    p9_genes: set[str]
    division: dict[str, str]
    acronym: dict[str, str]
    rows: pd.DataFrame


def structure_rows(
    set_table: pd.DataFrame, expr: dict[str, dict[str, float]], auto_ok: pd.Series
) -> pd.DataFrame:
    """Every structure of the adult table: used by the fit or not, and why.

    A declared structure is used when every subunit and marker gene has an ISH
    value there (ish.min_voxels voxels in at least one of its experiments) and every
    adult has autofluorescence above background (`auto_ok`, by structure).
    Columns: structure, acronym, division, used, reason, missing_genes.
    """
    needed = SUBUNITS + MARKERS
    rows = []
    for r in set_table.itertuples():
        missing = [g for g in needed if r.structure not in expr.get(g, {})]
        if not r.in_set:
            reason = f"not in the declared set: {r.reason}"
            missing = []
        elif missing:
            reason = "a subunit or marker gene has no ISH value there"
        elif not bool(auto_ok.get(r.structure, False)):
            reason = "autofluorescence not above background in every adult"
        else:
            reason = ""
        rows.append(
            dict(
                structure=r.structure,
                acronym=r.acronym,
                division=r.division,
                used=reason == "",
                reason=reason,
                missing_genes=" ".join(missing),
            )
        )
    return pd.DataFrame(rows)


def load_inputs() -> Inputs:
    """The adult profiles, the gene table's profiles and the structures of the fit."""
    per_mouse = profiles.load_per_mouse()
    set_table = load_structure_set()
    expr = gene_table.load_profiles()
    genes = gene_table.per_gene(gene_table.load_gene_table())
    role = dict(zip(genes["symbol"], genes["ontology_role"]))
    p9_genes = set(genes.loc[genes["p9_gene"], "symbol"])

    # autofluorescence above background in every adult, by structure
    auto = per_mouse.pivot(index="mouse", columns="structure", values="zref_auto")
    auto_ok = auto.reindex(index=ADULTS).notna().all(axis=0)
    rows = structure_rows(set_table, expr, auto_ok)
    used = sorted(rows.loc[rows["used"], "structure"])
    return Inputs(
        structures=used,
        nano=gene_ranking.adult_matrix(per_mouse, "zref_nano", used, ADULTS),
        auto=gene_ranking.adult_matrix(per_mouse, "zref_auto", used, ADULTS),
        expr=expr,
        role=role,
        p9_genes=p9_genes,
        division=dict(zip(set_table["structure"], set_table["division"])),
        acronym=dict(zip(set_table["structure"], set_table["acronym"])),
        rows=rows,
    )


# ===== Statistics =====


def composite(
    genes: Sequence[str], expr: dict[str, dict[str, float]], structures: list[str]
) -> np.ndarray:
    """One predictor from a set of genes: the mean of their rank profiles.

    Ranks rather than values because each Allen experiment carries its own
    arbitrary intensity scale, so raw numbers are not comparable between genes
    even when their orderings are.
    """
    return np.mean([rankdata([expr[g][s] for s in structures]) for g in genes], axis=0)


def first_pc(
    genes: Sequence[str], expr: dict[str, dict[str, float]], structures: list[str]
) -> tuple[np.ndarray, float]:
    """The dominant shared axis of a gene set, and the share of variance it carries.

    Used for the postsynaptic-density genes: rather than naming a handful of
    markers and arguing about the choice, take whatever those genes have most in
    common and call that the density axis.
    """
    m = np.array([rankdata([expr[g][s] for s in structures]) for g in genes])
    m = (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, keepdims=True)
    _, sv, vt = np.linalg.svd(m - m.mean(axis=0), full_matrices=False)
    pc = vt[0]

    # point it the same way as the genes
    if np.corrcoef(pc, m.mean(axis=0))[0, 1] < 0:
        pc = -pc
    return pc, float(sv[0] ** 2 / (sv**2).sum())


def residual(y: np.ndarray, predictors: Sequence[np.ndarray]) -> np.ndarray:
    """What is left of y after least squares on the predictors, plus an intercept."""
    design = np.column_stack(list(predictors) + [np.ones(len(y))])
    return y - design @ np.linalg.lstsq(design, y, rcond=None)[0]


def r_squared(y: np.ndarray, predictors: Sequence[np.ndarray]) -> float:
    """Variance explained on the data the fit was made from, so optimistic."""
    return float(1 - residual(y, predictors).var() / y.var())


def flexible(predictors: Sequence[np.ndarray]) -> list[np.ndarray]:
    """The same predictors, allowed to bend: x, x^2 and x^3 of each.

    A straight line through two rank variables assumes the relationship is not
    only monotone but evenly paced; control E of adult.beyond_controls shows that
    assumption is wrong here, so the quoted model bends.
    """
    return list(predictors) + [x**2 for x in predictors] + [x**3 for x in predictors]


def cv_predict(
    y: np.ndarray, predictors: Sequence[np.ndarray], folds: int = 5
) -> np.ndarray:
    """Each structure's prediction from a fit that never saw it.

    Structures are shuffled once with a fixed seed, cut into folds, and each fold
    predicted from a least-squares fit on the others.
    """
    n = len(y)
    order = np.random.default_rng(0).permutation(n)
    design = np.column_stack(list(predictors) + [np.ones(n)])
    predicted = np.empty(n)
    for k in range(folds):
        test = order[k::folds]
        train = np.setdiff1d(order, test)
        predicted[test] = (
            design[test] @ np.linalg.lstsq(design[train], y[train], rcond=None)[0]
        )
    return predicted


def cv_r2(y: np.ndarray, predictors: Sequence[np.ndarray], folds: int = 5) -> float:
    """Variance explained on structures the fit has never seen (cv_predict).

    The number to quote once a model has many terms: adding predictors always
    improves the in-sample fit, and only held-out structures can say whether it
    improved the prediction.
    """
    predicted = cv_predict(y, predictors, folds)
    return float(1 - np.var(y - predicted) / np.var(y))


def half_splits() -> list[tuple[list[int], list[int]]]:
    """Every way of cutting ten animals into two fives, each split counted once.

    A split is counted once by keeping only the halves that hold animal 0.
    """
    return [
        (list(p), [i for i in range(len(ADULTS)) if i not in p])
        for p in itertools.combinations(range(len(ADULTS)), BEYOND["half"])
        if 0 in p
    ]


def spearman_brown(r: float) -> float:
    """Reliability of a whole cohort from the agreement of its two halves; NaN at -1."""
    return 2 * r / (1 + r) if r > -1 else float("nan")


def replicates(agreement: float) -> bool:
    """Whether two half-cohort maps or leftovers agree at beyond_controls.replication."""
    return agreement >= BEYOND_CONTROLS["replication"]


def half_map(matrix: np.ndarray, rows: Sequence[int]) -> np.ndarray:
    """The ranked mean map of the adults in `rows` (matrix: adults x structures)."""
    return rankdata(matrix[list(rows)].mean(axis=0))


def full_map(matrix: np.ndarray) -> np.ndarray:
    """The ranked mean map of every adult."""
    return half_map(matrix, range(matrix.shape[0]))


def ceiling(
    matrix: np.ndarray, splits: list[tuple[list[int], list[int]]]
) -> tuple[list[float], float]:
    """The two half-maps' agreement over every split, and the ceiling.

    The ceiling is the Spearman-Brown reliability of the mean agreement: the share
    of the ten-adult map's variance that is reproducible.
    """
    agreement = [
        float(spearmanr(half_map(matrix, a), half_map(matrix, b)).statistic)
        for a, b in splits
    ]
    return agreement, spearman_brown(float(np.mean(agreement)))


def leftover_agreement(
    matrix: np.ndarray,
    predictors: Sequence[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
) -> list[float]:
    """How well the leftover of one half-cohort matches the other's, per split.

    Each half's ranked map is fitted to the predictors on its own, and the two
    leftovers are compared by Spearman.
    """
    return [
        float(
            spearmanr(
                residual(half_map(matrix, a), predictors),
                residual(half_map(matrix, b), predictors),
            ).statistic
        )
        for a, b in splits
    ]


def implied_replication(half_agreement: float, r2: float, terms: int, n: int) -> float:
    """The leftover's replication that the map's reliability and the fit alone imply.

    A half-map is signal (a share h, the half-against-half agreement) plus noise
    (1 - h). The fit removes a share r2 of the variance, all of it signal, and
    absorbs terms / n of the noise (terms counts the intercept). Two half leftovers
    then share the signal left, h - r2, out of h - r2 plus the noise left, so

        implied = (h - r2) / ((h - r2) + (1 - h) (1 - terms / n))

    A leftover that replicates near this value carries no evidence beyond the
    map's own reliability; only the ceiling and the calibration speak to it.
    """
    signal = half_agreement - r2
    noise = (1 - half_agreement) * (1 - terms / n)
    return float(signal / (signal + noise)) if signal > 0 else float("nan")


# ===== The model =====


def build_covariates(
    expr: dict[str, dict[str, float]],
    role: dict[str, str],
    auto: np.ndarray,
    structures: list[str],
) -> tuple[dict[str, np.ndarray], list[str], float]:
    """Every predictor by name, built in one place so every module builds them alike.

    Gria1 to Gria4 each ranked, their April composite (abundance_composite), the
    marker composite, psd_pc1 and autofluo (the ranked mean of `auto`, adults x
    structures). Returns the predictors, the postsynaptic-density genes behind
    psd_pc1 and the share of their variance it carries.
    """
    psd = sorted(
        g
        for g in expr
        if role.get(g, "") == PSD_ROLE and all(s in expr[g] for s in structures)
    )
    psd_pc1, share = first_pc(psd, expr, structures)
    out = {g: rankdata([expr[g][s] for s in structures]) for g in SUBUNITS}
    out["abundance_composite"] = composite(SUBUNITS, expr, structures)
    out["markers"] = composite(MARKERS, expr, structures)
    out["psd_pc1"] = psd_pc1
    out["autofluo"] = rankdata(auto.mean(axis=0))
    return out, psd, share


def predictors(
    covariates: dict[str, np.ndarray],
    groups: Sequence[str] = ORDER,
    composite_abundance: bool = False,
) -> list[np.ndarray]:
    """The straight predictors of the groups named, in the order of the budget.

    With `composite_abundance`, abundance is the April composite, one term.
    """
    names = []
    for group in groups:
        if group == "abundance" and composite_abundance:
            names.append("abundance_composite")
        else:
            names.extend(GROUPS[group])
    return [covariates[n] for n in names]


def model(
    covariates: dict[str, np.ndarray],
    groups: Sequence[str] = ORDER,
    composite_abundance: bool = False,
    bend: bool = True,
) -> list[np.ndarray]:
    """The columns of a model: the groups' predictors, bent unless `bend` is False.

    model(covariates) is the quoted model: the four subunits, the markers, psd_pc1
    and autofluorescence, each as x, x^2 and x^3.
    """
    xs = predictors(covariates, groups, composite_abundance)
    return flexible(xs) if bend else xs


def budget(
    y: np.ndarray,
    covariates: dict[str, np.ndarray],
    explainable: float,
    composite_abundance: bool = False,
) -> list[float]:
    """The cumulative share of the explainable variance: abundance, + density, + auto.

    Each step is the CV R2 of the bent model so far over the ceiling. The leftover
    is 1 minus the last.
    """
    return [
        cv_r2(y, model(covariates, ORDER[:k], composite_abundance)) / explainable
        for k in range(1, len(ORDER) + 1)
    ]


# ===== Drawing =====


def save(fig: plt.Figure, name: str) -> None:
    """Save a working figure as `name` in the output folder at 200 dpi, and close it."""
    path = OUT / name
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  -> {path}")


# ===== Steps =====


def figure_structures(rows: pd.DataFrame) -> None:
    """Draw fig0: the structures used per division, and those left out per reason."""
    fig, axes = plt.subplots(
        1, 2, figsize=(11.0, 3.9), gridspec_kw=dict(width_ratios=[1, 1.25])
    )
    counts = rows[rows["used"]].groupby("division").size().sort_values(ascending=False)
    axes[0].bar(
        range(len(counts)), counts.to_numpy(), color="0.65", edgecolor="0.25", lw=0.5
    )
    axes[0].set_xticks(range(len(counts)))
    axes[0].set_xticklabels(counts.index, fontsize=7, rotation=45, ha="right")
    axes[0].set_ylabel("structures used", fontsize=8)
    axes[0].set_title(
        f"what the fit runs on\n{int(rows['used'].sum())} grey-matter structures",
        fontsize=9,
    )
    tidy(axes[0])

    # the reasons, the declared set's own reasons pooled into one
    why = rows.loc[~rows["used"], "reason"].map(
        lambda r: "not in the declared set" if r.startswith("not in the declared") else r
    )
    why_counts = why.value_counts().sort_values()
    axes[1].barh(
        range(len(why_counts)),
        why_counts.to_numpy(),
        color=RED,
        edgecolor="0.25",
        linewidth=0.5,
        alpha=0.85,
    )
    axes[1].set_yticks(range(len(why_counts)))
    axes[1].set_yticklabels([w[:52] for w in why_counts.index], fontsize=7)
    axes[1].set_xlabel("structures left out", fontsize=8)
    axes[1].set_title(
        "and what it leaves out\ndecided by rule, before any fitting", fontsize=9
    )
    tidy(axes[1])
    fig.tight_layout()
    save(fig, "fig0_structures.png")


def step0_structures(inputs: Inputs) -> None:
    """Write which structures the fit runs on, and why the others are left out."""
    print("\nSTEP 0  which structures the fit runs on")
    rows = inputs.rows
    rows.to_csv(STRUCTURES_USED, index=False)
    out = rows[~rows["used"]]
    print(f"  {len(rows)} structures in the adult table, {int(rows['used'].sum())} used")
    for reason, n in out["reason"].value_counts().items():
        print(f"    {n:3d}  {reason}")
    lacking = defaultdict(int)
    for text in out["missing_genes"]:
        for g in text.split():
            lacking[g] += 1
    if lacking:
        print(
            "  genes without a value in a declared structure: "
            + ", ".join(f"{g} {n}" for g, n in sorted(lacking.items()))
        )
    figure_structures(rows)


def step1_ceiling(
    nano: np.ndarray, splits: list[tuple[list[int], list[int]]]
) -> tuple[float, list[float]]:
    """How reproducible the map itself is: nothing below can beat this."""
    print("\nSTEP 1  the ceiling: how much of this map is reproducible")
    agreement, explainable = ceiling(nano, splits)
    half = float(np.mean(agreement))
    print(f"  over {len(splits)} five-against-five splits of the ten adults:")
    print(f"    half-cohort against half-cohort   rho = {half:.3f}")
    print(f"    Spearman-Brown, all ten adults    {explainable:.3f}")
    print(f"  so {explainable:.1%} of the map's variance is reproducible, and every")
    print("  R2 below is read against that rather than against 1.")

    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    ax.hist(agreement, bins=25, color="0.7", edgecolor="0.35", linewidth=0.4)
    ax.axvline(half, color=RED, lw=1.8)
    ax.set_xlabel("Spearman between the two half-cohort maps", fontsize=8)
    ax.set_ylabel(f"splits of ten adults ({len(splits)})", fontsize=8)
    verdict = "is reproducible" if replicates(half) else "does not reproduce"
    ax.set_title(
        f"Step 1. the map {verdict}\n"
        f"half-cohorts agree at {half:.3f}; ceiling (Spearman-Brown) {explainable:.3f}",
        fontsize=9,
    )
    tidy(ax)
    fig.tight_layout()
    save(fig, "fig1_ceiling.png")
    return explainable, agreement


def partition_models(c: dict[str, np.ndarray]) -> list[tuple[str, str, list]]:
    """The models of variance_partition.csv as (key, label, columns); "model" is the
    quoted one."""
    return [
        ("gria1", "Gria1 alone", flexible([c["Gria1"]])),
        ("abundance", "abundance: Gria1-4, four terms", model(c, ("abundance",))),
        ("density", "density: markers and psd_pc1", model(c, ("density",))),
        ("autofluorescence", "autofluorescence alone", model(c, ("autofluorescence",))),
        ("abundance_density", "abundance + density", model(c, ("abundance", "density"))),
        (
            "straight",
            "abundance + density + autofluorescence, straight",
            model(c, bend=False),
        ),
        ("model", "abundance + density + autofluorescence", model(c)),
        (
            "april",
            "April model: Gria1-4 composite + density + autofluorescence",
            model(c, composite_abundance=True),
        ),
    ]


def partition_table(
    y: np.ndarray, models: list[tuple[str, str, list]], explainable: float
) -> pd.DataFrame:
    """variance_partition.csv: each model's terms, R2 held out and in-sample, shares."""
    rows = []
    for key, label, xs in models:
        cv = cv_r2(y, xs)
        rows.append(
            dict(
                key=key,
                model=label,
                terms=len(xs) + 1,
                bent=key != "straight",
                cv_r2=cv,
                in_sample_r2=r_squared(y, xs),
                share_of_ceiling=cv / explainable,
                left=1 - cv / explainable,
            )
        )
    return pd.DataFrame(rows)


def figure_covariates(table: pd.DataFrame, explainable: float) -> None:
    """Draw fig2: each model's cross-validated R2 against the ceiling."""
    fig, ax = plt.subplots(figsize=(7.8, 4.5))
    ax.barh(np.arange(len(table)), table["cv_r2"], color="0.65", edgecolor="0.25", lw=0.5)
    ax.axvline(explainable, color=RED, lw=1.8)
    ax.annotate(
        f"ceiling {explainable:.3f}\n(the map's own reliability)",
        (explainable, len(table) - 0.4),
        color=RED,
        fontsize=7.5,
        ha="right",
        va="top",
        xytext=(-6, 0),
        textcoords="offset points",
    )
    ax.set_yticks(np.arange(len(table)))
    ax.set_yticklabels(table["model"], fontsize=7.5)
    ax.set_xlim(min(0.0, float(table["cv_r2"].min()) - 0.02), 1.02)
    ax.set_xlabel(
        "variance of the map predicted on held-out structures (cross-validated R2)",
        fontsize=8,
    )
    quoted = float(table.set_index("key").loc["model", "cv_r2"])
    ax.set_title(
        f"Step 2. what receptor mRNA and synaptic density predict\n"
        f"the model, bent: {quoted / explainable:.1%} of the ceiling, "
        f"{1 - quoted / explainable:.1%} left",
        fontsize=9,
    )
    tidy(ax)
    fig.tight_layout()
    save(fig, "fig2_covariates.png")


def step2_budget(
    y: np.ndarray,
    covariates: dict[str, np.ndarray],
    psd: list[str],
    share: float,
    explainable: float,
) -> pd.DataFrame:
    """What receptor mRNA, synaptic density and autofluorescence predict."""
    print("\nSTEP 2  what receptor mRNA and synaptic density predict")
    print(
        f"  psd_pc1 is the first component of {len(psd)} postsynaptic-density genes, "
        f"carrying {share:.1%} of their variance"
    )
    print("  each predictor on its own (Spearman with the map):")
    for key in SUBUNITS + ("abundance_composite", "markers", "psd_pc1", "autofluo"):
        print(f"    {key:20s} rho {spearmanr(y, covariates[key]).statistic:+.3f}")

    table = partition_table(y, partition_models(covariates), explainable)
    table.to_csv(PARTITION, index=False)
    print("  the models, scored on structures the fit has not seen:")
    for r in table.itertuples():
        print(
            f"    {r.model:62s} CV R2 {r.cv_r2:+.3f}   in-sample {r.in_sample_r2:.3f}"
            f"   {r.share_of_ceiling:6.1%} of the ceiling"
        )
    steps = budget(y, covariates, explainable)
    print(
        f"  the budget: abundance {steps[0]:.1%}, + density {steps[1] - steps[0]:+.1%}, "
        f"+ autofluorescence {steps[2] - steps[1]:+.1%}; left {1 - steps[2]:.1%}"
    )
    figure_covariates(table, explainable)
    return table


def replication_table(
    nano: np.ndarray,
    splits: list[tuple[list[int], list[int]]],
    map_agreement: list[float],
    left_agreement: list[float],
) -> pd.DataFrame:
    """replication.csv: per split, its two halves and the two agreements."""
    return pd.DataFrame(
        dict(
            half_a=[" ".join(ADULTS[i] for i in a) for a, _ in splits],
            half_b=[" ".join(ADULTS[i] for i in b) for _, b in splits],
            map_agreement=map_agreement,
            leftover_agreement=left_agreement,
        )
    )


def step3_replication(
    nano: np.ndarray,
    xs: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    map_agreement: list[float],
) -> tuple[list[float], float]:
    """Whether the leftover replicates across mice, and what that is worth."""
    print("\nSTEP 3  does the leftover replicate across mice?")
    agreement = leftover_agreement(nano, xs, splits)
    y = full_map(nano)
    h = float(np.mean(map_agreement))
    implied = implied_replication(h, r_squared(y, xs), len(xs) + 1, len(y))
    replication_table(nano, splits, map_agreement, agreement).to_csv(
        REPLICATION, index=False
    )
    print(
        f"  half against half: the map {h:.3f}, its leftover {np.mean(agreement):.3f} "
        f"(lowest split {min(agreement):.3f})"
    )
    print(
        f"  the map's reliability and the fit alone imply {implied:.3f}: a leftover "
        "replicating near that carries no extra evidence"
    )
    return agreement, implied


def same_sign_share(
    nano: np.ndarray,
    xs: list[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    res: np.ndarray,
) -> np.ndarray:
    """Per structure, the share of the half-cohort leftovers with the cohort's sign.

    Every split gives two half-cohorts, each fitted on its own; a structure whose
    leftover has the same sign in all of them is above (or below) prediction in
    every five adults.
    """
    same = np.zeros(len(res))
    for a, b in splits:
        for half in (a, b):
            same += np.sign(residual(half_map(nano, half), xs)) == np.sign(res)
    return same / (2 * len(splits))


def residual_table(
    inputs: Inputs,
    y: np.ndarray,
    predicted: np.ndarray,
    res: np.ndarray,
    same: np.ndarray,
    per_adult: np.ndarray,
) -> pd.DataFrame:
    """residual_by_structure.csv: per structure, largest leftover first.

    Beside the cohort's leftover: the share of half-cohort leftovers with its sign
    (same_sign), and the leftover of each adult's own map (`per_adult`, adults x
    structures) as a mean, an SD and a t across the adults (t_adults), the
    reliability that the bars of figure 12 show in grey.
    """
    n = per_adult.shape[0]
    mean = per_adult.mean(axis=0)
    sd = per_adult.std(axis=0, ddof=1)
    s = inputs.structures
    table = pd.DataFrame(
        dict(
            structure=s,
            acronym=[inputs.acronym.get(x, "") for x in s],
            division=[inputs.division.get(x, "") for x in s],
            nano_rank=y,
            predicted_rank=predicted,
            residual=res,
            same_sign=same,
            adult_mean=mean,
            adult_sd=sd,
            t_adults=mean / (sd / np.sqrt(n)),
        )
    )
    return table.sort_values("residual", ascending=False, ignore_index=True)


def per_adult_leftovers(nano: np.ndarray, xs: list[np.ndarray]) -> np.ndarray:
    """Each adult's own leftover: its ranked map fitted to the model, adults x
    structures."""
    return np.array([residual(rankdata(row), xs) for row in nano])


def leftover_boot(nano: np.ndarray, xs: list[np.ndarray], seed: int = 0) -> np.ndarray:
    """Leftovers of resampled cohorts, n_boot_mice x structures.

    Each is the mean of the adults drawn with replacement, ranked and fitted to the
    model's columns `xs`.
    """
    boot = gene_ranking.bootstrap_maps(nano, seed=seed)
    return np.array([residual(rankdata(b), xs) for b in boot])


def leftover_genes(
    inputs: Inputs,
    res: np.ndarray,
    xs: list[np.ndarray],
    psd: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray, list[str]]:
    """Every gene against the leftover, with the leftover's own spatial null.

    The leftover is a map of its own, with its own smoothness, so a gene's rho with
    it is tested as in analysis 1: against the leftover's surrogates (Burt 2020, on
    the structures' one-hemisphere centroids). Genes in the model correlate with
    the leftover near zero by construction (in_model says which); a gene that
    follows it predicts what receptor mRNA and synaptic density leave. The gene sets
    of analysis 3 are read the same way. None of these tests was named before the
    leftover was seen; they describe it.

    Returns the gene table, the set table, the surrogates, every gene's null rho
    (genes x surrogates) and the genes in the order of its rows.
    """
    s = inputs.structures
    centroids = load_centroids().loc[s]
    d = spatial_null.distance_matrix(centroids)
    surr = spatial_null.surrogates(res, d, seed=0)
    boot = leftover_boot(inputs.nano, xs)
    vectors = gene_ranking.gene_vectors(inputs.expr, s)
    table, null = gene_ranking.rank_genes(res, surr, boot, vectors)
    table = gene_ranking.with_q_and_ranks(table, inputs.p9_genes)
    in_model = {g: "abundance" for g in SUBUNITS}
    in_model.update({g: "psd_pc1" for g in psd})
    in_model.update({g: "markers" for g in MARKERS})
    table["in_model"] = table["symbol"].map(in_model).fillna("")
    table["p9_gene"] = table["symbol"].isin(inputs.p9_genes)
    genes = gene_table.per_gene(gene_table.load_gene_table()).set_index("symbol")
    table["gene_sets"] = table["symbol"].map(genes["gene_sets"]).fillna("")
    table["reliability"] = table["symbol"].map(genes["reliability"])

    # the gene sets, as analysis 3 reads them on the map
    members = gene_sets.members_from_table(genes.reset_index())
    tested = list(vectors)
    rho = table.set_index("symbol")["rho"]
    sets = gene_sets.set_tests(members, {"leftover": rho}, {"leftover": null}, tested)
    table = table.sort_values("rho", ascending=False, ignore_index=True)
    return table, sets, surr, null, tested


def figure_residual(
    map_agreement: list[float],
    agreement: list[float],
    implied: float,
    table: pd.DataFrame,
) -> None:
    """Draw fig3: the half-cohort agreements, and the structures most off prediction."""
    fig, axes = plt.subplots(
        1, 2, figsize=(11.8, 4.5), gridspec_kw=dict(width_ratios=[1, 1.15])
    )
    bins = np.linspace(
        min(min(map_agreement), min(agreement)) - 0.005,
        max(max(map_agreement), max(agreement)) + 0.005,
        30,
    )
    axes[0].hist(
        map_agreement,
        bins=bins,
        color="0.6",
        edgecolor="0.3",
        lw=0.3,
        alpha=0.85,
        label="the map itself",
    )
    axes[0].hist(
        agreement,
        bins=bins,
        color=RED,
        edgecolor="0.3",
        lw=0.3,
        alpha=0.7,
        label="what is left of it",
    )
    axes[0].axvline(implied, color="0.1", lw=1, ls=(0, (3, 2)))
    axes[0].set_xlabel("half-cohort against half-cohort (Spearman)", fontsize=8)
    axes[0].set_ylabel(f"splits of ten adults ({len(agreement)})", fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False, loc="upper left")
    axes[0].set_title(
        f"Step 3. map {np.mean(map_agreement):.3f}, leftover {np.mean(agreement):.3f}\n"
        f"dashed: {implied:.3f}, what the ceiling and the fit alone imply",
        fontsize=9,
    )
    tidy(axes[0])

    show = pd.concat([table.head(8), table.tail(6)])
    axes[1].barh(
        np.arange(len(show)),
        show["residual"],
        color=[RED if v > 0 else DARK_BLUE for v in show["residual"]],
        edgecolor="0.25",
        linewidth=0.4,
    )
    axes[1].set_yticks(np.arange(len(show)))
    axes[1].set_yticklabels([s[:36] for s in show["structure"]], fontsize=6.8)
    axes[1].invert_yaxis()
    axes[1].axvline(0, color="0.3", lw=0.7)
    axes[1].set_xlabel("nano rank minus predicted rank", fontsize=8)
    axes[1].set_title("Step 4. where it is largest", fontsize=9)
    tidy(axes[1])
    fig.tight_layout()
    save(fig, "fig3_residual.png")


def step4_where(
    inputs: Inputs,
    xs: list[np.ndarray],
    psd: list[str],
    splits: list[tuple[list[int], list[int]]],
    map_agreement: list[float],
    agreement: list[float],
    implied: float,
) -> None:
    """Which structures carry the leftover, how steadily, and which genes follow it."""
    print("\nSTEP 4  where the leftover lives")
    y = full_map(inputs.nano)
    res = residual(y, xs)
    same = same_sign_share(inputs.nano, xs, splits, res)
    per_adult = per_adult_leftovers(inputs.nano, xs)
    table = residual_table(inputs, y, y - res, res, same, per_adult)
    table.to_csv(RESIDUALS, index=False)
    print(
        "  above prediction (ranks; share of half-cohorts with the same sign; t across"
        " the adults' own leftovers):"
    )
    for r in table.head(8).itertuples():
        print(
            f"    {r.structure[:48]:50s} {r.residual:+6.1f}   {r.same_sign:.2f}"
            f"   t {r.t_adults:+6.1f}"
        )
    print("  below:")
    for r in table.tail(6).itertuples():
        print(
            f"    {r.structure[:48]:50s} {r.residual:+6.1f}   {r.same_sign:.2f}"
            f"   t {r.t_adults:+6.1f}"
        )

    genes, sets, surr, null, tested = leftover_genes(inputs, res, xs, psd)
    genes.to_csv(LEFTOVER_GENES, index=False)
    sets.to_csv(LEFTOVER_SETS, index=False)
    np.savez(
        LEFTOVER_NULL,
        structures=np.array(inputs.structures),
        surrogates=surr.astype(np.float32),
        genes=np.array(tested),
        null_rho=null,
    )
    print(
        f"  every gene against the leftover ({len(genes)} genes, {len(surr)} "
        "surrogates of the leftover):"
    )
    for r in genes.head(8).itertuples():
        print(
            f"    {r.symbol:10s} rho {r.rho:+.3f}  spatial p {r.p_spatial:.4f}  "
            f"q {r.q_all:.3f}  {r.in_model}"
        )
    for g in ("Cacng8", "Gria1"):
        r = genes.set_index("symbol").loc[g]
        print(
            f"    {g}: rho {r['rho']:+.3f}, rank {int(r['rank_all'])} of {len(genes)}, "
            f"spatial p {r['p_spatial']:.4f}"
        )
    for r in sets.itertuples():
        p = f"p {r.p_spatial:.3f}" if r.tested else "not tested"
        print(
            f"    set {r.gene_set:20s} {r.n_genes:3d} genes, "
            f"median {r.median_rho:+.3f}, {p}"
        )
    figure_residual(map_agreement, agreement, implied, table)


def main() -> None:
    """Run the five steps on the ten adults and write their tables and figures."""
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs()
    step0_structures(inputs)
    print(f"  {len(inputs.structures)} structures, {len(ADULTS)} adults")

    # the ceiling, the budget, the leftover's replication, and where it lives
    splits = half_splits()
    explainable, map_agreement = step1_ceiling(inputs.nano, splits)
    covariates, psd, share = build_covariates(
        inputs.expr, inputs.role, inputs.auto, inputs.structures
    )
    y = full_map(inputs.nano)
    step2_budget(y, covariates, psd, share, explainable)
    xs = model(covariates)
    agreement, implied = step3_replication(inputs.nano, xs, splits, map_agreement)
    step4_where(inputs, xs, psd, splits, map_agreement, agreement, implied)
    print("\nNow run run_beyond_controls.py: it tries to break this seven ways.")
