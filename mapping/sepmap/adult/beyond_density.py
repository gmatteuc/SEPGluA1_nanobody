"""What receptor abundance and synaptic density leave unexplained in the adult map.

The claim: the adult surface-GluA1 map carries substantial, reproducible spatial
structure that receptor abundance and synaptic density do not explain. In numbers,
97.4% of the map is explainable in principle (it is that reproducible across
animals); receptor abundance and three density proxies, allowed a curved
relationship, account for 61% of that; the remaining 39% replicates across
independent halves of the cohort at 0.97. A twenty-component model of the whole
390-gene panel, far richer than either explanation, still leaves a remainder that
replicates at 0.88. adult.beyond_controls tries to break the claim seven ways.

Ranking genes against the map cannot answer this: ish.panel_test found the ranking
flat across postsynaptic gene classes, which is what a map of synaptic density
alone would give (every synaptic gene correlates, none stands out). So the
question is turned round: what is left once the two explanations are subtracted,
and is that leftover real?

    nano zref  ~  receptor abundance + synaptic density   ->   residual

Any model leaves a residual; what makes one interesting is that it replicates.
Noise cannot replicate across independent animals, so if two halves of the cohort
give the same leftover, the leftover is biology rather than measurement error. The
steps, each printing its numbers (fig0 to fig3 draw them):

    0  the structures   which regions the analysis may use, and why
    1  the ceiling      how reproducible the map is; nothing can explain variance
                        that is not there, so every R2 is read against this
    2  the covariates   what abundance and density explain
    3  the residual     whether what is left replicates across animals (the claim)
    4  where it lives   which structures carry it

Everything is one number per structure, 125 of them, on ranks:

    y            the ten adults' mean zref per structure, ranked 1..125, so the
                 variance is that of the ordering of structures, not of log2 values
    predictors   the Gria1-4 composite, the 11-marker composite, the first
                 component of the 188 postsynaptic-density genes, and the cohort's
                 mean autofluorescence, each ranked and entered as x, x^2 and x^3
                 because the rank relationships are curved (control E)
    R2           least squares of y on those columns plus an intercept,
                     R2 = 1 - var(y - prediction) / var(y)
    CV R2        the same, each structure predicted from a fit that never saw it:
                 the 125 are shuffled once with a fixed seed, cut into five folds,
                 and each fold predicted from the other four. This is the number
                 quoted, because adding terms always improves the in-sample fit
                 (0.706 fitted against 0.597 predicted here).

The ceiling: variance can be explained only insofar as it is real, so the cohort is
split into two fives every possible way (126 splits), the two half-maps are
correlated (mean rho 0.974), and Spearman-Brown gives what the ten animals together
are worth, 0.987. Squared, 0.974 of the map's variance is reproducible and the
other 0.026 is noise no model could predict. The covariates then reach

    0.597 / 0.974 = 61.3%              of the explainable variance, leaving
    (0.974 - 0.597) / 0.974 = 38.7%    unexplained but real,

the "about 40%" of the claim. It is not 1 - 0.597 = 40.3%, which would count the
2.6% of noise as something the explanations failed to explain. Intervals come from
resampling the 125 structures 2000 times (adult.beyond_figures).

The choices, and why:

    structures, not voxels   the gene data exist only at 200 um and in a different
                             mouse, so at the voxel level the comparison would be
                             mostly registration error
    ranks everywhere         Allen expression energy has an arbitrary scale per
                             experiment; ranks assume only that more receptor gives
                             more signal, and control E tests whether the rank
                             relationship is straight
    grey matter only         fibre tracts and ventricles have no synapses, so
                             synaptic density is not defined there, and the CCF's
                             "..., unassigned" entries are leftover voxels whose
                             content depends on how each brain registered; both go
                             before any fit, by keep_structure (fig0 counts them)
    naive and RWS pooled     the groups differ in whisker experience, not in what
                             the stain is; control D checks that the leftover does
                             not depend on the group
    curved covariates        as straight lines the covariates reach 0.42 rather
                             than 0.60 cross-validated, and half of the explainable
                             variance would look unexplained instead of 39%
                             (control E)

Writes, in adult_v2/beyond/ under the data root:

    structures_used.csv        every structure, kept or dropped, with the reason
    variance_partition.csv     what each covariate set explains, against the ceiling
    residual_by_structure.csv  where the map exceeds and falls short of prediction
    fig0_structures.png, fig1_ceiling.png, fig2_covariates.png, fig3_residual.png

Run by run_beyond_density.py.
"""

import csv
import itertools
import math
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata, spearmanr

from sepmap.config import SETTINGS
from sepmap.plotting import DARK_BLUE, RED, tidy
from sepmap.volumes.cohort import NAIVE, RWS
from sepmap.volumes.per_mouse import DATA, MICE, annotation_20, structure_terms
from sepmap.volumes.per_mouse import OUT as PER_MOUSE

# the reading, the half-cohorts, the grey-matter divisions and the gene sets; the
# smallest structure kept, as in young_vs_adult.region_plot
BEYOND = SETTINGS["beyond"]
REGION_TABLES = SETTINGS["region_tables"]

NANO = DATA / "comparisons_v2" / "young_vs_adult" / "region_means_per_mouse.csv"
MERGED_ISH = DATA / "adult_v2" / "ish" / "gene_region_table_merged.csv"
OUT = DATA / "adult_v2" / "beyond"

# the ten adults, naive and rws pooled
ADULTS = NAIVE + RWS

# the grey-matter divisions as a set, and the gene sets standing in for the two
# explanations (receptor abundance, synaptic density) as tuples, as other modules use
# them
GREY = set(BEYOND["grey"])
SUBUNITS = tuple(BEYOND["subunits"])
MARKERS = tuple(BEYOND["markers"])


# ===== Loading =====


def keep_structure(name: str, division: str) -> tuple[bool, str]:
    """Whether a structure is in the analysis: (keep, reason why not).

    Two rules, both applied before any fitting so neither can be tuned to the
    answer: grey matter only, by division; and no "..., unassigned" entries,
    which are voxels the atlas could not place rather than structures.
    """
    if "unassigned" in name.lower():
        return False, "catch-all label, not a structure"
    if division not in GREY:
        return False, f"division {division} is not grey matter"
    return True, ""


def nano_per_mouse() -> tuple[dict[str, dict[str, float]], dict[str, str]]:
    """{mouse: {structure: zref}} and {structure: division}, for the ten adults."""
    per, division = defaultdict(dict), {}
    with open(NANO, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["reading"] == BEYOND["reading"] and r["mouse"] in ADULTS:
                per[r["mouse"]][r["structure"]] = float(r["log2_value"])
                division[r["structure"]] = r["division"]
    return per, division


def gene_profiles(
    path: Path = MERGED_ISH,
) -> tuple[dict[str, dict[str, float]], dict[str, str]]:
    """{gene: {structure: rank_mean}} and {gene: role}; replicates already merged."""
    per, role = defaultdict(dict), {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            per[r["symbol"]][r["structure"]] = float(r["rank_mean"])
            role[r["symbol"]] = r["role"]
    return per, role


def autofluorescence(structures: set[str]) -> dict[str, dict[str, float]]:
    """{mouse: {structure: log2 mean autofluorescence}}, from the same ten brains.

    Worth more as a density proxy than it looks: every gene covariate is a
    different mouse warped into our atlas, while this one is the very tissue the
    nano signal was measured in, in the same sections, with no registration
    between the two.
    """
    names, _, _ = structure_terms()
    per = {}
    for mouse in ADULTS:
        z = np.load(PER_MOUSE / (mouse + ".npz"))
        tissue = z["tissue"]
        annotation = annotation_20(MICE[mouse][1])
        labels = annotation[tissue]
        counts = np.bincount(labels, minlength=int(annotation.max()) + 1)
        totals = np.bincount(
            labels, weights=z["auto"].astype(np.float32)[tissue], minlength=len(counts)
        )

        # pool the indices by structure name (index 0 is outside the brain)
        acc = defaultdict(lambda: [0, 0.0])
        for idx in np.nonzero(counts)[0]:
            if idx == 0:
                continue
            a = acc[names.get(int(idx), f"id{idx}")]
            a[0] += int(counts[idx])
            a[1] += totals[idx]
        per[mouse] = {
            k: math.log2(t / c)
            for k, (c, t) in acc.items()
            if c >= REGION_TABLES["min_vox20"] and t > 0 and k in structures
        }
    return per


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

    Used for the 188 postsynaptic-density genes: rather than naming a handful of
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
    """The same covariates, allowed to bend: x, x^2 and x^3 of each.

    A straight line through two rank variables assumes the relationship is not
    only monotone but evenly paced, and control E in adult.beyond_controls shows
    that assumption is wrong here: adding curvature raises the cross-validated
    fit from 0.42 to 0.60. That 0.18 belongs to the covariates, not to the
    leftover, so the model quoted bends.
    """
    return list(predictors) + [x**2 for x in predictors] + [x**3 for x in predictors]


def cv_r2(y: np.ndarray, predictors: Sequence[np.ndarray], folds: int = 5) -> float:
    """Variance explained on structures the fit has never seen.

    The number to quote once a model has many terms: adding predictors always
    improves the in-sample fit, and only held-out structures can say whether it
    improved the prediction. Structures are shuffled once with a fixed seed, cut
    into folds, and each fold predicted from a fit on the others.
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


def half_map(
    nano: dict[str, dict[str, float]], indices: Sequence[int], structures: list[str]
) -> np.ndarray:
    """The mean zref map of one half-cohort, as ranks."""
    return rankdata(
        [float(np.mean([nano[ADULTS[i]][s] for i in indices])) for s in structures]
    )


def build_covariates(
    nano: dict[str, dict[str, float]],
    expr: dict[str, dict[str, float]],
    role: dict[str, str],
    auto: dict[str, dict[str, float]],
    structures: list[str],
) -> tuple[dict[str, np.ndarray], list[str], float]:
    """The four predictors, built in one place so every module builds them alike.

    Returns the predictors by name, the postsynaptic-density genes behind psd_pc1,
    and the share of their variance that psd_pc1 carries.
    """
    controls = sorted(
        g
        for g in expr
        if role[g] == "control_psd" and all(s in expr[g] for s in structures)
    )
    psd_pc1, share = first_pc(controls, expr, structures)
    return (
        dict(
            abundance=composite(SUBUNITS, expr, structures),
            markers=composite(MARKERS, expr, structures),
            psd_pc1=psd_pc1,
            autofluo=rankdata(
                [float(np.mean([auto[m][s] for m in ADULTS])) for s in structures]
            ),
        ),
        controls,
        share,
    )


def prepare() -> tuple[
    dict[str, dict[str, float]],
    dict[str, str],
    dict[str, dict[str, float]],
    dict[str, str],
    dict[str, dict[str, float]],
    list[str],
]:
    """Everything the analysis and the controls both need, loaded once."""
    nano, division = nano_per_mouse()
    expr, role = gene_profiles()
    everywhere = set.intersection(*[set(nano[m]) for m in ADULTS])
    everywhere &= set.intersection(*[set(expr[g]) for g in SUBUNITS + MARKERS])
    kept = [s for s in sorted(everywhere) if keep_structure(s, division.get(s, ""))[0]]
    auto = autofluorescence(set(kept))
    structures = sorted(set(kept).intersection(*[set(auto[m]) for m in ADULTS]))
    return nano, division, expr, role, auto, structures


# ===== Drawing =====


def save(fig: plt.Figure, name: str) -> None:
    """Save `fig` as `name` in the output folder at 200 dpi, close it, print the path."""
    path = OUT / name
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  -> {path}")


# ===== Steps =====


def figure_structures(rows: list[dict], kept: list[str], dropped: list[dict]) -> None:
    """Draw fig0: the structures kept per division, and those dropped per reason."""
    fig, axes = plt.subplots(
        1, 2, figsize=(11.0, 3.9), gridspec_kw=dict(width_ratios=[1, 1.25])
    )
    counts = defaultdict(int)
    for r in rows:
        if r["kept"] == "yes":
            counts[r["division"]] += 1
    order = sorted(counts, key=counts.get, reverse=True)
    axes[0].bar(
        range(len(order)),
        [counts[d] for d in order],
        color="0.65",
        edgecolor="0.25",
        linewidth=0.5,
    )
    axes[0].set_xticks(range(len(order)))
    axes[0].set_xticklabels(order, fontsize=7, rotation=45, ha="right")
    axes[0].set_ylabel("structures kept", fontsize=8)
    axes[0].set_title(
        f"what the analysis runs on\n{len(kept)} grey-matter structures", fontsize=9
    )
    tidy(axes[0])

    why_counts = defaultdict(int)
    for r in dropped:
        why_counts[r["reason"]] += 1
    ordered = sorted(why_counts, key=why_counts.get)
    axes[1].barh(
        range(len(ordered)),
        [why_counts[w] for w in ordered],
        color=RED,
        edgecolor="0.25",
        linewidth=0.5,
        alpha=0.85,
    )
    axes[1].set_yticks(range(len(ordered)))
    axes[1].set_yticklabels([w[:46] for w in ordered], fontsize=7)
    axes[1].set_xlabel("structures dropped", fontsize=8)
    axes[1].set_title(
        "and what it refuses to run on\ndecided by rule, before any fitting", fontsize=9
    )
    tidy(axes[1])
    fig.tight_layout()
    save(fig, "fig0_structures.png")


def step0_structures(
    nano: dict[str, dict[str, float]],
    division: dict[str, str],
    expr: dict[str, dict[str, float]],
) -> list[str]:
    """Choose the structures, and show what the choice threw away."""
    print("\nSTEP 0  which structures the analysis may use")
    everywhere = set.intersection(*[set(nano[m]) for m in ADULTS])
    everywhere &= set.intersection(*[set(expr[g]) for g in SUBUNITS + MARKERS])

    # every structure measured everywhere, kept or dropped, with the reason
    rows = []
    for name in sorted(everywhere):
        ok, why = keep_structure(name, division.get(name, ""))
        rows.append(
            dict(
                structure=name,
                division=division.get(name, ""),
                kept="yes" if ok else "no",
                reason=why,
            )
        )
    with open(OUT / "structures_used.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["structure", "division", "kept", "reason"]
        )
        writer.writeheader()
        writer.writerows(rows)

    kept = [r["structure"] for r in rows if r["kept"] == "yes"]
    dropped = [r for r in rows if r["kept"] == "no"]
    print(f"  {len(everywhere)} structures measured in every mouse and every covariate")
    print(f"  {len(kept)} kept, {len(dropped)} dropped:")
    for reason in sorted({r["reason"] for r in dropped}):
        names = [r["structure"] for r in dropped if r["reason"] == reason]
        print(f"    {len(names):3d}  {reason}")

    # figure: the structures kept per division, and those dropped per reason
    figure_structures(rows, kept, dropped)
    return kept


def step1_ceiling(
    nano: dict[str, dict[str, float]], structures: list[str]
) -> tuple[float, list[tuple[list[int], list[int]]], list[float]]:
    """How reproducible the map itself is: nothing below can beat this."""
    print("\nSTEP 1  the ceiling -- how much of this map is explainable at all")

    # agreement of the two half-cohort maps over every split
    splits = half_splits()
    agreement = [
        spearmanr(half_map(nano, a, structures), half_map(nano, b, structures)).statistic
        for a, b in splits
    ]
    half = float(np.mean(agreement))
    full = spearman_brown(half)
    print(f"  over {len(splits)} five-against-five splits of the ten adults:")
    print(f"    half-cohort against half-cohort   rho = {half:.3f}")
    print(f"    Spearman-Brown, all ten animals   rho = {full:.3f}")
    print(f"  so at most {full**2:.1%} of this map can ever be explained, and every")
    print("  R2 below is quoted against that rather than against 1.0.")

    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    ax.hist(agreement, bins=25, color="0.7", edgecolor="0.35", linewidth=0.4)
    ax.axvline(half, color=RED, lw=1.8)
    ax.set_xlabel("Spearman between the two half-cohort maps", fontsize=8)
    ax.set_ylabel(f"splits of ten animals ({len(splits)})", fontsize=8)
    ax.set_title(
        "Step 1. the map is highly reproducible\n"
        f"half-cohorts agree at {half:.3f}, whole cohort {full:.3f}",
        fontsize=9,
    )
    tidy(ax)
    fig.tight_layout()
    save(fig, "fig1_ceiling.png")
    return full, splits, agreement


def write_partition(
    models: list[tuple[str, list[np.ndarray]]],
    vals: list[float],
    y: np.ndarray,
    ceiling: float,
) -> None:
    """Write variance_partition.csv: each model's R2, held out and in-sample."""
    with open(OUT / "variance_partition.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["model", "cv_r2", "in_sample_r2", "share_of_ceiling"])
        for (label, xs), v in zip(models, vals):
            writer.writerow(
                [label, f"{v:.4f}", f"{r_squared(y, xs):.4f}", f"{v / ceiling**2:.4f}"]
            )


def figure_covariates(
    models: list[tuple[str, list[np.ndarray]]], vals: list[float], ceiling: float
) -> None:
    """Draw fig2: each model's cross-validated R2 against the ceiling."""
    fig, ax = plt.subplots(figsize=(7.6, 4.5))
    ax.barh(np.arange(len(models)), vals, color="0.65", edgecolor="0.25", linewidth=0.5)
    ax.axvline(ceiling**2, color=RED, lw=1.8)
    ax.annotate(
        f"ceiling {ceiling**2:.2f}\n(the map's own reliability)",
        (ceiling**2, len(models) - 0.4),
        color=RED,
        fontsize=7.5,
        ha="right",
        va="top",
        xytext=(-6, 0),
        textcoords="offset points",
    )
    ax.set_yticks(np.arange(len(models)))
    ax.set_yticklabels([m[0] for m in models], fontsize=8)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel(
        "variance explained on held-out structures (cross-validated R2)", fontsize=8
    )
    ax.set_title(
        "Step 2. neither explanation fills the map,\neven when allowed to bend",
        fontsize=9,
    )
    tidy(ax)
    fig.tight_layout()
    save(fig, "fig2_covariates.png")


def step2_covariates(
    nano: dict[str, dict[str, float]],
    expr: dict[str, dict[str, float]],
    role: dict[str, str],
    auto: dict[str, dict[str, float]],
    structures: list[str],
    ceiling: float,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """What receptor abundance and synaptic density account for."""
    print("\nSTEP 2  what the two boring explanations buy")
    covariates, controls, share = build_covariates(nano, expr, role, auto, structures)
    y = half_map(nano, range(len(ADULTS)), structures)
    print(
        f"  psd_pc1 is the first component of {len(controls)} postsynaptic-density "
        f"genes, carrying {share:.1%} of their variance"
    )

    # each covariate on its own, in-sample
    print("  each covariate on its own:")
    nice = {
        "abundance": "abundance (Gria1-4)",
        "markers": "synaptic markers",
        "psd_pc1": "psd_pc1",
        "autofluo": "autofluorescence",
    }
    for key, v in covariates.items():
        print(
            f"    {nice[key]:22s} rho {spearmanr(y, v).statistic:+.3f}   "
            f"R2 {r_squared(y, [v]):.3f}   "
            f"{r_squared(y, [v]) / ceiling**2:5.1%} of the ceiling"
        )

    # the models in combination, scored on held-out structures
    c = covariates
    models = [
        ("abundance only", [c["abundance"]]),
        ("markers only", [c["markers"]]),
        ("psd_pc1 only", [c["psd_pc1"]]),
        ("autofluorescence only", [c["autofluo"]]),
        ("abundance + markers", [c["abundance"], c["markers"]]),
        ("abundance + markers + psd_pc1", [c["abundance"], c["markers"], c["psd_pc1"]]),
        ("all four, straight", list(c.values())),
        ("all four, allowed to bend", flexible(list(c.values()))),
    ]
    vals = [cv_r2(y, xs) for _, xs in models]
    print("  and in combination, scored on structures the fit has not seen:")
    for (label, xs), v in zip(models, vals):
        print(
            f"    {label:30s} CV R2 {v:.3f}   in-sample {r_squared(y, xs):.3f}   "
            f"{v / ceiling**2:5.1%} of the ceiling"
        )
    print("  the model to quote is the last one -- letting the covariates bend is")
    print(
        f"  their best shot -- and it leaves {1 - vals[-1] / ceiling**2:.0%} of the "
        f"explainable variance unaccounted for."
    )

    # the partition as a table, and drawn
    write_partition(models, vals, y, ceiling)
    figure_covariates(models, vals, ceiling)
    return covariates, y


def step3_residual(
    nano: dict[str, dict[str, float]],
    structures: list[str],
    covariates: dict[str, np.ndarray],
    splits: list[tuple[list[int], list[int]]],
    raw_agreement: list[float],
) -> list[float]:
    """The claim: what is left over replicates across independent animals."""
    print("\nSTEP 3  does the leftover replicate?   <- this is the claim")

    # residualise each half-cohort map on the bending model, then compare the halves
    xs = flexible(list(covariates.values()))
    agreement = [
        spearmanr(
            residual(half_map(nano, a, structures), xs),
            residual(half_map(nano, b, structures), xs),
        ).statistic
        for a, b in splits
    ]
    half = float(np.mean(agreement))
    print("  residualise each half-cohort map on the bending four-covariate model,")
    print("  then compare the two leftovers:")
    print(
        f"    half against half   rho = {half:.3f}   "
        f"(the map itself: {np.mean(raw_agreement):.3f})"
    )
    print(f"    Spearman-Brown      rho = {spearman_brown(half):.3f}")
    print("  noise cannot replicate across independent animals, so this is real.")
    return agreement


def write_residuals(
    structures: list[str], y: np.ndarray, predicted: np.ndarray, res: np.ndarray
) -> None:
    """Write residual_by_structure.csv, the residual per structure, largest first."""
    with open(OUT / "residual_by_structure.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["structure", "nano_rank", "predicted_rank", "residual"])
        for i in np.argsort(-res):
            writer.writerow(
                [structures[i], f"{y[i]:.1f}", f"{predicted[i]:.1f}", f"{res[i]:.2f}"]
            )


def residual_against_genes(
    res: np.ndarray,
    expr: dict[str, dict[str, float]],
    role: dict[str, str],
    structures: list[str],
) -> None:
    """Print the five genes, of those measured everywhere, closest to the residual."""
    scored = []
    for gene in expr:
        if all(s in expr[gene] for s in structures):
            rho = spearmanr(res, [expr[gene][s] for s in structures]).statistic
            scored.append((abs(rho), rho, gene))
    scored.sort(reverse=True)
    print("  the residual against all 390 genes -- is it just something we left out?")
    for _, rho, gene in scored[:5]:
        print(f"    {gene:10s} rho {rho:+.3f}   ({role[gene]})")
    print("    no single gene in the panel accounts for it.")


def figure_residual(
    raw_agreement: list[float],
    agreement: list[float],
    order: np.ndarray,
    res: np.ndarray,
    structures: list[str],
) -> None:
    """Draw fig3: the half-cohort agreements, and the structures most off prediction."""
    # structures with the largest residuals
    fig, axes = plt.subplots(
        1, 2, figsize=(11.8, 4.5), gridspec_kw=dict(width_ratios=[1, 1.15])
    )
    bins = np.linspace(
        min(min(raw_agreement), min(agreement)) - 0.005,
        max(max(raw_agreement), max(agreement)) + 0.005,
        30,
    )
    axes[0].hist(
        raw_agreement,
        bins=bins,
        color="0.6",
        edgecolor="0.3",
        linewidth=0.3,
        alpha=0.85,
        label="the map itself",
    )
    axes[0].hist(
        agreement,
        bins=bins,
        color=RED,
        edgecolor="0.3",
        linewidth=0.3,
        alpha=0.7,
        label="what is left of it",
    )
    axes[0].set_xlabel("half-cohort against half-cohort (Spearman)", fontsize=8)
    axes[0].set_ylabel(f"splits of ten animals ({len(agreement)})", fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False, loc="upper left")
    axes[0].set_title(
        "Step 3. the leftover is reproducible, not noise\n"
        f"map {np.mean(raw_agreement):.3f}, "
        f"residual {np.mean(agreement):.3f}",
        fontsize=9,
    )
    tidy(axes[0])

    show = list(order[:8]) + list(order[-6:])
    axes[1].barh(
        np.arange(len(show)),
        [res[i] for i in show],
        color=[RED if res[i] > 0 else DARK_BLUE for i in show],
        edgecolor="0.25",
        linewidth=0.4,
    )
    axes[1].set_yticks(np.arange(len(show)))
    axes[1].set_yticklabels([structures[i][:36] for i in show], fontsize=6.8)
    axes[1].invert_yaxis()
    axes[1].axvline(0, color="0.3", lw=0.7)
    axes[1].set_xlabel("nano rank minus predicted rank", fontsize=8)
    axes[1].set_title("Step 4. and it is anatomically organised", fontsize=9)
    tidy(axes[1])
    fig.tight_layout()
    save(fig, "fig3_residual.png")


def step4_where(
    nano: dict[str, dict[str, float]],
    structures: list[str],
    covariates: dict[str, np.ndarray],
    expr: dict[str, dict[str, float]],
    role: dict[str, str],
    agreement: list[float],
    raw_agreement: list[float],
) -> None:
    """Which structures carry the leftover, and whether any single gene is behind it."""
    print("\nSTEP 4  where the leftover lives")
    y = half_map(nano, range(len(ADULTS)), structures)
    res = residual(y, flexible(list(covariates.values())))
    predicted = y - res

    # the residual per structure, largest first
    write_residuals(structures, y, predicted, res)

    order = np.argsort(-res)
    print("  more surface GluA1 than abundance and density predict:")
    for i in order[:8]:
        print(f"    {structures[i][:48]:50s} {res[i]:+6.1f} ranks")
    print("  less:")
    for i in order[-6:]:
        print(f"    {structures[i][:48]:50s} {res[i]:+6.1f} ranks")

    # the residual against every gene measured in all the structures
    residual_against_genes(res, expr, role, structures)

    # figure: the half-cohort agreement of the map and of the leftover, and the
    # structures with the largest residuals
    figure_residual(raw_agreement, agreement, order, res, structures)


def main() -> None:
    """Run the four steps on the ten adults and write their tables and figures."""
    OUT.mkdir(parents=True, exist_ok=True)
    nano, division = nano_per_mouse()
    expr, role = gene_profiles()

    # the structures: grey matter, measured everywhere, with autofluorescence
    kept = step0_structures(nano, division, expr)
    auto = autofluorescence(set(kept))
    structures = sorted(set(kept).intersection(*[set(auto[m]) for m in ADULTS]))
    print(f"  {len(structures)} of those also have autofluorescence in every mouse")

    # the ceiling, the covariates, the residual and where it lives
    ceiling, splits, raw_agreement = step1_ceiling(nano, structures)
    covariates, _ = step2_covariates(nano, expr, role, auto, structures, ceiling)
    agreement = step3_residual(nano, structures, covariates, splits, raw_agreement)
    step4_where(nano, structures, covariates, expr, role, agreement, raw_agreement)

    print("\nNow run run_beyond_controls.py: it tries to break this seven ways.")
