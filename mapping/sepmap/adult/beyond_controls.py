"""Seven attempts to break the result of adult.beyond_density, and its check rows.

That module finds part of the adult map's reproducible pattern not predicted by
Gria1 expression or synapse density, with a leftover that replicates across mice.
Each control below is a specific way that could be an artefact, with the number that
says whether it is:

    A  a smooth spatial gradient?     an illumination artefact would look like
                                      one; anatomy would not. Positions are the
                                      structures' centroids in one hemisphere
                                      (sepmap.structures); averaged over both, every
                                      bilateral structure sat on the midline
    B  small or poorly covered        noisy means look like signal
       structures?
    C  one or two animals?            a single odd brain can carry a cohort mean
                                      of ten
    D  the whisker manipulation?      naive and RWS are pooled; if the leftover is
                                      the manipulation, the pooling hides it
    E  curvature the model misses?    a bent rank relationship fitted as a straight
                                      line leaves real structure in the residual:
                                      the held-out R2 of the straight model, of the
                                      model curved (squares and cubes, the check row
                                      curved) and bent to fifth powers. When curving
                                      buys more than beyond_controls.curvature_gain,
                                      the curved row's leftover, against its own
                                      floor, is the bound to read beside the main one
    F  the choice of predictors?      the strongest version: give the model the
                                      expression of every gene of the gene table
                                      measured in all the structures, as principal
                                      components (as many as cross-validation
                                      picks, inside each training fold), and ask
                                      whether what it leaves stands above what the
                                      same model leaves of a map made of those
                                      genes' expression (its own calibration
                                      floor). It is why the leftover is "not
                                      predicted by Gria1 expression or synapse
                                      density" and never "beyond gene expression"
    G  zref?                          the same test on the other readings of the
                                      map, and on zref with the 17-brain reference

E and F are cross-validated, because a flexible model always fits better on the
data it was fitted to; the question is whether it predicts better, and only
held-out structures can say. F picks its number of components inside each training
fold (nested), since the best of many held-out scores is itself optimistic.
Everything is imported from adult.beyond_density, so the two modules use the same
structures, predictors and arithmetic.

Beside the controls, the main model under other folds: one shuffling, single
shufflings, ten folds, leave one out and folds of spatial blocks, each with the
shares of Gria1 alone, density alone and both, and the share left. The check rows,
which change the model, are adult.beyond_checks'.

Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    controls.csv            one row per control, with its number and verdict
    gene_space.csv          control F: per number of components, R2 fitted and held out
    gene_space_calibration.csv  control F's model on the nano map and on maps made of
                            the genes' expression, each half of the experiments
    gene_space_summary.csv  control F's numbers: genes, components picked, nested
                            share of the ceiling, replication
    readings.csv            control G: per reading, R2 and the two replications
    folds.csv               the main model under its own folds and others: the
                            shares and the share left

Run by run_beyond_controls.py.
"""

from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from sepmap import structures
from sepmap.adult import beyond_calibration, profiles
from sepmap.adult.beyond_density import (
    BEYOND,
    OUT,
    Inputs,
    block_labels,
    ceiling,
    covariates_for,
    cv_r2,
    flexible,
    fold_labels,
    full_map,
    gene_matrix,
    half_splits,
    leftover_agreement,
    load_inputs,
    model_columns,
    model_genes,
    model_shares,
    per_adult_leftovers,
    predictors,
    r_squared,
    replicates,
    residual,
)
from sepmap.adult.profiles import ADULTS
from sepmap.config import SETTINGS
from sepmap.volumes.cohort import NAIVE, RWS
from sepmap.young_vs_adult.region_plot import REGION_MEANS

# the largest gene-space model of control F, and the threshold of each verdict
BEYOND_CONTROLS = SETTINGS["beyond_controls"]

CONTROLS = OUT / "controls.csv"
GENE_SPACE = OUT / "gene_space.csv"
GENE_SPACE_CALIBRATION = OUT / "gene_space_calibration.csv"
GENE_SPACE_SUMMARY = OUT / "gene_space_summary.csv"
READINGS = OUT / "readings.csv"
FOLDS = OUT / "folds.csv"

# the rows of each group in the adult matrices, naive then rws
NAIVE_ROWS = [ADULTS.index(m) for m in NAIVE]
RWS_ROWS = [ADULTS.index(m) for m in RWS]


# ===== Utilities =====


def replication(
    matrix: np.ndarray,
    columns: Sequence[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
) -> float:
    """How well the leftover of one half-cohort matches the other's, over every split."""
    return float(np.mean(leftover_agreement(matrix, columns, splits)))


def adult_pair_agreements(per_adult: np.ndarray) -> list[float]:
    """Spearman of every pair of adults' own leftovers (adults x structures)."""
    n = len(per_adult)
    return [
        float(spearmanr(per_adult[i], per_adult[j]).statistic)
        for i in range(n)
        for j in range(i + 1, n)
    ]


def group_leftovers(
    nano: np.ndarray, columns: list[np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    """The leftovers of the naive and of the RWS group's ranked mean maps."""
    naive = residual(rankdata(nano[NAIVE_ROWS].mean(axis=0)), columns)
    rws = residual(rankdata(nano[RWS_ROWS].mean(axis=0)), columns)
    return naive, rws


def mean_sizes(subset: list[str]) -> np.ndarray:
    """Each structure's mean volume over the adults, in 20 um voxels."""
    table = profiles.load_per_mouse()
    sizes = table[table["mouse"].isin(ADULTS)].groupby("structure")["n_vox20"].mean()
    return sizes.reindex(subset).to_numpy(float)


# ===== Controls =====


def control_a_space(
    leftover: np.ndarray,
    nano: np.ndarray,
    columns: list[np.ndarray],
    xyz: np.ndarray,
    splits: list[tuple[list[int], list[int]]],
) -> dict[str, str]:
    """Control A: whether the leftover is just a smooth gradient across the block.

    `xyz` holds each structure's centroid (AP, DV, ML) in mm. Returns the verdict
    row.
    """
    print("\nA  is it a smooth spatial gradient? (an illumination artefact)")

    # the residual against a quadratic in the three ranked positions, enough to
    # follow a ramp or a bowl across the block, as uneven illumination would give
    pos = [rankdata(xyz[:, i]) for i in range(3)]
    quad = pos + [p**2 for p in pos]
    smooth = r_squared(leftover, quad)
    for axis, p in zip("AP DV ML".split(), pos):
        print(f"   residual against {axis}: rho {spearmanr(leftover, p).statistic:+.3f}")
    print(f"   a smooth quadratic in all three axes explains R2 = {smooth:.3f} of it")

    # the replication with position added to the model, the ranked axes straight
    with_pos = replication(nano, columns + pos, splits)
    still = "still replicates" if replicates(with_pos) else "replicates only"
    print(f"   with position added to the model the leftover {still} at {with_pos:.3f}")

    # a gradient could explain the leftover if position explains more than
    # beyond_controls.gradient_r2 of it
    passed = smooth <= BEYOND_CONTROLS["gradient_r2"]
    if passed:
        verdict = "not a gradient; position is a weak predictor of it"
    else:
        verdict = "a gradient could explain it: look closer"
    print("   verdict: " + verdict)
    return dict(
        control="A spatial gradient",
        number=f"smooth R2 {smooth:.3f}, replication with position {with_pos:.3f}",
        verdict="pass" if passed else "does not pass",
    )


def control_b_size(leftover: np.ndarray, sizes: np.ndarray) -> dict[str, str]:
    """Control B: whether the leftover comes from small or poorly covered structures.

    `sizes` holds each structure's mean volume over the adults in 20 um voxels.
    Returns the verdict row.
    """
    print("\nB  is it small structures, where a mean is noisy?")
    ok = np.isfinite(sizes)
    rho = spearmanr(leftover[ok], np.log10(sizes[ok])).statistic

    # noise would make the leftover larger in the smaller structures, so compare
    # the median |residual| of the two halves by volume
    big = sizes[ok] >= np.median(sizes[ok])
    print(f"   residual against log structure volume: rho {rho:+.3f}")
    print(
        f"   |residual| in the larger half {np.median(np.abs(leftover[ok][big])):.1f} "
        f"ranks, smaller half {np.median(np.abs(leftover[ok][~big])):.1f}"
    )
    passed = abs(rho) <= BEYOND_CONTROLS["size_rho"]
    if passed:
        verdict = "size is not what the leftover is made of"
    else:
        verdict = "size drives it: look closer"
    print("   verdict: " + verdict)
    return dict(
        control="B structure size",
        number=f"rho with volume {rho:+.3f}",
        verdict="pass" if passed else "does not pass",
    )


def control_c_mice(nano: np.ndarray, columns: list[np.ndarray]) -> dict[str, str]:
    """Control C: whether one or two animals carry the leftover; its verdict row."""
    print("\nC  is it one or two animals?")
    per = per_adult_leftovers(nano, columns)
    pairs = adult_pair_agreements(per)

    # the adult whose leftover agrees least, on average, with the others
    mean_agreement = [
        np.mean(
            [spearmanr(per[i], per[j]).statistic for j in range(len(ADULTS)) if j != i]
        )
        for i in range(len(ADULTS))
    ]
    worst = int(np.argmin(mean_agreement))
    print(
        f"   each adult against each other: median rho {np.median(pairs):+.3f}, "
        f"range {min(pairs):+.3f} to {max(pairs):+.3f}"
    )
    print(
        f"   least typical animal: {ADULTS[worst]} at {mean_agreement[worst]:+.3f} "
        "mean agreement"
    )

    # the verdict rests on the worst pair, since an odd brain lowers every pair it
    # is in: below beyond_controls.pair_rho, one animal may carry the leftover
    passed = min(pairs) >= BEYOND_CONTROLS["pair_rho"]
    if passed:
        verdict = "every animal shows the same leftover"
    else:
        verdict = "one animal may be carrying it: look closer"
    print("   verdict: " + verdict)
    return dict(
        control="C single animals",
        number=f"pairwise median {np.median(pairs):+.3f}, min {min(pairs):+.3f}",
        verdict="pass" if passed else "does not pass",
    )


def control_d_groups(nano: np.ndarray, columns: list[np.ndarray]) -> dict[str, str]:
    """Control D: whether the leftover is the whisker manipulation; its verdict row."""
    print("\nD  is it the whisker manipulation? (naive and RWS are pooled)")
    naive, rws = group_leftovers(nano, columns)
    rho = spearmanr(naive, rws).statistic
    print(f"   leftover of the five naive against the five RWS: rho {rho:+.3f}")
    passed = rho >= BEYOND_CONTROLS["groups_rho"]
    if passed:
        verdict = "both groups give the same leftover, so it is not the manipulation"
    else:
        verdict = "the groups disagree: pooling is hiding something"
    print("   verdict: " + verdict)
    return dict(
        control="D naive vs RWS",
        number=f"rho {rho:+.3f}",
        verdict="pass" if passed else "does not pass",
    )


def curvature_scores(
    y: np.ndarray, covariates: dict[str, np.ndarray], terms: dict[str, tuple[str, ...]]
) -> tuple[float, float, float]:
    """The CV R2 of the model's predictors straight, to cubes and to fifth powers."""
    straight = predictors(covariates, terms)
    linear = cv_r2(y, straight)
    cubic = cv_r2(y, flexible(straight))
    fifth = flexible(straight) + [x**4 for x in straight] + [x**5 for x in straight]
    quintic = cv_r2(y, fifth)
    return linear, cubic, quintic


def control_e_curvature(
    y: np.ndarray, covariates: dict[str, np.ndarray], terms: dict[str, tuple[str, ...]]
) -> dict[str, str]:
    """Control E: whether the straight model misses curvature.

    The main model's predictors (`terms`) straight (the model), to cubes (the check
    row curved) and to fifth powers, each scored on held-out structures. If the
    cubes buy prediction, part of the leftover is curvature the straight model
    misses, and the curved row is the bound to read beside it; if fifth powers buy
    more still, even that row is not bent enough. Returns the verdict row.
    """
    print("\nE  does the straight model miss curvature?")
    linear, cubic, quintic = curvature_scores(y, covariates, terms)
    print(f"   cross-validated R2, straight              {linear:+.3f}   <- the model")
    print(f"   cross-validated R2, squares and cubes     {cubic:+.3f}   <- row curved")
    print(f"   cross-validated R2, up to fifth powers    {quintic:+.3f}")
    gain = BEYOND_CONTROLS["curvature_gain"]
    passed = cubic - linear <= gain
    if passed:
        verdict = "curving the terms buys nothing, so the straight model is adequate"
    elif quintic - cubic <= gain:
        verdict = (
            "curving the terms pays, so part of the leftover is curvature: read the "
            "check row curved, against its own floor, beside the main value"
        )
    else:
        verdict = (
            "curvature still pays past cubes: even the check row curved is not bent "
            "enough"
        )
    print("   verdict: " + verdict)
    return dict(
        control="E curvature",
        number=f"CV R2 {linear:.3f} straight, {cubic:.3f} cubic, {quintic:.3f} quintic",
        verdict="pass" if passed else "does not pass",
    )


def components(m: np.ndarray) -> np.ndarray:
    """The shared patterns of a set of genes across structures, strongest first.

    `m` is genes x structures of ranks; each gene is z-scored, the mean profile
    removed, and the rows of vt are returned (components x structures).
    """
    m = (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, keepdims=True)
    _, _, vt = np.linalg.svd(m - m.mean(axis=0), full_matrices=False)
    return vt


def best_k(y: np.ndarray, pcs: np.ndarray, labels: list[np.ndarray]) -> int:
    """The number of leading components with the best held-out R2 (folds `labels`)."""
    scores = [
        cv_r2(y, list(pcs[:k]), labels) for k in range(1, BEYOND_CONTROLS["max_pcs"] + 1)
    ]
    return int(np.argmax(scores)) + 1


def nested_cv_r2(
    y: np.ndarray, pcs: np.ndarray, seed: int = 0
) -> tuple[float, list[int]]:
    """Held-out R2 of the gene-space model with its size picked inside each fold.

    beyond_controls.nested_repeats shufflings of five outer folds; in each training
    part the number of components is picked by its own cross-validation
    (beyond_controls.nested_repeats inner shufflings), then fitted and the outer
    fold predicted. Returns the R2 and the numbers of components picked.
    """
    n = len(y)
    repeats = BEYOND_CONTROLS["nested_repeats"]
    outer = fold_labels(n, repeats=repeats, seed=seed)
    predicted = np.empty((repeats, n))
    picked = []
    for r, label in enumerate(outer):
        for fold in np.unique(label):
            test = label == fold
            train = ~test
            inner = fold_labels(int(train.sum()), repeats=repeats, seed=seed + 1 + r)
            k = best_k(y[train], pcs[:, train], inner)
            picked.append(k)
            design = np.column_stack(list(pcs[:k]) + [np.ones(n)])
            beta = np.linalg.lstsq(design[train], y[train], rcond=None)[0]
            predicted[r, test] = design[test] @ beta
    missed = np.var(y[None, :] - predicted, axis=1).mean()
    return float(1 - missed / np.var(y)), picked


def gene_space_calibration(inputs: Inputs, seed: int = 0) -> pd.DataFrame:
    """Control F's own floor: its model on maps made of the genes' expression.

    As adult.beyond_calibration, with control F's model: each gene's experiments
    in two halves; on the calibration structures, the components of the genes
    measured there in both halves, from each half. The nano map is read with each
    half's components; a known map is the model fitted to nano with one half's
    components (the number picked by cross-validation), analysed with the other
    half's, beyond_controls.f_draws draws of made-up animals each way. Columns:
    map, truth_from, predictors_from, draw, n_structures, n_genes, ceiling, cv_r2,
    left.
    """
    halves, _ = beyond_calibration.half_profiles()
    subset = beyond_calibration.calibration_structures(
        inputs.structures, halves, model_genes(inputs.terms)
    )
    positions = [inputs.structures.index(x) for x in subset]
    nano = inputs.nano[:, positions]
    genes = sorted(
        g
        for g in halves["A"]
        if all(x in halves[h][g] for h in beyond_calibration.HALVES for x in subset)
    )
    pcs = {
        h: components(gene_matrix(halves[h], subset, genes))
        for h in beyond_calibration.HALVES
    }
    splits = half_splits()
    y = full_map(nano)
    agreement, explainable = ceiling(nano, splits)
    h_nano = float(np.mean(agreement))
    rows = []
    for h in beyond_calibration.HALVES:
        cv, _ = nested_cv_r2(y, pcs[h], seed=seed)
        rows.append(
            dict(
                map=beyond_calibration.NANO,
                truth_from="",
                predictors_from=h,
                draw=-1,
                ceiling=explainable,
                cv_r2=cv,
                left=1 - cv / explainable,
            )
        )
    halves_named = beyond_calibration.HALVES
    children = np.random.SeedSequence(seed).spawn(len(halves_named))
    for h, child in zip(halves_named, children):
        other = beyond_calibration.other_half(h)
        k = best_k(y, pcs[h], fold_labels(len(y)))
        truth = y - residual(y, list(pcs[h][:k]))
        rng = np.random.default_rng(child)
        for draw in range(BEYOND_CONTROLS["f_draws"]):
            cohort = beyond_calibration.fake_cohort(truth, h_nano, rng)
            _, made_up = ceiling(cohort, splits)
            cv, _ = nested_cv_r2(full_map(cohort), pcs[other], seed=seed)
            rows.append(
                dict(
                    map="gene space",
                    truth_from=h,
                    predictors_from=other,
                    draw=draw,
                    ceiling=made_up,
                    cv_r2=cv,
                    left=1 - cv / made_up,
                )
            )
    out = pd.DataFrame(rows)
    out.insert(1, "n_structures", len(subset))
    out.insert(2, "n_genes", len(genes))
    return out


def gene_space_curve(
    y: np.ndarray, vt: np.ndarray, explainable: float, labels: list[np.ndarray]
) -> pd.DataFrame:
    """gene_space.csv: per number of leading components, R2 fitted and held out.

    `vt` holds the components (components x structures), 1 to
    beyond_controls.max_pcs of them in turn; share_of_ceiling is the held-out R2
    over the ceiling.
    """
    rows = []
    for k in range(1, BEYOND_CONTROLS["max_pcs"] + 1):
        pcs = [vt[i] for i in range(k)]
        cv = cv_r2(y, pcs, labels)
        rows.append(
            dict(
                n_components=k,
                r2=r_squared(y, pcs),
                cv_r2=cv,
                share_of_ceiling=cv / explainable,
            )
        )
    return pd.DataFrame(rows)


def control_f_gene_space(
    inputs: Inputs,
    y: np.ndarray,
    explainable: float,
    splits: list[tuple[list[int], list[int]]],
) -> tuple[dict[str, str], pd.DataFrame, pd.DataFrame, dict]:
    """Control F, the strongest: any combination of the gene table's genes may try.

    The genes measured in every structure are reduced to principal components. The
    curve of held-out R2 against the number of components (gene_space_curve) is
    drawn; the share quoted picks the number inside each training fold
    (nested_cv_r2). The verdict reads the model's own calibration
    (gene_space_calibration): the leftover passes when nano's share left, with
    either half's components, stands above every draw of a map made of the genes'
    expression. Returns the verdict row, the curve, the calibration and a summary
    (gene_space_summary.csv), k_most among it, the number picked most often.
    """
    s = inputs.structures
    genes = sorted(g for g in inputs.expr if all(x in inputs.expr[g] for x in s))
    print(
        "\nF  is it the choice of predictors? (give the model all "
        f"{len(genes)} genes measured in every structure)"
    )

    # the panel's shared patterns across structures, strongest first
    vt = components(gene_matrix(inputs.expr, s, genes))
    curve = gene_space_curve(y, vt, explainable, fold_labels(len(y)))

    # the share with the number of components picked inside each training fold, and
    # the leftover's replication with the number picked most often
    nested, picked = nested_cv_r2(y, vt)
    k_most = int(pd.Series(picked).mode().iloc[0])
    rep = replication(inputs.nano, [vt[i] for i in range(k_most)], splits)
    peak = int(curve.loc[curve["cv_r2"].idxmax(), "n_components"])
    print(
        f"   held-out R2 peaks at {peak} components "
        f"({curve['share_of_ceiling'].max():.0%} of the ceiling); picked inside each "
        f"fold ({min(picked)} to {max(picked)}, most often {k_most}): "
        f"{nested / explainable:.0%} of the ceiling, {1 - nested / explainable:.0%} left"
    )
    print(f"   the leftover with {k_most} components replicates at {rep:.3f}")

    # its own floor: the same model on maps made of the genes' expression
    calibration = gene_space_calibration(inputs)
    nano_left = calibration.loc[calibration["map"] == beyond_calibration.NANO, "left"]
    known = calibration.loc[calibration["map"] == "gene space", "left"]
    print(
        f"   on {int(calibration['n_structures'].iloc[0])} structures and "
        f"{int(calibration['n_genes'].iloc[0])} genes measured in both halves: nano "
        f"leaves {nano_left.min():.0%} to {nano_left.max():.0%}; a map made of the "
        f"genes' expression {known.min():.0%} to {known.max():.0%} "
        f"({len(known)} draws)"
    )
    passed = bool(nano_left.min() > known.max())
    if passed:
        verdict = (
            "even the gene table's components leave more than they leave of a map made "
            "of those genes"
        )
    else:
        verdict = "the gene table accounts for the map as far as Allen mismatch allows"
    print("   verdict: " + verdict)
    row = dict(
        control="F whole gene space",
        number=f"{len(genes)} genes, {k_most} components most often, nested CV R2 "
        f"{nested:.3f} ({nested / explainable:.0%} of ceiling); nano leaves "
        f"{nano_left.min():.0%} to {nano_left.max():.0%} against "
        f"{known.min():.0%} to {known.max():.0%}; replication {rep:.3f}",
        verdict="pass" if passed else "does not pass",
    )
    summary = dict(
        genes=len(genes),
        k_most=k_most,
        k_min=min(picked),
        k_max=max(picked),
        curve_peak=peak,
        nested_cv_r2=nested,
        share=nested / explainable,
        replication=rep,
    )
    return row, curve, calibration, summary


def fold_row(
    kind: str,
    key: str,
    label: str,
    y: np.ndarray,
    covariates: dict[str, np.ndarray],
    terms: dict[str, tuple[str, ...]],
    explainable: float,
    labels: list[np.ndarray] | None = None,
) -> dict:
    """One row of folds.csv: the main model's shares under folds `labels`.

    Gria1 alone, density alone, both, and the share left (beyond_density.model_shares),
    of the reproducible map, with folds `labels` (the main model's by default).
    """
    shares = model_shares(y, covariates, terms, explainable, labels)
    return dict(
        kind=kind,
        key=key,
        folds=label,
        n_structures=len(y),
        abundance=shares["abundance"],
        density_alone=shares["density"],
        model=shares["model"],
        left=shares["left"],
        lo=np.nan,
        hi=np.nan,
    )


def single_shufflings_row(
    y: np.ndarray,
    covariates: dict[str, np.ndarray],
    terms: dict[str, tuple[str, ...]],
    explainable: float,
) -> dict:
    """The spread a single shuffling of the folds would have, as one row of folds.csv.

    beyond_controls.n_single single shufflings, each seeded from one generator:
    their median row, with the 95% range of the share left (lo, hi).
    """
    n = len(y)
    rng = np.random.default_rng(0)
    single = []
    for _ in range(BEYOND_CONTROLS["n_single"]):
        labels = fold_labels(n, repeats=1, seed=int(rng.integers(1 << 30)))
        single.append(
            fold_row("folds", "", "", y, covariates, terms, explainable, labels)
        )
    table = pd.DataFrame(single)
    row = table.median(numeric_only=True).to_dict()
    row.update(
        kind="folds",
        key="single_shufflings",
        folds=f"single shufflings ({len(single)})",
        n_structures=n,
        lo=float(np.percentile(table["left"], 2.5)),
        hi=float(np.percentile(table["left"], 97.5)),
    )
    return row


def fold_rows(
    inputs: Inputs,
    covariates: dict[str, np.ndarray],
    y: np.ndarray,
    explainable: float,
    xyz: np.ndarray,
) -> pd.DataFrame:
    """folds.csv: the main model under its own folds and under others.

    Its folds (beyond.cv_repeats shufflings of five), one shuffling alone (the first
    of them), single shufflings (single_shufflings_row), ten folds, leave one out,
    and folds of spatial blocks of the centroids `xyz`.
    """
    n = len(y)
    terms = inputs.terms
    rows = [
        fold_row(
            "main",
            "main",
            f"the main model ({BEYOND['cv_repeats']} shufflings of five folds)",
            y,
            covariates,
            terms,
            explainable,
        ),
        fold_row(
            "folds",
            "one_shuffling",
            "one shuffling (the first of the main model's)",
            y,
            covariates,
            terms,
            explainable,
            fold_labels(n, repeats=1),
        ),
    ]
    rows.append(single_shufflings_row(y, covariates, terms, explainable))
    folds = [
        ("ten_folds", "ten folds", fold_labels(n, folds=10)),
        ("leave_one_out", "leave one out", [np.arange(n)]),
        ("spatial_blocks", "folds of spatial blocks", block_labels(xyz)),
    ]
    for key, label, labels in folds:
        rows.append(
            fold_row("folds", key, label, y, covariates, terms, explainable, labels)
        )
    return pd.DataFrame(rows)


def reading_matrices(inputs: Inputs) -> dict[str, np.ndarray]:
    """Every reading of the map as adults x structures, where every adult has it.

    zref is the map itself; the others come from the young-against-adult table of
    region_plot (zref_stored is zref with its 17-brain reference).
    """
    table = pd.read_csv(REGION_MEANS)
    table = table[table["mouse"].isin(ADULTS)]
    out = {"zref": inputs.nano}
    for reading in ("zref_stored", "cref", "subref", "ratio", "sepratio"):
        mine = table[table["reading"] == reading.replace("_stored", "")]
        wide = mine.pivot(index="mouse", columns="structure", values="log2_value")
        wide = wide.reindex(index=ADULTS, columns=inputs.structures)
        if wide.isna().any().any():
            print(f"   {reading:12s} not measured in every structure; skipped")
            continue
        out[reading] = wide.to_numpy(float)
    return out


def control_g_readings(
    inputs: Inputs, columns: list[np.ndarray], splits: list[tuple[list[int], list[int]]]
) -> tuple[dict[str, str], pd.DataFrame]:
    """Control G: whether any of this is specific to zref and its reference.

    The same predictors, ceiling and replication for every reading measured in
    every structure. Returns the verdict row and, per reading, the R2 of the model
    (in-sample), the map's replication and the leftover's.
    """
    print("\nG  is it zref? (the same test on every reading)")
    rows = []
    for reading, matrix in reading_matrices(inputs).items():
        y = full_map(matrix)
        agreement, _ = ceiling(matrix, splits)
        rows.append(
            dict(
                reading=reading,
                r2=r_squared(y, columns),
                map_replication=float(np.mean(agreement)),
                leftover_replication=replication(matrix, columns, splits),
            )
        )
        r = rows[-1]
        print(
            f"   {reading:12s} model R2 {r['r2']:.3f}; map replicates "
            f"{r['map_replication']:.3f}, leftover {r['leftover_replication']:.3f}"
        )
    table = pd.DataFrame(rows)
    passed = bool(
        (table["leftover_replication"] > BEYOND_CONTROLS["readings_replication"]).all()
    )
    if passed:
        verdict = "the leftover replicates under every reading"
    else:
        verdict = "some readings disagree: look closer"
    print("   verdict: " + verdict)
    return (
        dict(
            control="G reading choice",
            number="; ".join(
                f"{r.reading} {r.leftover_replication:.2f}" for r in table.itertuples()
            ),
            verdict="pass" if passed else "does not pass",
        ),
        table,
    )


def write_verdicts(verdicts: list[dict[str, str]]) -> None:
    """Write controls.csv and print the verdict."""
    pd.DataFrame(verdicts).to_csv(CONTROLS, index=False)
    print(f"\n-> {CONTROLS}")
    failed = [v["control"] for v in verdicts if v["verdict"] != "pass"]
    verdict = "every control passes" if not failed else f"look closer at {failed}"
    print("verdict: " + verdict)


def main() -> None:
    """Run the seven controls and the other folds, and write their tables."""
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs()
    s = inputs.structures
    covariates, _, _ = covariates_for(inputs)
    columns = model_columns(covariates, inputs.terms)
    splits = half_splits()
    y = full_map(inputs.nano)
    leftover = residual(y, columns)
    _, explainable = ceiling(inputs.nano, splits)
    print(
        f"{len(s)} structures, ceiling {explainable:.3f}, the model's leftover "
        f"replicates at {replication(inputs.nano, columns, splits):.3f}"
    )

    # the centroids in one hemisphere, and each structure's mean volume
    xyz = structures.centroid_xyz(s)
    sizes = mean_sizes(s)

    # the seven controls
    verdicts = [
        control_a_space(leftover, inputs.nano, columns, xyz, splits),
        control_b_size(leftover, sizes),
    ]
    vc = control_c_mice(inputs.nano, columns)
    vd = control_d_groups(inputs.nano, columns)
    ve = control_e_curvature(y, covariates, inputs.terms)
    vf, curve, f_calibration, f_summary = control_f_gene_space(
        inputs, y, explainable, splits
    )
    f_calibration.to_csv(GENE_SPACE_CALIBRATION, index=False)
    pd.DataFrame([f_summary]).to_csv(GENE_SPACE_SUMMARY, index=False)
    vg, readings = control_g_readings(inputs, columns, splits)
    verdicts += [vc, vd, ve, vf, vg]
    write_verdicts(verdicts)
    curve.to_csv(GENE_SPACE, index=False)
    readings.to_csv(READINGS, index=False)

    # the main model under other folds
    table = fold_rows(inputs, covariates, y, explainable, xyz)
    table.to_csv(FOLDS, index=False)
    print("\nThe main model under its folds and others (share of the reproducible map):")
    for r in table.itertuples():
        spread = f" ({r.lo:.1%} to {r.hi:.1%})" if np.isfinite(r.lo) else ""
        print(
            f"   {r.folds[:50]:52s} Gria1 {r.abundance:6.1%}, density "
            f"{r.density_alone:6.1%}, both {r.model:6.1%}, left {r.left:6.1%}{spread}"
        )
