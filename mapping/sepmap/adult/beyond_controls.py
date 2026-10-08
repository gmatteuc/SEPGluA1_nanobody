"""Seven attempts to break the result of adult.beyond_density.

That module finds part of the adult map's reproducible pattern not predicted by
receptor mRNA or synaptic density, with a leftover that replicates across mice. Each
control below is a specific way that could be an artefact, with the number that says
whether it is:

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
                                      line leaves real structure in the residual;
                                      this control is why the model bends, and it
                                      asks whether it bends enough
    F  the choice of predictors?      the strongest version: give the model the
                                      expression of every gene of the gene table
                                      measured in all the structures, as principal
                                      components (as many as cross-validation
                                      picks, inside each training fold), and ask
                                      whether what it leaves stands above what the
                                      same model leaves of a map made of those
                                      genes' expression (its own calibration
                                      floor). It is why the leftover is "not
                                      predicted by receptor mRNA or synaptic
                                      density" and never "beyond gene expression"
    G  zref?                          the same test on the other readings of the
                                      map, and on zref with the 17-brain reference

E and F are cross-validated, because a flexible model always fits better on the
data it was fitted to; the question is whether it predicts better, and only
held-out structures can say. F picks its number of components inside each training
fold (nested), since the best of many held-out scores is itself optimistic.
Everything is imported from adult.beyond_density, so the two modules use the same
structures, predictors and arithmetic.

Beside the controls, the variants: the quoted leftover under other folds (one
shuffling, single shufflings, ten folds, leave one out, spatial blocks) and on the
larger set of structures that the rule keeps once the markers measured by a single
Allen experiment are left out.

Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    controls.csv            one row per control, with its number and verdict
    gene_space.csv          control F: per number of components, R2 fitted and held out
    gene_space_calibration.csv  control F's model on the nano map and on maps made of
                            the genes' expression, each half of the experiments
    gene_space_summary.csv  control F's numbers: genes, components picked, nested
                            share of the ceiling, replication
    readings.csv            control G: per reading, R2 and the two replications
    variants.csv            the leftover under other folds and structures
    fig4_controls.png       A to D, the four artefact checks
    fig5_model_space.png    E and F, how much any model of this data can explain
    fig6_readings.png       G, the same test on every reading

Run by run_beyond_controls.py.
"""

from collections.abc import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from sepmap.adult import beyond_calibration as bc
from sepmap.adult import profiles
from sepmap.adult.beyond_density import (
    ADULTS,
    MARKERS,
    NAIVE_ROWS,
    OUT,
    RWS_ROWS,
    Inputs,
    block_labels,
    budget,
    build_covariates,
    ceiling,
    cv_r2,
    flexible,
    fold_labels,
    full_map,
    half_splits,
    leftover_agreement,
    load_inputs,
    model,
    per_adult_leftovers,
    predictors,
    r_squared,
    replicates,
    residual,
    save,
)
from sepmap.config import SETTINGS
from sepmap.ish import gene_table
from sepmap.plotting import RED, tidy
from sepmap.structures import load_centroids
from sepmap.young_vs_adult.region_plot import REGION_MEANS

# the largest gene-space model of control F, and the threshold of each verdict
BEYOND_CONTROLS = SETTINGS["beyond_controls"]

CONTROLS = OUT / "controls.csv"
GENE_SPACE = OUT / "gene_space.csv"
GENE_SPACE_CALIBRATION = OUT / "gene_space_calibration.csv"
GENE_SPACE_SUMMARY = OUT / "gene_space_summary.csv"
READINGS = OUT / "readings.csv"
VARIANTS = OUT / "variants.csv"

# the readings of control G as the figure names them: the map's own zref, the same
# with the young-against-adult tables' 17-brain reference, and the stored readings
READING_LABELS = {
    "zref": "zref\n(declared reference)",
    "zref_stored": "zref\n(17-brain reference)",
    "cref": "cref",
    "subref": "subref",
    "ratio": "nano / auto",
    "sepratio": "nano / SEP",
}


# ===== Utilities =====


def replication(
    matrix: np.ndarray,
    columns: Sequence[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
) -> float:
    """How well the leftover of one half-cohort matches the other's, over every split."""
    return float(np.mean(leftover_agreement(matrix, columns, splits)))


def gene_matrix(
    expr: dict[str, dict[str, float]], structures: list[str], genes: list[str]
) -> np.ndarray:
    """Rank profiles of many genes as one array, genes by structures."""
    # ranks, since each Allen experiment has its own arbitrary intensity scale
    return np.array([rankdata([expr[g][s] for s in structures]) for g in genes])


# ===== Controls =====


def control_a_space(
    res: np.ndarray,
    structures: list[str],
    nano: np.ndarray,
    xs: list[np.ndarray],
    coords: dict[str, np.ndarray],
    splits: list[tuple[list[int], list[int]]],
) -> dict[str, str] | None:
    """Control A: whether the leftover is just a smooth gradient across the block.

    `coords` holds each structure's centroid (AP, DV, ML) in mm. Returns the verdict
    row, or None when fewer than beyond_controls.min_centroids structures have one.
    """
    print("\nA  is it a smooth spatial gradient? (an illumination artefact)")

    # the structures with a centroid; on too few, a fit of six position terms would
    # explain much of the leftover by chance
    xyz = np.array([coords.get(s, np.full(3, np.nan)) for s in structures])
    ok = np.all(np.isfinite(xyz), axis=1)
    if ok.sum() < BEYOND_CONTROLS["min_centroids"]:
        print("   not enough centroids; skipped")
        return None

    # the residual against a quadratic in the three ranked positions, enough to
    # follow a ramp or a bowl across the block, as uneven illumination would give
    pos = [rankdata(xyz[ok, i]) for i in range(3)]
    quad = pos + [p**2 for p in pos]
    smooth = r_squared(res[ok], quad)
    for axis, p in zip("AP DV ML".split(), pos):
        print(f"   residual against {axis}: rho {spearmanr(res[ok], p).statistic:+.3f}")
    print(f"   a smooth quadratic in all three axes explains R2 = {smooth:.3f} of it")

    # the replication with position added to the model, the ranked axes straight
    if ok.all():
        with_pos = replication(nano, xs + [rankdata(xyz[:, i]) for i in range(3)], splits)
        still = "still replicates" if replicates(with_pos) else "replicates only"
        print(
            f"   with position added to the model the leftover {still} at {with_pos:.3f}"
        )
        number = f"smooth R2 {smooth:.3f}, replication with position {with_pos:.3f}"
    else:
        number = f"smooth R2 {smooth:.3f}"

    # a gradient could explain the leftover if position explains more than
    # beyond_controls.gradient_r2 of it
    passed = smooth <= BEYOND_CONTROLS["gradient_r2"]
    verdict = (
        "not a gradient; position is a weak predictor of it"
        if passed
        else "a gradient could explain it -- LOOK CLOSER"
    )
    print("   verdict: " + verdict)
    return dict(
        control="A spatial gradient",
        number=number,
        verdict="pass" if passed else "CHECK",
    )


def control_b_size(res: np.ndarray, sizes: np.ndarray) -> dict[str, str]:
    """Control B: whether the leftover comes from small or poorly covered structures.

    `sizes` holds each structure's mean volume over the adults in 20 um voxels.
    Returns the verdict row.
    """
    print("\nB  is it small structures, where a mean is noisy?")
    ok = np.isfinite(sizes)
    rho = spearmanr(res[ok], np.log10(sizes[ok])).statistic

    # noise would make the leftover larger in the smaller structures, so compare
    # the median |residual| of the two halves by volume
    big = sizes[ok] >= np.median(sizes[ok])
    print(f"   residual against log structure volume: rho {rho:+.3f}")
    print(
        f"   |residual| in the larger half {np.median(np.abs(res[ok][big])):.1f} ranks, "
        f"smaller half {np.median(np.abs(res[ok][~big])):.1f}"
    )
    passed = abs(rho) <= BEYOND_CONTROLS["size_rho"]
    verdict = (
        "size is not what the leftover is made of"
        if passed
        else "size drives it -- LOOK CLOSER"
    )
    print("   verdict: " + verdict)
    return dict(
        control="B structure size",
        number=f"rho with volume {rho:+.3f}",
        verdict="pass" if passed else "CHECK",
    )


def control_c_mice(
    nano: np.ndarray, xs: list[np.ndarray]
) -> tuple[dict[str, str], list[float]]:
    """Control C: whether the leftover is carried by one or two animals.

    Returns the verdict row and the agreement of every pair of adults.
    """
    print("\nC  is it one or two animals?")
    per = per_adult_leftovers(nano, xs)
    pairs = [
        float(spearmanr(per[i], per[j]).statistic)
        for i in range(len(ADULTS))
        for j in range(i + 1, len(ADULTS))
    ]

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
    verdict = (
        "every animal shows the same leftover"
        if passed
        else "one animal may be carrying it -- LOOK CLOSER"
    )
    print("   verdict: " + verdict)
    return (
        dict(
            control="C single animals",
            number=f"pairwise median {np.median(pairs):+.3f}, min {min(pairs):+.3f}",
            verdict="pass" if passed else "CHECK",
        ),
        pairs,
    )


def control_d_groups(
    nano: np.ndarray, xs: list[np.ndarray]
) -> tuple[dict[str, str], np.ndarray, np.ndarray]:
    """Control D: whether the leftover is the whisker manipulation, not the anatomy.

    Returns the verdict row and the leftovers of the naive and the RWS group.
    """
    print("\nD  is it the whisker manipulation? (naive and RWS are pooled)")
    naive = residual(rankdata(nano[NAIVE_ROWS].mean(axis=0)), xs)
    rws = residual(rankdata(nano[RWS_ROWS].mean(axis=0)), xs)
    rho = spearmanr(naive, rws).statistic
    print(f"   leftover of the five naive against the five RWS: rho {rho:+.3f}")
    passed = rho >= BEYOND_CONTROLS["groups_rho"]
    verdict = (
        "both groups give the same leftover, so it is not the manipulation"
        if passed
        else "the groups disagree -- pooling is hiding something"
    )
    print("   verdict: " + verdict)
    return (
        dict(
            control="D naive vs RWS",
            number=f"rho {rho:+.3f}",
            verdict="pass" if passed else "CHECK",
        ),
        naive,
        rws,
    )


def control_e_curvature(
    y: np.ndarray, covariates: dict[str, np.ndarray]
) -> tuple[dict[str, str], float, float, float]:
    """Control E: whether the bent model bends enough.

    The predictors straight, to cubes (the model) and to fifth powers, each scored
    on held-out structures. If fifth powers still buy prediction, part of the
    leftover is the model's failing. Returns the verdict row and the three CV R2.
    """
    print("\nE  is the model bent enough? (does more curvature keep paying?)")
    xs = predictors(covariates)
    linear = cv_r2(y, xs)
    cubic = cv_r2(y, flexible(xs))
    quintic = cv_r2(y, flexible(xs) + [x**4 for x in xs] + [x**5 for x in xs])
    print(f"   cross-validated R2, straight              {linear:+.3f}")
    print(f"   cross-validated R2, squares and cubes     {cubic:+.3f}   <- the model")
    print(f"   cross-validated R2, up to fifth powers    {quintic:+.3f}")
    passed = quintic - cubic <= BEYOND_CONTROLS["curvature_gain"]
    verdict = (
        "further bending buys nothing, so the model is adequate"
        if passed
        else "more curvature still pays -- the model is not bent enough"
    )
    print("   verdict: " + verdict)
    return (
        dict(
            control="E curvature",
            number=f"CV R2 {linear:.3f} straight, {cubic:.3f} cubic, "
            f"{quintic:.3f} quintic",
            verdict="pass" if passed else "CHECK",
        ),
        linear,
        cubic,
        quintic,
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
    halves_a, halves_b, _ = bc.experiment_halves(bc.per_experiment_profiles())
    halves = {"A": halves_a, "B": halves_b}
    structures = bc.calibration_structures(inputs.structures, halves)
    columns = [inputs.structures.index(x) for x in structures]
    nano = inputs.nano[:, columns]
    genes = sorted(
        g
        for g in halves_a
        if all(x in halves[h][g] for h in bc.HALVES for x in structures)
    )
    pcs = {h: components(gene_matrix(halves[h], structures, genes)) for h in bc.HALVES}
    splits = half_splits()
    y = full_map(nano)
    agreement, explainable = ceiling(nano, splits)
    h_nano = float(np.mean(agreement))
    rows = []
    for h in bc.HALVES:
        cv, _ = nested_cv_r2(y, pcs[h], seed=seed)
        rows.append(
            dict(
                map=bc.NANO,
                truth_from="",
                predictors_from=h,
                draw=-1,
                ceiling=explainable,
                cv_r2=cv,
                left=1 - cv / explainable,
            )
        )
    children = np.random.SeedSequence(seed).spawn(len(bc.HALVES))
    for h, child in zip(bc.HALVES, children):
        other = bc.HALVES[1 - bc.HALVES.index(h)]
        k = best_k(y, pcs[h], fold_labels(len(y)))
        truth = y - residual(y, list(pcs[h][:k]))
        rng = np.random.default_rng(child)
        for draw in range(BEYOND_CONTROLS["f_draws"]):
            cohort = bc.fake_cohort(truth, h_nano, rng)
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
    out.insert(1, "n_structures", len(structures))
    out.insert(2, "n_genes", len(genes))
    return out


def control_f_gene_space(
    inputs: Inputs,
    y: np.ndarray,
    explainable: float,
    splits: list[tuple[list[int], list[int]]],
) -> tuple[dict[str, str], pd.DataFrame, int, pd.DataFrame, dict]:
    """Control F, the strongest: any combination of the gene table's genes may try.

    The genes measured in every structure are reduced to principal components. The
    curve of held-out R2 against the number of components (1 to
    beyond_controls.max_pcs) is drawn; the share quoted picks the number inside
    each training fold (nested_cv_r2). The verdict reads the model's own
    calibration (gene_space_calibration): the leftover passes when nano's share
    left, with either half's components, stands above every draw of a map made of
    the genes' expression. Returns the verdict row, the curve (n_components, r2
    fitted, cv_r2 held out, share_of_ceiling), the number of components picked
    most often, the calibration and a summary (gene_space_summary.csv).
    """
    s = inputs.structures
    genes = sorted(g for g in inputs.expr if all(x in inputs.expr[g] for x in s))
    print(
        "\nF  is it the choice of predictors? (give the model all "
        f"{len(genes)} genes measured in every structure)"
    )

    # the panel's shared patterns across structures, strongest first
    vt = components(gene_matrix(inputs.expr, s, genes))
    labels = fold_labels(len(y))
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
    curve = pd.DataFrame(rows)

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
    nano_left = calibration.loc[calibration["map"] == bc.NANO, "left"]
    known = calibration.loc[calibration["map"] == "gene space", "left"]
    print(
        f"   on {int(calibration['n_structures'].iloc[0])} structures and "
        f"{int(calibration['n_genes'].iloc[0])} genes measured in both halves: nano "
        f"leaves {nano_left.min():.0%} to {nano_left.max():.0%}; a map made of the "
        f"genes' expression {known.min():.0%} to {known.max():.0%} "
        f"({len(known)} draws)"
    )
    passed = bool(nano_left.min() > known.max())
    verdict = (
        "even the gene table's components leave more than they leave of a map made "
        "of those genes"
        if passed
        else "the gene table accounts for the map as far as Allen mismatch allows"
    )
    print("   verdict: " + verdict)
    return (
        dict(
            control="F whole gene space",
            number=f"{len(genes)} genes, {k_most} components most often, nested CV R2 "
            f"{nested:.3f} ({nested / explainable:.0%} of ceiling); nano leaves "
            f"{nano_left.min():.0%} to {nano_left.max():.0%} against "
            f"{known.min():.0%} to {known.max():.0%}; replication {rep:.3f}",
            verdict="pass" if passed else "CHECK",
        ),
        curve,
        k_most,
        calibration,
        dict(
            genes=len(genes),
            k_most=k_most,
            k_min=min(picked),
            k_max=max(picked),
            curve_peak=int(curve.loc[curve["cv_r2"].idxmax(), "n_components"]),
            nested_cv_r2=nested,
            share=nested / explainable,
            replication=rep,
        ),
    )


def once_measured_markers() -> list[str]:
    """The marker genes of beyond.markers measured by a single usable Allen experiment."""
    genes = gene_table.per_gene(gene_table.load_gene_table()).set_index("symbol")
    return [g for g in MARKERS if int(genes.loc[g, "n_experiments_used"]) == 1]


def variants(
    inputs: Inputs,
    covariates: dict[str, np.ndarray],
    y: np.ndarray,
    explainable: float,
    coords: dict[str, np.ndarray],
) -> pd.DataFrame:
    """variants.csv: the quoted leftover under other folds and on other structures.

    Rows: the quoted folds (beyond.cv_repeats shufflings of five), the single
    shuffling of 26 September, single shufflings (beyond_controls.n_single of them:
    median and 95% range), ten folds, leave one out, folds of spatial blocks, and
    the model on the structures the rule keeps once the markers measured by one
    Allen experiment are left out (their ceiling recomputed there).
    """
    n = len(y)
    rng = np.random.default_rng(0)
    xyz = np.array([coords[x] for x in inputs.structures])
    folds = [
        ("quoted: 20 shufflings of five folds", None),
        ("one shuffling (the folds of 26 September)", fold_labels(n, repeats=1)),
        ("ten folds", fold_labels(n, folds=10)),
        ("leave one out", [np.arange(n)]),
        ("folds of spatial blocks", block_labels(xyz)),
    ]
    rows = []
    for label, labels in folds:
        left = 1 - budget(y, covariates, explainable, labels=labels)[-1]
        rows.append(dict(variant=label, n_structures=n, left=left, lo=np.nan, hi=np.nan))
    single = [
        1
        - budget(
            y,
            covariates,
            explainable,
            labels=fold_labels(n, repeats=1, seed=int(rng.integers(1 << 30))),
        )[-1]
        for _ in range(BEYOND_CONTROLS["n_single"])
    ]
    rows.insert(
        2,
        dict(
            variant=f"single shufflings ({len(single)})",
            n_structures=n,
            left=float(np.median(single)),
            lo=float(np.percentile(single, 2.5)),
            hi=float(np.percentile(single, 97.5)),
        ),
    )

    # the structures the rule keeps without the markers measured once
    once = once_measured_markers()
    markers = [g for g in MARKERS if g not in once]
    wider = load_inputs(markers)
    cov, _, _ = build_covariates(
        wider.expr, wider.role, wider.auto, wider.structures, markers
    )
    _, wider_ceiling = ceiling(wider.nano, half_splits())
    left = 1 - budget(full_map(wider.nano), cov, wider_ceiling)[-1]
    rows.append(
        dict(
            variant=f"markers measured once left out ({', '.join(once)})",
            n_structures=len(wider.structures),
            left=left,
            lo=np.nan,
            hi=np.nan,
        )
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
    inputs: Inputs, xs: list[np.ndarray], splits: list[tuple[list[int], list[int]]]
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
                r2=r_squared(y, xs),
                map_replication=float(np.mean(agreement)),
                leftover_replication=replication(matrix, xs, splits),
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
    verdict = (
        "the leftover replicates under every reading"
        if passed
        else "some readings disagree -- LOOK CLOSER"
    )
    print("   verdict: " + verdict)
    return (
        dict(
            control="G reading choice",
            number="; ".join(
                f"{r.reading} {r.leftover_replication:.2f}" for r in table.itertuples()
            ),
            verdict="pass" if passed else "CHECK",
        ),
        table,
    )


# ===== Figures =====


def figure_artefacts(
    res: np.ndarray,
    structures: list[str],
    sizes: np.ndarray,
    pairs: list[float],
    naive: np.ndarray,
    rws: np.ndarray,
    coords: dict[str, np.ndarray],
    passed: dict[str, bool],
) -> None:
    """Draw controls A to D: position, size, pairs of adults, naive against RWS.

    `passed` holds each control's verdict by letter (A is missing when skipped);
    each title says what its verdict says.
    """
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 3.9))

    # A: the leftover against anterior-posterior position
    xyz = np.array([coords.get(s, np.full(3, np.nan)) for s in structures])
    ok = np.all(np.isfinite(xyz), axis=1)
    axes[0].scatter(
        xyz[ok, 0], res[ok], s=12, facecolor="0.6", edgecolor="0.25", linewidth=0.3
    )
    axes[0].axhline(0, color="0.85", lw=0.7)
    axes[0].set_xlabel("structure centroid, anterior-posterior (mm)", fontsize=8)
    axes[0].set_ylabel("residual (ranks)", fontsize=8)
    if "A" not in passed:
        title = "A. not tested: too few centroids"
    elif passed["A"]:
        title = "A. not a front-to-back gradient"
    else:
        title = "A. a gradient could explain it"
    axes[0].set_title(title, fontsize=9)
    tidy(axes[0])

    # B: the leftover against structure volume
    good = np.isfinite(sizes)
    axes[1].scatter(
        np.log10(sizes[good]),
        res[good],
        s=12,
        facecolor="0.6",
        edgecolor="0.25",
        linewidth=0.3,
    )
    axes[1].axhline(0, color="0.85", lw=0.7)
    axes[1].set_xlabel("log10 structure volume (20 um voxels)", fontsize=8)
    axes[1].set_ylabel("residual (ranks)", fontsize=8)
    title = "B. not small-structure noise" if passed["B"] else "B. size may drive it"
    axes[1].set_title(title, fontsize=9)
    tidy(axes[1])

    # C: the agreement of every pair of adults, the median in red
    axes[2].hist(pairs, bins=20, color="0.7", edgecolor="0.35", linewidth=0.4)
    axes[2].axvline(float(np.median(pairs)), color=RED, lw=1.8)
    axes[2].set_xlabel("leftover of one adult against another (Spearman)", fontsize=8)
    axes[2].set_ylabel("pairs of adults", fontsize=8)
    title = "C. every animal shows it" if passed["C"] else "C. one animal may carry it"
    axes[2].set_title(title, fontsize=9)
    tidy(axes[2])

    # D: the naive leftover against the RWS one, with the identity line
    axes[3].scatter(naive, rws, s=12, facecolor="0.6", edgecolor="0.25", linewidth=0.3)
    lim = [min(naive.min(), rws.min()) - 3, max(naive.max(), rws.max()) + 3]
    axes[3].plot(lim, lim, color="0.75", ls="--", lw=0.8)
    axes[3].set_xlabel("leftover, five naive adults", fontsize=8)
    axes[3].set_ylabel("leftover, five RWS adults", fontsize=8)
    title = (
        "D. not the whisker manipulation" if passed["D"] else "D. naive and RWS disagree"
    )
    axes[3].set_title(f"{title}\nrho {spearmanr(naive, rws).statistic:+.2f}", fontsize=9)
    tidy(axes[3])

    not_ruled_out = [k for k in "ABCD" if not passed.get(k, False)]
    if not_ruled_out:
        title = (
            "Four ways the leftover could be an artefact; not ruled out: "
            + ", ".join(not_ruled_out)
        )
    else:
        title = "Four ways the leftover could be an artefact, and is not"
    fig.suptitle(title, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, "fig4_controls.png")


def figure_model_space(
    curve: pd.DataFrame,
    best_k: int,
    explainable: float,
    cubic: float,
    quintic: float,
    passed: dict[str, bool],
) -> None:
    """Draw controls F and E: the gene-space curve, and bending further."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ks = curve["n_components"]
    axes[0].plot(ks, curve["r2"], color="0.6", lw=1.5, label="fitted")
    axes[0].plot(ks, curve["cv_r2"], color=RED, lw=1.8, label="cross-validated")
    axes[0].axhline(explainable, color="0.3", ls="--", lw=1.2)
    axes[0].annotate(
        "ceiling", (ks.iloc[-1], explainable), fontsize=7.5, ha="right", va="bottom"
    )
    axes[0].axvline(best_k, color="0.4", ls=":", lw=1.0)
    axes[0].set_xlabel(
        "components of the expression of the genes measured in every structure",
        fontsize=8,
    )
    axes[0].set_ylabel("variance of the map explained", fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False)
    title = (
        "F. even the whole gene table falls short"
        if passed["F"]
        else "F. the whole gene table accounts for the map"
    )
    axes[0].set_title(
        f"{title}\nthe gap between the two lines is overfitting", fontsize=9
    )
    tidy(axes[0])

    axes[1].bar(
        [0, 1], [cubic, quintic], color=[RED, "0.65"], edgecolor="0.25", linewidth=0.5
    )
    axes[1].set_xticks([0, 1])
    axes[1].set_xticklabels(
        ["the model\n(squares, cubes)", "bent further\n(to fifth powers)"], fontsize=8
    )
    axes[1].set_ylabel("cross-validated R2", fontsize=8)
    title = (
        "E. and bending it further buys nothing"
        if passed["E"]
        else "E. and bending it further still pays"
    )
    axes[1].set_title(title, fontsize=9)
    tidy(axes[1])
    fig.tight_layout()
    save(fig, "fig5_model_space.png")


def figure_readings(table: pd.DataFrame, passed: bool) -> None:
    """Draw control G: per reading, the map's and the leftover's replication, and R2."""
    fig, ax = plt.subplots(figsize=(8.0, 4.0))
    x = np.arange(len(table))
    ax.bar(
        x - 0.22,
        table["map_replication"],
        width=0.2,
        color="0.55",
        edgecolor="0.25",
        linewidth=0.4,
        label="the map replicates",
    )
    ax.bar(
        x,
        table["leftover_replication"],
        width=0.2,
        color=RED,
        edgecolor="0.25",
        linewidth=0.4,
        label="the leftover replicates",
    )
    ax.bar(
        x + 0.22,
        table["r2"],
        width=0.2,
        color="0.8",
        edgecolor="0.25",
        linewidth=0.4,
        label="the model explains (R2)",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([READING_LABELS[r] for r in table["reading"]], fontsize=7.5)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    title = (
        "G. the same picture under every reading, not just zref"
        if passed
        else "G. not the same picture under every reading"
    )
    ax.set_title(title, fontsize=9)
    tidy(ax)
    fig.tight_layout()
    save(fig, "fig6_readings.png")


def mean_sizes(structures: list[str]) -> np.ndarray:
    """Each structure's mean volume over the adults, in 20 um voxels."""
    table = profiles.load_per_mouse()
    sizes = table[table["mouse"].isin(ADULTS)].groupby("structure")["n_vox20"].mean()
    return sizes.reindex(structures).to_numpy(float)


def write_verdicts(verdicts: list[dict[str, str] | None]) -> None:
    """Write controls.csv, a skipped control (None) left out, and print the verdict."""
    pd.DataFrame([v for v in verdicts if v]).to_csv(CONTROLS, index=False)
    print(f"\n-> {CONTROLS}")
    failed = [v["control"] for v in verdicts if v and v["verdict"] != "pass"]
    verdict = "every control passes" if not failed else f"look closer at {failed}"
    print("VERDICT: " + verdict)


def main() -> None:
    """Run the seven controls, write their verdicts and draw them."""
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs()
    s = inputs.structures
    covariates, _, _ = build_covariates(inputs.expr, inputs.role, inputs.auto, s)
    xs = model(covariates)
    splits = half_splits()
    y = full_map(inputs.nano)
    res = residual(y, xs)
    _, explainable = ceiling(inputs.nano, splits)
    print(
        f"{len(s)} structures, ceiling {explainable:.3f}, the model's leftover "
        f"replicates at {replication(inputs.nano, xs, splits):.3f}"
    )

    # the centroids in one hemisphere, and each structure's mean volume
    centroids = load_centroids()
    coords = {
        name: centroids.loc[name, ["ap_mm", "dv_mm", "ml_mm"]].to_numpy(float)
        for name in s
        if name in centroids.index
    }
    sizes = mean_sizes(s)

    verdicts = [
        control_a_space(res, s, inputs.nano, xs, coords, splits),
        control_b_size(res, sizes),
    ]
    vc, pairs = control_c_mice(inputs.nano, xs)
    vd, naive, rws = control_d_groups(inputs.nano, xs)
    ve, _, cubic, quintic = control_e_curvature(y, covariates)
    vf, curve, best_k, f_calibration, f_summary = control_f_gene_space(
        inputs, y, explainable, splits
    )
    f_calibration.to_csv(GENE_SPACE_CALIBRATION, index=False)
    pd.DataFrame([f_summary]).to_csv(GENE_SPACE_SUMMARY, index=False)
    vg, readings = control_g_readings(inputs, xs, splits)
    verdicts += [vc, vd, ve, vf, vg]
    write_verdicts(verdicts)
    curve.to_csv(GENE_SPACE, index=False)
    readings.to_csv(READINGS, index=False)

    passed = {v["control"][0]: v["verdict"] == "pass" for v in verdicts if v}
    figure_artefacts(res, s, sizes, pairs, naive, rws, coords, passed)
    figure_model_space(curve, best_k, explainable, cubic, quintic, passed)
    figure_readings(readings, passed["G"])

    # the leftover under other folds and on other structures
    table = variants(inputs, covariates, y, explainable, coords)
    table.to_csv(VARIANTS, index=False)
    print("\nVariants of the quoted leftover:")
    for r in table.itertuples():
        spread = f" ({r.lo:.1%} to {r.hi:.1%})" if np.isfinite(r.lo) else ""
        print(f"   {r.variant:60s} {r.n_structures:4d} structures  {r.left:6.1%}{spread}")
