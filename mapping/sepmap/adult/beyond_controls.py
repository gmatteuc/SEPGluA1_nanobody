"""Seven attempts to break the claim of adult.beyond_density.

That module finds that about 40% of the explainable variance of the adult
surface-GluA1 map is not accounted for by receptor abundance or synaptic density,
and that what remains replicates across independent animals. Each control below is
a specific way the claim could be wrong, with the number that says whether it is:

    A  a smooth spatial gradient?     an illumination artefact would look like
                                      one; anatomy would not
    B  small or poorly covered        noisy means look like signal
       structures?
    C  one or two animals?            a single odd brain can carry a cohort mean
                                      of ten
    D  the whisker manipulation?      naive and RWS are pooled; if the leftover is
                                      the manipulation, the pooling hides it
    E  curvature the model misses?    a bent rank relationship fitted as a straight
                                      line leaves real structure in the residual;
                                      this control is why the quoted model bends,
                                      and it asks whether it bends enough
    F  our choice of covariates?      the strongest version: give the model the
                                      expression of every panel gene measured in
                                      all the structures (253 of the 390) and see
                                      whether the leftover survives
    G  zref?                          the same test on every other reading

E and F are cross-validated, because a flexible model always fits better on the
data it was fitted to; the question is whether it predicts better, and only
held-out structures can say. Everything is imported from adult.beyond_density, so
the two modules use the same structures, covariates and arithmetic.

Writes, in adult_v2/beyond/ under the data root:

    controls.csv            one row per control, with its verdict
    fig4_controls.png       A to D, the four artefact checks
    fig5_model_space.png    E and F, how much any model of this data can explain
    fig6_readings.png       G, the same test on all five readings

Run by run_beyond_controls.py.
"""

import csv
from collections import defaultdict
from collections.abc import Sequence

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata, spearmanr

from sepmap.adult.beyond_density import (
    ADULTS,
    NAIVE,
    OUT,
    RWS,
    build_covariates,
    cv_r2,
    flexible,
    half_map,
    half_splits,
    prepare,
    r_squared,
    replicates,
    residual,
    save,
    spearman_brown,
)
from sepmap.config import SETTINGS
from sepmap.plotting import RED, tidy
from sepmap.volumes.per_mouse import annotation_20, structure_terms
from sepmap.young_vs_adult.region_plot import REGION_MEANS

# the largest gene-space model of control F, and the threshold of each verdict
BEYOND_CONTROLS = SETTINGS["beyond_controls"]


# ===== Utilities =====


def centroids(structures: list[str]) -> dict[str, np.ndarray]:
    """Mean (AP, DV, ML) position of each structure, in mm, from the CCF itself.

    Needed by control A: if the leftover were an imaging artefact it would vary
    smoothly with position in the block, so the first thing to ask of it is how
    much a smooth function of position can explain. A structure not in the atlas
    gets NaN; when none of them is, ValueError, since control A and the artefact
    figure both read every structure's centroid.
    """
    # the structure name of each annotation index, and the CCF annotation at 20 um
    names, _, _ = structure_terms()
    annotation = annotation_20("ccf")
    coords = {}

    # the annotation indices of each structure name (0 is outside the brain)
    per_name = defaultdict(list)
    for idx in np.unique(annotation):
        if idx == 0 or int(idx) not in names:
            continue
        per_name[names[int(idx)]].append(int(idx))

    # the indices of the structures asked for, each mapped back to its structure
    wanted = {s: per_name.get(s, []) for s in structures}
    flat = {i: s for s, ids in wanted.items() for i in ids}
    if not flat:
        raise ValueError(
            f"none of the {len(structures)} structures is in the CCF annotation, "
            "so they have no centroids"
        )

    # sum the voxel coordinates and count the voxels of each structure; 20 um voxels
    mask = np.isin(annotation, list(flat))
    ap, dv, ml = np.nonzero(mask)
    labels = annotation[mask]
    sums = defaultdict(lambda: np.zeros(4))
    for a, d, m, lab in zip(ap, dv, ml, labels):
        sums[flat[int(lab)]] += (a, d, m, 1)

    # the mean position in mm (a voxel is 0.02 mm), NaN for a structure with no voxels
    for s in structures:
        v = sums.get(s)
        coords[s] = (
            (v[:3] / v[3]) * 0.02 if v is not None and v[3] else np.full(3, np.nan)
        )
    return coords


def replication(
    nano: dict[str, dict[str, float]],
    structures: list[str],
    predictors: Sequence[np.ndarray],
    splits: list[tuple[list[int], list[int]]],
) -> float:
    """How well the leftover of one half-cohort matches the leftover of the other."""
    # each half's ranked map fitted to the predictors on its own, the Spearman of
    # the two leftovers, averaged over every split
    return float(
        np.mean(
            [
                spearmanr(
                    residual(half_map(nano, a, structures), predictors),
                    residual(half_map(nano, b, structures), predictors),
                ).statistic
                for a, b in splits
            ]
        )
    )


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
    y: np.ndarray,
    covariates: dict[str, np.ndarray],
    nano: dict[str, dict[str, float]],
    splits: list[tuple[list[int], list[int]]],
) -> dict[str, str] | None:
    """Control A: whether the leftover is just a smooth gradient across the block.

    Returns the verdict row, or None when fewer than beyond_controls.min_centroids
    structures have a centroid.
    """
    print("\nA  is it a smooth spatial gradient? (an illumination artefact)")

    # the structures with a centroid; on too few, a fit of six position terms would
    # explain much of the leftover by chance
    coords = centroids(structures)
    xyz = np.array([coords[s] for s in structures])
    ok = np.all(np.isfinite(xyz), axis=1)
    if ok.sum() < BEYOND_CONTROLS["min_centroids"]:
        print("   not enough centroids; skipped")
        return None

    # the residual against a quadratic in the three ranked positions, enough to
    # follow a ramp or a bowl across the block, as uneven illumination would give
    pos = [rankdata(xyz[ok, i]) for i in range(3)]
    quad = pos + [p**2 for p in pos]
    smooth = r_squared(res[ok], quad)

    # and against each axis on its own
    for axis, p in zip("AP DV ML".split(), pos):
        print(f"   residual against {axis}: rho {spearmanr(res[ok], p).statistic:+.3f}")
    print(f"   a smooth quadratic in all three axes explains R2 = {smooth:.3f} of it")

    # the replication with position added to the covariates: the bent covariates
    # and the three ranked axes, straight
    with_pos = replication(
        nano,
        structures,
        flexible(list(covariates.values()))
        + [rankdata([coords[s][i] for s in structures]) for i in range(3)],
        splits,
    )
    still = "still replicates" if replicates(with_pos) else "replicates only"
    print(
        f"   and with position added as a covariate the leftover {still} "
        f"at {with_pos:.3f}"
    )

    # the verdict, printed and as a row of controls.csv: a gradient could explain
    # the leftover if position explains more than beyond_controls.gradient_r2 of it
    if smooth > BEYOND_CONTROLS["gradient_r2"]:
        verdict = "a gradient could explain it -- LOOK CLOSER"
    else:
        verdict = "not a gradient; position is a weak predictor of it"
    print("   verdict: " + verdict)
    return dict(
        control="A spatial gradient",
        number=f"smooth R2 {smooth:.3f}, replication with position {with_pos:.3f}",
        verdict="pass" if smooth <= BEYOND_CONTROLS["gradient_r2"] else "CHECK",
    )


def control_b_size(
    res: np.ndarray, structures: list[str], nano_rows: dict[str, float]
) -> dict[str, str]:
    """Control B: whether the leftover comes from small or poorly covered structures.

    `nano_rows` holds each structure's mean volume in 20 um voxels. Returns the
    verdict row.
    """
    print("\nB  is it small structures, where a mean is noisy?")

    # the leftover against structure volume, structures with no volume left out
    size = np.array([nano_rows.get(s, np.nan) for s in structures])
    ok = np.isfinite(size)
    rho = spearmanr(res[ok], np.log10(size[ok])).statistic

    # noise would make the leftover larger in the smaller structures, so compare
    # the median |residual| of the two halves by volume
    big = size[ok] >= np.median(size[ok])
    print(f"   residual against log structure volume: rho {rho:+.3f}")
    print(
        f"   |residual| in the larger half {np.median(np.abs(res[ok][big])):.1f} ranks, "
        f"smaller half {np.median(np.abs(res[ok][~big])):.1f}"
    )

    # the verdict: size drives the leftover if |rho| exceeds beyond_controls.size_rho
    if abs(rho) > BEYOND_CONTROLS["size_rho"]:
        verdict = "size drives it -- LOOK CLOSER"
    else:
        verdict = "size is not what the leftover is made of"
    print("   verdict: " + verdict)
    return dict(
        control="B structure size",
        number=f"rho with volume {rho:+.3f}",
        verdict="pass" if abs(rho) <= BEYOND_CONTROLS["size_rho"] else "CHECK",
    )


def control_c_mice(
    nano: dict[str, dict[str, float]],
    structures: list[str],
    covariates: dict[str, np.ndarray],
) -> tuple[dict[str, str], dict[str, np.ndarray], list[float]]:
    """Control C: whether the leftover is carried by one or two animals.

    Returns the verdict row, each mouse's leftover, and the agreement of every
    pair of mice.
    """
    print("\nC  is it one or two animals?")

    # each mouse's own leftover: its ranked map fitted to the bent covariates
    xs = flexible(list(covariates.values()))
    per = {m: residual(rankdata([nano[m][s] for s in structures]), xs) for m in ADULTS}

    # the agreement of every pair of mice
    pairs = [
        spearmanr(per[a], per[b]).statistic
        for i, a in enumerate(ADULTS)
        for b in ADULTS[i + 1 :]
    ]

    # the mouse whose leftover agrees least, on average, with the others (the
    # first by name on a tie)
    worst_rho, worst = min(
        (
            float(
                np.mean([spearmanr(per[m], per[o]).statistic for o in ADULTS if o != m])
            ),
            m,
        )
        for m in ADULTS
    )
    print(
        f"   each mouse against each other mouse: median rho {np.median(pairs):+.3f}, "
        f"range {min(pairs):+.3f} to {max(pairs):+.3f}"
    )
    print(f"   least typical animal: {worst} at {worst_rho:+.3f} mean agreement")

    # the verdict rests on the worst pair, since an odd brain lowers every pair it
    # is in: below beyond_controls.pair_rho, one animal may carry the leftover
    if min(pairs) < BEYOND_CONTROLS["pair_rho"]:
        verdict = "one animal may be carrying it -- LOOK CLOSER"
    else:
        verdict = "every animal shows the same leftover"
    print("   verdict: " + verdict)
    return (
        dict(
            control="C single animals",
            number=f"pairwise median {np.median(pairs):+.3f}, min {min(pairs):+.3f}",
            verdict="pass" if min(pairs) >= BEYOND_CONTROLS["pair_rho"] else "CHECK",
        ),
        per,
        pairs,
    )


def control_d_groups(
    nano: dict[str, dict[str, float]],
    structures: list[str],
    covariates: dict[str, np.ndarray],
) -> tuple[dict[str, str], np.ndarray, np.ndarray]:
    """Control D: whether the leftover is the whisker manipulation, not the anatomy.

    Returns the verdict row and the leftovers of the naive and the RWS group.
    """
    print("\nD  is it the whisker manipulation? (naive and RWS are pooled)")

    # each group's leftover: its mean map, ranked, fitted to the bent covariates
    xs = flexible(list(covariates.values()))
    naive = residual(
        rankdata([float(np.mean([nano[m][s] for m in NAIVE])) for s in structures]), xs
    )
    rws = residual(
        rankdata([float(np.mean([nano[m][s] for m in RWS])) for s in structures]), xs
    )

    # how well the two groups' leftovers agree; below beyond_controls.groups_rho
    # the pooling hides a difference between them
    rho = spearmanr(naive, rws).statistic
    print(f"   leftover from the five naive against the five RWS: rho {rho:+.3f}")
    if rho < BEYOND_CONTROLS["groups_rho"]:
        verdict = "the groups disagree -- pooling is hiding something"
    else:
        verdict = "both groups give the same leftover, so it is not the manipulation"
    print("   verdict: " + verdict)
    return (
        dict(
            control="D naive vs RWS",
            number=f"rho {rho:+.3f}",
            verdict="pass" if rho >= BEYOND_CONTROLS["groups_rho"] else "CHECK",
        ),
        naive,
        rws,
    )


def control_e_curvature(
    y: np.ndarray, covariates: dict[str, np.ndarray]
) -> tuple[dict[str, str], float, float]:
    """Control E: whether the bending model bends enough.

    As straight lines the covariates reach a cross-validated 0.42, with squares and
    cubes 0.60, which is why adult.beyond_density quotes the bending model
    (`flexible`). If going on to fifth powers still buys prediction, the leftover is
    still partly the model's failing. Returns the verdict row and the
    cross-validated R2 of the cubic and the quintic model.
    """
    print("\nE  is the bending model bent enough? (does more curvature keep paying?)")

    # held-out R2 of the covariates straight, to cubes (the quoted model) and to
    # fifth powers; held out, because in-sample every added term fits better
    xs = list(covariates.values())
    linear = cv_r2(y, xs)
    cubic = cv_r2(y, flexible(xs))
    quintic = cv_r2(y, flexible(xs) + [x**4 for x in xs] + [x**5 for x in xs])
    print(f"   cross-validated R2, straight              {linear:+.3f}")
    print(f"   cross-validated R2, squares and cubes     {cubic:+.3f}   <- the model")
    print(f"   cross-validated R2, up to fifth powers    {quintic:+.3f}")

    # the verdict: the model is not bent enough if fifth powers gain more than
    # beyond_controls.curvature_gain of held-out R2
    gain = BEYOND_CONTROLS["curvature_gain"]
    if quintic - cubic > gain:
        verdict = "more curvature still pays -- the model is not bent enough"
    else:
        verdict = "further bending buys nothing, so the model is adequate"
    print("   verdict: " + verdict)
    return (
        dict(
            control="E curvature",
            number=f"CV R2 {linear:.3f} straight, {cubic:.3f} cubic, "
            f"{quintic:.3f} quintic",
            verdict="pass" if quintic - cubic <= gain else "CHECK",
        ),
        cubic,
        quintic,
    )


def control_f_gene_space(
    nano: dict[str, dict[str, float]],
    expr: dict[str, dict[str, float]],
    structures: list[str],
    y: np.ndarray,
    ceiling: float,
    splits: list[tuple[list[int], list[int]]],
) -> tuple[dict[str, str], list[dict[str, float]], int]:
    """Control F, the strongest: any combination of the panel's genes may try.

    The genes measured in every structure are reduced to principal components;
    models of 1 to beyond_controls.max_pcs components are scored by
    cross-validation, and the leftover of the best is tested for replication.
    Returns the verdict row, the curve (per number of components, n_components, r2
    fitted and cv_r2 cross-validated) and the best number of components.
    """
    # the panel genes measured in every one of the structures
    genes = sorted(g for g in expr if all(s in expr[g] for s in structures))
    print(
        "\nF  is it our choice of covariates? (give the model all "
        f"{len(genes)} genes measured in every structure)"
    )

    # each gene's rank profile z-scored and the mean profile removed; the rows of
    # vt are then the panel's shared patterns across structures, strongest first
    m = gene_matrix(expr, structures, genes)
    m = (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, keepdims=True)
    _, _, vt = np.linalg.svd(m - m.mean(axis=0), full_matrices=False)

    # models of the first 1 to beyond_controls.max_pcs components, each scored
    # in-sample and held out
    curve = []
    for k in range(1, BEYOND_CONTROLS["max_pcs"] + 1):
        pcs = [vt[i] for i in range(k)]
        curve.append(dict(n_components=k, r2=r_squared(y, pcs), cv_r2=cv_r2(y, pcs)))

    # the best model by held-out R2 (in-sample would always pick the most
    # components), and whether its leftover still replicates
    best = max(curve, key=lambda c: c["cv_r2"])
    best_k, best_cv = best["n_components"], best["cv_r2"]
    pcs = [vt[i] for i in range(best_k)]
    rep = replication(nano, structures, pcs, splits)
    print(
        f"   {len(genes)} genes reduced to components; the best model by "
        f"cross-validation uses {best_k}"
    )

    # its held-out R2 as a share of the explainable variance, the squared ceiling
    print(
        f"   in-sample R2 {r_squared(y, pcs):.3f}, cross-validated {best_cv:.3f}, "
        f"{best_cv / ceiling**2:.0%} of the ceiling"
    )
    still = "still replicates" if replicates(rep) else "replicates only"
    print(f"   and the leftover of THAT model {still} at {rep:.3f}")

    # the question is not whether a model this rich explains a lot (with 253 genes
    # it should) but whether it explains the map completely
    if rep < BEYOND_CONTROLS["replication"]:
        verdict = "the gene panel accounts for the map; the leftover is gone"
    else:
        verdict = "even the whole panel leaves a leftover that replicates"
    print("   verdict: " + verdict)
    return (
        dict(
            control="F whole gene space",
            number=f"{best_k} components, CV R2 {best_cv:.3f} "
            f"({best_cv / ceiling**2:.0%} of ceiling), replication {rep:.3f}",
            verdict="pass" if rep >= BEYOND_CONTROLS["replication"] else "CHECK",
        ),
        curve,
        best_k,
    )


def control_g_readings(
    expr: dict[str, dict[str, float]],
    role: dict[str, str],
    auto: dict[str, dict[str, float]],
    structures: list[str],
    splits: list[tuple[list[int], list[int]]],
) -> tuple[dict[str, str], list[dict]]:
    """Control G: whether any of this is specific to zref.

    The same covariates, ceiling and replication for each reading measured in
    every structure. Returns the verdict row and, per reading, a row of the R2 of
    the covariates (r2), the map's replication (map_replication) and the
    leftover's (leftover_replication).
    """
    print("\nG  is it zref? (the same test on every reading)")
    out = []
    for reading in ("zref", "cref", "subref", "ratio", "sepratio"):
        # the reading's log2 value per adult and structure
        per = defaultdict(dict)
        with open(REGION_MEANS, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r["reading"] == reading and r["mouse"] in ADULTS:
                    per[r["mouse"]][r["structure"]] = float(r["log2_value"])

        # the test needs the structures zref used, so a reading missing from any
        # of them in any mouse is skipped
        if not all(all(s in per[m] for s in structures) for m in ADULTS):
            print(f"   {reading:9s} not measured in every structure; skipped")
            continue

        # the bent covariates, and the reading's ranked mean map over the adults
        covariates, _, _ = build_covariates(per, expr, role, auto, structures)
        xs = flexible(list(covariates.values()))
        y = half_map(per, range(len(ADULTS)), structures)

        # how well the two half-cohort maps agree, then their leftovers, averaged
        # over the splits
        raw = float(
            np.mean(
                [
                    spearmanr(
                        half_map(per, a, structures), half_map(per, b, structures)
                    ).statistic
                    for a, b in splits
                ]
            )
        )
        rep = replication(per, structures, xs, splits)

        # one row per reading, R2 in-sample
        out.append(
            dict(
                reading=reading,
                r2=r_squared(y, xs),
                map_replication=raw,
                leftover_replication=rep,
            )
        )
        print(
            f"   {reading:9s} covariates explain R2 {out[-1]['r2']:.3f}; "
            f"map replicates {raw:.3f}, leftover {rep:.3f}"
        )

    # the verdict: the claim does not hang on zref if every reading's leftover
    # replicates above beyond_controls.readings_replication
    agree = all(
        r["leftover_replication"] > BEYOND_CONTROLS["readings_replication"] for r in out
    )
    if agree:
        verdict = "the leftover replicates under every reading"
    else:
        verdict = "some readings disagree -- LOOK CLOSER"
    print("   verdict: " + verdict)
    return (
        dict(
            control="G reading choice",
            number="; ".join(
                f"{r['reading']} {r['leftover_replication']:.2f}" for r in out
            ),
            verdict="pass" if agree else "CHECK",
        ),
        out,
    )


# ===== Figures =====


def figure_artefacts(
    res: np.ndarray,
    structures: list[str],
    sizes: np.ndarray,
    per_mouse: dict[str, np.ndarray],
    pairs: list[float],
    naive: np.ndarray,
    rws: np.ndarray,
    coords: dict[str, np.ndarray],
    passed: dict[str, bool],
) -> None:
    """Draw controls A to D: position, size, pairs of mice, naive against RWS.

    `passed` holds each control's verdict by letter (A is missing when skipped);
    each title says what its verdict says.
    """
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 3.9))

    # A: the leftover against anterior-posterior position, structures with a centroid
    xyz = np.array([coords[s] for s in structures])
    ok = np.all(np.isfinite(xyz), axis=1)
    axes[0].scatter(
        xyz[ok, 0], res[ok], s=12, facecolor="0.6", edgecolor="0.25", linewidth=0.3
    )
    axes[0].axhline(0, color="0.85", lw=0.7)
    axes[0].set_xlabel("structure centroid, anterior-posterior (mm)", fontsize=8)
    axes[0].set_ylabel("residual (ranks)", fontsize=8)

    # each title says what the control's verdict says
    if "A" not in passed:
        title = "A. not tested: too few centroids"
    elif passed["A"]:
        title = "A. not a front-to-back gradient"
    else:
        title = "A. a gradient could explain it"
    axes[0].set_title(title, fontsize=9)
    tidy(axes[0])

    # B: the leftover against structure volume, structures with a volume
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

    # C: the agreement of every pair of mice, the median in red
    axes[2].hist(pairs, bins=20, color="0.7", edgecolor="0.35", linewidth=0.4)
    axes[2].axvline(float(np.median(pairs)), color=RED, lw=1.8)
    axes[2].set_xlabel("leftover of one mouse against another (Spearman)", fontsize=8)
    axes[2].set_ylabel("pairs of animals", fontsize=8)
    title = "C. every animal shows it" if passed["C"] else "C. one animal may carry it"
    axes[2].set_title(title, fontsize=9)
    tidy(axes[2])

    # D: the naive leftover against the RWS one, with the identity line over both
    # ranges and a margin of three ranks
    axes[3].scatter(naive, rws, s=12, facecolor="0.6", edgecolor="0.25", linewidth=0.3)
    lim = [min(naive.min(), rws.min()) - 3, max(naive.max(), rws.max()) + 3]
    axes[3].plot(lim, lim, color="0.75", ls="--", lw=0.8)
    axes[3].set_xlabel("leftover, five naive animals", fontsize=8)
    axes[3].set_ylabel("leftover, five RWS animals", fontsize=8)
    title = (
        "D. not the whisker manipulation" if passed["D"] else "D. naive and RWS disagree"
    )
    axes[3].set_title(
        f"{title}\nrho {spearmanr(naive, rws).statistic:+.2f}",
        fontsize=9,
    )
    tidy(axes[3])

    # the overall title names the controls that did not rule their artefact out
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
    curve: list[dict[str, float]],
    best_k: int,
    ceiling: float,
    cubic: float,
    quintic: float,
    passed: dict[str, bool],
) -> None:
    """Draw controls F and E: the gene-space curve, and bending further.

    `passed` holds each control's verdict by letter; each title follows its own.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    # F: in-sample and held-out R2 against the number of components
    ks = [c["n_components"] for c in curve]
    axes[0].plot(ks, [c["r2"] for c in curve], color="0.6", lw=1.5, label="fitted")
    axes[0].plot(
        ks, [c["cv_r2"] for c in curve], color=RED, lw=1.8, label="cross-validated"
    )

    # the explainable variance (the squared ceiling), and the best model dotted
    axes[0].axhline(ceiling**2, color="0.3", ls="--", lw=1.2)
    axes[0].annotate(
        "ceiling", (ks[-1], ceiling**2), fontsize=7.5, ha="right", va="bottom"
    )
    axes[0].axvline(best_k, color="0.4", ls=":", lw=1.0)
    axes[0].set_xlabel(
        "components of the expression of the genes measured in every structure",
        fontsize=8,
    )
    axes[0].set_ylabel("variance of the map explained", fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False)

    # the title follows the verdict of F
    if passed["F"]:
        title = "F. even the whole panel falls short"
    else:
        title = "F. the whole panel accounts for the map"
    axes[0].set_title(
        f"{title}\nthe gap between the two lines is overfitting", fontsize=9
    )
    tidy(axes[0])

    # E: held-out R2 of the quoted model and of the model bent to fifth powers,
    # titled by the verdict of E
    axes[1].bar(
        [0, 1],
        [cubic, quintic],
        color=[RED, "0.65"],
        edgecolor="0.25",
        linewidth=0.5,
    )
    axes[1].set_xticks([0, 1])
    axes[1].set_xticklabels(
        ["the model\n(squares, cubes)", "bent further\n(to fifth powers)"], fontsize=8
    )
    axes[1].set_ylabel("cross-validated R2", fontsize=8)
    if passed["E"]:
        title = "E. and bending it further buys nothing"
    else:
        title = "E. and bending it further still pays"
    axes[1].set_title(title, fontsize=9)
    tidy(axes[1])

    fig.tight_layout()
    save(fig, "fig5_model_space.png")


def figure_readings(rows: list[dict], passed: bool) -> None:
    """Draw control G: per reading, the map's and the leftover's replication, and R2.

    The title follows the control's verdict, `passed`.
    """
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    labels = [r["reading"] for r in rows]
    x = np.arange(len(rows))

    # left of each reading: how well the two half-cohort maps agree
    ax.bar(
        x - 0.22,
        [r["map_replication"] for r in rows],
        width=0.2,
        color="0.55",
        edgecolor="0.25",
        linewidth=0.4,
        label="the map replicates",
    )

    # middle, in red: how well their leftovers agree, the claim itself
    ax.bar(
        x,
        [r["leftover_replication"] for r in rows],
        width=0.2,
        color=RED,
        edgecolor="0.25",
        linewidth=0.4,
        label="the leftover replicates",
    )

    # right: the in-sample R2 of the bent covariates
    ax.bar(
        x + 0.22,
        [r["r2"] for r in rows],
        width=0.2,
        color="0.8",
        edgecolor="0.25",
        linewidth=0.4,
        label="covariates explain (R2)",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")

    # the title follows the verdict of G
    if passed:
        title = "G. the same picture under every reading, not just zref"
    else:
        title = "G. not the same picture under every reading"
    ax.set_title(title, fontsize=9)
    tidy(ax)
    fig.tight_layout()
    save(fig, "fig6_readings.png")


def map_ceiling(
    nano: dict[str, dict[str, float]],
    structures: list[str],
    splits: list[tuple[list[int], list[int]]],
) -> float:
    """The map's reliability: Spearman-Brown of the mean half-cohort agreement."""
    # the two half-cohort maps' Spearman averaged over every split, then stepped up
    # to what the whole cohort is worth
    ceiling = spearman_brown(
        float(
            np.mean(
                [
                    spearmanr(
                        half_map(nano, a, structures), half_map(nano, b, structures)
                    ).statistic
                    for a, b in splits
                ]
            )
        )
    )
    return ceiling


def mean_sizes() -> dict[str, float]:
    """Mean structure volume over the adults, in 20 um voxels."""
    # each adult's voxel count per structure, from the zref rows, one per mouse
    # and structure (the table holds a row per reading)
    sizes = {}
    with open(REGION_MEANS, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["reading"] == "zref" and r["mouse"] in ADULTS:
                sizes.setdefault(r["structure"], []).append(float(r["n_vox20"]))
    size_mean = {k: float(np.mean(v)) for k, v in sizes.items()}
    return size_mean


def write_verdicts(verdicts: list[dict[str, str] | None]) -> None:
    """Write controls.csv, a skipped control (None) left out, and print the verdict."""
    # one row per control that ran
    path = OUT / "controls.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["control", "number", "verdict"])
        writer.writeheader()
        writer.writerows([v for v in verdicts if v])
    print(f"\n-> {path}")

    # the claim stands only if every control that ran passes
    failed = [v["control"] for v in verdicts if v and v["verdict"] != "pass"]
    if not failed:
        verdict = "every control passes; the claim stands as written."
    else:
        verdict = f"look closer at {failed}"
    print("VERDICT: " + verdict)


def main() -> None:
    """Run the seven controls, write their verdicts and draw them."""
    # the structures and the quoted model, as adult.beyond_density builds them
    nano, _, expr, role, auto, structures = prepare()

    # the ranked mean map of the ten adults and its leftover after the bent covariates
    covariates, _, _ = build_covariates(nano, expr, role, auto, structures)
    splits = half_splits()
    y = half_map(nano, range(len(ADULTS)), structures)
    res = residual(y, flexible(list(covariates.values())))

    # the ceiling and the quoted model's replication, as beyond_density finds them:
    # the baseline the controls are read against
    ceiling = map_ceiling(nano, structures, splits)
    quoted_replication = replication(
        nano, structures, flexible(list(covariates.values())), splits
    )
    print(
        f"{len(structures)} structures, ceiling {ceiling:.3f}, leftover of the "
        f"quoted model replicates at {quoted_replication:.3f}"
    )

    # mean structure volume over the adults, in 20 um voxels
    size_mean = mean_sizes()

    # the seven controls, each returning its verdict row (None when skipped) and
    # what its figure draws; the centroids are for the figure of A to D
    verdicts = []
    coords = centroids(structures)
    verdicts.append(control_a_space(res, structures, y, covariates, nano, splits))
    verdicts.append(control_b_size(res, structures, size_mean))
    vc, per_mouse, pairs = control_c_mice(nano, structures, covariates)
    verdicts.append(vc)
    vd, naive, rws = control_d_groups(nano, structures, covariates)
    verdicts.append(vd)
    ve, cubic, quintic = control_e_curvature(y, covariates)
    verdicts.append(ve)
    vf, curve, best_k = control_f_gene_space(nano, expr, structures, y, ceiling, splits)
    verdicts.append(vf)
    vg, reading_rows = control_g_readings(expr, role, auto, structures, splits)
    verdicts.append(vg)

    # the verdicts, and the figures
    write_verdicts(verdicts)
    sizes_arr = np.array([size_mean.get(s, np.nan) for s in structures])

    # each control's verdict by its letter, the first of its name; a skipped one
    # is missing
    passed = {v["control"][0]: v["verdict"] == "pass" for v in verdicts if v}
    figure_artefacts(
        res, structures, sizes_arr, per_mouse, pairs, naive, rws, coords, passed
    )
    figure_model_space(curve, best_k, ceiling, cubic, quintic, passed)
    figure_readings(reading_rows, passed["G"])
