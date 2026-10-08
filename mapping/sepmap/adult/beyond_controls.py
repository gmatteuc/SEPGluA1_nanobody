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
                                      components, and see whether a leftover
                                      survives. It is why the leftover is "not
                                      predicted by receptor mRNA or synaptic
                                      density" and never "beyond gene expression"
    G  zref?                          the same test on the other readings of the
                                      map, and on zref with the 17-brain reference

E and F are cross-validated, because a flexible model always fits better on the
data it was fitted to; the question is whether it predicts better, and only
held-out structures can say. Everything is imported from adult.beyond_density, so
the two modules use the same structures, predictors and arithmetic.

Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    controls.csv            one row per control, with its number and verdict
    gene_space.csv          control F: per number of components, R2 fitted and held out
    readings.csv            control G: per reading, R2 and the two replications
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

from sepmap.adult import profiles
from sepmap.adult.beyond_density import (
    ADULTS,
    NAIVE_ROWS,
    OUT,
    RWS_ROWS,
    Inputs,
    build_covariates,
    ceiling,
    cv_r2,
    flexible,
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
from sepmap.plotting import RED, tidy
from sepmap.structures import load_centroids
from sepmap.young_vs_adult.region_plot import REGION_MEANS

# the largest gene-space model of control F, and the threshold of each verdict
BEYOND_CONTROLS = SETTINGS["beyond_controls"]

CONTROLS = OUT / "controls.csv"
GENE_SPACE = OUT / "gene_space.csv"
READINGS = OUT / "readings.csv"

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


def control_f_gene_space(
    inputs: Inputs,
    y: np.ndarray,
    explainable: float,
    splits: list[tuple[list[int], list[int]]],
) -> tuple[dict[str, str], pd.DataFrame, int]:
    """Control F, the strongest: any combination of the gene table's genes may try.

    The genes measured in every structure are reduced to principal components;
    models of 1 to beyond_controls.max_pcs components are scored by cross-validation,
    and the leftover of the best is tested for replication. Returns the verdict
    row, the curve (n_components, r2 fitted, cv_r2 held out, share_of_ceiling) and
    the best number of components.
    """
    s = inputs.structures
    genes = sorted(g for g in inputs.expr if all(x in inputs.expr[g] for x in s))
    print(
        "\nF  is it the choice of predictors? (give the model all "
        f"{len(genes)} genes measured in every structure)"
    )

    # each gene's rank profile z-scored and the mean profile removed; the rows of
    # vt are then the panel's shared patterns across structures, strongest first
    m = gene_matrix(inputs.expr, s, genes)
    m = (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, keepdims=True)
    _, _, vt = np.linalg.svd(m - m.mean(axis=0), full_matrices=False)
    rows = []
    for k in range(1, BEYOND_CONTROLS["max_pcs"] + 1):
        pcs = [vt[i] for i in range(k)]
        cv = cv_r2(y, pcs)
        rows.append(
            dict(
                n_components=k,
                r2=r_squared(y, pcs),
                cv_r2=cv,
                share_of_ceiling=cv / explainable,
            )
        )
    curve = pd.DataFrame(rows)

    # the best model by held-out R2 (in-sample would always pick the most
    # components), and whether its leftover still replicates
    best = curve.loc[curve["cv_r2"].idxmax()]
    best_k = int(best["n_components"])
    rep = replication(inputs.nano, [vt[i] for i in range(best_k)], splits)
    print(
        f"   the best model by cross-validation uses {best_k} components: CV R2 "
        f"{best['cv_r2']:.3f}, {best['share_of_ceiling']:.0%} of the ceiling"
    )
    still = "still replicates" if replicates(rep) else "replicates only"
    print(f"   and the leftover of that model {still} at {rep:.3f}")

    # the question is not whether a model this rich predicts a lot (it should) but
    # whether it predicts the map completely
    passed = rep >= BEYOND_CONTROLS["replication"]
    verdict = (
        "even the whole gene table leaves a leftover that replicates"
        if passed
        else "the gene table accounts for the map; the leftover is gone"
    )
    print("   verdict: " + verdict)
    return (
        dict(
            control="F whole gene space",
            number=f"{len(genes)} genes, {best_k} components, CV R2 "
            f"{best['cv_r2']:.3f} ({best['share_of_ceiling']:.0%} of ceiling), "
            f"replication {rep:.3f}",
            verdict="pass" if passed else "CHECK",
        ),
        curve,
        best_k,
    )


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
    vf, curve, best_k = control_f_gene_space(inputs, y, explainable, splits)
    vg, readings = control_g_readings(inputs, xs, splits)
    verdicts += [vc, vd, ve, vf, vg]
    write_verdicts(verdicts)
    curve.to_csv(GENE_SPACE, index=False)
    readings.to_csv(READINGS, index=False)

    passed = {v["control"][0]: v["verdict"] == "pass" for v in verdicts if v}
    figure_artefacts(res, s, sizes, pairs, naive, rws, coords, passed)
    figure_model_space(curve, best_k, explainable, cubic, quintic, passed)
    figure_readings(readings, passed["G"])
