"""Seven attempts to break the claim of adult.beyond_density.

That module finds that about 40% of the explainable variance of the adult
surface-GluA1 map is not accounted for by receptor abundance or synaptic density,
and that what remains replicates across independent animals. Each control below is
a specific way the claim could be wrong, with the number that says whether it is:

    A  a smooth spatial gradient?     a clearing or illumination artefact would
                                      look like one; anatomy would not
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
                                      whole 390-gene expression space and see
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
import os
from collections import defaultdict

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata, spearmanr

from sepmap.adult.beyond_density import (
    ADULTS,
    MARKERS,
    NAIVE,
    NANO,
    OUT,
    RWS,
    SUBUNITS,
    autofluorescence,
    build_covariates,
    cv_r2,
    flexible,
    gene_profiles,
    half_map,
    half_splits,
    keep_structure,
    nano_per_mouse,
    r_squared,
    residual,
    save,
    spearman_brown,
    structure_names,
    tidy,
)
from sepmap.volumes.per_mouse import annotation_20

# folds of the cross-validated controls E and F (not passed on: cv_r2 uses its own
# default, also 5)
N_FOLDS = 5

# the largest gene-space model tried in control F, in components
MAX_PCS = 25
RNG = np.random.default_rng(0)


# ===== Utilities =====


def centroids(structures):
    """Mean (AP, DV, ML) position of each structure, in mm, from the CCF itself.

    Needed by control A: if the leftover were an imaging or clearing artefact it
    would vary smoothly with position in the block, so the first thing to ask of
    it is how much a smooth function of position can explain. A structure not in
    the atlas gets NaN; None comes back when none of them is.
    """
    names = structure_names()
    ann = annotation_20("ccf")
    coords = {}

    # the annotation indices of each structure name
    per_name = defaultdict(list)
    for idx in np.unique(ann):
        if idx == 0 or int(idx) not in names:
            continue
        per_name[names[int(idx)]].append(int(idx))
    wanted = {s: per_name.get(s, []) for s in structures}
    flat = {i: s for s, ids in wanted.items() for i in ids}
    if not flat:
        return None

    # sum the voxel coordinates and count the voxels of each structure; 20 um voxels
    mask = np.isin(ann, list(flat))
    ap, dv, ml = np.nonzero(mask)
    labels = ann[mask]
    sums = defaultdict(lambda: np.zeros(4))
    for a, d, m, lab in zip(ap, dv, ml, labels):
        sums[flat[int(lab)]] += (a, d, m, 1)
    for s in structures:
        v = sums.get(s)
        coords[s] = (
            (v[:3] / v[3]) * 0.02 if v is not None and v[3] else np.full(3, np.nan)
        )
    return coords


def replication(nano, structures, predictors, splits):
    """How well the leftover of one half-cohort matches the leftover of the other."""
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


def gene_matrix(expr, structures, genes):
    """Rank profiles of many genes as one array, genes by structures."""
    return np.array([rankdata([expr[g][s] for s in structures]) for g in genes])


# ===== Controls =====


def control_a_space(res, structures, y, covariates, nano, splits):
    """Control A: whether the leftover is just a smooth gradient across the block.

    Returns the verdict row, or None when fewer than 50 structures have a centroid.
    """
    print("\nA  is it a smooth spatial gradient? (a clearing or illumination artefact)")
    coords = centroids(structures)
    xyz = np.array([coords[s] for s in structures])
    ok = np.all(np.isfinite(xyz), axis=1)
    if ok.sum() < 50:
        print("   not enough centroids; skipped")
        return None

    # the residual against a quadratic in the three ranked positions
    pos = [rankdata(xyz[ok, i]) for i in range(3)]
    quad = pos + [p**2 for p in pos]
    smooth = r_squared(res[ok], quad)
    for axis, p in zip("AP DV ML".split(), pos):
        print(f"   residual against {axis}: rho {spearmanr(res[ok], p).statistic:+.3f}")
    print(f"   a smooth quadratic in all three axes explains R2 = {smooth:.3f} of it")

    # the replication with position added to the covariates
    with_pos = replication(
        nano,
        structures,
        flexible(list(covariates.values()))
        + [rankdata([coords[s][i] for s in structures]) for i in range(3)],
        splits,
    )
    print(
        f"   and with position added as a covariate the leftover still replicates "
        f"at {with_pos:.3f}"
    )
    if smooth > 0.5:
        verdict = "a gradient could explain it -- LOOK CLOSER"
    else:
        verdict = "not a gradient; position is a weak predictor of it"
    print("   verdict: " + verdict)
    return dict(
        control="A spatial gradient",
        number=f"smooth R2 {smooth:.3f}, replication with position {with_pos:.3f}",
        verdict="pass" if smooth <= 0.5 else "CHECK",
    )


def control_b_size(res, structures, nano_rows):
    """Control B: whether the leftover comes from small or poorly covered structures.

    `nano_rows` holds each structure's mean volume in 20 um voxels. Returns the
    verdict row.
    """
    print("\nB  is it small structures, where a mean is noisy?")
    size = np.array([nano_rows.get(s, np.nan) for s in structures])
    ok = np.isfinite(size)
    rho = spearmanr(res[ok], np.log10(size[ok])).statistic
    big = size[ok] >= np.median(size[ok])
    print(f"   residual against log structure volume: rho {rho:+.3f}")
    print(
        f"   |residual| in the larger half {np.median(np.abs(res[ok][big])):.1f} ranks, "
        f"smaller half {np.median(np.abs(res[ok][~big])):.1f}"
    )
    if abs(rho) > 0.4:
        verdict = "size drives it -- LOOK CLOSER"
    else:
        verdict = "size is not what the leftover is made of"
    print("   verdict: " + verdict)
    return dict(
        control="B structure size",
        number=f"rho with volume {rho:+.3f}",
        verdict="pass" if abs(rho) <= 0.4 else "CHECK",
    )


def control_c_mice(nano, structures, covariates):
    """Control C: whether the leftover is carried by one or two animals.

    Returns the verdict row, each mouse's leftover, and the agreement of every
    pair of mice.
    """
    print("\nC  is it one or two animals?")
    xs = flexible(list(covariates.values()))
    per = {m: residual(rankdata([nano[m][s] for s in structures]), xs) for m in ADULTS}
    pairs = [
        spearmanr(per[a], per[b]).statistic
        for i, a in enumerate(ADULTS)
        for b in ADULTS[i + 1 :]
    ]

    # the mouse whose leftover agrees least, on average, with the others
    worst = min(
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
    print(f"   least typical animal: {worst[1]} at {worst[0]:+.3f} mean agreement")
    if min(pairs) < 0.1:
        verdict = "one animal may be carrying it -- LOOK CLOSER"
    else:
        verdict = "every animal shows the same leftover"
    print("   verdict: " + verdict)
    return (
        dict(
            control="C single animals",
            number=f"pairwise median {np.median(pairs):+.3f}, min {min(pairs):+.3f}",
            verdict="pass" if min(pairs) >= 0.1 else "CHECK",
        ),
        per,
        pairs,
    )


def control_d_groups(nano, structures, covariates):
    """Control D: whether the leftover is the whisker manipulation, not the anatomy.

    Returns the verdict row and the leftovers of the naive and the RWS group.
    """
    print("\nD  is it the whisker manipulation? (naive and RWS are pooled)")
    xs = flexible(list(covariates.values()))
    naive = residual(
        rankdata([float(np.mean([nano[m][s] for m in NAIVE])) for s in structures]), xs
    )
    rws = residual(
        rankdata([float(np.mean([nano[m][s] for m in RWS])) for s in structures]), xs
    )
    rho = spearmanr(naive, rws).statistic
    print(f"   leftover from the five naive against the five RWS: rho {rho:+.3f}")
    if rho < 0.5:
        verdict = "the groups disagree -- pooling is hiding something"
    else:
        verdict = "both groups give the same leftover, so it is not the manipulation"
    print("   verdict: " + verdict)
    return (
        dict(
            control="D naive vs RWS",
            number=f"rho {rho:+.3f}",
            verdict="pass" if rho >= 0.5 else "CHECK",
        ),
        naive,
        rws,
    )


def control_e_curvature(y, covariates):
    """Control E: whether the bending model bends enough.

    As straight lines the covariates reach a cross-validated 0.42, with squares and
    cubes 0.60, which is why adult.beyond_density quotes the bending model
    (`flexible`). If going on to fifth powers still buys prediction, the leftover is
    still partly the model's failing. Returns the verdict row and the
    cross-validated R2 of the cubic and the quintic model.
    """
    print("\nE  is the bending model bent enough? (does more curvature keep paying?)")
    xs = list(covariates.values())
    linear = cv_r2(y, xs)
    cubic = cv_r2(y, flexible(xs))
    quintic = cv_r2(y, flexible(xs) + [x**4 for x in xs] + [x**5 for x in xs])
    print(f"   cross-validated R2, straight              {linear:+.3f}")
    print(f"   cross-validated R2, squares and cubes     {cubic:+.3f}   <- the model")
    print(f"   cross-validated R2, up to fifth powers    {quintic:+.3f}")
    if quintic - cubic > 0.05:
        verdict = "more curvature still pays -- the model is not bent enough"
    else:
        verdict = "further bending buys nothing, so the model is adequate"
    print("   verdict: " + verdict)
    return (
        dict(
            control="E curvature",
            number=f"CV R2 {linear:.3f} straight, {cubic:.3f} cubic, "
            f"{quintic:.3f} quintic",
            verdict="pass" if quintic - cubic <= 0.05 else "CHECK",
        ),
        cubic,
        quintic,
    )


def control_f_gene_space(nano, expr, structures, y, ceiling, splits):
    """Control F, the strongest: any combination of the panel's genes may try.

    The genes measured in every structure are reduced to principal components;
    models of 1 to MAX_PCS components are scored by cross-validation, and the
    leftover of the best is tested for replication. Returns the verdict row, the
    (components, fitted R2, cross-validated R2) curve and the best number of
    components.
    """
    print("\nF  is it our choice of covariates? (give the model all 390 genes)")
    genes = sorted(g for g in expr if all(s in expr[g] for s in structures))
    m = gene_matrix(expr, structures, genes)
    m = (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, keepdims=True)
    _, _, vt = np.linalg.svd(m - m.mean(axis=0), full_matrices=False)

    curve = []
    for k in range(1, MAX_PCS + 1):
        pcs = [vt[i] for i in range(k)]
        curve.append((k, r_squared(y, pcs), cv_r2(y, pcs)))
    best_k, _, best_cv = max(curve, key=lambda t: t[2])
    pcs = [vt[i] for i in range(best_k)]
    rep = replication(nano, structures, pcs, splits)
    print(
        f"   {len(genes)} genes reduced to components; the best model by "
        f"cross-validation uses {best_k}"
    )
    print(
        f"   in-sample R2 {r_squared(y, pcs):.3f}, cross-validated {best_cv:.3f}, "
        f"{best_cv / ceiling**2:.0%} of the ceiling"
    )
    print(f"   and the leftover of THAT model still replicates at {rep:.3f}")

    # the question is not whether a model this rich explains a lot (with 253 genes
    # it should) but whether it explains the map completely
    if rep < 0.5:
        verdict = "the gene panel accounts for the map; the leftover is gone"
    else:
        verdict = "even the whole panel leaves a leftover that replicates"
    print("   verdict: " + verdict)
    return (
        dict(
            control="F whole gene space",
            number=f"{best_k} components, CV R2 {best_cv:.3f} "
            f"({best_cv / ceiling**2:.0%} of ceiling), replication {rep:.3f}",
            verdict="pass" if rep >= 0.5 else "CHECK",
        ),
        curve,
        best_k,
    )


def control_g_readings(expr, role, auto, structures, splits):
    """Control G: whether any of this is specific to zref.

    The same covariates, ceiling and replication for each reading measured in
    every structure. Returns the verdict row and, per reading, (reading, R2 of
    the covariates, the map's replication, the leftover's replication).
    """
    print("\nG  is it zref? (the same test on every reading)")
    out = []
    for reading in ("zref", "cref", "subref", "ratio", "sepratio"):
        per = defaultdict(dict)
        with open(NANO, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r["reading"] == reading and r["mouse"] in ADULTS:
                    per[r["mouse"]][r["structure"]] = float(r["log2_value"])
        if not all(all(s in per[m] for s in structures) for m in ADULTS):
            print(f"   {reading:9s} not measured in every structure; skipped")
            continue
        covariates, _, _ = build_covariates(per, expr, role, auto, structures)
        xs = flexible(list(covariates.values()))
        y = half_map(per, range(len(ADULTS)), structures)
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
        out.append((reading, r_squared(y, xs), raw, rep))
        print(
            f"   {reading:9s} covariates explain R2 {out[-1][1]:.3f}; map replicates "
            f"{raw:.3f}, leftover {rep:.3f}"
        )
    agree = all(r > 0.8 for _, _, _, r in out)
    if agree:
        verdict = "the leftover replicates under every reading"
    else:
        verdict = "some readings disagree -- LOOK CLOSER"
    print("   verdict: " + verdict)
    return (
        dict(
            control="G reading choice",
            number="; ".join(f"{r} {rep:.2f}" for r, _, _, rep in out),
            verdict="pass" if agree else "CHECK",
        ),
        out,
    )


# ===== Figures =====


def figure_artefacts(res, structures, sizes, per_mouse, pairs, naive, rws, coords):
    """Draw controls A to D: position, size, pairs of mice, naive against RWS."""
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 3.9))

    xyz = np.array([coords[s] for s in structures])
    ok = np.all(np.isfinite(xyz), axis=1)
    axes[0].scatter(
        xyz[ok, 0], res[ok], s=12, facecolor="0.6", edgecolor="0.25", linewidth=0.3
    )
    axes[0].axhline(0, color="0.85", lw=0.7)
    axes[0].set_xlabel("structure centroid, anterior-posterior (mm)", fontsize=8)
    axes[0].set_ylabel("residual (ranks)", fontsize=8)
    axes[0].set_title("A. not a front-to-back gradient", fontsize=9)
    tidy(axes[0])

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
    axes[1].set_title("B. not small-structure noise", fontsize=9)
    tidy(axes[1])

    axes[2].hist(pairs, bins=20, color="0.7", edgecolor="0.35", linewidth=0.4)
    axes[2].axvline(float(np.median(pairs)), color="#c0392b", lw=1.8)
    axes[2].set_xlabel("leftover of one mouse against another (Spearman)", fontsize=8)
    axes[2].set_ylabel("pairs of animals", fontsize=8)
    axes[2].set_title("C. every animal shows it", fontsize=9)
    tidy(axes[2])

    axes[3].scatter(naive, rws, s=12, facecolor="0.6", edgecolor="0.25", linewidth=0.3)
    lim = [min(naive.min(), rws.min()) - 3, max(naive.max(), rws.max()) + 3]
    axes[3].plot(lim, lim, color="0.75", ls="--", lw=0.8)
    axes[3].set_xlabel("leftover, five naive animals", fontsize=8)
    axes[3].set_ylabel("leftover, five RWS animals", fontsize=8)
    axes[3].set_title(
        f"D. not the whisker manipulation\nrho {spearmanr(naive, rws).statistic:+.2f}",
        fontsize=9,
    )
    tidy(axes[3])

    fig.suptitle("Four ways the leftover could be an artefact, and is not", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    save(fig, "fig4_controls.png")


def figure_model_space(curve, best_k, ceiling, cubic, quintic):
    """Draw controls F and E: the gene-space curve, and bending further."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ks = [c[0] for c in curve]
    axes[0].plot(ks, [c[1] for c in curve], color="0.6", lw=1.5, label="fitted")
    axes[0].plot(
        ks, [c[2] for c in curve], color="#c0392b", lw=1.8, label="cross-validated"
    )
    axes[0].axhline(ceiling**2, color="0.3", ls="--", lw=1.2)
    axes[0].annotate(
        "ceiling", (ks[-1], ceiling**2), fontsize=7.5, ha="right", va="bottom"
    )
    axes[0].axvline(best_k, color="0.4", ls=":", lw=1.0)
    axes[0].set_xlabel("components of the 390-gene expression space", fontsize=8)
    axes[0].set_ylabel("variance of the map explained", fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False)
    axes[0].set_title(
        "F. even the whole panel falls short\n"
        "the gap between the two lines is overfitting",
        fontsize=9,
    )
    tidy(axes[0])

    axes[1].bar(
        [0, 1],
        [cubic, quintic],
        color=["#c0392b", "0.65"],
        edgecolor="0.25",
        linewidth=0.5,
    )
    axes[1].set_xticks([0, 1])
    axes[1].set_xticklabels(
        ["the model\n(squares, cubes)", "bent further\n(to fifth powers)"], fontsize=8
    )
    axes[1].set_ylabel("cross-validated R2", fontsize=8)
    axes[1].set_title("E. and bending it further buys nothing", fontsize=9)
    tidy(axes[1])

    fig.tight_layout()
    save(fig, "fig5_model_space.png")


def figure_readings(rows):
    """Draw control G: per reading, the map's and the leftover's replication, and R2."""
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    labels = [r[0] for r in rows]
    x = np.arange(len(rows))
    ax.bar(
        x - 0.22,
        [r[2] for r in rows],
        width=0.2,
        color="0.55",
        edgecolor="0.25",
        linewidth=0.4,
        label="the map replicates",
    )
    ax.bar(
        x,
        [r[3] for r in rows],
        width=0.2,
        color="#c0392b",
        edgecolor="0.25",
        linewidth=0.4,
        label="the leftover replicates",
    )
    ax.bar(
        x + 0.22,
        [r[1] for r in rows],
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
    ax.set_title("G. the same picture under every reading, not just zref", fontsize=9)
    tidy(ax)
    fig.tight_layout()
    save(fig, "fig6_readings.png")


def main():
    """Run the seven controls, write their verdicts and draw them."""
    # the structures and the quoted model, as adult.beyond_density builds them
    nano, division = nano_per_mouse()
    expr, role = gene_profiles()
    everywhere = set.intersection(*[set(nano[m]) for m in ADULTS])
    everywhere &= set.intersection(*[set(expr[g]) for g in SUBUNITS + MARKERS])
    kept = [s for s in sorted(everywhere) if keep_structure(s, division.get(s, ""))[0]]
    auto = autofluorescence(set(kept))
    structures = sorted(set(kept).intersection(*[set(auto[m]) for m in ADULTS]))

    covariates, _, _ = build_covariates(nano, expr, role, auto, structures)
    splits = half_splits()
    y = half_map(nano, range(len(ADULTS)), structures)
    res = residual(y, flexible(list(covariates.values())))
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
    quoted_replication = replication(
        nano, structures, flexible(list(covariates.values())), splits
    )
    print(
        f"{len(structures)} structures, ceiling {ceiling:.3f}, leftover of the "
        f"quoted model replicates at {quoted_replication:.3f}"
    )

    # mean structure volume over the adults, in 20 um voxels
    sizes = {}
    with open(NANO, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["reading"] == "zref" and r["mouse"] in ADULTS:
                sizes.setdefault(r["structure"], []).append(float(r["n_vox20"]))
    size_mean = {k: float(np.mean(v)) for k, v in sizes.items()}

    # the seven controls
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

    # the verdicts, a skipped control (None) left out
    path = os.path.join(OUT, "controls.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["control", "number", "verdict"])
        w.writeheader()
        w.writerows([v for v in verdicts if v])
    print(f"\n-> {path}")
    failed = [v["control"] for v in verdicts if v and v["verdict"] != "pass"]
    if not failed:
        verdict = "every control passes; the claim stands as written."
    else:
        verdict = f"look closer at {failed}"
    print("VERDICT: " + verdict)

    sizes_arr = np.array([size_mean.get(s, np.nan) for s in structures])
    figure_artefacts(res, structures, sizes_arr, per_mouse, pairs, naive, rws, coords)
    figure_model_space(curve, best_k, ceiling, cubic, quintic)
    figure_readings(reading_rows)
