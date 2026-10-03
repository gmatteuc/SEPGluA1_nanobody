"""The beyond-abundance regression shown as a regression: observed, predicted, residual.

adult.beyond_density and adult.beyond_figures report the regression's summary
numbers (R2, replication, controls); none of them shows the fit. This module does,
three ways:

    1  observed against predicted, one dot per structure, with the identity line.
       A good model puts the cloud on that line; the spread away from it is the
       leftover, drawn rather than summarised.
    2  residual against predicted, the standard diagnostic. A tilt or a fan here
       would mean the model is mis-specified rather than merely incomplete.
    3  the same three quantities painted back onto the brain (observed, predicted
       and residual maps), because the leftover is a spatial claim.

What is regressed on what:

    y            the ten adults' mean surface-GluA1 (zref) per structure, ranked
    predictors   four, each one value per structure, entered as x, x^2 and x^3:

      abundance  the four AMPA receptor subunit genes Gria1-4, the mean of their
                 rank profiles. The set is GO:0004971 intersected with
                 GO:0032281, less the delta receptors Grid1 and Grid2 (OVERRIDE in
                 ish.panel_build).
      markers    eleven canonical synaptic markers chosen by hand from the panel:
                 Syp, Syn1, Vamp2, Bsn and Syt1 presynaptic, Dlg4, Homer1,
                 Shank2, Shank3, Nlgn1 and Camk2a postsynaptic. A hand-made list
                 is arguable, hence the next predictor.
      psd_pc1    the first principal component of the 188 postsynaptic-density
                 genes (GO:0014069, less the subunit and localisation sets). No
                 gene is chosen individually: the component is whatever those
                 genes have most in common, and it is the stronger density
                 predictor of the two (0.26 against 0.15).
      autofluo   not a gene: the autofluorescence of the same ten brains per
                 structure, the only predictor measured in the tissue the map
                 comes from.

Writes, in adult_v2/beyond/for_sami/ under the data root:

    E_regression.png/.eps    the fit and its diagnostic
    F_maps.png/.eps          observed, predicted and residual on the brain
    regression_table.csv     every structure: observed, predicted, residual

Run by run_beyond_regression.py.
"""

import csv

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from sepmap.adult.beyond_density import (
    ADULTS,
    MARKERS,
    SUBUNITS,
    build_covariates,
    cv_r2,
    flexible,
    half_map,
    prepare,
    r_squared,
    residual,
)
from sepmap.adult.beyond_figures import FIGS, save
from sepmap.config import SETTINGS
from sepmap.plotting import DARK_GREY, MID_GREY, RED, tidy
from sepmap.volumes.per_mouse import annotation_20, structure_terms
from sepmap.volumes.to_ccf import CCF_CROP

# the coronal planes drawn, the floor of the colour scale (settings.toml says why), and
# the p at which the residuals' tilt or fan calls the model mis-specified
BEYOND_REGRESSION = SETTINGS["beyond_regression"]

# the planes as a tuple, in 20 um planes of the cropped CCF grid
PLANES = tuple(BEYOND_REGRESSION["planes"])


def paint(
    plane: int, value_by_name: dict[str, float], names: dict[int, str]
) -> np.ndarray:
    """One coronal slice with each structure filled by its value, NaN elsewhere."""
    labels = annotation_20("ccf")[plane]
    out = np.full(labels.shape, np.nan, np.float32)
    for idx in np.unique(labels):
        if idx == 0:
            continue
        v = value_by_name.get(names.get(int(idx), ""))
        if v is not None:
            out[labels == idx] = v
    return out


def draw_fit(
    ax: plt.Axes,
    observed: np.ndarray,
    predicted: np.ndarray,
    res: np.ndarray,
    structures: list[str],
    fitted_r2: float,
    cv: float,
) -> None:
    """Draw observed against predicted, the five largest residuals named."""
    lim = [
        min(observed.min(), predicted.min()) - 4,
        max(observed.max(), predicted.max()) + 4,
    ]
    ax.plot(lim, lim, color=MID_GREY, ls="--", lw=1.0, zorder=1)
    ax.scatter(
        predicted,
        observed,
        s=16,
        facecolor=DARK_GREY,
        edgecolor="0.2",
        linewidth=0.3,
        zorder=2,
    )
    worst = np.argsort(-np.abs(res))[:5]
    for i in worst:
        ax.annotate(
            structures[i][:24],
            (predicted[i], observed[i]),
            fontsize=6.5,
            color=RED,
            xytext=(4, 2),
            textcoords="offset points",
        )
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("predicted from abundance and density (rank)", fontsize=8.5)
    ax.set_ylabel("observed surface GluA1 (rank)", fontsize=8.5)
    ax.set_title(f"the fit\nR2 {fitted_r2:.2f} fitted, {cv:.2f} predicted", fontsize=9.5)
    tidy(ax)


def diagnostic(
    predicted: np.ndarray, res: np.ndarray
) -> tuple[float, float, float, float]:
    """The tilt and the fan of the residuals against the prediction.

    Returns Spearman's rho and p of the residual (tilt) and of its size (fan)
    against the prediction; a model that is incomplete rather than mis-specified
    has neither.
    """
    tilt = spearmanr(predicted, res)
    fan = spearmanr(predicted, np.abs(res))
    return (
        float(tilt.statistic),
        float(tilt.pvalue),
        float(fan.statistic),
        float(fan.pvalue),
    )


def draw_diagnostic(
    ax: plt.Axes, predicted: np.ndarray, res: np.ndarray, misspecified: bool
) -> None:
    """Draw the residual against the prediction, the standard diagnostic.

    The title follows `misspecified`, whether a tilt or a fan reached
    beyond_regression.diagnostic_p.
    """
    ax.axhline(0, color=MID_GREY, lw=1.0)
    ax.scatter(
        predicted,
        res,
        s=16,
        facecolor=DARK_GREY,
        edgecolor="0.2",
        linewidth=0.3,
        zorder=2,
    )
    ax.set_xlabel("predicted (rank)", fontsize=8.5)
    ax.set_ylabel("residual (ranks)", fontsize=8.5)
    if misspecified:
        p_max = BEYOND_REGRESSION["diagnostic_p"]
        title = (
            f"the diagnostic\na tilt or a fan (p < {p_max:g}), so the model\n"
            "may be mis-specified"
        )
    else:
        title = (
            "the diagnostic\nno tilt and no fan, so the model is\n"
            "incomplete rather than mis-specified"
        )
    ax.set_title(title, fontsize=9.5)
    tidy(ax)


def draw_leftover(ax: plt.Axes, res: np.ndarray) -> None:
    """Draw the distribution of the residuals, the leftover."""
    ax.hist(res, bins=26, color=MID_GREY, edgecolor="0.3", linewidth=0.4)
    ax.axvline(0, color=DARK_GREY, lw=1.4)
    ax.set_xlabel("residual (ranks)", fontsize=8.5)
    ax.set_ylabel("structures", fontsize=8.5)
    ax.set_title(
        f"the leftover\nspread {res.std():.1f} ranks over {len(res)} structures",
        fontsize=9.5,
    )
    tidy(ax)


def panel_e(
    observed: np.ndarray,
    predicted: np.ndarray,
    res: np.ndarray,
    structures: list[str],
    fitted_r2: float,
    cv: float,
    misspecified: bool,
) -> None:
    """Draw panel E: observed against predicted, the diagnostic, the residuals."""
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.3))

    # observed against predicted, the five largest residuals named; the residual
    # against the prediction; the residuals
    draw_fit(axes[0], observed, predicted, res, structures, fitted_r2, cv)
    draw_diagnostic(axes[1], predicted, res, misspecified)
    draw_leftover(axes[2], res)

    fig.suptitle(
        "E.  The regression behind the claim: surface GluA1 predicted from "
        "receptor abundance and synaptic density",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "E_regression")


def draw_map(
    fig: plt.Figure,
    axes: np.ndarray,
    r: int,
    c: int,
    plane: int,
    img: np.ndarray,
    title: str,
    cmap: str,
    span: float | None,
    structures: list[str],
) -> None:
    """Draw one map, plane `r` of column `c`, with the column's colour bar at the bottom.

    `span` is the symmetric limit of the residual, None for the ranks of the others.
    """
    ax = axes[r, c]
    if span is None:
        lo = 1 - BEYOND_REGRESSION["floor"] * (len(structures) - 1)
        im = ax.imshow(
            img, cmap=cmap, vmin=lo, vmax=len(structures), interpolation="nearest"
        )
    else:
        im = ax.imshow(img, cmap=cmap, vmin=-span, vmax=span, interpolation="nearest")

    # the image rasterised in the EPS, the text kept vector; NaN is
    # transparent, so the black face is the ground
    im.set_rasterized(True)
    ax.set_facecolor("black")
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)
    if r == 0:
        ax.set_title(title, fontsize=9)
    if c == 0:
        # the plane in the full CCF at 10 um, as the route's other figures give it
        ax.set_ylabel(f"CCF plane {2 * (plane + CCF_CROP[0])} / 10 um", fontsize=8)
    if r == len(PLANES) - 1:
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, orientation="horizontal")
        cb.ax.tick_params(labelsize=6.5)
        if span is None:
            label = "rank among structures"
        else:
            label = "observed minus predicted (ranks)"
        cb.set_label(label, fontsize=7)
        if span is None:
            # the floor sits below rank 1 so that nothing draws as black;
            # the ticks still stop at the real range
            cb.set_ticks([1, 25, 50, 75, 100, len(structures)])


def panel_f(
    observed: np.ndarray, predicted: np.ndarray, res: np.ndarray, structures: list[str]
) -> None:
    """Draw panel F: observed, predicted and residual maps on the PLANES."""
    names, _, _ = structure_terms()

    # (title, value per structure, colormap, symmetric limit or None for ranks)
    maps = [
        ("observed\nsurface GluA1", dict(zip(structures, observed)), "hot", None),
        (
            "predicted from abundance\nand synaptic density",
            dict(zip(structures, predicted)),
            "hot",
            None,
        ),
        (
            "what is left over",
            dict(zip(structures, res)),
            "RdBu_r",
            float(np.abs(res).max()),
        ),
    ]

    fig, axes = plt.subplots(len(PLANES), 3, figsize=(10.2, 3.2 * len(PLANES)))
    axes = np.atleast_2d(axes)
    for r, plane in enumerate(PLANES):
        for c, (title, values, cmap, span) in enumerate(maps):
            img = paint(plane, values, names)
            draw_map(fig, axes, r, c, plane, img, title, cmap, span, structures)

    fig.suptitle(
        "F.  The same three quantities on the brain.  Red in the third column is "
        "more surface GluA1 than\nabundance and density predict, blue is less.  "
        "Black is outside the brain, or a structure the\nanalysis excludes: "
        "fibre tracts, ventricles and unassigned voxels.",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "F_maps")


def main() -> None:
    """Fit the quoted model, write the table, and draw panels E and F."""
    FIGS.mkdir(parents=True, exist_ok=True)

    # the structures and the quoted model, as adult.beyond_density builds them
    nano, _, expr, role, auto, structures = prepare()

    # the fit
    covariates, controls, _ = build_covariates(nano, expr, role, auto, structures)
    model = flexible(list(covariates.values()))
    observed = half_map(nano, range(len(ADULTS)), structures)
    res = residual(observed, model)
    predicted = observed - res
    fitted, cv = r_squared(observed, model), cv_r2(observed, model)

    print(
        f"{len(structures)} structures; predictors: Gria1-4 ({len(SUBUNITS)} genes), "
        f"{len(MARKERS)} hand-picked markers, psd_pc1 of {len(controls)} genes, "
        f"autofluorescence"
    )
    print(
        f"  R2 {fitted:.3f} fitted, {cv:.3f} cross-validated; "
        f"observed against predicted rho {spearmanr(observed, predicted).statistic:.3f}"
    )
    tilt, p_tilt, fan, p_fan = diagnostic(predicted, res)
    misspecified = min(p_tilt, p_fan) < BEYOND_REGRESSION["diagnostic_p"]
    print(
        f"  residual spread {res.std():.1f} ranks; "
        f"residual against predicted rho {tilt:+.3f} (p {p_tilt:.2g}), "
        f"|residual| against predicted rho {fan:+.3f} (p {p_fan:.2g}) "
        f"(both should be ~0 if the model is not mis-specified)"
    )

    # every structure, largest residual first
    with open(FIGS / "regression_table.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["structure", "observed_rank", "predicted_rank", "residual"])
        for i in np.argsort(-res):
            writer.writerow(
                [
                    structures[i],
                    f"{observed[i]:.1f}",
                    f"{predicted[i]:.1f}",
                    f"{res[i]:.2f}",
                ]
            )

    panel_e(observed, predicted, res, structures, fitted, cv, misspecified)
    panel_f(observed, predicted, res, structures)
