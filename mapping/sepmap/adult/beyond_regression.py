"""The regression of analysis 4 shown as a regression: map, prediction, leftover.

adult.beyond_density reports the model's summary numbers (R2, replication) and
adult.beyond_controls its controls; this module shows the fit itself, three ways:

    1  the map against the prediction, one dot per structure, with the identity
       line. A good model puts the cloud on that line; the spread away from it is
       the leftover, drawn rather than summarised.
    2  the leftover against the prediction, the standard diagnostic. A tilt or a
       fan here would mean the model is mis-specified rather than incomplete.
    3  the same three quantities painted back onto the brain (map, prediction and
       leftover), because the leftover is a spatial claim.

What is regressed on what (adult.beyond_density has the reasons):

    y            the ten adults' mean zref per structure, ranked
    predictors   seven, one value per structure each, entered as x, x^2 and x^3:

      Gria1-4    the four AMPA receptor subunit genes, each its own term. The set is
                 GO:0004971 intersected with GO:0032281, less the delta receptors
                 Grid1 and Grid2 (ish.panel_build)
      markers    the synaptic marker genes of [beyond] markers, chosen by hand from
                 the panel, presynaptic and postsynaptic, as the mean of their ranks.
                 A hand-made list is arguable, hence the next predictor
      psd_pc1    the first principal component of the postsynaptic-density genes
                 of the ontology panel (role control_psd) measured in every
                 structure. No gene is chosen individually: the component is
                 whatever those genes have most in common
      autofluo   not a gene: the autofluorescence of the same ten brains per
                 structure, the only predictor measured in the tissue the map comes
                 from

Guided figure 12 (adult.beyond_figures) draws its maps from regression_table.csv.

Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    regression_table.csv    every structure: map, prediction, leftover (ranks)
    E_regression.png        the fit and its diagnostic (a working figure)
    F_maps.png              map, prediction and leftover on the brain (working)

Run by run_beyond_regression.py.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap.adult.beyond_density import (
    MARKERS,
    OUT,
    SUBUNITS,
    build_covariates,
    cv_r2,
    full_map,
    load_inputs,
    model,
    r_squared,
    residual,
    save,
)
from sepmap.config import SETTINGS
from sepmap.plotting import DARK_GREY, MID_GREY, RED, tidy
from sepmap.volumes.per_mouse import annotation_20, structure_terms
from sepmap.volumes.to_ccf import CCF_CROP

# the coronal planes drawn, the floor of the colour scale (settings.toml says why), and
# the p at which the residuals' tilt or fan calls the model mis-specified
BEYOND_REGRESSION = SETTINGS["beyond_regression"]

# the planes as a tuple, in 20 um planes of the cropped CCF grid
PLANES = tuple(BEYOND_REGRESSION["planes"])

REGRESSION = OUT / "regression_table.csv"


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
    ax.set_xlabel("predicted from receptor mRNA and density (rank)", fontsize=8.5)
    ax.set_ylabel("nano map (rank)", fontsize=8.5)
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
        "E.  The regression of analysis 4: the nano map predicted from receptor "
        "mRNA and synaptic density",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "E_regression.png")


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
        ("the nano map", dict(zip(structures, observed)), "hot", None),
        (
            "predicted from receptor mRNA\nand synaptic density",
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
        "a higher nano rank than\nreceptor mRNA and density predict, blue is lower.  "
        "Black is outside the brain, or a structure the\nfit does not use.",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "F_maps.png")


def regression_table(
    structures: list[str],
    acronym: dict[str, str],
    division: dict[str, str],
    observed: np.ndarray,
    predicted: np.ndarray,
    res: np.ndarray,
) -> pd.DataFrame:
    """regression_table.csv: every structure, the largest leftover first."""
    table = pd.DataFrame(
        dict(
            structure=structures,
            acronym=[acronym.get(s, "") for s in structures],
            division=[division.get(s, "") for s in structures],
            observed_rank=observed,
            predicted_rank=predicted,
            residual=res,
        )
    )
    return table.sort_values("residual", ascending=False, ignore_index=True)


def load_regression() -> pd.DataFrame:
    """The table run_beyond_regression wrote."""
    if not REGRESSION.exists():
        raise FileNotFoundError(
            f"{REGRESSION} not found: run run_beyond_regression.py first"
        )
    return pd.read_csv(REGRESSION, keep_default_na=False, na_values=[""])


def main() -> None:
    """Fit the model, write the table, and draw the working panels E and F."""
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs()
    structures = inputs.structures
    covariates, psd, _ = build_covariates(
        inputs.expr, inputs.role, inputs.auto, structures
    )
    xs = model(covariates)
    observed = full_map(inputs.nano)
    res = residual(observed, xs)
    predicted = observed - res
    fitted, cv = r_squared(observed, xs), cv_r2(observed, xs)
    print(
        f"{len(structures)} structures; predictors: {', '.join(SUBUNITS)} as separate "
        f"terms, {len(MARKERS)} marker genes, psd_pc1 of {len(psd)} genes, "
        "autofluorescence, each bent"
    )
    print(
        f"  R2 {fitted:.3f} fitted, {cv:.3f} cross-validated; "
        f"map against prediction rho {spearmanr(observed, predicted).statistic:.3f}"
    )
    tilt, p_tilt, fan, p_fan = diagnostic(predicted, res)
    misspecified = min(p_tilt, p_fan) < BEYOND_REGRESSION["diagnostic_p"]
    print(
        f"  residual spread {res.std():.1f} ranks; "
        f"residual against predicted rho {tilt:+.3f} (p {p_tilt:.2g}), "
        f"|residual| against predicted rho {fan:+.3f} (p {p_fan:.2g}) "
        "(both near 0 when the model is not mis-specified)"
    )
    regression_table(
        structures, inputs.acronym, inputs.division, observed, predicted, res
    ).to_csv(REGRESSION, index=False)
    print(f"  -> {REGRESSION}")
    panel_e(observed, predicted, res, structures, fitted, cv, misspecified)
    panel_f(observed, predicted, res, structures)
