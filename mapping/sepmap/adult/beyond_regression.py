"""The regression of analysis 4 shown as a regression: map, prediction, leftover.

adult.beyond_density reports the model's summary numbers (R2, replication) and
adult.beyond_controls its controls; this module writes the fit itself, for three
views of it (drawn by adult.plotting):

    1  the map against the prediction, one dot per structure, with the identity
       line. A good model puts the cloud on that line; the spread away from it is
       the leftover, drawn rather than summarised.
    2  the leftover against the prediction, the standard diagnostic. A tilt or a
       fan here would mean the model is mis-specified rather than incomplete.
    3  the same three quantities painted back onto the brain (map, prediction and
       leftover), because the leftover is a spatial claim.

What is regressed on what is the main model of adult.beyond_density, which has the
predictors and the reasons: the ten adults' mean zref per structure, ranked, on the
ranks of Gria1 and of synapse density (the mean rank of the density genes), two
straight terms.

Guided figure 04 (adult.beyond_figures) draws its maps from regression_table.csv.

Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    regression_table.csv    every structure: map, prediction, leftover (ranks)

Run by run_beyond_regression.py, which draws the working figures E_regression.png
(the fit and its diagnostic) and F_maps.png (the three maps on the brain).
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from sepmap.adult.beyond_density import (
    OUT,
    covariates_for,
    cv_r2,
    full_map,
    load_inputs,
    model_columns,
    r_squared,
    residual,
)
from sepmap.config import SETTINGS
from sepmap.volumes.per_mouse import annotation_20, structure_terms
from sepmap.volumes.to_ccf import CCF_CROP

# the coronal planes drawn, the floor of the colour scale (settings.toml says why), and
# the p at which the residuals' tilt or fan calls the model mis-specified
BEYOND_REGRESSION = SETTINGS["beyond_regression"]

# the planes as a tuple, in 20 um planes of the cropped CCF grid
PLANES = tuple(BEYOND_REGRESSION["planes"])

REGRESSION = OUT / "regression_table.csv"


def coronal_planes() -> list[dict]:
    """The PLANES as the maps draw them, one dict each.

    lab: the parcellation indices of the plane, (DV, ML); ccf_plane: its number in
    the full CCF at 10 um, as the route's other figures give it; names: the
    structure of each index.
    """
    names, _, _ = structure_terms()
    atlas = annotation_20("ccf")
    return [
        dict(lab=np.asarray(atlas[p]), ccf_plane=2 * (p + CCF_CROP[0]), names=names)
        for p in PLANES
    ]


def diagnostic(
    predicted: np.ndarray, leftover: np.ndarray
) -> tuple[float, float, float, float]:
    """The tilt and the fan of the residuals against the prediction.

    Returns Spearman's rho and p of the residual (tilt) and of its size (fan)
    against the prediction; a model that is incomplete rather than mis-specified
    has neither.
    """
    tilt = spearmanr(predicted, leftover)
    fan = spearmanr(predicted, np.abs(leftover))
    return (
        float(tilt.statistic),
        float(tilt.pvalue),
        float(fan.statistic),
        float(fan.pvalue),
    )


def regression_table(
    structures: list[str],
    acronym: dict[str, str],
    division: dict[str, str],
    observed: np.ndarray,
    predicted: np.ndarray,
    leftover: np.ndarray,
) -> pd.DataFrame:
    """regression_table.csv: every structure, the largest leftover first."""
    table = pd.DataFrame(
        dict(
            structure=structures,
            acronym=[acronym.get(s, "") for s in structures],
            division=[division.get(s, "") for s in structures],
            observed_rank=observed,
            predicted_rank=predicted,
            residual=leftover,
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


def main() -> dict:
    """Fit the main model, write the table, and print the diagnostic.

    Returns what the working figures draw: the map, the prediction and the leftover
    (ranks), the structures, the R2 fitted and held out, and whether the diagnostic
    calls the model mis-specified.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs()
    structures = inputs.structures
    covariates, _, _ = covariates_for(inputs)
    columns = model_columns(covariates, inputs.terms)
    observed = full_map(inputs.nano)
    leftover = residual(observed, columns)
    predicted = observed - leftover
    fitted, cv = r_squared(observed, columns), cv_r2(observed, columns)
    terms = [", ".join(inputs.terms[group]) for group in inputs.terms]
    print(f"{len(structures)} structures; predictors: {'; '.join(terms)}, straight")
    print(
        f"  R2 {fitted:.3f} fitted, {cv:.3f} cross-validated; "
        f"map against prediction rho {spearmanr(observed, predicted).statistic:.3f}"
    )
    tilt, p_tilt, fan, p_fan = diagnostic(predicted, leftover)
    misspecified = min(p_tilt, p_fan) < BEYOND_REGRESSION["diagnostic_p"]
    print(
        f"  residual spread {leftover.std():.1f} ranks; "
        f"residual against predicted rho {tilt:+.3f} (p {p_tilt:.2g}), "
        f"|residual| against predicted rho {fan:+.3f} (p {p_fan:.2g}) "
        "(both near 0 when the model is not mis-specified)"
    )
    regression_table(
        structures, inputs.acronym, inputs.division, observed, predicted, leftover
    ).to_csv(REGRESSION, index=False)
    print(f"  -> {REGRESSION}")
    return dict(
        observed=observed,
        predicted=predicted,
        leftover=leftover,
        structures=structures,
        fitted=fitted,
        cv=cv,
        misspecified=misspecified,
    )
