"""A spatial null for correlations across structures: surrogates with the same smoothness.

Two brain maps agree partly because both are smooth: neighbouring structures have
similar values, cortex and hippocampus are high, thalamus and hypothalamus low. A
map that knows only each structure's division already correlates strongly with many
genes (analysis 2, ish.divisions). Shuffling the structures destroys that
smoothness, so a permutation null is far too narrow and every rho looks significant
(the inflation Fulcher 2021 measured in mouse). The null here keeps the smoothness,
with variogram-matched surrogates
(Burt et al. 2020, NeuroImage 220:117038, the method of brainsmash's Base class,
written here in numpy because brainsmash needs scikit-learn and joblib, which
venv_atlas does not have):

    the map      its ranks: every test here is a rank correlation, so the
                 smoothness that matters is that of the ranks; and the values are
                 skewed (the hippocampus sits far above the rest), so their
                 variogram is set by a few large values, which the surrogates
                 match less well than they match the ranks
    distances    between the structures' centroids in one hemisphere, in mm
                 (structures.centroids)
    variogram    half the squared difference of every pair of structures, against
                 their distance, for the pairs closer than the pv-th percentile of
                 all distances (the short range, where the autocorrelation is),
                 smoothed with a Gaussian kernel onto nh distances; the kernel's
                 width is three steps of that grid
    surrogate    the map's values shuffled, then smoothed over each structure's
                 nearest neighbours (excluding itself) with an exponential kernel,
                 for each share delta of the structures; each smoothed map's
                 variogram is regressed on the map's, and the share that fits best
                 is kept, scaled by the slope, with noise of the intercept's size
                 added: sqrt|beta| x smoothed + sqrt|alpha| x N(0, 1)
    resample     the surrogate's values are then replaced by the map's ranks, in
                 the surrogate's order, so it differs from the map only in where
                 the values sit

The smoothing for one share is one fixed weight matrix, so thousands of surrogates
are a few matrix products. A gene is tested by correlating it with every surrogate
on the structures it shares with the map; its two-sided p counts the observed value
once,

    p = (1 + number of surrogates with |rho| >= |observed rho|) / (1 + surrogates)

The null is checked on maps whose answer is known (calibration): pairs of
independent Gaussian random fields on the same centroids, with an exponential
covariance fitted to the rising part of the map's own variogram (another generator
than the surrogates', so the null is not checked against itself), and fields tested
against the map itself with its surrogates. Independent maps should give p < 0.05
in 5% of tests. The surrogates match smooth maps best; a map rough at short range is
matched less well, and its p then errs on the conservative side.

Run by run_ish_spatial_null.py.
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import rankdata, spearmanr

from sepmap.config import SETTINGS
from sepmap.ish.numbers import numbers_frame, numbers_path
from sepmap.structures import TABLES

# the number of surrogates, the shares of neighbours, the variogram's range and
# grid, resampling, and the size of the calibration
SPATIAL_NULL = SETTINGS["spatial_null"]

SURROGATES = {
    "nano": TABLES / "surrogates_nano.npy",
    "auto": TABLES / "surrogates_auto.npy",
}
SURROGATE_STRUCTURES = TABLES / "surrogate_structures.csv"
VARIOGRAM = TABLES / "variogram.csv"
CALIBRATION = TABLES / "null_calibration.csv"
NUMBERS = numbers_path("spatial_null")

# the variogram's Gaussian kernel has an SD of its width over 2.68, so its weight
# falls to about 3% one width from its centre, and the width is three steps of the
# distance grid (Burt 2020)
KERNEL_FACTOR = 2.68
KERNEL_STEPS = 3.0

# surrogates drawn per matrix product; larger batches use more memory, not less time
BATCH = 500

# distances at which a variogram is drawn or fitted over its whole range
FIT_STEPS = 50

# the level at which the calibration counts false positives
ALPHA = 0.05

# the percentiles of a null distribution that a value must leave to pass at ALPHA
BAND = (2.5, 97.5)


# ===== Variograms =====


def distance_matrix(centroids: pd.DataFrame) -> np.ndarray:
    """Euclidean distances between the rows of (ap_mm, dv_mm, ml_mm), in mm."""
    xyz = centroids[["ap_mm", "dv_mm", "ml_mm"]].to_numpy(float)
    return np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(axis=-1))


def variogram_grid(d: np.ndarray, pv: float, nh: int) -> dict:
    """The pairs, distance grid and kernel weights of a smoothed variogram.

    Pairs closer than the pv-th percentile of all pairwise distances enter; the
    grid runs from the shortest to the longest of them in nh steps. Returns i, j
    (the pairs), u (their distances), h (the grid) and weights (nh x pairs, each
    row summing to 1).
    """
    i, j = np.triu_indices(d.shape[0], k=1)
    u = d[i, j]
    keep = u < np.percentile(u, pv)
    i, j, u = i[keep], j[keep], u[keep]
    h = np.linspace(u.min(), u.max(), nh)
    width = KERNEL_STEPS * (h[1] - h[0])
    w = np.exp(-np.square(KERNEL_FACTOR * np.abs(u[None, :] - h[:, None]) / width) / 2)
    return dict(i=i, j=j, u=u, h=h, weights=w / w.sum(axis=1, keepdims=True))


def smoothed_variogram(maps: np.ndarray, grid: dict) -> np.ndarray:
    """The smoothed variogram of each column of `maps` (structures x maps), nh x maps."""
    if maps.ndim == 1:
        maps = maps[:, None]
    half_sq = 0.5 * np.square(maps[grid["j"]] - maps[grid["i"]])
    return grid["weights"] @ half_sq


# ===== Surrogates =====


def neighbour_weights(d: np.ndarray, delta: float) -> np.ndarray:
    """Row-normalised exponential weights over each structure's nearest neighbours.

    The int(delta x n) nearest structures, the structure itself left out, weighted
    exp(-d / the farthest of them), as Burt 2020's exponential kernel; zero
    elsewhere.
    """
    n = d.shape[0]
    k = int(delta * n)
    order = np.argsort(d, axis=1)[:, 1 : k + 1]
    rows = np.arange(n)[:, None]
    near = d[rows, order]
    w = np.exp(-near / near.max(axis=1, keepdims=True))
    out = np.zeros((n, n))
    out[rows, order] = w / w.sum(axis=1, keepdims=True)
    return out


def fit_scale(target: np.ndarray, variograms: np.ndarray) -> tuple:
    """Least squares of the target variogram on each column: (alpha, beta, sse)."""
    v_mean = variograms.mean(axis=0)
    t_mean = target.mean()
    dv = variograms - v_mean
    beta = (dv * (target - t_mean)[:, None]).sum(axis=0) / np.square(dv).sum(axis=0)
    alpha = t_mean - beta * v_mean
    sse = np.square(target[:, None] - alpha - beta * variograms).sum(axis=0)
    return alpha, beta, sse


def surrogate_batch(
    x: np.ndarray,
    weights: list[np.ndarray],
    grid: dict,
    target: np.ndarray,
    n: int,
    rng: np.random.Generator,
    resample: bool,
) -> np.ndarray:
    """`n` surrogates of `x`, one per column, from the weight matrix of each share."""
    shuffled = np.column_stack([rng.permutation(x) for _ in range(n)])
    best_sse = np.full(n, np.inf)
    best = np.zeros_like(shuffled)
    best_alpha = np.zeros(n)
    best_beta = np.zeros(n)
    for w in weights:
        smoothed = w @ shuffled
        alpha, beta, sse = fit_scale(target, smoothed_variogram(smoothed, grid))
        better = sse < best_sse
        best_sse[better] = sse[better]
        best[:, better] = smoothed[:, better]
        best_alpha[better] = alpha[better]
        best_beta[better] = beta[better]
    noise = rng.standard_normal(best.shape)
    surr = np.sqrt(np.abs(best_beta)) * best + np.sqrt(np.abs(best_alpha)) * noise
    if resample:
        order = np.argsort(surr, axis=0)
        values = np.repeat(np.sort(x)[:, None], n, axis=1)
        np.put_along_axis(surr, order, values, axis=0)
    else:
        surr = surr - surr.mean(axis=0)
    return surr


def surrogates(
    x: np.ndarray,
    d: np.ndarray,
    n_surrogates: int | None = None,
    seed: int = 0,
    resample: bool | None = None,
) -> np.ndarray:
    """Variogram-matched surrogates of the ranks of `x` on distances `d`.

    Returns surrogates x structures, each a map of ranks (1 to n). The settings of
    [spatial_null] give the shares of neighbours, the variogram's range (pv) and
    grid (nh), and by default the number of surrogates and resampling.
    """
    if n_surrogates is None:
        n_surrogates = SPATIAL_NULL["n_surrogates"]
    if resample is None:
        resample = SPATIAL_NULL["resample"]
    x = np.asarray(x, dtype=float)
    if not np.isfinite(x).all():
        raise ValueError("the map has missing values; surrogates need a value everywhere")
    x = rankdata(x)
    rng = np.random.default_rng(seed)
    grid = variogram_grid(d, SPATIAL_NULL["pv"], SPATIAL_NULL["nh"])
    target = smoothed_variogram(x, grid)[:, 0]
    weights = [neighbour_weights(d, delta) for delta in SPATIAL_NULL["deltas"]]
    out = []
    for start in range(0, n_surrogates, BATCH):
        n = min(BATCH, n_surrogates - start)
        out.append(surrogate_batch(x, weights, grid, target, n, rng, resample).T)
    return np.vstack(out)


# ===== p values =====


def null_rho(surr: np.ndarray, gene: np.ndarray) -> np.ndarray:
    """Spearman of every surrogate (a row) with `gene`, over the same structures."""
    r_surr = rankdata(surr, axis=1)
    r_gene = rankdata(gene)
    r_surr = r_surr - r_surr.mean(axis=1, keepdims=True)
    r_gene = r_gene - r_gene.mean()
    num = r_surr @ r_gene
    den = np.sqrt(np.square(r_surr).sum(axis=1) * np.square(r_gene).sum())
    return num / den


def spatial_p(observed: float, null: np.ndarray) -> float:
    """Two-sided p of `observed` against a null sample, counting the observed once."""
    k = int(np.sum(np.abs(null) >= abs(observed)))
    return (k + 1) / (len(null) + 1)


# ===== Checks of the null =====


def exponential_variogram(h, nugget, sill, range_mm):
    """An exponential variogram: nugget + sill x (1 - exp(-h / range))."""
    return nugget + sill * (1 - np.exp(-h / range_mm))


def fit_exponential(x: np.ndarray, d: np.ndarray) -> tuple[float, float, float]:
    """(nugget, sill, range in mm) of an exponential variogram fitted to the map's.

    The map's ranks, scaled to 0-1, as the surrogates take them. The variogram of
    a map on a bounded brain rises, then falls again at the longest distances
    (structures at opposite ends are alike more often than chance), which no
    exponential describes; so the fit runs up to the distance where the map's
    smoothed variogram peaks, over every pair closer than that.
    """
    grid = variogram_grid(d, 100.0, FIT_STEPS)
    v = smoothed_variogram(rankdata(x) / len(x), grid)[:, 0]
    rising = grid["h"] <= grid["h"][np.argmax(v)]
    start = (float(v.min()), float(v.max() - v.min()), 1.0)
    bounds = ([0, 0, 0.05], [np.inf, np.inf, 50.0])
    params, _ = curve_fit(
        exponential_variogram, grid["h"][rising], v[rising], p0=start, bounds=bounds
    )
    return tuple(float(p) for p in params)


def random_fields(
    d: np.ndarray, params: tuple, n_maps: int, rng: np.random.Generator
) -> np.ndarray:
    """Gaussian random fields with an exponential covariance, maps x structures."""
    nugget, sill, range_mm = params
    cov = sill * np.exp(-d / range_mm) + 1e-9 * np.eye(d.shape[0])
    chol = np.linalg.cholesky(cov)
    smooth = rng.standard_normal((n_maps, d.shape[0])) @ chol.T
    return smooth + np.sqrt(nugget) * rng.standard_normal(smooth.shape)


def calibration_table(
    d: np.ndarray, params: tuple, map_values: np.ndarray, surr: np.ndarray, seed: int = 0
) -> pd.DataFrame:
    """Independent smooth maps tested as genes are, one row per test.

    Two designs: pairs of independent random fields, each pair tested with
    surrogates of its first map (`pair`); and random fields tested against the map
    itself with its own surrogates `surr` (`map`), the use the null is built for.
    Per test: rho, the ordinary Spearman p, the spatial p, and rho with the second
    map shuffled.
    """
    rng = np.random.default_rng(seed)
    n_maps = SPATIAL_NULL["n_calibration"]
    n_surr = SPATIAL_NULL["n_calibration_surrogates"]
    rows = []
    first = random_fields(d, params, n_maps, rng)
    second = random_fields(d, params, n_maps, rng)
    for k in range(n_maps):
        a, b = first[k], second[k]
        own = surrogates(a, d, n_surr, seed=seed + 1 + k)
        test = spearmanr(a, b)
        null = null_rho(own, b)
        rows.append(
            dict(
                design="pair",
                test=k,
                rho=float(test.statistic),
                p_ordinary=float(test.pvalue),
                p_spatial=spatial_p(test.statistic, null),
                rho_shuffled=float(spearmanr(a, rng.permutation(b)).statistic),
            )
        )
    for k in range(n_maps):
        b = second[k]
        test = spearmanr(map_values, b)
        rows.append(
            dict(
                design="map",
                test=k,
                rho=float(test.statistic),
                p_ordinary=float(test.pvalue),
                p_spatial=spatial_p(test.statistic, null_rho(surr, b)),
                rho_shuffled=float(spearmanr(map_values, rng.permutation(b)).statistic),
            )
        )
    return pd.DataFrame(rows)


def variogram_table(
    name: str, x: np.ndarray, d: np.ndarray, surr: np.ndarray, seed: int = 0
) -> pd.DataFrame:
    """The map's variogram over every distance, its surrogates' and a shuffled map's.

    One row per distance: the map's smoothed variogram, the median and the 5th and
    95th percentiles over the surrogates, the median over shuffled maps (as many
    as surrogates), and whether the distance is inside the matched range (pv).
    """
    rng = np.random.default_rng(seed)
    grid = variogram_grid(d, 100.0, FIT_STEPS)
    fitted_max = variogram_grid(d, SPATIAL_NULL["pv"], SPATIAL_NULL["nh"])["h"][-1]
    ranks = rankdata(x)
    shuffled = np.column_stack([rng.permutation(ranks) for _ in range(surr.shape[0])])
    v_map = smoothed_variogram(ranks, grid)[:, 0]
    v_surr = smoothed_variogram(surr.T, grid)
    v_shuffled = smoothed_variogram(shuffled, grid)
    return pd.DataFrame(
        dict(
            map=name,
            distance_mm=grid["h"],
            variogram=v_map,
            surrogate_median=np.median(v_surr, axis=1),
            surrogate_p5=np.percentile(v_surr, 5, axis=1),
            surrogate_p95=np.percentile(v_surr, 95, axis=1),
            shuffled_median=np.median(v_shuffled, axis=1),
            matched=grid["h"] <= fitted_max,
        )
    )


def false_positive_rates(calibration: pd.DataFrame) -> pd.DataFrame:
    """Per design, the share of tests below ALPHA by the ordinary and the spatial p."""
    rows = []
    for design, sub in calibration.groupby("design", sort=False):
        rows.append(
            dict(
                design=design,
                n_tests=len(sub),
                ordinary=float((sub["p_ordinary"] < ALPHA).mean()),
                spatial=float((sub["p_spatial"] < ALPHA).mean()),
                sd_rho=float(sub["rho"].std()),
                sd_rho_shuffled=float(sub["rho_shuffled"].std()),
            )
        )
    return pd.DataFrame(rows)


def matched_misfit(variogram: pd.DataFrame) -> dict[str, tuple[float, float]]:
    """Per map, the median and largest |surrogates - map| / map in the matched range.

    `variogram` is variogram_table's, one map or several.
    """
    out = {}
    for name, sub in variogram[variogram["matched"]].groupby("map"):
        rel = (sub["surrogate_median"] - sub["variogram"]).abs() / sub["variogram"]
        out[name] = (float(rel.median()), float(rel.max()))
    return out


def numbers_table(
    n_structures: int,
    params: tuple,
    misfit: dict[str, tuple[float, float]],
    rates: pd.DataFrame,
) -> pd.DataFrame:
    """numbers_spatial_null.csv: the numbers of this step that the text quotes."""
    rows = [
        ("null_structures", n_structures, "declared structures the surrogates cover"),
        ("null_surrogates", SPATIAL_NULL["n_surrogates"], "surrogates per map"),
        ("calibration_range_mm", round(params[2], 3), "range of the fitted fields, mm"),
        ("calibration_nugget", round(params[0], 4), "nugget of the fitted fields"),
    ]
    for name, (median, largest) in misfit.items():
        rows.append(
            (f"variogram_misfit_{name}", round(median, 3), "median, matched range")
        )
        rows.append(
            (f"variogram_misfit_max_{name}", round(largest, 3), "largest, matched range")
        )
    for _, r in rates.iterrows():
        what = f"{r['design']} design, {r['n_tests']} tests"
        rows.append((f"fpr_ordinary_{r['design']}", round(r["ordinary"], 4), what))
        rows.append((f"fpr_spatial_{r['design']}", round(r["spatial"], 4), what))
        rows.append((f"sd_rho_{r['design']}", round(r["sd_rho"], 3), "smooth maps"))
        rows.append(
            (f"sd_rho_shuffled_{r['design']}", round(r["sd_rho_shuffled"], 3), "shuffled")
        )
    return numbers_frame(rows)


def save_surrogates(surr: dict[str, np.ndarray], centroids: pd.DataFrame) -> None:
    """Write each map's surrogates (float32) and the structures of their columns.

    `centroids` is the declared structures' table in the order of the columns.
    """
    table = centroids.reset_index()[["structure", "ap_mm", "dv_mm", "ml_mm"]]
    table.to_csv(SURROGATE_STRUCTURES, index=False)
    for name, maps in surr.items():
        np.save(SURROGATES[name], maps.astype(np.float32))


def load_surrogates(
    map_name: str, declared: list[str] | None = None
) -> tuple[np.ndarray, list[str]]:
    """The surrogates of one map (nano or auto) and the structures of their columns.

    With `declared`, the run stops unless the surrogates were drawn on exactly those
    structures, in that order.
    """
    path = SURROGATES[map_name]
    if not path.exists():
        raise FileNotFoundError(f"{path} not found: run run_ish_spatial_null.py first")
    structures = pd.read_csv(SURROGATE_STRUCTURES)["structure"].tolist()
    surr = np.load(path)
    if surr.shape[1] != len(structures):
        raise ValueError(
            f"{path} has {surr.shape[1]} columns, {SURROGATE_STRUCTURES.name} lists "
            f"{len(structures)} structures; run run_ish_spatial_null.py --recompute"
        )
    if declared is not None and structures != declared:
        raise ValueError(
            "the surrogates were drawn on other structures than the declared set; "
            "run run_ish_spatial_null.py --recompute"
        )
    return surr, structures
