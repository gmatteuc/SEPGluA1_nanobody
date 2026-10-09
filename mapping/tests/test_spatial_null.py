"""The spatial null: surrogates keep a map's smoothness, and the p is calibrated."""

import numpy as np
import pytest
from scipy.stats import rankdata, spearmanr

from sepmap.ish import spatial_null

# surrogates per map in these tests, enough to place a p against 0.05 in a minute
N_SURROGATES = 500


def distances(seed=0, n=200, box_mm=10.0):
    """The distances in mm between `n` points uniform in a box."""
    rng = np.random.default_rng(seed)
    xyz = rng.uniform(0, box_mm, (n, 3))
    return np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(axis=-1))


def smooth_fields(d, n_maps, rng, length_mm=4.0, nugget=0.3):
    """Maps like the adult one: a large-scale smooth pattern plus structure-level noise.

    A Gaussian random field with a Gaussian covariance of length 4 mm, variance 1,
    plus independent noise of variance 0.3; a different generator from the
    surrogates', so the null is not checked against itself.
    """
    cov = np.exp(-np.square(d / length_mm) / 2) + 1e-6 * np.eye(d.shape[0])
    smooth = rng.standard_normal((n_maps, d.shape[0])) @ np.linalg.cholesky(cov).T
    return smooth + np.sqrt(nugget) * rng.standard_normal(smooth.shape)


def smooth_map(seed=0):
    """Distances between 200 points, and one smooth map on them."""
    d = distances(seed)
    return d, smooth_fields(d, 1, np.random.default_rng(seed + 1000))[0]


def matched_variograms(x, d, surr):
    """The map's rank variogram and the surrogates' mean, on the matched range."""
    grid = spatial_null.variogram_grid(
        d, spatial_null.SPATIAL_NULL["pv"], spatial_null.SPATIAL_NULL["nh"]
    )
    target = spatial_null.smoothed_variogram(rankdata(x), grid)[:, 0]
    mean = spatial_null.smoothed_variogram(surr.T, grid).mean(axis=1)
    return target, mean, grid


def test_surrogates_keep_the_variogram_and_a_shuffle_does_not():
    """The surrogates' variogram is close to the map's; a shuffled map's is not."""
    errors = []
    for seed in range(8):
        d, x = smooth_map(seed)
        surr = spatial_null.surrogates(x, d, N_SURROGATES, seed=seed)
        target, mean, grid = matched_variograms(x, d, surr)
        rng = np.random.default_rng(seed)
        shuffled = np.column_stack(
            [rng.permutation(rankdata(x)) for _ in range(N_SURROGATES)]
        )
        flat = spatial_null.smoothed_variogram(shuffled, grid).mean(axis=1)

        # the mean absolute difference over the matched range, as a share of the
        # map's mean level: single grid points at the shortest distances hold few
        # pairs and a variogram near zero, where a ratio point by point is noise
        level = target.mean()
        error = np.abs(mean - target).mean() / level
        error_shuffled = np.abs(flat - target).mean() / level
        errors.append(error)
        assert error < error_shuffled / 3

    # the smoothest of these random maps is matched less well (one in eight, 43%)
    assert np.median(errors) < 0.15


def test_resampled_surrogates_hold_the_map_ranks():
    """With resample, every surrogate is a permutation of the map's ranks."""
    d, x = smooth_map()
    surr = spatial_null.surrogates(x, d, 50, seed=0, resample=True)
    ranks = np.sort(rankdata(x))
    assert all(np.array_equal(np.sort(s), ranks) for s in surr)


def test_the_same_seed_gives_the_same_surrogates():
    """A seed fixes the surrogates; another seed changes them."""
    d, x = smooth_map()
    a = spatial_null.surrogates(x, d, 20, seed=3)
    b = spatial_null.surrogates(x, d, 20, seed=3)
    c = spatial_null.surrogates(x, d, 20, seed=4)
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)


def test_a_map_with_missing_values_is_refused():
    """Surrogates need a value in every structure."""
    d, x = smooth_map()
    x[5] = np.nan
    with pytest.raises(ValueError):
        spatial_null.surrogates(x, d, 10)


def test_a_planted_correlation_is_found():
    """A gene that is the map plus noise (rho about 0.5) gets p < 0.05 in 9 of 10 maps."""
    hits = 0
    for seed in range(10):
        d, x = smooth_map(seed=seed)
        rng = np.random.default_rng(100 + seed)
        gene = x + 1.5 * x.std() * rng.standard_normal(x.size)
        surr = spatial_null.surrogates(x, d, N_SURROGATES, seed=seed)
        rho = spearmanr(x, gene).statistic
        hits += spatial_null.spatial_p(rho, spatial_null.null_rho(surr, gene)) < 0.05
    assert hits >= 9


def test_independent_smooth_maps_are_not_called_related():
    """Of 100 pairs of unrelated smooth maps few pass the spatial p, many the other."""
    d = distances(7)
    rng = np.random.default_rng(7)
    first = smooth_fields(d, 100, rng)
    second = smooth_fields(d, 100, rng)
    n_spatial = n_ordinary = 0
    for k in range(100):
        test = spearmanr(first[k], second[k])
        surr = spatial_null.surrogates(first[k], d, N_SURROGATES, seed=k)
        n_spatial += (
            spatial_null.spatial_p(test.statistic, spatial_null.null_rho(surr, second[k]))
            < 0.05
        )
        n_ordinary += test.pvalue < 0.05
    assert n_spatial <= 10
    assert n_ordinary >= 25


def test_the_p_counts_the_observed_value_once():
    """No surrogate as extreme gives 1 / (n + 1); every one as extreme gives 1."""
    null = np.linspace(-0.3, 0.3, 99)
    assert spatial_null.spatial_p(0.9, null) == pytest.approx(1 / 100)
    assert spatial_null.spatial_p(-0.9, null) == pytest.approx(1 / 100)
    assert spatial_null.spatial_p(0.0, null) == pytest.approx(1.0)


def test_null_rho_is_spearman_per_surrogate():
    """The vectorised null equals scipy's Spearman for each surrogate."""
    rng = np.random.default_rng(0)
    surr = rng.standard_normal((5, 40))
    gene = rng.standard_normal(40)
    want = [spearmanr(s, gene).statistic for s in surr]
    assert np.allclose(spatial_null.null_rho(surr, gene), want)
