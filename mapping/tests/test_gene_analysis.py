"""Known-answer checks of analyses 1 to 3: gene ranking, divisions, gene sets."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from sepmap import config
from sepmap.ish import divisions, gene_ranking, gene_sets, panel_test

TABLES = Path(config.DATA) / "adult_v2" / "ish_analysis" / "tables"


def smooth_map(n=120, seed=0):
    """A map over n structures along a line, smooth, with its positions."""
    rng = np.random.default_rng(seed)
    x = np.sort(rng.uniform(0, 10, n))
    return x, np.sin(x) + 0.3 * rng.standard_normal(n)


def test_planted_gene_ranks_first():
    """A gene equal to the map plus a little noise is first, a random gene is not."""
    rng = np.random.default_rng(1)
    _, m = smooth_map()
    structures = [f"s{i}" for i in range(len(m))]
    profiles = {
        "planted": dict(zip(structures, m + 0.05 * rng.standard_normal(len(m)))),
        "random": dict(zip(structures, rng.standard_normal(len(m)))),
        "inverse": dict(zip(structures, -m)),
    }
    vectors = gene_ranking.gene_vectors(profiles, structures)
    surr = np.array([rng.permutation(m) for _ in range(200)])
    boot = np.array([m + 0.1 * rng.standard_normal(len(m)) for _ in range(50)])
    table, null = gene_ranking.rank_genes(m, surr, boot, vectors)
    table = gene_ranking.with_q_and_ranks(table, {"planted", "random", "inverse"})
    by = table.set_index("symbol")
    assert by.loc["planted", "rank_all"] == 1
    assert by.loc["inverse", "rank_all"] == 3
    assert by.loc["planted", "rho"] > 0.95
    assert null.shape == (3, 200)


def test_benjamini_hochberg_matches_hand_worked_example():
    """q of p = 0.01, 0.04, 0.03, 0.005 is 0.02, 0.04, 0.04, 0.02 (worked by hand)."""
    table = pd.DataFrame(
        dict(symbol=["a", "b", "c", "d"], rho=[0.4, 0.3, 0.2, 0.5])
    ).assign(p_spatial=[0.01, 0.04, 0.03, 0.005])
    out = gene_ranking.with_q_and_ranks(table, {"a", "b", "c", "d"})
    assert np.allclose(out["q_all"], [0.02, 0.04, 0.04, 0.02])
    assert np.allclose(out["q_p9"], out["q_all"])


def test_q_within_p9_uses_only_p9_genes():
    """A gene outside P9's panel gets no q_p9 and no rank_p9."""
    table = pd.DataFrame(dict(symbol=["a", "b"], rho=[0.4, 0.3], p_spatial=[0.01, 0.2]))
    out = gene_ranking.with_q_and_ranks(table, {"a"})
    assert np.isnan(out.loc[1, "q_p9"])
    assert np.isnan(out.loc[1, "rank_p9"])
    assert out.loc[0, "q_p9"] == pytest.approx(0.01)


def test_bootstrap_maps_are_means_of_resampled_adults():
    """A resampled map is a mean of the adults' rows, so it stays within their range."""
    rng = np.random.default_rng(2)
    adults = rng.standard_normal((10, 30))
    boot = gene_ranking.bootstrap_maps(adults, n=100, seed=0)
    assert boot.shape == (100, 30)
    assert (boot >= adults.min(axis=0) - 1e-12).all()
    assert (boot <= adults.max(axis=0) + 1e-12).all()


def test_gap_uses_only_structures_both_genes_have():
    """The gap is computed on the intersection of the two genes' structures."""
    rng = np.random.default_rng(3)
    _, m = smooth_map(80)
    structures = [f"s{i}" for i in range(80)]
    first = {s: v for s, v in zip(structures[:70], m[:70])}
    second = {s: v for s, v in zip(structures[20:], rng.standard_normal(60))}
    surr = np.array([rng.permutation(m) for _ in range(100)])
    row, null = gene_ranking.gap_row(m, surr, surr[:10], first, second, structures)
    assert row["n_structures"] == 50
    shared = list(range(20, 70))
    expected = (
        spearmanr(m[shared], m[shared]).statistic
        - spearmanr(m[shared], [second[structures[i]] for i in shared]).statistic
    )
    assert row["gap"] == pytest.approx(expected)
    assert null.shape == (100,)


def test_division_only_map_has_one_value_per_division():
    """Every structure of a division takes the division's median."""
    values = np.array([1.0, 2.0, 3.0, 10.0, 20.0])
    labels = np.array(["A", "A", "A", "B", "B"])
    out = divisions.division_only(values, labels)
    assert np.allclose(out, [2.0, 2.0, 2.0, 15.0, 15.0])


def test_shuffle_within_divisions_keeps_each_divisions_values():
    """A shuffled map holds, inside each division, a permutation of its own values."""
    rng = np.random.default_rng(4)
    values = np.arange(30, dtype=float)
    labels = np.array(["A"] * 10 + ["B"] * 12 + ["C"] * 8)
    parts = divisions.division_parts(np.arange(30), labels, min_structures=8)
    maps = divisions.shuffle_within(values, parts, 20, rng)
    for _, positions in parts:
        for row in maps:
            assert sorted(row[positions]) == sorted(values[positions])
    # the order inside a division changes, or the null would hold only the map itself
    assert (maps != values).any(axis=1).all()
    out = divisions.shuffled_within(values, values, parts, 50, rng)
    assert out.shape == (50,)


def test_within_rho_of_one_division_is_its_spearman():
    """With one division, the within rho is the plain Spearman inside it."""
    rng = np.random.default_rng(5)
    x = rng.standard_normal(20)
    y = x + rng.standard_normal(20)
    parts = [("A", np.arange(20))]
    assert divisions.within_rho(x, y, parts)[0] == pytest.approx(
        spearmanr(x, y).statistic
    )


def test_gene_sharing_only_the_contrast_between_divisions_has_no_within_rho():
    """A gene that is the division's level plus noise: high whole-brain rho, ~0 within."""
    rng = np.random.default_rng(6)
    labels = np.repeat(["A", "B", "C", "D"], 40)
    level = np.repeat([0.0, 3.0, 6.0, 9.0], 40)
    map_values = level + rng.standard_normal(160)
    gene = level + rng.standard_normal(160)
    parts = divisions.division_parts(np.arange(160), labels)
    whole = spearmanr(map_values, gene).statistic
    within = divisions.within_rho(map_values, gene, parts)[0]
    assert whole > 0.8
    # four divisions of 40 independent structures: the SD of the mean rho is ~0.08
    assert abs(within) < 0.25


def test_contrast_leaves_out_genes_on_both_sides():
    """A gene in a set of each side of a contrast is in neither, and is reported."""
    members = {name: [] for name in gene_sets.SET_ORDER}
    members["other postsynaptic"] = ["Dlg4", "Slc32a1"]
    members["presynaptic"] = ["Syp", "Slc32a1"]
    first, second, both = gene_sets.contrast_sides(
        members, "postsynaptic against presynaptic"
    )
    assert first == ["Dlg4"]
    assert second == ["Syp"]
    assert both == ["Slc32a1"]


def test_set_median_null_takes_the_median_over_the_set_per_surrogate():
    """The null of a set's median is, per surrogate, the median of its genes' rho."""
    null = np.array([[0.1, 0.5, -0.2], [0.3, -0.1, 0.0], [0.2, 0.2, 0.4]])
    rho = pd.Series({"a": 0.5, "b": 0.6, "c": 0.7})
    row, null_median = gene_sets.set_test(
        ["a", "b", "c"], rho, null, {"a": 0, "b": 1, "c": 2}
    )
    assert np.allclose(null_median, [0.2, 0.2, 0.0])
    assert row["median_rho"] == pytest.approx(0.6)
    # no surrogate median reaches 0.6, so p is the floor of 3 surrogates
    assert row["p_spatial"] == pytest.approx(0.25)


def test_null_partial_equals_panel_tests_partial():
    """The vectorised partial rho of a surrogate equals ish.panel_test.partial."""
    rng = np.random.default_rng(7)
    composite = rng.standard_normal(50)
    gene = composite + rng.standard_normal(50)
    maps = rng.standard_normal((5, 50))
    fast = gene_sets.null_partial(maps, gene, composite)
    slow = [panel_test.partial(m, gene, [composite]) for m in maps]
    assert np.allclose(fast, slow)


@pytest.mark.skipif(
    not (TABLES / "gene_ranking.csv").exists(), reason="data not connected"
)
def test_todays_ranking_is_complete_and_consistent():
    """Today's ranking: both maps, every gene once per map, P9's 100 genes, q in 0..1."""
    ranking = gene_ranking.load_ranking()
    assert set(ranking["map"]) == {"nano", "auto"}
    for _, mine in ranking.groupby("map"):
        assert mine["symbol"].is_unique
        assert int(mine["p9_gene"].sum()) == 100
        assert mine["q_all"].between(0, 1).all()
        assert (mine["p_spatial"] > 0).all()
