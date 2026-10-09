"""Known-answer checks of the overview step: April's headline, the numbers, the index."""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import f_oneway

from sepmap.ish import april_headline, numbers, overview
from sepmap.ish.figure_index import (
    DRAWN_BY,
    FIGURE_NUMBERS,
    PARTS,
    QUESTIONS,
    SUPPLEMENTS,
    figure_file,
)


def test_f_statistic_equals_scipy_column_by_column():
    """The vectorised F of every column equals scipy's one-way ANOVA F."""
    rng = np.random.default_rng(0)
    groups = np.array(["a"] * 5 + ["b"] * 7 + ["c"] * 4)
    values = rng.standard_normal((len(groups), 30))
    values[groups == "b"] += 0.8
    f = april_headline.f_statistic(values, groups)
    for k in range(values.shape[1]):
        expected = f_oneway(*[values[groups == g, k] for g in "abc"]).statistic
        assert f[k] == pytest.approx(expected, rel=1e-10)


def test_hand_split_separates_the_forebrain_auxiliary_genes():
    """Cacng8 goes to Aux forebrain, Cacng5 to Aux other, a scaffold gene stays put."""
    assert april_headline.headline_group("Cacng8", "auxiliary") == "aux_forebrain"
    assert april_headline.headline_group("Cacng5", "auxiliary") == "aux_other"
    assert april_headline.headline_group("Dlg2", "scaffold") == "scaffold"


def test_planted_group_difference_passes_the_null_and_shared_noise_does_not():
    """A group planted high passes its null; rho drawn as the surrogates are does not."""
    rng = np.random.default_rng(1)
    keys = [k for k, _ in april_headline.HEADLINE_ORDER]
    rows = []
    for i, key in enumerate(keys):
        for j in range(6):
            rows.append(dict(symbol=f"g{i}_{j}", group=key, rho_today=0.0))
    headline = pd.DataFrame(rows)
    planted = headline["group"] == keys[0]
    headline.loc[planted, "rho_today"] = 0.8
    headline.loc[~planted, "rho_today"] = 0.05 * rng.standard_normal((~planted).sum())
    genes = list(headline["symbol"])
    null = 0.2 * rng.standard_normal((len(genes), 999))
    groups, f_info = april_headline.headline_null(headline, null, genes)
    assert groups.loc[0, "p_spatial"] < 0.01
    assert f_info["p_spatial"] < 0.01

    # rho drawn as the surrogates' are: the groups differ only by chance
    unrelated = headline.assign(rho_today=0.2 * rng.standard_normal(len(headline)))
    _, f_unrelated = april_headline.headline_null(unrelated, null, genes)
    assert f_unrelated["p_spatial"] > 0.05


def test_group_with_too_few_genes_is_not_tested():
    """A group under ish_analysis.min_set_genes genes gets no p and no q."""
    rng = np.random.default_rng(2)
    keys = [k for k, _ in april_headline.HEADLINE_ORDER]
    rows = [dict(symbol=f"a{j}", group=keys[0], rho_today=0.3) for j in range(3)]
    rows += [dict(symbol=f"b{j}", group=keys[1], rho_today=0.1) for j in range(8)]
    headline = pd.DataFrame(rows)
    null = rng.standard_normal((len(headline), 199))
    groups, _ = april_headline.headline_null(headline, null, list(headline["symbol"]))
    small = groups.set_index("group").loc[keys[0]]
    assert not small["tested"]
    assert np.isnan(small["p_spatial"]) and np.isnan(small["q"])


def test_numbers_keep_their_step_and_their_text(tmp_path):
    """Two steps' numbers come back in one table, values as written, the index first."""
    pd.DataFrame(
        dict(name=["a", "b"], value=["0.120", "Cacng8"], what=["x", "y"])
    ).to_csv(tmp_path / "numbers_one.csv", index=False)
    pd.DataFrame(dict(name=["a"], value=["7"], what=["z"])).to_csv(
        tmp_path / "numbers_two.csv", index=False
    )
    gathered = numbers.gather_numbers(tmp_path)
    n = numbers.lookup(gathered)
    assert n["one.a"] == "0.120"
    assert n["one.b"] == "Cacng8"
    assert n["two.a"] == "7"
    assert len(numbers.text_lines(gathered)) == 3 + 2 + 2


def test_the_gap_in_words_follows_the_two_sided_test():
    """The index reads the Cacng8 - Gria1 gap by its two-sided p, the shares after."""
    n = {
        "gene_ranking.gap_p": "0.105",
        "gene_ranking.gap_equal_first_as_large": "0.0232",
        "gene_ranking.gap_equal_second_as_large": "0.0815",
    }
    words = overview.gap_words(n)
    assert words.startswith("inside the two-sided test fixed in advance")
    assert "p 0.105; 2.3% of those maps give a Cacng8 lead as large, 8.2%" in words
    assert "band" not in words
    n["gene_ranking.gap_p"] = "0.01"
    assert overview.gap_words(n).startswith("past the two-sided test")


def test_ordinal_words():
    """1st, 2nd, 3rd, 4th, 11th, 12th, 13th, 21st, 22nd."""
    words = [numbers.ordinal(k) for k in (1, 2, 3, 4, 11, 12, 13, 21, 22)]
    assert words == ["1st", "2nd", "3rd", "4th", "11th", "12th", "13th", "21st", "22nd"]


def test_every_guided_figure_has_a_question_a_part_and_a_script():
    """The walk names every figure once, main ones in order, detailed ones beneath."""
    walked = [key for _, keys in PARTS for key in keys]
    details = [d for key in walked for d, _ in SUPPLEMENTS.get(key, ())]
    assert sorted(walked + details) == sorted(FIGURE_NUMBERS)
    assert [FIGURE_NUMBERS[k] for k in walked] == sorted(
        FIGURE_NUMBERS[k] for k in walked
    )
    assert all("s" not in FIGURE_NUMBERS[k] for k in walked)
    for key in walked:
        for detail, _ in SUPPLEMENTS.get(key, ()):
            assert FIGURE_NUMBERS[detail].startswith(FIGURE_NUMBERS[key] + "s")
    assert set(QUESTIONS) == set(FIGURE_NUMBERS) == set(DRAWN_BY)


@pytest.mark.skipif(not numbers.ALL_NUMBERS.exists(), reason="data not connected")
def test_written_index_and_overview_fill_from_the_numbers():
    """Every number the index and the overview figure name is in the run's numbers."""
    written = pd.read_csv(numbers.ALL_NUMBERS, dtype=str)
    n = numbers.lookup(written)
    text = overview.figure_index(n)
    for key in FIGURE_NUMBERS:
        assert figure_file(key) in text
    content = overview.overview_content(n)
    assert [len(part["numbers"]) for part in content["parts"]] == [3, 3]
    assert sum(len(figures) for _, figures in overview.figure_map()) == 16
