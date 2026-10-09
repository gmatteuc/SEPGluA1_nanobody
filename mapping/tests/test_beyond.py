"""Known-answer checks of analysis 4 (beyond) and analysis 5 (the green channel)."""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import rankdata, spearmanr

from sepmap.adult import (
    beyond_calibration,
    beyond_checks,
    beyond_density,
    profiles,
    sep_channel_check,
)
from sepmap.adult.profiles import ADULTS

BEYOND = beyond_density.OUT
N_ADULTS = len(ADULTS)
SPLITS = beyond_density.half_splits()


def cohort(truth: np.ndarray, sd: float, seed: int = 0) -> np.ndarray:
    """Ten made-up adults: the truth plus independent noise of SD `sd`."""
    rng = np.random.default_rng(seed)
    return truth + rng.normal(0, sd, size=(N_ADULTS, len(truth)))


def predictors_and_truth(n: int = 300, seed: int = 0):
    """Two predictors (ranks) and a map that is a curved function of them."""
    rng = np.random.default_rng(seed)
    a, b = rng.standard_normal(n), rng.standard_normal(n)
    xs = [rankdata(a), rankdata(b)]
    truth = np.tanh(a) + 0.5 * b**2
    return xs, truth


def test_ceiling_is_the_reliability_not_its_square():
    """The ceiling equals the full map's squared agreement with the truth, not less."""
    rng = np.random.default_rng(1)
    truth = rng.standard_normal(500)
    adults = cohort(truth, sd=1.5)
    _, explainable = beyond_density.ceiling(adults, SPLITS)
    with_truth = spearmanr(beyond_density.full_map(adults), truth).statistic ** 2
    # 500 structures: the agreement is known to about +-0.02
    assert abs(explainable - with_truth) < 0.03
    assert abs(explainable**2 - with_truth) > 0.05


def test_spearman_brown_steps_up_a_half_agreement():
    """Two halves agreeing at 0.974 make a ten-adult map of reliability 0.987."""
    assert beyond_density.spearman_brown(0.974) == pytest.approx(0.98683, abs=1e-5)
    assert np.isnan(beyond_density.spearman_brown(-1.0))


def test_implied_replication_reproduces_the_hand_worked_value():
    """h 0.974, in-sample R2 0.701, 13 terms, 126 structures imply 0.922."""
    assert beyond_density.implied_replication(0.974, 0.701, 13, 126) == pytest.approx(
        0.922, abs=1e-3
    )
    assert np.isnan(beyond_density.implied_replication(0.9, 0.95, 13, 126))


def test_a_map_made_of_its_predictors_leaves_almost_nothing():
    """A curved function of the predictors, bent, is predicted to within a few %."""
    xs, truth = predictors_and_truth()
    adults = cohort(truth, sd=0.1)
    _, explainable = beyond_density.ceiling(adults, SPLITS)
    left = (
        1
        - beyond_density.cv_r2(
            beyond_density.full_map(adults), beyond_density.flexible(xs)
        )
        / explainable
    )
    assert left < 0.08


def test_a_leftover_of_pure_animal_noise_does_not_replicate():
    """When the predictors make the whole map, two halves' leftovers are unrelated."""
    xs, _ = predictors_and_truth()
    # a map inside the span of the bent predictors, so all the leftover is noise
    z = [(x - x.mean()) / x.std() for x in xs]
    truth = z[0] + 0.5 * z[1] ** 2
    adults = cohort(truth, sd=1.0)
    agreement = beyond_density.leftover_agreement(
        adults, beyond_density.flexible(xs), SPLITS[:20]
    )
    # 300 structures: two unrelated leftovers agree within about +-0.12
    assert abs(np.mean(agreement)) < 0.12


def test_a_planted_leftover_replicates_and_keeps_its_sign():
    """A part of the map the predictors cannot make replicates across half-cohorts."""
    xs, truth = predictors_and_truth()
    rng = np.random.default_rng(3)
    planted = 2.0 * rng.standard_normal(len(truth))
    adults = cohort(truth + planted, sd=0.5)
    model = beyond_density.flexible(xs)
    assert np.mean(beyond_density.leftover_agreement(adults, model, SPLITS[:20])) > 0.8
    res = beyond_density.residual(beyond_density.full_map(adults), model)
    same = beyond_density.same_sign_share(adults, model, SPLITS[:20], res)
    big = np.abs(res) > np.percentile(np.abs(res), 90)
    assert same[big].min() > 0.9


def test_cross_validation_does_not_reward_noise_predictors():
    """Twenty random predictors fit noise in-sample and predict nothing held out."""
    rng = np.random.default_rng(4)
    y = rng.standard_normal(120)
    xs = [rng.standard_normal(120) for _ in range(20)]
    assert beyond_density.r_squared(y, xs) > 0.1
    assert beyond_density.cv_r2(y, xs) < 0.05


def test_the_main_model_is_gria1_and_the_density_term_straight():
    """Two columns; six curved; five with the four subunits; three with autofluo."""
    names = beyond_density.SUBUNITS + ("density", "psd_pc1", "autofluo", "psd95")
    covariates = {k: np.arange(10.0) for k in names}
    main = beyond_density.model_terms()
    assert main == {"abundance": ("Gria1",), "density": ("density",)}
    assert len(beyond_density.model_columns(covariates, main)) == 2
    assert len(beyond_density.model_columns(covariates, main, bend=True)) == 6
    subunits = beyond_density.model_terms(abundance=beyond_density.SUBUNITS)
    assert len(beyond_density.model_columns(covariates, subunits)) == 5
    auto = beyond_density.model_terms(autofluorescence=True)
    assert len(beyond_density.model_columns(covariates, auto)) == 3
    assert len(beyond_density.model_columns(covariates, main, ("abundance",))) == 1


def test_a_model_without_its_predictor_on_these_structures_stops():
    """A term build_covariates could not make (a gene not measured) stops the run."""
    covariates = {"Gria1": np.arange(10.0)}
    with pytest.raises(ValueError, match="density"):
        beyond_density.model_columns(covariates, beyond_density.model_terms())


def test_the_model_genes_are_gria1_and_the_genes_of_its_composites():
    """The main model needs Gria1 and the chosen genes; psd_pc1 and PSD95 need none."""
    main = beyond_density.model_terms()
    assert beyond_density.model_genes(main) == ("Gria1",) + beyond_density.DENSITY_GENES
    first = beyond_density.model_terms(("first_proposal",))
    assert beyond_density.model_genes(first) == ("Gria1",) + beyond_density.FIRST_PROPOSAL
    extra = beyond_density.model_terms(("density", "psd_pc1"))
    assert beyond_density.model_genes(extra) == beyond_density.model_genes(main)
    assert beyond_density.model_genes(beyond_density.model_terms(("psd95",))) == (
        "Gria1",
    )


def test_the_four_parts_add_to_one_and_a_negative_one_is_kept():
    """Gria1 only, shared, density only and left from three shares; signs as computed."""
    parts = beyond_density.partition(dict(abundance=0.5, density=0.6, model=0.7))
    assert parts["gria1_only"] == pytest.approx(0.1)
    assert parts["density_only"] == pytest.approx(0.2)
    assert parts["shared"] == pytest.approx(0.4)
    assert parts["left"] == pytest.approx(0.3)
    assert sum(parts.values()) == pytest.approx(1.0)
    # a term that adds nothing held out costs a little: its own part below zero
    parts = beyond_density.partition(dict(abundance=0.5, density=0.02, model=0.49))
    assert parts["density_only"] == pytest.approx(-0.01)
    assert sum(parts.values()) == pytest.approx(1.0)


def test_a_map_of_gria1_alone_leaves_density_nothing_of_its_own():
    """When the map is Gria1 plus noise, density's own part is near zero."""
    rng = np.random.default_rng(15)
    n = 300
    gria1 = rng.standard_normal(n)
    density = 0.5 * gria1 + rng.standard_normal(n)
    adults = cohort(gria1, sd=0.5)
    covariates = {"Gria1": rankdata(gria1), "density": rankdata(density)}
    y = beyond_density.full_map(adults)
    _, explainable = beyond_density.ceiling(adults, SPLITS[:20])
    shares = beyond_density.model_shares(
        y, covariates, beyond_density.model_terms(), explainable
    )
    parts = beyond_density.partition(shares)
    # 300 structures: a held-out share is known to within a few points
    assert abs(parts["density_only"]) < 0.03
    assert parts["gria1_only"] > 0.2


def test_the_weights_recover_planted_coefficients():
    """A map 0.8 Gria1 + 0.4 density on ranks gives weights in that ratio."""
    rng = np.random.default_rng(16)
    n = 400
    a, b = rng.standard_normal(n), rng.standard_normal(n)
    covariates = {"Gria1": rankdata(a), "density": rankdata(b)}
    z = [(x - x.mean()) / x.std() for x in covariates.values()]
    y = 0.8 * z[0] + 0.4 * z[1] + 0.1 * rng.standard_normal(n)
    found = beyond_density.weights(y, covariates, beyond_density.model_terms())
    assert list(found) == ["Gria1", "density"]
    # the map's SD is about 0.9, so the z-scored weights are 0.8 / 0.9 and 0.4 / 0.9
    assert found["Gria1"] / found["density"] == pytest.approx(2.0, rel=0.1)
    assert found["Gria1"] == pytest.approx(
        0.8 / np.sqrt(0.8**2 + 0.4**2 + 0.01), rel=0.05
    )


def test_a_composite_is_a_rank_and_is_made_only_where_every_gene_is_measured():
    """psd95 and the density term are ranks; sap102 and a composite with a gap absent."""
    structures = ["a", "b", "c", "d"]
    genes = beyond_density.SUBUNITS + beyond_density.DENSITY_GENES + ("Psd1", "Psd2")
    rng = np.random.default_rng(14)
    expr = {g: dict(zip(structures, rng.standard_normal(4))) for g in genes}
    role = {"Psd1": beyond_density.PSD_ROLE, "Psd2": beyond_density.PSD_ROLE}
    synapses = pd.DataFrame(
        {"psd95": [0.3, 0.1, 0.4, 0.2], "sap102": [0.1, np.nan, 0.2, 0.3]},
        index=structures,
    )
    auto = rng.standard_normal((N_ADULTS, 4))
    covariates, psd, _ = beyond_density.build_covariates(
        expr, role, auto, structures, synapses=synapses
    )
    assert psd == ["Psd1", "Psd2"]
    assert list(covariates["psd95"]) == [3.0, 1.0, 4.0, 2.0]
    # ranks of four structures, ties averaged, add to 1 + 2 + 3 + 4
    assert covariates["density"].sum() == pytest.approx(10.0)
    assert "sap102" not in covariates
    assert "first_proposal" not in covariates


def test_the_genes_measured_once_are_those_of_the_model():
    """A density gene with one experiment is listed, a subunit outside the model not."""
    once = beyond_density.DENSITY_GENES[0]
    split = [
        g
        for g in beyond_density.SUBUNITS + beyond_density.DENSITY_GENES
        if g not in ("Gria4", once)
    ]
    assert beyond_calibration.measured_once(split, beyond_density.model_terms()) == [once]
    assert beyond_calibration.measured_once(
        split, beyond_density.model_terms(("psd95",))
    ) == ["psd95"]


def test_structure_rows_give_the_reason_a_structure_is_left_out():
    """A declared structure without a density gene's value is left out, and says why."""
    set_table = pd.DataFrame(
        dict(
            structure=["kept", "no_gene", "outside"],
            acronym=["K", "N", "O"],
            division=["TH", "TH", "fiber tracts"],
            in_set=[True, True, False],
            reason=["", "", "not grey matter"],
        )
    )
    genes = beyond_density.model_genes(beyond_density.model_terms())
    expr = {g: {"kept": 1.0, "no_gene": 1.0, "outside": 1.0} for g in genes}
    missing = beyond_density.DENSITY_GENES[1]
    del expr[missing]["no_gene"]
    rows = beyond_density.structure_rows(set_table, expr)
    by = rows.set_index("structure")
    assert bool(by.loc["kept", "used"])
    assert not by.loc["no_gene", "used"]
    assert by.loc["no_gene", "missing_genes"] == missing
    assert by.loc["no_gene", "reason"] == f"not measured: {missing}"
    assert by.loc["outside", "reason"].startswith("not in the declared set")


def test_experiment_halves_alternate_by_id_and_share_a_single_experiment():
    """Experiments 1, 2, 10 split into 1 and 10 against 2; one experiment is in both."""
    per = {
        "three": {
            "10": {"a": 3.0, "b": 1.0},
            "2": {"a": 1.0, "b": 3.0},
            "1": {"a": 2.0, "b": 1.0},
        },
        "one": {"5": {"a": 1.0, "b": 2.0}},
    }
    a, b, split = beyond_calibration.experiment_halves(per)
    assert split == ["three"]
    assert a["three"]["a"] > a["three"]["b"]
    assert b["three"]["a"] < b["three"]["b"]
    assert a["one"] == b["one"]


def test_made_up_adults_agree_between_halves_as_asked():
    """Noise sized from h gives two halves of five that agree at about h."""
    truth = np.random.default_rng(5).standard_normal(500)
    adults = beyond_calibration.fake_cohort(truth, 0.9, np.random.default_rng(6))
    agreement, _ = beyond_density.ceiling(adults, SPLITS[:30])
    # 500 structures: the mean agreement sits within about 0.02 of h
    assert abs(np.mean(agreement) - 0.9) < 0.03


def test_calibration_floor_grows_with_the_mismatch_of_the_predictors():
    """A known map leaves nothing with its own predictors, more with noisy copies."""
    xs, truth = predictors_and_truth()
    rng = np.random.default_rng(7)
    noisy = [rankdata(x + rng.normal(0, 60, len(x))) for x in xs]
    adults = beyond_calibration.fake_cohort(truth, 0.97, np.random.default_rng(8))
    exact = beyond_calibration.analyse(adults, beyond_density.flexible(xs), SPLITS[:20])
    mismatched = beyond_calibration.analyse(
        adults, beyond_density.flexible(noisy), SPLITS[:20]
    )
    assert exact["left"] < 0.08
    assert mismatched["left"] > exact["left"] + 0.1


def test_green_channel_rows_are_computed_on_the_structures_all_channels_have():
    """A structure missing one channel in an adult is left out of that adult's row."""
    rng = np.random.default_rng(9)
    structures = [f"s{i}" for i in range(60)]
    auto = rng.standard_normal(60)
    per = {}
    for mouse in ADULTS:
        sep = auto + 0.1 * rng.standard_normal(60)
        nano = rng.standard_normal(60)
        per[mouse] = {
            "sig": dict(zip(structures, nano)),
            "auto": dict(zip(structures, auto)),
            "sep": dict(zip(structures[1:], sep[1:])),
        }
    profile = dict(zip(structures, auto))
    rows = pd.DataFrame(sep_channel_check.channel_rows(per, profile))
    assert (rows["n_structures"] == 59).all()
    assert (rows["rho_sep_auto"] > 0.9).all()
    assert (rows["rho_sepresid_gria"].abs() < 0.5).all()


def test_the_first_shuffling_is_one_seeded_permutation_dealt_in_turn():
    """fold_labels' first shuffling deals a seeded permutation into the folds."""
    n = 126
    order = np.random.default_rng(0).permutation(n)
    label = beyond_density.fold_labels(n, repeats=3)[0]
    for k in range(5):
        assert set(np.nonzero(label == k)[0]) == set(order[k::5])


def test_spatial_blocks_keep_neighbours_in_one_fold():
    """Structures of one k-means block always share a fold, and every fold is used."""
    rng = np.random.default_rng(10)
    centres = rng.uniform(0, 10, (20, 3))
    xyz = np.repeat(centres, 6, axis=0) + 0.05 * rng.standard_normal((120, 3))
    for label in beyond_density.block_labels(xyz, n_blocks=20, repeats=3):
        assert len(np.unique(label)) == 5
        for i in range(20):
            assert len(np.unique(label[6 * i : 6 * i + 6])) == 1


def test_repeated_folds_score_a_known_model_as_one_shuffling_does_on_average():
    """The mean over shufflings sits inside the spread of single shufflings."""
    xs, truth = predictors_and_truth(n=150, seed=11)
    y = rankdata(truth + np.random.default_rng(12).normal(0, 0.5, len(truth)))
    single = [
        beyond_density.cv_r2(
            y,
            beyond_density.flexible(xs),
            beyond_density.fold_labels(len(y), repeats=1, seed=s),
        )
        for s in range(30)
    ]
    repeated = beyond_density.cv_r2(
        y, beyond_density.flexible(xs), beyond_density.fold_labels(len(y), repeats=30)
    )
    assert min(single) <= repeated <= max(single)


def test_jackknife_sd_of_a_mean_matches_its_standard_error():
    """For a mean, the delete-d jackknife SD is close to SD / sqrt(n)."""
    rng = np.random.default_rng(13)
    x = rng.standard_normal(100)
    d = 20
    values = [np.delete(x, rng.choice(100, d, replace=False)).mean() for _ in range(2000)]
    # 2000 subsamples: the jackknife SD is known to within about 3%
    assert beyond_calibration.jackknife_sd(values, 100, d) == pytest.approx(
        x.std(ddof=1) / 10, rel=0.08
    )


def test_nano_on_the_allen_grid_is_the_mean_of_each_voxels_tissue():
    """A 200 um voxel takes the mean of its tissue; under half tissue it has none."""
    values = np.zeros((20, 10, 10), dtype=np.float32)
    tissue = np.zeros((20, 10, 10), dtype=bool)
    values[:10] = 2.0
    values[:10, :5] = 4.0
    tissue[:10] = True
    tissue[10:, :4] = True
    means = profiles.grid_means(values, tissue, 0.5)
    assert means.shape == (2, 1, 1)
    assert means[0, 0, 0] == pytest.approx(3.0)
    assert np.isnan(means[1, 0, 0])
    with pytest.raises(ValueError):
        profiles.grid_means(values[:15], tissue[:15], 0.5)


@pytest.mark.skipif(
    not (BEYOND / "variance_partition.csv").exists(), reason="data not connected"
)
def test_written_tables_hold_the_main_model_its_check_rows_and_the_floor():
    """Two terms; four parts adding to one; every check row with a floor; no benchmark."""
    partition, parts, weights = beyond_density.load_partition()
    assert {"model", "abundance", "density", "curved"} <= set(partition.index)
    model = partition.loc["model"]
    assert model["left"] == pytest.approx(1 - model["share_of_ceiling"])
    assert model["terms"] == 3
    assert parts["share"].sum() == pytest.approx(1.0)
    assert parts.loc["left", "share"] == pytest.approx(model["left"])
    assert list(weights.index) == ["Gria1", "density"]
    folds = pd.read_csv(BEYOND / "folds.csv").set_index("key")
    assert folds.loc["main", "left"] == pytest.approx(model["left"])
    checks = beyond_checks.load_check_rows().set_index("key")
    assert list(checks.index) == list(beyond_density.BEYOND["checks"])
    assert checks["nano_minus_floor"].notna().all()
    used = pd.read_csv(BEYOND / "structures_used.csv")
    assert checks.loc["curved", "n_structures"] == int(used["used"].sum())
    calibration = beyond_calibration.load_calibration()
    assert set(calibration["map"]) == {
        beyond_calibration.NANO,
        beyond_calibration.FLOOR_MAP,
    }
