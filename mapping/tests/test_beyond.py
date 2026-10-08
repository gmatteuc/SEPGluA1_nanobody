"""Known-answer checks of analysis 4 (beyond) and analysis 5 (the green channel)."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import rankdata, spearmanr

from sepmap import config
from sepmap.adult import beyond_calibration as bc
from sepmap.adult import beyond_density as bd
from sepmap.adult import sep_channel_check as scc

BEYOND = Path(config.DATA) / "adult_v2" / "ish_analysis" / "beyond"
N_ADULTS = len(bd.ADULTS)
SPLITS = bd.half_splits()


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
    _, explainable = bd.ceiling(adults, SPLITS)
    with_truth = spearmanr(bd.full_map(adults), truth).statistic ** 2
    # 500 structures: the agreement is known to about +-0.02
    assert abs(explainable - with_truth) < 0.03
    assert abs(explainable**2 - with_truth) > 0.05


def test_spearman_brown_steps_up_a_half_agreement():
    """Two halves agreeing at 0.974 make a ten-adult map of reliability 0.987."""
    assert bd.spearman_brown(0.974) == pytest.approx(0.98683, abs=1e-5)
    assert np.isnan(bd.spearman_brown(-1.0))


def test_implied_replication_reproduces_the_hand_worked_value():
    """h 0.974, in-sample R2 0.701, 13 terms, 126 structures imply 0.922."""
    assert bd.implied_replication(0.974, 0.701, 13, 126) == pytest.approx(0.922, abs=1e-3)
    assert np.isnan(bd.implied_replication(0.9, 0.95, 13, 126))


def test_a_map_made_of_its_predictors_leaves_almost_nothing():
    """A curved function of the predictors, bent, is predicted to within a few %."""
    xs, truth = predictors_and_truth()
    adults = cohort(truth, sd=0.1)
    _, explainable = bd.ceiling(adults, SPLITS)
    left = 1 - bd.cv_r2(bd.full_map(adults), bd.flexible(xs)) / explainable
    assert left < 0.08


def test_a_leftover_of_pure_animal_noise_does_not_replicate():
    """When the predictors make the whole map, two halves' leftovers are unrelated."""
    xs, _ = predictors_and_truth()
    # a map inside the span of the bent predictors, so all the leftover is noise
    z = [(x - x.mean()) / x.std() for x in xs]
    truth = z[0] + 0.5 * z[1] ** 2
    adults = cohort(truth, sd=1.0)
    agreement = bd.leftover_agreement(adults, bd.flexible(xs), SPLITS[:20])
    # 300 structures: two unrelated leftovers agree within about +-0.12
    assert abs(np.mean(agreement)) < 0.12


def test_a_planted_leftover_replicates_and_keeps_its_sign():
    """A part of the map the predictors cannot make replicates across half-cohorts."""
    xs, truth = predictors_and_truth()
    rng = np.random.default_rng(3)
    planted = 2.0 * rng.standard_normal(len(truth))
    adults = cohort(truth + planted, sd=0.5)
    model = bd.flexible(xs)
    assert np.mean(bd.leftover_agreement(adults, model, SPLITS[:20])) > 0.8
    res = bd.residual(bd.full_map(adults), model)
    same = bd.same_sign_share(adults, model, SPLITS[:20], res)
    big = np.abs(res) > np.percentile(np.abs(res), 90)
    assert same[big].min() > 0.9


def test_cross_validation_does_not_reward_noise_predictors():
    """Twenty random predictors fit noise in-sample and predict nothing held out."""
    rng = np.random.default_rng(4)
    y = rng.standard_normal(120)
    xs = [rng.standard_normal(120) for _ in range(20)]
    assert bd.r_squared(y, xs) > 0.1
    assert bd.cv_r2(y, xs) < 0.05


def test_the_model_bends_the_four_subunits_as_separate_terms():
    """The quoted model has x, x^2 and x^3 of seven predictors; April's of four."""
    covariates = {k: np.arange(10.0) for k in bd.SUBUNITS}
    covariates.update(
        abundance_composite=np.arange(10.0),
        markers=np.arange(10.0),
        psd_pc1=np.arange(10.0),
        autofluo=np.arange(10.0),
    )
    assert len(bd.model(covariates)) == 3 * 7
    assert len(bd.model(covariates, composite_abundance=True)) == 3 * 4
    assert len(bd.model(covariates, ("abundance",), bend=False)) == 4


def test_structure_rows_give_the_reason_a_structure_is_left_out():
    """A declared structure without a marker gene's value is left out, and says why."""
    set_table = pd.DataFrame(
        dict(
            structure=["kept", "no_gene", "outside"],
            acronym=["K", "N", "O"],
            division=["TH", "TH", "fiber tracts"],
            in_set=[True, True, False],
            reason=["", "", "not grey matter"],
        )
    )
    genes = bd.SUBUNITS + bd.MARKERS
    expr = {g: {"kept": 1.0, "no_gene": 1.0, "outside": 1.0} for g in genes}
    del expr[bd.MARKERS[0]]["no_gene"]
    rows = bd.structure_rows(set_table, expr, pd.Series(True, index=set_table.structure))
    by = rows.set_index("structure")
    assert bool(by.loc["kept", "used"])
    assert not by.loc["no_gene", "used"]
    assert by.loc["no_gene", "missing_genes"] == bd.MARKERS[0]
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
    a, b, split = bc.experiment_halves(per)
    assert split == ["three"]
    assert a["three"]["a"] > a["three"]["b"]
    assert b["three"]["a"] < b["three"]["b"]
    assert a["one"] == b["one"]


def test_made_up_adults_agree_between_halves_as_asked():
    """Noise sized from h gives two halves of five that agree at about h."""
    truth = np.random.default_rng(5).standard_normal(500)
    adults = bc.fake_cohort(truth, 0.9, np.random.default_rng(6))
    agreement, _ = bd.ceiling(adults, SPLITS[:30])
    # 500 structures: the mean agreement sits within about 0.02 of h
    assert abs(np.mean(agreement) - 0.9) < 0.03


def test_calibration_floor_grows_with_the_mismatch_of_the_predictors():
    """A known map leaves nothing with its own predictors, more with noisy copies."""
    xs, truth = predictors_and_truth()
    rng = np.random.default_rng(7)
    noisy = [rankdata(x + rng.normal(0, 60, len(x))) for x in xs]
    adults = bc.fake_cohort(truth, 0.97, np.random.default_rng(8))
    exact = bc.analyse(adults, bd.flexible(xs), SPLITS[:20])
    mismatched = bc.analyse(adults, bd.flexible(noisy), SPLITS[:20])
    assert exact["left"] < 0.08
    assert mismatched["left"] > exact["left"] + 0.1


def test_green_channel_rows_are_computed_on_the_structures_all_channels_have():
    """A structure missing one channel in an adult is left out of that adult's row."""
    rng = np.random.default_rng(9)
    structures = [f"s{i}" for i in range(60)]
    auto = rng.standard_normal(60)
    per = {}
    for mouse in scc.ADULTS:
        sep = auto + 0.1 * rng.standard_normal(60)
        nano = rng.standard_normal(60)
        per[mouse] = {
            "sig": dict(zip(structures, nano)),
            "auto": dict(zip(structures, auto)),
            "sep": dict(zip(structures[1:], sep[1:])),
        }
    profile = dict(zip(structures, auto))
    rows = pd.DataFrame(scc.channel_rows(per, profile))
    assert (rows["n_structures"] == 59).all()
    assert (rows["rho_sep_auto"] > 0.9).all()
    assert (rows["rho_sepresid_gria"].abs() < 0.5).all()


@pytest.mark.skipif(
    not (BEYOND / "variance_partition.csv").exists(), reason="analysis 4 has not run"
)
def test_todays_tables_hold_the_model_the_april_row_and_the_calibration():
    """The partition has the quoted model and April's; the calibration both maps."""
    partition = pd.read_csv(BEYOND / "variance_partition.csv").set_index("key")
    assert {"model", "april", "gria1", "density"} <= set(partition.index)
    model = partition.loc["model"]
    assert model["left"] == pytest.approx(1 - model["share_of_ceiling"])
    assert model["terms"] == 22
    calibration = pd.read_csv(BEYOND / "calibration.csv")
    counts = calibration["map"].value_counts()
    assert counts[bc.ABUNDANCE_DENSITY] == counts[bc.GRIA1] > 0
