"""Known answers for the adult map by depth (adult.layers, layers_plotting) and the
close-up helpers it shares (young_vs_adult.closeup)."""

import numpy as np
import pandas as pd
import pytest

from sepmap import config
from sepmap.adult import layers, layers_plotting
from sepmap.young_vs_adult import closeup

TABLE = config.DATA / "adult_v2" / "layers" / "area_layers_per_mouse.csv"


def test_layer_token_reads_both_namings():
    """A layer named with the word and one named without it are both read."""
    assert layers.layer_token("Primary visual area, layer 2/3") == "2/3"
    assert (
        layers.layer_token("Primary somatosensory area, barrel field, layer 6b") == "6b"
    )
    assert layers.layer_token("Anterior cingulate area, ventral part, 6a") == "6a"
    assert layers.layer_token("Primary visual area") is None
    assert layers.layer_token("Field CA1") is None


def test_depth_groups_pool_layers_into_bands_and_skip_absent_ones():
    """Layers go to their band, 6a and 6b to L6; an area without layer 4 has no
    granular group, and a label outside the isocortex is in none."""
    layer_names = ["1", "2/3", "4", "5", "6a", "6b"]
    stru, divi, sub = {}, {}, {}
    for i, lay in enumerate(layer_names, start=1):
        stru[i], divi[i], sub[i] = "A", "Isocortex", f"Area A, layer {lay}"
    for i, lay in enumerate(["1", "2/3", "5", "6a", "6b"], start=7):
        stru[i], divi[i], sub[i] = "B", "Isocortex", f"Area B, layer {lay}"
    stru[12], divi[12], sub[12] = "MOB", "OLF", "Main olfactory bulb, layer 1"
    groups = layers.depth_groups(stru, divi, sub)
    assert groups[("band", "supragranular", "A")] == {1, 2}
    assert groups[("band", "granular", "A")] == {3}
    assert groups[("band", "infragranular", "A")] == {4, 5, 6}
    assert groups[("layer", "L6", "B")] == {10, 11}
    assert groups[("whole", "whole", "A")] == {1, 2, 3, 4, 5, 6}
    assert ("band", "granular", "B") not in groups
    assert ("layer", "L4", "B") not in groups
    assert all(12 not in ids for ids in groups.values())


def test_hierarchy_table_holds_the_published_scores():
    """The table holds 43 areas, Harris's 37 scores as published (Supplementary
    Table 9, Cre_conf CC+TC+CT), VISp lowest and ORBvl highest, and the six
    unscored areas last."""
    hier = layers.load_hierarchy()
    assert len(hier) == 43
    scored = hier[hier["hierarchy_score"].notna()]
    assert len(scored) == 37
    score = hier.set_index("acronym")["hierarchy_score"]
    assert score["VISp"] == pytest.approx(-0.420983039)
    assert score["ORBvl"] == pytest.approx(0.386146821)
    assert score["VISrl"] == pytest.approx(-0.054969795)
    assert score["SSp-bfd"] == pytest.approx(-0.157092165)
    assert hier.iloc[0]["acronym"] == "VISp"
    assert hier.iloc[36]["acronym"] == "ORBvl"
    assert list(hier["order"]) == list(range(43))
    assert set(hier.iloc[37:]["acronym"]) == {
        "AUDv",
        "SSp-un",
        "ECT",
        "GU",
        "PERI",
        "VISC",
    }
    assert set(hier["module"]) == {
        "Visual",
        "Auditory",
        "Somatomotor",
        "Medial",
        "Lateral",
        "Prefrontal",
    }


def test_check_hierarchy_refuses_a_missing_area():
    """A hierarchy table that misses an area of the atlas stops the run."""
    hier = pd.DataFrame({"acronym": ["VISp", "AUDp"]})
    layers.check_hierarchy(hier, {"VISp", "AUDp"})
    with pytest.raises(ValueError, match="missing"):
        layers.check_hierarchy(hier, {"VISp", "AUDp", "MOp"})


def test_exclusion_names_the_reason():
    """A cell without a value says whether it had no tissue, too little, or no signal."""
    assert layers.exclusion(500, 0.3) == ""
    assert layers.exclusion(0, None) == "no tissue in the sections"
    assert layers.exclusion(100, None).startswith("100 tissue voxels")
    assert layers.exclusion(5000, None) == "mean nano at or below background"


def test_split_half_recovers_a_shared_profile_and_not_noise():
    """Ten mice sharing a profile under noise agree between halves near 1, over the
    126 cuts of ten into two fives; ten mice of noise alone do not."""
    rng = np.random.default_rng(1)
    profile = np.linspace(-1, 1, 30)
    shared = profile + rng.normal(0, 0.2, (10, 30))
    r, whole, n_cuts = layers.split_half(shared)
    assert n_cuts == 126
    assert r > 0.9
    assert whole == pytest.approx(2 * r / (1 + r))
    noise = rng.normal(0, 1, (10, 30))
    r_noise, _, _ = layers.split_half(noise)

    # mean of 126 overlapping cuts of 30 random areas: well under the 0.9 above
    assert abs(r_noise) < 0.3


def _cells(rows):
    """A per-mouse table from (depth, area, mouse, zref or None) tuples."""
    out = []
    for depth, area, mouse, z in rows:
        out.append(
            dict(
                reading="zref",
                depth_kind="band",
                depth=depth,
                area=area,
                module="Visual",
                hierarchy_score=0.0,
                hierarchy_rank=1.0,
                mouse=mouse,
                group="naive",
                n_vox20=1000 if z is not None else 0,
                atlas_vox20=2000,
                coverage=0.5 if z is not None else 0.0,
                zref=np.nan if z is None else z,
                excluded=z is None,
                exclude_reason="" if z is not None else "no tissue in the sections",
            )
        )
    return pd.DataFrame(out, columns=layers.PER_MOUSE_COLUMNS)


def test_summary_table_gives_mean_sem_t_and_needs_enough_mice():
    """Mean, SD, SEM and t of the mice with a value; under min_mice, no statistics,
    and the missing mice named."""
    table = _cells(
        [("granular", "X", f"m{i}", v) for i, v in enumerate([1.0, 2.0, 3.0, 4.0])]
        + [("granular", "Y", "m0", 1.0), ("granular", "Y", "m1", 2.0)]
        + [("granular", "Y", "m2", None), ("granular", "Y", "m3", None)]
    )
    s = layers.summary_table(table, min_mice=3).set_index("area")
    sd = np.std([1, 2, 3, 4], ddof=1)
    assert s.at["X", "n_mice"] == 4
    assert s.at["X", "mean"] == pytest.approx(2.5)
    assert s.at["X", "sem"] == pytest.approx(sd / 2)
    assert s.at["X", "t"] == pytest.approx(2.5 / (sd / 2))
    assert s.at["Y", "n_mice"] == 2
    assert np.isnan(s.at["Y", "mean"])
    assert s.at["Y", "mice_missing"] == "m2;m3"


def test_contrast_is_supra_minus_infra_and_missing_with_either():
    """The contrast of a mouse is its supragranular minus infragranular value; a mouse
    missing one band has none, with that band's reason."""
    table = _cells(
        [
            ("supragranular", "X", "m0", 0.5),
            ("infragranular", "X", "m0", 0.2),
            ("supragranular", "X", "m1", 0.4),
            ("infragranular", "X", "m1", None),
        ]
    )
    out = layers.with_contrast(table)
    c = out[out["depth_kind"] == "contrast"].set_index("mouse")
    assert c.at["m0", "zref"] == pytest.approx(0.3)
    assert not c.at["m0", "excluded"]
    assert c.at["m1", "excluded"]
    assert c.at["m1", "exclude_reason"].startswith("infragranular")


def test_check_against_region_table_refuses_a_difference(tmp_path):
    """Whole-area cells that match the region table to its four decimals pass; a
    difference in the third decimal, or a cell with no row there, stops the run."""
    table = _cells([("whole", "X", "m0", 0.12344), ("whole", "X", "m1", -0.5)])
    table["depth_kind"] = "whole"
    region = pd.DataFrame(
        dict(
            reading="zref",
            division="Isocortex",
            acronym=["X", "X"],
            mouse=["m0", "m1"],
            log2_value=[0.1234, -0.5],
        )
    )
    path = tmp_path / "region_means_per_mouse.csv"
    region.to_csv(path, index=False)
    assert layers.check_against_region_table(table, path) == 2
    region.loc[1, "log2_value"] = -0.502
    region.to_csv(path, index=False)
    with pytest.raises(ValueError, match="differs"):
        layers.check_against_region_table(table, path)
    region.drop(index=1).to_csv(path, index=False)
    with pytest.raises(ValueError, match="no row"):
        layers.check_against_region_table(table, path)


def test_band_stats_mean_sd_and_threshold():
    """Per pixel the mean and SD over the mice with a value; under min_n, none."""
    maps = np.full((5, 1, 1, 2), np.nan)
    maps[:, 0, 0, 0] = [1, 2, 3, 4, 5]
    maps[:2, 0, 0, 1] = [1, 3]
    mean, sd, n = layers.band_stats(maps, min_n=3)
    assert mean[0, 0, 0] == pytest.approx(3.0)
    assert sd[0, 0, 0] == pytest.approx(np.sqrt(2.5))
    assert n[0, 0, 1] == 2
    assert np.isnan(mean[0, 0, 1]) and np.isnan(sd[0, 0, 1])


def test_smooth_within_keeps_a_constant_inside_its_mask():
    """Mask-normalised smoothing leaves a constant constant up to the tissue's edge,
    and nothing outside it."""
    m = np.zeros((20, 20, 20), np.float32)
    m[5:15, 5:15, 5:15] = 1
    v = 2 * m
    out = closeup.smooth_within(v, m, [3.0, 1.0, 1.0])
    assert np.allclose(out[m > 0], 2.0, atol=1e-5)
    assert np.all(out[m == 0] == 0)


def test_band_average_is_over_tissue_only():
    """A band's mean counts only the bins with tissue; a band without any is NaN."""
    sv = np.zeros((1, 2, 4), np.float32)
    sm = np.zeros((1, 2, 4), np.float32)
    sv[0, 0, :2], sm[0, 0, :2] = [2.0, 0.0], [1.0, 0.0]
    out = closeup.band_average(sv, sm, 0, 2)
    assert out[0, 0] == pytest.approx(2.0)
    assert np.isnan(out[0, 1])


def test_bars_grey_runs_from_light_grey_to_black():
    """0 and NaN give the lightest grey, t_max and beyond give black."""
    grey = layers_plotting.bars_grey(np.array([0.0, np.nan, 5.0, 10.0, 30.0]), 10.0)
    assert grey[:, 0] == pytest.approx([0.78, 0.78, 0.39, 0.0, 0.0])


def test_unscored_areas_sit_after_a_gap():
    """The scored areas fill the first positions; the unscored follow one step later."""
    x_of = layers_plotting.area_positions(layers.load_hierarchy())
    assert x_of["VISp"] == 0
    assert x_of["ORBvl"] == 36
    assert min(x_of[a] for a in ("AUDv", "SSp-un", "ECT", "GU", "PERI", "VISC")) == 38


@pytest.mark.skipif(not TABLE.is_file(), reason="data not connected")
def test_written_table_is_complete_and_consistent():
    """The written table: ten adults, 43 areas, no RSPd layer 4, a reason exactly where
    a value is missing, and each contrast the difference of its two bands."""
    t = pd.read_csv(TABLE, keep_default_na=False, na_values=[""])
    assert t["mouse"].nunique() == 10
    assert t["area"].nunique() == 43
    assert not ((t["area"] == "RSPd") & (t["depth"].isin(["granular", "L4"]))).any()
    excluded = t["excluded"].astype(bool)
    assert (t.loc[excluded, "exclude_reason"] != "").all()
    assert t.loc[~excluded, "exclude_reason"].isna().all()
    bands = t[t["depth_kind"] == "band"].pivot_table(
        index=["area", "mouse"], columns="depth", values="zref"
    )
    contrast = t[(t["depth_kind"] == "contrast") & ~excluded].set_index(["area", "mouse"])
    expected = bands["supragranular"] - bands["infragranular"]

    # both read back at four decimals, so they agree to two units of the last
    gap = (contrast["zref"] - expected.loc[contrast.index]).abs().max()
    assert gap < 2e-4
