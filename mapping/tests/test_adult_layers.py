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


def test_load_hierarchy_refuses_a_repeated_area(tmp_path):
    """A hierarchy table that lists an area twice stops the run."""
    path = tmp_path / "hierarchy.csv"
    pd.DataFrame(
        dict(
            acronym=["VISp", "AUDp", "VISp"],
            module=["Visual", "Auditory", "Visual"],
            hierarchy_score=[-0.4, -0.3, 0.1],
        )
    ).to_csv(path, index=False)
    with pytest.raises(ValueError, match="listed twice"):
        layers.load_hierarchy(path)


def test_check_hierarchy_refuses_a_missing_area():
    """A hierarchy table that misses an area of the atlas stops the run."""
    hier = pd.DataFrame({"acronym": ["VISp", "AUDp"]})
    layers.check_hierarchy(hier, {"VISp", "AUDp"})
    with pytest.raises(ValueError, match="missing"):
        layers.check_hierarchy(hier, {"VISp", "AUDp", "MOp"})


def test_exclusion_names_the_reason():
    """A cell without a value says whether it had no tissue, too little, or no signal;
    one with a value is left out only under the coverage it needs."""
    assert layers.exclusion(500, 0.3, 0.9, 0.25) == ""
    assert layers.exclusion(0, None, 0.0, 0.25) == "no tissue in the sections"
    assert layers.exclusion(100, None, 0.1, 0.25).startswith("100 tissue voxels")
    assert layers.exclusion(5000, None, 0.9, 0.25) == "mean nano at or below background"
    assert layers.exclusion(500, 0.3, 0.04, 0.25).startswith("4.0% of the area's")


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


def _cells(rows, kind="band"):
    """A per-mouse table from (depth, area, mouse, zref or None) tuples; a zref given
    as a string is a value left out for its coverage."""
    out = []
    for depth, area, mouse, z in rows:
        low = isinstance(z, str)
        reason = ""
        if z is None:
            reason = "no tissue in the sections"
        elif low:
            reason = "4.0% of the area's atlas voxels, under adult_layers.min_coverage"
        out.append(
            dict(
                reading="zref",
                depth_kind=kind,
                depth=depth,
                area=area,
                module="Visual",
                hierarchy_score=0.0,
                hierarchy_rank=1.0,
                mouse=mouse,
                group="naive",
                n_vox20=0 if z is None else 1000,
                atlas_vox20=25000 if low else 2000,
                coverage=0.0 if z is None else (0.04 if low else 0.5),
                zref=np.nan if z is None else float(z),
                excluded=reason != "",
                exclude_reason=reason,
            )
        )
    return pd.DataFrame(out, columns=layers.PER_MOUSE_COLUMNS)


def test_summary_table_gives_mean_sem_t_and_needs_enough_mice():
    """Mean, SD, SEM and t of the mice that count, a low-coverage value left out;
    under min_mice, no statistics, and the mice left out named."""
    table = _cells(
        [("granular", "X", f"m{i}", v) for i, v in enumerate([1.0, 2.0, 3.0, 4.0])]
        + [("granular", "X", "m4", "9.0")]
        + [("granular", "Y", "m0", 1.0), ("granular", "Y", "m1", 2.0)]
        + [("granular", "Y", "m2", None), ("granular", "Y", "m3", "3.0")]
    )
    s = layers.summary_table(table, min_mice=3).set_index("area")
    sd = np.std([1, 2, 3, 4], ddof=1)
    assert s.at["X", "n_mice"] == 4
    assert s.at["X", "mean"] == pytest.approx(2.5)
    assert s.at["X", "sem"] == pytest.approx(sd / 2)
    assert s.at["X", "t"] == pytest.approx(2.5 / (sd / 2))
    assert s.at["X", "mice_left_out"] == "m4"
    assert s.at["Y", "n_mice"] == 2
    assert np.isnan(s.at["Y", "mean"])
    assert s.at["Y", "mice_left_out"] == "m2;m3"


def test_contrasts_are_upper_minus_deeper_and_left_out_with_either():
    """A mouse's contrast is its upper minus its deeper cell; one missing either band
    has none, with that band's reason, and a low-coverage band keeps the difference
    but leaves the contrast out. L2/3 - L5 is taken from the layers."""
    bands = _cells(
        [
            ("supragranular", "X", "m0", 0.5),
            ("infragranular", "X", "m0", 0.2),
            ("supragranular", "X", "m1", 0.4),
            ("infragranular", "X", "m1", None),
            ("supragranular", "X", "m2", None),
            ("infragranular", "X", "m2", 0.1),
            ("supragranular", "X", "m3", "0.6"),
            ("infragranular", "X", "m3", 0.1),
        ]
    )
    layer_cells = _cells([("L2/3", "X", "m0", 0.3), ("L5", "X", "m0", 0.4)], kind="layer")
    out = layers.with_contrasts(pd.concat([bands, layer_cells], ignore_index=True))
    contrast = out[out["depth_kind"] == "contrast"]
    c = contrast[contrast["depth"] == "supragranular - infragranular"].set_index("mouse")
    assert c.at["m0", "zref"] == pytest.approx(0.3)
    assert not c.at["m0", "excluded"]
    assert c.at["m1", "excluded"] and np.isnan(c.at["m1", "zref"])
    assert c.at["m1", "exclude_reason"].startswith("infragranular: no tissue")
    assert c.at["m2", "excluded"] and np.isnan(c.at["m2", "zref"])
    assert c.at["m2", "exclude_reason"].startswith("supragranular: no tissue")
    assert c.at["m3", "excluded"]
    assert c.at["m3", "zref"] == pytest.approx(0.5)
    assert c.at["m3", "coverage"] == pytest.approx(0.04)
    lay = contrast[contrast["depth"] == "L2/3 - L5"].set_index("mouse")
    assert list(lay.index) == ["m0"]
    assert lay.at["m0", "zref"] == pytest.approx(-0.1)


def test_per_mouse_table_drops_undrawn_groups_and_marks_low_coverage():
    """A group with no voxel in the atlas gets no row; a cell under min_coverage keeps
    its zref but is left out, one under min_vox20 has none."""
    adults = layers.ADULTS
    key_ok = ("band", "infragranular", "A")
    key_low = ("band", "supragranular", "A")
    key_none = ("band", "granular", "A")
    groups = {key_ok: {1}, key_low: {2}, key_none: {3}}
    atlas_n = {key_ok: 1000, key_low: 10000, key_none: 0}

    # a cell holds voxels, mean sig, mean ratio, mean sepratio; sig 2 over a cortex
    # mean of 1, with median 0 and spread 1, is zref log2(2) = 1
    cells = {
        m: {key_ok: (800, 2.0, 1.0, 1.0), key_low: (800, 2.0, 1.0, 1.0)} for m in adults
    }
    cells[adults[0]][key_ok] = None
    n_vox = {m: {key_ok: 800, key_low: 800, key_none: 0} for m in adults}
    n_vox[adults[0]][key_ok] = 100
    norm = {m: (0.0, 1.0) for m in adults}
    refs = {m: {"cref": 1.0} for m in adults}
    hier = pd.DataFrame(
        dict(
            acronym=["A"], module=["Visual"], hierarchy_score=[0.0], hierarchy_rank=[1.0]
        )
    )
    t = layers.per_mouse_table(groups, cells, n_vox, norm, refs, hier, atlas_n, 0.25)
    assert set(t["depth"]) == {"supragranular", "infragranular"}
    assert len(t) == 2 * len(adults)
    ok = t[t["depth"] == "infragranular"].set_index("mouse")
    assert ok.loc[adults[1:], "zref"].tolist() == pytest.approx([1.0] * (len(adults) - 1))
    assert not ok.loc[adults[1:], "excluded"].any()
    assert ok.at[adults[0], "excluded"] and np.isnan(ok.at[adults[0], "zref"])
    assert ok.at[adults[0], "exclude_reason"].startswith("100 tissue voxels")
    low = t[t["depth"] == "supragranular"]
    assert low["excluded"].all()
    assert low["zref"].tolist() == pytest.approx([1.0] * len(adults))
    assert low["coverage"].tolist() == pytest.approx([0.08] * len(adults))


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

    # a row of the region table with no cell here stops the run too when it is one of
    # the table's mice, not when it is another brain's
    region.loc[1, "log2_value"] = -0.5
    extra = dict(reading="zref", division="Isocortex", acronym="Y", log2_value=0.2)
    other = pd.DataFrame([dict(extra, mouse="young0")])
    pd.concat([region, other]).to_csv(path, index=False)
    assert layers.check_against_region_table(table, path) == 2
    same = pd.DataFrame([dict(extra, mouse="m1")])
    pd.concat([region, same]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="no whole-area cell"):
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


def test_bars_grey_runs_from_black_to_light_grey_with_the_sem():
    """An SEM of 0 gives black; sem_max and beyond, or NaN, the lightest grey."""
    grey = layers_plotting.bars_grey(np.array([0.0, 0.05, 0.1, 0.3, np.nan]), 0.1)
    assert grey[:, 0] == pytest.approx([0.0, 0.39, 0.78, 0.78, 0.78])


def test_unscored_areas_sit_after_a_gap():
    """The scored areas fill the first positions; the unscored follow one step later."""
    x_of = layers_plotting.area_positions(layers.load_hierarchy())
    assert x_of["VISp"] == 0
    assert x_of["ORBvl"] == 36
    assert min(x_of[a] for a in ("AUDv", "SSp-un", "ECT", "GU", "PERI", "VISC")) == 38


@pytest.mark.skipif(not TABLE.is_file(), reason="data not connected")
def test_written_table_is_complete_and_consistent():
    """The written table: ten adults, 43 areas, no RSPd layer 4, a reason exactly where
    a cell is left out, a cell left out with a value only under min_coverage, and each
    contrast the difference of its two cells."""
    t = pd.read_csv(TABLE, keep_default_na=False, na_values=[""])
    assert t["mouse"].nunique() == 10
    assert t["area"].nunique() == 43
    assert not ((t["area"] == "RSPd") & (t["depth"].isin(["granular", "L4"]))).any()
    excluded = t["excluded"].astype(bool)
    assert (t.loc[excluded, "exclude_reason"] != "").all()
    assert t.loc[~excluded, "exclude_reason"].isna().all()
    assert t.loc[~excluded, "zref"].notna().all()
    low = excluded & t["zref"].notna() & (t["depth_kind"] != "contrast")
    assert (t.loc[low, "coverage"] < layers.ADULT_LAYERS["min_coverage"]).all()
    cells = t[t["depth_kind"].isin(["band", "layer"])].pivot_table(
        index=["area", "mouse"], columns="depth", values="zref"
    )

    # both read back at four decimals, so they agree to two units of the last
    for name, (_, upper), (_, deeper) in layers.CONTRASTS:
        contrast = t[(t["depth"] == name) & t["zref"].notna()]
        contrast = contrast.set_index(["area", "mouse"])
        expected = cells[upper] - cells[deeper]
        gap = (contrast["zref"] - expected.loc[contrast.index]).abs().max()
        assert gap < 2e-4
