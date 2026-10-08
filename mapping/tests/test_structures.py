"""The declared structure set, zref with a declared reference, and the centroids."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from sepmap import structures

ADULTS = [f"m{i}" for i in range(10)]


def per_mouse_table():
    """Ten adults: a grey structure in all, one in nine, a fibre tract, a catch-all."""
    rows = []
    for i, mouse in enumerate(ADULTS):
        cases = [
            ("CA1", "CA1", "HPF"),
            ("Ventral posteromedial nucleus", "VPM", "TH"),
            ("corpus callosum", "cc", "lfbs"),
            ("Primary somatosensory area, unassigned", "SSp-un", "Isocortex"),
        ]
        if i > 0:
            cases.append(("Pontine gray", "PG", "P"))
        for name, acronym, division in cases:
            rows.append(
                dict(
                    mouse=mouse,
                    structure=name,
                    acronym=acronym,
                    division=division,
                    nano_mean=10.0 + i,
                )
            )
    return pd.DataFrame(rows)


def test_structure_set_keeps_grey_matter_seen_in_every_adult():
    """Grey structures measured in all ten adults are kept, the rest dropped, with why."""
    table = structures.structure_set(per_mouse_table(), ADULTS).set_index("structure")
    assert set(table.index[table["in_set"]]) == {"CA1", "Ventral posteromedial nucleus"}
    assert table.loc["Pontine gray", "reason"] == "measured in fewer than 10 adults"
    assert table.loc["Pontine gray", "n_adults"] == 9
    assert "not grey matter" in table.loc["corpus callosum", "reason"]
    assert table.loc["Primary somatosensory area, unassigned", "reason"].startswith(
        "catch-all"
    )


def test_structure_below_background_in_one_adult_is_not_measured_there():
    """A nano mean at or below background in one adult takes the structure out."""
    table = per_mouse_table()
    first_ca1 = (table["mouse"] == "m0") & (table["structure"] == "CA1")
    table.loc[first_ca1, "nano_mean"] = -0.5
    result = structures.structure_set(table, ADULTS).set_index("structure")
    assert not result.loc["CA1", "in_set"]
    assert result.loc["CA1", "n_adults"] == 9


def test_zref_has_median_zero_and_spread_one_over_its_reference():
    """Over the reference set, zref has median 0 and p90 - p10 of 1, in every brain."""
    rng = np.random.default_rng(0)
    for _ in range(5):
        v = pd.Series(rng.normal(0, 2, 60), index=[f"s{i}" for i in range(60)])
        reference = [f"s{i}" for i in range(40)]
        z, _, _ = structures.zref(v, reference)
        p10, median, p90 = np.percentile(z[reference], [10, 50, 90])
        assert median == pytest.approx(0.0, abs=1e-12)
        assert p90 - p10 == pytest.approx(1.0, abs=1e-12)


def test_structure_outside_the_reference_changes_no_zref():
    """Adding a structure outside the reference set moves no other structure's zref."""
    rng = np.random.default_rng(1)
    v = pd.Series(rng.normal(0, 1, 30), index=[f"s{i}" for i in range(30)])
    reference = list(v.index)
    before, _, _ = structures.zref(v, reference)
    extra = pd.concat([v, pd.Series({"brainstem extreme": 9.0})])
    after, _, _ = structures.zref(extra, reference)
    assert np.allclose(after[reference], before[reference], atol=0, rtol=0)


def symmetric_labels():
    """A label volume symmetric about the midline, with one bilateral structure."""
    labels = np.zeros((10, 12, 40), dtype=np.int32)

    # structure 1 on both sides, 8 voxels from the midline (ML 20) each way
    labels[3:7, 4:8, 10:14] = 1
    labels[3:7, 4:8, 26:30] = 1
    return labels, ["", "bilateral"]


def test_one_hemisphere_centroid_is_off_the_midline():
    """A bilateral structure sits in its hemisphere with one, on the midline with both."""
    labels, names = symmetric_labels()
    one = structures.centroids(labels, names, ["bilateral"]).iloc[0]
    both = structures.centroids(labels, names, ["bilateral"], one_hemisphere=False)
    midline_mm = 19.5 * structures.VOXEL_MM
    assert one["ml_mm"] == pytest.approx(11.5 * structures.VOXEL_MM)
    assert both.iloc[0]["ml_mm"] == pytest.approx(midline_mm)
    assert one["n_voxels"] == 4 * 4 * 4


def test_border_voxels_peel_one_voxel_off_each_structure():
    """A 5 x 5 x 5 cube inside another label keeps its 3 x 3 x 3 core after erosion."""
    labels = np.ones((9, 9, 9), dtype=np.int32)
    labels[2:7, 2:7, 2:7] = 2
    border = structures.border_voxels(labels)
    core = (labels == 2) & ~border
    assert core.sum() == 27
    assert core[3:6, 3:6, 3:6].all()


@pytest.mark.skipif(
    not Path(structures.STRUCTURE_SET).exists(),
    reason="run_structure_set has not written the set on this data root",
)
def test_todays_set_has_204_structures_all_with_a_centroid():
    """The declared set of the 5 October tables: 204 structures, each with a centroid."""
    table = structures.load_structure_set()
    declared = set(table.loc[table["in_set"], "structure"])
    centroids = structures.load_centroids()
    assert len(declared) == 204
    assert set(centroids.index) == declared
    assert centroids["ap_mm"].notna().all()
