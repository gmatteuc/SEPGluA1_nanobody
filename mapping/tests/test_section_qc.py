"""Section QC of an ISH grid: dim sections flagged and set missing, nothing filled in."""

import numpy as np
import pytest

from sepmap.ish import section_qc


def gradient_grid(n_ap=30, n_dv=10, n_ml=12, seed=0):
    """A grid whose energy rises tenfold along AP, with 5% noise, and its brain mask."""
    rng = np.random.default_rng(seed)
    level = np.geomspace(1.0, 10.0, n_ap)[:, None, None]
    vol = level * (1 + 0.05 * rng.standard_normal((n_ap, n_dv, n_ml)))
    brain = np.ones(vol.shape, dtype=bool)
    return vol.astype(np.float32), brain


def statuses(vol, brain, axis, kept=()):
    """The status of every section of `vol` along `axis`."""
    table = section_qc.flag_sections(
        section_qc.section_profile(vol, brain, axis), set(kept)
    )
    return table.set_index("section")["status"]


def test_a_tenfold_dimmer_section_is_flagged_and_set_missing():
    """A section planted ten times dimmer than its neighbours is flagged, then NaN."""
    vol, brain = gradient_grid()
    vol[15] /= 10
    status = statuses(vol, brain, 0)
    assert status[15] == "flagged"
    assert (status.drop(15) == "ok").all()
    clean = section_qc.apply_flags(vol, [15], 0)
    assert np.isnan(clean[15]).all()


def test_a_run_of_three_failed_sections_is_flagged_whole():
    """Three dim sections in a row are each flagged: a good section lies on each side."""
    vol, brain = gradient_grid()
    vol[14:17] /= 20
    status = statuses(vol, brain, 0)
    assert (status[[14, 15, 16]] == "flagged").all()
    assert (status.drop([14, 15, 16]) == "ok").all()


def test_a_step_in_expression_is_not_flagged():
    """A section that matches the sections behind a large step is anatomy, not a flag."""
    vol, brain = gradient_grid()
    vol[15:] /= 50
    status = statuses(vol, brain, 0)
    assert "flagged" not in set(status)


def test_sections_among_faint_neighbours_are_not_judged():
    """Among neighbours below ish_qc.min_reference, a drop is noise, not a flag."""
    vol, brain = gradient_grid()
    vol *= 0.01 / vol.max()
    vol[15] /= 10
    table = section_qc.flag_sections(section_qc.section_profile(vol, brain, 0), set())
    status = table.set_index("section")["status"]
    assert (status == "faint neighbours").all()


def test_a_smooth_gradient_is_not_flagged():
    """A tenfold gradient across the brain is anatomy, and no section is flagged."""
    vol, brain = gradient_grid()
    assert (statuses(vol, brain, 0) == "ok").all()


def test_a_listed_true_absence_is_kept():
    """A dim section the exceptions list keeps is marked absence, not set missing."""
    vol, brain = gradient_grid()
    vol[4] *= 0.01
    status = statuses(vol, brain, 0, kept=[4])
    assert status[4] == "absence kept"
    assert "flagged" not in set(status)


def test_a_sagittal_series_is_judged_along_ml():
    """A dim plane along ML is flagged when the axis is ML, and not along AP."""
    # the gradient grid turned so that its 30 sections run along ML
    vol, _ = gradient_grid()
    vol = np.transpose(vol, (2, 1, 0)).copy()
    brain = np.ones(vol.shape, dtype=bool)
    vol[:, :, 20] /= 10
    along_ml = statuses(vol, brain, 2)
    along_ap = statuses(vol, brain, 0)
    assert along_ml[20] == "flagged"
    assert "flagged" not in set(along_ap)


def test_nothing_is_interpolated():
    """Every section not flagged keeps its values to the bit; missing voxels stay NaN."""
    vol, brain = gradient_grid()
    vol[3, 2, 2] = np.nan
    clean = section_qc.apply_flags(vol, [10, 11], 0)
    keep = [k for k in range(vol.shape[0]) if k not in (10, 11)]
    assert np.array_equal(clean[keep], vol[keep], equal_nan=True)
    assert np.isnan(clean[[10, 11]]).all()
    assert vol[10].min() > 0


def test_end_sections_with_little_brain_are_not_judged():
    """Sections with fewer in-brain voxels than ish_qc.min_plane_voxels are not judged."""
    # a first section holding 4 in-brain voxels, all near zero
    vol, brain = gradient_grid()
    brain[0] = False
    brain[0, :2, :2] = True
    vol[0] = 0.001
    status = statuses(vol, brain, 0)
    assert status[0] == "too little brain"
    assert (status.drop(0) == "ok").all()


def test_a_section_without_data_is_not_judged():
    """A section whose voxels are all missing is marked no data, not flagged."""
    vol, brain = gradient_grid()
    vol[8] = np.nan
    status = statuses(vol, brain, 0)
    assert status[8] == "no data"
    assert (status.drop(8) == "ok").all()


def test_exceptions_file_statuses(tmp_path):
    """Proposed and accepted exceptions keep their sections, rejected ones none."""
    path = tmp_path / "exceptions.csv"
    path.write_text(
        "symbol,experiment_id,sections,reason,status\n"
        "Tac1,1038,4 5,anterior,proposed\n"
        "Glra1,7,15,forebrain,accepted\n"
        "Nrgn,736,44,looked failed,rejected\n",
        encoding="utf-8",
    )
    exceptions = section_qc.load_exceptions(path)
    assert exceptions[("Tac1", "1038")]["sections"] == {4, 5}
    assert exceptions[("Glra1", "7")]["sections"] == {15}
    assert exceptions[("Nrgn", "736")]["sections"] == set()


def test_unknown_exception_status_is_refused(tmp_path):
    """A status other than proposed, accepted or rejected stops the run."""
    path = tmp_path / "exceptions.csv"
    path.write_text(
        "symbol,experiment_id,sections,reason,status\nTac1,1038,4,x,maybe\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        section_qc.load_exceptions(path)
