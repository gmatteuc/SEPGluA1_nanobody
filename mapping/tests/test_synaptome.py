"""The measured synapse density: reading the source, placing samples, per structure."""

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from sepmap.adult import synaptome as sy

# a small ontology: the root, a division, structure A with two layers, structure B, and
# a region P above structures C and D, with a layer P1 that no structure holds
ONTOLOGY = pd.DataFrame(
    dict(
        id=[1, 2, 10, 11, 12, 20, 30, 31, 32, 33],
        acronym=["root", "DIV", "A", "A1", "A2", "B", "P", "C", "D", "P1"],
        name=["root", "div", "a", "a1", "a2", "b", "p", "c", "d", "p1"],
        parent=[sy.NO_PARENT, 1, 2, 10, 10, 2, 2, 30, 30, 30],
    )
).set_index("id")
STRUCTURES = {10: "A", 20: "B", 31: "C", 32: "D"}
PARENT = ONTOLOGY["parent"].to_dict()
ACRONYM = ONTOLOGY["acronym"].to_dict()
CHILDREN = {}
for _i, _p in PARENT.items():
    CHILDREN.setdefault(_p, []).append(_i)


def samples(rows):
    """A sample table from (name, hemisphere, source_id, name_acronym, psd95) rows."""
    table = pd.DataFrame(
        rows, columns=["name", "hemisphere", "source_id", "name_acronym", "psd95"]
    )
    table["measured"] = True
    for measure in sy.MEASURES:
        if measure != "psd95":
            table[measure] = table["psd95"]
    return table


def test_git_blob_hash_is_the_one_git_gives():
    """The hash of a file's bytes equals git hash-object's (known answers)."""
    assert sy.git_blob_sha1(b"") == "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391"
    assert (
        sy.git_blob_sha1(b"hello world\n") == "3b18e512dba79e4c8300dd08aeb37f8e728b8dad"
    )


def write_xlsx(path):
    """A two-sheet workbook whose parts are listed out of order, with every cell kind."""
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rel_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    book = (
        f'<workbook {ns} xmlns:r="{rel_ns}"><sheets>'
        '<sheet name="first" sheetId="1" r:id="rId7"/>'
        '<sheet name="second" sheetId="2" r:id="rId3"/></sheets></workbook>'
    )
    rels = (
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f'<Relationship Id="rId3" Type="{rel_ns}/worksheet" '
        'Target="worksheets/sheet1.xml"/>'
        f'<Relationship Id="rId7" Type="{rel_ns}/worksheet" '
        'Target="/xl/worksheets/sheet2.xml"/></Relationships>'
    )
    strings = f"<sst {ns}><si><t>Type</t></si><si><t>VPM left</t></si></sst>"
    first = (
        f"<worksheet {ns}><sheetData>"
        '<row r="1"><c r="A1" t="s"><v>0</v></c>'
        '<c r="C1" t="inlineStr"><is><t>x</t></is></c>'
        '</row><row r="2"><c r="A2"><v>1</v></c><c r="B2"><v>0.25</v></c></row>'
        "</sheetData></worksheet>"
    )
    second = (
        f'<worksheet {ns}><sheetData><row r="1"><c r="A1" t="s"><v>1</v></c></row>'
        "</sheetData></worksheet>"
    )
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("xl/workbook.xml", book)
        z.writestr("xl/_rels/workbook.xml.rels", rels)
        z.writestr("xl/sharedStrings.xml", strings)
        z.writestr("xl/worksheets/sheet2.xml", first)
        z.writestr("xl/worksheets/sheet1.xml", second)


def test_xlsx_is_read_in_the_workbooks_order_with_every_cell_kind(tmp_path):
    """Sheets come in the workbook's order, whatever their part names; empty is None."""
    path = tmp_path / "book.xlsx"
    write_xlsx(path)
    first, second = sy.read_xlsx(path)
    assert first == [["Type", None, "x"], [1.0, 0.25, None]]
    assert second == [["VPM left"]]


def test_column_letters_count_from_zero():
    """A is column 0, Z 25, AA 26, ACV 775 (the last column of the density sheet)."""
    assert [sy.column_index(r) for r in ("A1", "Z9", "AA2", "ACV38")] == [0, 25, 26, 775]


def test_names_give_the_acronym_the_ontology_writes():
    """Main-list names give their region; cortical ones add the layer, as the source."""
    cases = {
        "'main_Thalamus_VPM left'": "VPM",
        "'main_Hippocampus_DG-mo inf left'": "DG-mo",
        "'main_Cerebellum_CUL4,5gr inf right'": "CUL4, 5gr",
        "'ctxall_Isocortex_layers_Layer2-3 left_intersect_Isocortex_AId left'": "AId2/3",
        "'ctxall_Isocortex_Mop left_intersect_Isocortex_layers_Layer6 left'": "Mop6a",
        "'ctxall_Cortex_layers_Layer4 ENTl left_intersect_Cortex_ENTl left'": "ENTl4",
    }
    for name, acronym in cases.items():
        assert sy.name_acronym(name) == acronym


def test_hemisphere_is_read_from_the_end_of_the_name():
    """A trailing space is ignored; a name without a hemisphere is refused."""
    assert sy.hemisphere("'main_Cerebellum_IP left '") == "left"
    assert sy.hemisphere("'main_Thalamus_VPM right'") == "right"
    with pytest.raises(ValueError):
        sy.hemisphere("'main_Thalamus_VPM'")


def test_measures_follow_the_subtype_settings():
    """PSD95 is 30 subtypes, SAP102 26, and the three kinds are 11, 7 and 19."""
    assert len(sy.MEASURES["psd95"]) == 30
    assert len(sy.MEASURES["sap102"]) == 26
    assert len(sy.MEASURES["all_puncta"]) == 37
    assert (len(sy.PSD95_ONLY), len(sy.SAP102_ONLY), len(sy.BOTH)) == (11, 7, 19)
    assert sy.PSD95_ONLY == tuple(range(11))


def test_a_sample_without_any_punctum_is_not_measured():
    """All 37 subtypes at zero is a sample not measured; densities are subtype means."""
    density = np.zeros((37, 3))
    density[:, 0] = 1.0
    density[list(sy.PSD95_ONLY), 1] = 1.0
    info = pd.DataFrame(
        dict(
            Region_list=["'main_Thalamus_VPM left'"] * 3,
            ara_id=[733.0, 733.0, np.nan],
            acronym=["VPM", "VPM", np.nan],
        )
    )
    table = sy.sample_table(density, info)
    assert list(table["measured"]) == [True, True, False]
    assert table.loc[1, "psd95"] == pytest.approx(11 / 30)
    assert table.loc[1, "psd95_only"] == 1.0
    assert table.loc[1, "sap102"] == 0.0


def test_check_source_refuses_subtypes_not_scaled_to_0_1():
    """A subtype that does not run from 0 to 1, or blocks out of order, stops the run."""
    density = np.tile(np.linspace(0, 1, 5), (37, 1))
    names = [f"s{i}" for i in range(5)]
    info = pd.DataFrame(dict(Region_list=names))
    order = np.arange(37)
    sy.check_source(density, names, info, order)
    unscaled = density.copy()
    unscaled[3] *= 0.5
    with pytest.raises(ValueError):
        sy.check_source(unscaled, names, info, order)
    mixed = order.copy()
    mixed[[10, 11]] = mixed[[11, 10]]
    with pytest.raises(ValueError):
        sy.check_source(density, names, info, mixed)


def test_a_layer_is_placed_in_its_structure_and_a_structure_in_itself():
    """A layer lies in its structure; a structure in itself."""
    assert sy.place(11, PARENT, CHILDREN, STRUCTURES, ACRONYM) == (10, -1, "")
    assert sy.place(20, PARENT, CHILDREN, STRUCTURES, ACRONYM) == (20, -1, "")


def test_a_region_above_structures_is_never_spread_onto_them():
    """A parent of two structures, and a layer of it, stay unplaced, standing for it."""
    sid, above, reason = sy.place(30, PARENT, CHILDREN, STRUCTURES, ACRONYM)
    assert (sid, above) == (-1, 30)
    assert reason.startswith(sy.REASON_ABOVE) and "holds 2" in reason
    sid, above, reason = sy.place(33, PARENT, CHILDREN, STRUCTURES, ACRONYM)
    assert (sid, above) == (-1, 30)
    assert reason.startswith(sy.REASON_NO_STRUCTURE)


def test_samples_are_matched_by_id_first_and_acronym_second():
    """The source's id wins; without one the acronym is looked up; else no match."""
    table = samples(
        [
            ("a", "left", 11.0, "B", 0.5),
            ("b", "left", np.nan, "B", 0.5),
            ("c", "left", np.nan, "X9", 0.5),
            ("d", "left", 30.0, "P", 0.5),
        ]
    )
    out = sy.match_samples(table, ONTOLOGY, STRUCTURES)
    assert list(out["matched_by"]) == ["id", "acronym", "", "id"]
    assert list(out["structure"]) == ["A", "B", "", ""]
    assert list(out["used"]) == [True, True, False, False]
    assert out.loc[2, "reason"].startswith(sy.REASON_NOT_CCF)
    assert out.loc[3, "above_id"] == 30


def test_a_sample_not_measured_is_dropped_even_with_an_id():
    """No punctum at all is not measured, whatever id the source gives."""
    table = samples([("a", "right", 11.0, "A1", 0.0)])
    table["measured"] = False
    out = sy.match_samples(table, ONTOLOGY, STRUCTURES)
    assert not out.loc[0, "used"]
    assert out.loc[0, "reason"].startswith(sy.REASON_NOT_MEASURED)


def placed(rows, voxels):
    """Per-structure densities of samples (name, hemisphere, id, acronym, psd95)."""
    out = sy.match_samples(samples(rows), ONTOLOGY, STRUCTURES)
    return sy.structure_density(out, voxels)


def test_units_are_weighted_by_their_voxels():
    """Layer A1 (3 voxels) at 1 and A2 (1 voxel) at 0 give A 0.75, all of A covered."""
    rows = [
        ("a1", "left", 11.0, "A1", 1.0),
        ("a1", "right", 11.0, "A1", 1.0),
        ("a2", "left", 12.0, "A2", 0.0),
    ]
    table = placed(rows, {10: 4, 11: 3, 12: 1})
    assert table.loc["A", "psd95"] == pytest.approx(0.75)
    assert table.loc["A", "weighting"] == "voxels"
    assert table.loc["A", "covered_share"] == pytest.approx(1.0)
    assert table.loc["A", "psd95_right"] == pytest.approx(1.0)
    assert table.loc["A", "psd95_left"] == pytest.approx(0.75)


def test_units_not_drawn_in_the_annotation_weigh_alike():
    """A unit without voxels makes the units weigh alike; the covered share is unknown."""
    rows = [("a1", "left", 11.0, "A1", 1.0), ("a2", "left", 12.0, "A2", 0.0)]
    table = placed(rows, {10: 4, 11: 3})
    assert table.loc["A", "psd95"] == pytest.approx(0.5)
    assert table.loc["A", "weighting"] == "equal"
    assert np.isnan(table.loc["A", "covered_share"])


def test_a_partly_sampled_structure_reports_its_covered_share():
    """Only A1 sampled: A takes its value, and 3 of A's 4 voxels are covered."""
    table = placed([("a1", "left", 11.0, "A1", 0.2)], {10: 4, 11: 3, 12: 1})
    assert table.loc["A", "psd95"] == pytest.approx(0.2)
    assert table.loc["A", "covered_share"] == pytest.approx(0.75)


def test_a_structure_under_a_sampled_region_stays_missing_with_its_reason():
    """C is never filled from P; the table says only a region above it was sampled."""
    rows = [("p", "left", 30.0, "P", 0.9), ("b", "left", 20.0, "B", 0.3)]
    out = sy.match_samples(samples(rows), ONTOLOGY, STRUCTURES)
    per_structure = sy.structure_density(out, {20: 5, 30: 9})
    set_table = pd.DataFrame(
        dict(
            structure=["B", "C"],
            acronym=["B", "C"],
            division=["DIV", "DIV"],
            in_set=[True, True],
        )
    )
    table = sy.density_table(
        set_table, ["B", "C"], per_structure, out, ONTOLOGY, STRUCTURES
    ).set_index("structure")
    assert table.loc["B", "measured"] and table.loc["B", "psd95"] == pytest.approx(0.3)
    assert not table.loc["C", "measured"] and np.isnan(table.loc["C", "psd95"])
    assert "region above it was sampled (P)" in table.loc["C", "reason"]


def test_the_coverage_rule_asks_for_80_percent_of_the_fit():
    """101 of 126 is enough, 100 is not (0.8 x 126 = 100.8)."""
    assert sy.in_main_model(101, 126)
    assert not sy.in_main_model(100, 126)


@pytest.mark.skipif(not sy.DENSITY_FILE.exists(), reason="synaptome not fetched")
def test_the_fetched_source_is_the_pinned_one():
    """The files match the commit; 37 subtypes by 775 samples, one sample empty."""
    sy.check_files()
    density, info = sy.load_source()
    assert density.shape == (37, 775)
    table = sy.sample_table(density, info)
    assert int((~table["measured"]).sum()) == 1


@pytest.mark.skipif(not sy.DENSITY.exists(), reason="run_synaptome has not run")
def test_todays_density_table():
    """Measured where no reason is given; densities in 0..1; the fit inside the set."""
    table = sy.load_density()
    assert (table["measured"] == (table["reason"] == "")).all()
    measured = table[table["measured"]]
    assert measured["psd95"].between(0, 1).all()
    assert table.loc[~table["measured"], "psd95"].isna().all()
    assert (table.loc[table["in_fit"], "in_set"]).all()
    assert Path(sy.COVERAGE).exists()
