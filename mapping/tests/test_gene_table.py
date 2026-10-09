"""The gene table: the union of the panels, exclusions, profiles, GO ancestry, sets."""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from sepmap.config import DATA
from sepmap.ish import gene_sets, gene_table, panel_fetch, regions

# the tables of run_ish_gene_table on this data root
HAVE_TABLES = gene_table.GENE_TABLE.exists() and gene_table.DOCUMENTATION.exists()


def panels():
    """Two small panels sharing the gene Gria1 and its experiment 11."""
    p9 = pd.DataFrame(
        dict(
            symbol=["Gria1", "Th"],
            experiment_id=["11", "12"],
            plane=["coronal", "coronal"],
        )
    )
    ontology = pd.DataFrame(
        dict(
            symbol=["Gria1", "Gria1", "Dlg4"],
            experiment_id=["11", "21", "22"],
            plane=["coronal", "sagittal", "sagittal"],
        )
    )
    return p9, ontology


def test_union_has_one_row_per_experiment_with_both_labels():
    """An experiment in both panels is one row, labelled as in both."""
    p9, ontology = panels()
    repair = [dict(symbol="Th", experiment_id="31", plane="sagittal")]
    table = gene_table.union_rows(p9, ontology, repair).set_index("experiment_id")
    assert sorted(table.index) == ["11", "12", "21", "22", "31"]
    assert table.loc["11", "p9_experiment"] and table.loc["11", "ontology_experiment"]
    assert not table.loc["21", "p9_experiment"] and table.loc["21", "ontology_experiment"]
    assert table.loc["31", "repair_experiment"]
    assert not table.loc["11", "repair_experiment"]


def test_unusable_grid_keeps_its_row_with_a_reason(tmp_path, monkeypatch):
    """A grid outside the reference box, or missing, is excluded with the reason why."""
    monkeypatch.setattr(regions, "ISH_DIR", tmp_path)
    monkeypatch.setattr(panel_fetch, "DEST", tmp_path)
    monkeypatch.setattr(regions, "FETCH_FAILURES", tmp_path / "none.csv")
    for eid, dims in (("11", "67 41 58"), ("21", "68 40 50")):
        (tmp_path / f"{eid}_energy.mhd").write_text(f"DimSize = {dims}\n")
        (tmp_path / f"{eid}_energy.raw").write_bytes(b"")
    p9, ontology = panels()
    table = gene_table.exclusion_reasons(gene_table.union_rows(p9, ontology, []))
    table = table.set_index("experiment_id")
    assert not table.loc["11", "excluded"]
    assert table.loc["21", "excluded"]
    assert "reference box" in table.loc["21", "exclude_reason"]
    assert table.loc["22", "exclude_reason"] == "grid not on disk"
    assert (table["exclude_reason"] != "").eq(table["excluded"]).all()


def test_profile_of_one_experiment_orders_structures_as_the_experiment():
    """A gene measured once has a merged profile in the order of its one experiment."""
    rng = np.random.default_rng(0)
    values = rng.gamma(2.0, 3.0, 40)
    profile = {f"s{i}": float(v) for i, v in enumerate(values)}
    rows = pd.DataFrame(gene_table.profile_rows({"Gria1": {"11": profile}}))
    merged = rows.set_index("structure")["rank_mean"]
    rho = spearmanr(merged[list(profile)], list(profile.values())).statistic
    assert rho == pytest.approx(1.0)
    assert merged.min() == 0.0 and merged.max() == 1.0


def test_go_ancestry_follows_is_a_and_part_of(tmp_path, monkeypatch):
    """A term under postsynaptic density, by is_a then part_of, has postsynapse above."""
    monkeypatch.setattr(gene_table, "CACHE", tmp_path)
    (tmp_path / "go-basic.obo").write_text(
        "format-version: 1.2\n"
        "data-version: releases/2026-01-01\n\n"
        "[Term]\nid: GO:0098794\nname: postsynapse\n\n"
        "[Term]\nid: GO:0099572\nname: postsynaptic specialization\n"
        "relationship: part_of GO:0098794 ! postsynapse\n\n"
        "[Term]\nid: GO:0014069\nname: postsynaptic density\n"
        "alt_id: GO:9999999\n"
        "is_a: GO:0099572 ! postsynaptic specialization\n\n"
        "[Typedef]\nid: part_of\nname: part of\n",
        encoding="utf-8",
    )
    parents, names, release = gene_table.load_obo(offline=True)
    assert release == "releases/2026-01-01"
    assert names["GO:0014069"] == "postsynaptic density"
    assert "GO:0098794" in gene_table.with_ancestors({"GO:0014069"}, parents)
    assert "GO:0098794" in gene_table.with_ancestors({"GO:9999999"}, parents)
    assert gene_table.with_ancestors({"GO:0098794"}, parents) == {"GO:0098794"}


def test_not_qualifiers_are_left_out():
    """An annotation with a NOT qualifier does not count, a dict branch is read too."""
    record = {
        "go": {
            "CC": [
                {"id": "GO:0014069", "qualifier": "located_in"},
                {"id": "GO:0098793", "qualifier": "NOT|located_in"},
            ],
            "MF": {"id": "GO:0004971", "qualifier": "enables"},
        }
    }
    assert gene_table.go_annotations(record, "CC") == {"GO:0014069"}
    assert gene_table.go_annotations(record, "MF") == {"GO:0004971"}


def test_gene_sets_follow_their_rules():
    """Panel roles, GO compartments with both excluded, and the marker lists."""
    genes = ["Gria1", "Cacng8", "Dlg4", "Syp", "Bsn", "Gad1", "Aqp4", "Vip", "Fos"]
    role = {"Gria1": "subunit", "Cacng8": "localisation", "Dlg4": "control_psd"}
    post, pre = gene_sets.POSTSYNAPSE, gene_sets.PRESYNAPSE
    components = {
        "Cacng8": {post},
        "Dlg4": {post},
        "Syp": {pre},
        "Bsn": {pre, post},
        "Gad1": {pre},
    }
    sets = gene_sets.set_members(genes, role, components)
    assert sets["subunits"] == ["Gria1"]
    assert sets["localisation"] == ["Cacng8"]
    assert sets["other postsynaptic"] == ["Dlg4"]
    assert sets["presynaptic"] == ["Gad1", "Syp"]
    assert sets["GABAergic markers"] == ["Gad1", "Vip"]
    assert sets["glia"] == ["Aqp4"]
    assert list(sets) == list(gene_sets.SET_ORDER)


@pytest.mark.skipif(not HAVE_TABLES, reason="data not connected")
def test_written_gene_table():
    """451 genes, a row per experiment, all of P9's genes, every exclusion explained."""
    table = gene_table.load_gene_table()
    doc = pd.read_csv(
        gene_table.DOCUMENTATION, encoding="utf-8-sig", keep_default_na=False
    )
    p9 = gene_table.read_panel("targets")
    assert table["symbol"].nunique() == 451
    assert not table.duplicated(["symbol", "experiment_id"]).any()
    assert not table["experiment_id"].duplicated().any()
    assert set(p9["symbol"]) <= set(table["symbol"])
    assert table.loc[table["excluded"], "exclude_reason"].notna().all()
    with_category = set(doc.loc[doc["p9_category"] != "", "symbol"])
    assert with_category == set(p9["symbol"])


@pytest.mark.skipif(
    not (DATA / "adult_v2" / "ish" / "gene_reliability.csv").exists(),
    reason="data not connected",
)
def test_reliability_without_qc_equals_that_of_5_october():
    """With no section set missing, 15 genes' reliabilities equal gene_reliability.csv."""
    old = pd.read_csv(DATA / "adult_v2" / "ish" / "gene_reliability.csv")
    genes = sorted(old.dropna(subset=["reliability"])["symbol"])[:15]
    ontology = gene_table.read_panel("ontology")
    mine = ontology[ontology["symbol"].isin(genes)].copy()
    mine["p9_experiment"] = False
    region = gene_table.region_table(mine, {})
    rel = gene_table.gene_reliability(gene_table.experiment_profiles(region))
    got = rel.set_index("symbol")["reliability"]
    want = old.set_index("symbol")["reliability"]

    # gene_reliability.csv holds four decimals
    assert np.allclose(got[genes], want[genes], atol=6e-5)
