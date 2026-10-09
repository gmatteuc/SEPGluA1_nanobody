"""Known-answer checks of the genes that follow the map and the tests named for the
leftover (ish.top_genes)."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control

from sepmap.adult import beyond_density as bd
from sepmap.ish import gene_sets, spatial_null, top_genes

TABLES = Path(top_genes.TOP_TABLE).parent


def test_the_genes_and_their_tiers():
    """Past BH, the named and the family's members the leftover holds; best first."""
    nano = pd.DataFrame(
        dict(
            symbol=["Cacng8", "Gria1", "Top", "Shisa6", "Cacng2", "Other"],
            q_all=[0.001, 0.5, 0.01, 0.6, 0.9, 0.5],
            rho=[0.8, 0.3, 0.7, 0.1, 0.05, 0.2],
        )
    )
    leftover = pd.DataFrame(dict(symbol=["Cacng8", "Gria1", "Top", "Shisa6", "Other"]))
    chosen = top_genes.chosen_genes(nano, leftover, q=0.05)
    assert list(chosen["symbol"]) == ["Cacng8", "Top", "Gria1", "Shisa6"]
    tier = dict(zip(chosen["symbol"], chosen["tier"]))
    assert tier["Cacng8"].startswith("1")
    assert tier["Shisa6"].startswith("2")
    assert tier["Top"].startswith("3") and tier["Gria1"].startswith("3")
    row = chosen.set_index("symbol").loc["Gria1"]
    assert row["named"] and not row["top"] and not row["family"]

    # Cacng2 is in the family, but the leftover table does not hold it
    assert "Cacng2" not in set(chosen["symbol"])


def test_every_family_member_is_listed_with_its_reason():
    """Tested, no usable experiment, too few structures, or not in the gene table."""
    genes = pd.DataFrame(
        dict(symbol=["Shisa6", "Cacng2", "Cacng3"], n_experiments_used=[1, 0, 2])
    )
    leftover = pd.DataFrame(dict(symbol=["Shisa6"]))
    members = top_genes.family_members(genes, leftover).set_index("symbol")
    assert len(members) == len(gene_sets.AMPA_FAMILY)
    assert "Gria1" not in members.index
    assert members.loc["Shisa6", "tested"]
    assert members.loc["Cacng2", "reason"] == "no usable Allen experiment"
    assert members.loc["Cacng3", "reason"].startswith("fewer than")
    assert members.loc["Prrt1", "reason"].startswith("not in the gene table")
    assert "partner subunit" in members.loc["Gria2", "sources"]
    assert "Schwenk 2012" in members.loc["Cacng8", "sources"]
    assert members.loc["Grid1", "sources"].startswith("GO:0032281")


def test_ranks_come_only_from_variants_holding_every_gene():
    """A variant of P9's genes alone ranks among fewer, so its rank is left out."""
    robust = pd.DataFrame(
        dict(
            variant=["primary"] * 3 + ["full"] * 3 + ["p9"] * 2,
            symbol=["A", "B", "C"] * 2 + ["A", "B"],
            rho=[0.5, 0.4, 0.1, 0.45, 0.35, 0.2, 0.3, 0.6],
            rank_all=[1, 2, 3, 1, 2, 3, 2, 1],
        )
    )
    ranges = top_genes.robustness_ranges(robust, ["B"]).iloc[0]
    assert ranges["n_variants"] == 2
    assert ranges["rho_variants_min"] == 0.35 and ranges["rho_variants_max"] == 0.6
    assert ranges["rank_variants_min"] == 2 and ranges["rank_variants_max"] == 2


def test_bh_within_the_family_is_over_the_family_alone():
    """q within the family equals BH on its p alone, and is smaller than over more."""
    p = pd.Series([0.001, 0.02, 0.2, 0.6], index=["a", "b", "c", "d"])
    q = top_genes.within_family_q(p)
    assert np.allclose(q.to_numpy(), false_discovery_control(p.to_numpy(), method="bh"))
    wider = false_discovery_control(np.concatenate([p.to_numpy(), np.full(20, 0.5)]))
    assert (q.to_numpy() <= wider[:4] + 1e-12).all()


def test_synaptic_go_terms_follow_ancestry_and_drop_not():
    """A term below GO:0045202 counts, a membrane term does not, a NOT never does."""
    documentation = pd.DataFrame(
        dict(
            symbol=["G"],
            name=["a gene"],
            ontology_role=["localisation"],
            ontology_go_terms=["GO:0000002"],
            gene_sets=["localisation"],
            p9_category=[""],
        )
    )
    cc = [
        dict(id="GO:0000001", qualifier="located_in"),
        dict(id="GO:0000003", qualifier="located_in"),
        dict(id="GO:0000004", qualifier="NOT located_in"),
    ]
    records = {"G": {"go": {"CC": cc}}}
    parents = {
        "GO:0000001": {"GO:0045202"},
        "GO:0000004": {"GO:0045202"},
        "GO:0000003": {"GO:0016020"},
    }
    names = {"GO:0000001": "postsynaptic density", "GO:0000002": "a panel term"}
    row = top_genes.annotation_table(["G"], documentation, records, parents, names)
    row = row.iloc[0]
    assert row["go_synapse_terms"] == "GO:0000001 postsynaptic density"
    assert row["go_panel_terms"] == "GO:0000002 a panel term"


def test_a_remainder_surrogate_keeps_the_genes_fit_on_the_model():
    """Each map's fit on the design is the gene's; its remainder is as large."""
    rng = np.random.default_rng(0)
    n = 60
    xyz = rng.uniform(0, 10, size=(n, 3))
    d = np.sqrt(((xyz[:, None] - xyz[None]) ** 2).sum(-1))
    design = np.column_stack([rng.standard_normal((n, 3)), np.ones(n)])
    gene = design[:, 0] + rng.standard_normal(n)
    maps = top_genes.remainder_surrogates(gene, design, d, 20, seed=0)
    hat = design @ np.linalg.pinv(design)
    assert np.allclose(maps @ hat.T, np.tile(hat @ gene, (20, 1)))
    rest = gene - hat @ gene
    assert np.allclose((maps - hat @ gene).std(axis=1), rest.std())


def synthetic_inputs(n: int = 80, seed: int = 0):
    """Ten made-up adults whose map is Gria1, synapse density and a planted pattern.

    Returns the inputs of analysis 4, the structures' centroids and the planted
    pattern; every predictor of the main model is a smooth field on random points.
    """
    rng = np.random.default_rng(seed)
    structures = [f"s{i:02d}" for i in range(n)]
    xyz = rng.uniform(0, 10, size=(n, 3))
    centroids = pd.DataFrame(xyz, index=structures, columns=["ap_mm", "dv_mm", "ml_mm"])
    d = spatial_null.distance_matrix(centroids)
    fields = spatial_null.random_fields(d, (0.05, 1.0, 3.0), 8, rng)
    gria1, density, auto, planted = fields[0], fields[1], fields[2], fields[3]
    profiles = {"Gria1": gria1}
    for k, gene in enumerate(("Gria2", "Gria3", "Gria4")):
        profiles[gene] = fields[4 + k]
    for gene in bd.MARKERS:
        profiles[gene] = density + 0.3 * rng.standard_normal(n)
    role = {}
    for k in range(3):
        profiles[f"Psd{k}"] = density + 0.3 * rng.standard_normal(n)
        role[f"Psd{k}"] = bd.PSD_ROLE
    profiles["Planted"] = planted + 0.2 * rng.standard_normal(n)
    profiles["Noise"] = fields[7]
    expr = {g: dict(zip(structures, v)) for g, v in profiles.items()}
    truth = gria1 + density + 1.2 * planted
    nano = truth + 0.3 * rng.standard_normal((len(bd.ADULTS), n))
    inputs = bd.Inputs(
        structures=structures,
        nano=nano,
        auto=auto + 0.3 * rng.standard_normal((len(bd.ADULTS), n)),
        synapses=pd.DataFrame(index=structures),
        expr=expr,
        role=role,
        p9_genes=set(),
        division={s: "Isocortex" for s in structures},
        acronym={s: s for s in structures},
        rows=pd.DataFrame(),
        terms=bd.model_terms(),
    )
    return inputs, centroids


def test_a_gene_that_is_the_leftover_takes_it_and_noise_does_not():
    """A planted pattern takes a share beyond both nulls; a smooth noise map less."""
    inputs, centroids = synthetic_inputs()
    planted, nulls = top_genes.added_share(
        inputs, "Planted", centroids, n_surrogates=40, seed=0
    )
    noise, _ = top_genes.added_share(inputs, "Noise", centroids, n_surrogates=40)
    assert set(nulls) == {"plain", "alike"}
    assert len(nulls["plain"]) == 40
    # the planted pattern carries about a third of the map's variance
    assert planted["taken"] > 0.1
    assert planted["p_taken"] <= 0.05 and planted["p_taken_alike"] <= 0.05
    assert noise["taken"] < planted["taken"] / 3
    assert planted["left_main"] == pytest.approx(noise["left_main"])
    assert planted["left_with_gene"] == pytest.approx(
        planted["left_main"] - planted["taken"]
    )


def test_a_gene_missing_structures_is_scored_on_those_it_has():
    """The main model is refitted on the structures the gene has, its ceiling too."""
    inputs, centroids = synthetic_inputs()
    for s in inputs.structures[:5]:
        del inputs.expr["Planted"][s]
    row, _ = top_genes.added_share(inputs, "Planted", centroids, n_surrogates=10)
    assert row["n_added"] == len(inputs.structures) - 5


def group_maps(family_mean: float, seed: int = 0):
    """rho of 20 family genes and 20 controls, and a null of 500 surrogates."""
    rng = np.random.default_rng(seed)
    family = [f"f{i}" for i in range(20)]
    controls = [f"c{i}" for i in range(20)]
    rho = pd.Series(
        np.concatenate([rng.normal(family_mean, 0.05, 20), rng.normal(0.0, 0.05, 20)]),
        index=family + controls,
    )
    null = rng.normal(0.0, 0.1, size=(40, 500))
    pairs = dict(zip(family, controls))
    return family, pairs, {"leftover": (rho, null, family + controls)}


def test_a_family_above_its_null_and_its_controls_passes_both():
    """A family at +0.3 against controls and a null at 0 passes both group tests."""
    family, pairs, maps = group_maps(0.3)
    tests, nulls = top_genes.group_tests(family, pairs, maps)
    assert (tests["p"] < 0.01).all()
    assert set(nulls["leftover"]) == {"spatial", "labels"}
    assert (tests["role"] == "the test").all()


def test_a_family_like_its_controls_does_not_pass_against_them():
    """A family drawn as its controls are gives a label p well above 0.05."""
    family, pairs, maps = group_maps(0.0, seed=3)
    tests, _ = top_genes.group_tests(family, pairs, maps)
    matched = tests[tests["test"].str.contains("controls")].iloc[0]
    assert matched["p"] > 0.05


@pytest.mark.skipif(
    not (TABLES / "top_genes.csv").is_file(), reason="top_genes.csv not written"
)
def test_todays_top_genes():
    """Today's table: Cacng8 and Gria1, every family member tested, the named tests."""
    table = top_genes.load_top_genes().set_index("symbol")
    assert {"Cacng8", "Gria1"} <= set(table.index)
    members = pd.read_csv(top_genes.FAMILY_TABLE)
    tested = set(members.loc[members["tested"].astype(str) == "True", "symbol"])
    assert tested == set(table.index[table["family"]])
    assert table.loc["Cacng8", "tier"].startswith("1")
    leftover = pd.read_csv(bd.LEFTOVER_GENES).set_index("symbol")
    assert table.loc["Cacng8", "leftover_rho"] == pytest.approx(
        leftover.loc["Cacng8", "rho"]
    )
    tests = pd.read_csv(top_genes.NAMED_TESTS, dtype={"tier": str})
    assert list(tests["tier"]) == ["1", "2", "2", "2", "2"]
    assert set(tests["map"]) == {"leftover", "nano"}
