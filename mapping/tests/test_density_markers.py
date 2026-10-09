"""The synapse-density genes: the pool, the exclusions, the choice and its halves."""

import gzip

import numpy as np
import pandas as pd
import pytest

from sepmap.adult import density_markers as dm
from sepmap.adult import synaptome
from sepmap.ish import gene_sets, gene_table
from sepmap.structures import STRUCTURE_SET, load_structure_set

# a small ontology: the postsynaptic density below the specialization by is_a, its
# membrane below it by part_of, the synapse above them; a child of an exclusion term
# by is_a; a plasticity term
OBO_TEXT = (
    "format-version: 1.2\n"
    "data-version: releases/2026-01-01\n\n"
    "[Term]\nid: GO:0045202\nname: synapse\n\n"
    "[Term]\nid: GO:0099572\nname: postsynaptic specialization\n"
    "relationship: part_of GO:0045202 ! synapse\n\n"
    "[Term]\nid: GO:0014069\nname: postsynaptic density\nalt_id: GO:9999999\n"
    "is_a: GO:0099572 ! postsynaptic specialization\n\n"
    "[Term]\nid: GO:0098839\nname: postsynaptic density membrane\n"
    "relationship: part_of GO:0014069 ! postsynaptic density\n\n"
    "[Term]\nid: GO:2000311\nname: regulation of AMPA receptor activity\n\n"
    "[Term]\nid: GO:2000969\nname: positive regulation of AMPA receptor activity\n"
    "is_a: GO:2000311 ! regulation of AMPA receptor activity\n\n"
    "[Term]\nid: GO:0048167\nname: regulation of synaptic plasticity\n\n"
    "[Typedef]\nid: part_of\nname: part of\n"
)

# twenty structures, every one measured for Gria1
STRUCTURES = [f"S{k:02d}" for k in range(20)]

# the genes of the synthetic gene table: (symbol, role, usable experiments, structures
# measured, GO terms); Edge is measured in exactly 90% of the structures
GENES = [
    ("Gria1", "subunit", 2, 20, ["GO:0014069"]),
    ("Homer1", "control_psd", 2, 20, ["GO:0098839"]),
    ("Dlg4", "control_psd", 2, 20, ["GO:0014069"]),
    ("Camk2a", "control_psd", 2, 20, ["GO:0014069", "GO:2000969"]),
    ("Loc1", "localisation", 2, 20, ["GO:0014069"]),
    ("Syn9", "", 2, 20, ["GO:0045202"]),
    ("One1", "control_psd", 1, 20, ["GO:0014069"]),
    ("Thin", "control_psd", 2, 17, ["GO:0014069"]),
    ("Edge", "control_psd", 2, 18, ["GO:0014069"]),
    ("Plast", "control_psd", 2, 20, ["GO:0014069", "GO:0048167"]),
    ("None1", "", 2, 20, []),
]


def synthetic_rule(tmp_path):
    """The candidate and exclusion tables of the synthetic gene table."""
    obo = tmp_path / "go-basic.obo"
    obo.write_text(OBO_TEXT, encoding="utf-8")
    parents, names, _ = gene_table.read_obo(obo)
    genes = pd.DataFrame(
        [dict(symbol=g, ontology_role=r, n_experiments_used=n) for g, r, n, _, _ in GENES]
    )
    profiles = {
        g: {s: float(k) for k, s in enumerate(STRUCTURES[:m])} for g, _, _, m, _ in GENES
    }
    terms = {g: {t: {"IDA"} for t in go} for g, _, _, _, go in GENES if go}
    return dm.candidate_table(genes, profiles, STRUCTURES, terms, parents, names)


def test_the_pool_takes_postsynaptic_genes_and_leaves_out_the_ampa_linked(tmp_path):
    """A gene below the density by part_of enters; subunit, family, set, term do not."""
    candidates, _ = synthetic_rule(tmp_path)
    table = candidates.set_index("symbol")
    assert set(table.index[table["in_pool"]]) == {"Homer1", "Edge", "Plast"}
    assert table.loc["Gria1", "exclude_rules"] == dm.RULE_SUBUNIT
    assert dm.RULE_FAMILY in table.loc["Dlg4", "exclude_rules"]
    assert table.loc["Camk2a", "exclude_rules"] == dm.RULE_GO
    assert table.loc["Loc1", "exclude_rules"] == dm.RULE_LOCALISATION

    # every gene excluded was otherwise eligible, so it is drawn in the ranking
    assert table.loc[["Gria1", "Dlg4", "Camk2a", "Loc1"], "eligible"].all()


def test_each_gene_out_of_the_pool_says_which_criterion_it_fails(tmp_path):
    """The first criterion failed, in the rule's order, is the reason given."""
    candidates, _ = synthetic_rule(tmp_path)
    reason = candidates.set_index("symbol")["reason"]
    assert reason["One1"].startswith("1 usable Allen experiments")
    assert reason["Thin"].startswith("measured in 17 of the 20")
    assert reason["None1"] == "no annotation in the mouse GAF"
    assert reason["Syn9"].startswith("not annotated to the postsynaptic")
    assert reason["Camk2a"].startswith("excluded: ")
    assert reason[["Homer1", "Edge", "Plast"]].eq("").all()


def test_an_exclusion_is_found_through_a_term_below_it(tmp_path):
    """Camk2a's positive regulation term lies below the excluding term by is_a."""
    _, excluded = synthetic_rule(tmp_path)
    row = excluded[(excluded["symbol"] == "Camk2a") & (excluded["rule"] == dm.RULE_GO)]
    assert len(row) == 1
    assert row["term"].iloc[0] == "GO:2000311 regulation of AMPA receptor activity"
    assert row["via"].iloc[0].startswith("GO:2000969")
    assert row["evidence"].iloc[0] == "IDA"


def test_plasticity_terms_are_listed_and_never_exclude(tmp_path):
    """A pool gene with a plasticity term keeps its place, the term recorded."""
    candidates, _ = synthetic_rule(tmp_path)
    plast = candidates.set_index("symbol").loc["Plast"]
    assert plast["in_pool"]
    assert "GO:0048167 regulation of synaptic plasticity" in plast["plasticity_terms"]


def test_descendants_follow_is_a_part_of_and_alt_ids(tmp_path):
    """Below the specialization: the density (is_a), its membrane (part_of), alt ids."""
    obo = tmp_path / "go-basic.obo"
    obo.write_text(OBO_TEXT, encoding="utf-8")
    parents, _, _ = gene_table.read_obo(obo)
    below = dm.descendants("GO:0099572", dm.children_of(parents))
    assert {"GO:0014069", "GO:0098839", "GO:9999999"} <= below
    assert "GO:0045202" not in below


def test_a_term_named_otherwise_stops_the_run():
    """An id whose release name differs is refused, and the term of the name is given."""
    names = {"GO:0000001": "one", "GO:0000002": "two"}
    dm.check_terms({"GO:0000001": "one"}, names)
    with pytest.raises(ValueError, match="GO:0000002"):
        dm.check_terms({"GO:0000001": "two"}, names)


def test_the_gaf_is_read_for_the_genes_asked_without_not(tmp_path):
    """Header skipped, other genes skipped, a NOT annotation counted and left out."""
    lines = [
        "!gaf-version: 2.2",
        "!date-generated: 2026-01-01",
        "MGI\tMGI:1\tHomer1\tlocated_in\tGO:0014069\tPMID:1\tIDA\t\tC\tname\t\tprotein"
        "\ttaxon:10090\t20260101\tMGI\t\t",
        "MGI\tMGI:1\tHomer1\tNOT|located_in\tGO:0045202\tPMID:1\tIDA\t\tC\tname\t\t"
        "protein\ttaxon:10090\t20260101\tMGI\t\t",
        "MGI\tMGI:2\tOther\tlocated_in\tGO:0014069\tPMID:1\tIDA\t\tC\tname\t\tprotein"
        "\ttaxon:10090\t20260101\tMGI\t\t",
    ]
    path = tmp_path / "mouse.gaf.gz"
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    annotations, n_not = dm.read_gaf(path, {"Homer1"})
    assert n_not == 1
    assert list(annotations["go_id"]) == ["GO:0014069"]
    assert dm.file_version("x.gaf.gz", path.read_bytes()) == "date-generated: 2026-01-01"


def planted_profiles(n, noise, seed=0):
    """PSD95 over n structures, and genes that are it plus noise of the given SDs."""
    rng = np.random.default_rng(seed)
    structures = [f"S{k:02d}" for k in range(n)]
    psd95 = pd.Series(np.arange(n, dtype=float), index=structures)
    profiles = {}
    for gene, sd in noise.items():
        values = psd95.to_numpy() + rng.normal(0, sd, n)
        profiles[gene] = dict(zip(structures, values))
    return structures, psd95, profiles


def eligible_table(pool, excluded=()):
    """A candidate table of eligible genes, some of them excluded."""
    rows = [
        dict(symbol=g, eligible=True, in_pool=True, excluded=False, exclude_rules="")
        for g in pool
    ]
    rows += [
        dict(
            symbol=g, eligible=True, in_pool=False, excluded=True, exclude_rules="GO term"
        )
        for g in excluded
    ]
    return pd.DataFrame(rows)


def test_each_excluded_gene_counts_once_under_its_first_reason():
    """Subunit before family before GO; a GO-only gene under its first listed term."""
    first, later = list(dm.EXCLUSION_TERMS)[:2]
    rows = [
        ("Gria2", True, dm.RULE_SUBUNIT, "Gria1 to Gria4"),
        ("Gria2", True, dm.RULE_FAMILY, "Schwenk et al. 2012"),
        ("Gria2", True, dm.RULE_GO, f"{first} {dm.EXCLUSION_TERMS[first]}"),
        ("Fam1", True, dm.RULE_FAMILY, "Schwenk et al. 2012"),
        ("Fam1", True, dm.RULE_GO, f"{later} {dm.EXCLUSION_TERMS[later]}"),
        ("Go1", True, dm.RULE_GO, f"{later} {dm.EXCLUSION_TERMS[later]}"),
        ("Go1", True, dm.RULE_GO, f"{first} {dm.EXCLUSION_TERMS[first]}"),
        ("Go2", True, dm.RULE_GO, f"{later} {dm.EXCLUSION_TERMS[later]}"),
        ("Rare", False, dm.RULE_GO, f"{first} {dm.EXCLUSION_TERMS[first]}"),
    ]
    excluded = pd.DataFrame(rows, columns=["symbol", "eligible", "rule", "term"])
    reasons = dm.exclusion_reasons(excluded).set_index("reason")
    assert reasons["n_genes"].sum() == 4
    assert reasons.loc[dm.RULE_SUBUNIT, "genes"] == "Gria2"
    assert reasons.loc[dm.RULE_FAMILY, "genes"] == "Fam1"
    assert reasons.loc[f"{first} {dm.EXCLUSION_TERMS[first]}", "genes"] == "Go1"
    assert reasons.loc[f"{later} {dm.EXCLUSION_TERMS[later]}", "genes"] == "Go2"


def test_the_choice_takes_the_pool_genes_highest_with_psd95():
    """An excluded gene above them all is passed over; the next three are chosen."""
    noise = dict(A=0.1, B=3.0, C=6.0, D=9.0, E=40.0)
    structures, psd95, profiles = planted_profiles(40, noise)
    agreement = dm.agreement_table(
        eligible_table(["B", "C", "D", "E"], ["A"]), profiles, psd95, structures
    )
    assert dm.chosen_genes(agreement) == ["B", "C", "D"]
    row = agreement.set_index("symbol").loc["A"]
    assert not row["chosen"] and np.isnan(row["rank"])
    assert agreement["rho"].is_monotonic_decreasing


def test_with_fewer_pool_genes_than_asked_the_term_takes_them_all():
    """Two genes in the pool give a term of two."""
    structures, psd95, profiles = planted_profiles(30, dict(B=3.0, C=6.0))
    agreement = dm.agreement_table(
        eligible_table(["B", "C"]), profiles, psd95, structures
    )
    assert dm.chosen_genes(agreement) == ["B", "C"]


def test_halves_choose_the_planted_genes_and_hold_out_their_agreement():
    """Three genes planted as PSD95 are chosen on nearly every half, agreeing held out."""
    noise = dict(B=3.0, C=3.0, D=3.0, E=60.0, F=60.0, G=60.0)
    structures, psd95, profiles = planted_profiles(40, noise)
    pool = sorted(noise)
    halves = dm.half_split(pool, profiles, psd95, structures, {})
    selection = dm.selection_table(halves, pool).set_index("symbol")["share"]
    assert (selection[["B", "C", "D"]] > 0.95).all()
    assert (selection[["E", "F", "G"]] < 0.05).all()

    # each half holds 20 structures; three copies of the map plus noise of SD 3 (a
    # tenth of its range) agree with it above 0.9 on any of them
    assert (halves["n_chosen"] == 20).all()
    assert halves["rho_chosen"].min() > 0.9


def test_the_comparison_reads_a_gene_list_and_a_map_alike():
    """A gene list is ranked on the structures; a map equal to PSD95 agrees at 1."""
    noise = dict(B=3.0, C=3.0, D=3.0, E=60.0)
    structures, psd95, profiles = planted_profiles(30, noise)
    fixed = {"list": ("E",), "map": psd95.copy()}
    members = {"list": "E", "map": "PSD95 itself"}
    division = {s: "X" if k % 2 else "Y" for k, s in enumerate(structures)}
    halves = dm.half_split(sorted(noise), profiles, psd95, structures, fixed, division)
    table = dm.comparison_table(
        ["B", "C", "D"], profiles, psd95, structures, halves, fixed, members, division
    ).set_index("composite")
    assert list(table.index) == [dm.CHOSEN, "list", "map"]
    assert table.loc["map", "rho_full_set"] == pytest.approx(1.0)
    assert table.loc["map", "rho_held_out_lo"] == pytest.approx(1.0)
    assert table.loc[dm.CHOSEN, "genes"] == "B C D"
    assert table.loc["list", "rho_full_set"] < table.loc[dm.CHOSEN, "rho_full_set"]


def test_inside_divisions_the_contrast_between_them_does_not_count():
    """A gene that follows only the divisions' levels agrees with PSD95 between them."""
    rng = np.random.default_rng(0)
    structures = [f"S{k:02d}" for k in range(40)]
    division = {s: "X" if k < 20 else "Y" for k, s in enumerate(structures)}
    level = np.array([0.0 if division[s] == "X" else 10.0 for s in structures])

    # PSD95: the divisions' levels and a pattern inside each; the gene: the levels and
    # noise unrelated to that pattern
    psd95 = pd.Series(level + rng.normal(0, 1, 40), index=structures)
    gene = pd.Series(level + rng.normal(0, 1, 40), index=structures)

    # two halves set apart, their order inside random, rank together at about 0.75
    assert dm.rho_with(gene, psd95)[1] > 0.6
    assert abs(dm.within_rho(gene, psd95, division)) < 0.35
    assert dm.within_rho(psd95, psd95, division) == pytest.approx(1.0)

    # with fewer structures in every division than a division needs, there is no value
    few = {s: f"D{k}" for k, s in enumerate(structures)}
    assert np.isnan(dm.within_rho(gene, psd95, few))


def test_a_rule_that_gives_other_genes_than_the_settings_stops_the_run():
    """The settings' genes pass in any order; another gene does not."""
    named = list(dm.DENSITY_MARKERS["chosen"])
    dm.check_chosen(list(reversed(named)))
    with pytest.raises(ValueError, match="changed"):
        dm.check_chosen(named[:-1] + ["Gria1"])


def test_the_model_structures_need_gria1_and_every_gene_chosen():
    """A structure without Gria1, or without one gene chosen, is out with its reason."""
    set_table = pd.DataFrame(
        dict(
            structure=["A", "B", "C", "D"],
            acronym=["a", "b", "c", "d"],
            division=["TH"] * 4,
            in_set=[True, True, True, False],
        )
    )
    profiles = {
        "Gria1": {"A": 1.0, "B": 1.0, "D": 1.0},
        "X": {"A": 1.0, "B": 1.0, "C": 1.0, "D": 1.0},
        "Y": {"A": 1.0, "C": 1.0, "D": 1.0},
    }
    psd95 = pd.Series({"A": 1.0, "B": np.nan})
    table = dm.structure_table(set_table, profiles, ["X", "Y"], psd95).set_index(
        "structure"
    )
    assert list(table.index) == ["A", "B", "C"]
    assert list(table["in_model"]) == [True, False, False]
    assert table.loc["B", "reason"] == "not measured: Y"
    assert table.loc["C", "reason"] == "Gria1 not measured"
    assert list(table["psd95_measured"]) == [True, False, False]


# the files the rule may read: the gene table, the genes' profiles, the declared set and
# the synaptome's density table; nothing that holds a nano value
RULE_INPUTS = (
    gene_table.GENE_TABLE.name,
    gene_table.PROFILES.name,
    STRUCTURE_SET.name,
    synaptome.DENSITY.name,
)
CONNECTED = (
    dm.GAF.exists()
    and gene_table.GENE_TABLE.exists()
    and gene_table.PROFILES.exists()
    and synaptome.DENSITY.exists()
)


@pytest.mark.skipif(not dm.GAF.exists(), reason="data not connected")
def test_the_fetched_go_files_are_the_pinned_release_and_name_the_rules_terms():
    """Both files carry the release's hash, and every term of the rule its name."""
    dm.check_files()
    _, names, release = gene_table.read_obo(dm.OBO)
    dm.check_terms(dm.POOL_TERMS, names)
    dm.check_terms(dm.EXCLUSION_TERMS, names)
    assert release.startswith("releases/")


@pytest.mark.skipif(not CONNECTED, reason="data not connected")
def test_the_rule_reads_no_nano_value_and_gives_the_settings_genes(monkeypatch):
    """Only the rule's inputs are read; Dlg4 and Camk2a are excluded, Homer1 is not."""
    read_csv = pd.read_csv
    read = []

    def guarded(path, *args, **kwargs):
        """pandas.read_csv for the rule's inputs only."""
        read.append(str(path))
        if not str(path).endswith(RULE_INPUTS):
            raise AssertionError(f"the rule read {path}, which it must not")
        return read_csv(path, *args, **kwargs)

    monkeypatch.setattr(pd, "read_csv", guarded)
    monkeypatch.setattr(np, "load", lambda *a, **k: pytest.fail("np.load was called"))
    genes = gene_table.per_gene(gene_table.load_gene_table())
    profiles = gene_table.load_profiles()
    set_table = load_structure_set()
    declared = sorted(set_table.loc[set_table["in_set"], "structure"])
    psd95 = synaptome.load_density()[synaptome.MEASURE]
    measured = [s for s in declared if np.isfinite(psd95.get(s, np.nan))]
    annotations, _ = dm.read_gaf(dm.GAF, set(genes["symbol"]))
    parents, names, _ = gene_table.read_obo(dm.OBO)
    terms = dm.annotated_terms(annotations)
    candidates, _ = dm.candidate_table(genes, profiles, declared, terms, parents, names)
    agreement = dm.agreement_table(candidates, profiles, psd95, measured)
    dm.check_chosen(dm.chosen_genes(agreement))
    assert len(read) == len(RULE_INPUTS)

    table = candidates.set_index("symbol")
    assert table.loc[["Dlg4", "Camk2a"], "eligible"].all()
    assert table.loc[["Dlg4", "Camk2a"], "excluded"].all()
    assert table.loc["Homer1", "in_pool"]


@pytest.mark.skipif(not dm.GAF.exists(), reason="data not connected")
def test_the_pinned_release_gives_the_family_its_go_genes():
    """GO:0032281 in release 2026-08-05's mouse GAF lists the family's GO genes."""
    genes = set()
    with gzip.open(dm.GAF, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            field = line.rstrip("\n").split("\t")
            if field[4] == "GO:0032281" and "NOT" not in field[3].split("|"):
                genes.add(field[2])
    assert genes == set(gene_sets.GO_AMPA_COMPLEX)
