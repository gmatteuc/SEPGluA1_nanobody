"""Analysis 3 of the ISH line: which kinds of genes match the map.

Python route, in run order (tools\\venv_atlas; run_closeup in tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid
     3. run_cohort             cohort mean, SD and n per voxel
     4. run_compare            young against adult: maps, table
        run_replot             redraws those maps (not in a full run)
     5. run_region_plot        region statistics, no warping
     6. run_region_groups      the same by system and by layer
     7. run_video              cohort videos
     8. run_video_compare      young beside adult, plane by plane
     9. run_closeup            close-ups and flatmaps (venv_flat)
    10. run_diagnostics        sheets that audit each step
    11. run_panel_build        the 390-gene ontology panel (network, cached)
    12. run_panel_fetch        its ISH grids (network, once)
    13. run_structure_set      A1: the declared structures, the
                               adult profiles; figure 01
    14. run_ish_section_qc     A2: the experiments of both panels
                               and the repair (network, once);
                               section QC; QC sheets
    15. run_ish_gene_table     A9: region means, the gene table,
                               merged profiles, gene sets,
                               documentation; figure 02
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the
                               map, the null, autofluorescence (A8),
                               the Cacng8 - Gria1 gap; figures 03 to 06
    18. run_ish_robustness     A3: the ranking under other choices;
                               figure 07
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figure 08, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation   <- this script
                               against matched controls; figures 09, 10
    21. run_beyond_density     analysis 4: what receptor mRNA and
                               synaptic density leave; the leftover
    22. run_beyond_controls    seven attempts to break it
    23. run_beyond_calibration the same model on maps whose answer
                               is known
    24. run_beyond_regression  the regression, per structure
    25. run_beyond_figures     figures 11 and 12
    26. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 13
    27. run_ish_overview       figures 00 and 14; the numbers for the text

Tests the gene sets fixed in advance (sepmap/ish/gene_sets.py) against the adult
map: each set's median rho against the medians its genes give with the map's
surrogates, the two contrasts named in advance, and the localisation genes against
expression-matched postsynaptic controls once the subunit composite is removed,
with its positive control and its spatial p (the method is in
sepmap/ish/gene_sets.py; every gene's rho and null come from run_ish_gene_ranking).
Writes, in adult_v2/ish_analysis/ under the data root:

    tables/gene_sets.csv               per set and gene, the gene's rho with both
                                       maps (the context group included)
    tables/set_tests.csv               per map and set: median rho, null band, p, q
    tables/contrasts.csv               the contrasts named in advance
    tables/localisation_test.csv       per control pool, localisation and control
                                       gene: plain and partial rho, reliability,
                                       expression, its match
    tables/localisation_summary.csv    every test of the localisation design, and
                                       the test and positive control with the
                                       control pool of 5 October
    tables/numbers_gene_sets.csv       the numbers of this step, for the text
    figures/09_gene_sets.png           the sets against the null
    figures/10_localisation.png        localisation against matched controls

    python run_ish_gene_sets.py
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sepmap import config, plotting, structures
from sepmap.adult import profiles
from sepmap.ish import gene_ranking, gene_sets, gene_table, spatial_null
from sepmap.ish import plotting as ish_plotting

ISH_ANALYSIS = config.SETTINGS["ish_analysis"]

OUT = config.DATA / "adult_v2" / "ish_analysis"


def numbers_table(members, tests, contrasts, summary):
    """The numbers of this step that the text quotes, one row each."""
    rows = []
    for name, symbols in members.items():
        rows.append((f"set_size_{name}", len(symbols), f"genes of {name} with a rho"))
    for _, r in tests.iterrows():
        key = f"{r['map']}_{r['gene_set']}"
        rows += [
            (f"set_median_{key}", round(r["median_rho"], 3), "the set's median rho"),
            (f"set_p_{key}", round(r["p_spatial"], 5), "its spatial p"),
            (f"set_q_{key}", round(r["q"], 5), "its q over the tested sets"),
        ]
    for _, r in contrasts.iterrows():
        key = r["contrast"].replace(" ", "_")
        rows += [
            (f"{key}_difference", round(r["difference"], 3), "difference of medians"),
            (f"{key}_p_spatial", round(r["p_spatial"], 5), "its spatial p"),
            (f"{key}_p_labels", round(r["p_labels"], 5), "labels permuted"),
            (f"{key}_n", f"{r['n_first']} against {r['n_second']}", "genes per side"),
        ]
    for _, r in summary.iterrows():
        key = r["test"].replace(",", "").replace("'", "").replace(" ", "_")
        rows += [
            (f"localisation_{key}_difference", round(r["difference"], 4), r["test"]),
            (f"localisation_{key}_p", round(r["p_labels"], 5), r["test"]),
            (f"localisation_{key}_detectable", round(r["detectable"], 4), r["test"]),
        ]
    main = summary.set_index("test").loc["matched controls"]
    rows.append(
        ("localisation_p_spatial", round(main["p_spatial"], 5), "matched, surrogates")
    )
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


def gene_levels(table):
    """{gene: reliability} and {gene: median energy}, from the usable experiments."""
    region = gene_table.load_region_table()
    used = set(table.loc[~table["excluded"], "experiment_id"])
    per = gene_table.experiment_profiles(region[region["experiment_id"].isin(used)])
    rel = gene_table.gene_reliability(per).set_index("symbol")
    return rel["reliability"].dropna().to_dict(), rel["median_energy"].to_dict()


def main():
    """Print the settings, then test the sets, the contrasts and localisation."""
    config.print_settings({})
    tables = OUT / "tables"
    figures = OUT / "figures"
    q = ISH_ANALYSIS["q"]

    # the sets, the context group, and every gene's rho and null with both maps
    table = gene_table.load_gene_table()
    genes = gene_table.per_gene(table)
    members = gene_sets.members_from_table(genes)
    context = gene_table.context_members(sorted(genes["symbol"]), offline=True)
    ranking = gene_ranking.load_ranking()
    rho, nulls = {}, {}
    for name in gene_ranking.MAPS:
        rho[name] = ranking[ranking["map"] == name].set_index("symbol")["rho"]
        nulls[name], null_genes = gene_ranking.load_null_rho(name)
    shown = dict(members)
    shown[gene_sets.CONTEXT_SET] = context
    member_rows = gene_sets.member_table(shown, ranking)
    member_rows.to_csv(gene_sets.MEMBERS, index=False)
    sizes = ", ".join(f"{k} {len(v)}" for k, v in shown.items())
    print(f"sets: {sizes}")

    # each set's median against the surrogates, and the contrasts named in advance
    tests = gene_sets.set_tests(members, rho, nulls, null_genes)
    tests.to_csv(gene_sets.SET_TESTS, index=False)
    contrasts, contrast_nulls = gene_sets.contrast_tests(
        members, rho["nano"], nulls["nano"], null_genes
    )
    contrasts.to_csv(gene_sets.CONTRAST_TESTS, index=False)
    nano = tests[(tests["map"] == "nano") & tests["tested"]]
    print(
        f"set tests: {int((nano['q'] < q).sum())} of {len(nano)} sets past the null at "
        f"q < {q}; contrasts: "
        + "; ".join(
            f"{r['contrast']} {r['difference']:+.3f} (p {r['p_spatial']:.4f})"
            for _, r in contrasts.iterrows()
        )
    )

    # localisation against expression-matched postsynaptic controls
    declared = structures.declared_structures()
    profile = profiles.load_profile()
    merged = gene_table.load_profiles()
    common, composite = gene_sets.subunit_composite(members["subunits"], merged, declared)
    reliability, level = gene_levels(table)
    map_common = profile.loc[common, "zref_nano"].to_numpy(float)
    sides = {g: "localisation" for g in members["localisation"]}
    sides.update({g: "control" for g in members["other postsynaptic"]})
    loc_table = gene_sets.partial_rows(
        sides, common, map_common, composite, merged, reliability, level
    )
    pairs = gene_sets.matched_controls(loc_table, level)
    loc_table["matched_to"] = loc_table["symbol"].map(pairs).fillna("")
    summary, label_nulls = gene_sets.localisation_tests(loc_table, pairs)

    # the same design with the control pool of 5 October, the ontology panel's
    # postsynaptic-density genes, which also hold genes GO puts on both sides
    panel_controls = genes.loc[genes["ontology_role"] == "control_psd", "symbol"]
    panel_sides = {g: "localisation" for g in members["localisation"]}
    panel_sides.update({g: "control" for g in panel_controls})
    panel_table = gene_sets.partial_rows(
        panel_sides, common, map_common, composite, merged, reliability, level
    )
    panel_pairs = gene_sets.matched_controls(panel_table, level)
    panel_table["matched_to"] = panel_table["symbol"].map(panel_pairs).fillna("")
    panel_summary, _ = gene_sets.localisation_tests(panel_table, panel_pairs)
    panel_summary = panel_summary[
        panel_summary["test"].isin(["matched controls", "positive control"])
    ].copy()
    panel_summary["test"] = "5 October's controls, " + panel_summary["test"]
    summary = pd.concat([summary, panel_summary], ignore_index=True)

    # the matched difference against the map's surrogates
    surr, listed = spatial_null.load_surrogates("nano")
    if listed != declared:
        raise ValueError(
            "the surrogates were drawn on other structures than the declared set; "
            "run run_ish_spatial_null.py --recompute"
        )
    columns = [listed.index(s) for s in common]
    p_spatial, _ = gene_sets.matched_spatial(
        loc_table, pairs, surr[:, columns], common, composite, merged
    )
    summary["p_spatial"] = np.where(
        summary["test"] == "matched controls", p_spatial, np.nan
    )
    loc_table.insert(1, "pool", "other postsynaptic (GO)")
    panel_table.insert(1, "pool", "5 October (control_psd)")
    pd.concat([loc_table, panel_table], ignore_index=True).to_csv(
        gene_sets.LOCALISATION, index=False
    )
    summary.to_csv(gene_sets.LOCALISATION_SUMMARY, index=False)
    main_row = summary.iloc[0]
    positive = summary.set_index("test").loc["positive control"]
    print(
        f"localisation: {len(common)} structures with every subunit; "
        f"{int(main_row['n_first'])} against {int(main_row['n_second'])} matched, "
        f"difference {main_row['difference']:+.3f}, p {main_row['p_labels']:.4f}, "
        f"spatial p {p_spatial:.4f}; positive control p "
        f"{positive['p_labels']:.4f}"
    )

    # the numbers for the text
    numbers = numbers_table(shown, tests, contrasts, summary)
    numbers.to_csv(tables / "numbers_gene_sets.csv", index=False)

    # figures 09 and 10
    rules = dict(gene_sets.GENE_SETS)
    rules[gene_sets.CONTEXT_SET] = gene_sets.CONTEXT_RULE
    n_surrogates = surr.shape[0]
    fig = ish_plotting.plot_gene_sets(
        member_rows,
        tests,
        contrasts,
        contrast_nulls,
        list(gene_sets.SET_ORDER),
        gene_sets.CONTEXT_SET,
        rules,
        q,
        n_surrogates,
        save=figures / ish_plotting.figure_file("gene_sets"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_localisation(
        loc_table,
        summary,
        label_nulls,
        pairs,
        p_spatial,
        n_surrogates,
        panel_table=panel_table,
        save=figures / ish_plotting.figure_file("localisation"),
    )
    plt.close(fig)
    print(f"figures: 09 and 10 in {figures}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="analysis 3: kinds of genes")
    parser.parse_args()
    main()
