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
    11. run_panel_build        390-gene ontology panel (network, cached)
    12. run_panel_fetch        its ISH grids (network, once)
    13. run_structure_set      A1: the declared structures, the
                               adult profiles; figure 01
    14. run_ish_section_qc     A2: the experiments of both panels
                               and the repair (network, once);
                               section QC; QC sheets
    15. run_ish_gene_table     A9: region means, the gene table,
                               merged profiles, gene sets,
                               documentation; figures 02, 02s
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the
                               map, the null, autofluorescence (A8),
                               the Cacng8 - Gria1 gap; figures 05,
                               06, 07s and 12, with their s
    18. run_ish_robustness     A3: the ranking under other choices;
                               figures 13, 13s
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figures 09, 09s, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation   <- this script
                               against matched controls; figures 10,
                               10s1, 10s2
    21. run_synaptome          the measured synapse density (network,
                               once): PSD95 puncta per structure,
                               its coverage of the fit; figure 14s1
    22. run_density_markers    the synapse-density genes, chosen by
                               PSD95 without the map (network,
                               once); figure 14s2
    23. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_calibration the same model on maps whose answer
                               is known
    26. run_beyond_regression  the regression, per structure
    27. run_beyond_figures     figures 03, 03s1, 03s2, 04, 11s1, 14
    28. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 07, 08, 11, 11s2,
                               sheets
    29. run_sep_channel_check  analysis 5: what the green channel
                               reports; figures 15, 15s
    30. run_ish_overview       figures 00, 16 and 16s, the figure
                               index; the numbers for the text

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
                                       control pool of 5 October (secondary)
    tables/localisation_power.csv      what the test finds when localisation genes
                                       do shape the map, per effect size
    tables/numbers_gene_sets.csv       the numbers of this step, for the text
    figures/10_gene_kinds.png          the sets, and localisation against
                                       matched controls
    figures/10s1_gene_sets.png         the sets in detail
    figures/10s2_localisation.png      the localisation test in detail

    python run_ish_gene_sets.py
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from sepmap import config, plotting, structures
from sepmap.ish import gene_ranking, gene_sets, gene_table
from sepmap.ish import plotting as ish_plotting
from sepmap.ish.figure_index import figure_file, figure_path

ISH_ANALYSIS = config.SETTINGS["ish_analysis"]


def draw_figures(member_rows, tests, contrasts, contrast_nulls, local):
    """Figure 10, and its detailed versions: the gene sets, the localisation test.

    `local` holds what gene_sets.localisation_design returns.
    """
    q = ISH_ANALYSIS["q"]
    rules = dict(gene_sets.GENE_SETS)
    rules[gene_sets.CONTEXT_SET] = gene_sets.CONTEXT_RULE
    _, term_names, _ = gene_table.load_obo(offline=True)
    rules = {k: gene_table.with_go_names(v, term_names) for k, v in rules.items()}
    n_surrogates = local["n_surrogates"]
    fig = ish_plotting.plot_gene_kinds(
        member_rows,
        tests,
        list(gene_sets.SET_ORDER),
        local["loc_table"],
        local["summary"],
        local["pairs"],
        local["power"],
        q,
        n_surrogates,
        save=figure_path("gene_kinds"),
    )
    plt.close(fig)
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
        save=figure_path("gene_sets"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_localisation(
        local["loc_table"],
        local["summary"],
        local["label_nulls"],
        local["pairs"],
        local["p_spatial"],
        n_surrogates,
        power=local["power"],
        panel_table=local["panel_table"],
        save=figure_path("localisation"),
    )
    plt.close(fig)
    drawn = [figure_file(k) for k in ("gene_kinds", "gene_sets", "localisation")]
    print(f"figures: {', '.join(drawn)} in {structures.FIGURES}")


def main():
    """Print the settings, then test the sets, the contrasts and localisation."""
    config.print_settings({})
    q = ISH_ANALYSIS["q"]

    # the sets, the context group, and every gene's rho and null with both maps
    table = gene_table.load_gene_table()
    genes = gene_table.per_gene(table)
    members = gene_sets.members_from_table(genes)
    context = gene_table.context_members(sorted(genes["symbol"]), offline=True)
    ranking = gene_ranking.load_ranking()
    rho, nulls = {}, {}
    for name in gene_ranking.MAPS:
        rho[name] = gene_ranking.map_rows(ranking, name)["rho"]
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
    named = "; ".join(
        f"{r['contrast']} {r['difference']:+.3f} (p {r['p_spatial']:.4f})"
        for _, r in contrasts.iterrows()
    )
    print(
        f"set tests: {int((nano['q'] < q).sum())} of {len(nano)} sets past the null at "
        f"q < {q}; contrasts: {named}"
    )

    # localisation against expression-matched postsynaptic controls
    merged = gene_table.load_profiles()
    reliability, level = gene_table.gene_levels(table)
    local = gene_sets.localisation_design(members, genes, merged, reliability, level)
    power = local["power"]
    power.to_csv(gene_sets.LOCALISATION_POWER, index=False)
    print(
        "power: label test p < 0.05 on "
        f"{power['power_labels'].iloc[0]:.1%} of maps with no effect; 80% found at a "
        f"matched difference of {gene_sets.detectable(power, 'power_labels'):+.3f} "
        f"(labels), {gene_sets.detectable(power, 'power_spatial'):+.3f} (spatial)"
    )
    pd.concat([local["loc_table"], local["panel_table"]], ignore_index=True).to_csv(
        gene_sets.LOCALISATION, index=False
    )
    summary = local["summary"]
    summary.to_csv(gene_sets.LOCALISATION_SUMMARY, index=False)
    by_test = summary.set_index("test")
    matched = by_test.loc["matched controls"]
    positive = by_test.loc["positive control"]
    print(
        f"localisation: {len(local['common'])} structures with every subunit; "
        f"{int(matched['n_first'])} against {int(matched['n_second'])} matched, "
        f"difference {matched['difference']:+.3f}, p {matched['p_labels']:.4f}, "
        f"spatial p {local['p_spatial']:.4f}; positive control p "
        f"{positive['p_labels']:.4f}"
    )

    # the numbers for the text
    numbers = gene_sets.numbers_table(shown, tests, contrasts, summary, power)
    numbers.to_csv(gene_sets.NUMBERS, index=False)

    # figure 10, and its detailed versions
    draw_figures(member_rows, tests, contrasts, contrast_nulls, local)


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="analysis 3: kinds of genes")
    parser.parse_args()
    main()
