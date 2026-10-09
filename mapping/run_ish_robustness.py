"""The gene ranking of the ISH line under other reasonable choices (A3).

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
    18. run_ish_robustness     A3: the ranking under other choices;   <- this script
                               figures 13, 13s
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figures 09, 09s, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation
                               against matched controls; figures 10,
                               10s1, 10s2
    21. run_synaptome          the measured synapse density (network,
                               once): PSD95 puncta per structure,
                               its coverage of the fit; figure 03s4
    22. run_density_markers    the synapse-density genes, chosen by
                               PSD95 without the map (network,
                               once); figure 03s3
    23. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_calibration the same model on maps whose answer
                               is known; the check rows, each
                               with its own floor
    26. run_beyond_regression  the regression, per structure
    27. run_beyond_figures     figures 03, 03s1, 03s2, 04, 11s1
    28. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 07, 08, 11, 11s2,
                               sheets
    29. run_sep_channel_check  analysis 5: what the green channel
                               reports; figures 14, 14s
    30. run_ish_overview       figures 00, 15 and 15s, the figure
                               index; the numbers for the text

Correlates every gene with the adult map again, changing one choice of the primary
ranking at a time (statistic, borders, reading, inputs, structure set) and once all
together as the route of 5 October ran, and compares each gene order with the
primary's (the method is in sepmap/ish/robustness.py). Writes, in
adult_v2/ish_analysis/ under the data root:

    tables/ranking_robustness.csv     per variant and gene: rho, structures, ranks
    tables/robustness_summary.csv     per variant: agreement with the primary over
                                      P9's genes and over all, Cacng8's and Gria1's
                                      rho and ranks, the gap between them
    tables/numbers_robustness.csv     the numbers of this step, for the text
    figures/13_robustness.png         the ranking under each choice (13s in
                                      detail)

    python run_ish_robustness.py

The row without section QC measures the flagged experiments again with nothing set
missing (about a minute).
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt

from sepmap import config, plotting, structures
from sepmap.adult import profiles
from sepmap.ish import gene_ranking, gene_sets, gene_table, robustness
from sepmap.ish.figure_index import figure_file, figure_path
from sepmap.ish.plotting import ranking as ranking_figures


def main():
    """Print the settings, then rank the genes under every variant; tables, figures."""
    config.print_settings({})

    # the inputs of the primary and of every variant
    declared = structures.declared_structures()
    set_table = structures.load_structure_set()
    profile = profiles.load_profile()
    table = gene_table.load_gene_table()
    genes = gene_table.per_gene(table)
    p9_genes = set(genes.loc[genes["p9_gene"], "symbol"])
    subunits = set(gene_sets.members_from_table(genes)["subunits"])
    inputs = robustness.variant_inputs(table, genes, profile, set_table, declared)
    print(f"variants: {len(inputs)}, {len(p9_genes)} of P9's genes among {len(genes)}")

    # every variant's rho per gene, and the summary: agreement, where the two genes
    # sit, the gap
    per_gene, summary = robustness.variant_tables(inputs, p9_genes)
    per_gene.to_csv(robustness.ROBUSTNESS, index=False)
    summary.to_csv(robustness.SUMMARY, index=False)
    for _, r in summary.iterrows():
        print(
            f"  {r['variant']:16s} agreement {r['agreement_p9']:.3f} over "
            f"{r['n_agreement_p9']} of P9's genes; Gria1 {r['rank_p9_Gria1']:.0f}, "
            f"Cacng8 {r['rank_p9_Cacng8']:.0f}, gap {r['gap']:+.3f}"
        )
    robustness.numbers_table(summary).to_csv(robustness.NUMBERS, index=False)

    # figure 13 and its detailed version, behind the gaps the band within which the
    # primary's two-sided test, fixed in advance, does not pass
    merged = gene_ranking.merged_gap(gene_ranking.load_gap())
    band = (-float(merged["equal_two_sided"]), float(merged["equal_two_sided"]))
    fig = ranking_figures.plot_robustness(summary, band, save=figure_path("robustness"))
    plt.close(fig)
    fig = ranking_figures.plot_robustness_detail(
        summary, per_gene, subunits, band, save=figure_path("robustness_detail")
    )
    plt.close(fig)
    drawn = [figure_file(k) for k in ("robustness", "robustness_detail")]
    print(f"figures: {', '.join(drawn)} in {structures.FIGURES}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="the gene ranking under other choices")
    parser.parse_args()
    main()
