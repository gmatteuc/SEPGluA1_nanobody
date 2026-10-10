"""Per-mouse region means and the young against adult statistics, with no warping.

Python route, in run order (tools\\venv_atlas; run_closeup and run_adult_layers in
tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid
     3. run_cohort             cohort mean, SD and n per voxel
     4. run_compare            young against adult: maps, table
        run_replot             redraws those maps (not in a full run)
     5. run_region_plot        region statistics, no warping        <- this script
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
    31. run_adult_layers       the adult map by depth (venv_flat)

Every brain is measured on the atlas of its own age. Per structure and
reading: one value per mouse, Welch and Mann-Whitney tests of young against
adult, and their Benjamini-Hochberg q. These are the numbers to quote; the
method is in sepmap/young_vs_adult/region_plot.py. Writes, in
comparisons_v2/young_vs_adult/ under the data root:

    region_means_per_mouse.csv    one row per reading, mouse and structure
    region_stats.csv              one row per reading and structure: the tests of
                                  young against adult, their q, the P20-only and
                                  naive-rws contrasts
    region_plot.png               a panel per reading, a dot per mouse

and an .eps beside the PNG.

    python run_region_plot.py

The readings are readings.in_force of settings.toml; V2_READINGS, a
comma-separated subset of them (ratio, sepratio, cref, subref, zref), replaces it
for one run and leaves the others out of this step. The run prints the readings in
force.
"""

import argparse

import matplotlib

from sepmap import config
from sepmap.volumes import cohort
from sepmap.young_vs_adult import region_plot


def main():
    """Print the readings in force, then measure the region statistics per mouse."""
    # settings in force
    readings = cohort.readings_in_force()
    config.print_settings({"readings": readings})

    # region means, tests and dot plot
    region_plot.main()


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="per-mouse region statistics")
    parser.parse_args()
    main()
