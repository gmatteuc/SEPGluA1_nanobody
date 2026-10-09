"""The three channel arms of the measurement argument, one region table per adult.

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
                               documentation; figure 02
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the
                               map, the null, autofluorescence (A8),
                               the Cacng8 - Gria1 gap; figures 05 to 07 and 15
    18. run_ish_robustness     A3: the ranking under other choices;
                               figure 16
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figure 10, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation
                               against matched controls; figures 08, 09
    21. run_synaptome          the measured synapse density (network,
                               once): PSD95 puncta per structure,
                               its coverage of the fit
    22. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    23. run_beyond_controls    seven attempts to break it
    24. run_beyond_calibration the same model on maps whose answer
                               is known
    25. run_beyond_regression  the regression, per structure
    26. run_beyond_figures     figures 03, 04 and 11
    27. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 12 to 14, sheets
    28. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 17
    29. run_ish_overview       figures 00 and 18, the figure index;
                               the numbers for the text

Not in the run order: no step of the ISH line reads its table, and it writes
adult_v2/arms/, the frozen run of 5 October, so a rerun would overwrite it. Kept
until Giulio decides whether it retires with run_ish_arms.py (docs/ISH_ANALYSIS.md,
Points to settle).

SEP/auto, nano/auto and nano/SEP per structure for the ten adults, checked
against run_region_plot's table for the two arms both compute. The method, and
why the SEP arms do not mean what their names promise, are in
sepmap/adult/arms.py. Writes, in adult_v2/arms/ under the data root:

    region_means_arms.csv    arm x mouse x structure, log2
    arms_consistency.png     the self-check: against young_vs_adult.region_plot,
                             and the log-space identity between the three arms

    python run_adult_arms.py
"""

import argparse

import matplotlib

from sepmap import config
from sepmap.adult import arms


def main():
    """Print the settings in force, then write the arms table and its self-check."""
    # settings in force
    config.print_settings({})

    # region table per arm and adult
    arms.main()


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="the channel arms per adult")
    parser.parse_args()
    main()
