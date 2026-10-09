"""The model of analysis 4 on maps whose answer is known, and its check rows.

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
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation
                               against matched controls; figures 10,
                               10s1, 10s2
    21. run_synaptome          the measured synapse density (network,
                               once): PSD95 puncta per structure,
                               its coverage of the fit; figure 14s1
    22. run_density_markers    the synapse-density genes, chosen by
                               PSD95 without the map (network,
                               once); figure 03s3
    23. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_calibration the same model on maps whose answer   <- this script
                               is known; the check rows, each
                               with its own floor
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

Each gene's Allen experiments split in two halves; a map made only of Gria1 and
synapse density, the main model fitted to nano with one half's predictors, given
ten made-up adults as noisy as ours and predicted from the other half: the floor.
The real map is read the same way, on the same structures, and their difference
resampled over structures. Then each check row of [beyond.checks], the main model
changed one thing at a time, on its own structures and with its own floor. The
methods are in sepmap/adult/beyond_calibration.py and beyond_checks.py. Writes, in
adult_v2/ish_analysis/beyond/ under the data root:

    calibration.csv             per folds (random, spatial blocks), map, direction
                                and draw: the ceiling, the CV R2, the share left
                                and the leftover's replication
    calibration_jackknife.csv   per subsample of the calibration structures: the
                                nano map's share left, the floor's, the difference
    check_rows.csv              per check row: its structures and terms, the
                                shares, the share left, its floor and nano minus
                                the floor with its interval
    check_rows_calibration.csv  every check row's calibration
    check_rows_jackknife.csv    every check row's paired jackknife

    python run_beyond_calibration.py
"""

import argparse

import matplotlib

from sepmap import config, plotting
from sepmap.adult import beyond_calibration, beyond_checks


def main():
    """Print the settings in force, then run the calibration and the check rows."""
    config.print_settings({})
    beyond_calibration.main()
    beyond_checks.main()


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(
        description="analysis 4 on maps whose answer is known, and its check rows"
    )
    parser.parse_args()
    main()
