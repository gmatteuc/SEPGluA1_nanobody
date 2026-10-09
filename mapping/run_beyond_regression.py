"""The regression of analysis 4, structure by structure.

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
    25. run_beyond_calibration the same model on maps whose answer
                               is known; the check rows, each
                               with its own floor
    26. run_beyond_regression  the regression, per structure        <- this script
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

The fit of run_beyond_density drawn as a fit: the map against the prediction,
the residual diagnostic, and the three maps on the brain; the method is in
sepmap/adult/beyond_regression.py. Writes, in adult_v2/ish_analysis/beyond/ under
the data root:

    regression_table.csv    every structure: map, prediction, leftover (ranks)
    E_regression.png        working figures: the fit and its diagnostic, and the
    F_maps.png              three maps on the brain

    python run_beyond_regression.py
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt

from sepmap import config, plotting
from sepmap.adult import beyond_regression
from sepmap.adult import plotting as adult_plotting

# the p at which the diagnostic calls the model mis-specified, and the floor of the
# rank maps' colour scale
BEYOND_REGRESSION = config.SETTINGS["beyond_regression"]


def main():
    """Print the settings in force, fit the main model, draw the fit and its maps."""
    config.print_settings({})
    found = beyond_regression.main()

    # the working figures: the fit and its diagnostic, the three maps on the brain
    out = beyond_regression.OUT
    fig = adult_plotting.plot_regression(
        found["observed"],
        found["predicted"],
        found["leftover"],
        found["structures"],
        found["fitted"],
        found["cv"],
        found["misspecified"],
        BEYOND_REGRESSION["diagnostic_p"],
        save=out / "E_regression.png",
    )
    plt.close(fig)
    fig = adult_plotting.plot_regression_maps(
        beyond_regression.coronal_planes(),
        found["observed"],
        found["predicted"],
        found["leftover"],
        found["structures"],
        BEYOND_REGRESSION["rank_floor"],
        save=out / "F_maps.png",
    )
    plt.close(fig)
    print(f"working figures: E_regression.png, F_maps.png in {out}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="the regression of analysis 4")
    parser.parse_args()
    main()
