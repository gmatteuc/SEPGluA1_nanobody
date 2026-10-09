"""Seven attempts to break the result of run_beyond_density, and its other folds.

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
                               once); figure 14s2
    23. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    24. run_beyond_controls    seven attempts to break it           <- this script
    25. run_beyond_calibration the same model on maps whose answer
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

Spatial gradient, structure size, single animals, the whisker manipulation,
curvature, the whole gene space and the reading; and the main model under other
folds; described in sepmap/adult/beyond_controls.py (the check rows, which change
the model, run with the calibration at step 25). Writes, in
adult_v2/ish_analysis/beyond/ under the data root:

    controls.csv            one row per control, with its number and verdict
    gene_space.csv          control F, per number of components
    gene_space_calibration.csv  control F's own floor
    gene_space_summary.csv  control F's numbers
    readings.csv            control G, per reading
    folds.csv               the main model under its folds and others: the
                            shares of Gria1, density and both, the share left
    fig4_controls.png       working figures: controls A to D, E and F, and G
    fig5_model_space.png
    fig6_readings.png

    python run_beyond_controls.py
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt

from sepmap import config, plotting
from sepmap.adult import beyond_controls
from sepmap.adult import plotting as adult_plotting


def main():
    """Print the settings in force, run the seven controls and the other folds, draw."""
    config.print_settings({})
    found = beyond_controls.main()

    # the working figures: A to D, E and F, G
    out = beyond_controls.OUT
    passed = found["passed"]
    figures = [
        adult_plotting.plot_artefacts(
            found["leftover"],
            found["xyz"],
            found["sizes"],
            found["pairs"],
            found["naive"],
            found["rws"],
            passed,
            save=out / "fig4_controls.png",
        ),
        adult_plotting.plot_model_space(
            found["curve"],
            found["k_most"],
            found["explainable"],
            (found["linear"], found["cubic"], found["quintic"]),
            passed,
            save=out / "fig5_model_space.png",
        ),
        adult_plotting.plot_readings(
            found["readings"], passed["G"], save=out / "fig6_readings.png"
        ),
    ]
    for fig in figures:
        plt.close(fig)
    print(f"\nworking figures: fig4_controls.png to fig6_readings.png in {out}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(
        description="seven controls of analysis 4, and its other folds"
    )
    parser.parse_args()
    main()
