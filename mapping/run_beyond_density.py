"""Analysis 4 of the ISH line: what Gria1 and synapse density leave of the map.

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
    23. run_beyond_density     analysis 4: what Gria1 and synapse   <- this script
                               density leave; the leftover
    24. run_beyond_controls    seven attempts to break it
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

The main model of [beyond], nano rank ~ Gria1 rank + synapse-density rank, two
straight terms: the ceiling (how reproducible the map is), what Gria1 alone, synapse
density alone and both predict on structures the fit has not seen, split into what
only Gria1 predicts, what the two share, what only density predicts and what is
left, and the weight of each term; whether the leftover replicates across mice,
where it lives, and every gene of the gene table against it with the leftover's
spatial null. Synapse density is the mean rank of the genes run_density_markers
chose ([density_markers] chosen); the method is in sepmap/adult/beyond_density.py.
Writes, in adult_v2/ish_analysis/beyond/ under the data root:

    structures_used.csv        every structure of the adult table, in the fit or not,
                               why, and whether PSD95 is measured there
    variance_partition.csv     what each model predicts, against the ceiling
    partition.csv              the four parts of the reproducible map
    weights.csv                the two weights, map and terms z-scored
    replication.csv            per split of the adults, the two agreements
    residual_by_structure.csv  per structure: map, prediction, leftover, steadiness
    leftover_genes.csv         every gene against the leftover, with its spatial p
    leftover_sets.csv          the gene sets of analysis 3 against the leftover
    leftover_null.npz          the leftover's surrogates, every gene's null rho
    fig0_structures.png        working figures: the structures of the fit, the
    fig1_ceiling.png           ceiling, what each model predicts, the leftover's
    fig2_covariates.png        replication and where it is largest
    fig3_residual.png

    python run_beyond_density.py
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from sepmap import config, plotting
from sepmap.adult import beyond_density
from sepmap.adult import plotting as adult_plotting


def main():
    """Print the settings in force, run the steps of analysis 4, draw them."""
    config.print_settings({})
    found = beyond_density.main()

    # the working figures, one per step
    out = beyond_density.OUT
    reproducible = beyond_density.replicates(float(np.mean(found["map_agreement"])))
    figures = [
        adult_plotting.plot_structures_used(
            found["structures_used"], save=out / "fig0_structures.png"
        ),
        adult_plotting.plot_ceiling(
            found["map_agreement"],
            found["explainable"],
            reproducible,
            save=out / "fig1_ceiling.png",
        ),
        adult_plotting.plot_covariates(
            found["partition"], found["explainable"], save=out / "fig2_covariates.png"
        ),
        adult_plotting.plot_residual(
            found["map_agreement"],
            found["agreement"],
            found["implied"],
            found["residuals"],
            save=out / "fig3_residual.png",
        ),
    ]
    for fig in figures:
        plt.close(fig)
    print(f"\nworking figures: fig0_structures.png to fig3_residual.png in {out}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(
        description="analysis 4: what Gria1 and synapse density leave"
    )
    parser.parse_args()
    main()
