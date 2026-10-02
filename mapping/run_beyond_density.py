"""What receptor abundance and synaptic density leave unexplained in the adult map.

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
    11. run_ish_regions        ISH per structure, 100-gene panel
    12. run_ish_compare        the adult map against each gene
    13. run_ish_words          annotation words of the ranking
    14. run_ish_roles          subunit against localisation genes
    15. run_adult_arms         the channel arms per adult
    16. run_ish_arms           the arms against the genes
    17. run_sep_channel_check  what the green channel reports
    18. run_panel_build        the 390-gene ontology panel
    19. run_panel_fetch        its ISH grids (network, once)
    20. run_ish_regions        ISH per structure, ontology panel
    21. run_ish_reliability    how reliable one ISH map is
    22. run_ish_panel_test     localisation against controls
    23. run_beyond_density     what abundance and density leave       <- this script
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_figures     the figures of that result
    26. run_beyond_regression  the regression, shown

The ceiling (how reproducible the map is), what the covariates explain, and
whether the residual replicates across independent halves of the cohort; the
argument is in sepmap/adult/beyond_density.py. Writes, in adult_v2/beyond/ under
the data root:

    structures_used.csv        every structure, kept or dropped, with the reason
    variance_partition.csv     what each covariate set explains, against the ceiling
    residual_by_structure.csv  where the map exceeds and falls short of prediction
    fig0_structures.png, fig1_ceiling.png, fig2_covariates.png, fig3_residual.png

    python run_beyond_density.py
"""

import argparse

import matplotlib

from sepmap import config
from sepmap.adult import beyond_density


def main():
    """Print the settings in force, then run the steps of the analysis."""
    # settings in force
    config.print_settings({})

    # ceiling, covariates, residual
    beyond_density.main()


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(
        description="what abundance and density leave unexplained"
    )
    parser.parse_args()
    main()
