"""What the green (SEP) channel reports in this tissue, from the raw channels.

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
    17. run_sep_channel_check  what the green channel reports         <- this script
    18. run_panel_build        the 390-gene ontology panel
    19. run_panel_fetch        its ISH grids (network, once)
    20. run_ish_regions        ISH per structure, ontology panel
    21. run_ish_reliability    how reliable one ISH map is
    22. run_ish_panel_test     localisation against controls
    23. run_beyond_density     what abundance and density leave
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_figures     the figures of that result
    26. run_beyond_regression  the regression, shown

Per adult: the dynamic range of each channel, what each tracks across
structures, and the SEP residual once autofluorescence is regressed out; the
method is in sepmap/adult/sep_channel_check.py. Writes, in adult_v2/arms/ under
the data root:

    sep_channel_check.csv    the per-mouse numbers behind the figure
    sep_channel_check.png    the figure

    python run_sep_channel_check.py
"""

import argparse

import matplotlib

from sepmap import config
from sepmap.adult import sep_channel_check


def main():
    """Print the settings in force, then measure what the green channel reports."""
    # settings in force
    config.print_settings({})

    # channel check
    sep_channel_check.main()


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="what the green channel reports")
    parser.parse_args()
    main()
