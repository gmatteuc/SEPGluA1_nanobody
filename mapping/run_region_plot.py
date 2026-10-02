"""Per-mouse region means and the young against adult statistics, with no warping.

Every brain is measured on the atlas of its own age. Per structure and
reading: one value per mouse, Welch and Mann-Whitney tests of young against
adult, and their Benjamini-Hochberg q. Writes region_means_per_mouse.csv,
region_stats.csv and region_plot.png into comparisons_v2/young_vs_adult/.
These are the numbers to quote; the method is in
sepmap/young_vs_adult/region_plot.py.

Python route, in run order (tools\\venv_atlas; run_closeup in tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid
     3. run_cohort             cohort mean, SD and n per voxel
     4. run_compare            young against adult: maps, table
        run_replot             redraws those maps (not in a full run)
     5. run_region_plot        region statistics, no warping          <- this script
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
    23. run_beyond_density     what abundance and density leave
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_figures     the figures of that result
    26. run_beyond_regression  the regression, shown

    python run_region_plot.py

V2_READINGS, a comma-separated subset of the readings (ratio, sepratio,
cref, subref, zref), leaves the others out of this step; the run prints the
readings in force.
"""

import argparse
import os

from sepmap import config
from sepmap.volumes import cohort
from sepmap.young_vs_adult import region_plot


def main():
    """Per-mouse region statistics."""
    # settings in force
    readings = " ".join(cohort.MODES)
    if os.environ.get("V2_READINGS", "").strip():
        readings += "  (from V2_READINGS)"
    config.print_settings({"readings": readings})

    # region means, tests and dot plot
    region_plot.main()


if __name__ == "__main__":
    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="per-mouse region statistics")
    parser.parse_args()
    main()
