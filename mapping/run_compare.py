"""Young against adult in the adult CCF: maps and the per-structure table.

Python route, in run order (tools\\venv_atlas; run_closeup in tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid
     3. run_cohort             cohort mean, SD and n per voxel
     4. run_compare            young against adult: maps, table       <- this script
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
    23. run_beyond_density     what abundance and density leave
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_figures     the figures of that result
    26. run_beyond_regression  the regression, shown

Folds the hemispheres of the cohort volumes and compares the young cohort with
the adults reading by reading (the method is in sepmap/young_vs_adult/compare.py).
Writes, in comparisons_v2/young_vs_adult/ under the data root:

    volumes_ccf20.npz     adult_*, young_*, log2_* maps, n maps, annot20 (AP, DV, ML half)
    slices_<reading>.png  dorsal up, midline right, no data in grey; an .eps beside it
    region_table.csv      per structure, every reading, both young groups
    cortex_table.txt      the cortical areas, also printed

    python run_compare.py

The readings are readings.in_force of settings.toml; V2_READINGS, a
comma-separated subset of them (ratio, sepratio, cref, subref, zref), replaces it
for one run and leaves the others out of this step. The run prints the readings in
force.
"""

import argparse

import matplotlib

from sepmap import config
from sepmap.volumes import cohort
from sepmap.young_vs_adult import compare


def main():
    """Print the readings in force, then compare and write the maps and tables."""
    # settings in force
    readings = cohort.readings_in_force()
    config.print_settings({"readings": readings})

    # maps and tables
    compare.main()


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="young against adult in the adult CCF")
    parser.parse_args()
    main()
