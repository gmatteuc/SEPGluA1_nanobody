"""Every brain onto the adult CCF grid at 20 um, one brain at a time.

Young brains are carried from the DeMBA atlas of their own age to the adult
CCF (brainglobe_ccf_translator, minutes per channel); adults, registered to the
CCF already, are only placed. Reads comparisons_v2/per_mouse/ and writes
comparisons_v2/per_mouse_ccf/<mouse>.npz; the method is in
sepmap/volumes/to_ccf.py.

Python route, in run order (tools\\venv_atlas; run_closeup in tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid      <- this script
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
    23. run_beyond_density     what abundance and density leave
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_figures     the figures of that result
    26. run_beyond_regression  the regression, shown

    python run_to_ccf.py [mouse ...]

With no mouse named, every brain of the cohort.
"""

import argparse

from sepmap import config
from sepmap.volumes import per_mouse, to_ccf


def main(mice):
    """Carry `mice`, or every brain when the list is empty, onto the CCF grid."""
    mice = mice or list(per_mouse.MICE)

    # settings in force
    config.print_settings({"mice": " ".join(mice)})

    # every brain on the adult CCF grid
    to_ccf.main(mice)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="every brain onto the adult CCF grid")
    parser.add_argument("mice", nargs="*", help="mice to process (default: every brain)")
    args = parser.parse_args()
    main(args.mice)
