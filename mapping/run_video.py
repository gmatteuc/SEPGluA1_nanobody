"""Cohort videos: the mean and its reliability t, plane by plane, per reading.

Python route, in run order (tools\\venv_atlas; run_closeup in tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid
     3. run_cohort             cohort mean, SD and n per voxel
     4. run_compare            young against adult: maps, table
        run_replot             redraws those maps (not in a full run)
     5. run_region_plot        region statistics, no warping
     6. run_region_groups      the same by system and by layer
     7. run_video              cohort videos                          <- this script
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

The layout and the colour ranges are described in sepmap/young_vs_adult/video.py.
Writes, in comparisons_v2/ccf/<cohort>/ under the data root:

    video_<reading>_<cohort>.mp4    the mean and its reliability t, plane by plane

    python run_video.py [cohort ...]

With no cohort named: young, adult, young_P20, naive and rws.

V2_READINGS, a comma-separated subset of the readings (ratio, sepratio,
cref, subref, zref), leaves the others out of this step; the run prints the
readings in force.
"""

import argparse
import os

import matplotlib

from sepmap import config
from sepmap.volumes import cohort
from sepmap.young_vs_adult import video

DEFAULT_COHORTS = ["young", "adult", "young_P20", "naive", "rws"]


def main(cohorts):
    """Print the settings in force, then write the videos of `cohorts`.

    `cohorts` empty means the default five, DEFAULT_COHORTS.
    """
    cohorts = cohorts or DEFAULT_COHORTS

    # settings in force
    readings = " ".join(cohort.MODES)
    if os.environ.get("V2_READINGS", "").strip():
        readings += "  (from V2_READINGS)"
    config.print_settings({"cohorts": " ".join(cohorts), "readings": readings})

    # one video per cohort and reading
    video.main(cohorts)


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")

    parser = argparse.ArgumentParser(description="cohort videos, plane by plane")
    parser.add_argument(
        "cohorts",
        nargs="*",
        help="cohorts to draw (default: young adult young_P20 naive rws)",
    )
    args = parser.parse_args()
    main(args.cohorts)
