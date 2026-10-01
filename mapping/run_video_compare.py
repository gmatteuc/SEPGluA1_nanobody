"""Young beside adult, plane by plane: a video per reading, or one plane as a still.

Three panels per frame: the young and the adult cohort mean on one colour
scale, and their difference. Writes into comparisons_v2/young_vs_adult/; the
layout is described in sepmap/young_vs_adult/video_compare.py.

Python route, in run order (tools\\venv_atlas; run_closeup in tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid
     3. run_cohort             cohort mean, SD and n per voxel
     4. run_compare            young against adult: maps, table
        run_replot             redraws those maps (not in a full run)
     5. run_region_plot        region statistics, no warping
     6. run_region_groups      the same by system and by layer
     7. run_video              cohort videos
     8. run_video_compare      young beside adult, plane by plane     <- this script
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

    python run_video_compare.py [reading ...] [--plane P] [--vmax V] [--dlim D]

With no reading named, every reading in force. --plane draws that CCF plane
(10 um numbering) as a PNG instead of the video; --vmax and --dlim set the
colour range of the two means and of the difference for this run only.
Readings and options can come in any order. V2_READINGS limits the readings
in force, as for run_cohort; the run prints them.
"""

import argparse
import os

from sepmap import config
from sepmap.volumes import cohort
from sepmap.young_vs_adult import video_compare


def main(readings, plane=None, vmax=None, dlim=None):
    """Young beside adult for `readings`, every reading in force when it is empty."""
    readings = readings or list(cohort.MODES)

    # settings in force
    in_force = " ".join(cohort.MODES)
    if os.environ.get("V2_READINGS", "").strip():
        in_force += "  (from V2_READINGS)"
    config.print_settings({
        "readings": " ".join(readings),
        "in force": in_force,
        "plane": "every plane, as a video" if plane is None else f"{plane}, as a still",
        "vmax": "per reading" if vmax is None else vmax,
        "dlim": "default" if dlim is None else dlim,
    })

    # young beside adult
    video_compare.main(readings, plane=plane, vmax=vmax, dlim=dlim)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="young beside adult, plane by plane")
    parser.add_argument("readings", nargs="*", help="readings (default: all in force)")
    parser.add_argument("--plane", type=int, help="one CCF plane (10 um) as a still")
    parser.add_argument("--vmax", type=float, help="colour range of the two means")
    parser.add_argument("--dlim", type=float, help="colour range of the difference")
    args = parser.parse_intermixed_args()
    main(args.readings, plane=args.plane, vmax=args.vmax, dlim=args.dlim)
