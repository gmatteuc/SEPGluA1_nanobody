"""One reading looked at closely: a coronal plane, its video, the cortical flatmaps.

All four views come from one set of prepared cohort volumes, with the colour
range tightened for cortex and a light smoothing that every title declares.
Writes into comparisons_v2/young_vs_adult/, or into a subfolder named after
--cmap. The flatmaps need ccf_streamlines, so this runs in tools\\venv_flat
and imports only sepmap.config and its own module; the method, and the
flatmap assets it needs in atlas_flatmap/, are described in
sepmap/young_vs_adult/closeup.py.

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
     9. run_closeup            close-ups and flatmaps (venv_flat)     <- this script
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

    python run_closeup.py [reading ...] [--plane 790] [--vmax V] [--dlim D]
        [--smooth 3,1,1] [--cmap NAME] [--no-video] [--no-flatmap]

With no reading named, zref. --plane is a CCF plane at 10 um (default 790,
where RL and AL are cut). --vmax and --dlim override the colour ranges of the
means and of the difference. --smooth is one sigma or three, comma
separated, in 20 um voxels along (AP, DV, ML); 0 turns the smoothing off.
--cmap draws the mean panels with another matplotlib colormap, into a
subfolder of that name. --no-video and --no-flatmap leave those out.
Readings and options can come in any order.
"""

import argparse

from sepmap import config
from sepmap.young_vs_adult import closeup


def smoothing(text):
    """The --smooth value as a list of sigmas: one number, or three comma separated."""
    return [float(x) for x in text.split(",")]


def main(readings, plane, vmax, dlim, smooth, want_video, want_flatmap, cmap_name):
    """Close-up views of `readings`, zref when the list is empty.

    `smooth` is a list of one or three sigmas, or None for the default.
    """
    readings = readings or ["zref"]
    sigmas = smooth if smooth is not None else list(closeup.SMOOTH)
    sigma = sigmas[0] if len(sigmas) == 1 else sigmas

    # settings in force
    config.print_settings({
        "readings": " ".join(readings),
        "plane": plane,
        "vmax": "per reading" if vmax is None else vmax,
        "dlim": "per reading" if dlim is None else dlim,
        "smooth": sigma,
        "cmap": "default" if cmap_name is None else cmap_name,
        "video": "yes" if want_video else "no",
        "flatmap": "yes" if want_flatmap else "no",
    })

    # close-up views
    closeup.main(readings, plane=plane, vmax=vmax, dlim=dlim, sigma=sigma,
                 want_video=want_video, want_flatmap=want_flatmap, cmap_name=cmap_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="one reading looked at closely")
    parser.add_argument("readings", nargs="*", help="readings to draw (default: zref)")
    parser.add_argument("--plane", type=int, default=closeup.PLANE,
                        help=f"CCF plane at 10 um (default {closeup.PLANE})")
    parser.add_argument("--vmax", type=float, help="colour range of the means")
    parser.add_argument("--dlim", type=float, help="colour range of the difference")
    parser.add_argument("--smooth", type=smoothing,
                        help="sigma in 20 um voxels: one, or three comma separated")
    parser.add_argument("--cmap", help="matplotlib colormap of the mean panels")
    parser.add_argument("--no-video", action="store_true", help="no video")
    parser.add_argument("--no-flatmap", action="store_true", help="no flatmaps")
    args = parser.parse_intermixed_args()
    main(args.readings, args.plane, args.vmax, args.dlim, args.smooth,
         want_video=not args.no_video, want_flatmap=not args.no_flatmap,
         cmap_name=args.cmap)
