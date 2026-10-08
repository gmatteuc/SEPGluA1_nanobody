"""One reading looked at closely: a coronal plane, its video, the cortical flatmaps.

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
     9. run_closeup            close-ups and flatmaps (venv_flat)   <- this script
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
                               documentation; figure 02
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the
                               map, the null, autofluorescence (A8),
                               the Cacng8 - Gria1 gap; figures 05 to 07 and 11
    18. run_ish_robustness     A3: the ranking under other choices;
                               figure 12
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figure 10, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation
                               against matched controls; figures 08, 09
    21. run_beyond_density     analysis 4: what receptor mRNA and
                               synaptic density leave; the leftover
    22. run_beyond_controls    seven attempts to break it
    23. run_beyond_calibration the same model on maps whose answer
                               is known
    24. run_beyond_regression  the regression, per structure
    25. run_beyond_figures     figures 03 and 04
    26. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 13
    27. run_ish_overview       figures 00 and 14, the figure index;
                               the numbers for the text

All four views come from one set of prepared cohort volumes, with the colour
range tightened for cortex and a light smoothing that every title declares. The
flatmaps need ccf_streamlines, so this runs in tools\\venv_flat and imports from
the package only sepmap.config and its own module, which needs only plotting and
hemispheres; the method, and the flatmap assets it needs in atlas_flatmap/, are
described in sepmap/young_vs_adult/closeup.py. Its defaults are the closeup table
of settings.toml. Writes, in
comparisons_v2/young_vs_adult/ under the data root, or in a subfolder named after
--cmap:

    detail_plane<P>_<reading>.png        one coronal plane: young, adult, difference
    detail_video_<reading>.mp4           the same three panels, plane by plane
    detail_flatmap_<reading>.png         the isocortex unrolled, full cortical depth
    detail_flatmap_layers_<reading>.png  the same by depth band

and an .eps beside each PNG.

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

import matplotlib

from sepmap import config
from sepmap.young_vs_adult import closeup

CLOSEUP = config.SETTINGS["closeup"]


def smoothing(text):
    """The --smooth value as a list of sigmas: one number, or three comma separated."""
    return [float(x) for x in text.split(",")]


def main(readings, plane, vmax, dlim, smooth, want_video, want_flatmap, cmap_name):
    """Print the settings in force, then draw the close-up views of `readings`.

    `readings` empty means zref; `smooth` is a list of one or three sigmas, or None
    for the default.
    """
    readings = readings or ["zref"]
    sigmas = smooth if smooth is not None else list(CLOSEUP["smooth"])
    sigma = sigmas[0] if len(sigmas) == 1 else sigmas

    # settings in force
    config.print_settings(
        {
            "readings": " ".join(readings),
            "plane": plane,
            "vmax": "per reading" if vmax is None else vmax,
            "dlim": "per reading" if dlim is None else dlim,
            "smooth": sigma,
            "cmap": "default" if cmap_name is None else cmap_name,
            "video": "yes" if want_video else "no",
            "flatmap": "yes" if want_flatmap else "no",
        }
    )

    # close-up views
    closeup.main(
        readings,
        plane=plane,
        vmax=vmax,
        dlim=dlim,
        sigma=sigma,
        want_video=want_video,
        want_flatmap=want_flatmap,
        cmap_name=cmap_name,
    )


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")

    parser = argparse.ArgumentParser(description="one reading looked at closely")
    parser.add_argument("readings", nargs="*", help="readings to draw (default: zref)")
    parser.add_argument(
        "--plane",
        type=int,
        default=CLOSEUP["plane"],
        help=f"CCF plane at 10 um (default {CLOSEUP['plane']})",
    )
    parser.add_argument("--vmax", type=float, help="colour range of the means")
    parser.add_argument("--dlim", type=float, help="colour range of the difference")
    parser.add_argument(
        "--smooth",
        type=smoothing,
        help="sigma in 20 um voxels: one, or three comma separated",
    )
    parser.add_argument("--cmap", help="matplotlib colormap of the mean panels")
    parser.add_argument("--no-video", action="store_true", help="no video")
    parser.add_argument("--no-flatmap", action="store_true", help="no flatmaps")
    args = parser.parse_intermixed_args()
    main(
        args.readings,
        args.plane,
        args.vmax,
        args.dlim,
        args.smooth,
        want_video=not args.no_video,
        want_flatmap=not args.no_flatmap,
        cmap_name=args.cmap,
    )
