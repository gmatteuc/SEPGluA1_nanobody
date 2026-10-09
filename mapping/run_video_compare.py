"""Young beside adult, plane by plane: a video per reading, or one plane as a still.

Python route, in run order (tools\\venv_atlas; run_closeup in tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid
     3. run_cohort             cohort mean, SD and n per voxel
     4. run_compare            young against adult: maps, table
        run_replot             redraws those maps (not in a full run)
     5. run_region_plot        region statistics, no warping
     6. run_region_groups      the same by system and by layer
     7. run_video              cohort videos
     8. run_video_compare      young beside adult, plane by plane   <- this script
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
                               documentation; figure 02
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the
                               map, the null, autofluorescence (A8),
                               the Cacng8 - Gria1 gap; figures 05 to 07 and 15
    18. run_ish_robustness     A3: the ranking under other choices;
                               figure 16
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figure 10, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation
                               against matched controls; figures 08, 09
    21. run_synaptome          the measured synapse density (network,
                               once): PSD95 puncta per structure,
                               its coverage of the fit
    22. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    23. run_beyond_controls    seven attempts to break it
    24. run_beyond_calibration the same model on maps whose answer
                               is known
    25. run_beyond_regression  the regression, per structure
    26. run_beyond_figures     figures 03, 04 and 11
    27. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 12 to 14, sheets
    28. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 17
    29. run_ish_overview       figures 00 and 18, the figure index;
                               the numbers for the text

Three panels per frame: the young and the adult cohort mean on one colour
scale, and their comparison; the layout is described in
sepmap/young_vs_adult/video_compare.py. Writes, in comparisons_v2/young_vs_adult/
under the data root:

    video_side_by_side_<reading>.mp4       a video per reading
    plane<P>_side_by_side_<reading>.png    with --plane, the still instead

with _vmax<V> before the extension when --vmax is given.

    python run_video_compare.py [reading ...] [--plane P] [--vmax V] [--dlim D]

With no reading named, every reading in force. --plane draws that CCF plane
(10 um numbering) as a PNG instead of the video; --vmax and --dlim set the
colour range of the two means and of the difference for this run only.
Readings and options can come in any order.

The readings in force, which a run with no reading named draws, are
readings.in_force of settings.toml; V2_READINGS, a comma-separated subset of them
(ratio, sepratio, cref, subref, zref), replaces it for one run. The run prints
them.
"""

import argparse

import matplotlib

from sepmap import config
from sepmap.volumes import cohort
from sepmap.young_vs_adult import video_compare


def main(readings, plane=None, vmax=None, dlim=None):
    """Print the settings in force, then draw young beside adult for `readings`.

    `readings` empty means every reading in force.
    """
    readings = readings or list(cohort.MODES)

    # settings in force
    in_force = cohort.readings_in_force()
    if plane is None:
        plane_text = "every plane, as a video"
    else:
        plane_text = f"{plane}, as a still"
    config.print_settings(
        {
            "readings": " ".join(readings),
            "in force": in_force,
            "plane": plane_text,
            "vmax": "per reading" if vmax is None else vmax,
            "dlim": "default" if dlim is None else dlim,
        }
    )

    # young beside adult
    video_compare.main(readings, plane=plane, vmax=vmax, dlim=dlim)


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")

    parser = argparse.ArgumentParser(description="young beside adult, plane by plane")
    parser.add_argument("readings", nargs="*", help="readings (default: all in force)")
    parser.add_argument("--plane", type=int, help="one CCF plane (10 um) as a still")
    parser.add_argument("--vmax", type=float, help="colour range of the two means")
    parser.add_argument("--dlim", type=float, help="colour range of the difference")
    args = parser.parse_intermixed_args()
    main(args.readings, plane=args.plane, vmax=args.vmax, dlim=args.dlim)
