"""Cohort videos: the mean and its reliability t, plane by plane, per reading.

Python route, in run order (tools\\venv_atlas; run_closeup in tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid
     3. run_cohort             cohort mean, SD and n per voxel
     4. run_compare            young against adult: maps, table
        run_replot             redraws those maps (not in a full run)
     5. run_region_plot        region statistics, no warping
     6. run_region_groups      the same by system and by layer
     7. run_video              cohort videos                        <- this script
     8. run_video_compare      young beside adult, plane by plane
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

Reads the cohort volumes run_cohort.py writes, the t from the _folded ones (each
brain's hemispheres averaged first); the layout and the colour ranges are
described in sepmap/young_vs_adult/video.py. Writes, in
comparisons_v2/ccf/<cohort>/ under the data root:

    video_<reading>_<cohort>.mp4    the mean and its reliability t, plane by plane

    python run_video.py [cohort ...]

With no cohort named: young, adult, young_P20, naive and rws.

The readings are readings.in_force of settings.toml; V2_READINGS, a
comma-separated subset of them (ratio, sepratio, cref, subref, zref), replaces it
for one run and leaves the others out of this step. The run prints the readings in
force.
"""

import argparse

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
    readings = cohort.readings_in_force()
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
