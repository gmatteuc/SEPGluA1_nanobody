"""Young against adult in the adult CCF: maps and the per-structure table.

Python route, in run order (tools\\venv_atlas; run_closeup in tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid
     3. run_cohort             cohort mean, SD and n per voxel
     4. run_compare            young against adult: maps, table     <- this script
        run_replot             redraws those maps (not in a full run)
     5. run_region_plot        region statistics, no warping
     6. run_region_groups      the same by system and by layer
     7. run_video              cohort videos
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
                               the Cacng8 - Gria1 gap; figures 05 to 07 and 12
    18. run_ish_robustness     A3: the ranking under other choices;
                               figure 13
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figure 10, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation
                               against matched controls; figures 08, 09
    21. run_synaptome          the measured synapse density (network,
                               once): PSD95 puncta per structure,
                               its coverage of the fit
    22. run_beyond_density     analysis 4: what receptor mRNA and
                               synaptic density leave; the leftover
    23. run_beyond_controls    seven attempts to break it
    24. run_beyond_calibration the same model on maps whose answer
                               is known
    25. run_beyond_regression  the regression, per structure
    26. run_beyond_figures     figures 03, 04 and 11
    27. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 14
    28. run_ish_overview       figures 00 and 15, the figure index;
                               the numbers for the text

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
