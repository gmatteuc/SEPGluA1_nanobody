"""Diagnostic sheets for every step of the route, so it can be audited by eye.

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
     9. run_closeup            close-ups and flatmaps (venv_flat)
    10. run_diagnostics        sheets that audit each step          <- this script
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

Per brain the tissue and level sheets (and the warp sheet of a young brain),
then the cohort-level sheets and the README index; the sheets are listed in
sepmap/diagnostics.py. Sheet 08 reads the background masks of
run_normalise_groups (MATLAB), so in a full run this comes after the plasticity
chain too. Writes, in comparisons_v2/processing_diagnostics/ under the data root:

    01_tissue_<mouse>.png    per brain: what was counted as tissue
    02_levels_<mouse>.png    per brain: where the background and threshold sit
    04_warp_<mouse>.png      per young brain: before and after DeMBA -> CCF
    03_, 05_ to 09_*.png     the cohort-level sheets
    README.md                the index: what to look for in each sheet

    python run_diagnostics.py [mouse ...]

With mouse names, only those brains' sheets are refreshed; with none, every
brain's sheets, the cohort-level sheets and the index.
"""

import argparse

import matplotlib

from sepmap import config, diagnostics


def main(mice):
    """Print the settings in force, then draw the diagnostic sheets.

    Only the sheets of `mice`, or every sheet when the list is empty.
    """
    # settings in force
    which = " ".join(mice) if mice else "every brain, and the cohort-level sheets"
    config.print_settings({"mice": which})

    # diagnostic sheets
    diagnostics.main(mice)


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")

    parser = argparse.ArgumentParser(description="diagnostic sheets of the route")
    parser.add_argument(
        "mice", nargs="*", help="brains whose sheets to refresh (default: every sheet)"
    )
    args = parser.parse_args()
    main(args.mice)
