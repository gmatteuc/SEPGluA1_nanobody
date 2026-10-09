"""How reliable one Allen ISH map is, from the genes measured more than once.

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

Not in the run order: it writes adult_v2/ish/, the frozen run of 5 October, so a
rerun would overwrite it; run_ish_gene_table.py (step 15) replaces it, and it moves
to archive/ with the next step of the refactor.

Reads the region table of one panel pass (by default ontology,
gene_region_table_panel.csv); the method is in sepmap/ish/reliability.py. Writes,
in adult_v2/ish/ under the data root:

    gene_reliability.csv            per gene: experiments, planes, reliability,
                                    median energy
    gene_region_table_merged.csv    one profile per gene, its experiments averaged
    ish_reliability.png             the reliability, and what it depends on

    python run_ish_reliability.py [--panel ontology|targets]

V2_ISH_PANEL and V2_ISH_TABLE are refused: the pass is chosen with --panel.
"""

import argparse
import os

import matplotlib

from sepmap import config
from sepmap.ish import reliability

# environment variables that are refused: the pass is chosen with --panel
OLD_VARIABLES = ("V2_ISH_PANEL", "V2_ISH_TABLE")


def main(panel_name):
    """Print the settings in force, then measure the reliability of one pass."""
    # settings in force
    table = reliability.ISH_PANELS[panel_name]["table"]
    config.print_settings({"panel": panel_name, "table": table})

    # reliability and merged profiles
    reliability.main(panel_name=panel_name)


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")

    parser = argparse.ArgumentParser(description="reliability of one Allen ISH map")
    parser.add_argument(
        "--panel",
        choices=list(reliability.ISH_PANELS),
        default=reliability.DEFAULT_PANEL,
        help=f"panel pass (default {reliability.DEFAULT_PANEL})",
    )
    args = parser.parse_args()

    # a run that still sets the old variables would silently get the other panel
    old = [name for name in OLD_VARIABLES if os.environ.get(name)]
    if old:
        parser.error(
            f"{', '.join(old)}: the panel is no longer chosen through the "
            f"environment; unset {'it' if len(old) == 1 else 'them'} and use "
            f"--panel {'|'.join(reliability.ISH_PANELS)}"
        )
    main(args.panel)
