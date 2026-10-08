"""Allen ISH expression per gene and structure, for one gene panel.

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
                               the Cacng8 - Gria1 gap; figures 05 to 07 and 12
    18. run_ish_robustness     A3: the ranking under other choices;
                               figure 13
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
    25. run_beyond_figures     figures 03, 04 and 11
    26. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 14
    27. run_ish_overview       figures 00 and 15, the figure index;
                               the numbers for the text

Not in the run order: it writes adult_v2/ish/, the frozen run of 5 October, so a
rerun would overwrite it; run_ish_gene_table.py (step 15) replaces it, and it moves
to archive/ with the next step of the refactor.

One mean per gene and structure, over the full structure and eroded by one
200 um voxel, keyed by structure name as the nano tables are (the method is in
sepmap/ish/regions.py). The panel, and the table written in adult_v2/ish/ under
the data root, are one of the passes of settings.toml ([ish_panels]):

    targets     gene_region_table.csv          the hand-written 100-gene panel
    ontology    gene_region_table_panel.csv    the 390-gene panel of run_panel_build

and beside the table, <table>_drops.csv lists the experiments dropped, and why.

    python run_ish_regions.py [gene ...] [--panel targets|ontology]

With gene symbols, only those genes of the panel. --panel defaults to targets.
V2_ISH_PANEL and V2_ISH_TABLE are refused: the pass is chosen with --panel.
"""

import argparse
import os

from sepmap import config
from sepmap.ish import regions

# environment variables that are refused: the pass is chosen with --panel
OLD_VARIABLES = ("V2_ISH_PANEL", "V2_ISH_TABLE")


def main(genes, panel_name):
    """Print the settings in force, then write the region table of one panel pass.

    With `genes`, only those genes of the panel.
    """
    # settings in force, with the panel file and the table of the pass
    panel, table = regions.panel_files(panel_name)
    config.print_settings(
        {
            "panel": f"{panel_name}  ({panel})",
            "table": table,
            "genes": " ".join(genes) if genes else "the whole panel",
        }
    )

    # one mean per gene and structure
    regions.main(genes or None, panel_name=panel_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ISH expression per gene and structure")
    parser.add_argument("genes", nargs="*", help="genes (default: the whole panel)")
    parser.add_argument(
        "--panel",
        choices=list(regions.ISH_PANELS),
        default=regions.DEFAULT_PANEL,
        help=f"panel pass (default {regions.DEFAULT_PANEL})",
    )
    args = parser.parse_intermixed_args()

    # a run that still sets the old variables would silently get the other panel
    old = [name for name in OLD_VARIABLES if os.environ.get(name)]
    if old:
        parser.error(
            f"{', '.join(old)}: the panel is no longer chosen through the "
            f"environment; unset {'it' if len(old) == 1 else 'them'} and use "
            f"--panel {'|'.join(regions.ISH_PANELS)}"
        )
    main(args.genes, args.panel)
