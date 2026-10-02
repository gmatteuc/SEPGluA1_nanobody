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
    21. run_ish_reliability    how reliable one ISH map is            <- this script
    22. run_ish_panel_test     localisation against controls
    23. run_beyond_density     what abundance and density leave
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_figures     the figures of that result
    26. run_beyond_regression  the regression, shown

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
