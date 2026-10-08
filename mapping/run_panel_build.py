"""The ontology-defined ISH gene panel: genes from GO terms, experiments from Allen.

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
    11. run_panel_build        390-gene ontology panel (network, cached)   <- this script
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

Builds the gene panel from Gene Ontology terms (the terms, and the one hand-made
call, are in sepmap/ish/panel_build.py), asking mygene.info and the Allen API and
caching every answer, so a re-run asks nothing twice. Writes, in adult_v2/panel/
under the data root:

    panel_v2.csv       one row per experiment: gene, role, id, plane, GO terms
    panel_genes.csv    one row per gene, with every term that claimed it
    cache/*.json       every API answer

    python run_panel_build.py
"""

import argparse

from sepmap import config
from sepmap.ish import panel_build


def main():
    """Print the settings in force, then build the panel."""
    # settings in force
    config.print_settings({})

    # the panel, from the ontology
    panel_build.main()


if __name__ == "__main__":
    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="the ontology-defined ISH gene panel")
    parser.parse_args()
    main()
