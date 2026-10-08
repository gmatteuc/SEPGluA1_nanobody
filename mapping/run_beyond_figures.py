"""The guided figures of analysis 4, with their intervals: figures 03, 04 and 11.

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
    21. run_synaptome          the measured synapse density (network,
                               once): PSD95 puncta per structure,
                               its coverage of the fit
    22. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    23. run_beyond_controls    seven attempts to break it
    24. run_beyond_calibration the same model on maps whose answer
                               is known
    25. run_beyond_regression  the regression, per structure
    26. run_beyond_figures     figures 03, 04 and 11                <- this script
    27. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 14
    28. run_ish_overview       figures 00 and 15, the figure index;
                               the numbers for the text

Jackknife intervals over structures, the noise band of the replication, the
numbers the text quotes, and the three figures; it reads the tables of steps 22 to
25. The method is in sepmap/adult/beyond_figures.py. Writes, under
adult_v2/ish_analysis/ in the data root:

    beyond/jackknife.csv                per subsample, the ceiling and the shares
    beyond/numbers_for_the_caption.txt  the figures' numbers as sentences
    tables/numbers_beyond.csv           the numbers of analysis 4, for the text
    figures/03_beyond_budget.png        what Gria1 and synapse density predict
    figures/04_beyond_where.png         where the leftover lives
    figures/11_leftover_genes.png       the genes and gene sets against the leftover

    python run_beyond_figures.py
"""

import argparse

import matplotlib

from sepmap import config, plotting
from sepmap.adult import beyond_figures


def main():
    """Print the settings in force, then run the intervals and figures 03, 04 and 11."""
    config.print_settings({})
    beyond_figures.main()


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="figures 03, 04 and 11 of analysis 4")
    parser.parse_args()
    main()
