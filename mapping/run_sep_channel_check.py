"""Analysis 5 of the ISH line: what the green channel reports; figure 14.

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
                               documentation; figures 02, 02s
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the
                               map, the null, autofluorescence (A8),
                               the Cacng8 - Gria1 gap; figures 05,
                               06, 07s and 12, with their s
    18. run_ish_robustness     A3: the ranking under other choices;
                               figures 13, 13s
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figures 09, 09s, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation
                               against matched controls; figures 10,
                               10s1, 10s2
    21. run_synaptome          the measured synapse density (network,
                               once): PSD95 puncta per structure,
                               its coverage of the fit; figure 03s4
    22. run_density_markers    the synapse-density genes, chosen by
                               PSD95 without the map (network,
                               once); figure 03s3
    23. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_calibration the same model on maps whose answer
                               is known; the check rows, each
                               with its own floor
    26. run_beyond_regression  the regression, per structure
    27. run_beyond_figures     figures 03, 03s1, 03s2, 04, 11s1
    28. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 07, 08, 11, 11s2,
                               sheets
    29. run_sep_channel_check  analysis 5: what the green channel   <- this script
                               reports; figures 14, 14s
    30. run_ish_overview       figures 00, 15 and 15s, the figure
                               index; the numbers for the text

Per adult, on the declared structures: the dynamic range of each raw channel, what
each follows across structures (the others, and Gria1 mRNA), and what is left of
the green channel once its autofluorescence part is regressed out; no ratio of
channels is taken. The method is in sepmap/adult/sep_channel_check.py. Writes, under
adult_v2/ish_analysis/ in the data root:

    green_channel/sep_channel_check.csv   per adult: the ranges and correlations
    green_channel/sep_channel_check.png   the working figure: the ranges, one adult,
                                          the Gria1 correlations
    tables/numbers_green_channel.csv      the numbers of analysis 5, for the text
    figures/14_green_channel.png          the guided figure (14s in detail)

    python run_sep_channel_check.py
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from sepmap import config, plotting, structures
from sepmap.adult import plotting as adult_plotting
from sepmap.adult import sep_channel_check
from sepmap.adult.profiles import ADULTS
from sepmap.ish import planes
from sepmap.ish.figure_index import figure_file, figure_path

ISH_FIGURES = config.SETTINGS["ish_figures"]
ISH = config.SETTINGS["ish"]


def main():
    """Print the settings, measure the channels per adult, draw figure 14."""
    config.print_settings({})
    per, rows = sep_channel_check.main()
    rows = pd.DataFrame(rows)
    n_structures = int(rows["n_structures"].iloc[0])
    numbers = sep_channel_check.numbers_table(rows, n_structures)
    numbers.to_csv(sep_channel_check.NUMBERS, index=False)

    # the working figure: the ranges, the first adult, the Gria1 correlations
    mouse = ADULTS[0]
    working = sep_channel_check.OUT / "sep_channel_check.png"
    fig = adult_plotting.plot_channel_check(
        rows, per[mouse], mouse, ISH["control_gene"], save=working
    )
    plt.close(fig)
    print(f"\nworking figure: {working}")

    # figure 14, and 14s with the first adult's raw channels on the plane of the
    # guided figures
    fig = adult_plotting.plot_green_channel(
        rows, n_structures, save=figure_path("green_channel")
    )
    plt.close(fig)
    plane = ISH_FIGURES["plane"]
    fig = adult_plotting.plot_green_channel_detail(
        rows,
        planes.channel_planes(mouse, plane),
        planes.label_plane(plane),
        mouse,
        plane,
        n_structures,
        ISH["control_gene"],
        save=figure_path("green_channel_detail"),
    )
    plt.close(fig)
    drawn = [figure_file(k) for k in ("green_channel", "green_channel_detail")]
    print(f"figures: {', '.join(drawn)} in {structures.FIGURES}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(
        description="analysis 5: what the green channel reports"
    )
    parser.parse_args()
    main()
