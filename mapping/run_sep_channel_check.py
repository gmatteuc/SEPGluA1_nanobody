"""Analysis 5 of the ISH line: what the green channel reports; figure 15.

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
                               its coverage of the fit; figure 14s
    22. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    23. run_beyond_controls    seven attempts to break it
    24. run_beyond_calibration the same model on maps whose answer
                               is known
    25. run_beyond_regression  the regression, per structure
    26. run_beyond_figures     figures 03, 03s, 04, 11s1 and 14
    27. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 07, 08, 11, 11s2,
                               sheets
    28. run_sep_channel_check  analysis 5: what the green channel   <- this script
                               reports; figures 15, 15s
    29. run_ish_overview       figures 00, 16 and 16s, the figure
                               index; the numbers for the text

Per adult, on the declared structures: the dynamic range of each raw channel, what
each follows across structures (the others, and Gria1 mRNA), and what is left of
the green channel once its autofluorescence part is regressed out; no ratio of
channels is taken. The method is in sepmap/adult/sep_channel_check.py. Writes, under
adult_v2/ish_analysis/ in the data root:

    green_channel/sep_channel_check.csv   per adult: the ranges and correlations
    green_channel/sep_channel_check.png   the working figure
    tables/numbers_green_channel.csv      the numbers of analysis 5, for the text
    figures/15_green_channel.png          the guided figure (15s in detail)

    python run_sep_channel_check.py
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from sepmap import config, plotting
from sepmap.adult import sep_channel_check
from sepmap.ish import planes
from sepmap.ish import plotting as ish_plotting
from sepmap.structures import TABLES

ISH_FIGURES = config.SETTINGS["ish_figures"]
ISH = config.SETTINGS["ish"]

FIGURES = config.DATA / "adult_v2" / "ish_analysis" / "figures"


def numbers_table(rows: pd.DataFrame, n_structures: int) -> pd.DataFrame:
    """The numbers of analysis 5 that the text quotes, one row each."""
    out = [
        ("adults", len(rows), "adults measured"),
        ("structures", n_structures, "declared structures"),
    ]
    for column in rows.columns:
        if column in ("mouse", "n_structures"):
            continue
        v = rows[column]
        out += [
            (
                f"{column}_mean",
                round(float(v.mean()), 3),
                f"{column}, mean of the adults",
            ),
            (f"{column}_min", round(float(v.min()), 3), f"{column}, lowest adult"),
            (f"{column}_max", round(float(v.max()), 3), f"{column}, highest adult"),
        ]
    return pd.DataFrame(out, columns=["name", "value", "what"], dtype=object)


def main():
    """Print the settings, measure the channels per adult, draw figure 15."""
    config.print_settings({})
    per, rows = sep_channel_check.main()
    rows = pd.DataFrame(rows)
    n_structures = int(rows["n_structures"].iloc[0])
    numbers_table(rows, n_structures).to_csv(
        TABLES / "numbers_green_channel.csv", index=False
    )

    # figure 15, and 15s with the first adult's raw channels on the plane of the
    # guided figures
    fig = ish_plotting.plot_green_channel(
        rows, n_structures, save=FIGURES / ish_plotting.figure_file("green_channel")
    )
    plt.close(fig)
    plane = ISH_FIGURES["plane"]
    mouse = sep_channel_check.ADULTS[0]
    fig = ish_plotting.plot_green_channel_detail(
        rows,
        planes.channel_planes(mouse, plane),
        planes.label_plane(plane),
        mouse,
        plane,
        n_structures,
        ISH["control_gene"],
        save=FIGURES / ish_plotting.figure_file("green_channel_detail"),
    )
    plt.close(fig)
    keys = ("green_channel", "green_channel_detail")
    drawn = [ish_plotting.figure_file(k) for k in keys]
    print(f"figures: {', '.join(drawn)} in {FIGURES}")


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
