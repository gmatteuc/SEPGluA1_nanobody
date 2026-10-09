"""The ISH line on one page, April's headline, the figure index and the numbers.

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
    29. run_ish_overview       figures 00 and 18, the figure index;   <- this script
                               the numbers for the text

Measures what is left of April's headline (the category violins and their ANOVA)
on today's inputs and against the map's surrogates, gathers the numbers every step
of the ISH line wrote into one table, and draws the overview and the index of the
figures (the method is in sepmap/ish/overview.py). Reads only tables; run it last.
Writes, in adult_v2/ish_analysis/ under the data root:

    tables/april_headline.csv          P9's genes: group, April's rho, today's rho
    tables/april_anova.csv             the group ANOVA under each choice
    tables/april_groups.csv            today's groups against the surrogates
    tables/numbers_overview.csv        the numbers of this step
    tables/numbers_for_the_text.csv    every step's numbers, the one source of
                                       docs/ISH_ANALYSIS.md (and a .txt to read)
    figures/00_overview.png            the question, the argument, which figure
                                       answers what
    figures/18_april_headline.png      April's headline, then and now
    figures/README.md                  the guided walk: each figure with its
                                       question, what to look at, what to take

    python run_ish_overview.py
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from sepmap import config, plotting
from sepmap.ish import gene_ranking, overview
from sepmap.ish import plotting as ish_plotting
from sepmap.structures import TABLES

OUT = config.DATA / "adult_v2" / "ish_analysis"


def main():
    """Print the settings, measure April's headline, gather the numbers, draw."""
    config.print_settings({})
    figures = OUT / "figures"

    # April's headline: P9's genes then and now, the ANOVA under each choice
    ranking = gene_ranking.load_ranking()
    headline = overview.headline_table(overview.load_p9_headline(), ranking)
    robustness_rows = pd.read_csv(TABLES / "ranking_robustness.csv")
    anova = overview.anova_table(headline, robustness_rows)
    headline.to_csv(overview.HEADLINE, index=False)
    anova.to_csv(overview.HEADLINE_ANOVA, index=False)
    first = anova.iloc[0]
    print(
        f"April: {first['n_genes']} genes, group ANOVA p {first['p']:.4f}; "
        f"{len(anova)} choices, p {anova['p'].min():.4f} to {anova['p'].max():.4f}"
    )

    # the same groups today, against the surrogates of the map
    null, null_genes = gene_ranking.load_null_rho("nano")
    groups, f_info = overview.headline_null(headline, null, null_genes)
    groups.to_csv(overview.HEADLINE_GROUPS, index=False)
    print(
        f"against the surrogates: F {f_info['f']:.2f} over {f_info['n_genes']} genes, "
        f"spatial p {f_info['p_spatial']:.4f}; {null.shape[1]} surrogates"
    )

    # the numbers of this step, then every step's in one table
    per_adult = pd.read_csv(gene_ranking.PER_ADULT)
    summary = pd.read_csv(TABLES / "robustness_summary.csv")
    mine = pd.concat(
        [
            overview.headline_numbers(headline, anova, groups, f_info),
            overview.leftover_extremes(),
            overview.ranking_extras(ranking, per_adult, summary, robustness_rows),
            overview.gene_extras(),
        ],
        ignore_index=True,
    )
    mine.to_csv(TABLES / "numbers_overview.csv", index=False)
    numbers = overview.gather_numbers()
    numbers.to_csv(overview.NUMBERS, index=False)
    with open(overview.NUMBERS_TXT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(overview.text_lines(numbers)) + "\n")
    n = overview.lookup(numbers)
    print(f"numbers: {len(numbers)} from {numbers['step'].nunique()} steps")

    # figure 18: April's headline; figure 00: the overview; the index
    fig = ish_plotting.plot_april_headline(
        headline,
        anova,
        groups,
        f_info,
        overview.HEADLINE_ORDER,
        null.shape[1],
        overview.ISH_ANALYSIS["q"],
        save=figures / ish_plotting.figure_file("april_headline"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_overview(
        overview.argument_text(n),
        overview.overview_rows(n),
        overview.ITEMS,
        save=figures / ish_plotting.figure_file("overview"),
    )
    plt.close(fig)
    with open(overview.INDEX, "w", encoding="utf-8") as fh:
        fh.write(overview.figure_index(n))
    missing = [
        key
        for key in ish_plotting.FIGURES
        if not (figures / ish_plotting.figure_file(key)).exists()
    ]
    drawn = [ish_plotting.figure_file(k) for k in ("overview", "april_headline")]
    print(f"figures: {', '.join(drawn)} and the index, in {figures}")
    if missing:
        print(f"  warning: the index names figures not drawn yet: {', '.join(missing)}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(
        description="the ISH line on one page, April's headline, the numbers"
    )
    parser.parse_args()
    main()
