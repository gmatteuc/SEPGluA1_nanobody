"""The ISH line on one page, April's headline, the figure index and the numbers.

Python route, in run order (tools\\venv_atlas; run_closeup and run_adult_layers in
tools\\venv_flat):
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
    29. run_sep_channel_check  analysis 5: what the green channel
                               reports; figures 14, 14s
    30. run_ish_overview       figures 00, 15 and 15s, the figure   <- this script
                               index; the numbers for the text
    31. run_adult_layers       the adult map by depth (venv_flat)

Measures what is left of April's headline (the category violins and their ANOVA)
on today's inputs and against the map's surrogates, gathers the numbers every step
of the ISH line wrote into one table, and draws the overview and the index of the
figures (the methods are in sepmap/ish/april_headline.py, sepmap/ish/numbers.py and
sepmap/ish/overview.py). Reads only tables; run it last.
Writes, in adult_v2/ish_analysis/ under the data root:

    tables/april_headline.csv          P9's genes: group, April's rho, today's rho
    tables/april_anova.csv             the group ANOVA under each choice
    tables/april_groups.csv            today's groups against the surrogates
    tables/numbers_overview.csv        the numbers of this step
    tables/numbers_for_the_text.csv    every step's numbers, the one source of
                                       docs/ISH_ANALYSIS.md (and a .txt to read)
    figures/00_overview.png            the question, the argument, which figure
                                       answers what
    figures/15_april_headline.png      April's headline, then and now (15s in
                                       detail)
    figures/README.md                  the guided walk: each figure with its
                                       question, what to look at, what to take

    python run_ish_overview.py
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from sepmap import config, plotting, structures
from sepmap.ish import april_headline, gene_ranking, numbers, overview, robustness
from sepmap.ish.figure_index import FIGURE_NUMBERS, figure_file, figure_path
from sepmap.ish.plotting import headline as headline_figures

ISH_ANALYSIS = config.SETTINGS["ish_analysis"]


def main():
    """Print the settings, measure April's headline, gather the numbers, draw."""
    config.print_settings({})
    q = ISH_ANALYSIS["q"]

    # April's headline: P9's genes then and now, the ANOVA under each choice
    ranking = gene_ranking.load_ranking()
    robustness_rows = robustness.load_robustness()
    headline = april_headline.headline_table(april_headline.load_p9_headline(), ranking)
    anova = april_headline.anova_table(headline, robustness_rows)
    headline.to_csv(april_headline.HEADLINE, index=False)
    anova.to_csv(april_headline.HEADLINE_ANOVA, index=False)
    first = anova.iloc[0]
    print(
        f"April: {first['n_genes']} genes, group ANOVA p {first['p']:.4f}; "
        f"{len(anova)} choices, p {anova['p'].min():.4f} to {anova['p'].max():.4f}"
    )

    # the same groups today, against the surrogates of the map
    null, null_genes = gene_ranking.load_null_rho("nano")
    groups, f_info = april_headline.headline_null(headline, null, null_genes)
    groups.to_csv(april_headline.HEADLINE_GROUPS, index=False)
    print(
        f"against the surrogates: F {f_info['f']:.2f} over {f_info['n_genes']} genes, "
        f"spatial p {f_info['p_spatial']:.4f}; {null.shape[1]} surrogates"
    )

    # the numbers of this step, then every step's in one table
    per_adult = pd.read_csv(gene_ranking.PER_ADULT)
    summary = robustness.load_summary()
    mine = pd.concat(
        [
            april_headline.headline_numbers(headline, anova, groups, f_info),
            overview.leftover_extremes(),
            overview.ranking_extras(ranking, per_adult, summary, robustness_rows),
            overview.gene_extras(),
        ],
        ignore_index=True,
    )
    mine.to_csv(numbers.numbers_path("overview"), index=False)
    gathered = numbers.gather_numbers()
    gathered.to_csv(numbers.ALL_NUMBERS, index=False)
    numbers.ALL_NUMBERS_TXT.write_text(
        "\n".join(numbers.text_lines(gathered)) + "\n", encoding="utf-8"
    )
    n = numbers.lookup(gathered)
    print(f"numbers: {len(gathered)} from {gathered['step'].nunique()} steps")

    # figure 15 and its detailed version: April's headline; figure 00: the
    # overview; the index
    fig = headline_figures.plot_april_headline(
        headline,
        anova,
        groups,
        f_info,
        null.shape[1],
        q,
        save=figure_path("april_headline"),
    )
    plt.close(fig)
    fig = headline_figures.plot_april_headline_detail(
        headline,
        anova,
        groups,
        f_info,
        april_headline.HEADLINE_ORDER,
        null.shape[1],
        q,
        save=figure_path("april_headline_detail"),
    )
    plt.close(fig)
    fig = headline_figures.plot_overview(
        overview.overview_content(n),
        overview.figure_map(),
        save=figure_path("overview"),
    )
    plt.close(fig)
    overview.INDEX.write_text(overview.figure_index(n), encoding="utf-8")
    missing = [key for key in FIGURE_NUMBERS if not figure_path(key).exists()]
    keys = ("overview", "april_headline", "april_headline_detail")
    drawn = [figure_file(k) for k in keys]
    print(f"figures: {', '.join(drawn)} and the index, in {structures.FIGURES}")
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
