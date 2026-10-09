"""Analysis 1 of the ISH line: each gene against the adult map, with the spatial null.

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
    17. run_ish_gene_ranking   analysis 1: each gene against the    <- this script
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
                               its coverage of the fit; figure 14s1
    22. run_density_markers    the synapse-density genes, chosen by
                               PSD95 without the map (network,
                               once); figure 14s2
    23. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_calibration the same model on maps whose answer
                               is known
    26. run_beyond_regression  the regression, per structure
    27. run_beyond_figures     figures 03, 03s1, 03s2, 04, 11s1, 14
    28. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 07, 08, 11, 11s2,
                               sheets
    29. run_sep_channel_check  analysis 5: what the green channel
                               reports; figures 15, 15s
    30. run_ish_overview       figures 00, 16 and 16s, the figure
                               index; the numbers for the text

Correlates every gene of the gene table with the adult nano map and with the
autofluorescence map of the same sections on the declared structures, tests each
rho against the map's surrogates, corrects within P9's genes and within all genes,
measures how much each rho depends on which adults were measured, and tests the
Cacng8 - Gria1 gap (the method is in sepmap/ish/gene_ranking.py). Writes, in
adult_v2/ish_analysis/ under the data root:

    tables/gene_ranking.csv            per map and gene: rho, structures, spatial p,
                                       q within P9's genes and within all, null band,
                                       spread over resampled adults, ranks, labels
    tables/gene_ranking_per_adult.csv  per adult and gene, rho with that adult's
                                       nano and autofluorescence maps
    tables/gap.csv                     the Cacng8 - Gria1 gap, merged and for each
                                       pairing of their experiments
    tables/null_rho.npz                every gene's rho with every surrogate, per map,
                                       and the gap of every surrogate
    tables/numbers_gene_ranking.csv    the numbers of this step, for the text
    figures/05_one_comparison.png      what one comparison is (05s in detail)
    figures/06_spatial_null.png        why a null, and that it works (06s in
                                       detail)
    figures/07s_gene_ranking.png       P9's genes against the map and its null,
                                       the Cacng8 - Gria1 gap
    figures/12_autofluorescence.png    the label or the tissue (12s in detail)

    python run_ish_gene_ranking.py
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from sepmap import config, plotting, structures
from sepmap.adult import profiles
from sepmap.ish import gene_ranking, gene_sets, gene_table, planes, spatial_null
from sepmap.ish import plotting as ish_plotting
from sepmap.ish.figure_index import figure_file, figure_path
from sepmap.volumes.per_mouse import structure_terms

ISH_ANALYSIS = config.SETTINGS["ish_analysis"]
ISH_FIGURES = config.SETTINGS["ish_figures"]

# surrogates of the nano map drawn on the plane of figure 06s
SURROGATES_SHOWN = 3


def draw_comparison(ranking, merged, table, profile, surr, nulls, vectors):
    """Figures 05 and 06 with their detailed versions: one comparison, the null."""
    declared = structures.declared_structures()
    set_table = structures.load_structure_set()
    n_surrogates = surr["nano"].shape[0]
    plane = ISH_FIGURES["plane"]
    names, _, _ = structure_terms()
    lab = planes.label_plane(plane)
    map_values = profile.loc[declared, "zref_nano"]
    nano = gene_ranking.map_rows(ranking, "nano")
    rows = planes.comparison_rows(
        nano, merged, table, lab, plane, gene_ranking.QUOTED_GENES
    )
    fig = ish_plotting.plot_one_comparison(
        map_values,
        rows,
        set_table,
        lab,
        names,
        plane,
        n_surrogates,
        save=figure_path("one_comparison"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_one_comparison_detail(
        map_values,
        planes.nano_plane(plane),
        rows,
        set_table,
        lab,
        names,
        plane,
        n_surrogates,
        save=figure_path("one_comparison_detail"),
    )
    plt.close(fig)

    # the spatial null, with surrogates on the plane and the two genes against it;
    # a surrogate holds the ranks 1 to n, drawn as 0 to 1
    variogram = pd.read_csv(spatial_null.VARIOGRAM)
    calibration = pd.read_csv(spatial_null.CALIBRATION)
    fig = ish_plotting.plot_spatial_null(
        variogram, calibration, n_surrogates, save=figure_path("spatial_null")
    )
    plt.close(fig)
    n = len(declared)
    map_ranks = dict(zip(declared, ish_plotting.ranks01(map_values.to_numpy())))
    surrogate_ranks = [
        dict(zip(declared, (surr["nano"][k] - 1) / (n - 1)))
        for k in range(SURROGATES_SHOWN)
    ]
    tested = list(vectors)
    shown = [
        (g, nano.loc[g, "rho"], nulls["nano"][tested.index(g)], nano.loc[g, "p_spatial"])
        for g in gene_ranking.GAP_GENES
    ]
    fig = ish_plotting.plot_spatial_null_detail(
        variogram,
        calibration,
        map_ranks,
        surrogate_ranks,
        lab,
        names,
        plane,
        n_surrogates,
        genes=shown,
        save=figure_path("spatial_null_detail"),
    )
    plt.close(fig)


def draw_ranking(ranking, gap, gap_null, per_adult, subunits, n_surrogates):
    """Figure 07s, P9's genes and the gap; figure 12 and 12s, autofluorescence."""
    q = ISH_ANALYSIS["q"]
    fig = ish_plotting.plot_gene_ranking(
        ranking,
        gap,
        gap_null,
        subunits,
        ISH_FIGURES["t_max"],
        q,
        n_surrogates,
        save=figure_path("gene_ranking"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_autofluorescence(
        ranking,
        per_adult,
        subunits,
        q,
        n_surrogates,
        save=figure_path("autofluorescence"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_autofluorescence_detail(
        ranking,
        per_adult,
        subunits,
        q,
        n_surrogates,
        save=figure_path("autofluorescence_detail"),
    )
    plt.close(fig)


def main():
    """Print the settings, then rank every gene against both maps; tables, figures."""
    config.print_settings({})
    q = ISH_ANALYSIS["q"]

    # the declared structures, the two maps per adult and as cohort means, surrogates
    declared = structures.declared_structures()
    profile = profiles.load_profile()
    per_mouse = profiles.load_per_mouse()
    surr = {
        name: spatial_null.load_surrogates(name, declared)[0]
        for name in gene_ranking.MAPS
    }
    n_surrogates = surr["nano"].shape[0]
    print(f"maps: nano and auto on {len(declared)} structures, {n_surrogates} surrogates")

    # the genes: merged profiles and labels
    table = gene_table.load_gene_table()
    genes = gene_table.per_gene(table)
    merged = gene_table.load_profiles()
    vectors = gene_ranking.gene_vectors(merged, declared)
    subunits = set(gene_sets.members_from_table(genes)["subunits"])
    print(f"genes: {len(vectors)} of {len(genes)} share enough declared structures")

    # every gene against each map: rho, spatial p, null band, spread over adults
    ranking, nulls, adults = gene_ranking.rank_against_maps(
        per_mouse, profile, declared, surr, vectors, genes
    )
    ranking.to_csv(gene_ranking.RANKING, index=False)
    counts = ", ".join(
        f"{name} {len(gene_ranking.passing(ranking, name, 'p9'))} of P9's and "
        f"{len(gene_ranking.passing(ranking, name, 'all'))} of all"
        for name in gene_ranking.MAPS
    )
    print(f"ranking: past the null at q < {q}: {counts}")

    # each adult's own maps against every gene
    per_adult = gene_ranking.per_adult_table(
        adults, vectors, profiles.ADULTS, profiles.GROUP
    )
    per_adult.to_csv(gene_ranking.PER_ADULT, index=False)
    print(f"per adult: {per_adult['mouse'].nunique()} adults x {len(vectors)} genes")

    # the Cacng8 - Gria1 gap, merged and per pairing of experiments
    gap, gap_null = gene_ranking.gap_with_nulls(
        table, merged, profile, declared, surr, adults["nano"]
    )
    gap.to_csv(gene_ranking.GAP, index=False)
    gene_ranking.save_null_rho(list(vectors), nulls, gap_null)
    first = gene_ranking.merged_gap(gap)
    print(
        f"gap: {first['gap']:+.3f} on {first['n_structures']} structures, p "
        f"{first['p_equal']:.4f} against maps related to both alike (c "
        f"{first['equal_weight']:.3f}), {first['p_spatial']:.4f} against unrelated "
        f"maps; {len(gap) - 1} pairings of experiments"
    )

    # the numbers for the text
    numbers = gene_ranking.numbers_table(ranking, gap, per_adult)
    numbers.to_csv(gene_ranking.NUMBERS, index=False)

    # figures 05, 06, 07s and 12, with their detailed versions
    draw_comparison(ranking, merged, table, profile, surr, nulls, vectors)
    draw_ranking(ranking, gap, gap_null, per_adult, subunits, n_surrogates)
    keys = (
        "one_comparison",
        "one_comparison_detail",
        "spatial_null",
        "spatial_null_detail",
        "gene_ranking",
        "autofluorescence",
        "autofluorescence_detail",
    )
    drawn = [figure_file(k) for k in keys]
    print(f"figures: {', '.join(drawn)} in {structures.FIGURES}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="analysis 1: each gene against the map")
    parser.parse_args()
    main()
