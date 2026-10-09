"""Analysis 2 of the ISH line: between divisions or within them (A6).

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
    19. run_ish_divisions      analysis 2 (A6): between or within   <- this script
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
    30. run_ish_overview       figures 00, 15 and 15s, the figure
                               index; the numbers for the text

Splits each gene's rho with the adult map into what the contrast between divisions
gives (rho with a map that knows only each structure's division) and what is left
inside divisions (Spearman within each division, averaged), and tests the within
rho against the map's surrogates and against shuffles inside divisions; both nulls
are checked on random smooth maps (the method is in sepmap/ish/divisions.py).
Writes, in adult_v2/ish_analysis/ under the data root:

    tables/within_division.csv          per gene: whole-brain, division-only and
                                        within rho, the within p of both nulls, q
    tables/within_division_detail.csv   per gene and division, rho inside it
    tables/within_calibration.csv       random smooth maps tested inside divisions
    tables/numbers_divisions.csv        the numbers of this step, for the text
    figures/09_between_within.png       whether genes follow the map inside
                                        divisions
    figures/09s_between_within_detail.png  the same in detail
    figures/genes/<gene>.png            with --sheets, one sheet per detail gene

    python run_ish_divisions.py [--sheets]

--sheets also draws the per-gene sheets (a few seconds each).
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sepmap import config, plotting, structures
from sepmap.adult import profiles
from sepmap.ish import (
    divisions,
    gene_ranking,
    gene_sets,
    gene_table,
    planes,
    spatial_null,
)
from sepmap.ish.figure_index import figure_file, figure_path
from sepmap.ish.plotting import sets as set_figures
from sepmap.volumes.per_mouse import structure_terms

ISH_ANALYSIS = config.SETTINGS["ish_analysis"]
ISH_FIGURES = config.SETTINGS["ish_figures"]


def draw_figures(within, detail, map_values, division, merged, calibration, subunits):
    """Figure 09 and its detailed version."""
    set_table = structures.load_structure_set()
    q = ISH_ANALYSIS["q"]
    coarse = pd.Series(
        divisions.division_only(map_values.to_numpy(float), division),
        index=map_values.index,
    )
    example = divisions.EXAMPLE_GENE
    min_structures = ISH_ANALYSIS["min_division_structures"]
    fig = set_figures.plot_between_within(
        within,
        map_values,
        coarse,
        merged[example],
        example,
        set_table,
        subunits,
        divisions.DETAIL_GENES,
        q,
        min_structures,
        save=figure_path("between_within"),
    )
    plt.close(fig)
    fig = set_figures.plot_between_within_detail(
        within,
        detail,
        map_values,
        coarse,
        merged[example],
        example,
        set_table,
        calibration,
        subunits,
        divisions.DETAIL_GENES,
        q,
        min_structures,
        save=figure_path("between_within_detail"),
    )
    plt.close(fig)
    drawn = [figure_file(k) for k in ("between_within", "between_within_detail")]
    print(f"figures: {', '.join(drawn)} in {structures.FIGURES}")


def draw_sheets(within, detail, map_values, merged, ranking, n_surrogates):
    """One sheet per gene of divisions.DETAIL_GENES, in figures/genes/."""
    set_table = structures.load_structure_set()
    plane = ISH_FIGURES["plane"]
    names, _, _ = structure_terms()
    lab = planes.label_plane(plane)
    folder = structures.FIGURES / "genes"
    folder.mkdir(parents=True, exist_ok=True)
    by_gene = within.set_index("symbol")
    nano_rows = gene_ranking.map_rows(ranking, "nano")
    for gene in divisions.DETAIL_GENES:
        fig = set_figures.plot_gene_sheet(
            gene,
            map_values,
            merged[gene],
            set_table,
            nano_rows.loc[gene],
            by_gene.loc[gene],
            detail[detail["symbol"] == gene],
            lab,
            names,
            plane,
            n_surrogates,
            save=folder / f"{gene}.png",
        )
        plt.close(fig)
    print(f"gene sheets: {len(divisions.DETAIL_GENES)} in {folder}")


def main(sheets):
    """Print the settings, split every gene's rho; tables, figures, gene sheets."""
    config.print_settings({"sheets": sheets})
    q = ISH_ANALYSIS["q"]

    # the declared structures, their divisions, the map and its surrogates
    declared = structures.declared_structures()
    set_table = structures.load_structure_set()
    division_of = dict(zip(set_table["structure"], set_table["division"]))
    division = np.array([division_of[s] for s in declared])
    profile = profiles.load_profile()
    map_values = profile.loc[declared, "zref_nano"]
    surr, _ = spatial_null.load_surrogates("nano", declared)
    print(f"map: {len(declared)} structures in {len(set(division))} divisions")

    # the genes
    genes = gene_table.per_gene(gene_table.load_gene_table())
    p9_genes = set(genes.loc[genes["p9_gene"], "symbol"])
    subunits = set(gene_sets.members_from_table(genes)["subunits"])
    merged = gene_table.load_profiles()
    vectors = gene_ranking.gene_vectors(merged, declared)

    # between and within, per gene, with the two nulls of the within rho; the
    # whole-brain rho must be analysis 1's
    within, detail = divisions.gene_rows(
        map_values.to_numpy(float), division, surr, vectors
    )
    within = divisions.with_q(within, p9_genes)
    gene_ranking.check_matches_ranking(
        within.set_index("symbol")["rho"], "the whole-brain rho"
    )
    within.to_csv(divisions.WITHIN, index=False)
    detail.to_csv(divisions.WITHIN_DETAIL, index=False)
    print(
        f"within: {int(within['rho_within'].notna().sum())} genes; median "
        f"{within['rho_within'].median():+.3f} against {within['rho'].median():+.3f} "
        f"whole-brain; past the spatial null at q < {q}: "
        f"{int((within['q_all_spatial'] < q).sum())} of all"
    )

    # the two nulls on random smooth maps with the map's smoothness
    calibration = divisions.within_calibration(
        map_values.to_numpy(float),
        division,
        surr,
        structures.load_centroids().loc[declared],
    )
    calibration.to_csv(divisions.CALIBRATION, index=False)
    rates = ", ".join(
        f"{null} {float((calibration[f'p_{null}'] < spatial_null.ALPHA).mean()):.1%}"
        for null in ("shuffle", "spatial")
    )
    print(
        f"calibration: {len(calibration)} random maps, p < {spatial_null.ALPHA} by "
        f"{rates}"
    )

    # the numbers for the text
    numbers = divisions.numbers_table(within, detail, calibration)
    numbers.to_csv(divisions.NUMBERS, index=False)

    # figure 09 and its detailed version, and the gene sheets
    draw_figures(within, detail, map_values, division, merged, calibration, subunits)
    if sheets:
        ranking = gene_ranking.load_ranking()
        draw_sheets(within, detail, map_values, merged, ranking, surr.shape[0])


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    parser = argparse.ArgumentParser(
        description="analysis 2: between or within divisions"
    )
    parser.add_argument(
        "--sheets", action="store_true", help="also draw the per-gene sheets"
    )
    args = parser.parse_args()
    main(args.sheets)
