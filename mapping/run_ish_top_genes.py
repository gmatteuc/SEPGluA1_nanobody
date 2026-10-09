"""The genes that follow the map most, characterised; the tests named for the leftover.

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
                               its coverage of the fit; figure 14s1
    22. run_density_markers    the synapse-density genes, chosen by
                               PSD95 without the map (network,
                               once); figure 14s2
    23. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_calibration the same model on maps whose answer
                               is known; the check rows, each
                               with its own floor
    26. run_beyond_regression  the regression, per structure
    27. run_beyond_figures     figures 03, 03s1, 03s2, 04, 11s1, 14
    28. run_ish_top_genes      the genes that follow the map most,   <- this script
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 07, 08, 11, 11s2,
                               sheets
    29. run_sep_channel_check  analysis 5: what the green channel
                               reports; figures 15, 15s
    30. run_ish_overview       figures 00, 16 and 16s, the figure
                               index; the numbers for the text

Describes the genes past the spatial null after BH in analysis 1, Gria1, Cacng8 and
every member of the AMPA receptor complex family on every axis of the ISH line, adds
each to the main model of analysis 4 against maps of its own smoothness, and runs the
tests named in advance for the leftover: Cacng8 alone, then the family as a group
and gene by gene, with two check rows of the group test and the tests again against
smoother nulls (the method is in sepmap/ish/top_genes.py). It reads the tables of
steps 15 to 23. Writes, in adult_v2/ish_analysis/ under the data root:

    tables/top_genes.csv              one row per gene characterised
    tables/named_tests.csv            the tests named for the leftover, their check
                                      rows, and the same read on the map
    tables/family_members.csv         the family: each member's sources, tested or
                                      why not
    tables/leftover_null_check.csv    the leftover's tests against smoother nulls
    tables/numbers_top_genes.csv      the numbers of this step, for the text
    figures/07_top_genes.png          the genes that follow the map
    figures/08_cacng8_gria1.png       Cacng8 against Gria1
    figures/11_leftover.png           the three tiers against the leftover
    figures/11s2_ampa_family.png      the family against the leftover in detail
    figures/top_genes/<gene>.png      with --sheets, one sheet per gene characterised

    python run_ish_top_genes.py [--sheets]

--sheets also draws the per-gene sheets (about a second each).
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sepmap import config, plotting, structures
from sepmap.adult import beyond_density, profiles
from sepmap.ish import divisions, gene_ranking, gene_sets, planes, robustness, top_genes
from sepmap.ish import plotting as ish_plotting
from sepmap.ish.figure_index import figure_file, figure_path
from sepmap.volumes.per_mouse import structure_terms

ISH_ANALYSIS = config.SETTINGS["ish_analysis"]
ISH_FIGURES = config.SETTINGS["ish_figures"]


def added_shares(inputs, symbols):
    """What each gene takes of the leftover, against maps of its smoothness."""
    centroids = structures.load_centroids()
    rows, nulls = [], {}
    for symbol in symbols:
        row, nulls[symbol] = top_genes.added_share(inputs, symbol, centroids)
        rows.append(dict(symbol=symbol, **row))
        print(
            f"  {symbol:8s} on {row['n_added']} structures takes "
            f"{row['taken']:+.3f} of the reproducible map; 95% of maps of its "
            f"smoothness take less than {row['taken_null_hi']:+.3f} (p "
            f"{row['p_taken']:.3f}), of maps alike to the model "
            f"{row['taken_null_alike_hi']:+.3f} (p {row['p_taken_alike']:.3f})",
            flush=True,
        )
    return pd.DataFrame(rows), nulls


def null_check(inputs, residual, null, leftover, family):
    """The leftover's tests against nulls smoother than its own surrogates."""
    check = top_genes.smooth_null_check(
        inputs,
        residual,
        null["surrogates"],
        leftover,
        family,
        structures.load_centroids(),
    )
    for r in check.itertuples():
        where = r.null if np.isnan(r.range_mm) else f"{r.null}, {r.range_mm:g} mm"
        line = f"null check, {where}: smoothness {r.neighbour_rho:.3f}"
        if np.isfinite(r.p_cacng8):
            line += (
                f"; Cacng8 p {r.p_cacng8:.4f}, the family's p {r.p_family:.3f}, "
                f"{int(r.genes_p05)} genes below p 0.05"
            )
        print(line)
    return check


def draw_figures(table, inputs, leftover, tests, group_nulls, members, residual, n):
    """Figures 07 and 08, and 11 and 11s2 on the main model's leftover.

    `table` is top_genes.csv, `residual` the leftover by structure and `n` holds the
    leftover's and the map's numbers of surrogates.
    """
    q = ISH_ANALYSIS["q"]
    set_table = structures.load_structure_set()
    map_values = profiles.load_profile().loc[
        structures.declared_structures(), "zref_nano"
    ]
    fig = ish_plotting.plot_top_genes(
        table, q, ISH_FIGURES["t_max"], save=figure_path("top_genes")
    )
    plt.close(fig)
    fig = ish_plotting.plot_cacng8_gria1(
        table.set_index("symbol"),
        map_values,
        inputs.expr,
        gene_ranking.merged_gap(gene_ranking.load_gap()),
        gene_ranking.load_gap_null()[0],
        set_table,
        n["map"],
        save=figure_path("cacng8_gria1"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_leftover(
        table,
        residual,
        inputs.expr[gene_sets.LEFTOVER_GENE],
        tests,
        group_nulls,
        leftover,
        set_table,
        q,
        n["leftover"],
        save=figure_path("leftover"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_ampa_family(
        table[table["family"]],
        tests,
        group_nulls,
        members,
        q,
        ISH_FIGURES["t_max"],
        n["leftover"],
        n["map"],
        save=figure_path("ampa_family"),
    )
    plt.close(fig)
    keys = ("top_genes", "cacng8_gria1", "leftover", "ampa_family")
    drawn = [figure_file(k) for k in keys]
    print(f"figures: {', '.join(drawn)} in {structures.FIGURES}")


def draw_sheets(table, inputs, residual, added_nulls, n):
    """One sheet per gene characterised, in figures/top_genes/."""
    plane = ISH_FIGURES["plane"]
    structure_names, _, _ = structure_terms()
    lab = planes.label_plane(plane)
    set_table = structures.load_structure_set()
    map_values = profiles.load_profile().loc[
        structures.declared_structures(), "zref_nano"
    ]
    folder = structures.FIGURES / "top_genes"
    folder.mkdir(parents=True, exist_ok=True)
    for _, row in table.iterrows():
        fig = ish_plotting.plot_top_gene_sheet(
            row,
            map_values,
            inputs.expr[row["symbol"]],
            residual,
            added_nulls[row["symbol"]],
            set_table,
            lab,
            structure_names,
            plane,
            n["leftover"],
            n["map"],
            save=folder / f"{row['symbol']}.png",
        )
        plt.close(fig)
    print(f"gene sheets: {len(table)} in {folder}")


def main(sheets):
    """Print the settings; characterise the genes, run the named tests, draw."""
    config.print_settings({"sheets": sheets})
    q = ISH_ANALYSIS["q"]

    # the tables of analyses 1 to 4, and the main model with its leftover's null
    nano = gene_ranking.load_ranking()
    nano = nano[nano["map"] == "nano"]
    within, _ = divisions.load_within()
    robust = robustness.load_robustness()
    nano_rho, nano_genes = gene_ranking.load_null_rho("nano")
    nano_null = dict(rho=nano_rho, genes=nano_genes)
    inputs = beyond_density.load_inputs()
    leftover = beyond_density.load_leftover_genes()
    null = beyond_density.load_leftover_null(inputs.structures)
    covariates, _, _ = beyond_density.covariates_for(inputs)
    y = beyond_density.full_map(inputs.nano)
    columns = beyond_density.model_columns(covariates, inputs.terms)
    residual = beyond_density.residual(y, columns)
    print(
        f"main model: {len(inputs.structures)} structures, "
        f"{' + '.join(inputs.terms['abundance'])} + the mean rank of "
        f"{', '.join(beyond_density.DENSITY_GENES)}; {len(leftover)} genes against "
        "its leftover"
    )

    # the genes, and the family with its matched controls
    chosen, members, family, pool, pairs, level = top_genes.choose_genes(
        nano, leftover, q
    )
    members.to_csv(top_genes.FAMILY_TABLE, index=False)
    print(
        f"genes: {int(chosen['top'].sum())} past the map's null after BH, "
        f"{len(family)} of the family's {len(members)} members tested, each matched "
        f"to one of {len(pool)} other postsynaptic genes"
    )

    # what kind of map each gene is, and what it takes of the leftover
    symbols = list(chosen["symbol"])
    annotation, release = top_genes.annotation(symbols)
    print(f"annotation: GO terms named from go-basic.obo, release {release}")
    added, added_nulls = added_shares(inputs, symbols)

    # the tests named for the leftover, and the same against smoother nulls
    tests, group_nulls = top_genes.named_tests(
        family, pairs, leftover, null, nano, nano_null, pool, level
    )
    tests.to_csv(top_genes.NAMED_TESTS, index=False)
    for r in tests.itertuples():
        print(f"tier {r.tier}, {r.map}: {r.test}: {r.first:+.3f}, p {r.p:.4f} ({r.role})")
    check = null_check(inputs, residual, null, leftover, family)
    check.to_csv(top_genes.NULL_CHECK, index=False)

    # the table, and the numbers for the text
    sources = members[["symbol", "sources"]].rename(columns={"sources": "family_sources"})
    table = top_genes.characterise(
        [
            chosen,
            annotation,
            sources,
            top_genes.ranking_columns(nano),
            top_genes.within_columns(within),
            top_genes.robustness_ranges(robust, symbols),
            top_genes.predictor_rhos(inputs, symbols),
            top_genes.leftover_columns(leftover),
            top_genes.family_q_table(family, leftover, nano),
            top_genes.control_columns(pairs, leftover, nano, level),
            added,
        ]
    )
    table.to_csv(top_genes.TOP_TABLE, index=False)
    numbers = top_genes.numbers_table(table, members, tests, check, q)
    numbers.to_csv(top_genes.NUMBERS, index=False)
    past = table[table["family"] & (table["leftover_q_family"] < q)]
    print(
        f"table: {len(table)} genes; family members past BH within the family against "
        f"the leftover: {len(past)}"
    )

    # figures 07, 08, 11 and 11s2, and the gene sheets
    leftover_map = pd.Series(residual, index=inputs.structures)
    n_surrogates = dict(leftover=int(null["rho"].shape[1]), map=int(nano_rho.shape[1]))
    draw_figures(
        table, inputs, leftover, tests, group_nulls, members, leftover_map, n_surrogates
    )
    if sheets:
        draw_sheets(table, inputs, leftover_map, added_nulls, n_surrogates)


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    parser = argparse.ArgumentParser(
        description="the genes that follow the map most; the tests named for the leftover"
    )
    parser.add_argument(
        "--sheets", action="store_true", help="also draw the per-gene sheets"
    )
    args = parser.parse_args()
    main(args.sheets)
