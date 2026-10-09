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
                               its coverage of the fit; figure 14s
    22. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    23. run_beyond_controls    seven attempts to break it
    24. run_beyond_calibration the same model on maps whose answer
                               is known
    25. run_beyond_regression  the regression, per structure
    26. run_beyond_figures     figures 03, 03s1, 03s2, 04, 11s1, 14
    27. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 07, 08, 11, 11s2,
                               sheets
    28. run_sep_channel_check  analysis 5: what the green channel
                               reports; figures 15, 15s
    29. run_ish_overview       figures 00, 16 and 16s, the figure
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
from sepmap.ish import divisions, gene_ranking, gene_table, planes, spatial_null
from sepmap.ish import plotting as ish_plotting
from sepmap.volumes.per_mouse import structure_terms

ISH_ANALYSIS = config.SETTINGS["ish_analysis"]
ISH_FIGURES = config.SETTINGS["ish_figures"]

OUT = config.DATA / "adult_v2" / "ish_analysis"


def numbers_table(within, detail, calibration):
    """The numbers of this step that the text quotes, one row each."""
    q = ISH_ANALYSIS["q"]
    have = within[within["rho_within"].notna()]
    p9 = within["p9_gene"]
    rows = [
        ("genes", len(within), "genes split into between and within"),
        ("genes_within", len(have), "genes in a division with enough structures"),
        (
            "agreement_division_only_all",
            round(
                divisions.agreement(
                    within, "rho_division_only", within["symbol"].notna()
                ),
                3,
            ),
            "Spearman over genes, whole-brain against division-only rho",
        ),
        (
            "agreement_division_only_p9",
            round(divisions.agreement(within, "rho_division_only", p9), 3),
            "the same over P9's genes",
        ),
        (
            "agreement_within_all",
            round(divisions.agreement(within, "rho_within", within["symbol"].notna()), 3),
            "Spearman over genes, whole-brain against within rho",
        ),
        (
            "agreement_within_p9",
            round(divisions.agreement(within, "rho_within", p9), 3),
            "the same over P9's genes",
        ),
        ("median_rho", round(within["rho"].median(), 3), "median whole-brain rho"),
        (
            "median_rho_within",
            round(have["rho_within"].median(), 3),
            "median within rho",
        ),
        (
            "median_rho_p9",
            round(within.loc[p9, "rho"].median(), 3),
            "median whole-brain rho, P9's genes",
        ),
        (
            "median_rho_within_p9",
            round(within.loc[p9, "rho_within"].median(), 3),
            "median within rho, P9's genes",
        ),
    ]
    for null in ("spatial", "shuffle"):
        rows += [
            (
                f"pass_all_{null}",
                int((within[f"q_all_{null}"] < q).sum()),
                f"genes past the within {null} null, BH within all, q < {q}",
            ),
            (
                f"pass_p9_{null}",
                int((within[f"q_p9_{null}"] < q).sum()),
                f"P9's genes past the within {null} null, BH within P9's",
            ),
            (
                f"calibration_{null}",
                round(float((calibration[f"p_{null}"] < 0.05).mean()), 4),
                f"random maps with p < 0.05 by the {null} null ({len(calibration)})",
            ),
        ]
    for gene in divisions.DETAIL_GENES:
        r = within[within["symbol"] == gene].iloc[0]
        rows += [
            (f"rho_division_only_{gene}", round(r["rho_division_only"], 3), gene),
            (f"rho_within_{gene}", round(r["rho_within"], 3), gene),
            (f"p_within_spatial_{gene}", round(r["p_within_spatial"], 5), gene),
            (f"p_within_shuffle_{gene}", round(r["p_within_shuffle"], 5), gene),
        ]
    divisions_used = sorted(set(detail["division"]))
    rows.append(("divisions_used", " ".join(divisions_used), "divisions entered"))
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


def within_calibration(map_values, division, surr, declared):
    """The two nulls of the within rho on random smooth maps with the map's
    smoothness (the spatial null's calibration fields)."""
    d = spatial_null.distance_matrix(structures.load_centroids().loc[declared])
    params = spatial_null.fit_exponential(map_values.to_numpy(float), d)
    fields = spatial_null.random_fields(
        d, params, ISH_ANALYSIS["n_calibration_within"], np.random.default_rng(1)
    )
    return divisions.calibration(map_values.to_numpy(float), division, surr, fields)


def draw_figures(within, detail, map_values, division, merged, calibration, subunits):
    """Figure 09 and its detailed version."""
    figures = OUT / "figures"
    set_table = structures.load_structure_set()
    q = ISH_ANALYSIS["q"]
    coarse = pd.Series(
        divisions.division_only(map_values.to_numpy(float), division),
        index=map_values.index,
    )
    min_structures = ISH_ANALYSIS["min_division_structures"]
    fig = ish_plotting.plot_between_within(
        within,
        map_values,
        coarse,
        merged["Cacng8"],
        "Cacng8",
        set_table,
        subunits,
        divisions.DETAIL_GENES,
        q,
        min_structures,
        save=figures / ish_plotting.figure_file("between_within"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_between_within_detail(
        within,
        detail,
        map_values,
        coarse,
        merged["Cacng8"],
        "Cacng8",
        set_table,
        calibration,
        subunits,
        divisions.DETAIL_GENES,
        q,
        min_structures,
        save=figures / ish_plotting.figure_file("between_within_detail"),
    )
    plt.close(fig)
    keys = ("between_within", "between_within_detail")
    drawn = [ish_plotting.figure_file(k) for k in keys]
    print(f"figures: {', '.join(drawn)} in {figures}")


def draw_sheets(within, detail, map_values, merged, ranking, n_surrogates):
    """One sheet per gene of divisions.DETAIL_GENES, in figures/genes/."""
    set_table = structures.load_structure_set()
    plane = ISH_FIGURES["plane"]
    names, _, _ = structure_terms()
    lab = planes.label_plane(plane)
    folder = OUT / "figures" / "genes"
    folder.mkdir(parents=True, exist_ok=True)
    by_gene = within.set_index("symbol")
    nano_rows = ranking[ranking["map"] == "nano"].set_index("symbol")
    for gene in divisions.DETAIL_GENES:
        fig = ish_plotting.plot_gene_sheet(
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
    tables = OUT / "tables"

    # the declared structures, their divisions, the map and its surrogates
    declared = structures.declared_structures()
    set_table = structures.load_structure_set()
    division_of = dict(zip(set_table["structure"], set_table["division"]))
    division = np.array([division_of[s] for s in declared])
    profile = profiles.load_profile()
    map_values = profile.loc[declared, "zref_nano"]
    surr, listed = spatial_null.load_surrogates("nano")
    if listed != declared:
        raise ValueError(
            "the surrogates were drawn on other structures than the declared set; "
            "run run_ish_spatial_null.py --recompute"
        )
    print(f"map: {len(declared)} structures in {len(set(division))} divisions")

    # the genes
    genes = gene_table.per_gene(gene_table.load_gene_table())
    p9_genes = set(genes.loc[genes["p9_gene"], "symbol"])
    subunits = {s for s, t in zip(genes["symbol"], genes["gene_sets"]) if "subunits" in t}
    merged = gene_table.load_profiles()
    vectors = gene_ranking.gene_vectors(merged, declared)

    # between and within, per gene, with the two nulls of the within rho
    within, detail = divisions.gene_rows(
        map_values.to_numpy(float), division, surr, vectors
    )
    within = divisions.with_q(within, p9_genes)
    ranking = gene_ranking.load_ranking()
    nano = ranking[ranking["map"] == "nano"].set_index("symbol")["rho"]
    largest = float((within.set_index("symbol")["rho"] - nano).abs().max())
    if largest > 1e-9:
        raise ValueError(
            f"the whole-brain rho differs from gene_ranking.csv by up to {largest:.2e}; "
            "run run_ish_gene_ranking.py first"
        )
    within.to_csv(divisions.WITHIN, index=False)
    detail.to_csv(divisions.WITHIN_DETAIL, index=False)
    q = ISH_ANALYSIS["q"]
    print(
        f"within: {int(within['rho_within'].notna().sum())} genes; median "
        f"{within['rho_within'].median():+.3f} against {within['rho'].median():+.3f} "
        f"whole-brain; past the spatial null at q < {q}: "
        f"{int((within['q_all_spatial'] < q).sum())} of all"
    )

    # the two nulls on random smooth maps with the map's smoothness
    calibration = within_calibration(map_values, division, surr, declared)
    calibration.to_csv(tables / "within_calibration.csv", index=False)
    rates = ", ".join(
        f"{null} {float((calibration[f'p_{null}'] < 0.05).mean()):.1%}"
        for null in ("shuffle", "spatial")
    )
    print(f"calibration: {len(calibration)} random maps, p < 0.05 by {rates}")

    # the numbers for the text
    numbers = numbers_table(within, detail, calibration)
    numbers.to_csv(tables / "numbers_divisions.csv", index=False)

    # figure 09 and its detailed version, and the gene sheets
    draw_figures(within, detail, map_values, division, merged, calibration, subunits)
    if sheets:
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
