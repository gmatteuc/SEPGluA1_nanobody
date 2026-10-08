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
                               documentation; figure 02
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the    <- this script
                               map, the null, autofluorescence (A8),
                               the Cacng8 - Gria1 gap; figures 05 to 07 and 12
    18. run_ish_robustness     A3: the ranking under other choices;
                               figure 13
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
    25. run_beyond_figures     figures 03, 04 and 11
    26. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 14
    27. run_ish_overview       figures 00 and 15, the figure index;
                               the numbers for the text

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
    figures/05_one_comparison.png      what one comparison is
    figures/06_spatial_null.png        why a null, and the null itself
    figures/07_gene_ranking.png        P9's genes against the map and its null
    figures/12_autofluorescence.png    the same on the autofluorescence map

    python run_ish_gene_ranking.py
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sepmap import config, plotting, structures
from sepmap.adult import profiles
from sepmap.ish import gene_ranking, gene_table, planes, spatial_null
from sepmap.ish import plotting as ish_plotting
from sepmap.volumes.per_mouse import structure_terms

ISH_ANALYSIS = config.SETTINGS["ish_analysis"]
ISH_FIGURES = config.SETTINGS["ish_figures"]

OUT = config.DATA / "adult_v2" / "ish_analysis"

# the genes of figure 05 beside the map: the top of P9's ranking, the subunit the
# label is on, and a glial gene as the control
COMPARISON_GENES = ("Cacng8", "Gria1", "Aqp4")


def numbers_table(ranking, gap, per_adult):
    """The numbers of this step that the text quotes, one row each."""
    q = ISH_ANALYSIS["q"]
    rows = []
    for name in gene_ranking.MAPS:
        mine = ranking[ranking["map"] == name]
        p9 = mine[mine["p9_gene"]]
        rows += [
            (f"{name}_genes", len(mine), f"genes correlated with the {name} map"),
            (f"{name}_p9_genes", len(p9), f"P9's genes correlated with the {name} map"),
            (
                f"{name}_pass_p9",
                int((p9["q_p9"] < q).sum()),
                f"P9's genes past the null, BH within P9's genes, q < {q}",
            ),
            (
                f"{name}_pass_all",
                int((mine["q_all"] < q).sum()),
                f"genes past the null, BH within all genes, q < {q}",
            ),
            (
                f"{name}_median_rho",
                round(float(mine["rho"].median()), 3),
                "median rho over all genes",
            ),
        ]
        for gene in ("Cacng8", "Gria1", "Aqp4"):
            r = mine[mine["symbol"] == gene].iloc[0]
            rows += [
                (f"{name}_rho_{gene}", round(r["rho"], 3), f"{gene}'s rho"),
                (f"{name}_p_{gene}", round(r["p_spatial"], 5), f"{gene}'s spatial p"),
                (f"{name}_q_p9_{gene}", round(r["q_p9"], 5), f"{gene}'s q within P9"),
                (f"{name}_rank_p9_{gene}", r["rank_p9"], f"{gene}'s rank among P9's"),
                (f"{name}_rank_all_{gene}", r["rank_all"], f"{gene}'s rank among all"),
            ]
    merged = gap[gap["kind"] == "merged profiles"].iloc[0]
    pairs = gap[gap["kind"] == "experiment pairing"]
    rows += [
        ("gap", round(merged["gap"], 3), "rho(Cacng8) - rho(Gria1), merged profiles"),
        ("gap_structures", int(merged["n_structures"]), "structures both genes have"),
        ("gap_p", round(merged["p_equal"], 5), "its p, maps related to both alike"),
        ("gap_equal_lo", round(merged["equal_lo"], 3), "2.5% of that null"),
        ("gap_equal_hi", round(merged["equal_hi"], 3), "97.5% of that null"),
        ("gap_p_unrelated", round(merged["p_spatial"], 5), "its p, unrelated maps"),
        ("gap_unrelated_lo", round(merged["null_lo"], 3), "2.5% of that null"),
        ("gap_unrelated_hi", round(merged["null_hi"], 3), "97.5% of that null"),
        ("gap_boot_lo", round(merged["boot_lo"], 3), "2.5% over resampled adults"),
        ("gap_boot_hi", round(merged["boot_hi"], 3), "97.5% over resampled adults"),
        ("gap_pairs_min", round(pairs["gap"].min(), 3), "smallest over pairings"),
        ("gap_pairs_max", round(pairs["gap"].max(), 3), "largest over pairings"),
        ("gap_pairs_p_min", round(pairs["p_equal"].min(), 5), "smallest p over pairings"),
        ("gap_pairs_p_max", round(pairs["p_equal"].max(), 5), "largest p over pairings"),
    ]
    for gene in ("Gria1", "Cacng8"):
        mine = per_adult[per_adult["symbol"] == gene]
        for name in gene_ranking.MAPS:
            rows += [
                (
                    f"per_adult_{name}_{gene}_min",
                    round(mine[f"rho_{name}"].min(), 3),
                    f"lowest rho of one adult's {name} map with {gene}",
                ),
                (
                    f"per_adult_{name}_{gene}_max",
                    round(mine[f"rho_{name}"].max(), 3),
                    f"highest rho of one adult's {name} map with {gene}",
                ),
            ]
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


def comparison_rows(ranking, merged, table, lab, plane):
    """The genes of figure 05: their ISH section on the plane and their rows."""
    nano = ranking[ranking["map"] == "nano"].set_index("symbol")
    n_p9 = int(nano["p9_gene"].sum())
    rows = []
    for gene in COMPARISON_GENES:
        used = table[(table["symbol"] == gene) & ~table["excluded"]]
        coronal = used[used["plane"] == "coronal"]
        own = coronal[coronal["p9_experiment"]]
        shown = (own if len(own) else coronal).iloc[0]
        flagged = [int(k) for k in str(shown["flagged_sections"]).split()]
        image = planes.ish_plane(
            shown["experiment_id"], "coronal", flagged, plane, lab.shape
        )
        caption = (
            f"Allen experiment {shown['experiment_id']}, one of the gene's "
            f"{len(used)} (200 um)"
        )
        if len(used) == 1:
            caption = f"Allen experiment {shown['experiment_id']}, its only one (200 um)"
        rows.append(
            dict(
                symbol=gene,
                image=image,
                caption=caption,
                profile=merged[gene],
                ranking=nano.loc[gene],
                n_p9=n_p9,
            )
        )
    return rows


def main():
    """Print the settings, then rank every gene against both maps; tables, figures."""
    config.print_settings({})
    tables = OUT / "tables"
    figures = OUT / "figures"

    # the declared structures, the two maps per adult and as cohort means, surrogates
    declared = structures.declared_structures()
    set_table = structures.load_structure_set()
    profile = profiles.load_profile()
    per_mouse = profiles.load_per_mouse()
    surr = {}
    for name in gene_ranking.MAPS:
        surr[name], listed = spatial_null.load_surrogates(name)
        if listed != declared:
            raise ValueError(
                "the surrogates were drawn on other structures than the declared set; "
                "run run_ish_spatial_null.py --recompute"
            )
    n_surrogates = surr["nano"].shape[0]
    print(f"maps: nano and auto on {len(declared)} structures, {n_surrogates} surrogates")

    # the genes: merged profiles and labels
    table = gene_table.load_gene_table()
    genes = gene_table.per_gene(table)
    merged = gene_table.load_profiles()
    vectors = gene_ranking.gene_vectors(merged, declared)
    subunits = {s for s, t in zip(genes["symbol"], genes["gene_sets"]) if "subunits" in t}
    print(f"genes: {len(vectors)} of {len(genes)} share enough declared structures")

    # every gene against each map: rho, spatial p, null band, spread over adults
    results, nulls, adults = {}, {}, {}
    for k, (name, column) in enumerate(gene_ranking.MAPS.items()):
        adults[name] = gene_ranking.adult_matrix(
            per_mouse, column, declared, profiles.ADULTS
        )
        boot = gene_ranking.bootstrap_maps(adults[name], seed=k)
        map_values = profile.loc[declared, column].to_numpy(float)
        results[name], nulls[name] = gene_ranking.rank_genes(
            map_values, surr[name], boot, vectors
        )
        print(f"  {name}: {len(results[name])} genes tested", flush=True)
    ranking = gene_ranking.ranking_table(results, genes)
    ranking.to_csv(gene_ranking.RANKING, index=False)
    q = ISH_ANALYSIS["q"]
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
    region = gene_table.load_region_table()
    used = set(table.loc[~table["excluded"], "experiment_id"])
    per_experiment = gene_table.experiment_profiles(
        region[region["experiment_id"].isin(used)]
    )
    boot = gene_ranking.bootstrap_maps(adults["nano"], seed=0)
    gap, gap_null = gene_ranking.gap_table(
        profile.loc[declared, "zref_nano"].to_numpy(float),
        surr["nano"],
        boot,
        merged,
        per_experiment,
        declared,
    )
    gap.to_csv(gene_ranking.GAP, index=False)
    np.savez(gene_ranking.NULL_RHO, genes=np.array(list(vectors)), gap=gap_null, **nulls)
    first = gap.iloc[0]
    print(
        f"gap: {first['gap']:+.3f} on {first['n_structures']} structures, p "
        f"{first['p_equal']:.4f} against maps related to both alike (c "
        f"{first['equal_weight']:.3f}), {first['p_spatial']:.4f} against unrelated "
        f"maps; {len(gap) - 1} pairings of experiments"
    )

    # the numbers for the text
    numbers = numbers_table(ranking, gap, per_adult)
    numbers.to_csv(tables / "numbers_gene_ranking.csv", index=False)

    # figure 05: one comparison
    plane = ISH_FIGURES["plane"]
    names, _, _ = structure_terms()
    lab = planes.label_plane(plane)
    map_values = profile.loc[declared, "zref_nano"]
    fig = ish_plotting.plot_one_comparison(
        map_values,
        planes.nano_plane(plane),
        comparison_rows(ranking, merged, table, lab, plane),
        set_table,
        lab,
        names,
        plane,
        n_surrogates,
        save=figures / ish_plotting.figure_file("one_comparison"),
    )
    plt.close(fig)

    # figure 06: the spatial null, with the two genes against it
    n = len(declared)
    map_ranks = dict(zip(declared, ish_plotting.ranks01(map_values.to_numpy())))
    surrogate_ranks = [
        dict(zip(declared, (surr["nano"][k] - 1) / (n - 1))) for k in range(3)
    ]
    nano = ranking[ranking["map"] == "nano"].set_index("symbol")
    tested = list(vectors)
    shown = [
        (g, nano.loc[g, "rho"], nulls["nano"][tested.index(g)], nano.loc[g, "p_spatial"])
        for g in gene_ranking.GAP_GENES
    ]
    fig = ish_plotting.plot_spatial_null(
        pd.read_csv(spatial_null.VARIOGRAM),
        pd.read_csv(spatial_null.CALIBRATION),
        map_ranks,
        surrogate_ranks,
        lab,
        names,
        plane,
        n_surrogates,
        genes=shown,
        save=figures / ish_plotting.figure_file("spatial_null"),
    )
    plt.close(fig)

    # figure 07: the ranking; figure 12: the autofluorescence map
    fig = ish_plotting.plot_gene_ranking(
        ranking,
        gap,
        gap_null,
        subunits,
        ISH_FIGURES["t_max"],
        q,
        n_surrogates,
        save=figures / ish_plotting.figure_file("gene_ranking"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_autofluorescence(
        ranking,
        per_adult,
        subunits,
        q,
        n_surrogates,
        save=figures / ish_plotting.figure_file("autofluorescence"),
    )
    plt.close(fig)
    print(f"figures: 05 to 07 and 12 in {figures}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="analysis 1: each gene against the map")
    parser.parse_args()
    main()
