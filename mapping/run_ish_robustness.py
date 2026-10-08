"""The gene ranking of the ISH line under other reasonable choices (A3).

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
    18. run_ish_robustness     A3: the ranking under other choices;   <- this script
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

Correlates every gene with the adult map again, changing one choice of the primary
ranking at a time (statistic, borders, reading, inputs, structure set) and once all
together as the route of 5 October ran, and compares each gene order with the
primary's (the method is in sepmap/ish/robustness.py). Writes, in
adult_v2/ish_analysis/ under the data root:

    tables/ranking_robustness.csv     per variant and gene: rho, structures, ranks
    tables/robustness_summary.csv     per variant: agreement with the primary over
                                      P9's genes and over all, Cacng8's and Gria1's
                                      rho and ranks, the gap between them
    tables/numbers_robustness.csv     the numbers of this step, for the text
    figures/13_robustness.png         the ranking under each choice

    python run_ish_robustness.py

The row without section QC measures the flagged experiments again with nothing set
missing (about a minute).
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sepmap import config, plotting, structures
from sepmap.adult import profiles
from sepmap.ish import gene_ranking, gene_table, regions, robustness, section_qc
from sepmap.ish import plotting as ish_plotting
from sepmap.volumes.per_mouse import structure_terms

OUT = config.DATA / "adult_v2" / "ish_analysis"


def numbers_table(summary):
    """The numbers of this step that the text quotes, one row each."""
    rows = []
    for _, r in summary.iterrows():
        name = r["variant"]
        rows += [
            (f"agreement_p9_{name}", round(r["agreement_p9"], 4), r["label"]),
            (f"agreement_all_{name}", round(r["agreement_all"], 4), r["label"]),
            (f"rank_p9_Gria1_{name}", r["rank_p9_Gria1"], r["label"]),
            (f"rank_p9_Cacng8_{name}", r["rank_p9_Cacng8"], r["label"]),
            (f"rho_Gria1_{name}", round(r["rho_Gria1"], 3), r["label"]),
            (f"rho_Cacng8_{name}", round(r["rho_Cacng8"], 3), r["label"]),
            (f"gap_{name}", round(r["gap"], 3), r["label"]),
        ]
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


def unflagged_rows(table, region):
    """The region table with the flagged experiments measured without section QC."""
    flags = section_qc.flagged_sections(section_qc.load_experiment_qc())
    flagged = table[table["experiment_id"].isin(flags) & ~table["excluded"]]
    ann = regions.annotation_200()
    names, _, _ = structure_terms()
    eroded = regions.eroded_annotation(ann)
    redone = pd.DataFrame(gene_table.region_rows(flagged, {}, ann, names, eroded))
    print(f"no QC: {len(flagged)} flagged experiments measured again", flush=True)
    return robustness.unflagged_region(region, redone)


def variant_inputs(table, genes, profile, set_table, declared):
    """{variant: (map, gene profiles, structures, statistic)} of every row."""
    region = gene_table.load_region_table()
    used = set(table.loc[~table["excluded"], "experiment_id"])
    region = region[region["experiment_id"].isin(used)]
    per_experiment = gene_table.experiment_profiles(region)
    eroded = robustness.merged_profiles(
        gene_table.experiment_profiles(region, "ish_mean_eroded")
    )
    no_qc = robustness.merged_profiles(
        gene_table.experiment_profiles(unflagged_rows(table, region))
    )
    p9_own = dict(zip(genes["symbol"], genes["p9_experiment_id"]))
    merged = gene_table.load_profiles()
    zref = profile["zref_nano"]
    zref_eroded = profile["zref_nano_eroded"]
    stored_zref = robustness.stored_map("zref")
    every = sorted(set_table["structure"])
    return {
        "primary": (zref, merged, declared, "spearman"),
        "pearson_log2": (
            zref,
            robustness.log2_profiles(per_experiment),
            declared,
            "pearson",
        ),
        "eroded_ish": (zref, eroded, declared, "spearman"),
        "eroded_nano": (zref_eroded, merged, declared, "spearman"),
        "eroded_both": (zref_eroded, eroded, declared, "spearman"),
        "ratio": (robustness.stored_map("ratio"), merged, declared, "spearman"),
        "stored_zref": (stored_zref, merged, declared, "spearman"),
        "no_qc": (zref, no_qc, declared, "spearman"),
        "p9_experiment": (
            zref,
            robustness.single_experiment(per_experiment, p9_own),
            declared,
            "spearman",
        ),
        "p9_divisions": (
            zref,
            merged,
            robustness.p9_division_structures(set_table),
            "spearman",
        ),
        "every_structure": (zref, merged, every, "spearman"),
        "before_build": (
            stored_zref,
            robustness.frozen_profiles(),
            sorted(stored_zref.index),
            "spearman",
        ),
    }


def main():
    """Print the settings, then rank the genes under every variant; tables, figure."""
    config.print_settings({})
    tables = OUT / "tables"

    # the inputs of the primary and of every variant
    declared = structures.declared_structures()
    set_table = structures.load_structure_set()
    profile = profiles.load_profile()
    table = gene_table.load_gene_table()
    genes = gene_table.per_gene(table)
    p9_genes = set(genes.loc[genes["p9_gene"], "symbol"])
    subunits = {s for s, t in zip(genes["symbol"], genes["gene_sets"]) if "subunits" in t}
    inputs = variant_inputs(table, genes, profile, set_table, declared)
    print(f"variants: {len(inputs)}, {len(p9_genes)} of P9's genes among {len(genes)}")

    # every variant's rho per gene; the primary must be analysis 1's ranking
    parts, gaps = [], {}
    for name, _, _ in robustness.VARIANTS:
        map_values, gene_profiles, structure_list, method = inputs[name]
        part = robustness.variant_rows(
            name, map_values, gene_profiles, structure_list, p9_genes, method
        )
        parts.append(part)
        if all(g in gene_profiles for g in gene_ranking.GAP_GENES):
            gaps[name] = robustness.gap_on_shared(
                map_values,
                gene_profiles[gene_ranking.GAP_GENES[0]],
                gene_profiles[gene_ranking.GAP_GENES[1]],
                structure_list,
                method,
            )
        else:
            gaps[name] = np.nan
    per_gene = pd.concat(parts, ignore_index=True)
    primary = per_gene[per_gene["variant"] == "primary"]
    ranking = gene_ranking.load_ranking()
    nano = ranking[ranking["map"] == "nano"].set_index("symbol")["rho"]
    largest = float((primary.set_index("symbol")["rho"] - nano).abs().max())
    if largest > 1e-9:
        raise ValueError(
            f"the primary row differs from gene_ranking.csv by up to {largest:.2e}; "
            "run run_ish_gene_ranking.py first"
        )
    per_gene.to_csv(robustness.ROBUSTNESS, index=False)

    # the summary: agreement, where the two genes sit, the gap
    summary = pd.DataFrame(
        [
            robustness.summary_row(
                v, per_gene[per_gene["variant"] == v[0]], primary, p9_genes, gaps[v[0]]
            )
            for v in robustness.VARIANTS
        ]
    )
    summary.to_csv(robustness.SUMMARY, index=False)
    for _, r in summary.iterrows():
        print(
            f"  {r['variant']:16s} agreement {r['agreement_p9']:.3f} over "
            f"{r['n_agreement_p9']} of P9's genes; Gria1 {r['rank_p9_Gria1']:.0f}, "
            f"Cacng8 {r['rank_p9_Cacng8']:.0f}, gap {r['gap']:+.3f}"
        )
    numbers_table(summary).to_csv(tables / "numbers_robustness.csv", index=False)

    # figure 13
    path = OUT / "figures" / ish_plotting.figure_file("robustness")
    gap = pd.read_csv(gene_ranking.GAP)
    merged = gap[gap["kind"] == "merged profiles"].iloc[0]
    band = (float(merged["equal_lo"]), float(merged["equal_hi"]))
    fig = ish_plotting.plot_robustness(summary, per_gene, subunits, band, save=path)
    plt.close(fig)
    print(f"figure: {path}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    # no options; parsing still gives the script its --help
    parser = argparse.ArgumentParser(description="the gene ranking under other choices")
    parser.parse_args()
    main()
