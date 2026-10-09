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
                               its coverage of the fit; figure 14s
    22. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    23. run_beyond_controls    seven attempts to break it
    24. run_beyond_calibration the same model on maps whose answer
                               is known
    25. run_beyond_regression  the regression, per structure
    26. run_beyond_figures     figures 03, 03s1, 03s2, 04, 11s1, 14
    27. run_ish_top_genes      the genes that follow the map most,   <- this script
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 07, 08, 11, 11s2,
                               sheets
    28. run_sep_channel_check  analysis 5: what the green channel
                               reports; figures 15, 15s
    29. run_ish_overview       figures 00, 16 and 16s, the figure
                               index; the numbers for the text

Describes the genes past the spatial null after BH in analysis 1, Gria1, Cacng8 and
every member of the AMPA receptor complex family on every axis of the ISH line, adds
each to the main model of analysis 4 against maps of its own smoothness, and runs the
tests named on 8 October for the leftover: Cacng8 alone, then the family as a group
and gene by gene, with two check rows of the group test and the tests again against
smoother nulls (the method is in sepmap/ish/top_genes.py). It reads the tables of
steps 15 to 22. Writes, in adult_v2/ish_analysis/ under the data root:

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
from sepmap.adult import beyond_density as bd
from sepmap.adult import profiles
from sepmap.ish import (
    divisions,
    gene_ranking,
    gene_sets,
    gene_table,
    planes,
    robustness,
    top_genes,
)
from sepmap.ish import plotting as ish_plotting
from sepmap.volumes.per_mouse import structure_terms

ISH_ANALYSIS = config.SETTINGS["ish_analysis"]
ISH_FIGURES = config.SETTINGS["ish_figures"]

OUT = config.DATA / "adult_v2" / "ish_analysis"


def load_leftover(inputs):
    """leftover_genes.csv and its null, checked against the main model's structures."""
    leftover = pd.read_csv(bd.LEFTOVER_GENES, keep_default_na=False, na_values=[""])
    leftover["in_model"] = leftover["in_model"].fillna("")
    with np.load(bd.LEFTOVER_NULL) as z:
        null = dict(
            structures=[str(s) for s in z["structures"]],
            surrogates=z["surrogates"].astype(float),
            rho=z["null_rho"],
            genes=[str(g) for g in z["genes"]],
        )
    if null["structures"] != inputs.structures:
        raise ValueError(
            f"{bd.LEFTOVER_NULL} was drawn on other structures than the main model's "
            f"{len(inputs.structures)}; run run_beyond_density.py first"
        )
    return leftover, null


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


def draw_figures(out, inputs, leftover, tests, group_nulls, members, residual, n):
    """Figures 07 and 08, and 11 and 11s2 on the main model's leftover.

    `n` holds the leftover's and the map's numbers of surrogates.
    """
    q = ISH_ANALYSIS["q"]
    figures = OUT / "figures"
    set_table = structures.load_structure_set()
    map_values = profiles.load_profile().loc[
        structures.declared_structures(), "zref_nano"
    ]
    fig = ish_plotting.plot_top_genes(
        out, q, ISH_FIGURES["t_max"], save=figures / ish_plotting.figure_file("top_genes")
    )
    plt.close(fig)
    gap = pd.read_csv(gene_ranking.GAP)
    gap_null, _ = gene_ranking.load_null_rho("gap")
    fig = ish_plotting.plot_cacng8_gria1(
        out.set_index("symbol"),
        map_values,
        inputs.expr,
        gap[gap["kind"] == "merged profiles"].iloc[0],
        gap_null[0],
        set_table,
        n["map"],
        save=figures / ish_plotting.figure_file("cacng8_gria1"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_leftover(
        out,
        residual,
        inputs.expr[gene_sets.LEFTOVER_GENE],
        tests,
        group_nulls,
        leftover,
        set_table,
        q,
        n["leftover"],
        save=figures / ish_plotting.figure_file("leftover"),
    )
    plt.close(fig)
    fig = ish_plotting.plot_ampa_family(
        out[out["family"]],
        tests,
        group_nulls,
        members,
        q,
        ISH_FIGURES["t_max"],
        n["leftover"],
        n["map"],
        save=figures / ish_plotting.figure_file("ampa_family"),
    )
    plt.close(fig)
    keys = ("top_genes", "cacng8_gria1", "leftover", "ampa_family")
    drawn = [ish_plotting.figure_file(k) for k in keys]
    print(f"figures: {', '.join(drawn)} in {figures}")


def draw_sheets(out, inputs, residual, added_nulls, n):
    """One sheet per gene characterised, in figures/top_genes/."""
    plane = ISH_FIGURES["plane"]
    structure_names, _, _ = structure_terms()
    lab = planes.label_plane(plane)
    set_table = structures.load_structure_set()
    map_values = profiles.load_profile().loc[
        structures.declared_structures(), "zref_nano"
    ]
    folder = OUT / "figures" / "top_genes"
    folder.mkdir(parents=True, exist_ok=True)
    for _, row in out.iterrows():
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
    print(f"gene sheets: {len(out)} in {folder}")


def choose_genes(nano, leftover, q):
    """The genes characterised, the family's members, and each member's control.

    Returns the chosen genes, family_members.csv's table, the members tested, the
    pool of other postsynaptic genes against the leftover, the pairs and every
    gene's median energy.
    """
    table = gene_table.load_gene_table()
    genes = gene_table.per_gene(table)
    chosen = top_genes.chosen_genes(nano, leftover, q)
    members = top_genes.family_members(genes, leftover)
    family = list(members.loc[members["tested"], "symbol"])
    _, level = gene_table.gene_levels(table)
    tested = set(leftover["symbol"])
    postsynaptic = gene_sets.members_from_table(genes)["other postsynaptic"]
    pool = [g for g in postsynaptic if g in tested]
    pairs = top_genes.family_controls(family, pool, level)
    return chosen, members, family, pool, pairs, level


def annotate(symbols):
    """What the gene table says each gene is, with GO's names of its terms."""
    documentation = pd.read_csv(
        gene_table.DOCUMENTATION, encoding="utf-8-sig", keep_default_na=False
    )
    records = gene_table.mygene_records(symbols, offline=True)
    parents, names, release = gene_table.load_obo(offline=True)
    print(f"annotation: GO terms named from go-basic.obo, release {release}")
    return top_genes.annotation_table(symbols, documentation, records, parents, names)


def named_tests(family, pairs, leftover, null, nano, nano_null, pool, level):
    """Tier 1, and tier 2 as a group with its check rows, on the leftover; the group
    tests on the map, as a description.

    `null` and `nano_null` hold each map's null rho and the genes of its rows.
    Returns named_tests.csv's table and the group tests' nulls.
    """
    maps = dict(
        leftover=(leftover.set_index("symbol")["rho"], null["rho"], null["genes"]),
        nano=(nano.set_index("symbol")["rho"], nano_null["rho"], nano_null["genes"]),
    )
    group, group_nulls = top_genes.group_tests(family, pairs, maps)
    checks = top_genes.check_tests(
        family, leftover, null["rho"], null["genes"], pool, level
    )
    tier_one = pd.DataFrame([top_genes.tier_one(leftover)])
    tests = pd.concat([tier_one, group, checks], ignore_index=True)
    for r in tests.itertuples():
        print(f"tier {r.tier}, {r.map}: {r.test}: {r.first:+.3f}, p {r.p:.4f} ({r.role})")
    return tests, group_nulls


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


def main(sheets):
    """Print the settings; characterise the genes, run the named tests, draw."""
    config.print_settings({"sheets": sheets})
    q = ISH_ANALYSIS["q"]

    # the tables of analyses 1 to 4, and the main model with its leftover's null
    nano = gene_ranking.load_ranking()
    nano = nano[nano["map"] == "nano"]
    within, _ = divisions.load_within()
    robust = pd.read_csv(robustness.ROBUSTNESS)
    nano_rho, nano_genes = gene_ranking.load_null_rho("nano")
    nano_null = dict(rho=nano_rho, genes=nano_genes)
    inputs = bd.load_inputs()
    leftover, null = load_leftover(inputs)
    covariates, _, _ = bd.covariates_for(inputs)
    y = bd.full_map(inputs.nano)
    residual = bd.residual(y, bd.model(covariates, inputs.terms))
    print(
        f"main model: {len(inputs.structures)} structures, density "
        f"{', '.join(inputs.terms['density'])}; {len(leftover)} genes against its "
        "leftover"
    )

    # the genes, and the family with its matched controls
    chosen, members, family, pool, pairs, level = choose_genes(nano, leftover, q)
    members.to_csv(top_genes.FAMILY_TABLE, index=False)
    print(
        f"genes: {int(chosen['top'].sum())} past the map's null after BH, "
        f"{len(family)} of the family's {len(members)} members tested, each matched "
        f"to one of {len(pool)} other postsynaptic genes"
    )

    # what kind of map each gene is, and what it takes of the leftover
    symbols = list(chosen["symbol"])
    annotation = annotate(symbols)
    added, added_nulls = added_shares(inputs, symbols)

    # the tests named for the leftover, and the same against smoother nulls
    tests, group_nulls = named_tests(
        family, pairs, leftover, null, nano, nano_null, pool, level
    )
    tests.to_csv(top_genes.NAMED_TESTS, index=False)
    check = null_check(inputs, residual, null, leftover, family)
    check.to_csv(top_genes.NULL_CHECK, index=False)

    # the table, and the numbers for the text
    sources = members[["symbol", "sources"]].rename(columns={"sources": "family_sources"})
    out = top_genes.characterise(
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
    for column in ("family_sources", "in_model", "matched_control"):
        out[column] = out[column].fillna("")
    out.to_csv(top_genes.TOP_TABLE, index=False)
    numbers = top_genes.numbers_table(out, members, tests, check, q)
    numbers.to_csv(top_genes.NUMBERS, index=False)
    past = out[out["family"] & (out["leftover_q_family"] < q)]
    print(
        f"table: {len(out)} genes; family members past BH within the family against "
        f"the leftover: {len(past)}"
    )

    # figures 07, 08, 11 and 11s2, and the gene sheets
    residual = pd.Series(residual, index=inputs.structures)
    n = dict(leftover=int(null["rho"].shape[1]), map=int(nano_rho.shape[1]))
    draw_figures(out, inputs, leftover, tests, group_nulls, members, residual, n)
    if sheets:
        draw_sheets(out, inputs, residual, added_nulls, n)


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
