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
    27. run_ish_top_genes      the genes that follow the map most,  <- this script
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 12 to 14, sheets
    28. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 17
    29. run_ish_overview       figures 00 and 18, the figure index;
                               the numbers for the text

Describes the genes past the spatial null after BH in analysis 1, Gria1, Cacng8 and
every member of the AMPA receptor complex family on every axis of the ISH line, adds
each to the main model of analysis 4 against maps of its own smoothness, and runs the
tests named on 8 October for the leftover: Cacng8 alone, then the family as a group
and gene by gene (the method is in sepmap/ish/top_genes.py). It reads the tables of
steps 15 to 22. Writes, in adult_v2/ish_analysis/ under the data root:

    tables/top_genes.csv              one row per gene characterised
    tables/named_tests.csv            the tests named for the leftover, and the
                                      same read on the map
    tables/family_members.csv         the family: each member's sources, tested or
                                      why not
    tables/numbers_top_genes.csv      the numbers of this step, for the text
    figures/12_top_genes.png          the genes that follow the map
    figures/13_cacng8_gria1.png       Cacng8 against Gria1
    figures/14_ampa_family.png        the AMPA receptor complex against the leftover
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

# the columns of a gene's row that the numbers file carries for Cacng8 and Gria1
GENE_NUMBERS = (
    ("rho", 3, "rho with the map"),
    ("p_spatial", 6, "its spatial p"),
    ("rank_all", 0, "its rank among all genes"),
    ("rho_within", 3, "mean rho inside divisions"),
    ("p_within", 6, "its spatial p"),
    ("rho_variants_min", 3, "lowest rho over the robustness variants"),
    ("rho_variants_max", 3, "highest rho over the robustness variants"),
    ("rank_variants_min", 0, "best rank over the variants holding every gene"),
    ("rank_variants_max", 0, "worst rank over the variants holding every gene"),
    ("reliability", 3, "agreement of its Allen experiments"),
    ("rho_Gria1", 3, "rho with Gria1, structures of the fit"),
    ("rho_markers", 3, "rho with the marker composite"),
    ("rho_psd_pc1", 3, "rho with psd_pc1"),
    ("rho_prediction", 3, "rho with the main model's prediction"),
    ("rho_psd95", 3, "rho with PSD95 punctum density, where measured"),
    ("leftover_rho", 3, "rho with the main model's leftover"),
    ("leftover_p", 6, "its spatial p, each surrogate through the same fit"),
    ("leftover_rank", 0, "its rank against the leftover"),
    ("left_main", 4, "share of the reproducible map the main model leaves there"),
    ("left_with_gene", 4, "the same with the gene added"),
    ("taken", 4, "share of the reproducible map the gene takes from the leftover"),
    ("taken_null_hi", 4, "95% of plain surrogates of the gene take less than this"),
    ("p_taken", 6, "share of those maps that take as much, one-sided"),
    ("taken_null_alike_hi", 4, "95% of maps alike to the model take less than this"),
    ("p_taken_alike", 6, "share of those maps that take as much, one-sided"),
)


def numbers_table(table, members, tests, pairs):
    """The numbers of this step that the text quotes, one row each."""
    q = ISH_ANALYSIS["q"]
    top = table[table["top"]]
    free = top[top["in_model"] != "abundance"]
    follow = free[free["leftover_p"] < 0.05]
    take = free[free["p_taken"] < 0.05]
    take_alike = free[free["p_taken_alike"] < 0.05]
    untested = members[~members["tested"]]
    rows = [
        ("characterised", len(table), "genes characterised"),
        ("top_genes", len(top), f"genes past the map's null, BH over all, q < {q}"),
        ("top_list", " ".join(top["symbol"]), "those genes, best first"),
        (
            "top_within",
            int((top["q_within"] < q).sum()),
            f"of them, past the null inside divisions, q < {q}",
        ),
        (
            "top_follow_leftover",
            len(follow),
            "of them, not a term of the model, rho with the leftover at p < 0.05",
        ),
        ("top_follow_leftover_list", " ".join(follow["symbol"]), "those genes"),
        (
            "top_take",
            len(take),
            "of them, not a term of the model, taking more than 95% of its "
            "plain surrogates",
        ),
        ("top_take_list", " ".join(take["symbol"]), "those genes"),
        (
            "top_take_alike",
            len(take_alike),
            "of them, not a term of the model, taking more than 95% of maps alike to it",
        ),
        ("top_take_alike_list", " ".join(take_alike["symbol"]), "those genes"),
        (
            "added_surrogates",
            config.SETTINGS["top_genes"]["n_added_surrogates"],
            "maps of a gene's smoothness added in its place",
        ),
    ]
    by_gene = table.set_index("symbol")
    for gene in top_genes.NAMED_GENES:
        for column, digits, what in GENE_NUMBERS:
            if column not in by_gene:
                continue
            value = by_gene.loc[gene, column]
            # adding 0.0 turns a rounded -0.0 into 0.0
            value = int(value) if digits == 0 else round(float(value), digits) + 0.0
            rows.append((f"{column}_{gene}", value, f"{gene}: {what}"))
    rows += [
        ("family_listed", len(members), "members of the AMPA receptor complex family"),
        ("family_tested", int(members["tested"].sum()), "of them tested"),
        ("family_not_tested", " ".join(untested["symbol"]), "not tested"),
        ("family_controls", len(set(pairs.values())), "matched postsynaptic controls"),
    ]
    for map_name in ("leftover", "nano"):
        mine = tests[(tests["map"] == map_name) & (tests["tier"] == "2")]
        spatial = mine[mine["test"].str.contains("surrogates")].iloc[0]
        matched = mine[mine["test"].str.contains("controls")].iloc[0]
        key = f"family_{map_name}"
        rows += [
            (f"{key}_median", round(spatial["first"], 3), "the family's median rho"),
            (f"{key}_null_median", round(spatial["second"], 3), "over the surrogates"),
            (f"{key}_null_lo", round(spatial["null_lo"], 3), "2.5% of that null"),
            (f"{key}_null_hi", round(spatial["null_hi"], 3), "97.5% of that null"),
            (f"{key}_p", round(spatial["p"], 6), "the family's spatial p"),
            (
                f"{key}_controls_median",
                round(matched["second"], 3),
                "the matched controls' median rho",
            ),
            (
                f"{key}_controls_difference",
                round(matched["difference"], 3),
                "family minus controls, medians",
            ),
            (
                f"{key}_controls_critical",
                round(matched["critical"], 3),
                "the difference that would give p < 0.05",
            ),
            (f"{key}_controls_p", round(matched["p"], 6), "label-permutation p"),
        ]
        column = "leftover_q_family" if map_name == "leftover" else "q_family"
        past = table[table["family"] & (table[column] < q)]
        rows += [
            (
                f"{key}_pass_within",
                len(past),
                f"members past BH within the family, q < {q}",
            ),
            (f"{key}_pass_within_list", " ".join(past["symbol"]), "those members"),
        ]
    return pd.DataFrame(rows, columns=["name", "value", "what"], dtype=object)


def control_columns(pairs, leftover, nano, level):
    """Per family member, its matched control, both levels and the control's rho."""
    left_rho = leftover.set_index("symbol")["rho"]
    map_rho = nano.set_index("symbol")["rho"]
    rows = [
        dict(
            symbol=gene,
            median_energy=level.get(gene, np.nan),
            matched_control=control,
            control_median_energy=level.get(control, np.nan),
            control_rho=map_rho[control],
            control_leftover_rho=left_rho[control],
        )
        for gene, control in pairs.items()
    ]
    return pd.DataFrame(rows)


def main(sheets):
    """Print the settings; characterise the genes, run the named tests, draw."""
    config.print_settings({"sheets": sheets})
    figures = OUT / "figures"
    q = ISH_ANALYSIS["q"]

    # the tables of analyses 1 to 4, and the main model with its leftover's null
    ranking = gene_ranking.load_ranking()
    nano = ranking[ranking["map"] == "nano"]
    within, _ = divisions.load_within()
    robust = pd.read_csv(robustness.ROBUSTNESS)
    leftover = pd.read_csv(bd.LEFTOVER_GENES, keep_default_na=False, na_values=[""])
    leftover["in_model"] = leftover["in_model"].fillna("")
    inputs = bd.load_inputs()
    with np.load(bd.LEFTOVER_NULL) as z:
        null_structures = [str(s) for s in z["structures"]]
        leftover_null = z["null_rho"]
        leftover_genes = [str(g) for g in z["genes"]]
    if null_structures != inputs.structures:
        raise ValueError(
            f"{bd.LEFTOVER_NULL} was drawn on other structures than the main model's "
            f"{len(inputs.structures)}; run run_beyond_density.py first"
        )
    print(
        f"main model: {len(inputs.structures)} structures, density "
        f"{', '.join(inputs.terms['density'])}; {len(leftover)} genes against its "
        "leftover"
    )

    # the genes, and the family with its matched controls
    table = gene_table.load_gene_table()
    genes = gene_table.per_gene(table)
    chosen = top_genes.chosen_genes(nano, leftover, q)
    members = top_genes.family_members(genes, leftover)
    members.to_csv(top_genes.FAMILY_TABLE, index=False)
    family = list(members.loc[members["tested"], "symbol"])
    _, level = gene_table.gene_levels(table)
    in_leftover = set(leftover["symbol"])
    pool = [
        g
        for g in gene_sets.members_from_table(genes)["other postsynaptic"]
        if g in in_leftover
    ]
    pairs = top_genes.family_controls(family, pool, level)
    print(
        f"genes: {int(chosen['top'].sum())} past the map's null after BH, "
        f"{len(family)} of the family's {len(members)} members tested, each matched "
        f"to one of {len(pool)} other postsynaptic genes"
    )

    # what kind of map each gene is
    symbols = list(chosen["symbol"])
    documentation = pd.read_csv(
        gene_table.DOCUMENTATION, encoding="utf-8-sig", keep_default_na=False
    )
    records = gene_table.mygene_records(symbols, offline=True)
    parents, names, release = gene_table.load_obo(offline=True)
    annotation = top_genes.annotation_table(
        symbols, documentation, records, parents, names
    )
    ranges = top_genes.robustness_ranges(robust, symbols)
    rhos = top_genes.predictor_rhos(inputs, symbols)
    print(f"annotation: GO terms named from go-basic.obo, release {release}")

    # what each gene takes of the leftover, against maps of its smoothness
    centroids = structures.load_centroids()
    added_rows, added_nulls = [], {}
    for symbol in symbols:
        row, added_nulls[symbol] = top_genes.added_share(inputs, symbol, centroids)
        added_rows.append(dict(symbol=symbol, **row))
        print(
            f"  {symbol:8s} on {row['n_added']} structures takes "
            f"{row['taken']:+.3f} of the reproducible map; 95% of maps of its "
            f"smoothness take less than {row['taken_null_hi']:+.3f} (p "
            f"{row['p_taken']:.3f}), of maps alike to the model "
            f"{row['taken_null_alike_hi']:+.3f} (p {row['p_taken_alike']:.3f})",
            flush=True,
        )
    added = pd.DataFrame(added_rows)

    # tier 1, and tier 2 as a group on the leftover and, as a description, the map
    nano_null, nano_genes = gene_ranking.load_null_rho("nano")
    maps = dict(
        leftover=(leftover.set_index("symbol")["rho"], leftover_null, leftover_genes),
        nano=(nano.set_index("symbol")["rho"], nano_null, nano_genes),
    )
    group, group_nulls = top_genes.group_tests(family, pairs, maps)
    tests = pd.concat(
        [pd.DataFrame([top_genes.tier_one(leftover)]), group], ignore_index=True
    )
    tests.to_csv(top_genes.NAMED_TESTS, index=False)
    for r in tests.itertuples():
        print(f"tier {r.tier}, {r.map}: {r.test}: {r.first:+.3f}, p {r.p:.4f} ({r.role})")

    # tier 2 gene by gene, BH within the family
    p_leftover = leftover.set_index("symbol")["p_spatial"].reindex(family)
    p_map = nano.set_index("symbol")["p_spatial"].reindex(family)
    family_q = pd.DataFrame(
        dict(
            symbol=family,
            leftover_q_family=top_genes.within_family_q(p_leftover).to_numpy(),
            q_family=top_genes.within_family_q(p_map).to_numpy(),
        )
    )

    # the table
    sources = members[["symbol", "sources"]].rename(columns={"sources": "family_sources"})
    out = top_genes.characterise(
        [
            chosen,
            annotation,
            sources,
            top_genes.ranking_columns(nano),
            top_genes.within_columns(within),
            ranges,
            rhos,
            top_genes.leftover_columns(leftover),
            family_q,
            control_columns(pairs, leftover, nano, level),
            added,
        ]
    )
    for column in ("family_sources", "in_model", "matched_control"):
        out[column] = out[column].fillna("")
    out.to_csv(top_genes.TOP_TABLE, index=False)
    past = out[out["family"] & (out["leftover_q_family"] < q)]
    print(
        f"table: {len(out)} genes; family members past BH within the family against "
        f"the leftover: {len(past)}"
    )

    # the numbers for the text
    numbers = numbers_table(out, members, tests, pairs)
    numbers.to_csv(top_genes.NUMBERS, index=False)

    # figure 12, the genes that follow the map
    n_surrogates = int(leftover_null.shape[1])
    fig = ish_plotting.plot_top_genes(
        out,
        q,
        ISH_FIGURES["t_max"],
        n_surrogates,
        save=figures / ish_plotting.figure_file("top_genes"),
    )
    plt.close(fig)

    # figure 13, Cacng8 against Gria1, with the leftover of the main model
    set_table = structures.load_structure_set()
    declared = structures.declared_structures()
    map_values = profiles.load_profile().loc[declared, "zref_nano"]
    covariates, _, _ = bd.covariates_for(inputs)
    y = bd.full_map(inputs.nano)
    residual = pd.Series(
        bd.residual(y, bd.model(covariates, inputs.terms)), index=inputs.structures
    )
    gap = pd.read_csv(gene_ranking.GAP)
    gap_null, _ = gene_ranking.load_null_rho("gap")
    n_map_surrogates = int(nano_null.shape[1])
    fig = ish_plotting.plot_cacng8_gria1(
        out.set_index("symbol"),
        map_values,
        inputs.expr,
        residual,
        gap[gap["kind"] == "merged profiles"].iloc[0],
        gap_null[0],
        set_table,
        n_surrogates,
        n_map_surrogates,
        save=figures / ish_plotting.figure_file("cacng8_gria1"),
    )
    plt.close(fig)

    # figure 14, the AMPA receptor complex against the leftover
    fig = ish_plotting.plot_ampa_family(
        out[out["family"]],
        tests,
        group_nulls,
        members,
        q,
        ISH_FIGURES["t_max"],
        n_surrogates,
        n_map_surrogates,
        save=figures / ish_plotting.figure_file("ampa_family"),
    )
    plt.close(fig)
    drawn = [ish_plotting.figure_file(k) for k in ("top_genes", "cacng8_gria1")]
    drawn.append(ish_plotting.figure_file("ampa_family"))
    print(f"figures: {', '.join(drawn)} in {figures}")

    # the gene sheets
    if not sheets:
        return
    plane = ISH_FIGURES["plane"]
    structure_names, _, _ = structure_terms()
    lab = planes.label_plane(plane)
    folder = figures / "top_genes"
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
            n_surrogates,
            n_map_surrogates,
            save=folder / f"{row['symbol']}.png",
        )
        plt.close(fig)
    print(f"gene sheets: {len(out)} in {folder}")


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
