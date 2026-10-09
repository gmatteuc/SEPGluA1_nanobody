"""The synapse-density genes of analysis 4, chosen by measured PSD95 without the map.

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
    22. run_density_markers    the synapse-density genes, chosen by   <- this script
                               PSD95 without the map (network,
                               once); figure 03s3
    23. run_beyond_density     analysis 4: what Gria1 and synapse
                               density leave; the leftover
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_calibration the same model on maps whose answer
                               is known; the check rows, each
                               with its own floor
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

Downloads the mouse GO annotation (MGI's GAF) and go-basic.obo of one GO release from
the GO Consortium's archive into reference/go/<release>/ under the data root (only
what is missing, each file checked against its SHA-256), then applies the rule of
[density_markers]: the pool of postsynaptic genes, the AMPA-linked genes left out,
the pool ranked by Spearman with the measured PSD95 punctum density, the highest
chosen, and the choice repeated on random halves of the structures. It reads the gene
table, the declared set and the synaptome's density table, and no nano value; the
model of analysis 4 takes the genes from settings.toml, and the run stops at its end
if the rule chose others. The method is in sepmap/adult/density_markers.py. Writes,
under the data root:

    reference/go/<release>/           the GAF, go-basic.obo and fetch_log.txt
                                      (release, date, hashes, each file's version)
    adult_v2/ish_analysis/density_markers/
        candidates.csv                every gene of the gene table against the rule:
                                      each criterion, its postsynaptic and plasticity
                                      terms, in the pool or why not
        excluded.csv                  every reason a gene is excluded: the rule, the
                                      term, the gene's own term below it, its evidence
        agreement.csv                 the eligible genes by Spearman with PSD95, the
                                      pool's ranks, the genes chosen
        structures.csv                the declared structures: Gria1, the genes chosen
                                      and PSD95 measured or not; in the model or why not
        halves.csv                    per random half, the genes chosen and each
                                      composite's agreement on the other half
        selection.csv                 how often each pool gene is chosen
        comparison.csv                each composite's agreement with PSD95, held out
                                      and on the full set
    adult_v2/ish_analysis/tables/numbers_density_markers.csv   the numbers for the text
    adult_v2/ish_analysis/figures/03s3_density_markers.png     the choice (and .eps)

    python run_density_markers.py [--offline]

--offline stops instead of downloading when a file is missing.
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from sepmap import config, plotting
from sepmap.adult import density_markers, synaptome
from sepmap.adult import plotting as adult_plotting
from sepmap.ish import gene_table
from sepmap.ish import plotting as ish_plotting
from sepmap.ish.figure_index import figure_path
from sepmap.structures import load_structure_set

# the number of genes chosen and of random halves
DENSITY_MARKERS = config.SETTINGS["density_markers"]


def main(offline):
    """Print the settings, fetch the annotation, apply the rule, write and draw."""
    config.print_settings({"offline": offline})
    density_markers.OUT.mkdir(parents=True, exist_ok=True)

    # the GO release's files, downloaded once and checked against their hashes
    fetched = density_markers.fetch(offline)
    releases = density_markers.versions()
    print(
        f"GO release {density_markers.GO_RELEASE}: {len(fetched)} of "
        f"{len(density_markers.FILES)} files downloaded into {density_markers.REFERENCE}"
    )
    for name, version in releases.items():
        print(f"  {name}: {version}")

    # the inputs: the gene table, the declared set and PSD95 density; no nano value
    genes = gene_table.per_gene(gene_table.load_gene_table())
    profiles = gene_table.load_profiles()
    set_table = load_structure_set()
    declared = sorted(set_table.loc[set_table["in_set"], "structure"])
    psd95 = synaptome.load_density()[synaptome.MEASURE]
    measured = [s for s in declared if np.isfinite(psd95.get(s, np.nan))]
    print(
        f"inputs: {len(genes)} genes, {len(declared)} declared structures, PSD95 "
        f"measured in {len(measured)} of them"
    )

    # the annotation of the genes, and the rule's terms checked against their names
    annotations, n_not = density_markers.read_gaf(
        density_markers.GAF, set(genes["symbol"])
    )
    density_markers.check_gaf(annotations)
    parents, names, _ = gene_table.read_obo(density_markers.OBO)
    density_markers.check_terms(density_markers.POOL_TERMS, names)
    density_markers.check_terms(density_markers.EXCLUSION_TERMS, names)
    terms = density_markers.annotated_terms(annotations)
    print(
        f"annotation: {len(annotations)} annotations of {len(terms)} genes, {n_not} with "
        f"NOT left out; the {len(density_markers.POOL_TERMS)} pool and "
        f"{len(density_markers.EXCLUSION_TERMS)} exclusion terms carry their names"
    )

    # the pool, and the genes excluded from it
    candidates, excluded = density_markers.candidate_table(
        genes, profiles, declared, terms, parents, names
    )
    candidates.to_csv(density_markers.CANDIDATES, index=False)
    excluded.to_csv(density_markers.EXCLUDED, index=False)
    eligible = candidates[candidates["eligible"]]
    out = eligible[eligible["excluded"]]
    pool = list(candidates.loc[candidates["in_pool"], "symbol"])
    print(
        f"pool: {len(eligible)} genes meet the rule, {len(out)} of them excluded as "
        f"AMPA-linked, {len(pool)} in the pool"
    )
    for r in out.itertuples():
        print(f"  excluded {r.symbol}: {r.exclude_rules}")
    plastic = candidates[candidates["in_pool"] & (candidates["plasticity_terms"] != "")]
    print(
        f"  {len(plastic)} pool genes carry a plasticity term, kept: "
        + " ".join(plastic["symbol"])
    )

    # the choice: the pool ranked by agreement with PSD95
    agreement = density_markers.agreement_table(candidates, profiles, psd95, measured)
    agreement.to_csv(density_markers.AGREEMENT, index=False)
    chosen = density_markers.chosen_genes(agreement)
    picked = agreement[agreement["chosen"]]
    print(
        "chosen: "
        + ", ".join(f"{r.symbol} {r.rho:+.3f}" for r in picked.itertuples())
        + f" (Spearman with PSD95 over {int(picked['n_structures'].min())} structures)"
    )
    if len(chosen) < DENSITY_MARKERS["n_markers"]:
        print(
            f"  warning: only {len(chosen)} genes in the pool, fewer than the "
            f"{DENSITY_MARKERS['n_markers']} asked; the term takes them all"
        )

    # the declared structures the model can run on
    structures = density_markers.structure_table(set_table, profiles, chosen, psd95)
    structures.to_csv(density_markers.STRUCTURES, index=False)
    n_model = int(structures["in_model"].sum())
    print(f"structures: {n_model} of {len(declared)} declared have Gria1 and the genes")
    for reason, part in structures[~structures["in_model"]].groupby("reason"):
        divisions = part["division"].value_counts().to_dict()
        print(f"  {len(part):3d}  {reason}: {divisions}")

    # the choice repeated on random halves, and the composites compared
    fixed, members = density_markers.fixed_composites(genes, profiles, measured)
    halves = density_markers.half_split(pool, profiles, psd95, measured, fixed)
    selection = density_markers.selection_table(halves, pool)
    comparison = density_markers.comparison_table(
        chosen, profiles, psd95, measured, halves, fixed, members
    )
    halves.to_csv(density_markers.HALVES, index=False)
    selection.to_csv(density_markers.SELECTION, index=False)
    comparison.to_csv(density_markers.COMPARISON, index=False)
    top = selection.head(DENSITY_MARKERS["n_markers"] + 3)
    print(
        f"halves: {len(halves)}; chosen most often "
        + ", ".join(f"{r.symbol} {r.share:.0%}" for r in top.itertuples())
    )
    for r in comparison.itertuples():
        print(
            f"  {r.composite:15s} held out {r.rho_held_out_median:+.3f} "
            f"({r.rho_held_out_lo:+.3f} to {r.rho_held_out_hi:+.3f}), full set "
            f"{r.rho_full_set:+.3f} on {r.n_structures}"
        )

    # the numbers for the text
    counts = dict(
        declared=len(declared),
        gria1=int(structures["gria1_measured"].sum()),
        psd95=len(measured),
        model=n_model,
    )
    numbers = density_markers.numbers_table(
        candidates, excluded, agreement, selection, comparison, counts, releases
    )
    numbers.to_csv(density_markers.NUMBERS, index=False)

    # figure 03s3, with the genes left out grouped by their first reason
    reasons = density_markers.exclusion_reasons(excluded)
    fig = adult_plotting.plot_density_markers(
        agreement,
        density_markers.composite_on(chosen, profiles, measured),
        psd95.reindex(measured),
        ish_plotting.group_of(set_table),
        selection,
        comparison,
        reasons,
        DENSITY_MARKERS["n_halves"],
        density_markers.GO_RELEASE,
        save=figure_path("density_markers"),
    )
    plt.close(fig)
    print(
        f"figure: {figure_path('density_markers')}; {len(reasons)} reasons for the "
        f"{int(reasons['n_genes'].sum())} genes left out"
    )

    # the genes chosen against those the model reads
    density_markers.check_chosen(chosen)
    print("check passed: the genes chosen are those of [density_markers] chosen")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    parser = argparse.ArgumentParser(
        description="the synapse-density genes, chosen by measured PSD95 without the map"
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="stop instead of downloading a missing file",
    )
    args = parser.parse_args()
    main(args.offline)
