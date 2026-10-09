"""The measured synapse density: PSD95 puncta per structure, and how much of the fit.

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
    21. run_synaptome          the measured synapse density (network,   <- this script
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

Downloads the PSD95 and SAP102 punctum densities of Zhu et al. 2018, as Hansen et al.
share them, from their repository at a pinned commit into reference/synaptome/ under
the data root (only what is missing, each file checked against the commit's hash),
places each of the 775 samples in a structure of the adult table through the CCF
ontology, and takes per structure the PSD95 density analysis 4 uses, with three
variants. Then counts how many structures of the declared set and of analysis 4's
fit it covers, which decides whether it replaces the mRNA density terms in the main
model (beyond.min_psd95_coverage), and how it agrees with those terms, with Gria1,
with the nano map and across the two hemispheres. The method is in
sepmap/adult/synaptome.py. Writes, under the data root:

    reference/synaptome/              the source's files and fetch_log.txt (source,
                                      commit, citation, hashes, date)
    adult_v2/ish_analysis/synaptome/
        samples.csv                   per sample: name, hemisphere, ids, densities,
                                      its unit and structure, or why it has none
        density.csv                   per structure of the adult table: in the
                                      declared set and in the fit, measured or why
                                      not, its units, weights and densities
        coverage.csv                  per division, the declared and fitted
                                      structures and how many are measured
        agreement.csv                 Spearman of each density with the mRNA density
                                      terms, Gria1, the nano and autofluorescence maps
    adult_v2/ish_analysis/tables/numbers_synaptome.csv   the numbers for the text
    adult_v2/ish_analysis/figures/14s1_synaptome_detail.png
                                      the coverage, the two hemispheres and the
                                      agreement (and .eps); figure 14, with the check
                                      rows of analysis 4, is drawn by step 27

    python run_synaptome.py [--offline]

--offline stops instead of downloading when a file is missing.
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from sepmap import config, plotting
from sepmap.adult import beyond_density, profiles, synaptome
from sepmap.adult import plotting as adult_plotting
from sepmap.ish import gene_table
from sepmap.ish.figure_index import figure_path
from sepmap.structures import load_structure_set
from sepmap.volumes.per_mouse import annotation_20

# the share of the fit the measured density must cover to enter the main model
BEYOND = config.SETTINGS["beyond"]
ISH = config.SETTINGS["ish"]


def agreement_terms(inputs, set_table):
    """The maps each density is compared with, on the fit and on the declared set.

    On the structures of the fit, the density terms, Gria1, autofluorescence and
    the map as the fit builds them; on the declared set, Gria1's profile and the
    two maps. Returns both and the declared structures.
    """
    covariates, _, _ = beyond_density.build_covariates(
        inputs.expr, inputs.role, inputs.auto, inputs.structures
    )
    gene = ISH["control_gene"]
    fit_terms = {
        name: pd.Series(covariates[key], index=inputs.structures)
        for name, key in (
            ("markers", "markers"),
            ("psd_pc1", "psd_pc1"),
            (gene, gene),
            ("autofluorescence", "autofluo"),
        )
    }
    fit_terms["nano"] = pd.Series(inputs.nano.mean(axis=0), index=inputs.structures)
    profile = profiles.load_profile()
    set_terms = {
        gene: pd.Series(gene_table.load_profiles()[gene]),
        "autofluorescence": profile["zref_auto"],
        "nano": profile["zref_nano"],
    }
    declared = sorted(set_table.loc[set_table["in_set"], "structure"])
    return fit_terms, set_terms, declared


def main(offline):
    """Print the settings, then fetch, place and measure, and write the tables."""
    config.print_settings({"offline": offline})
    synaptome.OUT.mkdir(parents=True, exist_ok=True)

    # the source's files, downloaded once and checked against the commit
    fetched = synaptome.fetch(offline)
    print(
        f"source: {len(synaptome.FILES)} files of {synaptome.REPOSITORY} at "
        f"{synaptome.COMMIT[:7]} in {synaptome.REFERENCE}, {len(fetched)} downloaded"
    )

    # each sample placed in a structure through the CCF ontology
    samples, ontology, structures, shape = synaptome.placed_samples()
    samples.to_csv(synaptome.SAMPLES, index=False)
    print(
        f"source read: {shape[0]} subtypes x {shape[1]} samples, each subtype divided "
        f"by its largest value; {int((~samples['measured']).sum())} sample with no "
        "punctum at all"
    )
    by = samples["matched_by"].value_counts().to_dict()
    print(
        f"samples placed: {int(samples['used'].sum())} of {len(samples)} in "
        f"{samples.loc[samples['used'], 'structure'].nunique()} structures "
        f"({by.get('id', 0)} with the source's Allen id, {by.get('acronym', 0)} found "
        "by acronym)"
    )
    left_out = samples.loc[~samples["used"], "reason"].str.split(": ", n=1, expand=True)
    for kind, part in left_out.groupby(0, sort=False):
        print(f"  {len(part):3d}  {kind}: " + "; ".join(sorted(set(part[1]))))

    # the density per structure, weighted by the units' voxels in the 20 um CCF
    parent = ontology["parent"].to_dict()
    voxels = synaptome.id_voxels(annotation_20("ccf"), parent)
    per_structure = synaptome.structure_density(samples, voxels, parent)
    set_table = load_structure_set()
    inputs = beyond_density.load_inputs(measured=False)
    table = synaptome.density_table(
        set_table, inputs.structures, per_structure, samples, ontology, structures
    )
    table.to_csv(synaptome.DENSITY, index=False)
    n_equal = int((per_structure["weighting"] == "equal").sum())
    print(
        f"density: {len(per_structure)} structures with a sample, "
        f"{int(table['measured'].sum())} of them in the adult table; {n_equal} with "
        "units weighted alike (a unit not drawn in the CCF)"
    )

    # coverage of the declared set and of the fit, and the rule
    coverage = synaptome.coverage_table(table)
    coverage.to_csv(synaptome.COVERAGE, index=False)
    n_fit = len(inputs.structures)
    n_fit_measured = int(coverage["fit_measured"].sum())
    n_set = int(coverage["declared"].sum())
    n_set_measured = int(coverage["declared_measured"].sum())
    if synaptome.in_main_model(n_fit_measured, n_fit):
        verdict = "replaces the mRNA density terms"
    else:
        verdict = "is a variant only"
    print(
        f"coverage: fit {n_fit_measured} of {n_fit} ({n_fit_measured / n_fit:.0%}), "
        f"declared {n_set_measured} of {n_set}; the rule asks "
        f"{BEYOND['min_psd95_coverage']:.0%}, so PSD95 {verdict}"
    )

    # agreement with the terms of the fit, and on the declared set
    fit_terms, set_terms, declared = agreement_terms(inputs, set_table)
    agreement = synaptome.agreement_table(
        table, fit_terms, set_terms, inputs.structures, declared
    )
    agreement.to_csv(synaptome.AGREEMENT, index=False)
    hemispheres = synaptome.hemisphere_agreement(table, inputs.structures)
    measured = agreement[
        (agreement["density"] == synaptome.MEASURE) & (agreement["structures"] == "fit")
    ]
    print(
        f"agreement of {synaptome.MEASURE} on the fit (Spearman): "
        + ", ".join(f"{r.term} {r.rho:+.2f}" for r in measured.itertuples())
        + f"; left against right {hemispheres[1]:+.2f} (n = {hemispheres[0]})"
    )

    # the numbers for the text
    numbers = synaptome.numbers_table(samples, table, agreement, hemispheres)
    numbers.to_csv(synaptome.NUMBERS, index=False)

    # figure 14s1; figure 14 needs analysis 4's check rows and is drawn at step 27
    fig = adult_plotting.plot_synaptome_detail(
        table,
        coverage,
        agreement,
        fit_terms["markers"],
        BEYOND["min_psd95_coverage"],
        save=figure_path("synaptome_detail"),
    )
    plt.close(fig)
    print(f"figure: {figure_path('synaptome_detail')}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    parser = argparse.ArgumentParser(
        description="the measured synapse density: PSD95 puncta per structure"
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="stop instead of downloading a missing file",
    )
    args = parser.parse_args()
    main(args.offline)
