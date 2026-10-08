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
                               documentation; figure 02
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the
                               map, the null, autofluorescence (A8),
                               the Cacng8 - Gria1 gap; figures 05 to 07 and 12
    18. run_ish_robustness     A3: the ranking under other choices;
                               figure 13
    19. run_ish_divisions      analysis 2 (A6): between or within
                               divisions; figure 10, gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation
                               against matched controls; figures 08, 09
    21. run_synaptome          the measured synapse density (network,   <- this script
                               once): PSD95 puncta per structure,
                               its coverage of the fit
    22. run_beyond_density     analysis 4: what receptor mRNA and
                               synaptic density leave; the leftover
    23. run_beyond_controls    seven attempts to break it
    24. run_beyond_calibration the same model on maps whose answer
                               is known
    25. run_beyond_regression  the regression, per structure
    26. run_beyond_figures     figures 03, 04 and 11
    27. run_sep_channel_check  analysis 5: what the green channel
                               reports; figure 14
    28. run_ish_overview       figures 00 and 15, the figure index;
                               the numbers for the text

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
        feasibility.png               the coverage and the agreement (and .eps)
    adult_v2/ish_analysis/tables/numbers_synaptome.csv   the numbers for the text

    python run_synaptome.py [--offline]

--offline stops instead of downloading when a file is missing.
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

from sepmap import config, plotting
from sepmap.adult import beyond_density, profiles, synaptome
from sepmap.ish import gene_table
from sepmap.ish import plotting as ish_plotting
from sepmap.structures import load_structure_set
from sepmap.volumes.per_mouse import annotation_20


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

    # the densities and the samples, checked
    density, info = synaptome.load_source()
    samples = synaptome.sample_table(density, info)
    print(
        f"source read: {density.shape[0]} subtypes x {density.shape[1]} samples, every "
        f"subtype scaled 0..1; {int((~samples['measured']).sum())} sample with no "
        "punctum at all"
    )

    # each sample placed in a structure through the CCF ontology
    ontology = synaptome.load_ontology()
    structures = synaptome.structure_ids()
    samples = synaptome.match_samples(samples, ontology, structures)
    samples.to_csv(synaptome.SAMPLES, index=False)
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
    per_structure = synaptome.structure_density(samples, voxels)
    set_table = load_structure_set()
    inputs = beyond_density.load_inputs()
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
    n_fit, n_fit_measured = len(inputs.structures), int(coverage["fit_measured"].sum())
    n_set, n_set_measured = (
        int(coverage["declared"].sum()),
        int(coverage["declared_measured"].sum()),
    )
    replaces = synaptome.in_main_model(n_fit_measured, n_fit)
    print(
        f"coverage: fit {n_fit_measured} of {n_fit} ({n_fit_measured / n_fit:.0%}), "
        f"declared {n_set_measured} of {n_set}; the rule asks "
        f"{synaptome.BEYOND['min_psd95_coverage']:.0%}, so PSD95 "
        + ("replaces the mRNA density terms" if replaces else "is a variant only")
    )

    # agreement with the terms of the fit, and on the declared set
    covariates, _, _ = beyond_density.build_covariates(
        inputs.expr, inputs.role, inputs.auto, inputs.structures
    )
    fit_terms = {
        name: pd.Series(covariates[key], index=inputs.structures)
        for name, key in (
            ("markers", "markers"),
            ("psd_pc1", "psd_pc1"),
            ("Gria1", "Gria1"),
            ("autofluorescence", "autofluo"),
        )
    }
    fit_terms["nano"] = pd.Series(inputs.nano.mean(axis=0), index=inputs.structures)
    profile = profiles.load_profile()
    declared = sorted(set_table.loc[set_table["in_set"], "structure"])
    set_terms = {
        "Gria1": pd.Series(gene_table.load_profiles()["Gria1"]),
        "autofluorescence": profile["zref_auto"],
        "nano": profile["zref_nano"],
    }
    rows = synaptome.agreement_rows(table, fit_terms, inputs.structures, "fit")
    rows += synaptome.measure_rows(table, inputs.structures, "fit")
    rows += synaptome.agreement_rows(table, set_terms, declared, "declared")
    agreement = pd.DataFrame(rows)
    agreement.to_csv(synaptome.AGREEMENT, index=False)
    hemispheres = synaptome.hemisphere_agreement(table, inputs.structures)
    psd95 = agreement[
        (agreement["density"] == synaptome.MEASURE) & (agreement["structures"] == "fit")
    ]
    print(
        "agreement of psd95 on the fit (Spearman): "
        + ", ".join(f"{r.term} {r.rho:+.2f}" for r in psd95.itertuples())
        + f"; left against right {hemispheres[1]:+.2f} (n = {hemispheres[0]})"
    )

    # the numbers for the text
    numbers = synaptome.numbers_table(samples, table, agreement, hemispheres)
    numbers.to_csv(synaptome.NUMBERS, index=False)

    # the diagnostic figure
    path = synaptome.OUT / "feasibility.png"
    fig = ish_plotting.plot_synaptome(
        table,
        coverage,
        agreement,
        fit_terms["markers"],
        synaptome.BEYOND["min_psd95_coverage"],
        save=path,
    )
    plt.close(fig)
    print(f"figure: {path}")


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
