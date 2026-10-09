"""The declared structure set of the ISH line, and the adult profiles it reads.

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
    13. run_structure_set      A1: the declared structures, the     <- this script
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
    27. run_ish_top_genes      the genes that follow the map most,
                               characterised; Cacng8 and the AMPA
                               receptor complex family against the
                               leftover; figures 07, 08, 11, 11s2,
                               sheets
    28. run_sep_channel_check  analysis 5: what the green channel
                               reports; figures 15, 15s
    29. run_ish_overview       figures 00, 16 and 16s, the figure
                               index; the numbers for the text

Measures the three channels of every adult per structure, plain and eroded by one
20 um voxel, declares the structures every comparison of the ISH line uses (grey
matter measured in all ten adults, settings [structures]), and takes zref with that
set as reference (the methods are in sepmap/structures.py and
sepmap/adult/profiles.py). Writes, in adult_v2/ish_analysis/ under the data root:

    tables/adult_per_mouse.csv       per adult and structure: voxels, the channel
                                     means plain and eroded (the cache), cref, zref
    tables/structure_set.csv         every structure of the adult table, in the
                                     set or not, and why
    tables/zref_reference.csv        per adult and channel, the zero and spread of
                                     zref, and the stored ones of region_plot
    tables/adult_profile.csv         per structure, the mean over the adults
    tables/centroids.csv             the declared structures' centroids, one
                                     hemisphere, in mm
    tables/numbers_structure_set.csv the numbers of this step, for the text
    figures/01_structures.png        which structures, and why (and .eps)

    python run_structure_set.py [--recompute]

--recompute measures the channel means again from the per-mouse files (about ten
minutes) instead of reading them from adult_per_mouse.csv.
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from sepmap import config, plotting, structures
from sepmap.adult import profiles
from sepmap.ish import planes
from sepmap.ish import plotting as ish_plotting
from sepmap.ish.figure_index import figure_path
from sepmap.volumes.per_mouse import annotation_20, structure_terms

ISH_FIGURES = config.SETTINGS["ish_figures"]
STRUCTURES = config.SETTINGS["structures"]


def main(recompute):
    """Print the settings, then declare the set, take zref and write the tables."""
    config.print_settings({"recompute": recompute})
    structures.TABLES.mkdir(parents=True, exist_ok=True)
    structures.FIGURES.mkdir(parents=True, exist_ok=True)

    # the CCF as structure codes, its interior and the ontology
    names, acronyms, divisions = structure_terms()
    ann = annotation_20("ccf")
    labels, label_names = structures.name_volume(ann, names)
    meta = profiles.structure_meta(names, acronyms, divisions)

    # the channel means of every adult, plain and eroded, measured or from the cache
    if recompute or not profiles.PER_MOUSE_TABLE.exists():
        print("channel means of the ten adults, plain and eroded:", flush=True)
        interior = ~structures.border_voxels(labels) & (labels > 0)
        means = profiles.channel_table(labels, interior, label_names, meta)
        source = "measured"
    else:
        means = profiles.load_channel_table()
        source = f"read from {profiles.PER_MOUSE_TABLE.name}"
    n_plain = int(means["eroded_is_plain"].sum())
    print(
        f"channel means: {means['mouse'].nunique()} adults, {len(means)} rows, "
        f"{n_plain} eroded below the smallest size and kept plain ({source})"
    )

    # the declared set
    set_table = structures.structure_set(means, profiles.ADULTS)
    set_table.to_csv(structures.STRUCTURE_SET, index=False)
    declared = sorted(set_table.loc[set_table["in_set"], "structure"])
    by_division = set_table[set_table["in_set"]]["division"].value_counts().to_dict()
    print(f"structure set: {len(declared)} declared of {len(set_table)} ({by_division})")

    # zref with the declared reference, per adult, and each adult's reference
    per_mouse, references = profiles.with_zref(means, declared)
    per_mouse.to_csv(profiles.PER_MOUSE_TABLE, index=False)
    reference = profiles.reference_table(references)
    reference.to_csv(profiles.REFERENCE, index=False)
    nano = reference[reference["channel"] == "nano"]
    shift = (nano["median"] - nano["median_stored"]).abs()
    print(f"zref reference: zero moved by {shift.min():.3f} to {shift.max():.3f} log2")

    # the adult profile
    profile = profiles.adult_profile(per_mouse, set_table)
    profile.to_csv(profiles.PROFILE, index=False)
    profile = profiles.load_profile()
    agreement = profiles.stored_agreement(profile)
    print(
        f"adult profile: {len(profile)} structures, {int(profile['in_set'].sum())} "
        f"in the set; order against the stored zref rho {agreement:.5f}"
    )

    # centroids of the declared structures, one hemisphere
    centroids = structures.centroids(labels, label_names, declared)
    centroids.to_csv(structures.CENTROIDS, index=False)
    n_missing = int(centroids["ap_mm"].isna().sum())
    print(
        f"centroids: {len(centroids) - n_missing} of {len(declared)}, ML up to "
        f"{centroids['ml_mm'].max():.2f} mm"
    )

    # the numbers for the text
    numbers = profiles.numbers_table(set_table, reference, centroids, agreement)
    numbers.to_csv(profiles.NUMBERS, index=False)

    # figure 01
    plane = ISH_FIGURES["plane"]
    fig = ish_plotting.plot_structures(
        set_table,
        profile,
        reference,
        agreement,
        np.asarray(ann[planes.crop_index(plane)]),
        names,
        plane,
        STRUCTURES["min_adults"],
        save=figure_path("structures"),
    )
    plt.close(fig)
    print(f"figure: {figure_path('structures')}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")
    plotting.set_style()

    parser = argparse.ArgumentParser(description="the declared structure set (A1)")
    parser.add_argument(
        "--recompute",
        action="store_true",
        help="measure the channel means again (about ten minutes)",
    )
    args = parser.parse_args()
    main(args.recompute)
