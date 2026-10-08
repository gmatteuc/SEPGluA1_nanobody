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
    11. run_panel_build        the 390-gene ontology panel (network, cached)
    12. run_panel_fetch        its ISH grids (network, once)
    13. run_structure_set      A1: the declared structures, the adult    <- this script
                               profiles; figure 01
    14. run_ish_section_qc     A2: section QC of every experiment; QC sheets
    15. run_ish_gene_table     A9: region means, the gene table, merged profiles,
                               gene sets, documentation; the panel repair
                               (network); figure 02
    16. run_ish_spatial_null   A7: surrogate maps and their checks
    17. run_ish_gene_ranking   analysis 1: each gene against the map, the null,
                               autofluorescence (A8), the Cacng8 - Gria1 gap;
                               figures 03 to 06
    18. run_ish_robustness     A3: the ranking under other choices
    19. run_ish_divisions      analysis 2 (A6): between or within divisions;
                               per-gene sheets
    20. run_ish_gene_sets      analysis 3: kinds of genes; localisation against
                               matched controls
    21. run_beyond_density     analysis 4: what abundance and density leave
    22. run_beyond_controls    seven attempts to break it
    23. run_beyond_calibration the same model on maps whose answer is known
    24. run_beyond_regression  the regression, per structure
    25. run_beyond_figures     figures 11 and 12
    26. run_sep_channel_check  analysis 5: what the green channel reports
    27. run_ish_overview       figures 00 and 14; the numbers for the text

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
import pandas as pd
from scipy.stats import spearmanr

from sepmap import config, plotting, structures
from sepmap.adult import profiles
from sepmap.ish import plotting as ish_plotting
from sepmap.volumes.per_mouse import annotation_20, structure_terms

ISH_FIGURES = config.SETTINGS["ish_figures"]
STRUCTURES = config.SETTINGS["structures"]

OUT = config.DATA / "adult_v2" / "ish_analysis"

# the first 10 um CCF plane of the adult crop (volumes.per_mouse.atlas_grid)
CROP_START = 180


def numbers_table(set_table, reference, centroids, agreement):
    """The numbers of this step that the text quotes, one row each."""
    nano = reference[reference["channel"] == "nano"]
    shift = (nano["median"] - nano["median_stored"]).abs()
    rows = [
        ("structures_in_table", len(set_table), "structures in the adult table"),
        (
            "structures_all_adults",
            int((set_table["n_adults"] >= STRUCTURES["min_adults"]).sum()),
            "measured in all the adults the rule asks for",
        ),
        ("structures_declared", int(set_table["in_set"].sum()), "the declared set"),
        (
            "zref_zero_shift_min",
            round(shift.min(), 4),
            "smallest shift of an adult's zero",
        ),
        (
            "zref_zero_shift_max",
            round(shift.max(), 4),
            "largest shift of an adult's zero",
        ),
        (
            "cohort_map_agreement",
            round(agreement, 5),
            "Spearman of the cohort map before and after A1, declared set",
        ),
        (
            "centroids_ml_max_mm",
            round(float(centroids["ml_mm"].max()), 3),
            "largest ML centroid, mm (one hemisphere: below 5.7)",
        ),
    ]
    in_set = set_table[set_table["in_set"]]
    for division, n in in_set["division"].value_counts().items():
        rows.append(
            (f"declared_{division}", int(n), f"declared structures in {division}")
        )
    return pd.DataFrame(rows, columns=["name", "value", "what"])


def stored_agreement(profile):
    """Spearman of the 10 adults' mean stored zref with the new one, over the set."""
    table = pd.read_csv(profiles.REGION_MEANS)
    stored = table[(table["reading"] == "zref") & table["mouse"].isin(profiles.ADULTS)]
    stored = stored.groupby("structure")["log2_value"].mean()
    in_set = profile[profile["in_set"]]
    return float(spearmanr(stored.reindex(in_set.index), in_set["zref_nano"]).statistic)


def main(recompute):
    """Print the settings, then declare the set, take zref and write the tables."""
    config.print_settings({"recompute": recompute})
    tables = OUT / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    # the CCF as structure codes, its interior and the ontology
    names, acro, divi = structure_terms()
    ann = annotation_20("ccf")
    labels, label_names = structures.name_volume(ann, names)
    meta = profiles.structure_meta(names, acro, divi)

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
    agreement = stored_agreement(profile)
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
    numbers = numbers_table(set_table, reference, centroids, agreement)
    numbers.to_csv(tables / "numbers_structure_set.csv", index=False)

    # figure 01
    figures = OUT / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    plane = ISH_FIGURES["plane"]
    lab = ann[(plane - CROP_START) // 2]
    fig = ish_plotting.plot_structures(
        set_table,
        profile,
        reference,
        agreement,
        np.asarray(lab),
        names,
        plane,
        STRUCTURES["min_adults"],
        save=figures / "01_structures.png",
    )
    plt.close(fig)
    print(f"figure: {figures / '01_structures.png'}")


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
