"""The adult map by cortical depth: flatmaps, every area by band, laminar profiles.

Python route, in run order (tools\\venv_atlas; run_closeup and run_adult_layers in
tools\\venv_flat):
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
    11. run_ish_regions        ISH per structure, 100-gene panel
    12. run_ish_compare        the adult map against each gene
    13. run_ish_words          annotation words of the ranking
    14. run_ish_roles          subunit against localisation genes
    15. run_adult_arms         the channel arms per adult
    16. run_ish_arms           the arms against the genes
    17. run_sep_channel_check  what the green channel reports
    18. run_panel_build        the 390-gene ontology panel
    19. run_panel_fetch        its ISH grids (network, once)
    20. run_ish_regions        ISH per structure, ontology panel
    21. run_ish_reliability    how reliable one ISH map is
    22. run_ish_panel_test     localisation against controls
    23. run_beyond_density     what abundance and density leave
    24. run_beyond_controls    seven attempts to break it
    25. run_beyond_figures     the figures of that result
    26. run_beyond_regression  the regression, shown
    27. run_adult_layers       the adult map by depth (venv_flat)     <- this script

Every isocortex area of the ten adults, measured per mouse in three depth bands
and five layers, ordered along the cortical hierarchy of Harris et al. 2019, beside
the flatmaps of the young against adult close-up drawn for the adults alone. It
reads the per-mouse files of steps 1 and 2 and checks itself against step 5's
region_means_per_mouse.csv; the method is in sepmap/adult/layers.py. The flatmaps
need ccf_streamlines, so this runs in tools\\venv_flat; with --no-flatmap the
tables and figures 02 and 03 also run in tools\\venv_atlas. Writes, in
adult_v2/layers/ under the data root:

    area_layers_per_mouse.csv        one row per adult, area and depth, kept or
                                     missing with the reason
    area_layers_summary.csv          per area and depth: n, mean, SD, SEM, t
    depth_summary.csv                per depth: reliability, rho with the hierarchy
    01_flatmaps_by_band_<smooth>.png mean and between-mouse SD per depth band
    02_areas_by_band.png             every area per band, a dot per adult
    03_laminar_profiles.png          eleven areas, five layers, a line per adult
    flatmap_bands.npz                the cached band maps of every adult

and an .eps beside each PNG.

    python run_adult_layers.py [--no-flatmap] [--reproject]

--no-flatmap leaves the flatmaps out; --reproject projects every adult again
instead of reading the cached band maps.
"""

import argparse

import matplotlib
import matplotlib.pyplot as plt

from sepmap import config
from sepmap.adult import layers, layers_plotting
from sepmap.volumes import cohort
from sepmap.young_vs_adult import region_groups, region_plot

ADULT_LAYERS = config.SETTINGS["adult_layers"]
CLOSEUP = config.SETTINGS["closeup"]
YOUNG_VS_ADULT = config.SETTINGS["young_vs_adult"]


def smoothing_tag(sigma):
    """The smoothing as a file name carries it: _smooth3x1x1, or _nosmooth."""
    if not any(sigma):
        return "_nosmooth"
    return "_smooth" + "x".join(f"{s:g}" for s in sigma)


def main(want_flatmap, reproject):
    """Print the settings in force, then measure, write and draw the map by depth."""
    # settings in force
    config.print_settings(
        {
            "reading": "zref",
            "adults": len(layers.ADULTS),
            "hierarchy": ADULT_LAYERS["hierarchy"],
            "flatmap": "yes" if want_flatmap else "no",
            "reproject": "yes" if reproject else "no",
        }
    )
    layers.OUT.mkdir(parents=True, exist_ok=True)

    # the labels of every isocortex area at every depth, and the hierarchy
    stru, divi, sub = region_groups.load_parcellation_terms()
    groups = layers.depth_groups(stru, divi, sub)
    hier = layers.load_hierarchy()
    layers.check_hierarchy(hier, {area for _, _, area in groups})
    atlas_n = layers.atlas_counts(groups)
    undrawn = [f"{area} {depth}" for (_, depth, area), n in atlas_n.items() if n == 0]
    n_scored = int(hier["hierarchy_score"].notna().sum())
    print(
        f"{len(hier)} isocortex areas, {n_scored} scored by Harris 2019; "
        f"{len(groups)} area-depth groups, left out with no voxel in the atlas: "
        f"{', '.join(undrawn) or 'none'}"
    )

    # every brain of the region tables, for the zref median and spread; the adults' cells
    print("measuring every brain on its own atlas")
    cells, n_vox, refs, struct_mean = layers.measure(groups, stru, divi)
    norm = region_plot.range_match(layers.ROUTE_MICE, struct_mean, refs)

    # the per-mouse table, checked against the region table
    table = layers.per_mouse_table(groups, cells, n_vox, norm, refs, hier, atlas_n)
    table = layers.with_contrast(table)
    n_checked = layers.check_against_region_table(table)
    print(f"whole-area zref equals region_means_per_mouse.csv in all {n_checked} cells")

    # group statistics, and per depth the reliability and the hierarchy rho
    summary = layers.summary_table(table)
    depths = layers.depth_summary(table, summary)
    table.round({"coverage": 4, "zref": 4}).to_csv(
        layers.OUT / "area_layers_per_mouse.csv", index=False
    )
    summary.round({"mean": 4, "sd": 4, "sem": 4, "t": 3}).to_csv(
        layers.OUT / "area_layers_summary.csv", index=False
    )
    depths.round(4).to_csv(layers.OUT / "depth_summary.csv", index=False)
    missing = table[table["excluded"] & (table["depth_kind"] != "whole")]
    print(
        f"{len(table)} cells written, {len(missing)} of them missing (whole-area cells "
        f"aside), in {missing['area'].nunique()} areas: "
        + ", ".join(f"{a} {n}" for a, n in missing.groupby("area").size().items())
    )
    for _, r in depths.iterrows():
        print(
            f"  {r['depth']:30s} half rho {r['half_rho']:.3f}"
            f"  all ten {r['whole_rho']:.3f}  ({r['n_areas_complete']} areas)"
            f"   hierarchy rho {r['hierarchy_rho']:+.2f}  p {r['hierarchy_p']:.3f}"
            f"  ({r['n_areas_scored']} areas)"
        )

    # figures 02 and 03
    n_adults = len(layers.ADULTS)
    fig = layers_plotting.plot_area_bands(
        table, summary, depths, hier, n_adults, save=layers.OUT / "02_areas_by_band.png"
    )
    plt.close(fig)
    fig = layers_plotting.plot_laminar_profiles(
        table, summary, hier, n_adults, save=layers.OUT / "03_laminar_profiles.png"
    )
    plt.close(fig)
    print("wrote 02_areas_by_band.png and 03_laminar_profiles.png")
    if not want_flatmap:
        return

    # each adult's band maps, their mean and SD per pixel, and figure 01
    sigma = list(CLOSEUP["smooth"])
    maps = layers.band_maps(sigma, recompute=reproject)
    min_n = YOUNG_VS_ADULT["min_n_adult"]
    mean, sd, n = layers.band_stats(maps, min_n)
    name = f"01_flatmaps_by_band{smoothing_tag(sigma)}.png"
    fig = layers_plotting.plot_band_flatmaps(
        mean, sd, n_adults, min_n, sigma, save=layers.OUT / name
    )
    plt.close(fig)
    n_pixels = [int((nb >= min_n).sum()) for nb in n]
    print(f"wrote {name}; pixels with at least {min_n} adults per band: {n_pixels}")


if __name__ == "__main__":
    # figures go to files, never to a window
    matplotlib.use("Agg")

    parser = argparse.ArgumentParser(description="the adult map by cortical depth")
    parser.add_argument("--no-flatmap", action="store_true", help="no flatmaps")
    parser.add_argument(
        "--reproject",
        action="store_true",
        help="project every adult again instead of reading the cached band maps",
    )
    args = parser.parse_args()

    # the reading this step draws must be among those in force
    if "zref" not in cohort.MODES:
        parser.error("this step reads zref, which V2_READINGS leaves out")
    main(want_flatmap=not args.no_flatmap, reproject=args.reproject)
