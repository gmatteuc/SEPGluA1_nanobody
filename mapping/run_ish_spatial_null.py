"""The spatial null of the ISH line: surrogate maps with the adult map's smoothness.

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
    16. run_ish_spatial_null   A7: surrogate maps and their checks   <- this script
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

Draws variogram-matched surrogates (Burt 2020) of the adult nano map and of the
autofluorescence map on the declared structures, and checks them: the variogram of
the surrogates against the map's, and the false positives of the spatial p on
independent random smooth maps (the method is in sepmap/ish/spatial_null.py).
Writes tables only; figure 06 is drawn by run_ish_gene_ranking, whose last panel
needs the genes. In adult_v2/ish_analysis/tables/ under the data root:

    surrogates_nano.npy, surrogates_auto.npy  surrogates x structures, as ranks
    surrogate_structures.csv                  the structures of their columns
    variogram.csv                             per map and distance: the map's
                                              variogram, the surrogates' and a
                                              shuffled map's
    null_calibration.csv                      one row per test of an independent
                                              random map: rho, ordinary and
                                              spatial p
    numbers_spatial_null.csv                  the numbers of this step, for the text

    python run_ish_spatial_null.py [--recompute]

--recompute draws the surrogates again instead of reading them (a few seconds
each map; the calibration is always run, a few minutes).
"""

import argparse

import pandas as pd

from sepmap import config, structures
from sepmap.adult import profiles
from sepmap.ish import spatial_null

# the maps, their profile columns, and the seed of each one's surrogates
MAPS = {"nano": ("zref_nano", 0), "auto": ("zref_auto", 1)}


def main(recompute):
    """Print the settings, draw or read the surrogates, check them, write the tables."""
    config.print_settings({"recompute": recompute})

    # the declared structures, their centroids and the two maps
    declared = structures.declared_structures()
    profile = profiles.load_profile().loc[declared]
    centroids = structures.load_centroids().loc[declared]
    d = spatial_null.distance_matrix(centroids)
    print(f"maps: {len(declared)} declared structures, distances up to {d.max():.2f} mm")

    # the surrogates of each map, drawn or read
    surr = {}
    if recompute or not all(p.exists() for p in spatial_null.SURROGATES.values()):
        for name, (column, seed) in MAPS.items():
            surr[name] = spatial_null.surrogates(profile[column].to_numpy(), d, seed=seed)
        spatial_null.save_surrogates(surr, centroids)
        source = "drawn"
    else:
        for name in MAPS:
            surr[name], _ = spatial_null.load_surrogates(name, declared)
        source = "read"
    shapes = ", ".join(f"{k} {v.shape[0]} x {v.shape[1]}" for k, v in surr.items())
    print(f"surrogates: {shapes} ({source})")

    # the variogram of each map, its surrogates' and a shuffled map's
    variogram = pd.concat(
        [
            spatial_null.variogram_table(name, profile[column].to_numpy(), d, surr[name])
            for name, (column, _) in MAPS.items()
        ],
        ignore_index=True,
    )
    variogram.to_csv(spatial_null.VARIOGRAM, index=False)
    misfit = spatial_null.matched_misfit(variogram)
    text = ", ".join(
        f"{k} median {m:.0%}, largest {x:.0%}" for k, (m, x) in misfit.items()
    )
    print(f"variogram: surrogates against the map inside the matched range: {text}")

    # the calibration: independent random fields with the map's smoothness
    nano = profile["zref_nano"].to_numpy()
    params = spatial_null.fit_exponential(nano, d)
    calibration = spatial_null.calibration_table(d, params, nano, surr["nano"])
    calibration.to_csv(spatial_null.CALIBRATION, index=False)
    rates = spatial_null.false_positive_rates(calibration)
    for _, r in rates.iterrows():
        print(
            f"calibration ({r['design']}, {r['n_tests']} tests): p < "
            f"{spatial_null.ALPHA} in "
            f"{r['ordinary']:.1%} by the ordinary p, {r['spatial']:.1%} by the spatial p"
        )

    # the numbers for the text
    numbers = spatial_null.numbers_table(len(declared), params, misfit, rates)
    numbers.to_csv(spatial_null.NUMBERS, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="the spatial null of the ISH line (A7)")
    parser.add_argument(
        "--recompute", action="store_true", help="draw the surrogates again"
    )
    args = parser.parse_args()
    main(args.recompute)
