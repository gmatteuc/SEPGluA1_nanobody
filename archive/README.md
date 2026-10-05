# Archive

Retired code, kept until Giulio has checked what replaced it (decision L2 of
the refactor). Nothing outside this folder calls into it, and
`sep_setup_paths` does not put it on the MATLAB path. The tag `refactor-start`
keeps every file, so deleting one loses nothing.

| file | what it was | replaced by | can go |
|---|---|---|---|
| `P6_analyze_group_averages_and_normalize.m`, `P7_analyze_group_differences.m` | the first normalisation and group-difference scripts of the plasticity comparison | P6bis and P7bis, now `group_comparison/run_normalise_groups.m` and `run_group_differences.m` | with the rest of the folder |
| `plot_abs_slice.m`, `plot_diff_slice.m`, `write_diff_video.m` | the slice figures and videos of P6 and P7 | the figures and videos of `run_group_differences` | with P6 and P7 |
| `lr_sum_and_diff.m` | left-right sum and difference of a TIFF stack, not called | `common/compute_lr_stats.m` | with the rest of the folder |
| `scratch.m` | working fragments of the plasticity comparison, run by hand in a workspace that held its stacks: mask smoothing, normalisation checks, region quantification, left-right videos; not called | `group_comparison/` (`run_normalise_groups`, `run_group_differences`) | with the rest of the folder |
| `compare_young_vs_adult_lrsum.py`, `replot_young_vs_adult_lrsum.py`, `region_means_raw_per_mouse.py`, `region_ratio_young_vs_adult.py`, `plot_region_ratio_young_vs_adult.py` | the first young-against-adult comparison (September 2026, three P20 brains), whose outputs are in `<data>\comparisons\young_P20_vs_adult_nano\` | the Python route, `mapping/` (`sepmap/volumes/` and `sepmap/young_vs_adult/`) | with the rest of the folder |
| `extractAxioscanImages.m` | an image extractor that came with BioformatsImage, not called | `preprocessing/run_extract_and_center` | with the rest of the folder |

The rest of the folder goes once Giulio has checked the Python route against
the old code, run from the tag `refactor-start` in a check tree of its own
([ROADMAP.md](../docs/ROADMAP.md), section 4).

Deleted on 5 October, and still at the tag `refactor-start`:

- `landmark_refine`, the image matcher behind the control-point GUI's old `r`
  key, retired in L3: its Python package and LoFTR weights, its MATLAB
  wrappers, its setup script and its requirements;
- `bk/`, old backup copies of early helpers. Its one note still in use, the
  elastix install, is now in the [README](../README.md#setup).

## Do not run

- The five `.py` scripts write to `D:\sep_histology\data\comparisons\young_P20_vs_adult_nano\`,
  a path written into each of them, with none of the guards of `get_paths.m`
  and `mapping/sepmap/config.py`: run from a copy of the code, they would
  still write into the production data. `compare_young_vs_adult_lrsum.py`
  creates that folder as soon as it is imported.
- `extractAxioscanImages.m` writes an `extracted_cy3\` folder next to the
  `.czi` files it reads: run on a raw folder, it writes beside the raw data.
