# Shared code

Functions more than one pipeline uses: the cohort table, reading volumes,
left-right statistics, background pixels and the colours. Nothing here runs on
its own. The paths are at the code root: `get_paths` (where the code and the
data are, `SEP_DATA_ROOT` included) and `sep_setup_paths` (the MATLAB path,
once per session).

| file | what it does | used by |
|---|---|---|
| `cohort.csv` | the cohort table: one row per mouse, with its group, age in days, the subfolder of its raw files on the share, and its cohort and order in the Python route | `get_cohort`; `mapping/sepmap/config.py` |
| `get_cohort.m` | every mouse of the table, or the mice of some groups, or named mice, each with its folder under the data root; `get_cohort('verify')` checks the adults' order | every MATLAB driver that takes mice |
| `get_cohort_spec.m` | a cohort name (`'naive'`, `'young_P20'`) into its mice, its stack files and its atlas | `run_normalise_groups`, `collect_by_group`, `run_characterize_distribution` |
| `loadVolume.m` | read a multi-page TIFF, or a folder of TIFFs, into an H x W x n array | `run_residual_correction`, `run_register_to_atlas`, `run_collect_by_group` |
| `compute_lr_stats.m` | left minus mirrored right, and left plus mirrored right, of a volume | `run_group_differences`, `run_characterize_distribution`, `run_compare_with_allen_ish` |
| `select_background_pixels.m` | the background pixels of one slice image, at the knee of its percentile curve | `run_nano_equalisation`, `run_normalise_groups`, `select_reference_pixels` |
| `sep_palette.m` | the project's colours and colormaps, by name (what each is for: `docs/STYLE.md`, Figures) | the figures of every MATLAB pipeline |
| `get_color2color_colormap.m` | a diverging colormap from one colour through white to another | `sep_palette`, the difference maps, `run_characterize_distribution`, `run_compare_with_allen_ish` |

## The cohort table

- A mouse's name is the name of its raw folder on the share and of its
  folder under the data root (`<data>\<group>\<name>\`). A young mouse's age
  is the one in that name (MG904_SepGluA_P22 is P22); the adults' ages are not
  recorded and stay empty (NaN in MATLAB).
- The first 17 rows, the adults, keep their order. `run_normalise_groups`
  and `run_group_differences` select adults by position within their group
  (the behaviour mice of `run_group_differences` by name, `behavior_mice`),
  and `get_cohort('verify')` stops a run when the order has changed. New mice
  go at the end.
- `share_subdir` is used only by `run_copy_raw_data`: the local copy always
  puts the `.czi` files at the root of the mouse's folder.
- `mapping_cohort` and `mapping_order` are read only by the Python route: the
  cohort a brain enters there (empty: not taken) and the order it stacks its
  brains in, which is not the MATLAB order. Each route keeps the order its
  results were produced in: a cohort mean summed in another order changes in
  its last bits.
- `../tests/test_backward_compat` checks the adults' order and folders, that
  the young cohort `get_cohort` returns is the table's young rows in their
  order (the test reads the table, so a new row needs no change to it), and
  that each young brain's age matches its name.

## Notes

- `select_background_pixels` marks a slice as all background when it is
  nearly empty or its centre of mass lies outside the middle third of the
  image; with no knee found it falls back to the percentile `p_max`, with a
  warning.
- `compute_lr_stats` pairs plane i with plane ML + 1 - i; for an odd width the
  midline plane is left out (the registered grids are 1140 planes wide).
