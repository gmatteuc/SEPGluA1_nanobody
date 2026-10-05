# Adult analyses in MATLAB

The earlier analyses of the adult map: how the signal of a channel is
distributed across the pooled adult brain, the nano map against the Allen in
situ hybridisation (ISH) maps, and nano against autofluorescence structure by
structure. They keep running until the Python route (`../mapping/`) answers
the same questions with its additions A1 to A5, then move to `archive/`
([`../docs/ROADMAP.md`](../docs/ROADMAP.md), sections 3 and 4).

## Run order

| step | script | what it does | replaced by |
|---|---|---|---|
| 1 | `run_characterize_distribution` | pool the normalised volumes of some cohorts into one; per voxel, the cohort mean of the hemisphere sum (left plus mirrored right), its reliability t (mean / SEM) and its z-score over the brain; per structure (the finest structures of nine divisions), a distance-weighted or eroded mean, the enrichment call, the mean of each mouse and the SEM; videos and bar charts | A1, A4 |
| 2 | `run_compare_with_allen_ish` | for each gene of `<data>\gene_targets.csv`, the Allen ISH grid with its failed sections repaired, resized onto the CCF and folded like the map; Spearman and Pearson correlations with the map across structures, a voxel Pearson; the genes ranked | A1 to A3 |
| 3 | `run_compare_nano_with_autofluorescence` | per structure and per division, a one-sided paired signed-rank test of nano above autofluorescence across mice, on the means and on the contrast of the z-scored channels | A5 |

Step 1 reads the normalised stacks that steps 1 and 2 of
[`../group_comparison/`](../group_comparison/README.md) write
(`<group>\<channel>_4d_normalized.mat`). It runs once per channel, with
`channel = 'nano'`, then `'auto'` for step 3; steps 2 and 3 read its outputs,
in either order. Run `sep_setup_paths` first, once per MATLAB session; the
header of each script gives its settings, inputs and outputs.

The runs behind the outputs of record:

1. `run_characterize_distribution` with `groups_to_merge = {'naive', 'rws'}`,
   `apply_smoothing = false`, `region_agg_method = 'distweight'` and
   `compute_per_mouse_sem = true`, once with `channel = 'nano'` (the bar
   charts of April and May 2026, rerun on 4 September 2026) and once with
   `'auto'`;
2. `run_compare_with_allen_ish` with `channel = 'nano'`, on the 100-gene
   panel (21 and 22 April 2026);
3. `run_compare_nano_with_autofluorescence` with `agg_method = 'distweight'`,
   `apply_bonferroni = false` and `alpha = 0.05` (May 2026).

`SEP_MERGE_SPECS` in the environment replaces the cohorts of step 1
(`young_P20` for the P20 brains on their own atlas, which wrote
`merged_young_P20_nano\`, superseded by the Python route).

## Outputs

Under `<data>\comparisons\`, each figure as `.fig` and `.png`:

- `merged_<cohorts>_<channel>\` (step 1): the cohort map
  `mean_lr_sum_*.mat`, the region table `Region_MeanSum_Table_*.csv` and
  `.mat`, the per-mouse means `per_mouse_region_means_*.mat`, the region
  masks `roi_masks_*.mat` and label positions `label_centroids_*.mat`; the
  bar charts by division and across divisions, raw and z-scored, with and
  without `_withSEM`; the videos `lr_sum_<channel>_*`, `lr_sum_threshold_*`
  and `lr_sum_zscore_*`; the distance-weight sheets
  `Diagnostic_DistWeight_*`
- `merged_naive_rws_<channel>_vs_ish_<gene>_nosmooth\` (step 2, one per gene):
  `gene_result_<gene>.mat`, `ish_lr_sum_<gene>.mat`,
  `Region_NanoVsISH_Table_<gene>.csv`, the scatter, the paired bars, the
  comparison video and two diagnostic sheets; and
  `merged_naive_rws_<channel>_vs_ish_summary_nosmooth\`:
  `gene_panel_summary.csv` and `.mat`, the correlation bars and violins
- `nano_vs_auto\` (step 3): the p values per structure and per division
  (`NanoVsAuto_*_signrank_*.csv`), the paired and delta-z bars
  (`Region_NanoVsAuto_*`, `Macro_NanoVsAuto_*`) and the video of the voxelwise
  contrast (`Contrast_video_NanoMinusAuto_z_*.mp4`)

Step 2 downloads each Allen grid once, into `<data>\atlas_ish\`.

The outputs of record stay where they are, frozen, as what the additions are
compared with: `merged_naive_rws_nano\`, `merged_naive_rws_auto\`,
`merged_naive_rws_vs_ish_*` and `nano_vs_auto\`. Step 2's folders have no
channel in their names: they predate it, and a rerun writes
`merged_naive_rws_nano_vs_ish_*` beside them. `../mapping/run_ish_compare.py`
reads their `gene_panel_summary.csv`. Which figure is where:
[`../docs/FIGURES.md`](../docs/FIGURES.md).

## Where the code is

Each script holds all its code, its local functions at the end; unlike the
other MATLAB pipelines, they are not split into a driver and a function in
`pipeline/`, and they keep the style they were written in, since they
retire. Shared: `../common/` (`get_cohort_spec`, `compute_lr_stats`,
`get_color2color_colormap`), `../atlas/` (`get_atlas_crop`,
`get_allen_region_mask`), `write_lr_video` in `../group_comparison/pipeline/`;
`plot_violinplot`, in this folder, draws step 2's violins.

## Notes

- These outputs carry known defects, which the additions are not to
  reproduce ([`../docs/ROADMAP.md`](../docs/ROADMAP.md), section 4): the slab
  artefact of the `abs()` in step 1's adult outputs; step 2's ISH grid
  stretched onto the atlas box (a 1.5 to 2.5% scale error), and its section
  repair, which "repaired" a true absence of expression (Slc17a6) and copied
  the last good section over the posterior planes (Cacng8); the
  autofluorescence of steps 1 and 3, from an `auto_4d.mat` of the
  registration before the one behind the nano stack, which
  `run_collect_by_group` no longer writes; step 3's Bonferroni, which cannot
  reject at n = 10; step 1's reliability bars, which start at white. They
  are not quoted.
- Zero in the z-scored bars and the enrichment call is the brain's voxel
  mean; A1 and A4 move it to the brain's median structure (S2 of
  `../docs/ROADMAP.md`).
- What the Python route already covers of the 88 outputs and tests of these
  scripts, item by item: [`../docs/REFACTOR_COVERAGE.md`](../docs/REFACTOR_COVERAGE.md).
- Steps 1 and 2 put the CCF folder on the path themselves (`addpath`), which
  elsewhere only `get_atlas` does: start a registration in a fresh session
  after them.
- `force_recompute_masks` in step 1 recomputes the label positions, the
  region table and the per-mouse means; the region masks, once built, are
  reused.
