# Production settings for the reference run

Written on 30 Sep 2026 (step 0.6 of [REFACTOR_PLAN.md](REFACTOR_PLAN.md)) for
the reference run of step 3: the old code, at the tag `refactor-start`, run in
the check trees on G:. It lists, driver by driver, the settings that run uses.
Every later run of the new code uses the same values, so old and new are
compared under identical settings. The cohort tree ran on the plan's small
reference set (Verification design), two mice per adult group, not on the
whole groups the plasticity tables below assume; everything else is as listed.

Names are those since the refactor, with the old name in parentheses (the full
table: [refactor_name_map.csv](../refactor_name_map.csv)). Three things changed
how a setting is given:

- the constants of the Python route are keys of `mapping/settings.toml`, named
  below as `table.key`;
- its two ISH passes are chosen with `--panel targets|ontology`, and
  `V2_ISH_PANEL` and `V2_ISH_TABLE` are refused;
- `v2_paths.py` is `mapping/sepmap/config.py`.

How to read the tables:

- committed: the value in commit `04c0484`, the last one before step 0. Step 0
  changes no setting value. It adds `SEP_DATA_ROOT` and the data-root guards,
  the guard in the register driver's `align` mode, and lets `V2_ISH_PANEL` be
  a path relative to the data root.
- reference: the value the reference run uses, and why.
- hidden state: what a run depends on besides its settings (caches reused
  when present, files read from earlier runs, globs, index-based selections).
- Points that needed a decision are marked TO CONFIRM and collected at the
  end, as drafted on 30 Sep.

## Common to every run

| item | reference |
|---|---|
| code | cohort tree: worktree at `refactor-start` in `G:\sep_refactor\ref\code`; registration checks: `G:\sep_refactor\reg\old\code` |
| MATLAB | R2024b, a fresh session per driver, `restoredefaultpath` then `sep_setup_paths`, started by `tools\run_matlab_detached.ps1` with `-CodeDir`, `-DataRoot` and `-LogDir` |
| toolboxes | Image Processing, Statistics, Parallel Computing (`run_nano_equalisation` uses `parfor`) |
| elastix | 5.1.0, installed per user and on the PATH (`run_register_to_atlas` mode `register`, `run_add_sep_channel`) |
| settings with no environment override | set by editing the settings block in the worktree copy, the edit saved as a patch beside the outputs. TO CONFIRM |

### Environment variables

The data root is the folder `data` beside the code folder
(`D:\sep_histology\data` on the lab computer), or the folder `SEP_DATA_ROOT`
names. A copy of the code is refused the production data.

| variable | read by | committed default (unset) | reference |
|---|---|---|---|
| `SEP_DATA_ROOT` | `get_paths.m`, `mapping/sepmap/config.py`, `atlas/build_demba_atlas.py` (from step 0) | `<code>\..\data` | the check tree's data folder, in every run |
| `SEP_COHORT_SPECS` | `run_normalise_groups` (P6bis) | `{'rws','naive'}` | `rws,naive,behavior` (S6) |
| `P4BIS_MICE` | `run_add_sep_channel` (P4bis) | groups young, naive, rws | `MG914_SepGluA_P28`, then `CGF027_Gria1` (registration tree only) |
| `SEP_MERGE_SPECS` | `run_characterize_distribution` (P8) | its own list | unset: it is not run |
| `V2_READINGS` | `sepmap/volumes/cohort.py`, and through it every step that draws or tabulates | `readings.in_force`, all five readings | unset |
| `V2_ISH_PANEL` | old `v2_ish_regions` only; refused since step 5 | `data\gene_targets.csv` (100 genes) | pass 1 unset; pass 2 `adult_v2\panel\panel_v2.csv` (relative, inside the tree's data) |
| `V2_ISH_TABLE` | old `v2_ish_regions`, `v2_ish_reliability` only; refused since step 5 | `gene_region_table.csv` in `v2_ish_regions`, `gene_region_table_panel.csv` in `v2_ish_reliability` | pass 1 unset; pass 2 `gene_region_table_panel.csv` |
| `AUTO_ANNOTATION_PYTHON` | `registration/pipeline/auto_annotate.m` | the engine's `.venv` (`registration\auto_annotation\.venv` after the merge) | `D:\sep_histology\code\auto_annotation\.venv\Scripts\python.exe` (a worktree has no `.venv`) |

In the old code both ISH variables had to be unset for every script outside
pass 2: with `V2_ISH_TABLE` left set, pass 1's `v2_ish_regions` would
overwrite the 390-gene table. The new code takes the pass from `--panel`
(`[ish_panels]` in `settings.toml`) and stops if either variable is set.

### Interpreters

| environment | path | Python | used by |
|---|---|---|---|
| analysis | `D:\sep_histology\code\tools\venv_atlas` | 3.12.7 | every `mapping/run_*.py` except `run_closeup`; holds the DeMBA-to-CCF deformation fields `run_to_ccf` needs (downloaded on first use if missing) |
| flatmaps | `D:\sep_histology\code\tools\venv_flat` | 3.12.7 | `run_closeup` only (ccf_streamlines) |
| automatic annotation | `D:\sep_histology\code\auto_annotation\.venv` (`registration\auto_annotation\.venv` after the merge) | 3.12.7 | `run_register_to_atlas` mode `autoannotate` (torch, GPU) |

The environments are ignored by git, so the worktrees have none. The
reference run calls `main`'s interpreters on the worktree's scripts
(`<python> <tree>\code\v2_x.py` for the old code, `mapping\run_x.py` for the
new), with `SEP_DATA_ROOT` set. TO CONFIRM for the two analysis environments
(the plan states it for the engine only).

## Run order

### Cohort tree

1. `run_collect_by_group` (P5) for rws, naive, behavior.
2. `run_normalise_groups` (P6bis) for rws, naive, behavior.
3. `run_group_differences` (P7bis) naive against rws, then naive against
   behavior.
4. The Python route, in the order of the headers of `mapping/run_*.py`:

| # | run script (old module) | stage |
|---|---|---|
| 1 | `run_per_mouse` (`v2_per_mouse`) | volumes |
| 2 | `run_to_ccf` (`v2_to_ccf`) | volumes |
| 3 | `run_cohort` (`v2_cohort`) | volumes |
| 4 | `run_compare` (`v2_compare`) | young against adult |
| 5 | `run_region_plot` (`v2_region_plot`) | young against adult |
| 6 | `run_region_groups` (`v2_region_groups`) | young against adult |
| 7 | `run_video` (`v2_video`) | young against adult |
| 8 | `run_video_compare` (`v2_video_compare`) | young against adult |
| 9 | `run_closeup` (`v2_inspect`; venv_flat) | young against adult |
| 10 | `run_diagnostics` (`v2_diagnostics`) | diagnostics |
| 11 | `run_ish_regions --panel targets` (`v2_ish_regions`) | ISH pass 1 (100 genes) |
| 12 | `run_ish_compare` (`v2_ish_compare`) | ISH pass 1 |
| 13 | `run_ish_words` (`v2_ish_words`) | ISH pass 1 |
| 14 | `run_ish_roles` (`v2_ish_roles`) | ISH pass 1 |
| 15 | `run_adult_arms` (`v2_adult_arms`) | adult arms |
| 16 | `run_ish_arms` (`v2_ish_arms`) | adult arms |
| 17 | `run_sep_channel_check` (`v2_sep_channel_check`) | adult arms |
| 18 | `run_panel_build` (`v2_panel_build`) | ISH pass 2 (390 genes) |
| 19 | `run_panel_fetch` (`v2_panel_fetch`) | ISH pass 2 |
| 20 | `run_ish_regions --panel ontology` (`v2_ish_regions` with both ISH variables set) | ISH pass 2 |
| 21 | `run_ish_reliability --panel ontology` (`v2_ish_reliability`) | ISH pass 2 |
| 22 | `run_ish_panel_test` (`v2_ish_panel_test`) | ISH pass 2 |
| 23 | `run_beyond_density` (`v2_beyond_density`) | beyond abundance |
| 24 | `run_beyond_controls` (`v2_beyond_controls`) | beyond abundance |
| 25 | `run_beyond_figures` (`v2_beyond_figures`) | beyond abundance |
| 26 | `run_beyond_regression` (`v2_beyond_regression`) | beyond abundance |

Not run: `run_replot` (`v2_replot`; it only redraws `run_compare`'s figures
from the volumes `run_compare` saved).

What fixes the order:

- `run_diagnostics` sheet 08 reads the background masks of
  `run_normalise_groups`, and sheet 07 reads `region_table.csv`
  (`run_compare`) and `region_stats.csv` (`run_region_plot`).
- `young_vs_adult\region_means_per_mouse.csv` (`run_region_plot`) is read by
  `run_ish_compare`, `run_ish_roles`, `run_ish_panel_test`, `run_adult_arms`
  (a consistency check) and the four `run_beyond_*`.
- `adult_v2\ish\gene_region_table.csv` (pass 1) is read by `run_ish_compare`,
  `run_ish_roles`, `run_ish_arms` and `run_sep_channel_check`.
- `run_ish_words` reads `gene_correlations.csv` (`run_ish_compare`);
  `run_ish_arms` reads `arms\region_means_arms.csv` (`run_adult_arms`).
- `run_beyond_density` and `run_beyond_controls` read
  `gene_region_table_merged.csv` (pass 2, `run_ish_reliability`);
  `run_beyond_figures` draws panel D only if `controls.csv`
  (`run_beyond_controls`) exists.

Within these constraints the order is free; the production run of 24 to 26
Sep 2026 differs only in running `v2_ish_roles` after the arms.

### Registration tree

Inputs are restored from `reg\old` before each stage, so no stage reads
another's outputs.

1. MG914_SepGluA_P28: `run_extract_and_center` (P1),
   `run_residual_correction` (P2), `run_nano_equalisation` (P2bis).
2. MG914_SepGluA_P28 and CGF027_Gria1: `run_register_to_atlas` (P4) mode
   `register`, then `run_add_sep_channel` (P4bis) on the same mouse.
3. MG914_SepGluA_P28: `run_register_to_atlas` mode `autoannotate`.

The restore matters most between stages 1 and 2: `run_extract_and_center`
rewrites `sliceinfo.mat`, and in the old code `register` and
`run_add_sep_channel` read the `px_atlas` inside it to size the registered
volume (see below; since fix 26 of step 8 they set the 10 um grid
themselves).

## Plasticity comparison (MATLAB)

### run_collect_by_group.m (P5_collect_data_by_group.m)

| setting | committed | reference | why |
|---|---|---|---|
| `mousetypes_list` | `{'young'}` | `{'rws','naive','behavior'}` | the three adult groups `run_normalise_groups` and `run_group_differences` read (S6) |
| `age_filter` | `[20]` | `[]` | the adults have no recorded age; `[]` writes plain `nano_4d.mat` |
| `skip_missing` | `true` | `true` | all 17 adults are registered (26 Nov 2025): the whole groups give 5, 5 and 7 mice |

The old script also had `correction_type`, which it never used, and an inline
`min_num_contrib = 3` feeding averages it no longer saved; step 8 removed both
with the rest of its dead code (`nano_4d.mat` unchanged).

- Reads: `<group>\<mouse>\lightsuite\volume_registered\chan02_NANO.tiff` of
  every mouse in the group (the old script also read `chan03_AUTO.tiff` and
  `chan05_MASK.tiff` for the averages it did not save).
- Writes: `<group>\nano_4d.mat`; `<group>\collected_mice.mat` (new for the
  adults, whose stacks predate that file).
- Hidden state: only `nano_4d` is saved, so the `auto_4d.mat`, `mask_4d.mat`
  and averages of 19 Nov 2025 stay on disk. They are not copied into the
  tree, and nothing in the reference run reads them. A mouse whose tiff is
  missing is skipped with only a printed note, and the index-based selection
  of `run_normalise_groups` would then pick the wrong mice.

### run_normalise_groups.m (P6bis_analyze_group_averages_and_normalize.m)

| setting | committed | reference | why |
|---|---|---|---|
| `mice`, `mousetypes` | the 17 adults, registry order | unchanged | index-based selection depends on the order |
| `mousetypes_list` | `{'rws','naive','behavior'}` | unchanged | |
| `selected_mice_idx_list` | rws `1:5`, naive `1:5`, behavior `[1,3,4,5]` | unchanged | behavior = MG705, MG709, MG716, MG718 (S6) |
| `plot_verification_video` | `false` | unchanged | |
| `channel` | `'nano'` | unchanged | |
| `cohort_specs` | `{'rws','naive'}` | `SEP_COHORT_SPECS=rws,naive,behavior` | the naive against behavior comparison needs behavior normalised by the same code (S6) |

Fixed in code, unchanged: pooled planes `1:5:900` (times `ap_scale`, 1 for
the adults); diagnostic planes 450, 500, 550; `plot_limit` 5000; 150
histogram bins; cortex mask Isocortex layers 1 to 5 minus retrosplenial,
anterior cingulate and prelimbic; the background percentile schedule by AP
plane (six bands, from 0/95 to 20/70); robust `fitlm` of each mouse against
the pooled median (polyfit if it fails, identity below 100 pixels);
`rng(42)` for the plotted subsample only; raw zeros set to NaN after the fit.

- Reads: `<group>\nano_4d.mat` (this run's `run_collect_by_group`);
  `data\atlas` through `get_atlas_crop('ccf')`.
- Writes: `<group>\nano_4d_normalized.mat`,
  `<group>\nano_4d_normalized_bkgmask.mat`,
  `<group>\Background_mask_diagnostics_trace_<group>_nano.fig/.png`,
  `<group>\global_diagnostics\normalization_checks_nano\`.
- Hidden state: an adult group is selected by position in the fourth
  dimension of `nano_4d.mat`, which assumes `run_collect_by_group` stacked the
  whole group in registry order; `collected_mice.mat` is not consulted for
  the adults.

### run_group_differences.m (P7bis_analyze_group_differences.m)

Run twice, with only `exp_type` changed.

| setting | committed | reference | why |
|---|---|---|---|
| `ctrl_type` | `'naive'` | `'naive'` | |
| `exp_type` | `'behavior'` | run 1 `'rws'`, run 2 `'behavior'` | both comparisons Sami approved (S6) |
| `selected_mice_idx_list` | as `run_normalise_groups` | unchanged | used only for the mouse names in figures; the data are whatever `run_normalise_groups` saved |
| `behavior_subset` (inline in the old script; `behavior_mice` since fix 11 of step 8) | `[1,2,3,4]` of the saved volume | unchanged | `run_normalise_groups` saves four, so all four; the old comment naming three is corrected. Since step 8 the setting names the mice, `behavior_mice = {'MG705_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1'}`, matched to the names saved with the volume; a name not saved stops the run |
| `generate_diff_videos` | `true` | `true` | the approved folders hold every video |
| `generate_individual_diff_videos` | `true` | `true` | |
| `generate_t_scored_videos` | `true` | `true` | |
| `generate_surprise_videos` | `true` | `true` | |
| `generate_rolling_videos` | `true` | `true` | |
| `generate_signed_diff_videos` | `true` | `true` | |
| `perform_area_based_analysis_fine` | `false` | `false` | |
| `perform_area_based_analysis_coarse` | `false` | `false` | |
| `apply_smoothing` | `true` | `true` | |
| `smooth_sigma` | `5.0` | `5.0` | voxels, 3D Gaussian over each mouse's tissue; the old code set the voxels outside the tissue to 0 after it. Since fix 23 of step 8 the Gaussian is normalised by the smoothed tissue mask and the voxels outside the tissue stay NaN, left out of every mean, SEM and t |
| `min_mice_per_group` (new, 4 Oct 2026) | `3` | none: every mouse counted at every voxel (2, fixed in the code, from fix 23 of step 8) | the fewest mice with a value a voxel needs in each group to get a t and a surprise. With 2, voxels at tissue edges gave \|t\| up to 73 (naive against rws); Giulio's decision of 4 Oct 2026 |
| `bar_measure` (new, 6 Oct 2026) | `'cluster'` | none: the bars were the summed surprise, then (4 Oct) the share | the measure the region bars draw; all five (`share`, `sum`, `q99`, `topvol`, `cluster`) are in the table and the comparison figure whatever the choice |
| `n_permutations` (new, 6 Oct 2026) | `'all'` | none | splits of the pooled mice for the region permutation test: all of them (252 for 5 against 5, 126 for 5 against 4), or a number for a random subset with seed 0 |
| `cluster_p`, `cluster_connectivity` (new, 6 Oct 2026) | `0.01`, `18` | none | a voxel joins a cluster from p < 0.01 of the rolling median; voxels touch by a face or an edge |
| `topvol_mm3` (new, 6 Oct 2026) | `0.1` | none | the volume over which the top-volume measure averages a region's most surprising voxels (100,000 voxels of 10 um) |
| `region_quantile` (new, 6 Oct 2026) | `0.99` | none | the quantile measure |
| `permutation_workers` (new, 6 Oct 2026) | `16` | none | thread workers of the permutation test; changes the time, not the result |
| `channel` | `'nano'` | `'nano'` | |

Fixed in code, unchanged:

| what | value |
|---|---|
| atlas | `data\atlas\annotation_10.nii.gz`, AP crop `[180 1079]` (CCF only) |
| experimental onto control fit | planes `200:700`, first-order polyfit between the group means of each mouse's mean tissue intensity per plane |
| shared scale | planes `300:500` |
| group videos | LR diff and sum `[0 2.5]`, group difference `[-2.5 2.5]` |
| individual videos | `[0 1.5]` and `[0 10]`; signed `[-0.75 0.75]` |
| t videos | `[-6 6]`; surprise video `[0 8]`; surprise-masked video at p < 0.05 |
| slab figure | plane 565 +/- 10, median over the slab, alpha rising to full at p < 0.01, t `[-6 6]`; since 4 Oct 2026 the median is over the voxels with a t and a pixel is shown only where plane 565 has a t (before, wherever a plane of the slab had tissue of both groups) |
| individual slab figure | plane 565 +/- 10, `[0 1.5]` and `[0 10]`, no atlas overlay |
| rolling videos | +/- 10 planes, p < 0.01; individual `[0 2]` and `[0 10]`; since 4 Oct 2026 the t video shows each frame on its central plane's voxels with a t, as the slab figure |
| regional surprise bars | rolling median +/- 10 planes, surprise summed above p < 0.01 (since 4 Oct 2026 the bar is the fraction of the region's voxels with a t that is at p < 0.01, the sum kept in the table), left hemisphere, 56 listed regions ("Parafascicular nucleus" twice, and "Mediodorsal nucleus of the thalamus" resolving to the intermediodorsal nucleus); 55 since fix 9 of step 8, the second parafascicular removed and the mediodorsal nucleus measured; since 4 Oct 2026 over the voxels with a t only (the rolling median no longer fills a voxel without one from its neighbours), and 71 regions: the 43 isocortical areas at the atlas's structure level (S1 as its seven subfields; ACA, ORB, AI and RSP as their parts) and the 28 subcortical regions of before, by acronym, each region of the list taken out of the region holding it (HPF without SUB, HY without STN and ZI, MBmot without SCm), so no voxel counts in two bars |

- Reads: `naive\` and `<exp>\` `nano_4d_normalized.mat` and
  `nano_4d_normalized_bkgmask.mat` (this run's `run_normalise_groups`);
  `data\atlas`.
- Writes: `data\comparisons\naive_vs_<exp>_nano\`. No `.mat`: the values are
  in the `.fig` files, and the regional bars' in
  `Region_Surprise_DiffSum_<comp>.csv` (since 4 Oct 2026).
- Hidden state: `perform_area_based_analysis_fine` would read
  `WholeBrain_TMap_*.mat` from an earlier run (off). In the old script `t_lim`
  is set inside the t-video block and used by the surprise block, so a quick
  check must turn `generate_t_scored_videos` and `generate_surprise_videos`
  off together (or all six video flags).

Against the outputs Sami approved (S6). The approved folders are
`comparisons\naive_vs_rws\` and `naive_vs_behavior\`, without `_nano`, from
before the channel was added to the comparison tag. Every name that carries the
tag gains `_nano` after it (`Slab_Avg_565_naive_vs_rws_surpmask` becomes
`Slab_Avg_565_naive_vs_rws_nano_surpmask`, and likewise
`Normalization_Profiles_LR_*`, `Region_Surprise_Bar_DiffSum_*` and the
group-difference videos). Names without the tag are unchanged
(`Indiv_Slab_Avg_565_<group>`, `Individual_*_<group>.mp4`,
`lr_diff_sum_nano_<group>.mp4`). The `bk\` subfolders are not part of the
reference.

The inputs of the approved run appear to survive: the naive and rws
normalised volumes of 27 and 26 Nov 2025 as `nano_4d_normalized_bk.mat` (five
mice each, the older format without a `channel` field), and behavior's of 5
Dec 2025 as the current `behavior\nano_4d_normalized.mat` and its mask (MG705,
MG709, MG716, MG718). The naive and rws masks of that time were overwritten on
4 Sep and 1 May 2026. The fresh `run_normalise_groups` differs from them at
least by the raw-zero fix (4 Sep 2026). A `run_group_differences` run on those
files would separate a difference due to the code from one due to the
renormalisation. TO CONFIRM.

## Registration checks (MATLAB)

### run_register_to_atlas.m (P4_register_to_atlas.m)

The working copy differs from the commit in its run settings (MG914,
`annotate`, `demba_p28`); they stay uncommitted.

| setting | committed | working copy | reference | why |
|---|---|---|---|---|
| `groups_to_process` | `{'young'}` | `{'young'}` | MG914 `{'young'}`, CGF027 `{'naive'}` | ignored once `mice_to_process` is set |
| `mice_to_process` | `{'MG904_SepGluA_P22'}` | `{'MG914_SepGluA_P28'}` | `{'MG914_SepGluA_P28'}` or `{'CGF027_Gria1'}`, one per run | |
| `run_mode` | `'register'` | `'annotate'` | `'register'`; MG914 also `'autoannotate'` | stages 2 and 3 |
| `atlas_extent_slices` | `15` | `15` | `15` | no effect: both mice have control points, so each keeps the margin in its `regopts.mat` and nothing is written back |
| `atlas_key` | `'demba_p22'` | `'demba_p28'` | MG914 `'demba_p28'`, CGF027 `'ccf'` | the driver refuses an atlas that does not match the age or the mouse's `px_atlas` |
| `correction_type` | `'slicewise'` | same | unchanged | used by `align` only; the correction files hold their own, which replaces it |
| `use_equalized_nano` | `1` | same | unchanged | used by `align` only |
| `allow_image_only_registration` | `false` | same | `false` | both mice have control points on every slice |

- Reads, `register`: in `lightsuite\`, `regopts.mat`, the file matching
  `*tform.mat` (`atlas2histology_tform.mat`), `cutting_angle_data.mat`, the
  file matching `*20um.tif`, `sliceinfo.mat`, `volume_aligned\`;
  `local_settings.txt` (mouse folder first, then `lightsuite\`) for the
  `px_atlas` check; the atlas through `which('average_template_10.nii.gz')`.
- Writes, `register`: `elastix_forward\`, `elastix_reverse\`,
  `transform_params.mat`, `volume_registered\`.
- Reads, `autoannotate`: `regopts.mat`, `plane_anchors.mat`,
  `auto_atlas_planes.mat`, `volume_for_inspection.tiff`, the weights in
  `registration\auto_annotation\weights\` (v1.0, tracked). Writes
  `auto_proposal_controlpoints.mat` and `auto_proposal_info.mat` (with date,
  device and model version).
- Hidden state: the three globs pick up any stray file with a matching name
  (today each folder holds exactly one). In the old code the registered grid
  came from `sliceinfo.mat`'s `px_atlas`, not from the atlas: MG914's still
  says 10, written by `run_extract_and_center` on 12 Aug under the pre-DeMBA
  settings, and that is how every young brain was registered. Since fix 26 of
  step 8, 'register' sets the 10 um grid itself
  (`registration/pipeline/registered_grid_um.m`), stops before registering
  when a brain's `px_register` is not 20, and only prints a line for a
  `px_atlas` other than 10. Elastix's sampler is unseeded, so tolerances come
  from the old-against-old pair. The saved proposal ran on `cuda`; its anchors
  and planes are the ones in `plane_anchors.mat` today.

### run_add_sep_channel.m (P4bis_add_sep_channel.m)

| setting | committed | reference | why |
|---|---|---|---|
| `groups_to_process` | `{'young','naive','rws'}` | unchanged | overridden by `P4BIS_MICE` |
| `mice_to_process` | `{}` | `P4BIS_MICE=MG914_SepGluA_P28`, then `CGF027_Gria1` | stage 2, after `register` on the same mouse |
| `overwrite` | `false` | `false` | the tree gets no `volume_aligned_sep\` or `volume_registered_sep\` as input, otherwise the mouse is skipped (CGF027 has both in production). TO CONFIRM |
| `recovery_levels` | `[4 2 1]` | unchanged | |
| `min_slice_corr` | `0.99` | unchanged | |

- Reads: `volume_centered\` DAPI and EGFP (first file whose name contains
  `dapi` or `egfp`), `volume_aligned\chan01_DAPI.tiff`, `sliceinfo.mat`,
  `volume_for_ordering_processing_decisions.txt`, `transform_params.mat` and
  `elastix_reverse\`, `volume_registered\chan01_DAPI.tiff` (for the check).
- Writes: `volume_aligned_sep\`, `volume_registered_sep\`, `sep_work\`
  (removed at the end), `comparisons_v2\processing_diagnostics\sep_channel\<mouse>.png`
  and `.txt` (the text carries the date).
- Hidden state: the elastix folder named in `transform_params.mat` is used as
  it is when it exists. In a `transform_params.mat` copied from production it
  names the production folder on D:, so `run_add_sep_channel` must run after
  `register` in the same tree, which rewrites it. The check compares against
  the `volume_registered\` of that same `register` run. MG914 has no
  production SEP output: it is compared old against old only.

### run_extract_and_center.m (P1_extract_and_center_data.m)

| setting | committed | reference | why |
|---|---|---|---|
| `groups_to_process` | `{'young'}` | unchanged | |
| `mice_to_process` | MG909, MG910, MG911, MG912, MG913, MG914 | `{'MG914_SepGluA_P28'}` | stage 1; extraction is per mouse |
| `atlas_key` | `'ccf'` | unchanged | only put on the path |

- Reads: `young\MG914_SepGluA_P28\*.czi` (seven files), `local_settings.txt`
  (mouse folder first), `lightsuite\volume_for_ordering_processing_decisions.txt`.
- Writes: `lightsuite\volume_centered\chan01_DAPI`, `chan02_Cy5`,
  `chan03_Cy3`, `chan04_EGFP` `.tiff`, `volume_for_ordering.tiff`,
  `sliceinfo.mat`, `volume_ordered.tiff`.
- Hidden state: `local_settings.txt` now holds the DeMBA values (`px_atlas`
  20, `atlasaplims [64 560]`, changed 29 Sep); production's `sliceinfo.mat`
  was written on 12 Aug with 10 and `[180 1079]`. The extraction settings
  (`px_process` 5, `px_register` 20, `slicethickness` 150, `regchan` dapi) are
  the same, so only those two fields and the absolute paths should differ. The
  point subsampling in `extractBrainLimits3` is unseeded.

### run_residual_correction.m (P2_residual_correction_analysis.m)

| setting | committed | reference | why |
|---|---|---|---|
| `groups_to_process` | `{'young'}` | unchanged | |
| `mice_to_process` | the same six | `{'MG914_SepGluA_P28'}` | per mouse, independent of the others |
| `atlas_key` | `'ccf'` | unchanged | only put on the path; the old script also loaded the annotation without using it (step 8 removed the load) |
| `doPlotBkg` | `true` | unchanged | |
| `savePlotBkg` | `true` | unchanged | |
| `saveRatioMap` | `false` | unchanged | |

Fixed in code: reference pixels with `range_frac` 0.20 and a window of 75
(first nine slices) or 50; bisquare `robustfit`; both corrections,
`slicewise` and `global`, computed and saved (L4 archives `global` later);
`randn` only jitters dots in a figure.

- Reads: `volume_centered\chan02_Cy5.tiff` (nano), `chan03_Cy3.tiff` (auto).
- Writes: `correction_output\corrected_volume_{slicewise,global}.mat`,
  `scaled_auto_volume_{slicewise,global}.mat`,
  `scaled_difference_video_{slicewise,global}.mp4`, `diagnostic_plots\`.

### run_nano_equalisation.m (P2bis_nano_correction_analysis.m)

| setting | committed | reference | why |
|---|---|---|---|
| `groups_to_process` | `{'young'}` | unchanged | |
| `mice_to_process` | the same six | the same six. TO CONFIRM | see hidden state |
| `atlas_key` | `'ccf'` | unchanged | only put on the path |
| `save_results` | `true` | unchanged | |
| `base_output_dir` | `data\intensity_diagnostics` | unchanged | |

Fixed in code: background percentiles 15/75 for the first nine slices, then
15/50; moving median of the slice medians over 5 slices; `parfor`.

- Reads: `volume_centered\chan02_Cy5.tiff` of every listed mouse.
- Writes: `correction_output\equalized_volume.mat` per mouse;
  `intensity_diagnostics\Intensity_Stats_<time>.mat`, `Plot_*_<time>.png`,
  `Video_<mouse>_<time>.mp4` (names carry the run's clock time).
- Hidden state: the mice are stacked into one array padded to the largest
  height and width, and the padding is filled with the slice's mode before
  the background is selected. A brain's equalisation therefore depends on
  which brains run with it. Production ran on these six on 12 Aug; MG914 is
  the tallest but MG911 is wider (2472 against 2259 px), so MG914 alone may not
  reproduce its production `equalized_volume.mat`. Running the six needs the
  other five's `chan02_Cy5.tiff` in the tree.

## Python route

No setting changes for the reference run: the committed values are the
reference values. They were constants in the `v2_*.py` modules; they are now
keys of `mapping/settings.toml` (`table.key` below, same values). Values
marked "in code" are fixed by the method and stay in their module. Only the
environment (above) and the options below are set; every run script runs with
no option except `run_closeup` and the `--panel` of the two ISH passes.

### Volumes and young against adult

| run script | options | settings in force | reads | writes |
|---|---|---|---|---|
| `run_per_mouse` | none (all 17 in `MICE`) | `tissue.mad_k` 4; 10 adults on `ccf`, 5 P20 on `demba_p20`, MG911 `demba_p16`, MG904 `demba_p22` | `volume_registered\` nano and auto, `volume_registered_sep\chan02_SEP.tiff`, the atlases | `comparisons_v2\per_mouse\<mouse>.npz` |
| `run_to_ccf` | none | in code: DeMBA canvas 705 x 400 x 570, CCF 660 x 400 x 570, adult crop planes 90 to 540 | per-mouse files | `per_mouse_ccf\<mouse>.npz` |
| `run_cohort` | none | `readings.ratio_clip` 20, `readings.log2_floor` 0.02, `readings.in_force` (five readings), `readings.signed` zref; cohorts young, young_P20, young_P16, young_P22, naive, rws, adult | both per-mouse folders | `ccf\<cohort>\*_{mean,sd,n}.npy`, `mice.txt`; the cache `per_mouse\<mouse>_scalars.npz`; since fix 24 of step 8 also `ccf\<cohort>\*_folded_{mean,sd,n}.npy` (each brain's hemispheres averaged first) |
| `run_compare` | none | `young_vs_adult.smooth` 1 voxel (log2 map only), `.min_n_young` 2, `.min_n_adult` 5, `.min_table_vox20` 100; young and young_P20 | cohort volumes | `young_vs_adult\volumes_ccf20.npz`, `region_table.csv`, `cortex_table.txt`, `slices_<reading>.png/.eps` |
| `run_region_plot` | none | `region_tables.min_vox20` 250, `region_plot.min_young` 2 and `.min_adult` 4, zref over the structures shared by all 17 brains, BH within each reading | per-mouse files | `region_means_per_mouse.csv`, `region_stats.csv`, `region_plot.png/.eps` |
| `run_region_groups` | none | `region_groups.min_young` 3 and `.min_adult` 5; systems, divisions and layers in code | per-mouse files | `group_stats.csv`, `group_plot`, `laminar_plot` |
| `run_video` | none: young, adult, young_P20, naive, rws | `videos.mean_vmax` per reading, `videos.t_pct` 95, `videos.min_n` per cohort, `videos.fps` 12 | cohort volumes; since fix 24 the t from the `_folded` ones | `ccf\<cohort>\video_<reading>_<cohort>.mp4` (25 videos) |
| `run_video_compare` | none: all five readings, no still | `videos.mean_vmax`, `videos.log2_lim` 1.5, `videos.fps` 12 | cohort volumes | `young_vs_adult\video_side_by_side_<reading>.mp4` |
| `run_closeup` (venv_flat) | none, then `--cmap jet`, then `--cmap turbo`. TO CONFIRM | zref, `closeup.plane` 790, `closeup.vmax` 0.9, `closeup.dlim` 0.5, `closeup.smooth` 3 x 1 x 1 voxels, video and flatmaps | cohort volumes, `data\atlas_flatmap\` | `young_vs_adult\detail_*_zref.*`, and `jet\`, `turbo\` |
| `run_diagnostics` | none (the cohort sheets and the index need no options) | `tissue.mad_k` (sheets 01 and 02); the rest in code | per-mouse files, cohort volumes, `region_table.csv`, `region_stats.csv`, `naive\nano_4d_normalized_bkgmask.mat` (CGF027), `young\nano_4d_normalized_bkgmask_P20.mat` (MG903) | `processing_diagnostics\` sheets 01 to 09 and `README.md` |

### ISH, both passes

| run script | pass | settings in force | reads | writes (under `adult_v2\`) |
|---|---|---|---|---|
| `run_ish_regions --panel targets` | 1 | `ish_regions.min_voxels` 3; in code: 200 um grid, missing = -1, reference box 67 x 41 x 58, one-voxel erosion | `data\gene_targets.csv`, `data\atlas_ish\` | `ish\gene_region_table.csv`, `ish\gene_region_table_drops.csv` |
| `run_ish_compare` | 1 | `ish.min_voxels` 10, `ish.min_structures` 50, `ish.min_genes_ranking` 20, five readings | pass-1 table, `region_means_per_mouse.csv`, the frozen `comparisons\merged_naive_rws_vs_ish_summary_nosmooth\gene_panel_summary.csv` of `run_compare_with_allen_ish` (P9) | `ish\gene_correlations.csv`, `ish_old_vs_new.png` |
| `run_ish_words` | 1 | `ish_words.min_genes` 5, `.max_share` 0.8, `.n_boot` 2000, seed 0 | `gene_correlations.csv`, the mygene cache `ish\annotation\` (95 files) | `feature_enrichment.csv`, `ish_word_enrichment.png` |
| `run_ish_roles` | 1 | `ish.reading` zref, the curated roles (in code), seed 0; since step 8 `ish_roles.evidence_p` 0.05 (the figure title's wording only) | pass-1 table, `region_means_per_mouse.csv` | `gene_roles.csv`, `role_summary.csv`, `ish_roles.png` |
| `run_panel_build` | 2 | `ish_panel_build.max_hits` and `.allen_rows`; the GO terms and the Grid1/Grid2 override in code | `gene_targets.csv`, the API cache `panel\cache\` (436 files) | `panel\panel_v2.csv`, `panel_genes.csv` |
| `run_panel_fetch` | 2 | in code: reference box 67 x 41 x 58, timeout 180 s, 2 retries | `panel_v2.csv`, `data\atlas_ish\` | missing grids; `panel\fetch_failures.csv` |
| `run_ish_regions --panel ontology` | 2 | as pass 1 | `adult_v2\panel\panel_v2.csv`, `data\atlas_ish\` | `ish\gene_region_table_panel.csv`, `gene_region_table_panel_drops.csv` |
| `run_ish_reliability --panel ontology` | 2 | `ish.min_voxels` 10, `ish.min_structures_pair` 50, seed 0 | `gene_region_table_panel.csv` | `gene_reliability.csv`, `gene_region_table_merged.csv`, `ish_reliability.png` |
| `run_ish_panel_test` | 2 | `ish.reading` zref, `ish_panel_test.min_structures` 80, `.min_reliability` 0.3, `.n_perm` 20000, seed 0 | merged table, reliability, `region_means_per_mouse.csv` | `panel_test.csv`, `ish_panel_test.png` |

### Adult arms and beyond abundance

| run script | settings in force | reads | writes (under `adult_v2\`) |
|---|---|---|---|
| `run_adult_arms` | `region_tables.min_vox20` 250; arms sepauto, ratio, sepratio | per-mouse files, `region_means_per_mouse.csv` | `arms\region_means_arms.csv`, `arms_consistency.png` |
| `run_ish_arms` | `ish.control_gene` Gria1, `ish.min_structures` 50, seed 0 | `region_means_arms.csv`, pass-1 table | `arms\arm_gene_correlations.csv`, `arms_vs_genes.png` |
| `run_sep_channel_check` | `region_tables.min_vox20` 250, seed 0 | per-mouse files, pass-1 table | `arms\sep_channel_check.csv/.png` |
| `run_beyond_density` | `beyond.reading` zref, `region_tables.min_vox20` 250, `beyond.half` 5, `beyond.grey` (the grey-matter rule), seed 0; since step 8 `beyond_controls.replication` 0.5 also decides whether its titles and printed lines say the map and the leftover replicate (so does `run_beyond_figures`) | `region_means_per_mouse.csv`, `gene_region_table_merged.csv` | `beyond\structures_used.csv`, `variance_partition.csv`, `residual_by_structure.csv`, `fig0` to `fig3` |
| `run_beyond_controls` | 5 folds, `beyond_controls.max_pcs` 25, the controls' bounds in `[beyond_controls]`, seed 0 | as above | `beyond\controls.csv`, `fig4` to `fig6` |
| `run_beyond_figures` | `beyond_figures.n_boot` 2000, `.n_perm` 10000, `.boot_splits` 20, seed 0 | as above, `controls.csv` | `beyond\for_sami\A` to `D`, `numbers_for_the_caption.txt` |
| `run_beyond_regression` | `beyond_regression.planes` 215, 265, 315; `.floor` 0.12; since step 8 `.diagnostic_p` 0.05 (panel E's wording only) | as above, per-mouse files | `beyond\for_sami\E_regression`, `F_maps`, `regression_table.csv` |

### Hidden state in the Python route

- `run_cohort` reuses `per_mouse\<mouse>_scalars.npz` when its recorded
  modification time matches the per-mouse file. A copy that keeps modification
  times would reuse it: the caches are deleted before each run (plan), and
  before the first production run after step 8, whose `subref` fix changes
  the subcortex mean they hold.
- Since fix 24 of step 8, `run_video` stops, naming `run_cohort.py`, when a
  cohort has no `_folded` files.
- `run_per_mouse` includes SEP only where
  `volume_registered_sep\chan02_SEP.tiff` exists; with it missing,
  `run_cohort` stops (sepratio needs it). All 17 have it (23 and 24 Sep 2026).
- `run_panel_build` and `run_ish_words` call the network only for what their
  caches lack. `run_panel_fetch` skips grids on disk, but retries on every run
  the two that failed in production (Gria1 526 and Negr1 695, "no energy.mhd in
  the zip"), so it always calls the network. TO CONFIRM.
- The old `v2_ish_compare` skipped the old-against-new check without a word
  when the summary of `run_compare_with_allen_ish` (P9) was missing;
  `run_ish_compare` stops instead (step 8).
  `run_adult_arms` skips its consistency check when
  `region_means_per_mouse.csv` is missing.
- `run_diagnostics` sheet 08 reads `young\nano_4d_normalized_bkgmask_P20.mat`
  (`run_normalise_groups` for young_P20, 4 Sep 2026, three mice), which the
  reference run does not regenerate. TO CONFIRM.
- A figure open in a viewer is written as `*_new.png` instead.

## Expected differences from the production outputs

Not refactor errors; the output comparison maps or excludes them.

- `run_group_differences`: the `_nano` in folder and file names (above);
  inputs renormalised since 5 Dec 2025.
- `run_collect_by_group`: new `collected_mice.mat` for the three adult groups.
- `run_extract_and_center` (MG914): `px_atlas`, `atlasaplims` and the paths
  inside `sliceinfo.mat`.
- `run_nano_equalisation`: the clock time in `intensity_diagnostics\` names.
- `run_register_to_atlas` mode `autoannotate`: date and device in
  `auto_proposal_info.mat`; `run_add_sep_channel`: the date line in the
  `.txt`.
- ISH pass 1 writes `gene_region_table_drops.csv`; production holds the older
  name `gene_drops.csv`.
- `young_vs_adult\barrel_only\` (27 Sep) is written by no code of the
  repository: by `D:\sep_histology\sandbox_barrel\group_plot_barrel_only.py`
  (item 10 below).
- `processing_diagnostics\README.md` is regenerated.

## TO CONFIRM

As drafted on 30 Sep, before the reference run, each with its outcome.

1. Settings with no environment override (the groups and age filter of
   `run_collect_by_group`, `exp_type` of `run_group_differences`, the mouse
   lists and modes of `run_extract_and_center`, `run_residual_correction`,
   `run_nano_equalisation` and `run_register_to_atlas`): edited in the
   worktree copy and saved as a patch beside the outputs, with the same
   values carried into the branch's drivers. Outcome: the settings were
   changed in driver copies (`tools/sep_run_driver_copy`), kept in each check
   tree's `driver_copies\`.
2. The two analysis environments: run from `main` (`tools\venv_atlas`,
   `tools\venv_flat`) on the worktree's scripts. Outcome: as proposed.
3. S6: also run today's `run_group_differences` on the normalised volumes of
   Nov and Dec 2025 (`naive` and `rws` `nano_4d_normalized_bk.mat`,
   `behavior`'s current file), to separate a code difference from the
   renormalisation. Outcome: run on 1 Oct; it reproduced the approved
   figures (REFACTOR_PLAN.md, Progress, "1 Oct, S6, informative").
4. `run_diagnostics` sheet 08: copy `young\nano_4d_normalized_bkgmask_P20.mat`
   into the tree as a frozen input (proposed), or add `run_collect_by_group`
   with `age_filter [20]` and `run_normalise_groups` for `young_P20` to the
   reference run (it would now stack five P20 brains instead of three).
   Outcome: copied as a frozen input, as proposed.
5. `run_panel_fetch`: run it with the network (it retries the two failures),
   or leave it out and copy `fetch_failures.csv` as a frozen input (proposed:
   leave it out; the grids it would fetch are already on disk). Outcome: left
   out, as proposed.
6. `run_closeup`: include the `--cmap jet` and `--cmap turbo` sets, as in
   production, or only the default set. Outcome: only the default set.
7. `run_nano_equalisation` stage 1: the six mice of the production run
   (proposed), or MG914 alone and accept a difference. Outcome: MG914 alone,
   old code against new on the same input, in `G:\sep_refactor\pre\`.
8. `run_add_sep_channel`: leave `volume_aligned_sep\` and
   `volume_registered_sep\` out of the registration tree's inputs (proposed),
   or set `overwrite = true`. Outcome: left out, as proposed.
9. `run_register_to_atlas` mode `autoannotate`: on the GPU, as the saved
   proposal was made, when no other process holds it. Outcome: on the GPU.
10. `young_vs_adult\barrel_only\`: which code wrote it, and whether it is
    excluded from the comparison. Outcome: written by
    `D:\sep_histology\sandbox_barrel\group_plot_barrel_only.py`, outside the
    repository, so no run of either code rewrites it; ROADMAP.md, section 8,
    proposes an option of `run_region_groups` for it.
