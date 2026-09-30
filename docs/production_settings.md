# Production settings for the reference run

Draft for Giulio's sign-off, 30 Sep 2026 (step 0.6 of
[REFACTOR_PLAN.md](REFACTOR_PLAN.md)). It lists, driver by driver, the settings
the reference run of step 3 uses: the old code, at the tag `refactor-start`,
run in the check trees on G:. Every later run of the new code uses the same
values, so old and new are compared under identical settings.

How to read the tables:

- committed: the value in commit `04c0484`, the last one before step 0. Step 0
  changes no setting value. It adds `SEP_DATA_ROOT` and the data-root guards,
  the guard in P4 `align`, and lets `V2_ISH_PANEL` be a path relative to the
  data root.
- reference: the value the reference run uses, and why.
- hidden state: what a run depends on besides its settings (caches reused
  when present, files read from earlier runs, globs, index-based selections).
- Points that need a decision are marked TO CONFIRM and collected at the end.

## Common to every run

| item | reference |
|---|---|
| code | cohort tree: worktree at `refactor-start` in `G:\sep_refactor\ref\code`; registration checks: `G:\sep_refactor\reg\old\code` |
| MATLAB | R2024b, a fresh session per driver, `restoredefaultpath` then `sep_setup_paths`, started by `tools\run_matlab_detached.ps1` with `-CodeDir`, `-DataRoot` and `-LogDir` |
| toolboxes | Image Processing, Statistics, Parallel Computing (P2bis uses `parfor`) |
| elastix | 5.1.0, installed per user and on the PATH (P4 `register`, P4bis) |
| settings with no environment override | set by editing the settings block in the worktree copy, the edit saved as a patch beside the outputs. TO CONFIRM |

### Environment variables

| variable | read by | committed default (unset) | reference |
|---|---|---|---|
| `SEP_DATA_ROOT` | `get_paths.m`, `v2_paths.py`, `build_demba_atlas.py` (from step 0) | `<code>\..\data` | the check tree's data folder, in every run |
| `SEP_COHORT_SPECS` | P6bis | `{'rws','naive'}` | `rws,naive,behavior` (S6) |
| `P4BIS_MICE` | P4bis | groups young, naive, rws | `MG914_SepGluA_P28`, then `CGF027_Gria1` (registration tree only) |
| `SEP_MERGE_SPECS` | P8 | P8's own list | unset: P8 is not run |
| `V2_READINGS` | `v2_cohort`, and through it every v2 script that draws or tabulates | all five readings | unset |
| `V2_ISH_PANEL` | `v2_ish_regions` | `data\gene_targets.csv` (100 genes) | pass 1 unset; pass 2 `adult_v2\panel\panel_v2.csv` (relative, inside the tree's data) |
| `V2_ISH_TABLE` | `v2_ish_regions`, `v2_ish_reliability` | `gene_region_table.csv` in `v2_ish_regions`, `gene_region_table_panel.csv` in `v2_ish_reliability` | pass 1 unset; pass 2 `gene_region_table_panel.csv` |
| `AUTO_ANNOTATION_PYTHON` | `auto_annotate.m` | `auto_annotation\.venv` beside the code | `D:\sep_histology\code\auto_annotation\.venv\Scripts\python.exe` (a worktree has no `.venv`) |
| `LANDMARK_REFINE_PYTHON` | `landmark_refine*.m` (P4 `annotate`) | the landmark_refine venv | unset: `annotate` is not run |

Both ISH variables must be unset for every script outside pass 2. With
`V2_ISH_TABLE` left set, pass 1's `v2_ish_regions` would overwrite the
390-gene table.

### Interpreters

| environment | path | Python | used by |
|---|---|---|---|
| analysis | `D:\sep_histology\code\tools\venv_atlas` | 3.12.7 | every `v2_*.py` except `v2_inspect`; holds the DeMBA-to-CCF deformation fields `v2_to_ccf` needs (downloaded on first use if missing) |
| flatmaps | `D:\sep_histology\code\tools\venv_flat` | 3.12.7 | `v2_inspect` only (ccf_streamlines) |
| automatic annotation | `D:\sep_histology\code\auto_annotation\.venv` | 3.12.7 | P4 `autoannotate` (torch, GPU) |

The environments are ignored by git, so the worktrees have none. The
reference run calls `main`'s interpreters on the worktree's scripts
(`<python> <tree>\code\v2_x.py`), with `SEP_DATA_ROOT` set. TO CONFIRM for the
two analysis environments (the plan states it for the engine only).

## Run order

### Cohort tree

1. P5 for rws, naive, behavior.
2. P6bis for rws, naive, behavior.
3. P7bis naive against rws, then P7bis naive against behavior.
4. The Python route:

| # | module | stage |
|---|---|---|
| 1 | `v2_per_mouse` | volumes |
| 2 | `v2_to_ccf` | volumes |
| 3 | `v2_cohort` | volumes |
| 4 | `v2_compare` | young against adult |
| 5 | `v2_region_plot` | young against adult |
| 6 | `v2_region_groups` | young against adult |
| 7 | `v2_video` | young against adult |
| 8 | `v2_video_compare` | young against adult |
| 9 | `v2_inspect` (venv_flat) | young against adult |
| 10 | `v2_diagnostics` | diagnostics |
| 11 | `v2_ish_regions` | ISH pass 1 (100 genes) |
| 12 | `v2_ish_compare` | ISH pass 1 |
| 13 | `v2_ish_words` | ISH pass 1 |
| 14 | `v2_ish_roles` | ISH pass 1 |
| 15 | `v2_adult_arms` | adult arms |
| 16 | `v2_ish_arms` | adult arms |
| 17 | `v2_sep_channel_check` | adult arms |
| 18 | `v2_panel_build` | ISH pass 2 (390 genes) |
| 19 | `v2_panel_fetch` | ISH pass 2 |
| 20 | `v2_ish_regions` with both ISH variables set | ISH pass 2 |
| 21 | `v2_ish_reliability` | ISH pass 2 |
| 22 | `v2_ish_panel_test` | ISH pass 2 |
| 23 | `v2_beyond_density` | beyond abundance |
| 24 | `v2_beyond_controls` | beyond abundance |
| 25 | `v2_beyond_figures` | beyond abundance |
| 26 | `v2_beyond_regression` | beyond abundance |

Not run: `v2_replot` (it only redraws `v2_compare`'s figures from the volumes
`v2_compare` saved).

What fixes the order:

- `v2_diagnostics` sheet 08 reads P6bis's background masks, and sheet 07 reads
  `region_table.csv` (`v2_compare`) and `region_stats.csv` (`v2_region_plot`).
- `young_vs_adult\region_means_per_mouse.csv` (`v2_region_plot`) is read by
  `v2_ish_compare`, `v2_ish_roles`, `v2_ish_panel_test`, `v2_adult_arms` (a
  consistency check) and the four `v2_beyond_*`.
- `adult_v2\ish\gene_region_table.csv` (pass 1) is read by `v2_ish_compare`,
  `v2_ish_roles`, `v2_ish_arms` and `v2_sep_channel_check`.
- `v2_ish_words` reads `gene_correlations.csv` (`v2_ish_compare`); `v2_ish_arms`
  reads `arms\region_means_arms.csv` (`v2_adult_arms`).
- `v2_beyond_density` and `v2_beyond_controls` read
  `gene_region_table_merged.csv` (pass 2, `v2_ish_reliability`);
  `v2_beyond_figures` draws panel D only if `controls.csv`
  (`v2_beyond_controls`) exists.

Within these constraints the order is free; the production run of 24 to 26
Sep 2026 differs only in running `v2_ish_roles` after the arms.

### Registration tree

Inputs are restored from `reg\old` before each stage, so no stage reads
another's outputs.

1. MG914_SepGluA_P28: P1, P2, P2bis.
2. MG914_SepGluA_P28 and CGF027_Gria1: P4 `register`, then P4bis on the same
   mouse.
3. MG914_SepGluA_P28: P4 `autoannotate`.

The restore matters most between stages 1 and 2: P1 rewrites `sliceinfo.mat`,
and P4 `register` and P4bis read the `px_atlas` inside it to size the
registered volume (see P1 and P4).

## Plasticity comparison (MATLAB)

### P5_collect_data_by_group.m

| setting | committed | reference | why |
|---|---|---|---|
| `mousetypes_list` | `{'young'}` | `{'rws','naive','behavior'}` | the three adult groups P6bis and P7bis read (S6) |
| `age_filter` | `[20]` | `[]` | the adults have no age; `[]` writes plain `nano_4d.mat` |
| `skip_missing` | `true` | `true` | all 17 adults are registered (26 Nov 2025): P5 must print 5, 5 and 7 mice |
| `correction_type` | `'slicewise'` | unchanged | not used by P5 |
| `min_num_contrib` (inline) | `3` | unchanged | feeds `avg_mask`, which is no longer saved |

- Reads: `<group>\<mouse>\lightsuite\volume_registered\chan02_NANO.tiff`,
  `chan03_AUTO.tiff`, `chan05_MASK.tiff` of every mouse in the group.
- Writes: `<group>\nano_4d.mat`; `<group>\collected_mice.mat` (new for the
  adults, whose stacks predate that file).
- Hidden state: only `nano_4d` is saved; the other saves are commented out,
  so the `auto_4d.mat`, `mask_4d.mat` and averages of 19 Nov 2025 would stay
  on disk. They are not copied into the tree, and nothing in the reference run
  reads them. A mouse whose tiff is missing is skipped with only a printed
  note, and P6bis's index-based selection would then pick the wrong mice.

### P6bis_analyze_group_averages_and_normalize.m

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

- Reads: `<group>\nano_4d.mat` (this run's P5); `data\atlas` through
  `get_atlas_crop('ccf')`.
- Writes: `<group>\nano_4d_normalized.mat`,
  `<group>\nano_4d_normalized_bkgmask.mat`,
  `<group>\Background_mask_diagnostics_trace_<group>_nano.fig/.png`,
  `<group>\global_diagnostics\normalization_checks_nano\`.
- Hidden state: an adult group is selected by position in the fourth
  dimension of `nano_4d.mat`, which assumes P5 stacked the whole group in
  registry order; `collected_mice.mat` is not consulted for the adults.

### P7bis_analyze_group_differences.m

Run twice, with only `exp_type` changed.

| setting | committed | reference | why |
|---|---|---|---|
| `ctrl_type` | `'naive'` | `'naive'` | |
| `exp_type` | `'behavior'` | run 1 `'rws'`, run 2 `'behavior'` | both comparisons Sami approved (S6) |
| `selected_mice_idx_list` | as P6bis | unchanged | used only for the mouse names in figures; the data are whatever P6bis saved |
| behavior subselection (inline) | `[1,2,3,4]` of the saved volume | unchanged | P6bis saves four, so all four; the comment naming three is corrected in Fixes |
| `generate_diff_videos` | `true` | `true` | the approved folders hold every video |
| `generate_individual_diff_videos` | `true` | `true` | |
| `generate_t_scored_videos` | `true` | `true` | |
| `generate_surprise_videos` | `true` | `true` | |
| `generate_rolling_videos` | `true` | `true` | |
| `generate_signed_diff_videos` | `true` | `true` | |
| `perform_area_based_analysis_fine` | `false` | `false` | |
| `perform_area_based_analysis_coarse` | `false` | `false` | |
| `apply_smoothing` | `true` | `true` | |
| `smooth_sigma` | `5.0` | `5.0` | voxels, NaN-tolerant 3D Gaussian over tissue |
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
| slab figure | plane 565 +/- 10, median over the slab, alpha rising to full at p < 0.01, t `[-6 6]` |
| individual slab figure | plane 565 +/- 10, `[0 1.5]` and `[0 10]`, no atlas overlay |
| rolling videos | +/- 10 planes, p < 0.01; individual `[0 2]` and `[0 10]` |
| regional surprise bars | rolling median +/- 10 planes, surprise summed above p < 0.01, left hemisphere, 56 listed regions ("Parafascicular nucleus" twice) |

- Reads: `naive\` and `<exp>\` `nano_4d_normalized.mat` and
  `nano_4d_normalized_bkgmask.mat` (this run's P6bis); `data\atlas`.
- Writes: `data\comparisons\naive_vs_<exp>_nano\`. No `.mat`: the values are
  in the `.fig` files.
- Hidden state: `perform_area_based_analysis_fine` would read
  `WholeBrain_TMap_*.mat` from an earlier run (off). `t_lim` is set inside the
  t-video block and used by the surprise block, so a quick check must turn
  `generate_t_scored_videos` and `generate_surprise_videos` off together (or
  all six video flags).

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
4 Sep and 1 May 2026. The fresh P6bis differs from them at least by the
raw-zero fix (4 Sep 2026). A P7bis run on those files would separate a
difference due to the code from one due to the renormalisation. TO CONFIRM.

## Registration checks (MATLAB)

### P4_register_to_atlas.m

The working copy differs from the commit in its run settings (MG914,
`annotate`, `demba_p28`); they stay uncommitted.

| setting | committed | working copy | reference | why |
|---|---|---|---|---|
| `groups_to_process` | `{'young'}` | `{'young'}` | MG914 `{'young'}`, CGF027 `{'naive'}` | ignored once `mice_to_process` is set |
| `mice_to_process` | `{'MG904_SepGluA_P22'}` | `{'MG914_SepGluA_P28'}` | `{'MG914_SepGluA_P28'}` or `{'CGF027_Gria1'}`, one per run | |
| `run_mode` | `'register'` | `'annotate'` | `'register'`; MG914 also `'autoannotate'` | stages 2 and 3 |
| `atlas_extent_slices` | `15` | `15` | `15` | no effect: both mice have control points, so each keeps the margin in its `regopts.mat` and nothing is written back |
| `atlas_key` | `'demba_p22'` | `'demba_p28'` | MG914 `'demba_p28'`, CGF027 `'ccf'` | P4 refuses an atlas that does not match the age or the mouse's `px_atlas` |
| `correction_type` | `'slicewise'` | same | unchanged | used by `align` only |
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
  `auto_annotation\weights\` (v1.0, tracked). Writes
  `auto_proposal_controlpoints.mat` and `auto_proposal_info.mat` (with date,
  device and model version).
- Hidden state: the three globs pick up any stray file with a matching name
  (today each folder holds exactly one). The registered grid comes from
  `sliceinfo.mat`'s `px_atlas`, not from the atlas: MG914's still says 10,
  written by P1 on 12 Aug under the pre-DeMBA settings, and that is how every
  young brain was registered. Elastix's sampler is unseeded, so tolerances come
  from the old-against-old pair. The saved proposal ran on `cuda`; its anchors
  and planes are the ones in `plane_anchors.mat` today.

### P4bis_add_sep_channel.m

| setting | committed | reference | why |
|---|---|---|---|
| `groups_to_process` | `{'young','naive','rws'}` | unchanged | overridden by `P4BIS_MICE` |
| `mice_to_process` | `{}` | `P4BIS_MICE=MG914_SepGluA_P28`, then `CGF027_Gria1` | stage 2, after P4 `register` on the same mouse |
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
  names the production folder on D:, so P4bis must run after P4 `register` in
  the same tree, which rewrites it. The check compares against the
  `volume_registered\` of that same `register` run. MG914 has no production
  SEP output: it is compared old against old only.

### P1_extract_and_center_data.m

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

### P2_residual_correction_analysis.m

| setting | committed | reference | why |
|---|---|---|---|
| `groups_to_process` | `{'young'}` | unchanged | |
| `mice_to_process` | the same six | `{'MG914_SepGluA_P28'}` | per mouse, independent of the others |
| `atlas_key` | `'ccf'` | unchanged | the annotation is loaded and not used |
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

### P2bis_nano_correction_analysis.m

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

No constant in any `v2_*.py` changes for the reference run: the committed
values are the reference values. Only the environment (above) and the
arguments below are set. Every module runs with no arguments except
`v2_inspect`.

### Volumes and young against adult

| module | arguments | settings in force | reads | writes |
|---|---|---|---|---|
| `v2_per_mouse` | none (all 17 in `MICE`) | `MAD_K` 4; 10 adults on `ccf`, 5 P20 on `demba_p20`, MG911 `demba_p16`, MG904 `demba_p22` | `volume_registered\` nano and auto, `volume_registered_sep\chan02_SEP.tiff`, the atlases | `comparisons_v2\per_mouse\<mouse>.npz` |
| `v2_to_ccf` | none | DeMBA canvas 705 x 400 x 570, CCF 660 x 400 x 570, adult crop planes 90 to 540 | per-mouse files | `per_mouse_ccf\<mouse>.npz` |
| `v2_cohort` | none | `RATIO_CLIP` 20, `Z_FLOOR` 0.02, five readings, `SIGNED_READINGS` zref; cohorts young, young_P20, young_P16, young_P22, naive, rws, adult | both per-mouse folders | `ccf\<cohort>\*_{mean,sd,n}.npy`, `mice.txt`; the cache `per_mouse\<mouse>_scalars.npz` |
| `v2_compare` | none | `SMOOTH` 1 voxel (log2 map only), `MIN_N_YOUNG` 2, `MIN_N_ADULT` 5, young and young_P20 | cohort volumes | `young_vs_adult\volumes_ccf20.npz`, `region_table.csv`, `cortex_table.txt`, `slices_<reading>.png/.eps` |
| `v2_region_plot` | none | `MIN_VOX` 250, zref over the structures shared by all 17 brains, BH within each reading | per-mouse files | `region_means_per_mouse.csv`, `region_stats.csv`, `region_plot.png/.eps` |
| `v2_region_groups` | none | systems, divisions and layers as in the file | per-mouse files | `group_stats.csv`, `group_plot`, `laminar_plot` |
| `v2_video` | none: young, adult, young_P20, naive, rws | `MEAN_VMAX` per reading, `T_PCT` 95, `MIN_N` per cohort, 12 fps | cohort volumes | `ccf\<cohort>\video_<reading>_<cohort>.mp4` (25 videos) |
| `v2_video_compare` | none: all five readings, no still | `MEAN_VMAX`, `LOG2_LIM` 1.5, 12 fps | cohort volumes | `young_vs_adult\video_side_by_side_<reading>.mp4` |
| `v2_inspect` (venv_flat) | none, then `--cmap jet`, then `--cmap turbo`. TO CONFIRM | zref, plane 790, `VMAX` 0.9, `DLIM` 0.5, smoothing 3 x 1 x 1 voxels, video and flatmaps | cohort volumes, `data\atlas_flatmap\` | `young_vs_adult\detail_*_zref.*`, and `jet\`, `turbo\` |
| `v2_diagnostics` | none (the cohort sheets and the index need no arguments) | as in the file | per-mouse files, cohort volumes, `region_table.csv`, `region_stats.csv`, `naive\nano_4d_normalized_bkgmask.mat` (CGF027), `young\nano_4d_normalized_bkgmask_P20.mat` (MG903) | `processing_diagnostics\` sheets 01 to 09 and `README.md` |

### ISH, both passes

| module | pass | settings in force | reads | writes (under `adult_v2\`) |
|---|---|---|---|---|
| `v2_ish_regions` | 1 | 200 um grid, missing = -1, `MIN_VOXELS` 3, reference box 67 x 41 x 58, one-voxel erosion | `data\gene_targets.csv`, `data\atlas_ish\` | `ish\gene_region_table.csv`, `ish\gene_region_table_drops.csv` |
| `v2_ish_compare` | 1 | `MIN_ISH_VOXELS` 10, `MIN_STRUCTURES` 50, `MIN_GENES` 20, five readings | pass-1 table, `region_means_per_mouse.csv`, P9's frozen `comparisons\merged_naive_rws_vs_ish_summary_nosmooth\gene_panel_summary.csv` | `ish\gene_correlations.csv`, `ish_old_vs_new.png` |
| `v2_ish_words` | 1 | `MIN_GENES` 5, `MAX_SHARE` 0.8, 2000 bootstraps, seed 0 | `gene_correlations.csv`, the mygene cache `ish\annotation\` (95 files) | `feature_enrichment.csv`, `ish_word_enrichment.png` |
| `v2_ish_roles` | 1 | zref, the curated roles, seed 0 | pass-1 table, `region_means_per_mouse.csv` | `gene_roles.csv`, `role_summary.csv`, `ish_roles.png` |
| `v2_panel_build` | 2 | the GO terms and the Grid1/Grid2 override as in the file | `gene_targets.csv`, the API cache `panel\cache\` (436 files) | `panel\panel_v2.csv`, `panel_genes.csv` |
| `v2_panel_fetch` | 2 | reference box 67 x 41 x 58, timeout 180 s, 2 retries | `panel_v2.csv`, `data\atlas_ish\` | missing grids; `panel\fetch_failures.csv` |
| `v2_ish_regions` | 2 | as pass 1 | `adult_v2\panel\panel_v2.csv`, `data\atlas_ish\` | `ish\gene_region_table_panel.csv`, `gene_region_table_panel_drops.csv` |
| `v2_ish_reliability` | 2 | `MIN_ISH_VOXELS` 10, `MIN_SHARED` 50, seed 0 | `gene_region_table_panel.csv` | `gene_reliability.csv`, `gene_region_table_merged.csv`, `ish_reliability.png` |
| `v2_ish_panel_test` | 2 | zref, `MIN_STRUCTURES` 80, `MIN_RELIABILITY` 0.3, 20000 permutations, seed 0 | merged table, reliability, `region_means_per_mouse.csv` | `panel_test.csv`, `ish_panel_test.png` |

### Adult arms and beyond abundance

| module | settings in force | reads | writes (under `adult_v2\`) |
|---|---|---|---|
| `v2_adult_arms` | `MIN_VOX` 250; arms sepauto, ratio, sepratio | per-mouse files, `region_means_per_mouse.csv` | `arms\region_means_arms.csv`, `arms_consistency.png` |
| `v2_ish_arms` | control Gria1, `MIN_STRUCTURES` 50, seed 0 | `region_means_arms.csv`, pass-1 table | `arms\arm_gene_correlations.csv`, `arms_vs_genes.png` |
| `v2_sep_channel_check` | `MIN_VOX` 250, seed 0 | per-mouse files, pass-1 table | `arms\sep_channel_check.csv/.png` |
| `v2_beyond_density` | zref, `MIN_VOX` 250, halves of 5, the grey-matter rule, seed 0 | `region_means_per_mouse.csv`, `gene_region_table_merged.csv` | `beyond\structures_used.csv`, `variance_partition.csv`, `residual_by_structure.csv`, `fig0` to `fig3` |
| `v2_beyond_controls` | 5 folds, up to 25 components, seed 0 | as above | `beyond\controls.csv`, `fig4` to `fig6` |
| `v2_beyond_figures` | 2000 bootstraps, 10000 permutations, 20 splits, seed 0 | as above, `controls.csv` | `beyond\for_sami\A` to `D`, `numbers_for_the_caption.txt` |
| `v2_beyond_regression` | planes 215, 265, 315; floor 0.12 | as above, per-mouse files | `beyond\for_sami\E_regression`, `F_maps`, `regression_table.csv` |

### Hidden state in the Python route

- `v2_cohort` reuses `per_mouse\<mouse>_scalars.npz` when its recorded
  modification time matches the per-mouse file. A copy that keeps modification
  times would reuse it: the caches are deleted before each run (plan).
- `v2_per_mouse` includes SEP only where `volume_registered_sep\chan02_SEP.tiff`
  exists; with it missing, `v2_cohort` stops (sepratio needs it). All 17 have
  it (23 and 24 Sep 2026).
- `v2_panel_build` and `v2_ish_words` call the network only for what their
  caches lack. `v2_panel_fetch` skips grids on disk, but retries on every run
  the two that failed in production (Gria1 526 and Negr1 695, "no energy.mhd in
  the zip"), so it always calls the network. TO CONFIRM.
- `v2_ish_compare` skips the old-against-new check without a word when P9's
  summary is missing, and `v2_adult_arms` skips its consistency check when
  `region_means_per_mouse.csv` is missing.
- `v2_diagnostics` sheet 08 reads `young\nano_4d_normalized_bkgmask_P20.mat`
  (P6bis for young_P20, 4 Sep 2026, three mice), which the reference run does
  not regenerate. TO CONFIRM.
- A figure open in a viewer is written as `*_new.png` instead.

## Expected differences from the production outputs

Not refactor errors; the output comparison maps or excludes them.

- P7bis: the `_nano` in folder and file names (above); inputs renormalised
  since 5 Dec 2025.
- P5: new `collected_mice.mat` for the three adult groups.
- P1 (MG914): `px_atlas`, `atlasaplims` and the paths inside `sliceinfo.mat`.
- P2bis: the clock time in `intensity_diagnostics\` names.
- P4 `autoannotate`: date and device in `auto_proposal_info.mat`; P4bis: the
  date line in the `.txt`.
- ISH pass 1 writes `gene_region_table_drops.csv`; production holds the older
  name `gene_drops.csv`.
- `young_vs_adult\barrel_only\` (27 Sep) has no code that writes it. TO
  CONFIRM.
- `processing_diagnostics\README.md` is regenerated.

## TO CONFIRM

1. Settings with no environment override (P5 groups and age filter, P7bis
   `exp_type`, the mouse lists and modes of P1, P2, P2bis and P4): edited in
   the worktree copy and saved as a patch beside the outputs, with the same
   values carried into the branch's drivers.
2. The two analysis environments: run from `main` (`tools\venv_atlas`,
   `tools\venv_flat`) on the worktree's scripts.
3. S6: also run today's P7bis on the normalised volumes of Nov and Dec 2025
   (`naive` and `rws` `nano_4d_normalized_bk.mat`, `behavior`'s current file),
   to separate a code difference from the renormalisation.
4. `v2_diagnostics` sheet 08: copy `young\nano_4d_normalized_bkgmask_P20.mat`
   into the tree as a frozen input (proposed), or add P5 with `age_filter [20]`
   and P6bis for `young_P20` to the reference run (it would now stack five P20
   brains instead of three).
5. `v2_panel_fetch`: run it with the network (it retries the two failures),
   or leave it out and copy `fetch_failures.csv` as a frozen input (proposed:
   leave it out; the grids it would fetch are already on disk).
6. `v2_inspect`: include the `--cmap jet` and `--cmap turbo` sets, as in
   production, or only the default set.
7. P2bis stage 1: the six mice of the production run (proposed), or MG914
   alone and accept a difference.
8. P4bis: leave `volume_aligned_sep\` and `volume_registered_sep\` out of the
   registration tree's inputs (proposed), or set `overwrite = true`.
9. P4 `autoannotate`: on the GPU, as the saved proposal was made, when no other
   process holds it.
10. `young_vs_adult\barrel_only\`: which code wrote it, and whether it is
    excluded from the comparison.
