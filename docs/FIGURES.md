# Figures

Every figure that has been shown or used so far: which run script makes it,
from which inputs, where the output lands, and which settings matter. The
results themselves are in [SCIENTIFIC_CONTEXT.md](SCIENTIFIC_CONTEXT.md); what
is still to come, in [ROADMAP.md](ROADMAP.md).

How to read this page:

- `<data>` is the data root: `get_paths().data` in MATLAB, `config.DATA` in
  the Python route, the folder `data` beside the code folder
  (`D:\sep_histology\data` on the analysis computer). `SEP_DATA_ROOT` moves
  it, for checks on a copy only.
- Scripts are named as after the refactor; the name a figure was made under,
  when different, is in brackets. The full table of old and new names is
  [refactor_name_map.csv](refactor_name_map.csv).
- Python run scripts run in `tools\venv_atlas`, except `run_closeup.py`, which
  runs in `tools\venv_flat`. Their run order is in the header of every
  `mapping/run_*.py`. Their settings are in `mapping/settings.toml`, by table
  (`[closeup]`, `[videos]` and so on); each run prints the settings in force.
- MATLAB drivers set their settings in their `%% Settings` block; they run
  after `sep_setup_paths`, once per session.
- Code states: the tag `grant-2026-09` is the code behind the grant figures;
  `refactor-start` is the code before the refactor. Outputs carry their file
  dates, which tell which code made them.
- A figure open in a viewer while a Python run saves it is written as
  `<name>_new.png` instead.
- "Changed since" lists what a later fix changes in a figure made before
  3 October 2026. A rerun with the merged code gives the new version
  ([ROADMAP.md](ROADMAP.md), section 2).

---

## The grant: SNSF Weave, Figure 3

Panels b to d show young against adult (line 4), on the readings of the
Python route as of 25 to 27 September 2026, code at `grant-2026-09`. The
young cohort is the seven brains of `young` (five at P20, MG911 at P16,
MG904 at P22), the adults the ten of `adult` (five naive, five RWS). The
grant's caption calls the pups "aged between P16 and P20": MG904 is P22,
grouped with the P20 brains as Sami El-Boustani asked.

| panel | what it shows | made by | output, under `<data>\comparisons_v2\young_vs_adult\` |
|---|---|---|---|
| 3a | surface GluA1, SEP fluorescence and Hoechst in one section, with a confocal close-up | not made by this code | |
| 3b | young and adult maps at the level of RL and AL | `mapping/run_closeup.py` (`v2_inspect.py`) | `detail_plane790_zref.png/.eps` |
| 3c | the young-minus-adult difference on the flattened cortex | `mapping/run_closeup.py` | `detail_flatmap_layers_zref.png/.eps`, one panel per depth band (which band, to confirm, below); the full-depth version is `detail_flatmap_zref.png` |
| 3d | V1, RL+AL, the other higher visual areas, S1, S2 and medial prefrontal cortex, a dot per mouse, Mann-Whitney stars | `mapping/run_region_groups.py` (`v2_region_groups.py`) | `group_plot.png/.eps`, values in `group_stats.csv`; or the barrel-field variant in `barrel_only\` (to confirm, below) |

**Inputs.** For 3b and 3c, the cohort volumes `<data>\comparisons_v2\ccf\young\`
and `ccf\adult\` (`run_cohort.py`, after `run_per_mouse.py` and
`run_to_ccf.py`), the CCF annotation, and for the flatmaps the
ccf_streamlines assets in `<data>\atlas_flatmap\` (listed in
`sepmap/young_vs_adult/closeup.py`). For 3d, the per-brain files of
`<data>\comparisons_v2\per_mouse\` measured on each brain's own atlas, with no
warping.

**Settings that matter.**

- Reading `zref` for the maps; the plane, 790 on the 10 um CCF grid
  (`[closeup] plane`, or `--plane`); the colour ranges `zref` 0.9 for the means
  and 0.5 for the difference (`[closeup] vmax`, `dlim`, or `--vmax`,
  `--dlim`); a light smoothing of 3 x 1 x 1 voxels of 20 um along AP, DV and
  ML (`[closeup] smooth`, `--smooth`), declared in every title and invisible
  to the statistics.
- The mean panels use `PuOr_r` for `zref`, the difference `RdBu_r`. `--cmap
  jet` and `--cmap turbo` (asked for by Sami El-Boustani) draw the mean panels
  with that colormap into the subfolders `jet\` and `turbo\`; the difference
  stays blue-red.
- The depth bands of the layer flatmap (`BANDS` in `closeup.py`): layers 1 and
  2/3, layer 4, layers 5 and 6.
- 3d: the systems (`SYSTEMS` in `sepmap/young_vs_adult/region_groups.py`,
  each naming its members), at least 3 young and 5 adult brains with a value
  (`[region_groups]`), the uncorrected Mann-Whitney p for the stars (the
  column `mannwhitney_p`; BH q in `mannwhitney_q_BH`).
- The outlines of RL and AL in green in the grant are not drawn by the code
  (its flatmaps outline VISp, VISrl, VISal and SSp-bfd brighter).

**To confirm.** Which colormap set (default, `jet\` or `turbo\`) the grant's
3b and 3c used. Which depth band of the layer flatmap 3c shows: the caption
names no layer. Whether 3d's S1 is all of primary somatosensory cortex
(`group_plot.png`) or the barrel field alone (`barrel_only\group_plot.png`,
27 September, written by `D:\sep_histology\sandbox_barrel\group_plot_barrel_only.py`,
a script outside this repository that reruns `v2_region_groups` with S1
narrowed to SSp-bfd). Which reading 3d shows (presumably `zref`, as the maps).
The grant's "mPFC": the nearest group in the code is `frontal` (ACA, ORB, PL,
ILA, FRP, DP).

**Changed since.** A1 changes the reference of `zref` and so every value in
these panels. The fixes of 3 October leave the `zref` close-ups unchanged; in
`group_plot.png` and `group_stats.csv` they change the `subref` panel (its
reference now leaves out fibre tracts and ventricles) and, slightly, the
`sepratio` panel (MG897's missing voxels no longer count as zero).

## The plasticity comparison Sami El-Boustani approved (5 December 2025)

Naive against RWS, and naive against behaviour (line 1). The approved
outputs are in `<data>\comparisons\naive_vs_rws\` and
`<data>\comparisons\naive_vs_behavior\` (a copy on the lab share is identical,
file by file). Its headline is the increase of the nanobody signal in S1 after
RWS: the hemisphere-sum t map of the slab, masked by surprise at p < 0.01.

| figure | what it shows |
|---|---|
| `Slab_Avg_565_<comp>_surpmask.png/.fig` | the group t maps of the left-right difference and sum, median over the slab, opaque where the surprise reaches p < 0.01 |
| `Indiv_Slab_Avg_565_<group>.png/.fig` | every mouse's left-right difference and sum in the same slab |
| `Region_Surprise_Bar_DiffSum_<comp>.png/.fig` | the surprise summed per region, for the difference and the sum |
| `Normalization_Profiles_LR_<comp>.png/.fig` | the experimental group's plane profile aligned onto the control group's |
| `lr_diff_sum_nano_<group>.mp4`, `lr_diff_sum_nano_groupdiff_<comp>.mp4`, `t_*`, `surp_*`, `*_surpmask*.mp4`, `Individual_*_<group>.mp4` | the same, plane by plane |

`<comp>` is `naive_vs_rws` or `naive_vs_behavior`; the `.fig` files hold the
values (no `.mat` is written).

**Made by** `group_comparison/run_group_differences.m`
(`P7bis_analyze_group_differences.m`), once per comparison, after
`run_collect_by_group.m` (`P5`) and `run_normalise_groups.m` (`P6bis`) for the
three groups. The code of 5 December 2025 is not tagged; on 1 October 2026
today's code reproduced these figures from the same inputs (Progress in
[REFACTOR_PLAN.md](REFACTOR_PLAN.md)).

**Inputs.** `<data>\<group>\nano_4d_normalized.mat` and
`nano_4d_normalized_bkgmask.mat` of naive and of the experimental group. The
approved run read the normalised volumes of 26 and 27 November 2025 (naive and
rws, kept as `nano_4d_normalized_bk.mat`) and behaviour's of 5 December 2025.
Their mice: naive CGF027, CGF028, CGF033, CGF034, CGF035; RWS MG691, MG692,
MG693, MG736, MG737; behaviour MG705, MG709, MG716, MG718.

**Settings that matter.**

- `run_group_differences`: `ctrl_type = 'naive'`, `exp_type` (`'rws'` or
  `'behavior'`), `behavior_subset` (all four behaviour mice saved),
  `channel = 'nano'`, `apply_smoothing = true` with `smooth_sigma = 5` (a 3D
  Gaussian over each mouse's tissue, in 10 um voxels, after which the voxels
  outside the tissue are set to 0; fix 23 leaves them out), every video flag
  on, both region analyses off.
- `run_normalise_groups`: the mice kept per group (`selected_mice_idx_list`:
  rws and naive all five, behaviour four of seven), and `cohort_specs`, or
  `SEP_COHORT_SPECS=rws,naive,behavior` so behaviour is normalised by the same
  code.
- Fixed in the code: the slab is plane 565 of the cropped volume and 10 planes
  on each side (the crop starts at plane 180 of the 10 um CCF annotation); the
  experimental group is aligned onto the control group by a line between their
  mean plane profiles over planes 200 to 700, then one shared scale from planes
  300 to 500; t limits -6 to 6; the regional bars sum the surprise above
  p < 0.01 over a rolling median of plus or minus 10 planes, on the left
  hemisphere, over a fixed list of regions.

**Names.** Today's code puts the channel in the comparison tag, so it writes
`naive_vs_rws_nano\` with `_nano` after the tag in every name that carries it
(`Slab_Avg_565_naive_vs_rws_nano_surpmask`). It never overwrites the approved
folders.

**Changed since.** After smoothing, the voxels outside the tissue are left
out instead of set to 0 (fix 23, accepted on 3 October): they no longer
enter the group means and the SEM's n as zeros, which changes the values of
every figure. In the regional bars, the parafascicular nucleus was listed
twice, and the bar labelled mediodorsal nucleus summed the intermediodorsal
nucleus, whose name contains the one asked for; the duplicate goes, and the
bar sums the mediodorsal nucleus. "Hemishpere" becomes "Hemisphere" in the
individual videos' titles.

## The adult map

How the nanobody signal is distributed across the adult brain (line 2): the
ten naive and RWS adults pooled.

### Python route

| output | made by | what it shows |
|---|---|---|
| `<data>\comparisons_v2\ccf\adult\<reading>_mean.npy`, `_sd.npy`, `_n.npy` | `mapping/run_cohort.py` (`v2_cohort.py`) | the cohort mean, SD and brain count per voxel, on the CCF at 20 um |
| `<data>\comparisons_v2\ccf\adult\video_<reading>_adult.mp4` | `mapping/run_video.py` (`v2_video.py`) | the mean and its reliability t, plane by plane |
| `<data>\comparisons_v2\young_vs_adult\region_means_per_mouse.csv` | `mapping/run_region_plot.py` (`v2_region_plot.py`) | one value per reading, mouse and structure; the table the adult and ISH analyses read |

The same videos exist for `naive` and `rws` alone. Settings: the readings
(`[readings] in_force`, or `V2_READINGS`); the tissue mask (`[tissue]`); the
smallest structure, 250 voxels of 20 um (`[region_tables] min_vox20`); for the
videos the colour range per reading (`[videos] mean_vmax`), the t panel's
upper end (`t_pct`) and the brains a voxel needs (`min_n`).

**Changed since.** The videos' reliability t will count each mouse once over
its two hemispheres (fix of 3 October; today it takes the larger of the two
hemispheres' counts). The `subref` volumes and video change with the `subref`
fix. A1 changes every `zref` value.

### P8 and P10, until A4 and A5 replace them

The bar charts shown to Sami El-Boustani in April and May 2026, rerun on
4 September 2026, are in `<data>\comparisons\merged_naive_rws_nano\`:

- `Region_MeanSum_BarByMacro_merged_naive_rws_nano_nosmooth_withSEM.png` and
  the `Zscore` version: every region by division, per-mouse mean and SEM,
  bars shaded by reliability; `BarAcrossDivi` the divisions; without
  `_withSEM` the cohort means alone; `Region_MeanSum_Table_*.csv` the values;
- `lr_sum_nano_*.mp4`, `lr_sum_zscore_*.mp4`, `lr_sum_threshold_*.mp4`: the
  pooled map plane by plane, raw, z-scored, and with the enrichment contour.

Made by `P8_characterize_merged_distribution.m` from the normalised volumes
of naive and rws (`run_normalise_groups`). Settings: `groups_to_merge =
{'naive', 'rws'}`, `channel`, `apply_smoothing = false`, `region_agg_method =
'distweight'` with `dist_weight_power = 4`, planes 100 to 700, the enrichment
thresholds (`raw_threshold = 1.5`, `zscore_threshold = 0`),
`compute_per_mouse_sem = true`. The same run with `channel = 'auto'` wrote
`merged_naive_rws_auto\`.

The autofluorescence control, `<data>\comparisons\nano_vs_auto\` (May 2026):
`Region_NanoVsAuto_BarByMacro_*`, `Region_NanoVsAuto_DeltaZ_BarByMacro_*`,
`Macro_NanoVsAuto_Paired_*`, `Macro_NanoVsAuto_DeltaZ_*`, their signed-rank
tables, and `Contrast_video_NanoMinusAuto_z_nosmooth.mp4`. Made by
`P10_compare_nano_vs_auto.m` from P8's per-mouse caches of both channels.
Settings: `agg_method = 'distweight'`, `apply_bonferroni = false`,
`alpha = 0.05`, one-sided paired signed-rank tests.

These outputs carry the known defects of the retiring route
([ROADMAP.md](ROADMAP.md), section 4): the slab artefact of P8's `abs()`, the
stale autofluorescence of P8's `auto` run and of P10, reliability bars that
start at white. They are not quoted.

## What the map measures

### Beyond abundance and density

The figures prepared for Sami El-Boustani on 26 September 2026, in
`<data>\adult_v2\beyond\for_sami\`:

| figure | made by | what it shows |
|---|---|---|
| `A_what_explains.png/.eps` | `mapping/run_beyond_figures.py` (`v2_beyond_figures.py`) | the map's explainable variance and what each explanation predicts of it |
| `B_leftover_real.png/.eps` | `run_beyond_figures.py` | the leftover replicating across independent halves of the cohort, against a noise null |
| `C_where.png/.eps` | `run_beyond_figures.py` | the structures where the map most exceeds and falls short of the prediction |
| `D_controls.png/.eps` | `run_beyond_figures.py` | the seven controls |
| `numbers_for_the_caption.txt` | `run_beyond_figures.py` | every number of A to D as a sentence |
| `E_regression.png/.eps`, `F_maps.png/.eps`, `regression_table.csv` | `mapping/run_beyond_regression.py` (`v2_beyond_regression.py`) | observed against predicted, the residual diagnostic, and observed, predicted and residual on three coronal planes |

Their working figures are in `<data>\adult_v2\beyond\`: `fig0_structures` to
`fig3_residual` and the tables `structures_used.csv`, `variance_partition.csv`,
`residual_by_structure.csv` (`mapping/run_beyond_density.py`), and
`fig4_controls` to `fig6_readings` with `controls.csv`
(`mapping/run_beyond_controls.py`).

**Inputs.** `region_means_per_mouse.csv` (`run_region_plot`) and
`<data>\adult_v2\ish\gene_region_table_merged.csv` (`run_ish_reliability`, the
390-gene pass); panel D reads `controls.csv`, so `run_beyond_controls` comes
before `run_beyond_figures`; panel F also reads the per-brain files.

**Settings that matter.** `[beyond]`: the reading (`zref`), animals per half
(5), the grey-matter divisions kept, the gene sets standing for abundance
(Gria1 to Gria4) and synaptic density (the markers); `[beyond_controls]`, each
control's threshold; `[beyond_figures]`, bootstrap and permutation counts;
`[beyond_regression]`, the planes drawn (215, 265 and 315 of the cropped 20 um
grid) and the colour floor. Random seeds are fixed (0).

**Changed since.** Panel F's rows are labelled with their CCF planes (they
read 4.3, 5.3 and 6.3 mm). `fig5_model_space` says how many genes the model
used (it said all 390). Panel A's EPS draws the ceiling band light grey under
the bars (it was solid over them). The printed conclusions are reworded to
follow the numbers. A1 to A3 do not move this result.

### The ISH comparison

No presentation of the ISH work exists. The outputs of record are those of
P9's final run, and the Python route's.

**P9 (21 and 22 April 2026), until A1 to A3 replace it.**
`<data>\comparisons\merged_naive_rws_vs_ish_summary_nosmooth\`:
`correlation_barchart_4metrics.png`, `violin_spearman_dw.png`,
`violin_spearman_ero.png`, `violin_voxel_pearson.png`,
`correlation_eroded_vs_distweight.png`, and the ranking `gene_panel_summary.csv`;
one folder per gene, `merged_naive_rws_vs_ish_<gene>_nosmooth\`, with
`Scatter_NanoVsISH_<gene>.png`, the paired bars, the comparison video and
two diagnostic sheets. Made by `P9_compare_nano_vs_allen_ish.m` from P8's
cohort cache, `<data>\gene_targets.csv` (the 100-gene panel) and the Allen
grids in `<data>\atlas_ish\`. Settings: the nine divisions, erosion radius 3,
distance-weight power 4, planes 100 to 700. These folders have no channel in
their names: they predate it, and today's P9 would write
`merged_naive_rws_nano_vs_ish_*`. They carry P9's defects (the stretched
grid, the section repair).

**Python route (25 and 26 September 2026).** In `<data>\adult_v2\ish\`
unless stated:

| figure | made by | what it shows |
|---|---|---|
| `ish_old_vs_new.png` | `mapping/run_ish_compare.py` | each gene's rho against P9's (reads P9's `gene_panel_summary.csv`) |
| `ish_word_enrichment.png` | `mapping/run_ish_words.py` | annotation words and GO terms at the top of the ranking |
| `ish_roles.png` | `mapping/run_ish_roles.py` | subunit against localisation genes, with the exact permutation |
| `ish_reliability.png` | `mapping/run_ish_reliability.py --panel ontology` | how reliable one Allen ISH map is, from genes measured twice |
| `ish_panel_test.png` | `mapping/run_ish_panel_test.py` | localisation genes against expression-matched controls, and the positive control |
| `adult_v2\arms\arms_vs_genes.png` | `mapping/run_ish_arms.py` | the channel arms against the genes, plain and with Gria1 partialled out |
| `adult_v2\arms\arms_consistency.png` | `mapping/run_adult_arms.py` | the arms' self-check against `run_region_plot`'s table |

The tables beside them hold the values (`gene_correlations.csv`,
`feature_enrichment.csv`, `role_summary.csv`, `gene_reliability.csv`,
`panel_test.csv`, `arm_gene_correlations.csv`).

Run order: the 100-gene pass (`run_ish_regions.py --panel targets`, then
compare, words, roles, the arms), then the 390-gene pass (`run_panel_build`,
`run_panel_fetch`, which needs the network, `run_ish_regions.py --panel
ontology`, `run_ish_reliability`, `run_ish_panel_test`). Settings: the
reading (`[ish] reading`, `zref`) and the gene for total receptor
(`control_gene`, Gria1); `[ish_words]`, `[ish_panel_test]` (20,000
permutations, the reliability a gene needs); `PYTHONHASHSEED` for
`run_ish_words` until its fix. The structure set is today every structure the
cohort measures, which is why S5 holds the numbers back until A1 to A3.

**Changed since.** The words' bootstrap intervals no longer depend on the
order of Python sets; each permutation test of the panel test gets its own
random generator (its p moves within Monte Carlo error); the arms' self-check
gets a margin for floating point and draws its figure before it stops; the
table of dropped experiments gives the real reason for Gria1 and Negr1 (no
`energy.mhd` in the downloaded file).

### What the green channel reports

`<data>\adult_v2\arms\sep_channel_check.png` and `sep_channel_check.csv`,
made by `mapping/run_sep_channel_check.py` (`v2_sep_channel_check.py`) from the
per-brain files of the ten adults and the 100-gene table: each channel's
dynamic range, what each tracks across structures, and the SEP residual once
autofluorescence is regressed out. The figure behind the finding that the
green channel is mostly autofluorescence.

## Young against adult: figures and videos

All in `<data>\comparisons_v2\young_vs_adult\` unless stated, from the run of
24 and 25 September 2026 (seven young brains, ten adults), every reading in
force (`ratio`, `sepratio`, `cref`, `subref`, `zref`).

| output | made by | what it shows |
|---|---|---|
| `slices_<reading>.png/.eps` | `mapping/run_compare.py` (`v2_compare.py`); `run_replot.py` redraws them from `volumes_ccf20.npz` | adult, young and their difference on six coronal planes, dorsal up |
| `region_table.csv`, `cortex_table.txt`, `volumes_ccf20.npz` | `run_compare.py` | per structure from the voxel maps; everything the slice figures are drawn from |
| `region_plot.png/.eps`, `region_stats.csv` | `mapping/run_region_plot.py` | a dot per mouse per structure, and the tests: the numbers to quote |
| `group_plot.png/.eps`, `laminar_plot.png/.eps`, `group_stats.csv` | `mapping/run_region_groups.py` | by system and by layer within each cortical system |
| `barrel_only\group_plot.png`, `laminar_plot.png`, `group_stats.csv` | `sandbox_barrel\group_plot_barrel_only.py`, outside the repository | the same with S1 narrowed to the barrel field |
| `<data>\comparisons_v2\ccf\<cohort>\video_<reading>_<cohort>.mp4` | `mapping/run_video.py` | a cohort's mean and its reliability t, plane by plane (young, adult, young_P20, naive, rws) |
| `video_side_by_side_<reading>.mp4` | `mapping/run_video_compare.py` | young, adult and their comparison, plane by plane |
| `detail_plane790_zref`, `detail_video_zref.mp4`, `detail_flatmap_zref`, `detail_flatmap_layers_zref` (and in `jet\`, `turbo\`) | `mapping/run_closeup.py` | `zref` close up: one plane, its video, the flattened cortex through its depth and by band |

**Inputs.** The maps and videos read the cohort volumes of `run_cohort`,
which reads the brains carried to the CCF by `run_to_ccf` (young brains
warped from the DeMBA atlas of their age, adults only placed). The tables
(`region_plot`, `region_groups`) read the per-brain files of `run_per_mouse`
on each brain's own atlas, with no warping. Quote numbers from the tables and
use the maps for the pattern: maps and tables do not give the same number
(`<data>\comparisons_v2\README.md`, "Reading the numbers off the figures").

**Settings that matter.**

- The cohorts (`COHORTS` in `sepmap/volumes/cohort.py`, then the shared
  cohort table `common\cohort.csv` of decision Y4): `young` pools P16, P20
  and P22.
- `[young_vs_adult]`: brains a voxel needs to be compared (2 young, 5 adult),
  the smoothing of the log2 map (1 voxel of 20 um); `[region_plot]` and
  `[region_groups]`: brains a structure or group needs to be tested.
- `zref`'s reference: today the structures shared by all seventeen brains for
  the tables, every structure of the brain for the maps (A1 makes them one
  set).
- `[videos]`: frames per second, the colour range per reading (`mean_vmax`),
  `log2_lim` for the side-by-side comparison, `min_n` per cohort.
  `run_video_compare.py --plane P` writes a still
  (`plane<P>_side_by_side_<reading>.png`) instead of the video; `--vmax V`
  adds `_vmax<V>` to the name, `--dlim` changes nothing in the name.
- `run_closeup.py`: as for the grant panels above.

**Changed since.**

- `slices_<reading>.png`: the panel titles gave plane numbers 180 too high
  (526 to 1074 for 346 to 894); the maps were right.
- `subref` everywhere: its reference now leaves out fibre tracts and
  ventricles (and the unassigned labels).
- `sepratio` young means: MG897's missing voxels no longer count as zero in
  the cohort mean, and the region sums become NaN-aware.
- `region_plot.png`: the `ratio` panel's title says what its zero means.
- The videos' reliability t: each mouse counted once over both hemispheres.
- A1 moves every `zref` value; the laminar contrast per mouse is added to
  `run_region_groups`' outputs (step 9).

## Diagnostic sheets

### Python route

`<data>\comparisons_v2\processing_diagnostics\`, made by
`mapping/run_diagnostics.py` (`v2_diagnostics.py`), with a `README.md` index
saying what to look for in each sheet:

| sheet | per | checks |
|---|---|---|
| `01_tissue_<mouse>.png` | brain | what was counted as tissue |
| `02_levels_<mouse>.png` | brain | where the background and the threshold sit |
| `03_coverage.png` | cohort | the planes each brain covers |
| `04_warp_<mouse>.png` | young brain | before and after the DeMBA-to-CCF warp |
| `05_cohort_n.png` | cohort | how many brains behind each voxel |
| `06_scaling.png` | cohort | backgrounds and cortex means across brains |
| `07_route_agreement.png` | cohort | the warped voxel route against the unwarped region route |
| `08_mask_vs_p6bis.png` | two brains | the tissue mask against the plasticity chain's background mask |
| `09_denominators.png` | cohort | the two reference channels against each other |

Inputs: the per-brain files, the cohort volumes, `region_table.csv` and
`region_stats.csv`; sheet 08 reads the background masks of
`run_normalise_groups` (`naive\nano_4d_normalized_bkgmask.mat` for CGF027,
`young\nano_4d_normalized_bkgmask_P20.mat` for MG903), so it runs after the
plasticity chain. `run_diagnostics.py <mouse>` refreshes one brain's sheets.
Settings: `[tissue]` (the mask's threshold, `mad_k`, appears in the titles).
Changed since: sheet 01's title says what the tissue mask needs.

`processing_diagnostics\sep_channel\<mouse>.png` and `.txt`, one per brain,
are written by `registration/run_add_sep_channel.m` (`P4bis_add_sep_channel.m`):
every slice of the SEP channel recovered against the aligned DAPI, and the
registered DAPI of that run on top of the one `register` wrote. Settings:
`recovery_levels`, `min_slice_corr` (0.99). Changed since: the panel title
names `run_register_to_atlas` without TeX subscripts.

### Preprocessing, registration and normalisation (MATLAB)

Under each brain's folder, `<data>\<group>\<mouse>\lightsuite\`:

- `correction_output\diagnostic_plots\reference_pix_analysis_slice_<n>.png`
  and `reference_pix_selection_slice_J_<n>.png`, the fit of nano on
  autofluorescence per slice, and `correction_output\scaled_difference_video_<type>.mp4`:
  `preprocessing/run_residual_correction.m` (`P2`). Settings: `doPlotBkg`,
  `savePlotBkg`, `saveRatioMap`. A slice with fewer than two reference pixels
  draws its overlay with the previous slice's slope (question 22,
  [ROADMAP.md](ROADMAP.md)).
- `slice_order_montage.png`, every slice as a tile in its curated order:
  the slice-order editor of `preprocessing/run_order_slices.m` (`P1bis`),
  saved next to the decisions file. The README's image is one of these.
- `<mouse>_dim<d>_initial_registration.png`, the initial 3D fit to the atlas:
  `registration/run_register_to_atlas.m`, mode `align`.

Across brains:

- `<data>\intensity_diagnostics\`, the slice-to-slice equalisation of the nano
  channel (`Plot_*_<time>.png`, `Video_<mouse>_<time>.mp4`, the statistics):
  `preprocessing/run_nano_equalisation.m` (`P2bis`); the names carry the run's
  clock time.
- `<data>\<group>\Background_mask_diagnostics_trace_<group>_<channel>.png` and
  `<group>\global_diagnostics\normalization_checks_<channel>\`:
  `group_comparison/run_normalise_groups.m`. Changed since: the background
  trace takes its y range from every mouse, not the last one.
- `<data>\young\registration_qc\` (`ATLAS_PARAMETERS.md`,
  `atlas_crop_and_ap_mapping.png`, `atlas_comparison_adult_vs_p20.png`,
  `atlas_region_comparison.png`, `demba_to_allen_transform_qc.png`): the atlas
  checks in `atlas/qc/` (`atlas_diagnostics`, `compare_atlases_montage`,
  `compare_atlas_regions`, `check_demba_to_allen`), run by hand.
  `registration/qc/check_registration_error.m` prints each aligned brain's
  fit error.

## Superseded outputs, not to quote

- `<data>\comparisons\young_P20_vs_adult_nano\`: the first young-against-adult
  attempt (4 September 2026, three P20 brains), made by five Python scripts now
  in `archive/` (`compare_young_vs_adult_lrsum.py` and the others); replaced by
  the Python route.
- `<data>\comparisons\merged_young_P20_nano\`: P8 run on the P20 brains
  (`SEP_MERGE_SPECS=young_P20`); replaced by the Python route, and later by
  A4's young cohort.
- `<data>\comparisons\merged_naive_rws\`: P8's earlier nano run, before the
  channel entered its folder names; `merged_naive_rws_nano\` is the later one.
- The `Diagnostic_DistWeight_*` sheets in P8's folders: dropped by decision
  (29 September), not rebuilt.
