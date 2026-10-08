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
- "Changed since" lists what a later fix or decision changed in a figure made
  before 3 October 2026. The rerun of 5 October 2026 wrote the new versions
  at the same paths; the versions shown before are kept in
  `<data>\backup_before_rerun_2026-10-05\` for the Python route, and in the
  approved folders for the plasticity comparison ([ROADMAP.md](ROADMAP.md),
  section 2).

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
these panels. The fixes of 3 October leave the `zref` close-ups and every
`zref` and `cref` value unchanged; in `group_plot.png` and `group_stats.csv`
they change the `subref` panel and its title (the reference is now TH, HY,
PAL, MB, P and MY: fibre tracts, ventricles and unassigned labels left out).
MG897's missing `sepratio` voxels, now left out of its region means, move one
marker of `group_plot.eps` by a thousandth of a point; the PNG is unchanged.
Decision 4 of 4 October carries the young brains to the CCF without
darkening their tissue edges, so the young cortical surface reads brighter in
the close-ups and flatmaps (the young outer cortical shell, 0.548 to 0.599 of
the interior), and takes each brain's backgrounds over its imaged voxels
only, which moves the differences of `group_stats.csv` by at most 0.01
(RL+AL `zref` +0.232 to +0.235).

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
| `Region_Surprise_Bar_DiffSum_<comp>.png/.fig` | the surprise summed per region, for the difference and the sum (since 4 October 2026 the fraction of each region's voxels with a t at p < 0.01, the sum in `Region_Surprise_DiffSum_<comp>.csv`) |
| `Region_Surprise_Bar_<measure>_<comp>.png/.fig` | since 6 October 2026: the regions scored by `bar_measure` (the mass of the heaviest cluster at p < 0.01, Giulio's choice of 7 October 2026), signed, shaded by the corrected p of an exact label permutation, starred below 0.05: the test of a search over all 71 regions. Since 7 October 2026 also a dagger where only the uncorrected p is below 0.05, the uncorrected p in bold beside the bar of each region named in advance (`a_priori_regions`, the barrel field), whose test it is, and two lines under the title saying so; a dagger on any other region is no test. Names in magenta: the regions expected to change |
| `Region_Measures_<comp>.png/.fig` | since 6 October 2026: the regions in the first ten of any of the five measures (share, sum, q99, topvol, cluster), with their rank and corrected p under each; the values in `Region_Surprise_DiffSum_<comp>.csv`, the permutation null in `Region_Permutation_Null_<comp>.mat` |
| `Normalization_Profiles_LR_<comp>.png/.fig` | the experimental group's plane profile aligned onto the control group's |
| `lr_diff_sum_nano_<group>.mp4`, `lr_diff_sum_nano_groupdiff_<comp>.mp4`, `t_*`, `surp_*`, `*_surpmask*.mp4`, `Individual_*_<group>.mp4` | the same, plane by plane |

`<comp>` is `naive_vs_rws` or `naive_vs_behavior`; the `.fig` files hold the
values (no `.mat` is written).

**Made by** `group_comparison/run_group_differences.m`
(`P7bis_analyze_group_differences.m`), once per comparison, after
`run_collect_by_group.m` (`P5`) and `run_normalise_groups.m` (`P6bis`) for the
three groups. The code of 5 December 2025 is not tagged; on 1 October 2026
today's code reproduced these figures from the same inputs (Progress in
[history/REFACTOR_PLAN.md](history/REFACTOR_PLAN.md)).

**Inputs.** `<data>\<group>\nano_4d_normalized.mat` and
`nano_4d_normalized_bkgmask.mat` of naive and of the experimental group. The
approved run read the normalised volumes of 26 and 27 November 2025 (naive and
rws, kept as `nano_4d_normalized_bk.mat`) and behaviour's of 5 December 2025.
Their mice: naive CGF027, CGF028, CGF033, CGF034, CGF035; RWS MG691, MG692,
MG693, MG736, MG737; behaviour MG705, MG709, MG716, MG718.

**Settings that matter.**

- `run_group_differences`: `ctrl_type = 'naive'`, `exp_type` (`'rws'` or
  `'behavior'`), `behavior_mice` (the four behaviour mice saved, by name;
  `behavior_subset`, by position, before step 8), `channel = 'nano'`,
  `apply_smoothing = true` with `smooth_sigma = 5` (a 3D Gaussian over each
  mouse's tissue, in 10 um voxels, normalised by the smoothed tissue mask;
  outside the tissue a mouse has no value since fix 23, where the approved run
  set those voxels to 0), every video flag on, both region analyses off.
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
folders. The merged code wrote `naive_vs_rws_nano\` and
`naive_vs_behavior_nano\` on 5 October 2026, from production's current
normalised volumes, which are not the approved run's inputs
([ROADMAP.md](ROADMAP.md), section 2).

**Changed since.** Each mouse's voxels outside its tissue are left out
instead of set to 0 after smoothing (fix 23, accepted on 3 October): every
mean, SEM and t is taken over the mice with tissue at that voxel, and the t
and surprise videos show only the voxels that have a t. Values inside the
tissue are unchanged; the t maps, the surprise and the bars change wherever
some mouse lacks tissue, and in the individual slab figures 0.2 to 0.5% of
the pixels shown differ. Naive against behaviour changes over the whole slab,
since MG709 has no tissue there and no longer counts as a fourth mouse. The
S1 result before and after, on the approved run's inputs, is in
[ROADMAP.md](ROADMAP.md), section 1. In the regional bars, the
parafascicular nucleus was listed twice, and the bar labelled mediodorsal
nucleus summed the intermediodorsal nucleus (88,205 voxels), whose name
contains the one asked for; the duplicate goes, and the bar sums the
mediodorsal nucleus (691,695 voxels). "Hemishpere" becomes "Hemisphere" in
the individual videos' titles. Since 4 October a voxel gets a t only where
each group has at least three mice, the slab figures and bars use only the
voxels with a t of their own, and each regional bar is the share of one of 71
atlas regions' voxels at p < 0.01, no voxel in two bars (decisions 1 and 2,
[ROADMAP.md](ROADMAP.md), section 1).

## The adult map

How the nanobody signal is distributed across the adult brain (line 2): the
ten naive and RWS adults pooled.

### Python route

| output | made by | what it shows |
|---|---|---|
| `<data>\comparisons_v2\ccf\adult\<reading>_mean.npy`, `_sd.npy`, `_n.npy` | `mapping/run_cohort.py` (`v2_cohort.py`) | the cohort mean, SD and brain count per voxel, on the CCF at 20 um; `<reading>_folded_mean.npy`, `_sd`, `_n` the same with each brain's hemispheres averaged first, on the left half |
| `<data>\comparisons_v2\ccf\adult\video_<reading>_adult.mp4` | `mapping/run_video.py` (`v2_video.py`) | the mean and its reliability t, plane by plane |
| `<data>\comparisons_v2\young_vs_adult\region_means_per_mouse.csv` | `mapping/run_region_plot.py` (`v2_region_plot.py`) | one value per reading, mouse and structure; the table the adult and ISH analyses read |

The same videos exist for `naive` and `rws` alone. Settings: the readings
(`[readings] in_force`, or `V2_READINGS`); the tissue mask (`[tissue]`); the
smallest structure, 250 voxels of 20 um (`[region_tables] min_vox20`); for the
videos the colour range per reading (`[videos] mean_vmax`), the t panel's
upper end (`t_pct`) and the brains a voxel needs (`min_n`).

**Changed since.** The videos' reliability t is taken over each brain's two
hemispheres averaged first, one value per brain, from the `_folded` files
(fix 24; before, the SEM's n was the larger of the two hemispheres' counts).
The absolute t rises in 89 to 96% of the voxels shown, by a median factor of
1.10 to 1.14 for `cref` and `zref` (1.02 to 1.06 for `ratio`, 1.00 to 1.04 for
`sepratio`, 1.09 to 1.11 for `subref`), and the t panel's upper end with it
(adult `cref` 22.9 to 28.0, `zref` 14.8 to 16.9); the mean panel is
unchanged. The `subref` volumes and video change with the `subref` fix.
Decision 4 takes each brain's backgrounds over its imaged voxels only, which
moves the readings slightly everywhere. A1 changes every `zref` value.

### The MATLAB distribution and autofluorescence control, until A4 and A5 replace them

The bar charts shown to Sami El-Boustani in April and May 2026, rerun on
4 September 2026, are in `<data>\comparisons\merged_naive_rws_nano\`:

- `Region_MeanSum_BarByMacro_merged_naive_rws_nano_nosmooth_withSEM.png` and
  the `Zscore` version: every region by division, per-mouse mean and SEM,
  bars shaded by reliability; `BarAcrossDivi` the divisions; without
  `_withSEM` the cohort means alone; `Region_MeanSum_Table_*.csv` the values;
- `lr_sum_nano_*.mp4`, `lr_sum_zscore_*.mp4`, `lr_sum_threshold_*.mp4`: the
  pooled map plane by plane, raw, z-scored, and with the enrichment contour.

Made by `adult_matlab/run_characterize_distribution.m`
(`P8_characterize_merged_distribution.m`) from the normalised volumes of
naive and rws (`run_normalise_groups`). Settings: `groups_to_merge =
{'naive', 'rws'}`, `channel`, `apply_smoothing = false`, `region_agg_method =
'distweight'` with `dist_weight_power = 4`, planes 100 to 700, the enrichment
thresholds (`raw_threshold = 1.5`, `zscore_threshold = 0`),
`compute_per_mouse_sem = true`. The same run with `channel = 'auto'` wrote
`merged_naive_rws_auto\`.

The autofluorescence control, `<data>\comparisons\nano_vs_auto\` (May 2026):
`Region_NanoVsAuto_BarByMacro_*`, `Region_NanoVsAuto_DeltaZ_BarByMacro_*`,
`Macro_NanoVsAuto_Paired_*`, `Macro_NanoVsAuto_DeltaZ_*`, their signed-rank
tables, and `Contrast_video_NanoMinusAuto_z_nosmooth.mp4`. Made by
`adult_matlab/run_compare_nano_with_autofluorescence.m`
(`P10_compare_nano_vs_auto.m`) from the per-mouse caches of both channels
that `run_characterize_distribution` wrote.
Settings: `agg_method = 'distweight'`, `apply_bonferroni = false`,
`alpha = 0.05`, one-sided paired signed-rank tests.

These outputs carry the known defects of the retiring route
([ROADMAP.md](ROADMAP.md), section 4): the slab artefact of the `abs()` in
`run_characterize_distribution`, the stale autofluorescence of its `auto` run
and of `run_compare_nano_with_autofluorescence`, reliability bars that start
at white. They are not quoted.

## What the map measures

### The ISH analysis (October 2026)

The guided figures of the ISH line, built on 8 October 2026 and run on a full
copy of the production inputs: what the adult map is, read against the Allen
ISH maps. They are numbered in the order of the argument, in
`<data>\adult_v2\ish_analysis\figures\`, PNG and EPS at 150 dpi;
`figures\README.md`, written by `run_ish_overview.py`, walks through them with
the numbers of the run. The story and every number:
[ISH_ANALYSIS.md](ISH_ANALYSIS.md).

| figure | made by | what it shows |
|---|---|---|
| `00_overview.png` | `mapping/run_ish_overview.py` | the question, the argument in two parts with this run's numbers and where each part stands, and its limit; one row per step with the numbers of the run and what stays open; the A-items and their state |
| `01_structures.png` | `mapping/run_structure_set.py` | the declared structures (A1): the rule as a funnel, kept and left out per division, the map on plane 700, how far each adult's `zref` moves |
| `02_genes.png` | `mapping/run_ish_gene_table.py` | the two panels and their union, experiments per gene, the sections set missing in P9's experiments, reliability, what was left out and the repair |
| `03_beyond_budget.png` | `mapping/run_beyond_figures.py` | part 1: the map against Gria1, against synaptic density and against the model; the variance budget with the calibration floor and the one-Gria1-experiment benchmark; the calibration; the leftover half against half; the leftover under other folds and structures |
| `04_beyond_where.png` | `run_beyond_figures.py` | where the leftover sits on three planes, and the structures furthest from prediction |
| `05_one_comparison.png` | `mapping/run_ish_gene_ranking.py` | what one gene's rho is: nano, Cacng8, Gria1 and Aqp4 on plane 700, as measured, as ranks, and the scatter of ranks |
| `06_spatial_null.png` | `run_ish_gene_ranking.py` | why a null (unrelated smooth maps correlate), the surrogates' variogram against the map's, three surrogates, the false-positive rates, Cacng8 and Gria1 against their nulls |
| `07_gene_ranking.png` | `run_ish_gene_ranking.py` | P9's genes ranked, each bar against its null band, grey by its steadiness across adults, autofluorescence's rho beside it; the Cacng8 - Gria1 gap against maps related to both alike; the two genes with the label and with the tissue |
| `08_gene_sets.png` | `mapping/run_ish_gene_sets.py` | the gene sets fixed in advance against their null bands; the two contrasts named in advance and the genes of the presynaptic set |
| `09_localisation.png` | `run_ish_gene_sets.py` | localisation genes against expression-matched controls once the subunit composite is removed, the positive controls, what the test can find, every test of the design |
| `10_between_within.png` | `mapping/run_ish_divisions.py` | rho against a division-only map, the mean rho inside divisions, five genes division by division, which within null holds |
| `11_leftover_genes.png` | `run_beyond_figures.py` | every gene and gene set against the leftover of part 1, against surrogates put through the same fit |
| `12_autofluorescence.png` | `run_ish_gene_ranking.py` | every gene against the autofluorescence map of the same sections, adult by adult for Gria1 and Cacng8, how many genes each map passes at three thresholds, the genes past each null |
| `13_robustness.png` | `mapping/run_ish_robustness.py` | the ranking under eleven variants: gene order, Cacng8's and Gria1's ranks, the gap against its null |
| `14_green_channel.png` | `mapping/run_sep_channel_check.py` | what the green channel reports: three channels, one adult's raw planes, which channel follows which, their ranges, against Gria1 |
| `15_april_headline.png` | `run_ish_overview.py` | April's category violins beside today's, gene by gene, the ANOVA p under each choice, today's groups against the null |
| `qc\00_flagged.png`, `qc\<gene>_<experiment>.png` | `mapping/run_ish_section_qc.py --sheets` | every experiment with a section set missing, kept as a true absence or kept at a step; one sheet per experiment: the section profile with its flags, the orientation check |
| `genes\<gene>.png` | `run_ish_divisions.py --sheets` | Cacng8, Gria1, Grm5, Dlg2 and Aqp4: rank maps, the scatter with a fitted line per division |

Working figures beside the tables: `beyond\fig0_structures` to
`fig6_readings`, `E_regression`, `F_maps` (steps 21, 22 and 24) and
`green_channel\sep_channel_check.png` (step 26).

**Inputs.** The per-brain files and `region_means_per_mouse.csv` of the route's
steps 1 to 5 (the stored `cref`, `zref` and `ratio` rows), the ontology panel
and its grids of steps 11 and 12 (`adult_v2\panel\`, `<data>\atlas_ish\`),
`mapping/gene_targets.csv` and `mapping/ish_section_exceptions.csv`, the
frozen P9 table `<data>\comparisons\merged_naive_rws_vs_ish_summary_nosmooth\gene_panel_summary.csv`
(figure 15), and the caches under `ish_analysis\cache\` (Allen experiment
lists, mygene records, `go-basic.obo`). Steps 13 to 27 run in order.

**Settings that matter.** `[structures] min_adults` (10) and `grey`;
`[ish_qc]` (three neighbours each side; a section below 0.2 of the median of
its neighbours and of the brightest on each side is set missing, not judged where
the neighbours read below 0.1); `[ish_analysis] q` (0.05),
`min_division_structures` (8), `min_set_genes` (5); `[spatial_null]`
(10,000 surrogates, the variogram matched to the 25th percentile of the
distances, 2,000 calibration maps); `[beyond]` (20 shufflings of the folds), `[beyond_calibration]` (the
jackknife), `[beyond_figures]`; `[ish_figures] plane` (700) and `t_max` (40, the grey of a
gene's bar). Random seeds are fixed.

### The ISH comparison (April 2026)

**`adult_matlab/run_compare_with_allen_ish.m` (`P9_compare_nano_vs_allen_ish.m`;
21 and 22 April 2026), until the MATLAB adult scripts retire.**
`<data>\comparisons\merged_naive_rws_vs_ish_summary_nosmooth\`:
`correlation_barchart_4metrics.png`, `violin_spearman_dw.png`,
`violin_spearman_ero.png`, `violin_voxel_pearson.png`,
`correlation_eroded_vs_distweight.png`, and the ranking `gene_panel_summary.csv`;
one folder per gene, `merged_naive_rws_vs_ish_<gene>_nosmooth\`, with
`Scatter_NanoVsISH_<gene>.png`, the paired bars, the comparison video and
two diagnostic sheets. Made from the cohort cache of
`run_characterize_distribution`, `<data>\gene_targets.csv` (the 100-gene
panel) and the Allen grids in `<data>\atlas_ish\`. Settings: the nine
divisions, erosion radius 3, distance-weight power 4, planes 100 to 700.
These folders have no channel in their names: they predate it, and a rerun
writes `merged_naive_rws_nano_vs_ish_*`. They carry the script's defects
(the stretched grid, the section repair). `violin_spearman_ero.png` is the
headline Sami El-Boustani saw; figure 15 of the ISH analysis shows what is
left of it.

## Young against adult: figures and videos

All in `<data>\comparisons_v2\young_vs_adult\` unless stated, first made on
24 and 25 September 2026 and rerun on 5 October (seven young brains, ten
adults), every reading in force (`ratio`, `sepratio`, `cref`, `subref`,
`zref`). The barrel-field variant was not rerun (its script no longer runs,
[ROADMAP.md](ROADMAP.md), section 8).

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

- The cohorts: the `mapping_cohort` column of the cohort table
  `common\cohort.csv` (decision Y4), grouped by `COHORTS` in
  `sepmap/volumes/cohort.py`; `young` pools P16, P20 and P22.
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
- `subref` everywhere, maps, tables, figures and video: its reference now
  leaves out fibre tracts, ventricles and the unassigned labels, and its
  titles say so (each brain's `subref` moves by one constant; the numbers are
  in [ROADMAP.md](ROADMAP.md), section 1).
- `sepratio` young means: MG897's 183 missing voxels no longer count as zero
  in the cohort mean (the young mean there rises by a factor of 7/6, n 7 to
  6), and the region tables leave its 117 missing voxels out of its
  nucleus accumbens mean (one value of `region_means_per_mouse.csv` and three
  of `region_stats.csv` move in the fourth decimal). The P20-only maps
  `log2_alt_*` in `volumes_ccf20.npz`, drawn by no figure, are smoothed
  without counting a missing value as zero.
- `region_plot.png`: the `ratio` panel's title says what its zero means
  ("0 = equally bright").
- The videos' reliability t: each brain counted once, its hemispheres
  averaged first (the adult map, above).
- Decision 4 of 4 October: each brain's backgrounds over its imaged voxels
  only, and the young brains carried to the CCF without darkening their
  tissue edges. In the maps, videos and close-ups the young cortical surface
  reads brighter. In `group_stats.csv` the differences move by at most 0.01;
  per structure most move by a few thousandths, and a few small structures
  at the tissue surface by more (the supramammillary nucleus, measured in
  four or five young brains, by 0.25 to 0.57 log2 depending on the
  reading). The `subref` structures at Welch q < 0.05 are 83 (104 before
  fix 1, 81 with fix 1 alone). Four `zref` structures that sat on the
  Benjamini-Hochberg boundary leave q < 0.05, so that count is not quoted.
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
Changed since: sheet 01's title says what the tissue mask needs; the titles of
sheets 06, 08 and 09 state what the numbers show (fix 27). Since decision 4
the sheets draw unimaged voxels white instead of the no-data grey
([ROADMAP.md](ROADMAP.md), section 8).

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
  [ROADMAP.md](ROADMAP.md)). Changed since: the ratio video
  (`ratio_map_video.mp4`, only with `saveRatioMap`, off in production) takes
  the palette's difference map, blue below 1 and red above (fix 27).
- `slice_order_montage.png`, every slice as a tile in its curated order:
  the slice-order editor of `preprocessing/run_order_slices.m` (`P1bis`),
  saved next to the decisions file. The README's image is one of these.
- `<mouse>_dim<d>_initial_registration.png`, the initial 3D fit to the atlas:
  `registration/run_register_to_atlas.m`, mode `align`.

Across brains:

- `<data>\intensity_diagnostics\`, the slice-to-slice equalisation of the nano
  channel (`Plot_*_<time>.png`, `Video_<mouse>_<time>.mp4`, the statistics):
  `preprocessing/run_nano_equalisation.m` (`P2bis`); the names carry the run's
  clock time. Changed since: the saved inter-quartile ranges are positive
  (they were saved negated; fix 8).
- `<data>\<group>\Background_mask_diagnostics_trace_<group>_<channel>.png` and
  `<group>\global_diagnostics\normalization_checks_<channel>\`:
  `group_comparison/run_normalise_groups.m`. Changed since: the background
  trace takes its y range from every mouse, not the last one.
- `<data>\young\registration_qc\` (`ATLAS_PARAMETERS.md`,
  `atlas_crop_and_ap_mapping.png`, `atlas_comparison_adult_vs_p20.png`,
  `atlas_region_comparison.png`, `demba_to_allen_transform_qc.png`): the atlas
  checks in `atlas/qc/` (`atlas_diagnostics`, `compare_atlases_montage`,
  `compare_atlas_regions`, `check_demba_to_allen`, which reads the annotation
  `demba_to_allen.py` beside it writes), run by hand.
  `registration/qc/check_registration_error.m` prints each aligned brain's
  fit error.

## Superseded outputs, not to quote

- `<data>\comparisons\young_P20_vs_adult_nano\`: the first young-against-adult
  attempt (4 September 2026, three P20 brains), made by five Python scripts now
  in `archive/` (`compare_young_vs_adult_lrsum.py` and the others); replaced by
  the Python route.
- `<data>\comparisons\merged_young_P20_nano\`: `run_characterize_distribution`
  (P8) run on the P20 brains (`SEP_MERGE_SPECS=young_P20`); replaced by the
  Python route, and later by A4's young cohort.
- `<data>\comparisons\merged_naive_rws\`: the earlier nano run of
  `run_characterize_distribution`, before the channel entered its folder
  names; `merged_naive_rws_nano\` is the later one.
- The `Diagnostic_DistWeight_*` sheets in the folders of
  `run_characterize_distribution`: dropped by decision (29 September), not
  rebuilt.
- `<data>\adult_v2\beyond\` and its `for_sami\` panels A to F with
  `numbers_for_the_caption.txt` (`run_beyond_*`, 26 September, rerun on 5
  October): replaced by figures 03, 04 and 11 of the ISH analysis and
  `<data>\adult_v2\ish_analysis\beyond\`. What changed: the four subunits
  enter as separate predictors (36% left becomes 27%), the ceiling is
  Spearman-Brown's value and not its square (97.4% becomes 98.7%), `zref`
  takes the declared reference, the covariates come after section QC, control
  A uses one-hemisphere centroids, the held-out R2 is averaged over 20
  shufflings of the folds, the range over structures comes from a jackknife
  (the bootstrap let a structure sit in a training and a test fold), and a
  calibration with the same model gives the floor the leftover is set beside,
  compared on the same structures.
- `<data>\adult_v2\ish\`: `ish_old_vs_new.png`, `ish_word_enrichment.png`,
  `ish_roles.png`, `ish_reliability.png`, `ish_panel_test.png` and their tables
  (`run_ish_regions`, `run_ish_compare`, `run_ish_words`, `run_ish_roles`,
  `run_ish_reliability`, `run_ish_panel_test`; 25 and 26 September, rerun on
  5 October): replaced by figures 02, 05 to 13 and 15. What changed: the
  declared structure set, section QC, one gene table merging every usable
  experiment, and a spatial null for every correlation; the word and roles
  tests are retired, their questions asked by the gene sets fixed in advance
  (figure 08) and the localisation test (figure 09).
- `<data>\adult_v2\arms\arms_vs_genes.png` and `arm_gene_correlations.csv`
  (`run_ish_arms`): retired (decision 5 of the ISH discussion, [history/ISH_DISCUSSION.md](history/ISH_DISCUSSION.md)): the green
  channel is not total receptor, so the tests of the channel ratios against
  genes have no premise. `arms_consistency.png` and `region_means_arms.csv`
  (`run_adult_arms`) stay, frozen; no step reads them.
- `<data>\adult_v2\arms\sep_channel_check.png` and `.csv`
  (`run_sep_channel_check`, 26 September, rerun on 5 October): replaced by
  figure 14 and `<data>\adult_v2\ish_analysis\green_channel\`, on the declared
  structures and the merged Gria1 profile, with neither "total receptor" nor
  "surface fraction" in its labels.
