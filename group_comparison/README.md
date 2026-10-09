# Plasticity comparison

Where an experience changes the nanobody signal: naive mice against mice
after rhythmic whisker stimulation (RWS) or after behavior, voxel by voxel,
on the registered adult brains (Allen CCF). The line of work is paused, not
closed: the code stays runnable with identical results, so it can be resumed
with more animals ([`../docs/ADDING_DATA.md`](../docs/ADDING_DATA.md), step 5).

## Run order

| step | script | what it does |
|---|---|---|
| 1 | `run_collect_by_group` | stack the registered nano volumes of a group's mice into one 4D array (AP x DV x ML x mouse), in cohort-table order |
| 2 | `run_normalise_groups` | per mouse and plane, a background mask; each mouse fitted onto the group's median cortex with a robust line, then applied to the whole volume |
| 3 | `run_group_differences` | each mouse folded onto the left hemisphere (L - R and L + R; the group means take \|L - R\|); the experimental group aligned onto the control group's profile; tissue smoothed in 3D, each mouse NaN outside its tissue; group means, Welch t and surprise (-log10 p) maps over the mice with tissue at each voxel, a t only where each group has at least `min_mice_per_group` of them (3), slab figures, regional bars, five region measures with an exact label permutation test, videos |

The comparison as it was approved (5 December 2025) runs:

1. `run_collect_by_group` with `mousetypes_list = {'rws', 'naive', 'behavior'}`
   and `age_filter = []`;
2. `run_normalise_groups` with `SEP_COHORT_SPECS=rws,naive,behavior` in the
   environment;
3. `run_group_differences` twice, `exp_type = 'rws'`, then `'behavior'`
   (`ctrl_type = 'naive'`).

Its mice: naive CGF027, CGF028, CGF033, CGF034, CGF035; RWS MG691, MG692,
MG693, MG736, MG737; behavior MG705, MG709, MG716, MG718, four of the seven
registered. The selection is written at the top of steps 2 and 3: by position
within each group (`selected_mice_idx_list`), and in step 3 the behavior mice
to analyse by name (`behavior_mice`).

Steps 1 and 2 also serve the young brains: with `mousetypes_list = {'young'}`
and `age_filter = [20]`, then `SEP_COHORT_SPECS=young_P20`, they give the
P20 stack and its background masks, which sheet 08 of
`../mapping/run_diagnostics.py` compares with the Python route's tissue mask.

## Outputs

- `<data>\<group>\nano_4d<tag>.mat` and `collected_mice<tag>.mat` (step 1;
  `<tag>` is empty for a whole group, `_P20` for `age_filter = [20]`)
- `<data>\<group>\nano_4d_normalized<tag>.mat`,
  `nano_4d_normalized_bkgmask<tag>.mat`,
  `Background_mask_diagnostics_trace_<group>_nano.fig/.png` and
  `global_diagnostics\normalization_checks_nano\` (step 2)
- `<data>\comparisons\<ctrl>_vs_<exp>_<channel>\` (step 3): the profile
  alignment, the slab figures around plane 565 (group t maps masked by
  surprise, individual mice), the regional surprise bars and the videos, each
  figure as `.fig` and `.png`. No `.mat`: the values are in the `.fig` files,
  and the regional bars' in `Region_Surprise_DiffSum_<comp>.csv` (per region:
  its atlas division, its voxels in the left hemisphere, the regions of the
  list taken out of it, and for L - R and L + R its voxels with a t, those at
  p < 0.01, their fraction, which is the bar, and their summed surprise, the
  bar before 4 October 2026). Since 6 October 2026 also
  `Region_Surprise_Bar_<bar_measure>_<comp>` (the bars of the measure chosen,
  the cluster mass, shaded by their corrected permutation p and starred where
  it is below 0.05; a dagger where only the uncorrected p is below 0.05; the
  uncorrected p of the regions named in advance in bold beside their bars;
  the names of the regions expected to change in magenta; two lines under the
  title say which p tests what), `Region_Measures_<comp>` (each
  region's rank and corrected p under the five measures), the measures'
  columns of the table (per map and measure the score with its sign, p and
  corrected p; the heaviest cluster's voxels and peak; the top volume's
  voxels) and `Region_Permutation_Null_<comp>.mat` (every split's scores). With
  `n_permutations` a number, these four carry `_perm<n>` in their names.

The approved outputs are in `<data>\comparisons\naive_vs_rws\` and
`naive_vs_behavior\`, named before the channel joined the comparison tag.
Today's code writes `naive_vs_rws_nano\` and `naive_vs_behavior_nano\`
beside them, last on 5 October 2026 from the normalised stacks then in the
group folders (naive of 4 September 2026, RWS of 1 May 2026, behavior of
5 December 2025), not those of the approved run. Naive and RWS were
normalised months apart: a clean final version reruns step 2 on the three
groups together, then step 3 ([`../docs/ROADMAP.md`](../docs/ROADMAP.md),
section 2).

## Per-mouse values

`run_per_mouse_values`, after step 3, gives one value per mouse in regions
named in advance, the barrel field (SSp-bfd) and, as a cortical control, the
primary visual area (VISp), to show which mice carry an effect. Each mouse is
taken as step 3 takes it, and in each region, over its own voxels:

| value | what it is |
|---|---|
| `ai` | asymmetry index, mean \|L - R\| over mean (L + R), on the maps of step 3 (folded, smoothed, the experimental group aligned) |
| `ai_raw` | the same on the collected stack (`nano_4d.mat`) less the mouse's off-tissue level, smoothed the same way: no normalisation, no alignment |
| `signed`, `signed_raw` | mean (L - R) over mean (L + R): above 0, the left hemisphere higher |
| `sum_rel`, `sum_rel_raw` | mean L + R in the region over the mouse's own mean over the isocortex |
| `loo_ai`, `loo_ai_raw` | the AI in the heaviest positive L - R cluster of SSp-bfd that the comparison finds without the mouse |

- The off-tissue level is the median of the raw stack over the mouse's
  background voxels (step 2's mask) outside the atlas brain: the slide around
  the section, about 500 raw where tissue is 1,200 to 2,700. Each mouse's
  normalised volume is checked to be its raw volume through its line of step
  2, to the bit.
- The alignment's intercept enters the experimental mice's L + R on the maps
  of step 3, and so their `ai`, `signed` and `sum_rel`; the raw values do not
  have it.
- Leave-one-out: without each mouse in turn, the comparison is redone as step
  3 does it (the alignment refitted on the others, the t with three mice per
  group, the surprise, the rolling median, p < 0.01, 18-connected, within
  SSp-bfd), by `region_permutation_test` on the observed split alone, and the
  mouse is read in the cluster the others give: no mouse is read in a cluster
  its own data helped define. The one thing not redone is step 2, which fits
  each mouse onto its group's median cortex: a scale per mouse, which cannot
  move a cluster. Fold 0 keeps every mouse and must give step 3's cluster; it
  is checked against step 3's table.
- Each value is compared between the groups by the exact permutation of the
  difference of the group means (252 splits for 5 and 5, 126 for 5 and 4),
  two-sided and one-sided for the experimental group higher (not for the
  signed values, whose side the experiment does not name), with the Welch t
  and Hedges' g beside it.
- Outputs, in the comparison's folder, `<tag>` the comparison and the
  smoothing (`naive_vs_rws_nano_smooth5`): `Per_Mouse_Values_<tag>` (.csv,
  one row per mouse, and the figure), `Per_Mouse_Stats_<tag>.csv`,
  `Per_Mouse_LOO_<tag>.csv` (one row per fold) and the cache of the mice's
  maps, `Per_Mouse_Maps_<tag>.mat` (`force_recompute_mice` redoes it).

## Where the code is

Each driver sets its settings and calls one function in `pipeline/`
(`collect_by_group`, `normalise_groups`, `group_differences`,
`per_mouse_region_values`); the video
writers (`write_lr_*`, `open_lr_video`, `set_lr_colormap`) and the atlas
outlines (`lr_atlas_boundaries`) sit beside them, and so do the parts of
step 3 that other code reuses: each mouse's tissue and its smoothing
(`tissue_only`), its plane profile (`plane_tissue_means`), the alignment of
the two groups (`align_exp_to_ctrl`), the regions of the bars
(`surprise_regions`, `surprise_region_acronyms`) and their permutation test
(`region_permutation_test`). Shared: `../common/`
(`get_cohort`, `get_cohort_spec`, `compute_lr_stats`, `sep_palette`) and
`../atlas/` (`get_atlas_crop`, `get_allen_region_mask`).

## Notes

- The region measures (since 6 October 2026, `region_permutation_test`).
  The share of significant voxels, the bar since 4 October, dilutes a focal
  bump inside a large region (the RWS bump inside SSp-bfd), and the summed
  surprise grows with the region's size. Each region is therefore scored five
  ways on the same rolling median, here of the surprise signed by the group
  difference, positive and negative effects apart, the larger with its sign:
  the share; the summed surprise at p < 0.01; the 99th percentile; the mean of
  its most surprising voxels over 0.1 mm^3 (`topvol`, 100,000 voxels of 10 um,
  all of them in a smaller region); and the mass of its heaviest cluster at
  p < 0.01 (18-connected, within the region's own voxels). Each score is
  tested by relabelling the mice: every split of the pooled mice into groups
  of the original sizes (252 for 5 against 5, 126 for 5 against 4), each with
  the t, the three-mice rule, the surprise and the rolling median computed
  again; a region's p is the share of splits reaching its |score|, the
  observed one included (so at least 1/126, and 2/252 for equal groups, whose
  splits come in mirror pairs), and its corrected p the share whose largest
  |score| over the regions reaches it. The splits relabel each mouse's folded
  map after the alignment of the two groups: the groups were normalised on
  scales of their own (step 2), and refitting the alignment line between mixed
  groups would put two scales into each group. The sum of the measures is
  taken on the signed median, one sign at a time, so it is not the table's
  `surprise_sum` (the unsigned median, both signs together), which stays
  beside it for reference. A corrected p compares a region's score with the
  largest score of any region under relabelling, and the measures differ in
  which regions give that largest score: for the sum, the top volume and the
  cluster, the ten largest regions (HPF, CP, MBmot, HY, MOs, PIR, MOp, SSs,
  VISp, SSp-bfd) give 54 to 98% of the splits' maxima in both comparisons;
  for the share and the 99th percentile, 1 to 16%, the rest coming from the
  small regions, whose share and percentile vary most.
  Giulio chose the cluster mass for the bars on 7 October 2026, from the
  comparison figure of both comparisons. The share's p, in the table, is the
  chance level of the share bars asked for in April (Sami El-Boustani).
- What the two p of the bars test. The corrected p is the test of a search
  over all 71 regions. A region named in advance (`a_priori_regions`: since
  7 October 2026 the barrel field, since RWS stimulates the whiskers and the
  behavior task uses one whisker) needs no search: it is one test, and its
  uncorrected p is that test, written in bold beside its bar. For any other
  region the uncorrected p is not a test, and its dagger is no result: in
  behavior L + R, 50 of the 59 regions with a score are below 0.05
  uncorrected. The table's `a_priori` column marks the regions named in
  advance. The barrel field's cluster mass, higher in the experimental group
  in all four (production run of 7 October 2026, `comparisons\naive_vs_<exp>_nano\`;
  in brackets, the run on the December 2025 inputs):

  | comparison | map | rank | uncorrected p | corrected p |
  |---|---|---|---|---|
  | naive vs RWS | L - R | 1st of 71 | 0.040 (0.032) | 0.635 (0.635) |
  | naive vs RWS | L + R | 5th | 0.056 (0.048) | 0.667 (0.683) |
  | naive vs behavior | L - R | 6th | 0.15 (0.119) | 0.944 (0.929) |
  | naive vs behavior | L + R | 4th | 0.0079 (0.0079) | 0.016 (0.016) |

  The result quoted for RWS (Giulio, 8 October 2026) is the production run's
  targeted test of the barrel field, the prediction the experiment was
  designed on: the heaviest cluster of the L - R map, first of all regions,
  p = 0.040 (L + R 0.056). The search over all regions does not show it.
  After behavior, the barrel field's L + R holds under both p, its L - R
  under neither. The two runs differ slightly because production's naive
  stack was normalised again in September 2026.

- The headline of the approved comparison is a small increase of the
  nanobody signal in S1 after RWS, in the hemisphere-sum t map of the slab
  around plane 565, masked at p < 0.01, uncorrected, with five mice per
  group. On 1 October 2026 today's step 3, run on the normalised volumes of
  that time (`nano_4d_normalized_bk.mat` in `naive\` and `rws\`), reproduced
  it: the slab t maps and surprise masks correlate with the approved ones at
  0.997 to 0.998 for RWS and 0.988 to 0.999 for behavior
  ([`../docs/history/REFACTOR_PLAN.md`](../docs/history/REFACTOR_PLAN.md),
  Progress). This checks the computation on the same mice; it is not a
  replication in new animals.
- Step 3 used to smooth each mouse's tissue and then set the voxels outside
  it to 0, so those zeros entered the group means and the SEM's n. Since
  fix 23 of step 8 they are left out. On the same inputs the RWS S1 increase
  holds (the difference map unchanged in S1, the sum map's significant area
  there a third smaller), while naive against behavior changes over the
  whole slab: MG709 has no tissue in the slab around plane 565 and no longer
  counts there as a fourth mouse (`../docs/ROADMAP.md`, section 1). Voxels
  left with two mice in a group, at tissue edges, gave |t| up to 73; since
  Giulio's decision of 4 October 2026 a voxel gets a t only with at least
  three mice in each group (`min_mice_per_group`).
- Step 2 selects the adults by position: step 1 stacks a group in table
  order, and step 2 picks mice by their place in that stack; step 3 names
  them by the same positions, except the behavior mice it analyses, which it
  takes by name from those step 2 saved (a name not saved stops the run). The
  adults' rows of the cohort table keep their order (`get_cohort('verify')`).
- This route normalises each mouse onto its group and the experimental group
  onto the control group; the Python route (`../mapping/`) scales each brain
  on its own. A result that uses both must make them comparable first.
- Only the nano stack is written. The `auto_4d.mat`, `mask_4d.mat` and average
  files in the adults' group folders are older (19 November 2025) and come
  from an earlier registration: not to be read (`collect_by_group`'s help).
  The autofluorescence per brain is in the registered tiffs, which the
  Python route reads.
- `../adult_matlab/run_characterize_distribution`, kept until the Python
  route replaces it, reads the normalised stacks of step 2; the other two
  scripts of `../adult_matlab/` read its outputs.
