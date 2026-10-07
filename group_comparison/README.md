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
  shaded by their corrected permutation p, a dagger where only the
  uncorrected p is below 0.05, the uncorrected p of the regions named in
  advance beside their bars), `Region_Measures_<comp>` (each
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

## Where the code is

Each driver sets its settings and calls one function in `pipeline/`
(`collect_by_group`, `normalise_groups`, `group_differences`); the video
writers (`write_lr_*`, `open_lr_video`, `set_lr_colormap`) and the atlas
outlines (`lr_atlas_boundaries`) sit beside them. Shared: `../common/`
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
  The corrected p is the test of a search over all 71 regions. A region
  named in advance (`a_priori_regions`: since 7 October 2026 the barrel
  field, since RWS stimulates the whiskers and the behavior task uses one
  whisker) needs no search: it is one test, and its uncorrected p is that
  test. For any other region the uncorrected p is not a test, and its
  dagger is no result. In the run on the December 2025 inputs the barrel
  field is first by cluster mass in RWS L - R, at p = 0.032 uncorrected
  and 0.635 corrected. The table's `a_priori` column marks these regions.

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
