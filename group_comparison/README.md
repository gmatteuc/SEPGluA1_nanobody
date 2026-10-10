# Plasticity comparison

Where an experience changes the nanobody signal: naive mice against mice
after rhythmic whisker stimulation (RWS) or after behavior, voxel by voxel,
on the registered adult brains (Allen CCF). The line of work is paused, not
closed: the code stays runnable with identical results, so it can be resumed
with more animals ([`../docs/ADDING_DATA.md`](../docs/ADDING_DATA.md), step 5).

## Run order

| step | script | what it does |
|---|---|---|
| 1 | `run_collect_by_group` | stack the registered nano volumes of a group's mice into one 4D array (AP x DV x ML x mouse), in cohort-table order; the autofluorescence too with `channels = {'nano', 'auto'}` |
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
  `<tag>` is empty for a whole group, `_P20` for `age_filter = [20]`), and
  `auto_4d<tag>.mat` with `channels = {'auto'}`; each stack with its mice and
  the registered files they were read from
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

In short (9 October 2026). Giulio asked for a per-mouse measure and test that
is step 3's barrel-field test seen from the mice, as permissive as it. The
selection-matched test below is that as far as the choice goes: the same
cluster search, with the same ten mice, under the same 252 splits. It does
not agree with step 3: after RWS its p is 0.35, one-sided (step 3's cluster
mass over the same splits: 0.040 with either sign, the p quoted; 0.020 with
the positive sign alone). It is not harsher on the mice. Step 3 scores each
split's cluster by its mass, which grows with the cluster's extent, and the
extent is what is rare under relabelling (5 of the 252 splits as heavy in
the positive direction). A mean over the cluster cannot carry the extent, and
every split's cluster separates its own groups by construction: the observed
difference, +0.084, is about the null's median, +0.080. No per-mouse view of
the asymmetry separates the groups:

| view | chosen, with which mice | null | how strict | p after RWS |
|---|---|---|---|---|
| selection-matched (`sm_ai_raw`) | step 3's cluster, all ten mice; every mouse read in it | the search redone in each of the 252 splits | as step 3 about the selection; every mouse in sample | 0.35 one-sided, 0.60 either direction |
| leave-one-out (`loo_ai_raw`) | a cluster of the other nine mice; each mouse read out of sample | the leave-one-out redone in each split | stricter about circularity by design | 0.13 one-sided |
| whole region (`ai_raw_ssp_bfd`) | nothing chosen | the values permuted | no selection | 0.25 one-sided |
| step 3's cluster mass, for reference | step 3's cluster, all ten mice | the same 252 splits | the test quoted | 0.020 positive sign, 0.040 either sign |

After behavior none separates the groups either (selection-matched 0.53
either direction, 0.36 behavior higher; step 3's mass 0.15 either sign). In
the RWS cluster the autofluorescence is higher after RWS too, by about half
the nano channel's difference (results below).

`run_per_mouse_values`, after step 3, gives one value per mouse in regions
named in advance, the barrel field (SSp-bfd) and, as a cortical control, the
primary visual area (VISp), to show which mice carry an effect. Each mouse is
taken as step 3 takes it. Its first values are read in the barrel-field
cluster of step 3's test, under the selection-matched test below, which has
the same selection and splits as that test, with the mice's mean asymmetry
in place of the cluster's mass; the others over a whole region, over the
mouse's own voxels:

| value | what it is |
|---|---|
| `sm_ai_raw`, `sm_ai`, `sm_ai_auto` | the AI (below) in the heaviest cluster of SSp-bfd where \|L - R\| is higher in the experimental group, as step 3 finds it with all the mice: on the raw stack, on the maps of step 3, and on the autofluorescence |
| `ai` | asymmetry index, mean \|L - R\| over mean (L + R), on the maps of step 3 (folded, smoothed, the experimental group aligned) |
| `ai_raw` | the same on the collected stack (`nano_4d.mat`) less the mouse's off-tissue level, smoothed the same way: no normalisation, no alignment |
| `sum_rel`, `sum_rel_raw` | mean L + R in the region over the mouse's own mean over the isocortex |
| `loo_ai`, `loo_ai_raw` | the AI in the heaviest cluster of SSp-bfd where \|L - R\| is higher in the experimental group (a positive t), as the comparison finds it without the mouse |

There is no signed value. Left and right are not certain for every brain, so
the stimulated side is not known mouse by mouse, and a signed L - R would
carry arbitrary signs; the test takes \|L - R\| and L + R for the same
reason. The whisker stimulated was C2 (possibly B2), the same in every RWS
mouse (Giulio, 9 October 2026).

### The selection-matched test

The main per-mouse test, fixed on 9 October 2026 before any of its numbers
(Giulio: the voxelwise test seen from the mice, as permissive as it; the
leave-one-out below is stricter). Its p is reported whatever it is; no other
variant is tried.

- The cluster. With the groups as they are, the heaviest cluster of SSp-bfd
  where \|L - R\| is higher in the experimental group (a positive t), found
  with all the mice as step 3 and `region_permutation_test` find it: the t
  where each group has three mice with a value, its surprise, the median over
  +/- 10 planes, p < 0.01, 18-connected, within the region's voxels of the
  left hemisphere. It is the cluster whose mass is the barrel field's score
  in the region test (p 0.040 after RWS, production run of 7 October). Each
  mouse's asymmetry index is read in it, mean \|L - R\| over mean (L + R)
  over the cluster's voxels where the mouse has a value: on the raw stack
  less the mouse's off-tissue level (`sm_ai_raw`, the value to read: no
  normalisation, no alignment), on the maps of step 3 (`sm_ai`) and on the
  autofluorescence (`sm_ai_auto`). The statistic: the experimental group's
  mean minus the control group's.
- Its null. Under each of the 252 splits of the ten mice into five and five
  (126 for five and four after behavior), the cluster search is redone with
  the split's groups (its experimental group in the experimental group's
  place, a positive t), every mouse is read in that split's cluster, and the
  same statistic is taken. The mice keep the alignment of their true groups,
  as the region test keeps it under its splits. The one-sided p is the share
  of the splits whose statistic reaches the observed one, the observed split
  included; after RWS it is the test, its direction named in advance (RWS
  potentiates the stimulated barrels' synapses, Gambino et al. 2014). The
  two-sided p is the search in either direction, as the region test's score
  takes the larger of its two signs: per split, the larger of the statistic
  in its positive cluster and, in its negative cluster (where its control
  group's \|L - R\| is higher), the control group's mean minus the
  experimental group's. After behavior the direction is carried over from
  RWS and the two-sided p comes first. A split whose search finds no
  cluster, or that leaves a group without a mouse with a value, gives a
  statistic of 0, which does not reach a positive observed one; how many
  there are is given, and the p with those splits left out beside it.
- Why it is valid. Each value is read in a cluster chosen with the same mice
  and the same labels, so the observed statistic is inflated by the
  selection. The selection is redone the same way in every split, so the
  null is inflated alike, and the p compares like with like. Shuffling the
  observed values instead leaves the selection out of the null: that p is
  anti-conservative and is kept in the statistics table only, so labelled.
- Why it matches the voxelwise test. The same cluster search, the same
  splits, the same relabelling of the mice's maps as the barrel field's
  cluster mass in step 3's region test; only the number taken from each
  split's cluster differs (the mice's asymmetry, not the cluster's summed
  surprise). The run checks it: its first split must give step 3's cluster,
  and the cluster mass over its splits must give step 3's p.
- The leave-one-out is stricter: each mouse is read in a cluster found
  without it, from four mice against five, so no value carries its own
  mouse's selection. It is kept, with its full test, in a figure of its own.
- The autofluorescence (review finding F12). A misregistration of the
  surface between the hemispheres would show as \|L - R\| in every channel.
  Each mouse's autofluorescence AI is read in the same clusters, the
  observed one and every split's, and tested the same way. The channel is
  stacked from the registered volumes as the nano channel is (step 1 with
  `channels = {'auto'}`; production's `auto_4d.mat` of 19 November 2025
  predates the registration and mixes two of them, and is refused), checked
  voxel for voxel against the nano stack (the voxels a section reached must
  agree), less its off-tissue level (the same voxels as the nano channel's)
  and smoothed in the nano channel's tissue, as the raw stack is.

Added after the results (9 October 2026, evening), the design above left as
it was written. In the event the leave-one-out gave the lower p (results
below): stricter about circularity by design, not in its p. The region
test's p 0.040 is that of either sign, the larger of the two signs' masses;
with the positive sign alone, the like of this test's one-sided p, it is
0.020 over the same splits (not a test of step 3, whose p stays 0.040). The
test maps' AI (`sm_ai`) carries the alignment's intercept, which stays with
the true experimental mice under every split, so it is less exchangeable
under relabelling than the raw AI, which is the value to read.

- The off-tissue level is the median of the raw stack over the mouse's
  background voxels (step 2's mask) outside the atlas brain: the slide around
  the section, about 500 raw where tissue is 1,200 to 2,700. Each mouse's
  normalised volume is checked to be its raw volume through its line of step
  2, to the bit.
- The alignment's intercept enters the experimental mice's L + R on the maps
  of step 3, and so their `ai` and `sum_rel`; the raw values do not have it.
- An index over a mean L + R at or below 0 is NaN. On the maps of step 3 the
  zero is the normalisation's, not the absence of signal: a mouse's
  off-tissue level lands between about -2,400 (CGF033) and +500 (CGF027) on
  the control group's scale. So `ai_raw` is the index to read as a ratio,
  and `ai` the index as the test sees it.
- Leave-one-out: without each mouse in turn, the comparison is redone as step
  3 does it (the alignment refitted on the others, the t with three mice per
  group, the surprise, the rolling median, p < 0.01, 18-connected, within
  SSp-bfd), by `region_permutation_test` on the observed split alone, and the
  mouse is read in the cluster the others give: no mouse is read in a cluster
  its own data helped define. Fold 0 keeps every mouse and must give step 3's
  cluster; it is checked against step 3's table. Step 2 is not redone: it
  fits every mouse of a group onto the group's median cortex, so without one
  mouse the others' lines would change too. To first order they change by
  one line common to the group, which leaves the fold's t of \|L - R\| as it
  is (the fold's alignment is refitted, and a scale common to both groups
  does not change a t); what is left is second order.
- Each region's value is compared between the groups by the exact
  permutation of the difference of the group means (252 splits for 5 and 5,
  126 for 5 and 4), one-sided (the experimental group higher) and two-sided,
  with the Welch t and Hedges' g beside it; a cluster's value by its own test,
  which redoes the cluster under every split (the selection-matched test
  above, the leave-one-out's below). For RWS the one-sided test is the one
  named before any number: RWS potentiates the stimulated barrels' synapses
  and brings AMPA receptors to their surface (Gambino et al. 2014). For
  behavior that direction is carried over from RWS, not named on its own, so
  its two-sided p comes first (`direction_named` in the driver).
- The leave-one-out's two values get a second, full test (`loo_relabel`):
  the leave-one-out is redone under each of the 252 splits, every fold's
  cluster found again with the split's groups (the fold keeping the alignment
  of its true groups, as the region test keeps the alignment under its
  splits), and the p of the difference of the group means is taken over
  those values. A split under which no fold of a group finds a cluster gives
  no difference and is left out (9 of 252 after RWS, 15 of 126 after
  behavior). This p is the leave-one-out's test, in the table and in the
  leave-one-out's figure; its two-sided p is the share of the splits whose
  \|difference\| reaches the observed one, every fold's cluster on the split's
  experimental-higher side (not the selection-matched test's either
  direction, which searches both sides). The plain permutation of the
  values, kept in `Per_Mouse_Stats_<tag>.csv` only and named
  anti-conservative, holds each fold's cluster as the
  true groups gave it: free of circularity, since no value comes from a
  cluster its own mouse helped define (leave-one-subject-out; Esterman et al.
  2010), but it treats the ten values as exchangeable, which they are not:
  each fold's cluster is chosen with the other mice's labels, so mice of one
  group that share a pattern of asymmetry lift one another's values whatever
  the pattern's cause. After RWS it gives 0.020 where the full test gives
  0.13.
- Outputs, in the comparison's folder, `<tag>` the comparison and the
  smoothing (`naive_vs_rws_nano_smooth5`): `Per_Mouse_Values_<tag>` (.csv,
  one row per mouse, and the figure: the selection-matched values, the
  test's splits under both its statistic and step 3's, then the regions),
  `Per_Mouse_Stats_<tag>.csv` (per value the p of its test, the column
  `test` saying which, and for the selection-matched values the null's
  median difference and the cluster the either-direction difference comes
  from; for the cluster values the shuffled p, the Welch p and Hedges' g,
  which ignore the selection, in columns named anti-conservative and
  inflated, the plain columns empty), `Per_Mouse_LOO_<tag>`
  (the leave-one-out's figure; .csv, one row per fold, with where its cluster
  sits; .mat, the folds' clusters, voxel by voxel), and the caches of the
  mice's maps, `Per_Mouse_Maps_<tag>.mat`, of their autofluorescence,
  `Per_Mouse_Auto_<tag>.mat`, of the selection-matched test,
  `Per_Mouse_Selection_<tag>.mat` (every split's clusters, their mass and
  every mouse's values in them), and of the leave-one-out redone under every
  split, `Per_Mouse_LOO_Relabelled_<tag>.mat` (`force_recompute_mice` redoes
  them all). The maps' cache is read only if each group's normalised stack
  still holds the mice and the lines of step 2 it was made from, the others
  only if made from those maps and settings. For ten mice about 25 minutes
  for the maps, 25 for the autofluorescence, 6 for the selection-matched test
  and an hour for the leave-one-out under every split (25 minutes for 5
  against 4), a few minutes from the caches.

Results, 9 October 2026, on copies of the stacks of the production run of
7 October. The cluster of all the mice is step 3's barrel-field cluster (its
voxel count and mass, which is what step 3's table holds; step 3 keeps no
voxel lists), and every split's cluster mass, both signs, is that of step
3's saved null (`Region_Permutation_Null_<comp>.mat`, production's, read
only), to the bit, in both comparisons: the selection-matched test's splits
are step 3's, and the cluster mass over them gives step 3's p (0.040 after
RWS, 0.15 after behavior). After RWS (5 and 5) every split's positive
cluster is its mirror image's negative one; 5 against 4 has no mirror
splits. The leave-one-out's and the regions' values are those of the
earlier run of the day, to the bit.

The selection-matched test: the difference of the group means, experimental
minus control, with the null's median over the splits that find a cluster
in brackets; then the p named first and the other, and in brackets the same
two with the splits that find no cluster left out. The anti-conservative
shuffle of the values, the Welch p and Hedges' g are in
`Per_Mouse_Stats_<tag>.csv` only. Last row, for reference: step 3's
statistic, the cluster's mass, over the same splits.

| value | after RWS (one-sided, either direction) | after behavior (either direction, one-sided) |
|---|---|---|
| `sm_ai_raw` | +0.084 (+0.080): 0.35, 0.60 (0.46, 0.61) | +0.072 (+0.080): 0.53, 0.36 (0.56, 0.56); either direction +0.090, in the naive-higher cluster |
| `sm_ai` | +0.131 (+0.123): 0.37, 0.63 (0.48, 0.65) | +0.054 (+0.105): 0.67, 0.52 (0.71, 0.81); either direction +0.106, in the naive-higher cluster |
| `sm_ai_auto` | +0.048 (+0.036): 0.29, 0.50 (0.38, 0.52) | +0.047 (+0.025): 0.40, 0.18 (0.43, 0.29) |
| splits without a cluster | 60 of 252 (8 in neither direction) | 46 of 126 (6 in neither direction) |
| step 3's cluster mass | 0.020, 0.040 (step 3's p) | 0.15 (step 3's p), 0.056 |

- After RWS the cluster of all the mice is 3,662 voxels (0.0037 mm^3, mass
  8,775), planes 554 to 570 (CCF 733 to 749, the posterior barrel field),
  3.2 mm from the midline, 95% in layer 2/3 and 5% in layer 1. On the raw
  stack every RWS mouse but MG693 is above every naive mouse in it (MG736
  0.173, MG691 0.154, MG737 0.145, MG692 0.122, MG693 0.077; naive 0.028 to
  0.080, CGF035 the highest). The difference, +0.084, is reached by 89 of
  the 252 splits: one-sided p 0.35. It is about the null's median, +0.080,
  over the 192 splits that find a cluster. The per-mouse view does not hold
  what the region test holds (p 0.020 with the positive sign alone, 0.040
  with either sign, step 3's p), and the null shows why: what is rare under
  relabelling is the cluster's extent, not the asymmetry of its mice. Only
  6 of the 252 splits, the observed one among them, give a cluster as large
  in the positive direction (5 as heavy; with either sign 10 as heavy,
  step 3's p 0.040), but the other splits' clusters, smaller (median 188
  voxels, quartiles 36 and 658), separate their labelled groups as much:
  the 89 splits that reach the observed difference have clusters of 111
  voxels at the median, and of the 30 splits with a cluster of 1,000 voxels
  or more, 10 reach it; of the 5 as heavy, 2, the observed one among them.
  Over the splits the smaller
  clusters separate their groups the more (Spearman -0.25 between a split's
  cluster size and its difference). A mean over the cluster does not grow
  with its size; the cluster's mass does. The figure's top right panel
  shows both statistics over the same splits. The cluster where naive is the
  more asymmetric (2,194 voxels, layer 2/3, planes 451 to 470) gives naive
  minus RWS +0.082, as large.
- The autofluorescence in the same cluster is also higher after RWS, by
  about half the nano channel's difference (+0.048; RWS 0.122, naive 0.073,
  with CGF033 at 0.145, below MG692 alone), and is not separated by the
  same test (p 0.29). Mouse by mouse the two channels' AI correlate weakly
  there (r 0.32), and within the groups not at all (r -0.33, each group's
  mean taken out): the pooled r is the groups' difference. Since the nano
  channel itself shows no difference beyond the selection, the control
  neither clears the cluster of a misregistration of the surface nor points
  to one. Over the null's splits, though, a split's nano and
  autofluorescence differences correlate (r 0.66; 0.52 after behavior): the
  clusters of \|L - R\| in the nano channel pick up in part an asymmetry
  that the autofluorescence shares. The autofluorescence stacks reach the
  nano stack's voxels (Dice of the voxels a section reached 0.99998 to
  0.999997 in all 14 mice), and the nano tiffs beside them are the nano
  stack, voxel for voxel, in all 17 adults.
- After behavior the cluster of all the mice is 1,489 voxels, planes 471 to
  486 (CCF 650 to 665), 3.0 mm from the midline, all in layer 6a. Every
  behavior mouse is above every naive mouse in it on the raw stack (0.080 to
  0.146 against 0.010 to 0.061), +0.072, yet the splits' own clusters
  separate their groups as much: behavior higher, the direction carried
  over from RWS, p 0.36, the null's median +0.080. The either-direction p,
  0.53, is that of the reverse cluster, 64 voxels in layer 1, where naive
  minus behavior is +0.090, the larger of the two differences here; step 3's
  score takes the behavior-higher cluster (mass 3,454 against 143). The
  autofluorescence is higher after behavior too (+0.047, p 0.40 either
  direction), and there it follows the nano channel mouse by mouse (r 0.75,
  and 0.74 within the groups; MG705 and MG716 high in both): in this
  cluster part of the nano asymmetry goes with the autofluorescence's.
- The leave-one-out (below) gives the lower p after RWS, 0.13 against 0.35,
  though each mouse is read out of sample. A reading of why: a cluster chosen
  on noise separates the mice that chose it as well as a real one does, but
  not a mouse left out of the choice. Read in sample, the noise-chosen
  clusters that fill the null keep their full differences; read out of
  sample, they lose more of them than a real cluster would.

The whole regions and the leave-one-out. Each cell: the p named first (after
RWS one-sided, RWS higher; after behavior two-sided), then, for a whole
region, Hedges' g; for the leave-one-out values the p of the leave-one-out
redone under every split alone (their g, which ignores the selection as the
shuffle of the values does, is in `Per_Mouse_Stats_<tag>.csv` only, named
inflated). After behavior VISp has five and three mice; no p is corrected
across the values:

| value | RWS, SSp-bfd | RWS, VISp | behavior, SSp-bfd | behavior, VISp |
|---|---|---|---|---|
| `ai_raw` | 0.25, 0.47 | 0.37, 0.24 | 0.087, 1.10 | 0.16, 0.98 |
| `ai` | 0.14, 0.69 | 0.31, 0.30 | 0.98, 0.02 | 0.57, 0.49 |
| `sum_rel_raw` | 0.048, 1.09 | 0.044, 1.03 | 0.040, 1.65 | 0.036, 1.54 |
| `sum_rel` | 0.12, 0.74 | 0.044, 1.11 | 0.024, 2.02 | 0.089, 1.38 |
| `loo_ai_raw` | 0.13 | | 0.56 | |
| `loo_ai` | 0.19 | | 0.68 | |

- After RWS the asymmetry of the whole barrel field does not separate the
  groups. In the leave-one-out cluster four of the five RWS mice are above
  every naive mouse on the raw stack (MG736 0.175, MG737 0.148, MG691
  0.135, MG692 0.099; the highest naive mouse, CGF035, 0.081), and MG693
  (0.043) is among them. With every fold redone under each split, 32 of the
  243 splits reach the observed difference on the raw stack (p 0.13; 0.19
  on the maps of step 3). With the groups swapped the same procedure finds
  the naive-higher cluster above and a larger difference (raw, +0.074
  against +0.066; on the maps of step 3 every naive mouse is above every RWS
  mouse). The leave-one-out finds group-coherent spots of asymmetry in
  either direction, and the RWS-higher one does not stand out among them;
  the two groups are also two cohorts, normalised months apart.
- The RWS-higher cluster sits at the same planes in every fold, 553 to 570,
  in the region's dorsal part, midway across it, 0.001 to 0.006 mm^3. Its
  share of the cluster of all the mice runs from 0.44 (without MG693, a
  larger cluster centred 8 columns more medial) to 0.99. 87 to 100% of it is
  in layer 2/3, the rest mostly in layer 1 (3% in layer 4 without CGF034):
  just under the pia, where a left-right misregistration of the surface
  would show as |L - R|.
- L + R relative to the isocortex is higher after RWS in VISp as much as in
  the barrel field: it does not single out the barrel field. VISp is
  covered at less than half in CGF034, CGF035 and MG691.
- After behavior (alignment slope 3.12) L + R relative to the isocortex is
  higher in both regions, the asymmetry nowhere, and the leave-one-out
  clusters move from fold to fold (layers 5 and 6a, one of 85 voxels in
  layer 1; none of the voxels of the cluster of all mice in six of the nine
  folds). VISp is missing in MG709 and covered at less than half in the
  other three behavior mice.

## Each mouse's part in the barrel-field cluster

In short (10 October 2026). Giulio asked whether the RWS result is driven by
a bump in one mouse, which one, and why one mouse's bump would be enough for
the test while the mice's means do not separate. `run_mouse_influence`, after
`run_per_mouse_values`, answers from its caches, on copies of the stacks of the
production run of 7 October. Nothing in it was chosen after seeing the values;
the raw stack is the reading to look at.

- It is not one mouse. The cluster of all ten mice (3,662 voxels) is where
  several mice line up. Without MG691, MG692, MG736 or MG737 it keeps 28 to
  38% of its voxels, without CGF033, the naive mouse lowest there, 24%; its p
  over every split of the nine mice left (126) is then 0.056, 0.048, 0.10,
  0.087 and 0.20. Without MG693, the RWS mouse lowest there, it grows to 5,508
  voxels (p 0.008). Without CGF027, CGF028, CGF034 or CGF035 it keeps 67 to
  90% (p 0.024 to 0.063).
- Each of those mice is asymmetric there in its own way
  (`Influence_Maps_<tag>_raw`, `Influence_Profiles_<tag>`). MG736 in one
  section: the cluster's 17 planes are the planes between two of its section
  edges, its asymmetry runs across the whole window along ML, and its
  autofluorescence does not have it (0.082 in the cluster, nano 0.173). MG737
  in a superficial patch over that section and the next behind it. MG691 over
  a third of the barrel field (above every naive mouse at 34% of the region's
  voxels, 17% by chance), most of all 200 to 450 um in front of the cluster,
  and its autofluorescence is as asymmetric (0.139 in the cluster, nano
  0.154). MG692 moderately along 500 um of AP, less than its autofluorescence
  (0.122, autofluorescence 0.176). MG693's own bump is 400 um behind the
  cluster. All five naive mice are low there: the highest of them reaches
  0.089 in the cluster's planes, and 0.10 to 0.20 within 600 um of it. The
  cluster is no mouse's peak: each mouse's highest asymmetry in the barrel
  field is 600 to 1,290 um away, in nine of the ten at an edge of the region
  (layer 6a or 6b by the white matter, layer 1 by the pia).
- Why a few mice are enough for step 3's test while their means are not: the
  cluster's mass sums the surprise of the contiguous voxels where the Welch t
  passes p < 0.01, so it grows with the extent of the place where the RWS
  mice are high and the naive mice low and alike. The same search in any
  split of the mice finds the place where its labelled groups are most apart,
  and they are apart there about as much (+0.080 at the null's median, +0.084
  observed): what is rare under relabelling is the extent (5 of 252 splits
  as heavy), and the extent is that of the overlap of four RWS mice over one
  naive group that happens to be low. Take out any one of the mice that make
  the overlap and two thirds or more of it go.
- Whether the RWS mice agree on the place more than relabelled groups do: on
  the raw stack, no. The largest patch where four of the five RWS mice are
  above every naive mouse is 13,605 voxels (2,690 of them in the cluster), as
  large as under relabelling (median 13,546, p 0.50); where all five are, 1,120
  voxels away from the cluster (p 0.74). On the test's maps, where all five RWS
  mice are above every naive mouse at 83 to 99% of the cluster's voxels, the
  five-of-five patch is 5,400 voxels, 2,865 in the cluster, against a median
  of 1,835 (p 0.079).
- Its size: 0.0037 mm^3, 170 um along AP, 220 um along ML, 102 to 311 um
  below the pia (95% layer 2/3, 5% layer 1), about two thirds of a barrel
  column's width (300 um, the mouse's C2 column, Lefort et al. 2009) and a
  fifth of the column's layer 2/3. Along AP it is about one section: the
  sections of these brains are mostly 14 to 18 planes apart, the cluster's
  17 planes are one section of MG736 and of MG737, and the rolling median of
  step 3 (11 of 21 planes) is passed by one section of 15 planes.

### How it is done

- Folds. Fold 0 keeps every mouse; then each mouse is left out in turn and the
  comparison is redone as step 3 does it (the alignment refitted on the mice
  kept; a t where each group has three mice with a value; its surprise; the
  median over +/- 10 planes; p < 0.01; 18-connected within SSp-bfd), by
  `region_permutation_test` under every split of the mice kept into groups of
  their sizes (252 for 5 and 5, 126 for 5 and 4). Per fold: the heaviest
  cluster where |L - R| is higher after RWS, its voxels, mass and overlap with
  the cluster of all the mice; its p, the share of the splits whose positive
  cluster is as heavy; and step 3's p of either sign. Fold 0 gives step 3's
  cluster and p (3,662 voxels, mass 8,774.95, p 0.0397), and every fold's
  cluster is the leave-one-out's of `run_per_mouse_values`, voxel for voxel.
  As there, step 2 is not redone without the mouse.
- Readings. raw: |L - R| / (L + R) of the collected stack less the mouse's
  off-tissue level, smoothed as step 3 (`run_per_mouse_values`' raw stack);
  auto: the same on the autofluorescence; test: |L - R| on the test's maps,
  what step 3's t compares (their L + R has the normalisation's zero, so no
  ratio). Each mouse's value over the cluster is its selection-matched value
  of `run_per_mouse_values`, checked to the bit.
- Maps and profiles. In a window around the cluster (its bounding box and 600
  um on each side): a coronal view, the ratio of the means over the cluster's
  planes; a view from above, the ratio of the means over the cluster's depth
  below the pia (102 to 311 um) in each column; a profile along AP, the ratio
  of the means over the cluster's coronal footprint in each plane. One colour
  scale for the raw and the autofluorescence asymmetry. The depth below the
  pia is the distance to the atlas brain's outer edge, its unlabelled holes
  filled. The sections: the registered volume holds each section over the
  planes nearest it, so the collected stack, unsmoothed, changes from one plane
  to the next only between two sections; the profiles draw that change.
- Counts. At each of the 2,751,666 voxels of SSp-bfd where all ten mice have a
  value (of 3,143,890; all of the cluster's), the number of RWS mice above
  every naive mouse, 0 to 5; for each
  level k, the largest connected patch with at least k; the same under each
  of the 252 splits, its p the share reaching the observed patch. Every level
  is given. On exchangeable mice a voxel has all five with probability 1/252,
  four or more with 1/42, so the patch asks whether those voxels lie
  together.
- Size and peaks. The cluster's extent from the first voxel's edge to the
  last's; its widths along the surface (whose normal is the depth's mean
  gradient over the cluster) and its thickness through it; against a barrel
  column of about 300 um, its layer 2/3 from 128 to 418 um below the pia
  (Lefort et al. 2009, from slices; the CCF's layers need not have the same
  depths). Each mouse's peak raw asymmetry over the complete voxels of the
  region, and each RWS mouse's peak excess over the highest naive mouse.

### Results

Each mouse left out (`Influence_Folds_<tag>`): the cluster's voxels, the
share of the cluster of all the mice it keeps, its mass, its p over every
split of the mice kept, and step 3's p of either sign (the larger sign naive
higher where marked).

| left out | voxels | share kept | mass | p | either sign |
|---|---|---|---|---|---|
| none | 3,662 | 1.00 | 8,775 | 0.020 | 0.040 |
| CGF027 | 3,031 | 0.82 | 7,166 | 0.040 | 0.079 |
| CGF028 | 2,638 | 0.67 | 6,008 | 0.063 | 0.095 |
| CGF033 | 959 | 0.24 | 2,226 | 0.20 | 0.024, naive higher |
| CGF034 | 3,652 | 0.82 | 8,440 | 0.032 | 0.079 |
| CGF035 | 3,908 | 0.90 | 9,315 | 0.024 | 0.032 |
| MG691 | 1,133 | 0.28 | 2,615 | 0.056 | 0.12, naive higher |
| MG692 | 1,878 | 0.38 | 4,102 | 0.048 | 0.087 |
| MG693 | 5,508 | 0.66 | 14,565 | 0.008 | 0.016 |
| MG736 | 1,216 | 0.33 | 2,731 | 0.10 | 0.079, naive higher |
| MG737 | 1,132 | 0.31 | 2,524 | 0.087 | 0.17 |

Every fold's cluster sits at planes 553 to 570, its centre within 85 um of
that of the cluster of all the mice. Without CGF033, the naive mouse lowest in
the RWS cluster (raw AI 0.028), the naive-higher cluster becomes the heavier
one, and is itself as rare as step 3's cluster (mass 12,696, p 0.024).

Each mouse in the cluster (`Influence_Mice_<tag>.csv`): its asymmetry, and the
share of the cluster's voxels where it is above every mouse of the other group
(chance 0.17 on exchangeable mice), on the raw stack, the test's maps and the
autofluorescence; its asymmetry along AP through the cluster, over the
cluster's planes (`Influence_Profiles_<tag>.csv`).

| mouse | raw AI | auto AI | above, raw | above, test | above, auto | in the cluster's planes |
|---|---|---|---|---|---|---|
| CGF027 | 0.055 | 0.043 | 0.03 | 0 | 0 | 0.052 |
| CGF028 | 0.030 | 0.070 | 0 | 0 | 0 | 0.031 |
| CGF033 | 0.028 | 0.145 | 0 | 0 | 0.13 | 0.033 |
| CGF034 | 0.057 | 0.074 | 0 | 0 | 0 | 0.058 |
| CGF035 | 0.080 | 0.035 | 0 | 0 | 0 | 0.079 |
| MG691 | 0.154 | 0.139 | 0.90 | 0.99 | 0.47 | 0.150 |
| MG692 | 0.122 | 0.176 | 0.72 | 0.98 | 0.81 | 0.121 |
| MG693 | 0.077 | 0.095 | 0.22 | 0.83 | 0.04 | 0.078 |
| MG736 | 0.173 | 0.082 | 0.94 | 0.99 | 0 | 0.160 |
| MG737 | 0.145 | 0.117 | 0.85 | 0.99 | 0.14 | 0.141 |

Over the whole region the same share is near or below chance for every RWS
mouse but MG691 (0.34; MG692 to MG737 0.08 to 0.16), and two naive mice are
above every RWS mouse more often than chance (CGF027 0.29, CGF035 0.22).

The counts (`Influence_Consistency_<tag>`): the largest patch where at least k
RWS mice are above every naive mouse, its voxels in the cluster, the null's
median and 95th percentile over the 252 splits, and the p.

| reading | k | patch | in the cluster | null median | null 95% | p |
|---|---|---|---|---|---|---|
| raw | 3 | 130,480 | 3,110 | 110,330 | 327,550 | 0.41 |
| raw | 4 | 13,605 | 2,690 | 13,546 | 52,920 | 0.50 |
| raw | 5 | 1,120 | 0 | 1,794 | 6,294 | 0.74 |
| test | 3 | 144,400 | 3,662 | 91,826 | 372,640 | 0.34 |
| test | 4 | 32,753 | 3,632 | 12,600 | 57,989 | 0.13 |
| test | 5 | 5,400 | 2,865 | 1,835 | 6,818 | 0.079 |
| auto | 3 | 205,670 | 0 | 86,784 | 357,410 | 0.17 |
| auto | 4 | 23,555 | 0 | 16,316 | 66,997 | 0.35 |
| auto | 5 | 4,317 | 0 | 2,656 | 9,346 | 0.26 |

At k = 1 and 2 the patches span about half and a fifth of the region in every
split (raw p 0.41 and 0.58, test 0.27 and 0.26, auto 0.087 and 0.075). On the
raw stack 482 of
the cluster's voxels have all five RWS mice above every naive mouse, in
pieces; the five-of-five patch of the raw stack lies elsewhere.

The cluster (`Influence_Cluster_<tag>.csv`): 3,662 voxels, 0.0037 mm^3; planes
554 to 570 (CCF 733 to 749), 3.2 mm from the midline; 170 um along AP, 180
along DV, 220 along ML; along the surface 169 um (AP) by 186 um, 225 um thick
through it (the surface faces up and 31 degrees laterally there); 102 to 311
um below the pia, median 197; 95% layer 2/3, 5% layer 1. Against a barrel
column: its widest is 0.62 of the column's 300 um, its volume 0.18 of the
column's layer 2/3 (0.020 mm^3), its depth the bottom of layer 1, layer 2 and
the top of layer 3 by Lefort et al.'s depths (L1 to 128 um, L2 to 269, L3 to
418). The
CCF has no barrels, so which barrel it is over is not known; the whisker
stimulated was C2 (possibly B2).

- The sections. In the collected stacks the sections are mostly 14 to 18
  planes apart (140 to 180 um). The cluster's planes are one section of
  MG736 (its edges at 556 and 570, its peak at plane 563, the section's
  middle) and of MG737 (552 to 570); MG691's and MG692's elevation spans
  several sections, MG693's is two sections behind. MG736's asymmetry has the
  shape a difference between the two halves of one section would have (one
  section, across the whole window); that its autofluorescence lacks it says
  the difference is in the nano signal of that section, whether from the
  receptor, the staining or the imaging.
- Peaks. Each mouse's highest raw asymmetry over the region's complete voxels
  is 0.43 to 0.86, 600 to 1,290 um from the cluster: in layer 6a or 6b, 850 to
  1,130 um deep, in six mice; in layer 1 in three (MG692, MG693, CGF034); in
  layer 2/3 in CGF028. The RWS mice's largest excess over the highest naive
  mouse is in the same places (630 to 1,290 um away).
- Outputs, in the comparison's folder, `<tag>` the comparison and the
  smoothing: `Influence_Folds_<tag>` (.csv, one row per fold, and the figure:
  voxels, mass and p with each mouse left out), `Influence_Maps_<tag>_raw`,
  `_auto` and `_test` (each mouse's coronal view and view from above, the
  cluster outlined), `Influence_Profiles_<tag>` (.csv, one row per plane, and
  the figure: each mouse along AP with its section edges),
  `Influence_Consistency_<tag>` (.csv, one row per reading and level, and the
  figure: the counts' views, the patches against every split, each mouse's
  share of the cluster), `Influence_Mice_<tag>.csv`,
  `Influence_Cluster_<tag>.csv`, and the cache of the folds and the counts
  under every split, `Influence_<tag>.mat` (read when made from the same maps,
  autofluorescence and settings; `force_recompute` redoes it). For ten mice
  about 25 minutes for the folds and 3 for the counts, a few minutes from the
  cache. Not run after behavior.

## Where the code is

Each driver sets its settings and calls one function in `pipeline/`
(`collect_by_group`, `normalise_groups`, `group_differences`,
`per_mouse_region_values`, `mouse_influence`); the video
writers (`write_lr_*`, `open_lr_video`, `set_lr_colormap`) and the atlas
outlines (`lr_atlas_boundaries`) sit beside them, and so do the parts of
step 3 that other code reuses: each mouse's tissue and its smoothing
(`tissue_only`), its plane profile (`plane_tissue_means`), the alignment of
the two groups (`align_exp_to_ctrl`), the regions of the bars
(`surprise_regions`, `surprise_region_acronyms`) and their permutation test
(`region_permutation_test`); and the parts the per-mouse steps share: the box
around the clusters' region (`region_box_masks`), a mouse's box on the test's
scale (`aligned_box_lr`), the mice's |L - R| on the band and its candidate
voxels (`band_stack`, `band_geometry`), the check of the maps' cache against
the normalised stacks (`check_maps_cache`) and the mice's places in a
collected stack (`collected_index`). Shared: `../common/`
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
- Step 1 stacks the nano channel, and the autofluorescence when asked
  (`channels = {'auto'}`); each stack it writes holds its mice and the
  registered files they were read from (`collected_mice`, `source_files`).
  The `auto_4d.mat`, `mask_4d.mat` and average files in the adults'
  production group folders are older (19 November 2025) and come from an
  earlier registration: not to be read (`collect_by_group`'s help);
  `run_per_mouse_values` refuses an `auto_4d.mat` without `source_files`.
  The Python route reads the autofluorescence from the registered tiffs.
- `../adult_matlab/run_characterize_distribution`, kept until the Python
  route replaces it, reads the normalised stacks of step 2; the other two
  scripts of `../adult_matlab/` read its outputs.
