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
  8,775), planes 554 to 570 (CCF index 732 to 748, the posterior barrel field),
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
  486 (CCF index 649 to 664), 3.0 mm from the midline, all in layer 6a. Every
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

- Which mice. Not one RWS mouse's bump: four RWS mice, MG691, MG692, MG736
  and MG737, are high in the same patch of upper layer 2/3, about 200 um
  across, at the back of the barrel field, where the five naive mice are low
  and alike. Without any one of the four, 28 to 38% of the cluster of all
  ten mice (3,662 voxels) is left, 30 to 42% with the alignment held; its p
  over every split of the nine mice left (126) is then 0.048 to 0.10. MG693,
  the fifth, is low there, and without it the cluster grows (5,508 voxels).
- One naive mouse counts too, through the alignment, not through its
  asymmetry. Without CGF033 only 24% of the cluster is left, but with the
  alignment held 83% is, as without the other naive mice (57 to 76%). Step 3
  puts the RWS mice on the naive mice's scale with the line fitted from the
  RWS group's mean tissue profile along AP onto the naive group's; on the
  test's maps its slope multiplies every RWS mouse's |L - R| (the intercept
  cancels in L - R, the common factor in the t). After step 2, CGF033's tissue
  profile is 1.6 to 2.3 times as bright as any other mouse's and varies 2.1
  to 3.8 times as much along AP (mean 6,718 against 2,982 to 4,126, standard
  deviation along AP 3,926 against 1,047 to 1,887), and it raises the slope
  from 1.11 without it to 1.49 with it. Why step 2, which fits each mouse's
  cortex onto its group's median mouse, leaves CGF033's whole-plane profile
  this far from the others is not looked into here.
- Why a few mice are enough for step 3's test while their means are not. The
  cluster's mass sums the surprise of the contiguous voxels where the Welch t
  passes p < 0.01, so it grows with the extent of the place where the RWS mice
  are high and the naive mice low and alike. The same search in any split of
  the mice finds the place where its labelled groups are most apart, and they
  are apart there about as much (+0.080 at the null's median, +0.084
  observed): what is rare under relabelling is the extent (5 of 252 splits as
  heavy). And the extent depends on the slope, which the splits keep fixed:
  with all ten mice, the RWS-higher cluster grows from 549 voxels (p 0.23) at
  slope 1, each group on its own scale of step 2, to 3,662 (p 0.020) at the
  fitted 1.49 and 5,165 (p 0.016) at 1.7, while the naive-higher cluster
  shrinks; the RWS-higher one is the heavier only above a slope between 1.3
  and 1.4, and at slope 1 the naive-higher cluster is as rare as a split can be
  (p 0.004, either sign 0.008). The slope's spread without one mouse (1.11 to
  1.67, jackknife standard error 0.45) spans that crossing, so step 3's p
  0.040 holds at the fitted slope, not across the slope's own uncertainty.
  One figure shows it, `Rescaling_Test_<tag>` (`run_rescaling_figure`, from
  this step's tables and the scale-free test's): each sign's p and cluster
  size against the slope, the fitted slope with its standard error, the slope
  without CGF033, and beside them the scale-free test's p (0.32).
- Whether the RWS mice agree on the place more than relabelled groups do: no.
  On the raw stack, which has no slope, the largest patch where four of the
  five RWS mice are above every naive mouse is as large as under relabelling
  (13,605 voxels, median 13,546, p 0.50), and where all five are, 1,120
  voxels away from the cluster (p 0.74). On the test's maps the five-of-five
  patch looks borderline (5,400 voxels, 2,865 in the cluster, p 0.079), but
  that is the slope again: it lifts all five RWS mice above every naive mouse
  across the whole region, at 1.61 times the voxels expected on exchangeable
  mice, and at slope 1 the patch's p is 0.83 (0.78 for four of five).
- In the autofluorescence the RWS mice are the more asymmetric across the
  barrel field, a tendency rather than a result: at every level the voxels
  where the RWS mice are above every naive mouse are 1.2 to 2.2 times those
  expected (p 0.08 to 0.15 under relabelling), but not in the cluster.
- Each of the four is asymmetric there in its own way
  (`Influence_Maps_<tag>_raw`, `Influence_Profiles_<tag>`). MG736 in about
  one section (its section edges within two planes of the cluster's), across
  the whole window along ML, and its autofluorescence does not have it (0.082
  in the cluster, nano 0.173). MG737 in a superficial patch over that section
  and the next behind it. MG691 at 0.154 in the cluster (0.12 to 0.16 along
  its planes), about as its autofluorescence (0.139); its higher values 250
  to 460 um in front of the cluster's centre are the edges of a gap where the
  footprint has no tissue (planes 517 to 536), as high in the
  autofluorescence. MG692 moderately along 500 um of AP, less than its
  autofluorescence (0.122, autofluorescence 0.176). MG693's own bump is 400
  um behind the cluster. The naive mice are at 0.028 to 0.080
  in the cluster.
- Its size: 0.0037 mm^3, 170 um along AP, 220 um along ML, 102 to 311 um
  below the pia (95% layer 2/3, 5% layer 1), about two thirds of a barrel
  column's width (300 um, the mouse's C2 column, Lefort et al. 2009) and a
  fifth of the column's layer 2/3. Along AP it is about one section: the
  sections of these brains are mostly 14 to 18 planes apart, and the
  cluster's 17 planes lie within a plane or two of one section of MG736, of
  MG737 and of the naive CGF027.
- Open: why CGF033's profile is where it is after step 2. Step 3's cluster
  test on a reading without the slope, the raw asymmetry index, is the
  scale-free test below (10 October 2026): there the barrel-field cluster is
  not unusual (one-sided p 0.32).

### How it is done

- Folds. Fold 0 keeps every mouse; then each mouse is left out in turn and the
  comparison is redone as step 3 does it (the alignment refitted on the mice
  kept; a t where each group has three mice with a value; its surprise; the
  median over +/- 10 planes; p < 0.01; 18-connected within SSp-bfd), by
  `region_permutation_test` under every split of the mice kept into groups of
  their sizes (252 for 5 and 5, 126 for 5 and 4). Per fold: the heaviest
  cluster where |L - R| is higher after RWS, its voxels, mass and overlap with
  the cluster of all the mice; its p, the share of the splits whose positive
  cluster is as heavy; the p of the naive-higher cluster the same way; and
  step 3's p of either sign. Fold 0 gives step 3's cluster and p (3,662
  voxels, mass 8,774.95, p 0.0397), and every fold's cluster is the
  leave-one-out's of `run_per_mouse_values`, voxel for voxel. As there, step 2
  is not redone without the mouse.
- Folds with the alignment held. The same folds with the alignment of all
  ten mice (slope 1.4919) instead of the refitted one: what a mouse's own
  values do, apart from its part in the slope. Fold 0 is checked to be the
  same either way.
- Slopes. The comparison of all ten mice redone with the alignment's slope at
  1 to 1.7 in steps of 0.1 (`sweep_slopes`: from each group on its own scale
  of step 2 to past the slopes refitted without one mouse, 1.11 to 1.67),
  under every split, both signs. Step 3's splits relabel the mice after the
  alignment, so they keep its slope (Notes, the region measures); the slopes
  show what that fixed slope does. The jackknife standard error of the slope,
  from the folds' refitted slopes s_i, is sqrt(9/10 * sum((s_i - mean(s))^2)).
- Readings. raw: |L - R| / (L + R) of the collected stack less the mouse's
  off-tissue level, smoothed as step 3 (`run_per_mouse_values`' raw stack);
  auto: the same on the autofluorescence; test: |L - R| on the test's maps,
  what step 3's t compares (their L + R has the normalisation's zero, so no
  ratio); unscaled, for the counts: the test's maps at slope 1. Each mouse's
  value over the cluster is its selection-matched value of
  `run_per_mouse_values`, checked to the bit.
- Maps and profiles. In a window around the cluster (its bounding box and 600
  um on each side): a coronal view, the ratio of the means over the cluster's
  planes; a view from above, the ratio of the means over the cluster's depth
  below the pia (102 to 311 um) in each column; a profile along AP, the ratio
  of the means over the cluster's coronal footprint in each plane, with the
  share of the footprint where the mouse has tissue. One colour scale for the
  raw and the autofluorescence asymmetry. The depth below the pia is the
  distance to the atlas brain's outer edge, its unlabelled holes filled. The
  sections: the registered volume holds each section over the planes nearest
  it, so the collected stack, unsmoothed, changes from one plane to the next
  only between two sections; the profiles draw that change.
- Counts. At each of the 2,751,666 voxels of SSp-bfd where all ten mice have a
  value (of 3,143,890; all of the cluster's), the number of RWS mice above
  every naive mouse, 0 to 5. On exchangeable mice a voxel has at least k with
  probability C(5, k) / C(10, k) (1/252 for all five, 1/42 for four or more),
  which is also the share averaged over every split (checked in the run). For
  each level k: the voxels with at least k against that expectation, which a
  group more asymmetric everywhere exceeds; and their largest connected patch,
  which asks whether they lie together. Each against the same under each of
  the 252 splits, its p the share reaching the observed one. Every level is
  given.
- Size and peaks. The cluster's extent from the first voxel's edge to the
  last's; its widths along the surface (whose normal is the depth's mean
  gradient over the cluster) and its thickness through it; against a barrel
  column of about 300 um, its layer 2/3 from 128 to 418 um below the pia
  (Lefort et al. 2009, from slices; the CCF's layers need not have the same
  depths). Its planes also as Allen's CCF index, counted from 0. Each
  mouse's peak raw asymmetry over the complete voxels of the region, and each
  RWS mouse's peak excess over the highest naive mouse: one voxel's maximum of
  a ratio, which lands where L + R is low, so it says only where the cluster
  is not.

### Results

Each mouse left out (`Influence_Folds_<tag>`): the fold's refitted slope; the
cluster's voxels, the share of the cluster of all the mice it keeps, its p
over every split of the mice kept, and step 3's p of either sign (* where the
larger sign is naive higher); then, with the alignment held at 1.49, the
cluster's voxels, share kept and p.

| left out | slope | voxels | kept | p | either sign | held: voxels | kept | p |
|---|---|---|---|---|---|---|---|---|
| none | 1.49 | 3,662 | 1.00 | 0.020 | 0.040 | 3,662 | 1.00 | 0.020 |
| CGF027 | 1.67 | 3,031 | 0.82 | 0.040 | 0.079 | 2,097 | 0.57 | 0.048 |
| CGF028 | 1.49 | 2,638 | 0.67 | 0.063 | 0.095 | 2,638 | 0.67 | 0.063 |
| CGF033 | 1.11 | 959 | 0.24 | 0.20 | 0.024 * | 4,221 | 0.83 | 0.024 |
| CGF034 | 1.57 | 3,652 | 0.82 | 0.032 | 0.079 | 2,492 | 0.66 | 0.032 |
| CGF035 | 1.62 | 3,908 | 0.90 | 0.024 | 0.032 | 3,001 | 0.76 | 0.024 |
| MG691 | 1.44 | 1,133 | 0.28 | 0.056 | 0.12 * | 1,270 | 0.32 | 0.048 |
| MG692 | 1.42 | 1,878 | 0.38 | 0.048 | 0.087 | 2,210 | 0.42 | 0.048 |
| MG693 | 1.62 | 5,508 | 0.66 | 0.008 | 0.016 | 4,671 | 0.62 | 0.008 |
| MG736 | 1.49 | 1,216 | 0.33 | 0.10 | 0.079 * | 1,231 | 0.33 | 0.10 |
| MG737 | 1.50 | 1,132 | 0.31 | 0.087 | 0.17 | 1,108 | 0.30 | 0.087 |

Every fold's cluster sits at planes 553 to 570, its centre within 85 um of
that of the cluster of all the mice. Without CGF033 and the slope refitted,
the naive-higher cluster becomes the heavier one, and is itself as rare as
step 3's cluster (mass 12,696, p 0.024); with the slope held it stays the
lighter (7,609 against 9,832).

All ten mice at each slope of the alignment (`Influence_Slopes_<tag>.csv`,
and the folds' figure): each sign's heaviest cluster, its p over the 252
splits, and step 3's p of either sign.

| slope | RWS higher: voxels | mass | p | naive higher: mass | p | either sign |
|---|---|---|---|---|---|---|
| 1.0 | 549 | 1,254 | 0.23 | 10,731 | 0.004 | 0.008 * |
| 1.1 | 953 | 2,181 | 0.13 | 9,489 | 0.016 | 0.032 * |
| 1.2 | 1,634 | 3,766 | 0.071 | 8,361 | 0.028 | 0.056 * |
| 1.3 | 2,185 | 5,159 | 0.056 | 7,347 | 0.032 | 0.063 * |
| 1.4 | 2,842 | 6,799 | 0.040 | 6,259 | 0.044 | 0.079 |
| 1.49, fitted | 3,662 | 8,775 | 0.020 | 5,381 | 0.048 | 0.040 |
| 1.6 | 4,552 | 10,958 | 0.016 | 4,304 | 0.056 | 0.032 |
| 1.7 | 5,165 | 12,538 | 0.016 | 3,439 | 0.071 | 0.032 |

The RWS-higher cluster keeps its place at every slope (planes 553 to 570; at
1.5 and above it holds all of the fitted slope's cluster, below it lies
within it).

Each mouse (`Influence_Mice_<tag>.csv`): its tissue profile over the
alignment's planes (mean and standard deviation along AP); its asymmetry in
the cluster, raw and autofluorescence; and the share of the region's voxels
where it is above every mouse of the other group (chance 0.17 on exchangeable
mice), raw and autofluorescence, and of the cluster's on the raw stack (no
chance level: the cluster was chosen on these mice); its asymmetry along AP
through the cluster in `Influence_Profiles_<tag>.csv`.

| mouse | tissue mean | tissue sd | raw AI | auto AI | above, region, raw | above, region, auto | above, cluster, raw |
|---|---|---|---|---|---|---|---|
| CGF027 | 2,982 | 1,047 | 0.055 | 0.043 | 0.29 | 0.11 | 0.03 |
| CGF028 | 4,126 | 1,887 | 0.030 | 0.070 | 0.05 | 0.06 | 0 |
| CGF033 | 6,718 | 3,926 | 0.028 | 0.145 | 0.10 | 0.09 | 0 |
| CGF034 | 3,212 | 1,513 | 0.057 | 0.074 | 0.15 | 0.09 | 0 |
| CGF035 | 3,011 | 1,274 | 0.080 | 0.035 | 0.22 | 0.19 | 0 |
| MG691 | 3,539 | 1,075 | 0.154 | 0.139 | 0.34 | 0.19 | 0.90 |
| MG692 | 3,310 | 1,058 | 0.122 | 0.176 | 0.11 | 0.31 | 0.72 |
| MG693 | 3,977 | 1,695 | 0.077 | 0.095 | 0.08 | 0.16 | 0.22 |
| MG736 | 3,657 | 1,249 | 0.173 | 0.082 | 0.14 | 0.21 | 0.94 |
| MG737 | 3,981 | 1,313 | 0.145 | 0.117 | 0.16 | 0.26 | 0.85 |

The counts (`Influence_Consistency_<tag>`): at each level k, the voxels where
at least k RWS mice are above every naive mouse over those expected on
exchangeable mice, with its p; and the largest patch of them, its voxels in
the cluster, the null's median over the 252 splits, and its p.

| reading | k | voxels / expected | p | patch | in the cluster | null median | p |
|---|---|---|---|---|---|---|---|
| raw | 3 | 0.92 | 0.52 | 130,482 | 3,110 | 110,330 | 0.41 |
| raw | 4 | 0.97 | 0.46 | 13,605 | 2,690 | 13,546 | 0.50 |
| raw | 5 | 0.57 | 0.69 | 1,120 | 0 | 1,794 | 0.74 |
| test | 3 | 1.20 | 0.30 | 144,398 | 3,662 | 91,827 | 0.34 |
| test | 4 | 1.33 | 0.23 | 32,753 | 3,632 | 12,600 | 0.13 |
| test | 5 | 1.61 | 0.16 | 5,400 | 2,865 | 1,835 | 0.079 |
| slope 1 | 3 | 0.36 | 0.91 | 25,115 | 2,849 | 82,552 | 0.85 |
| slope 1 | 4 | 0.39 | 0.85 | 5,659 | 1,814 | 12,148 | 0.78 |
| slope 1 | 5 | 0.31 | 0.84 | 837 | 0 | 1,901 | 0.83 |
| auto | 3 | 1.60 | 0.12 | 205,670 | 0 | 86,784 | 0.17 |
| auto | 4 | 1.59 | 0.15 | 23,555 | 0 | 16,316 | 0.35 |
| auto | 5 | 2.20 | 0.083 | 4,317 | 0 | 2,656 | 0.26 |

At k = 1 and 2 the patches span about half and a fifth of the region in every
split; the autofluorescence's voxels are 1.23 and 1.48 times those expected
there (p 0.087 and 0.083; patches p 0.087 and 0.075), the raw stack's 1.03 and
0.96. On the raw stack 482 of the cluster's voxels have all five RWS mice
above every naive mouse, in pieces; the five-of-five patch of the raw stack
lies elsewhere.

The cluster (`Influence_Cluster_<tag>.csv`): 3,662 voxels, 0.0037 mm^3; planes
554 to 570 of the volumes (Allen's CCF index 732 to 748, counted from 0),
3.2 mm from the midline; 170 um along AP, 180 along DV, 220 along ML; along
the surface 169 um (AP) by 186 um, 225 um thick through it (the surface faces
up and 31 degrees laterally there); 102 to 311 um below the pia, median 197;
95% layer 2/3, 5% layer 1. Against a barrel column: its widest is 0.62 of the
column's 300 um, its volume 0.18 of the column's layer 2/3 (0.020 mm^3), its
depth the bottom of layer 1, layer 2 and the top of layer 3 by Lefort et al.'s
depths (L1 to 128 um, L2 to 269, L3 to 418). The CCF has no barrels, so which
barrel it is over is not known; the whisker stimulated was C2 (possibly B2).

- The sections. In the collected stacks the sections are mostly 14 to 18
  planes apart (140 to 180 um). The plane-to-plane change, spread over one or
  two planes at each edge, puts one section of MG736 at about planes 556 to
  569, one of MG737 at 553 to 570 and one of the naive CGF027 at 554 to 569:
  the cluster's 17 planes (554 to 570) are about one section, not exactly any
  mouse's. MG736's asymmetry, highest at plane 563, the section's middle, has
  the shape a difference between the two halves of one section would have
  (one section, across the whole window); that its autofluorescence lacks it
  says the difference is in the nano signal of that section, whether from
  the receptor, the staining or the imaging. MG691's and MG692's elevation
  spans several sections, MG693's is two sections behind.
- MG691's gap. Over the cluster's footprint MG691 has no tissue at planes 517
  to 536 (the profiles' opaque grey, the grey patch of its view from above),
  and only 18% at plane 516 and 1 to 72% at 537 to 539. Its raw asymmetry is
  0.31 and 0.33 on the two edges, its autofluorescence 0.25 and 0.34, and
  falls to 0.16 by the cluster: the asymmetry of the tissue's edge, in both
  channels.
- Peaks. Each mouse's highest raw asymmetry over the region's complete voxels
  is 0.43 to 0.86, 600 to 1,290 um from the cluster: at the region's edge in
  seven mice, in layer 6b by the white matter (1,047 to 1,126 um deep:
  CGF027, CGF035, MG691, MG736) or in layer 1 by the pia (44 to 94 um:
  CGF034, MG692, MG693); in upper layer 6a in two (850 and 880 um: MG737,
  CGF033); in layer 2/3 in CGF028. The RWS mice's largest excess over the
  highest naive mouse is in the same places (630 to 1,290 um away). Single
  voxels of a ratio, these land where L + R is low and say only that the
  cluster is no mouse's highest point.
- Outputs, in the comparison's folder, `<tag>` the comparison and the
  smoothing: `Influence_Folds_<tag>` (.csv, one row per fold, refitted and
  held, and the figure: voxels, mass and p with each mouse left out, each sign
  at each slope, and the notes), `Influence_Slopes_<tag>.csv` (one row per
  slope), `Influence_Maps_<tag>_raw`, `_auto` and `_test` (each mouse's
  coronal view and view from above, the cluster outlined),
  `Influence_Profiles_<tag>` (.csv, one row per plane, and the figure: each
  mouse along AP with its section edges and where it has no tissue),
  `Influence_Consistency_<tag>` (.csv, one row per reading and level, and the
  figure: the counts' views, the patches against every split, each mouse's
  share of the region and of the cluster), `Influence_Mice_<tag>.csv`,
  `Influence_Cluster_<tag>.csv`, and the cache of the folds, the slopes and
  the counts under every split, `Influence_<tag>.mat` (read when made from the
  same maps, autofluorescence and settings, a missing part computed and
  added; `force_recompute` redoes it). For ten mice without a pool about 25
  minutes for the folds, as much for the held folds, 20 for the slopes and 4
  for the counts; with a thread pool of 8 (`n_workers`), 6, 4 and 4 minutes
  for the held folds, the slopes and the counts; a few minutes from the
  cache. Not run after behavior.

## The scale-free test

In short (10 October 2026). On a reading no scale can change, the raw
asymmetry index, step 3's test does not find the barrel-field result after
RWS: SSp-bfd's heaviest RWS-higher cluster is 332 voxels on the edge of step
3's patch, as heavy as in 80 of the 252 splits (one-sided p 0.32, the test
named in advance; either sign 0.56). No region passes correction, after RWS
or after behavior. This agrees with the slope sweep of the section above,
which is what shows that step 3's p 0.040 rests on the alignment's slope; the
index alone does not single out the slope, since it also leaves out each
mouse's line of step 2 and divides each voxel by its own L + R.

The design below was fixed on 10 October 2026, before any of its numbers
(`run_scale_free_test`). Its p are reported whatever they are; no other
variant is tried.

Why. Step 3's barrel-field result after RWS depends on the alignment of the
RWS group onto the naive group (the section above). The fitted slope, 1.49,
set mostly by CGF033, multiplies every RWS mouse's |L - R| on the test's maps.
At slope 1 the RWS-higher cluster is small (p 0.23) and the naive-higher one
is the heavier (p 0.004); the signs cross between 1.3 and 1.4, inside the
slope's own spread (jackknife standard error 0.45). The decisive check is the
same test on a reading no scale can change.

- The reading. Per mouse and voxel, the asymmetry index
  AI = |L - R| / (L + R) on the collected stack, the mouse's off-tissue level
  subtracted, as `run_per_mouse_values` reads its raw stack. Each mouse's
  tissue is smoothed before the ratio exactly as step 3 smooths it
  (`tissue_only`: NaN-aware, the same Gaussian), then folded
  (`compute_lr_stats`). A voxel is missing where L + R is not above 0 or
  outside the mouse's tissue. No normalisation of step 2 and no alignment of
  the groups: a factor on a mouse or on a group divides out of the ratio. The
  folded grid, the regions and the masks are step 3's.
- The statistics. Step 3's region test on these maps, run by step 3's own
  code (`permutation_stacks`, `region_permutation_test`, `plot_measure_bars`):
  the Welch t per voxel, RWS against naive, a t only where each group has at
  least three mice with a value; its surprise -log10 p; the median over +/- 10
  planes; clusters at p < 0.01, 18-connected, within each region's own voxels
  of the left hemisphere. Each region is scored by its heaviest cluster's mass
  with its sign, RWS higher positive, the heavier of the two signs. The exact
  permutation over all 252 splits of the ten mice, everything recomputed in
  each; a region's p is the share of the splits reaching its |score|, its
  corrected p the share whose largest |score| over the regions reaches it.
- The test named in advance. SSp-bfd, one-sided, RWS higher (Gambino et al.
  2014): the mass of its heaviest RWS-higher cluster against the RWS-higher
  masses of all 252 splits, the observed split included. The p of either sign
  beside it. Also where SSp-bfd's heaviest RWS-higher and naive-higher
  clusters sit (voxels, planes as Allen's CCF index) and how they overlap
  step 3's cluster (3,662 voxels, CCF 732 to 748).
- After behavior. The same over the 126 splits of the four behavior mice of
  step 3, the p of either sign first, the direction carried over from RWS
  beside it; the overlap with step 3's behavior cluster (1,489 voxels).
- What it can say. An RWS-higher cluster that holds on the index does not come
  from the slope. One that does not leaves step 3's p 0.040 resting on the
  alignment. (Added at the review, 10 October 2026: the index differs from
  step 3's maps in more than the slope, since it also leaves out each mouse's
  line of step 2, divides each voxel by its own L + R and takes its zero from
  the off-tissue level. So a cluster that does not hold on the index agrees
  with the slope sweep above, which isolates the slope, but does not itself
  isolate it.)
- Checks in the run. Every mouse's index in the box around SSp-bfd is that of
  `run_per_mouse_values`' raw stack, to the bit, with the same off-tissue
  level. Step 3's cluster is the leave-one-out's fold 0 of
  `run_per_mouse_values` (step 3 keeps no voxels), checked against step 3's
  table. Since the review, also the selection-matched test's observed
  cluster (`Per_Mouse_Selection_<tag>.mat`), voxel for voxel; its mass and
  its p of either sign over that test's splits, which are step 3's, are
  checked against step 3's table, and give step 3's one-sided p.
- Outputs, in the comparison's folder, `<tag>` the comparison and the
  smoothing: `Scale_Free_Regions_<tag>.csv` (one row per region: each sign's
  heaviest cluster, its mass, voxels, peak, own one-sided p and planes; the
  score, p and corrected p), `Scale_Free_Clusters_<tag>.csv` (SSp-bfd's two
  clusters and step 3's: where each sits, the overlap; `p_one_sided`, the
  share of the splits whose cluster of the row's sign is as heavy, and
  `p_either_sign`, the region's p, for the heavier sign only; a sign without
  a cluster has mass 0, which every split reaches, so its p is 1),
  `Scale_Free_Bars_<tag>` (the bars in step 3's style, the test named in
  advance under the title, and in bold beside SSp-bfd's bar the p of that
  test: one-sided after RWS, of either sign after behavior),
  `Scale_Free_TMap_<tag>` (the index's t in SSp-bfd at each cluster's planes
  and from above, the clusters outlined), and the cache of the test,
  `Scale_Free_<tag>.mat` (read when made from the same stacks and settings;
  `force_recompute` redoes it).

### Results

10 October 2026, on copies of the stacks of the production run of 7 October
(the development root of the two sections above). Both checks held: every
mouse's index was that of the raw stack, to the bit, and the leave-one-out's
fold 0 was step 3's cluster (3,662 voxels, mass 8,774.95 after RWS; 1,489
voxels after behavior). At the review's rerun from the cache (10 October
2026), fold 0 was also the selection-matched test's cluster, voxel for voxel,
with step 3's mass and p of either sign (one-sided, 0.020 after RWS and 0.056
after behavior); the index's clusters and every p were those of the first
run. For ten mice the maps took 12 minutes and the 252
splits 18 on 16 thread workers (14 and 16 minutes for behavior's 126); a
couple of minutes from the cache.

After RWS, SSp-bfd (252 splits; planes as Allen's CCF index; the one-sided p
is that of the row's sign, the p of either sign the region's, given for the
heavier sign):

| cluster | voxels | mass | planes | one-sided p | p of either sign | in step 3's cluster |
|---|---|---|---|---|---|---|
| index, RWS higher | 332 | 727 | 735 to 747 | 0.32 (80 of 252 splits) | 0.56 | 74 voxels: 22% of it, 2% of step 3's |
| index, naive higher | none | 0 | | 1 (no cluster) | | |
| step 3, RWS higher (test maps) | 3,662 | 8,775 | 732 to 748 | 0.020 | 0.040 | |

- The test named in advance does not hold: the RWS-higher cluster's mass is
  reached by 80 of the 252 splits (p 0.32; the splits' median 155, quartiles
  9 and 1,022, 37 splits without one). With either sign the score is +727,
  p 0.56, corrected p 1.
- Where. The index's cluster lies at step 3's planes and depth, on the upper
  and medial edge of its patch (centre 3.09 mm from the midline against
  3.18, 1.03 mm from the top of the volume against 1.09). On the index the
  t is positive over the whole patch (mean 3.3 over step 3's voxels, which
  were chosen with these same mice, so it is no test), but the heaviest
  cluster where its rolled surprise reaches p < 0.01 is a tenth of step 3's
  in size and holds 74 of its voxels. The direction is kept, the extent that
  step 3's test rewards is not (`Scale_Free_TMap_<tag>`).
- The naive-higher cluster that is the heavier at slope 1 on the test's maps
  (p 0.004) has no counterpart on the index: nowhere in SSp-bfd does the
  surprise's median over +/- 10 planes, on which the clusters are drawn,
  reach p < 0.01 on the naive-higher side (its lowest, p 0.013). Single
  voxels do: with their own Welch df, 8,490 of the region's 3.06 million
  voxels with a t have a naive-higher t at p < 0.01 (4,803 RWS-higher), in
  patches of up to 2,533 voxels, but not over the median's 21 planes. Which
  of the index's differences from the test's maps removes the cluster (the
  slope, each mouse's line of step 2, the division by each voxel's L + R) is
  not tested here.
- Every region: none at a corrected p < 0.05 (the smallest 0.83). The
  heaviest scores are naive-higher in MOp (8,841, p 0.17), SSp-tr (7,166,
  p 0.016) and GENd (6,975, p 0.040), RWS-higher in ACAv (4,375, p 0.024)
  and VISp (3,239, p 0.13). Four regions have an uncorrected p below 0.05
  (SSp-tr, GENd, ACAv, ILA), none named in advance, so none is a test.

After behavior, SSp-bfd (126 splits):

| cluster | voxels | mass | planes | one-sided p | p of either sign | in step 3's cluster |
|---|---|---|---|---|---|---|
| index, behavior higher | 1,271 | 2,761 | 663 to 676 | 0.17 (21 of 126 splits) | 0.24 | none |
| index, naive higher | 254 | 574 | 617 to 624 | 0.37 | | none |
| step 3, behavior higher (test maps) | 1,489 | 3,454 | 649 to 664 | 0.056 | 0.15 | |

- The p named first, either sign: +2,761, p 0.24, corrected p 1. Behavior
  higher, the direction carried over from RWS: p 0.17.
- The index's behavior-higher cluster lies just behind step 3's and higher
  in the cortex (1.66 mm from the top of the volume against 2.31), with no
  voxel in common.
- Every region: none at a corrected p < 0.05 (the smallest 0.66); six with an
  uncorrected p below 0.05 (SSp-n, RSPv, VISrl, MD, PF, CL), none named in
  advance.

What it means for the RWS claim. Step 3's barrel-field cluster is rare under
relabelling only through its extent. On a reading free of any scale the RWS
mice are still the more asymmetric over that patch (chosen with them), but
the patch where the difference reaches p < 0.01 is no larger than relabelled
groups give. With these ten mice the index shows no RWS-specific change of
the nanobody signal's asymmetry in the barrel field: an absence of evidence,
not evidence of no change. That the extent comes from the slope of the
alignment is what the slope sweep of the influence section shows (at slope 1
the RWS-higher cluster has 549 voxels, p 0.23; step 3's p 0.040 holds only at
the fitted slope); the index agrees with it but does not isolate the slope,
since it also leaves out each mouse's line of step 2, divides each voxel by
its own L + R and takes its zero from the off-tissue level. The index shows
no change after behavior either, where step 3's own test did not pass.

## Where the code is

Each driver sets its settings and calls one function in `pipeline/`
(`collect_by_group`, `normalise_groups`, `group_differences`,
`per_mouse_region_values`, `mouse_influence`, `scale_free_test`,
`plot_rescaling_test`); the video
writers (`write_lr_*`, `open_lr_video`, `set_lr_colormap`) and the atlas
outlines (`lr_atlas_boundaries`) sit beside them, and so do the parts of
step 3 that other code reuses: each mouse's tissue and its smoothing
(`tissue_only`), its plane profile (`plane_tissue_means`), the alignment of
the two groups (`align_exp_to_ctrl`), the regions of the bars
(`surprise_regions`, `surprise_region_acronyms`), the mice's maps on the
voxels a split can give a t (`permutation_stacks`), their permutation test
(`region_permutation_test`) and the bars of a region measure
(`plot_measure_bars`, with `p_shade_index`, `add_p_colorbar`,
`highlight_surprise_regions` and `expected_regions`); and the parts the
per-mouse steps share: a mouse's off-tissue level (`off_tissue_level`), the
box around the clusters' region (`region_box_masks`), a mouse's box on the
test's scale (`aligned_box_lr`), the mice's |L - R| on the band and its
candidate voxels (`band_stack`, `band_geometry`), the check of the maps' cache
against the normalised stacks (`check_maps_cache`) and the mice's places in a
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

  A test of the rescaling (10 October 2026; the sections on each mouse's
  influence and on the scale-free test, below) shows that this p depends on
  the line that aligns the RWS group onto the naive one. Its slope, 1.49, is
  set mostly by one bright naive mouse (CGF033). Without the rescaling the
  barrel field gives p 0.23, and on the asymmetry ratio |L - R| / (L + R),
  which no scale changes, p 0.32. Giulio (10 October) keeps p 0.040 as the
  result to quote, with this dependence stated beside it.

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
