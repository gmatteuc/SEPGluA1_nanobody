# Roadmap

What comes next, in order, and why. Written on 3 and 4 October 2026, when
the refactor of [REFACTOR_PLAN.md](REFACTOR_PLAN.md) had done its steps 0 to 8
on the branch `refactor` and the documents and the merge were next. What each
line of work found is in [SCIENTIFIC_CONTEXT.md](SCIENTIFIC_CONTEXT.md); which
script makes each figure, in [FIGURES.md](FIGURES.md).

Paths below are relative to the code root or to the data root, as
`get_paths.m` and `mapping/sepmap/config.py` define them: the code folder and
the `data` folder beside it (`D:\sep_histology\code` and
`D:\sep_histology\data` on the analysis computer). `SEP_DATA_ROOT` moves the
data root, for checks on a copy only.

---

## Sequence

```
1.  Merge the refactor (step 11)                       [next]
2.  Rerun the production outputs with the merged code  [after 1]
3.  A1 to A5 in the Python route (step 9)              [after 1; A1 first]
4.  P8, P9 and P10 to archive/                         [once A1 to A5 answer their questions]
5.  The remaining young brains                         [beside 3; into the Python route after A1]
      MG914's SEP channel, then MG896, MG906, MG895, then MG907, MG908
6.  Young against adult across ages, P16 to P36        [once 5 is done]
7.  A6 to A10                                          [after 4]
8.  Code questions held for later                      [any time]
9.  LightSuite: pull requests upstream (Z2), then the swap (Z1)
10. The automatic annotation's validation sandbox into this repository
Open throughout: comparable normalisation of the two routes; the open scientific questions
```

Rationale: A1 is load-bearing. It fixes the reference of `zref`, which today
is the set of structures shared by every brain in the run, young brains
included, so each new young brain moves every adult's `zref`, and with it the
young-against-adult differences behind the grant figures. With A1 done first,
the six young brains add data without moving the numbers already quoted.
Their registration (MATLAB, with sittings at the annotation GUI) does not
depend on A1 and can run beside it.

---

## 1. Merge the refactor

Steps 0 to 8 are done (Progress in [REFACTOR_PLAN.md](REFACTOR_PLAN.md)).
The fixes Giulio accepted on 3 October (yes to 1 to 9 and 11 to 21, no to
10; the new fixes 23 to 27) are on `refactor`, each in its own commit, and
were rechecked against the reference: the plasticity chain, P2 and P2bis on
MG914, the whole Python route, and fix 23 on the inputs of the approved
run. Every output that changed is traced to its fix. Left before the merge:

- **The registration check of step 8** (MG914, 'register' then the SEP
  channel), which stopped on 3 October when drive C: ran out of space. To
  rerun before the merge, with room on C: for MATLAB's temporary files.
- **Decisions from step 8** (below): fix 23's edge voxels, the region list of
  the surprise bars, and the cohort videos' mean panel.
- The documents of step 10 and the final report.
- The merge itself, as the plan's "Merging" lists it: every MATLAB session
  closed; the uncommitted run settings of `P4_register_to_atlas.m` on `main`
  written down, stashed and re-applied to `registration/run_register_to_atlas.m`;
  the engine's environment moved by hand to `registration\auto_annotation\.venv`
  (git does not move ignored files) and its self-test run before the first
  `annotate`, otherwise the GUI opens without the automatic layer; the nested
  duplicate `matlab_elastix-master` deleted.

### The fixes that change results

Measured on the check trees on 4 October; the production outputs change only
at the rerun after the merge (section 2).

- **Subref's reference (fix 1).** It is now exactly TH, HY, PAL, MB, P and
  MY: fibre tracts, ventricles and the unassigned labels left out. Each
  brain's `subref` moves by one constant, +0.19 to +0.45 log2 in the adults
  (mean +0.34) and +0.02 to +0.32 in the young brains (mean +0.13). The
  young-against-adult `subref` difference moves by -0.21 log2 (-0.19 from
  the fibre tracts and ventricles, -0.02 from the unassigned labels), and its
  structures at Welch q < 0.05 go from 104 to 81 (Mann-Whitney 97 to 75).
  `ratio`, `sepratio`, `cref` and `zref` are unchanged.
- **Missing values left out (fixes 3 and 25).** Only MG897 has voxels without
  a `sepratio` value: 183 in the CCF maps, and in its own atlas 117, all in
  the nucleus accumbens. At the 183 the young `sepratio` mean rises by a
  factor of 7/6 (5/4 for `young_P20`) and its n falls from 7 to 6 (5 to 4);
  in the region tables one value moves in the fourth decimal. The P20-only
  contrast maps of `volumes_ccf20.npz` (`log2_alt_*`, drawn by no figure)
  are now smoothed without counting a missing value as zero.
- **The plasticity comparison's smoothing (fix 23).** Each mouse is NaN
  outside its tissue, the smoothing is normalised by the smoothed tissue
  mask, and every mean, SEM and t is taken over the mice with tissue on both
  sides of the voxel. Values inside the tissue are unchanged; t changes
  wherever some mouse lacks tissue. Rerun on the inputs of the approved run,
  as the S6 run of 1 October was (the check's tables are in
  `G:\sep_refactor\s6\scratch\step8b\`), naive against RWS, barrel field at
  slab 565:

  | | approved run | before the fix | after |
  |---|---|---|---|
  | difference map (L - R), pixels at p < 0.01 | 271 | 292 | 292, the same pixels |
  | sum map (L + R), pixels at p < 0.01 | 532 | 437 | 287, all among the 437 |
  | surprise bar, L - R | 14,155 (rank 1 of 36) | 14,562 (1 of 36) | 15,041 (1 of 42) |
  | surprise bar, L + R | 89,140 (rank 14 of 38) | 84,521 (13 of 38) | 81,773 (5 of 35) |
  | per mouse, L + R, naive and RWS: mean (SEM), Welch p, 5 and 5 mice | 1.389 (0.143) and 1.766 (0.161), p 0.118 | p 0.120 | 1.414 (0.135) and 1.776 (0.157), p 0.119 |

  - The RWS S1 increase holds. The difference map is unchanged in S1; the
    sum map's significant area there shrinks by about a third. The barrel
    field rises in the sum bars because large bars made of zeros counted as
    data disappear (midbrain, motor related 5.22M to 159; VISp 0.96M to
    4,208; HPF 0.91M to 14,493; RSP 588k to 18,180). Per mouse, the barrel
    field is 26% higher in RWS, not significant with five against five,
    before and after; no mouse's own value moves by more than 0.006.
  - Naive against behaviour changes most. MG709 has no tissue at slab 565,
    but the old code counted it as a fourth behaviour mouse. With the three
    real ones, the sum map is significant over about a quarter of the slab
    (76,879 pixels against 7,888; in the barrel field 2,944 against 1), and
    the difference map has 286 significant pixels, mostly positive, against
    14, all negative. Per mouse, the barrel field's L + R is 55% higher after
    behaviour (naive 1.372, SEM 0.131; behaviour 2.119, SEM 0.133; 5 and 3
    mice; Welch p 0.009), before and after the fix. That figure needs
    re-reading.
  - New with the fix: at tissue edges, voxels left with two or three mice
    per group give very large t (9 pixels of the RWS slab above |t| = 10, up
    to 73), and probably the new L - R bars at tissue edges (olfactory
    tubercle, orbital, VISp, HPF, midbrain motor), not checked directly.

### Decisions open from step 8

- **Fix 23's edge voxels.** A minimum number of mice per voxel (for example
  3 per group), and/or the surprise bars and the slab opacity summed only
  over voxels with a t of their own (today the rolling median of +/- 10
  planes fills a voxel without a t from its neighbours). Either moves the
  December comparison again.
- **The surprise bars' region list** (55 names, each an exact atlas name,
  none twice since fix 9). Nested regions count voxels twice: SUB inside HPF
  (4.9% of HPF), STN and ZI inside HY (1.3%, 11.8%), SCm inside MBmot
  (26.5%); keep both, drop one, or take the smaller out of the larger? The
  bars are sums, so large regions win over small nuclei (HPF 21.3M voxels,
  STN 0.1M): sum, mean, or fraction of voxels over the threshold? Was the
  second parafascicular entry meant to be another region (PVT, SSp-n, SSp-m,
  SSp-un, VISpm, VISpl and AUDpo are absent)?
- **The cohort videos' mean panel** (fix 24). The t is now over each brain's
  hemispheres averaged first; the mean panel still folds the cohort map, so
  it matches `run_compare` and `run_video_compare`. Taking the per-brain
  average there too is one line, and would move the mean by more than 1% in
  12 to 25% of voxels. Each frame's header still gives the larger side's
  count of mice, not the t's n.

## 2. Rerun the production outputs

The refactor's checks ran on copies, so after the merge every output under
the data root is still the old code's. The fixes change some of them (listed
per figure in [FIGURES.md](FIGURES.md)), so both routes are rerun once on the
production data:

- **Python route**, in the run order of the `mapping/run_*.py` headers, with
  the `<mouse>_scalars.npz` caches in `comparisons_v2\per_mouse\` deleted
  first: they are reused when their source file's date matches, and would
  keep the subcortex mean from before the subref fix. This rewrites
  `comparisons_v2\` and `adult_v2\`, the grant figures included.
  `run_cohort` now also writes the `<reading>_folded_{mean,sd,n}.npy` files
  the videos' t needs (about 26 GB more on D:, and 55 minutes instead of 36),
  and `run_video` stops until it has. The data
  tree as of 30 September is in the snapshot `G:\sep_histology_snapshot_2026-09-29`;
  copying `comparisons_v2\young_vs_adult\` aside under a dated name before
  the rerun keeps the grant state next to the new one.
- **Plasticity comparison.** The approved figures of 5 December 2025 are in
  `comparisons\naive_vs_rws\` and `comparisons\naive_vs_behavior\`. Today's
  code writes to `naive_vs_rws_nano\` and `naive_vs_behavior_nano\`, so it
  cannot overwrite them. To decide: rerun the whole chain (collect, normalise,
  differences) on today's code, or only the differences on the preserved
  inputs of the approved run (`nano_4d_normalized_bk.mat` of naive and rws,
  behavior's file of 5 December), as the informative S6 run did. The second
  isolates the effect of the fixes; the first is what a new animal would go
  through. Settle fix 23's edge voxels first (section 1): so far the fixed
  group-difference step has run only on the S6 copy on G:.

## 3. The scientific additions A1 to A5 (step 9)

The Python route does not yet answer everything P8, P9 and P10 answer. Of
their 88 outputs and tests, 28 are covered (usually more rigorously), 26
partly, 19 not at all; 4 were dropped by decision and 11 are file handling
([REFACTOR_COVERAGE.md](REFACTOR_COVERAGE.md), which also holds the
specification of each addition; where it differs from the plan, the plan
applies).

| | addition | the question | replaces | what it moves |
|---|---|---|---|---|
| A1 | one declared structure set, written once to a table every analysis reads, and the `zref` reference taken from it | which structures enter, and a `zref` that depends only on the brain and a fixed list | P8's z-score and region set | every `zref` value, the young-against-adult differences of the grant figures included; the enrichment calls |
| A2 | ISH quality control per section: failed sections flagged and set to missing, never interpolated; a reviewed list of genuine regional absence; the effect on the ranking reported | are the ISH structure means corrupted by failed or dim sections | P9's section repair | the gene ranks (an unreviewed scan flags 14 of 95 genes; Gria1 is clean) |
| A3 | the gene ranking on A1's set and A2's tables, through one shared adult profile with a minimum-mice rule; ranking, roles, panel test, arms and words rerun; robustness to the statistic and to borders; a sensitivity table under three structure sets | which genes predict the adult map, on structures the cohort actually measures | P9's region set, metric comparison and headline ranking | Gria1's rank and the ISH numbers held back under S5 |
| A4 | the adult distribution: every structure by division, per-mouse mean and SEM, reliability, the enrichment call, an eroded mean beside the plain one; for the young cohort on its own too | how the signal is distributed across the adult brain, and how reproducibly | P8's bar charts, tables and enrichment | the enrichment list (zero becomes the brain's median structure) |
| A5 | autofluorescence as a parallel control per structure and division, and its own distribution | is the nano pattern its own signal, structure by structure | P10, and P8's `auto` run | P10's lists (the question is kept, not the output) |

The decisions they rest on (S1 to S5 in the plan):

- **S1, the structure set.** Grey-matter structures seen in all ten adults,
  with every dropped structure and the reason. The rule keeps 234 of 280
  structures; today 22 structures of the adult table are seen in fewer than
  five adults, mostly pons, medulla and midbrain, and their values do not
  replicate across mice. P9's nine-division set is reported alongside.
- **S2, what "enriched" means.** Per structure, an exact signed-rank test of
  `zref` against the brain's median structure across mice, with
  Benjamini-Hochberg across structures beside the uncorrected p; the spatial
  null (A7) is the stronger version. At n = 10 the smallest two-sided exact p
  is 1/512. Zero now means the median structure, not the voxel mean as in P8:
  to agree with Sami El-Boustani before the figures are shown. If BH leaves
  nothing standing, the correction is reconsidered together (a declared set
  of structures from the grant: V1, S1, RL, AL, LI and the control regions;
  tests at division level; the spatial null), not applied blindly.
- **S3, nano against autofluorescence.** Two-sided tests on centred
  contrasts of the per-mouse log2 nano and autofluorescence structure means,
  never the stored `ratio` reading, whose level depends on exposure; each
  contrast shown beside the two values; BH within each family at structure
  level, Bonferroni over the nine divisions.
- **S5, which ISH numbers to quote.** None of Gria1's rank, the powered panel
  test or the roles permutation until A1 to A3 have run. Cacng8 first under
  every structure set tried holds already.

Also in step 9, because no output holds the numbers yet: the per-mouse
supragranular-minus-infragranular contrast of each cortical system goes into
`run_region_groups`' outputs. SCIENTIFIC_CONTEXT states the result from a
read-only check of 1 October and quotes no value until `run_region_groups`
writes it.

Each addition is committed with its before and after written down, and the
numbers marked provisional in [SCIENTIFIC_CONTEXT.md](SCIENTIFIC_CONTEXT.md)
and [adult_ish_design.md](adult_ish_design.md) are corrected in the same
step.

## 4. Retiring P8, P9 and P10

A MATLAB analysis retires when the Python route answers the same scientific
question, not necessarily with the same output. Once A1 to A5 are in:

- `P8_characterize_merged_distribution.m`, `P9_compare_nano_vs_allen_ish.m`,
  `P10_compare_nano_vs_auto.m` and `plot_violinplot.m` (used by P9 only) move
  to `archive/`; the override `SEP_MERGE_SPECS` goes with P8.
- Their outputs in `comparisons\` stay where they are, frozen, as what the
  additions were compared with (`merged_naive_rws_nano\`,
  `merged_naive_rws_auto\`, the `merged_naive_rws_vs_ish_*` folders,
  `nano_vs_auto\`). P9's `gene_panel_summary.csv` stays an input of
  `run_ish_compare` (the old-against-new ranking check).
- Giulio checks the Python route against the old code, run from the tag
  `refactor-start` in its own check tree; then `archive/` is deleted (L2), the
  tag keeping the files.

Differences from P8 to P10 caused by their own defects are expected, not
errors of the additions:

- P8's adult outputs carry an artefact of the `abs()` over slabs: P6bis was
  not rerun after the raw-zero fix.
- P9 stretched the 200 um ISH grid onto the atlas box (a 1.5 to 2.5% scale
  error); its missing-data value was interpolated and clamped to 0; the
  off-reference Olig2 and Calb2 grids were stretched.
- P9's section repair "repaired" a true absence of expression (Slc17a6) and
  copied the last good section over posterior planes (Cacng8).
- P8's and P10's autofluorescence came from a stale `auto_4d.mat` that the
  collect step no longer writes; P10's inputs mix caches of two dates.
- P10's Bonferroni over structures cannot reject at n = 10; P9's category
  ANOVA was anticonservative.
- P8's reliability bars run from white, so the least reliable structures are
  white bars on white; A4 does not copy the scale.

## 5. The remaining young brains

Six young brains are not yet in the analyses. Where each stood on 3 October,
read from its folder under `data\young\`:

| brain | age | atlas | reached | next |
|---|---|---|---|---|
| MG914 | P28 | `demba_p28` | registered (automatic annotation, reviewed, 30 September) | `run_add_sep_channel`, then the Python route |
| MG896 | P28 | `demba_p28` | aligned (slice order, residual correction and equalisation done) | cutting angle, anchors, `autoannotate`, review, `register`, SEP channel |
| MG906 | P32 | `demba_p32` | aligned | as MG896 |
| MG895 | P36 | `demba_p36` | aligned | as MG896 |
| MG907 | P36 | `demba_p36` | extracted; residual correction and equalisation done | slice order, `align`, then as MG896 |
| MG908 | P32 | `demba_p32` | extracted; residual correction and equalisation done | as MG907 |

The order is MG896, MG906, MG895, then MG907 and MG908 (the plan), so the
first three give one brain at each new age. The steps of a new brain are in
[ADDING_DATA.md](ADDING_DATA.md). MG912 (P20) is aligned but was dropped for
poor slice quality (`comparisons_v2\README.md`). The DeMBA atlases of P28,
P32 and P36 are built (`data\atlas_demba_p28\`, `_p32\`, `_p36\`).

Things to watch with these brains:

- **The registered grid.** Since fix 26 of step 8, `register` sets it
  itself: every brain is sampled on the 10 um grid of
  `registration/pipeline/registered_grid_um.m`, the grid the Python route
  expects, whatever `px_atlas` the brain's `sliceinfo.mat` holds. It stops
  before registering any brain whose `px_register` there is not 20 (LightSuite
  places the slices on a grid of half `px_register`), and only prints a line
  for a `px_atlas` other than 10. Every brain registered so far has
  `px_atlas` 10 and `px_register` 20, so nothing changes for them.
- **Never re-align an annotated brain**: control points follow slice
  positions. The `align` mode refuses a brain that has them.
- **A clean test of the automatic annotation.** Its models (v1.0,
  `registration/auto_annotation/weights/VERSION.txt`) were trained on 24
  brains: adults and young brains of P16 to P22. None is older than P22 and
  younger than adult, so these brains are its first test on unseen ages. Keep
  what each review changed (`annotation_provenance.mat` records how many saved
  pairs are the proposal's, unchanged), and, as planned with the sandbox,
  annotate one of them by hand as well to compare the two registrations.
  Once the six are reviewed, retraining with them would give the models
  older young brains (section 10).

Into the Python route:

- Each brain already has its row in the cohort table (`common\cohort.csv`,
  one table read by both routes, decision Y4), with its age, from which its
  atlas follows; it enters the Python route when its `mapping_cohort` and
  `mapping_order` are set ([ADDING_DATA.md](ADDING_DATA.md), step 4).
- New cohorts per age (`young_P28`, `young_P32`, `young_P36`), and their
  `min_n` in `[videos]` of `mapping/settings.toml` (a cohort of one brain
  needs 1, as `young_P16` and `young_P22` have).
- `young` pools P16, P20 and P22 today (`sepmap/volumes/cohort.py`), the
  group of the grant figures. Decide before the first run with the new
  brains what `young` holds. Recommended: leave it as the P16 to P22 group,
  so the grant comparison stays reproducible, and compare the new ages
  through their own cohorts and the age analysis of section 6, never by
  pooling P16 with P36.
- Check the per-brain sheets of `run_diagnostics`, and the warp sheet
  (`04_warp_<mouse>.png`) of each new age in particular: CCF Translator has
  not carried a P28 to P36 brain in this project before.

## 6. Young against adult across ages

Two ages cannot show a window ([SCIENTIFIC_CONTEXT.md](SCIENTIFIC_CONTEXT.md),
line of work 4, "A critical period cannot be claimed yet"). With P28, P32 and
P36 added, the question becomes a trajectory:

- **What the new ages decide.** If V1's supragranular excess peaks near P28
  and falls by P36 while S1 has already declined, that is a window moving
  across modalities. The same decline everywhere at the same rate is
  maturation, which the later development of the superficial layers alone
  would also give.
- **The grant's prediction.** Associative visuo-tactile areas (RL, AL, LI)
  have a delayed or prolonged window of elevated surface GluA1 compared with
  primary sensory cortex.
- **What to build.** Per structure and per system, the per-mouse values
  against age (P16 to P36, then adult), on the readings that compare across
  sessions (`cref`, `subref`, `zref`), with the laminar contrast of step 9;
  age as a continuous variable, since there are one or two brains per age.
- **Limits that stay.** Young and adult brains were imaged in different
  sessions, and autofluorescence rises with age, so `ratio` stays a bound.
  Each age has its own atlas, and the young and adult atlases differ most
  next to white matter (`young\registration_qc\`), so results there need
  care. The
  grant's plasticity-potential index (surface GluA1 normalised for area,
  layer, age, section quality, staining batch and imaging conditions) is not
  computed: staining batch and section quality are not recorded in the
  cohort table.

## 7. A6 to A10, after P8 to P10 retire

| | addition | priority |
|---|---|---|
| A6 | within-division gene agreement, per-gene structure scatters | after the refactor |
| A7 | spatial null: surrogate maps with the same spatial smoothness, for the gene correlations and for enrichment (Fulcher 2021) | after the refactor |
| A8 | the gene panel against the autofluorescence map | after the refactor |
| A9 | gene documentation table for the supplement (Sami El-Boustani's request of 28 April) | after the refactor |
| A10 | autofluorescence and nano-minus-autofluorescence videos, threshold contours, per-gene videos for a few genes of interest | optional (S4) |

Outside both routes and still open: a DAPI control, run through the same
chain as nano (also asked for on 28 April).

## 8. Code questions held for later

### Question 22: P2 slices with too few reference pixels

`run_residual_correction` fits nano on autofluorescence per slice, over that
slice's reference pixels. With fewer than two, there is no fit: the slope is
stored as NaN, and the slice's diagnostic overlay is drawn with the previous
slice's slope (on slice 1 it would stop with an error). In production this
happened on MG904 slice 40, MG911 slice 50 and MG912 slices 7, 13 and 14
(`slice_data.slope` in each brain's
`lightsuite\correction_output\corrected_volume_slicewise.mat`; the overlays are
`diagnostic_plots\reference_pix_analysis_slice_<n>.png` beside it). MG904 and
MG911 are in the analysed cohort.

To discuss (Giulio, 3 October: a guard that rescues such a slice). Options:
fall back to the brain's mean fit (the `global` correction) or to the
neighbouring slices' fit, flag the slice in its figure and in the log, and
draw its overlay without a slope. Before choosing, check what the
`slicewise` correction gives such a slice (its scaled autofluorescence has no
fit) and what reaches its registered AUTO section.

### Output file names that record their settings

The style rule is that a file name carries every setting that changes what
the file holds ([STYLE.md](STYLE.md), Figures), and that existing names stay.
Today two run scripts break it: `run_video_compare --dlim` writes the same
name as the default (only `--vmax` adds `_vmax<V>`), and `run_closeup
--vmax`, `--dlim` and `--smooth` write the default names too (only `--cmap`
goes into a subfolder). Giulio (3 October): names should record their
settings, as long as they do not become long.

Recommended: a suffix only for a setting that differs from its default in
`settings.toml` (as `_vmax<V>` does), so a default run keeps today's names and
nothing that reads existing outputs breaks. Renaming outputs that already
exist would break their readers, and is not proposed.

### Replacing or retiring the artifact annotator

`run_annotate_artifacts` opens `ArtifactAnnotator`, a GUI to outline
artefacts slice by slice. At `align`, its mask joins the background mask of
`run_residual_correction` in the MASK channel; the NANO and AUTO channels are
not masked by it. Three brains have one (MG691, MG692, MG693:
`correction_output\artifact_mask_volume_slicewise.mat`), and no analysis reads
the MASK channel today: the plasticity chain stopped reading it with the
dead-code fix of step 8, and the Python route takes its tissue mask from the
autofluorescence. The GUI also fails when Esc is pressed again while its close
dialog is open (bug list of the plan). Giulio would rather replace it with
something better or retire it.

To decide between:

- **Retire.** `run_annotate_artifacts` leaves the run order, its code goes
  to `archive/`, and `align` writes the MASK channel from the background mask
  alone, as it does today for every brain but three.
- **Replace.** Artefacts (tears, folds, bubbles, out-of-focus tiles) flagged
  automatically per section, for example where the autofluorescence or DAPI
  departs from the neighbouring sections, reviewed in a montage like the
  slice-order editor's, and written as a mask the Python route reads beside
  its tissue mask. Worth it only if an artefact is shown to move a result.

### The barrel-field variant of the grouped figure

On 27 September Sami El-Boustani asked whether the somatosensory effect
survives when S1 means the barrel field alone. The answer,
`comparisons_v2\young_vs_adult\barrel_only\`, was written by a script outside
this repository (`D:\sep_histology\sandbox_barrel\group_plot_barrel_only.py`),
which imports `v2_region_groups` from the code root: it stops working with the
merge. Recommended: an option of `run_region_groups` that narrows primary
somatosensory cortex to SSp-bfd and writes into `barrel_only\`, as the
script did.

### Smaller items

- `write_lr_video` takes its region boundaries from the medio-lateral
  gradient only (a design choice to confirm).
- `select_background_pixels`' first knee search window differs from its
  fallbacks' (perhaps intended).
- Whether "same render" in the output comparison should count pixels or
  colour channels.
- Three small MATLAB findings of the style pass, none changing an output:
  `load_allen_atlas` and `group_differences` clear variables that do not
  exist, an `exist` check in `plot_individual_slabs` is always true, and
  `group_t_and_surprise` computes t maps that nothing reads.
- The plasticity comparison has two lists of highlighted regions: six in the
  regional surprise bars, RL included, five in the coarse t-score bars (off by
  default), without RL. One list would do.
- No check drives the engine's `sections` mode (the `U` key).
- Three places in the Python route still count a missing value as zero, or
  mix counts, found with fix 25 and left as they are because a fix would move
  values where data exist: `per_mouse.block_means` averages unsectioned
  10 um voxels in as zero, and those means also set the backgrounds;
  `to_ccf` warps a young brain linearly over the zeros outside its tissue,
  which darkens tissue edges; `compare` and `video_compare` apply the
  brain-count minimum with `cref_n` for every reading (no voxel is affected
  today).
- zref's median and spread count "brain, unassigned" as one structure among
  several hundred (`mouse_scalars`, `region_plot`); A1's declared set settles
  it.

## 9. LightSuite: upstream contributions (Z2), then the swap (Z1)

Our copy of LightSuite is upstream commit `2f16206` (29 October 2025) plus
our changes, each described in
[third_party/LightSuite/PATCHES.md](../third_party/LightSuite/PATCHES.md): two
registration fixes still missing upstream (LS1, LS2; LS2 is required for the
20 um DeMBA atlases), GUI fixes and review aids (LS3 to LS6), and the plugin
hook through which our automatic annotation enters the GUI (LS7).

- **Z2, first.** LS1 to LS7 offered upstream as separate pull requests,
  after rebasing on upstream (its GUI now has its own numbering and binds
  backspace to "delete last point"; its slice-order editor was rewritten).
  Also worth asking for: saving the per-slice transforms, an output-folder
  option, relative elastix paths, a seed in the align step, and the
  registration parameters (extent factor, CPD iterations, B-spline bins and
  samples) as options with upstream's defaults, which the swap needs. The
  mechanics: a fork holding the series, vendored with `git subtree`. Sami
  El-Boustani's parked wishes for the control-point GUI (contrast per channel,
  a way past a failed in-plane pre-alignment, moving an annotated slice's
  plane with the wheel) belong to the same contribution; confirm them with
  him before building.
- **Z2, then Z1.** Upstream was 168 commits ahead on 30 September. Its
  changes alter results: B-spline settings, the align extent (plane indices
  shift by 30, so re-aligning an annotated brain invalidates its points), CPD
  iterations, and an extraction that no longer crops and edits the curated
  decisions file. So: never rerun extraction or align on an existing brain;
  measure the noise floor of the current code first; register the reference
  mice with upstream from their existing align outputs and control points
  and compare voxel by voxel, per region and on the headline numbers;
  revalidate the automatic annotation, whose models learned the current
  geometry; record which version processed which cohort. PATCHES.md has the
  full list.
- Upstream's README announces a Python version (PyLightSuite), not public on
  30 September. Worth watching before investing in Z1.

## 10. The automatic annotation's validation sandbox

The engine is in this repository (`registration/auto_annotation/`, merged on
29 September). The work behind it, how the models were built and how they
were validated, is in a separate repository, `SEPGluA1_autoannotation`
(`D:\sep_histology\sandbox_autoannotate`): the training and validation code,
the lab notebook `LOG.md`, and `AUDIT.md`, which maps every claim to its
evidence. By decision it stays separate until after the refactor; its scripts
and `AUDIT.md` point at the tag `refactor-start` of this repository.

Reasons to merge it: retraining with the new brains (section 5) needs the
training code next to the engine it feeds; and a reader of the engine's
README should reach the evidence without a second repository. What a merge
would take:

- code and notes only (the dataset, caches, models and galleries stay on
  disk, as the sandbox's own ignore rules already keep them out of git), in a
  folder beside the engine;
- its scripts importing the engine's `core.py` instead of their own copies of
  the algorithms, where they still differ;
- its CUDA environment frozen as a `tools\requirements_<env>.txt`, like the
  others;
- `AUDIT.md`'s paths rewritten to the merged layout.

## 11. Comparable normalisation of the two routes

The two routes share the preprocessing: `run_nano_equalisation` evens the
nano channel from slice to slice within each brain, `run_residual_correction`
scales the autofluorescence onto the nano slice by slice, and both routes read
the registered volumes. After that they part:

| | plasticity comparison (MATLAB) | Python route |
|---|---|---|
| tissue | background masks from the nano channel, percentile bounds that change along AP | autofluorescence above its off-tissue median, nano non-zero, inside the atlas brain |
| background | absorbed by the intercept of the per-mouse fit | each channel minus its off-tissue median |
| per brain | a robust line (slope and intercept) fitted onto the group's median cortex | one number per brain (`cref`: the isocortex mean; `zref`: median and spread over structures) |
| between groups | the experimental group's profile aligned onto the control group's, then one shared scale | none: no brain is fitted to another |
| smoothing | 3D Gaussian, sigma 5 voxels, before the comparison | none on the tables; light, declared smoothing on maps only |
| hemispheres | left plus right and left minus right per mouse | left and mirrored right averaged |

Both are kept as they are for now (decision of 30 September). They answer
different questions, and an intercept per mouse does not preserve ratios
between regions, so numbers from the two cannot be put side by side. If
analyses of both routes go into one paper, either:

- **run the plasticity question through the Python readings**: the naive and
  RWS cohorts are already there (the column `naive_minus_rws_log2` of
  `region_stats.csv`), the behaviour mice would need `run_per_mouse` and
  `run_to_ccf`, and a voxelwise group test with the hemispheres kept apart
  would be added; or
- **keep both and say why**, reporting the Python readings as a check of the
  MATLAB result.

The check costs little and is informative already. In the per-structure table,
which takes both hemispheres together, the sign of the barrel field's naive
minus RWS difference depends on the reading: RWS is higher under `cref`
(-0.0749 log2) and `zref` (-0.0032), naive under `ratio` (+0.2356), `subref`
(+0.0546 with fix 1's reference, in the rerun on a copy; +0.0722 before it)
and `sepratio` (+0.0620), and no test is computed
(`comparisons_v2\young_vs_adult\region_stats.csv`, SSp-bfd,
`naive_minus_rws_log2`). That is no consistent whole-structure shift. The
MATLAB result is voxelwise, in one slab, so the two do not contradict each
other, and the side of the stimulation decides which hemispheric reading fits
(section 12).

## 12. Open scientific questions

**Plasticity localisation (paused).** The S1 increase after RWS is in the
approved figures, and the informative S6 run reproduced the computation on the
same inputs: slab t maps and surprise masks at 0.997 to 0.998 for RWS (0.988
to 0.999 for behaviour), individual maps at 1.000000, regional bars at 0.994 to
0.998 ([REFACTOR_PLAN.md](REFACTOR_PLAN.md), Progress, 1 October). That checks
the computation, not the effect. With the smoothing fixed (fix 23, section 1)
the increase is still there, smaller in the sum map; per mouse the barrel
field is 26% higher after RWS, p 0.12 with five against five. What would make
it convincing:

- more animals per group, the number fixed in advance;
- the protocol written into the cohort table: the RWS mice had one RWS
  session and were then perfused; the behaviour mice were perfused on their
  first day above chance in a single-whisker detection task (Giulio,
  2 October); which whiskers, which side, and the delay before perfusion are
  not recorded in the repository, and the side decides whether the
  hemispheric sum or difference is the right reading;
- a correction beyond the voxelwise p < 0.01 (uncorrected, on 3D-smoothed
  maps, five mice per group): a cluster-level permutation, or a test of
  declared structures (S1 barrel field, VPM, the posterior thalamic complex);
- the comparison repeated in new animals;
- what the naive-against-behaviour comparison shows, written down. The
  sources record only that it reproduced; with fix 23, which leaves out
  MG709 (no tissue at the slab), it shows a broad increase of the hemisphere
  sum over about a quarter of the slab, the barrel field included
  (section 1).

**What the label reports.** The tissue was not permeabilised (no or very
little detergent; Giulio, 2 October), so the nanobody reaches the surface
receptors, as the grant describes. Still missing: a specificity control (a
knockout or no-primary stain), a total-receptor measure on a subset of the
same brains (a total-GluA1 antibody, or autoradiography), since the green
SEP channel turned out to be mostly autofluorescence and cannot serve, and the
DAPI control. The beyond-abundance leftover is consistent with a surface
signal and does not show it: translation, turnover, subunit composition or
nanobody access would land there too, and only a total-receptor channel
separates them.

**The ISH comparison.** Provisional until A1 to A3. Its p values are
anticonservative (genes in a set are co-expressed); the spatial null (A7) is
the fix. Gene sets fixed before looking (SynGO, Koopmans 2019) and a
background panel of a few thousand Allen genes would turn the ranking into a
test (`adult_ish_design.md`, "Later").

**Young against adult.** A critical period cannot be claimed from P16 to P22
against adults (section 6). Also open:

- MG904 is two to four times brighter in nano than the other young brains,
  an outlier in both channel ratios and in neither internal reading
  (`comparisons_v2\README.md`); quote `ratio` and `sepratio` young means only
  after looking at the per-mouse dots.
- The grant's Figure 3b caption gives the pups as "aged between P16 and P20";
  MG904 is P22 (a mouse's age is the one in its raw folder's name) and was
  grouped with the P20 brains in the figure, as Sami El-Boustani asked. Any
  later text says P16 to P22.
- The argument that S1's lead speaks against a critical period rests on S1's
  layer 4 window closing early; no paper in `data\ref_papers\` covers it, so
  it needs a reference before it is used.

**The adults' age.** Not recorded (NaN in the cohort registry); well above P60
(Giulio, 2 October). Fill it in when found: the Allen ISH mice are 56 days
old, and the age series of section 6 ends at "adult".

---

## Decisions taken

For the reasons, see the plan and the documents named.

- **Atlases** (2 September): young brains to the DeMBA atlas of their own age,
  adults to the CCF; region statistics on each brain's own atlas, maps carried
  to the CCF with CCF Translator.
- **Routes** (29 and 30 September): Python is the main analysis route (adult
  map, ISH, young against adult); MATLAB keeps preprocessing, registration and
  the plasticity comparison, with results identical.
- **The plasticity reference** (S6): the outputs approved on 5 December 2025,
  with their mice; the line is paused, not closed.
- **The SEP channel** (26 September): mostly autofluorescence in these
  sections; `sepratio` is not a surface fraction.
- **The beyond-abundance result** (26 September): about 39% of the map's
  explainable variance is not explained by receptor abundance or synaptic
  density, and that leftover replicates across independent halves of the
  cohort at 0.934 (`adult_v2\beyond\for_sami\numbers_for_the_caption.txt`,
  `adult_v2\beyond\variance_partition.csv`, `controls.csv`). It uses all ten
  adults and the grey-matter rule already, so A1 to A3 do not move it.
- **The ISH tests**: the ranking is descriptive; the powered panel test is
  negative with a working positive control; numbers held back under S5.
- **Normalisation**: both routes kept as they are (section 11).
- **Missing values** (3 October): a mean leaves a missing value out instead
  of counting it as zero, wherever the code allows (fixes 3, 23 and 25 of
  step 8), even where that changes the approved December 2025 comparison.
- **The automatic annotation** stays on the GPU (3 October): speed matters
  more than proposals that repeat to the last pixel; the proposals are
  reviewed by hand anyway.
- **MG904** is P22, grouped with the P20 brains in the grant figure.
- **Code**: retired code to `archive/` until checked (L2); landmark_refine
  retired, the cheap point carry-over kept (L3); no licence file until
  publication (L6); one cohort table for both languages (Y4).

## Cross-references

- [REFACTOR_PLAN.md](REFACTOR_PLAN.md): the refactor, its decisions, bug list
  and progress.
- [REFACTOR_COVERAGE.md](REFACTOR_COVERAGE.md): P8 to P10 item by item, and
  the specification of A1 to A10.
- [SCIENTIFIC_CONTEXT.md](SCIENTIFIC_CONTEXT.md): the questions, the results,
  the grant, the papers.
- [FIGURES.md](FIGURES.md): which script makes each figure.
- [ADDING_DATA.md](ADDING_DATA.md): a new brain, end to end.
- [adult_ish_design.md](adult_ish_design.md): the adult and ISH analyses, as
  designed and run in September 2026.
- [third_party/LightSuite/PATCHES.md](../third_party/LightSuite/PATCHES.md):
  our LightSuite changes, Z1 and Z2.
- `comparisons_v2\README.md` under the data root: the Python route's
  readings, checks and traps when reading numbers off its figures.
