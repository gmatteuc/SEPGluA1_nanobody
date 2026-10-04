# Repository refactor plan

Draft 6, 30 Sep 2026, agreed with Giulio (checks made proportionate). Step 0 committed on 30 Sep (6a29c37).

The project is to be handed over, and this refactor consolidates it: code
another scientist can pick up and run, documents that say what it is for, and
an analysis route whose results can be defended. The grant work is done (last
figures: the cortical flatmaps), and the automatic annotation is merged and
has passed its first new brain (MG914, 30 Sep).

Draft 1 was an inventory and a layout. A read-only study of the code (30 Sep),
and an independent review of the resulting draft, changed several things:

- the Python route does not yet cover everything P8, P9 and P10 answer: five
  additions are required before they retire (item by item, with the
  specification of each addition, in [REFACTOR_COVERAGE.md](REFACTOR_COVERAGE.md));
- some quoted ISH results do not hold as quoted (see S5);
- the first verification design could have written into the production data
  tree and had no clean baseline; this draft replaces it;
- our LightSuite changes are not GUI-only: two are registration bug fixes;
- the Python half of the style guide has a full draft
  ([STYLE_python_draft.md](STYLE_python_draft.md)), derived from the Genedata
  and MaxWell projects;
- the documents follow the two-photon imaging repository (the dendrites
  repository, `D:\dendrites\code`): README with the science up front,
  scientific context, adding data, roadmap, a README per folder.

## The project in four lines of work

The documents (step 10) are written around these, and they decide what is
important in the code:

1. **Plasticity localisation** (the original aim). Surface GluA1 after an
   experience (rhythmic whisker stimulation, RWS, as in Gambino et al. 2014; or behavior) against naive
   mice, to find where plasticity lands in cortex or elsewhere. There were
   effects, but not robust enough to be sure, and more animals would be a
   large investment, so the line is paused, not closed. MATLAB: collect by
   group, normalise, group differences (today P5, P6bis, P7bis). It stays
   runnable, with identical results, and documented well enough to be resumed
   with more animals.
2. **Adult distribution.** How the nanobody signal (in principle surface AMPA
   receptors containing GluA1) is distributed across the whole brain, and how
   reproducible that is across mice. The adult cohort is the 10 naive and RWS
   mice pooled (checked: the part of the map that abundance and density do
   not explain agrees between the five naive and the five RWS mice, rho
   +0.88); behavior adults are not characterised here, only in the plasticity
   comparison. Was P8, and P10 for the autofluorescence control. Becomes
   Python.
3. **What the map measures.** Against the Allen in situ hybridisation maps:
   the map is not explained by receptor abundance or synaptic density alone,
   which points to surface GluA1. Was P9; the Python route extends it (390-gene
   panel, gene reliability, the beyond-abundance regression and its controls).
4. **Young against adult** (grant, September 2026). Young mice (P16 to P36)
   against adults, looking for differences in cortex that could correspond to
   critical periods. Python.

Preprocessing and registration (MATLAB, with LightSuite and the automatic
annotation) serve all four.

## Rules

- Results stay identical unless a change is intended, listed and checked. A
  bug found on the way is listed here and fixed in the Fixes step (8), one
  commit each, with a before/after check.
- Pure moves (`git mv`) are committed separately from content changes, so git
  keeps each file's history. A script that becomes a function is moved with
  `git mv` and a new short driver is written beside it. No case-only renames
  (`core.ignorecase` is true). The style-pass commits go into
  `.git-blame-ignore-revs`.
- Renaming follows an explicit old-to-new table, never a pattern on P plus
  digits: ages (P20), atlas keys (`demba_p20`), cohort tags (`young_P20`),
  data file names (`nano_4d_P20.mat`) and P7bis's P99 metric keep theirs. The
  messages that name a script (about 60) are updated in the same commit.
- The work happens on a branch whose code folder sits inside its own check
  tree on G: (see Verification design), so `main` stays usable. During the
  refactor `main` gets only uncommitted run-setting edits; any fix needed on
  `main` is logged and ported to the branch by hand.
- **No check ever writes into the production data tree.** Checks run on copies
  (check trees), in fresh MATLAB sessions, and print the code root, the data
  root, the commit and `which` of the path functions they resolved.
- Raw `.czi` files, on `S:` and under `data\<group>\<mouse>\` on D:, are never
  written, moved or deleted. Check trees hold copies. Caches are deleted only
  inside check trees. The G: snapshot is never used as a check tree.
- Commits are authored by Giulio alone, and code and documents credit no
  other author.
- LightSuite stays third-party code. Only the files we already changed are
  changed, and every change is listed in `third_party/LightSuite/PATCHES.md`.

## Decisions already taken

- Runnable scripts are named `run_...` with a short header explaining the
  flow; the P-numbers go. Code is grouped into folders by pipeline, drivers at
  the top of each, the functions they call in a `pipeline/` subfolder.
- Python is the main analysis route: adult distribution, ISH, young against
  adult. MATLAB stays for preprocessing, registration, and the plasticity
  comparison (naive against RWS, naive against behavior).
- A MATLAB analysis retires only when the Python route answers the same
  scientific question, not necessarily with the same output. If it does not,
  the question is added to Python first.
- Dropped on 29 Sep, not rebuilt in Python: P8's per-mouse videos, P9's
  voxel-level ISH correlation, and P8's distance-weight diagnostic figure. The
  border guard that figure justified is not dropped (A3, A4).
- Giulio checks the Python route against the old code after the refactor.
  That check runs the old code from the tag `refactor-start` in its own check
  tree, never the 29 Sep snapshot's code as it is: its Python writes to
  `D:\sep_histology\data` and its sandbox scripts import
  `D:\sep_histology\code`.
- The Python style follows the Genedata and MaxWell projects. READMEs have no
  section on how the code was written.
- The validation sandbox (`SEPGluA1_autoannotation`) stays a separate
  repository for now; merging it in is a possible later step, after this
  refactor.
- LightSuite: our current copy stays. The patches are restructured during the
  refactor (steps 4 and 6). Later, as a separate track, they are offered
  upstream, and the copy is replaced by the latest upstream plus our patches,
  with a rerun check.
- Order: the refactor comes before the remaining young mice (MG896, MG906,
  MG895, then MG907 and MG908), which then go through the refactored code.

Decided on 30 Sep (Giulio):

- **L1.** The Python route is the folder `mapping/` (entry points, settings,
  tests) with the package `mapping/sepmap/` inside it.
- **L2.** Retired code goes to `archive/` until Giulio's check, then is
  deleted (the tag keeps it).
- **L3.** landmark_refine (the image matching behind the GUI's r key) is
  retired. The cheap point carry-over stays (t, take the neighbour's points;
  p, carry forward as provisional), so a brain can still be annotated on a
  computer without a GPU. The retirement is one change in step 6: P4
  `annotate` stops starting the worker (its stale "(r, t)" comment goes), the
  files go to `archive/`, the r key's code leaves LightSuite's GUI, and
  PATCHES.md records it.
- **L4.** P2's `global` correction variant is archived; `slicewise` stays.
- **L5.** The reference adult is CGF027_Gria1 (naive, CCF, hand-annotated,
  with its SEP channel); with the lighter checks it is part of the small
  reference set, and the registration check uses MG914 alone.
- **L6.** No licence file for now. (The grant commits to releasing the analysis
  code on publication: the licence, and the GPL obligations of our modified
  LightSuite, come back then.)
- **L7.** The literature goes into `docs/SCIENTIFIC_CONTEXT.md` as in the imaging
  repository: grouped by the question each paper bears on, what it shows and
  what it means here. The PDFs stay in `data\ref_papers\`.
- **Y1.** No comment above a function repeating its docstring.
- **Y2.** Double quotes everywhere, existing code included, for consistency.
- **Y4.** One cohort table, read by MATLAB and Python alike. Done in step 8,
  as a fix with a before/after check, not in step 5: the Python route stacks
  its brains in its own order (young P20, P16, P22, naive, RWS) and the
  MATLAB registry in the legacy order, so a shared table either keeps a
  per-route order or changes cohort means in their last bits.
- **Y5.** The two Python verification tools are restyled.
- **Y6.** The MATLAB half of the style guide is the imaging repository's
  `docs/STYLE.md`, adapted (paths through `get_paths` and `SEP_DATA_ROOT`, the
  cohort registry, the hot colormap, grey reliability bars).
- **Y7.** The style pass works from the guide's quotes only.
- **Z1, Z2** (the LightSuite swap and the upstream contributions): later.

## Decisions on the science (30 Sep, Giulio)

- **S6. The plasticity comparison as Sami approved it.** The reference is the
  output of 5 Dec 2025, approved by Sami: `naive_vs_rws` and
  `naive_vs_behavior` in `S:\ElboustaniLab\#SHARE\Processed\SEP-GluA1_Project\comparisons\comparisons\`
  (the copies in `data\comparisons\` are identical, file by file). Its headline
  is the increase of the nanobody signal in S1 after RWS (hemisphere-sum t map,
  slab 565, surprise-masked at p < 0.01). Its mice: naive CGF027, CGF028,
  CGF033, CGF034, CGF035; RWS MG691, MG692, MG693, MG736, MG737; behavior
  MG705, MG709, MG716, MG718 (of the 7 in the registry; MG709 has no tissue at
  slab 565, which is probably why P7bis's comment speaks of three mice). This
  is exactly what the code selects today, so the settings stay as they are,
  and the comment is corrected (Fixes). Whether today's code still shows the
  S1 increase is checked once, for information: today's P7bis on the
  preserved inputs of that run (the normalised volumes of 26 and 27 Nov 2025,
  `nano_4d_normalized_bk.mat`), the result noted in `docs/ROADMAP.md`. It is
  not a gate: the old code changed a lot, and reproducing old outputs
  exactly is not the goal.
- **S1 to S5** as recommended below, with two amendments: Benjamini-Hochberg
  is reported beside the uncorrected p values, and if it leaves nothing
  standing the correction is reconsidered together (for example a declared set
  of structures of interest from the grant: V1, S1, RL, AL, LI and the control
  regions; tests at division level; the spatial null) rather than applied
  blindly; the per-gene videos stay optional, drawn for a few genes of
  interest (A10).
- **A6 to A9** wait until after the refactor.
- **Normalisation across the routes.** The Python route normalises each brain
  on its own (background subtraction, then one scale per brain: cref, zref);
  the plasticity comparison fits each brain to its group's median and the
  experimental group onto the control group. Both are kept as they are for
  now. If analyses of both routes end up in the same paper, their
  preprocessing and normalisation must be made comparable, or the difference
  discussed: this is an open item in `docs/ROADMAP.md`.

### The recommendations S1 to S5, as agreed

- **S1. The declared structure set.** Which structures enter the adult
  distribution, the ISH gene ranking, the nano-against-autofluorescence
  summaries and the zref reference. Today the gene ranking averages each
  structure over however many adults have it: 22 structures of the adult table
  are seen in fewer than 5 of the 10 adults (15 of them in the gene ranking,
  mostly pons, medulla and midbrain), and their values do not replicate across
  mice. zref is referenced to the structures shared by all 17 brains, young
  included, so adding young brains moves every adult's zref. **(rec.)**
  grey-matter structures (the rule the beyond-abundance analysis already uses)
  seen in all 10 adults, written once to a table that every analysis reads,
  with every dropped structure and the reason (number of mice, white matter,
  catch-all label). The all-10 rule keeps 234 of 280 structures; some drop
  only because of the 250-voxel minimum per mouse. P9's nine-division set is
  reported alongside.
- **S2. What "enriched" means.** P8 called a region enriched at LR-sum >= 1.5
  or z >= 0 (zero at the voxel mean); Sami asked (28 Apr) for a threshold of 2
  on the LR-sum and, in the same feedback, for a permutation null. **(rec.)**
  per structure, an exact signed-rank test of zref against the brain's median
  structure across mice, Benjamini-Hochberg across structures, with the spatial
  null (A7) as the stronger version. Zero now means the median structure, not
  the voxel mean, so "enriched" changes meaning: to agree with Sami. At n = 10
  the smallest two-sided exact p is 1/512.
- **S3. The autofluorescence control per structure (P10's question).**
  **(rec.)** two-sided tests on centred contrasts. Inputs are the per-mouse
  log2 nano and log2 autofluorescence structure means, never the stored
  `ratio` reading, whose level depends on exposure. Each contrast is shown
  beside the nano and autofluorescence values (6 of P10's 21 delta-z hits were
  low-autofluorescence periventricular structures with nano below the
  median). BH within each contrast family at structure level (Bonferroni
  cannot reject at n = 10); Bonferroni over the 9 divisions. This keeps P10's
  question, not its lists: P10 was one-sided, on a scale dominated by the
  hippocampus, with inputs from three dates.
- **S4. The per-gene delta-z video and P8's threshold videos:** optional or
  dropped? And does the 29 Sep drop of the "distance-weight diagnostic" cover
  P9's eroded-against-distance-weighted ranking comparison? **(rec.)** drop
  the per-gene video; the threshold contours survive only inside the
  optional A10; P9's comparison is covered by A3's robustness check.
- **S5. Which ISH numbers to quote.** None of the following until A1 to A3
  have run: Gria1's rank (22 in old P9; 10 in the Python route as is; 17 once
  each structure must be seen in at least 5 adults); the powered panel test
  (p = 0.74) and the roles permutation, which run on the same unrestricted
  structure set; "control genes cluster near zero", which holds for the
  control_inhib and control_glia medians but not gene by gene (Chrm1, Htr3a
  and Drd1 sit at +0.66 to +0.67). What holds already: Cacng8 first under
  every structure set tried (rho 0.77 to 0.82). The beyond-abundance result is
  not affected (it already uses all adults and the grey-matter rule).

## Target layout

```
README.md               the project, the four lines of work, setup, pipelines, outputs
ruff.toml               Python lint and format settings
sep_setup_paths.m       adds an explicit list of code folders, LightSuite once and
                        below ours; never genpath of the root, never data or atlas
                        folders; refuses to run when a script the branch moved is
                        back at the root under its old name (an editor tab saved
                        after the merge would do that). A minimal version, which
                        only adds the code root, is committed in step 0
get_paths.m             stays exactly one level below the project root (it derives
                        the root from its own location); gains SEP_DATA_ROOT and
                        the production-data guard

common/                 cohort registry, volume reading, left/right statistics, colours
atlas/                  atlas lookup, building DeMBA atlases, atlas QC
preprocessing/          raw copy, extraction, slice order, residual correction,
                        nano equalisation, artifact annotation
registration/           registration to the atlas, the SEP channel, QC, and the
                        automatic annotation: MATLAB side, the GUI plugin
                        (annotation_gui/), the Python engine (auto_annotation/)
group_comparison/       plasticity localisation: collect by group, normalise,
                        group differences
mapping/                the Python route (L1): run_*.py in run order, settings.toml,
                        README.md, tests/, and the package:
    sepmap/             config.py (paths and cohort; imports no analysis code),
                        plotting.py, volumes/, young_vs_adult/, adult/, ish/,
                        structures.py (the declared set, added in step 9)
tests/                  MATLAB tests, sep_run_all_tests
tools/                  verification tools (MATLAB ones with a sep_ prefix),
                        detached-run helper, requirements_<env>.txt
third_party/            README.md (licences, do-not-run list), LightSuite (+ UPSTREAM,
                        PATCHES.md), matlab_elastix, yamlmatlab, BioformatsImage
archive/                retired code until checked (L2)
docs/                   STYLE.md, scientific context (with the literature), adding
                        data, roadmap, figures, production settings, plans
assets/                 images for the README
```

Names that also exist in the dendrites repository (`setup_paths`,
`project_paths`, `compare_outputs`, `check_code_identity`, `run_driver_copy`,
...) get a `sep_` prefix here: both projects run in the same MATLAB, and a
shadowed `run_driver_copy` would redirect through the other project's variable
and write into real data. A test (step 4) lists every function on the project
path, third_party included, and the dendrites repository's function names when
its checkout is present; it fails on duplicate names (case-insensitive) and
runs `which -all` on each name to catch shadowed MATLAB and toolbox functions.

Output folders and file names under `data\` do not change, with two
exceptions: the detached runner names its logs after the script, and
`v2_diagnostics` regenerates `processing_diagnostics\README.md` with the script
names in it (an expected difference in the output comparison). Data files are
not edited; a table of old and new script names, with the tag
`refactor-start`, goes into `docs/ADDING_DATA.md`.

## Inventory

**live**: part of a pipeline in use; **tool**: run by hand when needed;
**QC**: checks, kept; **retired**: to `archive/`.

### Preprocessing (MATLAB)

| now | status | becomes | notes |
|---|---|---|---|
| P0_copy_raw_data.m | live | preprocessing/run_copy_raw_data.m | reads S:, writes D: only |
| P1_extract_and_center_data.m | live | preprocessing/run_extract_and_center.m | |
| P1bis_order_slices.m | live, manual | preprocessing/run_order_slices.m | opens SliceOrderEditor |
| P2_residual_correction_analysis.m | live | preprocessing/run_residual_correction.m | 498 lines, split |
| P2bis_nano_correction_analysis.m | live | preprocessing/run_nano_equalisation.m | 543 lines; P4 reads its equalized_volume.mat (use_equalized_nano = 1) |
| P3_annotate_artifacts.m | live, manual | preprocessing/run_annotate_artifacts.m | |
| ArtifactAnnotator.m, select_reference_pixels.m | live | preprocessing/pipeline/ | |
| select_background_pixels.m | live | common/ | P2bis, P6bis |
| make_ordering_volume.m, make_atlas_reference_sheet.m | tool | preprocessing/ | |
| BioformatsImage/explore_czi_G.m | tool (ours, inside a vendored folder) | preprocessing/ | dumps .czi channel metadata; its path to S: becomes a setting |
| BioformatsImage/extractAxioscanImages.m | not called | archive/ | writes folders next to the .czi files it reads: never to be run on raw folders |

### Registration

| now | status | becomes | notes |
|---|---|---|---|
| P4_register_to_atlas.m | live | registration/run_register_to_atlas.m | five modes stay; the settings block stays one small block at the top (edited daily) |
| P4bis_add_sep_channel.m | live | registration/run_add_sep_channel.m | |
| auto_annotate.m | live | registration/pipeline/ | finds the engine through get_paths, not its own location |
| auto_annotation/, setup_auto_annotation.ps1 | live | registration/auto_annotation/ (setup.ps1 inside) | the setup script's paths are fixed in the same commit |
| the automatic-annotation part of LightSuite's GUI | live | registration/annotation_gui/ | moved out as a plugin in step 6 (see LightSuite) |
| recompute_backvalues.m | live | registration/pipeline/ | |
| remap_control_points.m | tool | registration/ | carries points across a reorder |
| check_registration_error.m, verify_demba_setup.m | QC | registration/qc/ | verify_demba_setup checks behaviour, not P4's source text (it fails today) |
| landmark_refine*, setup_landmark_refine.ps1 | L3 | archive/ in step 6 | stays in place until then, so P4 annotate keeps working |

### Atlas

| now | status | becomes |
|---|---|---|
| get_atlas.m, get_atlas_crop.m, cohort_atlas_key.m, get_allen_region_mask.m | live | atlas/ |
| build_demba_atlas.py | tool | atlas/ |
| atlas_diagnostics.m, check_demba_to_allen.m, compare_atlas_regions.m, compare_atlases_montage.m | QC | atlas/qc/ |

### Shared

| now | status | becomes |
|---|---|---|
| get_paths.m | live | stays at the root |
| get_cohort.m, get_cohort_spec.m, loadVolume.m, compute_lr_stats.m, get_color2color_colormap.m | live | common/ |
| test_backward_compat.m | QC | tests/ (its hard-coded paths and its addpath of `main`'s code are fixed in step 4, before it is ever run from the branch) |
| docs/adult_ish_design.md | doc | docs/, folded into SCIENTIFIC_CONTEXT and ROADMAP |

### Plasticity localisation (MATLAB, kept)

| now | becomes | notes |
|---|---|---|
| P5_collect_data_by_group.m | group_comparison/run_collect_by_group.m | |
| P6bis_analyze_group_averages_and_normalize.m | group_comparison/run_normalise_groups.m | 841 lines, split |
| P7bis_analyze_group_differences.m | group_comparison/run_group_differences.m | 1699 lines, split; the comparison (naive against rws or behavior) is one setting |
| write_lr_video*.m, write_lr_indiv_*.m (6 files) | group_comparison/pipeline/ | |

### Python route

In step 4 the `v2_*.py` files move unchanged into `mapping/`, flat, so their
imports keep working and each still runs as `python v2_x.py`. In step 5 they
become the package: flat imports become package imports; the 11 modules that
define the data root (the other 15 import it from `v2_per_mouse`) read it from
`config.py`; every module loses its `__main__` block and is called from a
`run_*.py` next to the package. That is a code change, checked by rerunning.

| now | becomes |
|---|---|
| v2_per_mouse, v2_to_ccf, v2_cohort | sepmap/volumes/ |
| v2_compare, v2_replot, v2_region_plot, v2_region_groups, v2_video, v2_video_compare | sepmap/young_vs_adult/ |
| v2_inspect | sepmap/young_vs_adult/closeup.py (not `inspect`, which shadows the standard library). It runs in the flatmap environment, so it imports only `config.py`; the constants it copies today (SIGNED_READINGS, LOG2_FLOOR, MIN_N) move into `settings.toml` |
| v2_beyond_density, v2_beyond_controls, v2_beyond_regression, v2_beyond_figures, v2_sep_channel_check, v2_adult_arms | sepmap/adult/ |
| v2_panel_build, v2_panel_fetch, v2_ish_regions, v2_ish_reliability, v2_ish_compare, v2_ish_words, v2_ish_roles, v2_ish_panel_test, v2_ish_arms | sepmap/ish/ |
| v2_diagnostics | sepmap/diagnostics.py |

The ISH chain runs twice today (the 100-gene panel, then the 390-gene panel),
and only environment variables say so (`V2_ISH_PANEL`, `V2_ISH_TABLE`; the
first is an absolute production path). Both passes become explicit steps in
the run order. `V2_READINGS` and the MATLAB overrides `P4BIS_MICE` and
`SEP_COHORT_SPECS` are kept, or fail loudly under an old name
(`SEP_MERGE_SPECS` retires with P8); every run prints the settings in force.

### Waiting, then retired

P8, P9 and P10 stay where they are, untouched (not rerun, moved or restyled),
until the additions (step 9) are done; then they go to `archive/`. Their frozen
outputs in `data\comparisons\` (and on the G: snapshot) are what the additions
are compared with. The path guard checks only for the scripts the branch has
moved.

| now | why |
|---|---|
| P8, P9, P10, plot_violinplot.m (P9 only) | once A1 to A5 are in the Python route |
| P6_analyze_group_averages_and_normalize.m, P7_analyze_group_differences.m, plot_abs_slice.m, plot_diff_slice.m, write_diff_video.m | replaced by P6bis/P7bis; to archive in step 4 |
| compare_young_vs_adult_lrsum.py, replot_young_vs_adult_lrsum.py, region_means_raw_per_mouse.py, region_ratio_young_vs_adult.py, plot_region_ratio_young_vs_adult.py | the first young-against-adult attempt; to archive in step 4 (they still write `D:\sep_histology\data` literally, with no guard: never run from a check tree) |
| lr_sum_and_diff.m, scratch.m, bk/ | not called; to archive in step 4 (the elastix recipe of bk/LightSuite.txt is copied into the README in the same step) |

### Inputs and couplings to keep

- `data\gene_targets.csv` (the hand-written 100-gene panel) is not in the
  repository, and `*.csv` is ignored: it becomes a versioned, documented input.
- `v2_ish_compare` reads P9's `gene_panel_summary.csv`
  (`data\comparisons\merged_naive_rws_vs_ish_summary_nosmooth\`, from an older
  P9 whose output folder had no channel tag) for the old-against-new ranking
  check. It stays a frozen file, copied into every check tree. Today
  `v2_ish_compare` skips that check silently when the file is missing: it
  should fail instead (Fixes). Olig2 and Calb2 are left out when the agreement
  is quoted (their old grids were stretched off the reference).
- `v2_diagnostics` sheet 08 reads P6bis's background masks, so it runs after
  P6bis in every run order.
- The GUI's `*` outlier mark repeats P4's `report_suspect_pairs` thresholds
  (20 px, 30 px); after the plugin move, both read one definition.
- The sandbox repository and its audit document call the main code by path:
  they are pointed at the tag `refactor-start`. `sandbox_barrel`, the MATLAB
  bookmarks and the project's working notes are updated to the new names after
  the merge.

## Python route: what it must gain before P8, P9 and P10 retire

Item by item, with the specification of each addition, in
[REFACTOR_COVERAGE.md](REFACTOR_COVERAGE.md); where a specification there
differs from this plan, the plan applies. Of 88 outputs and tests of P8, P9
and P10, 28 are covered (usually more rigorously: explicit background, a tissue
mask from autofluorescence, half-cohort reproducibility, the 390-gene panel,
gene reliability, the powered and controlled tests), 26 partly, 19 not at all;
4 were dropped by decision and 11 are file handling. One missing item is the
dependence of the gene ranking on the structure set, the reason for S1 and S5.

| | addition | replaces | priority |
|---|---|---|---|
| A1 | the declared structure set and zref reference (S1): one table, used by tables and maps alike | P8's z-scoring and region set | required, first |
| A2 | ISH per-section quality control: failed sections flagged and set to missing, never interpolated; a reviewed exceptions list for genuine regional absence; the effect on the ranking reported. An unreviewed first scan flags 14 of 95 genes, two of them (Tac1, Glra1) likely genuine absence; Gria1 is clean | P9's section repair | required |
| A3 | the gene ranking on A1's set and A2's tables (one shared `adult_profile` with the minimum-mice rule), rerun of the ranking, roles, panel test, arms and words; its robustness to the statistic (Spearman, Pearson) and to borders on both sides (ISH full against eroded, nano plain against eroded); a sensitivity table under three structure sets | P9's region set, metric comparison, headline ranking | required |
| A4 | adult distribution: every structure by division, per-mouse mean and SEM, reliability, the enrichment call (S2), an eroded nano mean beside the plain one (calls that change are flagged); also for the young cohort on its own | P8's bar charts, tables, enrichment | required |
| A5 | autofluorescence as a parallel control per structure and division (S3), and the autofluorescence distribution itself | all of P10, P8's `auto` run | required |
| A6 | within-division gene agreement, per-gene structure scatters | P9's per-gene figures | after the refactor |
| A7 | spatial null (surrogate maps with the same spatial smoothness) for the gene correlations and for enrichment | Sami's permutation request, never built | after the refactor |
| A8 | the gene panel against the autofluorescence map | P9's `auto` option, never run | after the refactor |
| A9 | gene documentation table for the supplement | Sami's request of 28 Apr | after the refactor |
| A10 | autofluorescence and nano-minus-autofluorescence videos, threshold contours, per-gene videos for a few genes of interest | P8/P9/P10 videos | optional (S4) |

These change numbers on purpose (A1 moves every zref, so the young-against-
adult differences behind the grant figures; A1 to A3 move the gene ranks), so
they come after the moves and the identity checks, each with its before and
after written down. `docs/FIGURES.md` records each grant figure against the
tag that made it, and the report gives the new values beside the old.

Outside this refactor, and still open: a DAPI control (in neither route).

### Known defects of the retiring route

Differences from P8 to P10 caused by these are expected, not refactor errors.
The list goes into `docs/ROADMAP.md` and the final report.

- P8's adult outputs carry an artefact of the `abs()` over slabs (P6bis was not
  rerun after the raw-zero fix).
- P9 stretched the 200 um ISH grid onto the atlas box (a 1.5 to 2.5% scale
  error, up to about 0.17 mm); its missing-data value (-1) was interpolated and
  clamped to 0; the off-reference Olig2 and Calb2 grids were stretched.
- P9's section repair "repaired" true absence of expression (Slc17a6) and
  copied the last good section over posterior planes (Cacng8, sections 48 to
  66).
- P8's and P10's autofluorescence came from a stale `auto_4d.mat` (19 Nov) that
  P5 no longer writes; P10's inputs mix caches from 2 May and 4 Sep.
- P10's Bonferroni over structures cannot reject at n = 10; P9's category
  ANOVA was anticonservative.

## Verification design

Proportionate, not exhaustive (Giulio, 30 Sep): the MATLAB code is old, it
changed a lot over time, and its runs are long. The checks prove that the
refactor changes nothing with as little running as possible, and comparisons
with historical outputs inform rather than decide.

### Two questions, kept apart

1. **Does the refactor change anything?** The new code must behave as the code
   at `refactor-start` does. Most steps (moves, comments, layout) are proven
   without running anything, by comparing the code's parse trees. Runs are
   needed only where the structure of the code changes (scripts into
   functions, the big files split, the Python package), and they use a small
   reference set.
2. **Do today's results still look like the old ones?** Informative only: one
   run of today's P7bis on the preserved inputs of the approved December 2025
   comparison (S6), noted in `docs/ROADMAP.md`. Old outputs, old names and old
   defects are not chased.

### Where checks run

- **Two small check trees on G:**, each with its own code copy inside it, next
  to its data copy, so the data root a copy derives from its own location is
  the copy's:

  ```
  G:\sep_refactor\ref\code      worktree at tag refactor-start (the old code)
  G:\sep_refactor\ref\data      the small input set and the reference outputs,
                                read-only after step 3
  G:\sep_refactor\check\code    the refactor branch
  G:\sep_refactor\check\data    inputs restored from ref before each run
  ```

  Tens of GB in all.
- **Guards**, committed on `main` in step 0 so the tag carries them:
  `get_paths.m` and the Python data root refuse the production data (and any
  folder inside it) when the code folder is not the production one, a data
  root inside the code folder, and any data root inside the G: snapshot.
  `SEP_DATA_ROOT` overrides the derived root, for runs that need it.
- **Inputs are listed, not guessed.** Each tree's input list is built from what
  the drivers read, and a missing input stops the run. Also copied: the frozen
  P9 summary, `gene_targets.csv`, and the downloaded caches (Allen, ontology,
  gene annotations), byte-identical. The network steps (`v2_panel_build`,
  `v2_panel_fetch`) run from cache only.

### The small reference set

- **Plasticity chain** (P5, P6bis, P7bis): 2 naive (CGF027, CGF028), 2 RWS
  (MG691, MG692) and 2 behavior mice (MG705, MG716); P7bis once per
  comparison; videos off (the t-score and surprise videos off together, since
  one sets a limit the other uses).
- **Python route:** the per-mouse steps (`v2_per_mouse`, `v2_to_ccf`,
  `v2_cohort` and what reads their outputs) on 2 adults and 2 young brains of
  two ages (CGF027, MG691, MG903 at P20, MG911 at P16), so the atlas handling
  of each age is exercised; the table-level analyses (region statistics, both
  ISH passes, the beyond-abundance analysis, the SEP-channel check) on the
  existing full tables, since they are fast.
- **Registration:** P4 `register` and P4bis on a copy of MG914's folder, run
  twice with the old code to measure elastix's unseeded variability; P4
  `autoannotate` on MG914 once the GPU is free. Align is never rerun (the
  step 0 guard refuses).
- **Preprocessing** (P1 to P3): the moves are proven by code identity; if P2 or
  P2bis is restructured, it is run on MG914 alone, old against new on the same
  input (P2bis's result depends on which mice run together, so it will not
  match production, and does not need to).
- The settings of each run are in `docs/production_settings.md`; settings
  without an environment variable are edited in the check copy, and the edit
  is saved as a patch beside the outputs.

### What is compared with what

- **Pass criterion.** Where the old code is deterministic, new must equal old:
  tables and `.mat` values exactly, figures up to renderer noise. Where it is
  not, new against old must differ no more than old against old. Measured in
  step 4: the registration of MG914 came out identical voxel for voxel in two
  old runs and the new one, despite elastix's unseeded sampler; the automatic
  proposal on the GPU is not reproducible (see the bug list), and two
  Python outputs depend on string hashing or file dates (see the step-4
  check below).
- **Quick and full checks.** A quick check after each commit (the plasticity
  chain on the small set without videos, one brain through the Python
  per-mouse steps); the full small-set check at the end of each step, a few
  hours.
- **Caches.** Computed caches (`*_scalars.npz`, the ISH tables
  `gene_region_table*.csv`) are deleted in every tree before each run; file
  dates are checked after it, so a check that rewrote nothing cannot pass.
- **Order of a run:** P5, P6bis, P7bis, then the Python route (sheet 08 of the
  diagnostics reads P6bis's masks).

### Sessions, tools, environments

- **Sessions and paths.** Fresh MATLAB, `restoredefaultpath`, then
  `sep_setup_paths`; the explicit folder list keeps `archive/`, `tests/`,
  `tools/venv_*`, `.git` and worktree folders off the path.
- **The detached runner** (`tools/run_matlab_detached.ps1`) requires
  `-CodeDir`, `-DataRoot` and `-LogDir` with no defaults, sets `SEP_DATA_ROOT`
  for the MATLAB it starts, refuses a log folder under the production data
  when the data root is a copy, and runs `restoredefaultpath; cd(CodeDir);
  sep_setup_paths; <script>` after printing the commit and the data root.
- **Verification tools.** The MATLAB ones come from the imaging repository's
  `tools/`, copied with a `sep_` prefix (the shadowing risk is MATLAB's); the
  two Python ones (`check_code_identity.py`, `compare_outputs.py`) come from
  the refactor template and keep their names. All are adapted: the redirect
  goes through `SEP_DATA_ROOT` (inputs and outputs), and the identity check
  takes an old-to-new name map, so a moved file cannot pass uncompared.
- **Code identity** is compared through the parse tree, never byte by byte
  (line endings differ between the working tree and a checkout). The baseline
  for the style pass is taken after the moves and splits.
- **Environments** are frozen as they are (`pip freeze` of each, pip always
  called as `python -m pip`), never recreated during the refactor, and never
  given an editable install. `venv_atlas` holds 4.7 GB of DeMBA-to-CCF
  deformation fields downloaded in September: they are copied, with a file
  list and hashes, to `data\atlas\ccf_translator_fields\` (the snapshot
  refresh takes them to G:), and the
  README says how to restore them after a reinstall. pytest and ruff go into a
  separate `tools\venv_dev`, so the analysis environments do not change. Runs
  from a check tree set `AUTO_ANNOTATION_PYTHON` to `main`'s engine
  interpreter.
- **Merging.** Giulio stops the landmark_refine worker and closes every MATLAB
  session first (an open editor tab saved after the merge would recreate the
  old P4 at the root). The uncommitted run settings are written down, stashed,
  and re-applied to `run_register_to_atlas.m`. After the merge, the ignored
  engine environment is moved by hand to
  `registration\auto_annotation\.venv` (git moves only tracked files) and its
  self-test run; the ignored leftovers at the old paths are removed with
  Giulio's agreement.

## Steps

### 0. Safety commit on `main`

1. `SEP_DATA_ROOT` and the production-data guard in `get_paths.m` and in every
   Python file that defines the data root (11 `v2_*.py` modules and
   `build_demba_atlas.py`); the Python default found by walking up from the
   file to the folder that holds `get_paths.m`, then its parent plus `data`,
   so it survives the files moving one level down in step 4;
   `V2_ISH_PANEL` relative to the data root. Checked by comparing `get_paths()`
   and `DATA` with the variable unset against their old values; no driver runs
   on production.
2. A guard in P4 `align` that refuses when control points already exist, as
   `angle` already does. Only the guard is committed; the run settings in P4
   stay uncommitted.
3. `pip freeze` of the four environments into `tools/`; the deformation fields
   copied with their hashes; `tools\venv_dev` with pytest and ruff.
4. `.gitignore`: anchor `data/`, `output/` and `comparisons/` to the root;
   re-include `third_party/**` (the vendored re-include rules stop matching
   after the move), `assets/` and `docs/` (test fixtures and the versioned CSV
   inputs get their re-include lines when they first appear); the
   local-settings line for the editor folder moves to
   `.git/info/exclude`. `.gitattributes` marks `*.pt`, `*.ckpt`, `*.nii.gz` and
   `*.jar` as binary.
5. A minimal `sep_setup_paths.m` (adds the code root, as today's drivers
   expect), so the tag carries it and the runner works on the old code; the
   detached runner and the verification tools, as above.
6. `docs/production_settings.md`, drafted for Giulio's sign-off.
7. Commit this plan, `REFACTOR_COVERAGE.md` and the style draft.

### 1. Decisions

L1 to L7, Y1 to Y7 and S1 to S6 are decided. `docs/STYLE.md` is assembled: the shared
rules and the MATLAB half (Y6), the Python half from the draft, and the
style-pass rules (how comment-only and code changes are made and checked).
The reference examples are named at the start of step 7.

### 2. Freeze

Refresh the G: snapshot (only what changed since 29 Sep is copied). Tag `main`
as `refactor-start` (`before-refactor` stays). Tag the code state behind the
grant figures as `grant-2026-09`, on the last commit before step 0 (its
Python route is unchanged since 26 Sep, before the last grant figures;
Giulio confirms). `refactor-start` differs from it only by step 0. Create the two worktrees in their check
trees.

### 3. Reference capture

Build the two small trees, run the old code on the small reference set as
described in Verification design (MG914's registration twice), keep every
output and log, record run times. Run it when the machine is free: another
heavy process would not change the results, but both would slow down.

Then the informative run of S6: today's P7bis on the preserved December 2025
inputs, compared by eye and by the `.fig` data with the approved figures; the
outcome goes into `docs/ROADMAP.md` and does not hold up step 4.

### 4. Organise

One commit of pure moves (`git mv` into the layout, the Python files flat into
`mapping/`, the files retired now into `archive/`, LightSuite into
`third_party/LightSuite/`). Then the content commit: `sep_setup_paths`, the
drivers' own `addpath(genpath(...))` calls removed, `get_paths`' LightSuite
path, the rewritten `.gitignore`, the engine lookup and setup paths, the path
test with duplicate names, `test_backward_compat.m`'s paths, and the elastix
recipe of `bk/LightSuite.txt` copied into the README before `bk/` is
archived. LightSuite gets
its `UPSTREAM` file (repository, commit `2f16206`, 29 Oct 2025), `PATCHES.md`
describing today's differences, a one-line modification notice in each
modified file (GPL-3.0 section 5a; a comment, so code-identical), and
`third_party/README.md` with the licences and the do-not-run list. Full check.

### 5. Scripts into functions, one pipeline at a time

Common and atlas, preprocessing, registration, plasticity comparison, then the
Python route (first sub-step: the package, its imports and `config.py`). Each
driver becomes header, settings and a few calls; the body moves unchanged into
`pipeline/` or the package. `checkcode` on every new function (a script turned
function can silently pick up a built-in, as P4 once did with `settings`).
Quick check per commit, full check per pipeline.

### 6. Split the big files

P7bis (1699 lines), P6bis (841), P2bis (543), P4 (534), P2 (498), the longest
Python `main()` functions (v2_region_plot 215 lines, v2_region_groups 190).
One file at a time, same checks.

The annotation GUI: the automatic-annotation layer moves out of LightSuite's
file into `registration/annotation_gui/` as a plugin behind one generic hook
(LS7), and the r key goes (L3). Checked with the sandbox's GUI drive scripts
pointed at the branch, the engine's self-test, and a checklist Giulio runs on
MG914 in the check tree (keys a, j, k, K, u, U, save; with nothing
edited, `atlas2histology_tform.mat`, `plane_anchors.mat`,
`auto_atlas_planes.mat` and `annotation_provenance.mat` come out identical),
plus one opening with no engine installed.

### 7. Style pass

The reference examples are named from the first files restyled and approved.
Then about eight groups of files, each restyled against `docs/STYLE.md`.
Comment and layout changes are checked automatically for code identity; small
code changes by rerunning. Every possible bug noticed goes to the list below.
A final check that no file credits an author other than Giulio.

### 8. Fixes

Each bug of the list below gets its own commit and a before/after check on the
check tree. Fixes that would change production numbers wait for Giulio's
decision.

### 9. Scientific additions

A1, A2, A3, A4, A5 in that order (A10 if wanted), each with its before and
after; the
quoted numbers in the documents are corrected. Then P8 to P10 go to
`archive/`.

### 10. Tests and documents

Tests, where mistakes are easy and silent: cohort registry, atlas lookup,
left/right statistics, the path test (from step 4), the Python route's region
means on synthetic volumes (plant a known pattern, recover it), the automatic
annotation's self-test. Skipped, not failed, when data is not connected.

Documents, on the model of the imaging repository (Python READMEs follow the
sections of the Genedata and MaxWell projects, MATLAB ones the imaging
repository's):

- **README.md**: an image with a caption (a coronal zref plane or a flatmap),
  the project in one paragraph and the four lines of work, setup (MATLAB
  R2024b, R2022b at least for `fitgeotform2d`, with the Image Processing,
  Statistics and Parallel Computing toolboxes; elastix 5.1.0 on the path,
  installed per user (the recipe copied in step 4); npy-matlab; the
  BrainGlobe atlas cache in `~/.brainglobe`; the three Python environments
  (analysis, flatmaps, automatic annotation), plus `tools\venv_dev` for tests
  and linting), the pipelines table, where outputs go, conventions.
- **docs/SCIENTIFIC_CONTEXT.md**: the question; the four lines of work and what
  each found; the grant this serves (the imaging repository describes this
  project as part of Aim 2 of the SNSF Weave grant: the two documents cite each
  other); related literature grouped by the question each paper bears on,
  with what it shows and what it means here (L7; the PDFs stay in
  `data\ref_papers\`, whose `PAPER_SUMMARIES.md` is the start); a
  question-to-code table; the terms
  used in the code (readings, nano, auto, DIFF, LR-sum, the atlases, ages).
- **docs/ADDING_DATA.md**: a new mouse end to end, manual steps included (slice
  order, artifacts, anchors, review), for young mice and for more plasticity
  animals; the data layout and channel names; the old-to-new script names; the
  traps: channel names change between stages; control points follow slice
  positions, so a reorder needs `remap_control_points`; never re-align an
  annotated brain; never put an atlas folder on the path (LightSuite finds the
  atlas by `which('average_template_10.nii.gz')`, the 20 um DeMBA volumes are
  renamed `*_10` and need `px_atlas = 20`, so a second atlas folder or a wrong
  `px_atlas` halves the AP scale without an error); a slice still orange at
  save has no control points and registers from images alone; which
  LightSuite version processed which cohort.
- **docs/ROADMAP.md**: open questions and the decisions taken (the
  normalisation of the two routes and how they would be made comparable for
  one paper, atlas choice,
  what the SEP channel is, the beyond-abundance result, the ISH tests, the
  known defects of the retiring route, the paused plasticity line and what
  would make it convincing).
- **docs/FIGURES.md**: every figure or number that went into the grant or will
  go into a paper, with the script, settings and tag that make it.
- One README per pipeline folder (run order, outputs, where the code is),
  and for `tests/`, `tools/` and `third_party/`.

Sources for the scientific context: the grant and the papers in
`data\ref_papers\` (added 30 Sep); for the plasticity comparison, Giulio's
account (a small increase of the nanobody signal in S1 after RWS, as expected
for plasticity localised in S1) and its outputs in `data\comparisons\`; for
the ISH work, which was never presented, the output figures of the final P9
run.

### 11. Report and merge

What changed, how each step was checked, the bug list, the open decisions.
Merge as described under Merging. Then the remaining young mice through the
refactored code, and MG914's SEP channel and Python route with them.

## LightSuite

**What we have.** At our first commit, 151 of our 153 LightSuite files were
identical to upstream commit `2f16206` (29 Oct 2025). Two edits predate that
commit: the demo script's paths (reverted to upstream in step 4) and a PNG
export in `alignSliceVolume` that upstream has since fixed its own way (dropped
at the swap). Upstream is 168 commits ahead and active. Its README says a
Python version "is being developed at PyLightSuite"; that repository was not
public on 30 Sep.

**Our changes**, as a patch series, each small enough to be one upstream pull
request:

| | change | kind | effect on results |
|---|---|---|---|
| LS1 | `registerSlicesToAtlas`: 0x4 placeholders when there are no control points | crash fix | none today (the branch ran only in the 2 Sep image-only trials, since redone) |
| LS2 | `alignSliceVolume`: atlas resolution from `px_atlas` instead of a hard-coded 10 um | registration fix | required for the 20 um DeMBA atlases (the AP scale was off by 2); none for the adults |
| LS3 | GUI fixes: clicks on points, dragging, partial-pair title, no empty control-point file on save | GUI | none |
| LS4 | GUI: atlas-plane prediction without backward extrapolation, order-conflict warning | GUI | none when every slice has points |
| LS5 | GUI: editing, carry-forward and t, numbered labels with flags, controls window, the `*` outlier mark | GUI | none on registration code; a slice left orange is saved empty and registers from images alone |
| LS6 | SliceOrderEditor montage | GUI | none |
| LS7 | GUI: one generic plugin hook, off by default (new, written in step 6) | GUI | none |
| LS8 | GUI: the r key's landmark_refine proposals | project-specific | removed in step 6 (L3); t and p stay |

Both registration bugs are still in upstream. In step 6 the automatic-
annotation layer (about 550 lines today inside LightSuite's GUI file: anchors,
proposals, accept and re-propose keys, provenance) moves into our tree as a
plugin that P4 passes in when the engine is installed. Without the plugin, the
GUI behaves as our patched GUI (as upstream once LS3 to LS5 are merged). The
project thresholds (the suggested anchor count, the outlier rule) become
plugin settings.

**Upstream pull requests** (Z2): LS1 to LS7, after rebasing on upstream (its
GUI now has its own numbering and binds backspace to "delete last point";
its SliceOrderEditor was rewritten). Wish-list requests from our workarounds:
saving the per-slice transforms, an output-folder option, relative elastix
paths, a seed in the align step, the demo's settings bug, and extentfactor, the
CPD iteration count and the B-spline bins and samples as options with
upstream's defaults, so a cohort can be registered with the old parameters
(needed for Z1's first option). The patch mechanics (a fork with the series,
vendored with `git subtree`) and the full swap plan from the study go into
`PATCHES.md` in step 4.

**The swap** (Z1), as its own step after the refactor. Changes that alter
results: upstream's B-spline settings (every registration), the align extent
(6 to 10: plane indices shift by 30, so re-aligning an annotated brain
invalidates its points), a lower CPD iteration count, and extraction that no
longer crops and edits the curated decisions file (which would invalidate
artifact masks and P2 outputs). Also: P1bis must use upstream's SliceOrderEditor
(ours drops its crop columns on save); P4bis gains the crop step or refuses;
`cpwt` stays 0.2 (upstream's demo uses 0.4); the path starts at
`<LightSuite>/src`; our matlab_elastix is checked against the fork upstream
now asks for. So: never re-run extraction or align on an existing brain;
measure the noise floor of the current code first (register only); register
the reference mice with upstream from their existing align outputs and control
points and compare voxelwise, per region, and on the headline numbers;
re-validate the automatic annotation, whose models were trained on the current
geometry; record which version processed which cohort.

## Bugs and oddities (fixed in step 8, each on its own, unless another step is named)

- P7bis takes the behavior mice by index `[1 2 3 4]` from whatever P6bis
  saved, while its comment names three mice: select the same four by name
  and correct the comment (S6: no change in results).
- P6bis writes diagnostics to the current folder when its diagnostics folder
  variable is missing (to become an error).
- P4 `register` (and P4bis) take the grid of the registered volume from
  `sliceinfo.mat`, which P1 writes from `local_settings.txt`. Every young brain
  was extracted in August with `px_atlas = 10`, so all of them, the ones still
  to register included, come out on the same 10 um-sampled grid that the
  Python route expects (it averages 2x2x2 blocks down to 20 um). A P1 run now
  would read `px_atlas = 20` and give a 20 um grid without any warning. To fix
  before any new brain is extracted: P4 sets the output grid explicitly and
  checks it against the rest of the cohort; `ADDING_DATA.md` names the trap.
- `verify_demba_setup.m` checks P4's source text for `atlas_key =
  'demba_p20';` and fails today (to check behaviour instead).
- `test_backward_compat.m` hard-codes 10 `D:\` paths (two more in comments)
  and adds `main`'s code folder to the path (fixed in step 4, before use).
- `v2_ish_compare` skips the old-against-new check silently when P9's summary
  is missing (to fail instead).
- `v2_video.MIN_N` has no `young_P22` key (KeyError; fixed here, not in A10).
- `v2_ish_words` gives slightly different bootstrap intervals (`gap_lo`,
  `gap_hi`) from run to run: one random generator is shared across features in
  the order they were collected from Python sets, which changes with each
  run's string hashing (step 4 check: medians 0.004 apart, at most 0.19; every
  feature, p and q identical). Iterate the features sorted.
- The automatic annotation's matcher is not reproducible on the GPU: two runs
  of the same code on MG914 moved the proposed points by a median of 0.07 px,
  at most 13 px (atlas landmarks identical). Ask PyTorch for deterministic
  algorithms in the engine; the proposals are reviewed by hand anyway.
- `run_group_differences` (P7bis): with `perform_area_based_analysis_coarse`
  on and `perform_area_based_analysis_fine` off, the coarse block reads
  `half_width` before anything sets it. Both are off in production. Its
  comment beside the behavior subset still says "subselect 3 of ... 4".
- **First of step 8, a safety bug.** `run_order_slices` 'apply':
  LightSuite's `generateReordedVolume` takes the decisions file,
  `volume_for_ordering.tiff` and its output `volume_ordered.tiff` from the
  absolute `procpath` and `volorder` that P1 stored in `sliceinfo.mat`, not
  from the cohort folder. On a copied tree it reads and overwrites the
  original folder, which the data-root guard does not see (or fails if that
  drive is absent); if the decisions file is missing there it keeps the
  original order without a word. Fix: set `sliceinfo.procpath` and
  `sliceinfo.volorder` from the cohort folder before the call, as
  `register_to_atlas` already does for `opts.procpath`. Until then 'apply'
  is never run on a copied mouse folder (its help says so).
- `compare.draw_figures` (young against adult) titles each slice
  `plane {zc * 2 + CCF_AP0}`, which reads 180 planes too high (526 to 1074
  instead of 346 to 894 at 10 um): every figure shown with those titles
  carries wrong plane numbers; the maps themselves are right.
- P8 colours its reliability bars with the first 200 levels of
  `flipud(gray(256))`; its comment says this avoids white at the low end, but
  level 1 is white and the dark end is cut, so the least reliable structures
  are white bars on white. P8 retires; A4 must not copy the scale.
- `ArtifactAnnotator` (P3's GUI, unchanged since before the refactor) breaks
  when Esc is pressed again while its close dialog is open or the window is
  slow: the queued key presses run the close handler on a deleted figure
  (`uiresume(src)` on an invalid object), and a second launch can delete its
  window before `uiwait`. Seen in the step 6 hand check (2 Oct). Not fixed:
  P3 has not been used in production, and Giulio would rather replace the
  annotator with something better or retire it (ROADMAP).
- The annotate mode's "the automatic annotation is not installed" line is
  easy to miss among the start-up messages (Giulio, 2 Oct): make it a
  `warning` or a banner.
- `run_nano_equalisation` (P2bis) stops when `save_results` is false:
  `timestamp` is set only in the save branch and the first video needs it.
  Its two videos per mouse cannot be switched off. It also saves the
  inter-quartile range with the wrong sign (25th minus 75th percentile) as
  `stats_intensity_iqr_*`; nothing reads it.
- `run_residual_correction` (P2) stops when `doPlotBkg` is false
  (`select_reference_pixels` then never assigns its figure output), and
  loads the atlas annotation without using it.
- `explore_czi_G` reads `globalMeta` before setting it; the error is caught,
  so the first file reports no scene information and later files the
  previous file's.
- `add_sep_channel.m:279`: since step 4 its panel title names
  `run_register_to_atlas`, and the TeX interpreter draws the underscores as
  subscripts; give the title `'Interpreter', 'none'`. Also in the register
  code, an `annotated{end+1}` keeps the script form `%#ok<SAGROW>`, which a
  function does not honour (style pass).
- `v2_adult_arms` stops at its self-check before drawing
  `arms_consistency.png`, so the figure that would show a drift is missing
  exactly when the check fails (the table is written).
- `v2_adult_arms`' consistency check against `v2_region_plot` compares values
  stored to 4 decimals with `<=` half the last digit (5.0e-05) and no margin for
  floating-point error, so it fails when a difference lands exactly on the
  bound (5.000e-05 in the step-3 reference run); the table is written before
  it stops.
- `v2_video`'s reliability t uses n = max(nL, nR), which overstates n when the
  two hemispheres come from different mice.
- The `ratio` reading's level depends on exposure: the zero line in
  `region_plot.png` and the ratio video's scale carry no "nano equals
  autofluorescence" meaning; only the order within a mouse does (a caption
  fix at least).
- `adult_profile` is written three times in the ISH scripts, with no
  minimum-mice rule (replaced in A3).
- `v2_ish_regions`' docstring describes an orientation test that is not in the
  code.
- LightSuite's copy holds scripts that write to `S:` (`compare_mice.m`,
  `protocol_manuscript_generate_plots.m`) and a machine-specific demo: kept,
  marked do-not-run in `third_party/README.md`.
- A nested duplicate `matlab_elastix-master/matlab_elastix-master/` sits on
  disk (ignored): Giulio deletes it before the merge, or a stray copy is left
  on `main`.
- elastix 5.1.0 is installed per user and recorded only in
  `bk/LightSuite.txt` (copied into the README in step 4).

The style pass will add to this list.

## Progress

- **30 Sep, step 0** (`6a29c37`, `c2e4cdc`, `8f16f14` on `main`): the data-root
  variable and its guards in both languages, the P4 align guard, the
  verification tools, frozen environments, the plan documents.
- **30 Sep, step 2**: snapshot refreshed (2,156 files, 7.2 GB, nothing
  deleted); tags `refactor-start` (`c2e4cdc`) and `grant-2026-09` (`04c0484`)
  pushed; check trees `G:\sep_refactor\ref` and `G:\sep_refactor\check`.
- **1 Oct, step 3**: the reference run of the old code on the small set, all
  stages. Old code bugs met on the way: `v2_adult_arms`' self-check (in the
  bug list).
- **1 Oct, S6, informative**: today's P7bis on the inputs of the approved
  December 2025 figures reproduces them (slab t maps and surprise masks
  correlate 0.997 to 0.998 for RWS, 0.988 to 0.999 for behavior; individual
  maps 1.000000; regional bars 0.994 to 0.998). The residue comes from the
  background masks, which were regenerated since. The S1 increase is there.
- **1 Oct, step 4** (branch `refactor`: `a50bc66` pure moves, `3a735f0` paths
  and references, `e37d6bc`, `684ce64`): passed. Code identity shows only the
  intended files changed; the path test finds no clash. Run on the check tree
  and compared with the reference:
  - plasticity chain: every `.mat` and every P7bis output identical; 11
    diagnostic PNGs differ by 1 or 2 anti-aliasing pixels;
  - Python route: 817 of 825 files identical, and the same step fails
    (`v2_adult_arms`); the rest is expected or old-code non-determinism: the
    `*_scalars.npz` caches store their source's file date, the diagnostics
    index and sheet 08's title name the renamed scripts, `v2_ish_words`'
    bootstrap bounds move with string hashing (bug list);
  - registration of MG914: identical voxel for voxel (old, old and new), its
    transform file identical; the automatic proposal moves by a median of
    0.06 px between old and new, less than between two runs of the old code
    on the GPU (0.07 px median, 13 px at most; bug list).
- **1 Oct, step 5, plasticity chain** (`ca5b2f0` pure moves, `a984279`):
  `run_collect_by_group`, `run_normalise_groups` and `run_group_differences`
  keep their settings and call `collect_by_group`, `normalise_groups` and
  `group_differences` in `group_comparison/pipeline/`, whose bodies are the old
  scripts' bodies (parse trees identical); P7bis's behavior subset is the
  setting `behavior_subset`. Fresh run on the check tree against the
  reference: all 12 `.mat` files and every P7bis output identical, 7
  diagnostic PNGs differ by 1 or 2 anti-aliasing pixels; `sep_test_path`
  passes. The Python stage now pins `PYTHONHASHSEED=0`, and the reference's
  `v2_ish_words` outputs were regenerated with it by the old code (identical
  on two runs; the unseeded originals are kept in `G:\sep_refactor\ref_unseeded`).
- **1 Oct, step 5, Python route** (`add8135` pure moves, `f473cb5`): the
  `v2_*.py` scripts are the package `mapping/sepmap/` (`volumes`,
  `young_vs_adult`, `adult`, `ish`, `diagnostics`, `config`), with one
  `mapping/run_<step>.py` per step and `mapping/settings.toml`. In it for
  now: the constants two modules kept in step by hand (all copies were equal)
  and the two ISH passes, chosen with `--panel targets|ontology` instead of
  `V2_ISH_PANEL`/`V2_ISH_TABLE` (refused if set). Each run prints the
  settings in force. Reviewed from three sides before the run (one minor
  finding). Full route on the check tree against the reference, 120 min: the
  same step fails (`run_adult_arms`, the known self-check); 259 output files
  identical; the 8 that differ were each checked to differ only where a
  script is named (diagnostics README, sheet 08 and `slices_sepratio`
  titles) or in the scalars cache's source date. Every other parameter is
  still a constant in its module, and `matplotlib.use` is still in the
  modules (style pass).
- **1 Oct, step 5, registration** (`89a24c4` pure moves, `1d240a0`):
  `run_register_to_atlas` and `run_add_sep_channel` keep their settings and
  call `register_to_atlas` and `add_sep_channel` in `registration/pipeline/`
  (parse trees identical to the old bodies). Four variables of the align
  branch come from `load` without an output; none is a function on the path,
  and the function's help says so for step 6. `verify_demba_setup` reads the
  moved settings code. Reviewed for the GUI modes, which cannot run headless:
  their windows keep their own state, nothing reads the old script's
  variables. MG914 on the check tree against the reference: registered
  volumes identical page for page, `transform_params.mat` the same, the
  automatic proposal within the old-against-old spread (median 0.06 px).
  Still to try by hand, after step 6 changes the GUI: annotate mode on MG914
  in the check tree (the step 6 checklist).
- **1 Oct, step 5, preprocessing** (`6ed2472` pure moves, `15c2f2e`): the
  six drivers (`run_copy_raw_data`, `run_extract_and_center`,
  `run_order_slices`, `run_residual_correction`, `run_nano_equalisation`,
  `run_annotate_artifacts`) keep their settings and call one function each
  in `preprocessing/pipeline/` (parse trees identical to the old bodies);
  `explore_czi_G`'s folder is a setting at its top. P2 and P2bis run on
  MG914 alone, old code against new on identical inputs in
  `G:\sep_refactor\pre\{ref,check}`: all 7 `.mat` files, the 4 videos and 94
  of 102 figures identical, 8 per-slice PNGs differ by 1 or 2 anti-aliasing
  pixels. P0 and P1 by code identity and review; the two GUIs (slice order,
  artifacts) reviewed in the code: each blocks until its window closes and
  keeps its own state. The QC scripts and the tools `make_ordering_volume`
  and `make_atlas_reference_sheet` stay scripts (hand-run audits; headers
  in the style pass). To try by hand with the step 6 checklist: both GUIs on
  `G:\sep_refactor\pre\gui\data` (no `sliceinfo.mat` there on purpose, so
  'apply' cannot run).
- **1 Oct, step 6, the big files** (`24c41cf`, `e88475b`, `d67d8a2`,
  `5e67ade`, `0a8fce9`, `510bc6a`, `586430f`): `group_differences`,
  `normalise_groups`, `residual_correction`, `nano_equalisation` and
  `register_to_atlas` are a short main function of `%%` steps calling local
  functions with explicit inputs and outputs; `region_plot` and
  `region_groups` have a short `main()`. Every old statement is still there
  once, unchanged (checked per statement). Split in four worktrees, reviewed
  (no findings), checked one at a time on the check trees: plasticity 135
  files, every `.mat`, `.fig` and P7bis output identical, 10 PNGs differ by
  1 to 3 anti-aliasing pixels; P2/P2bis on MG914 identical apart from 9 such
  PNGs; MG914's registration identical page for page, the proposal within
  the GPU spread; the Python steps and their readers, 72 files identical
  (3 EPS differ in their creation date). Not run: the angle, annotate and
  align modes (their code moved into `load_regopts`, `set_cutting_angle`,
  `annotate_control_points`, `refuse_annotated_mice`, `bridge_preprocessing`
  and `align_slices`, statements unchanged), for the hand check.
- **1 Oct, step 6, the annotation GUI** (`303c566`, `d8356d6`, `b83da33`,
  `dfef3b5`): the r key and landmark_refine are retired (L3, LS8; files in
  `archive/`, P4 annotate no longer starts the worker). The automatic
  annotation left LightSuite's GUI file (-505 lines) for
  `registration/annotation_gui/auto_annotation_plugin.m`, behind one generic
  hook (LS7: `value = plugin(event, gui_fig, gui_data, value, info)`, events
  open, key, planes, title, labels, window, edit, save; without a plugin every
  call hands its value back), with its settings in `annotation_settings.m`.
  `annotate` passes the plugin when `auto_annotate('check')` finds the
  engine, otherwise it says so and opens the plain GUI. Checked by driving
  the old and new GUI headless through the sandbox's drive scripts plus
  three new ones (hand annotation, reopening a review, save unchanged), with
  and without the engine: the same files, apart from click times and the
  GPU proposal's usual spread. Still by hand: `G:\sep_refactor\gui_check\HAND_CHECK.md`.
  **For the merge (step 11):** move `code\auto_annotation\.venv` to
  `registration\auto_annotation\.venv` (git does not move ignored files) and
  run its self-test before the first `annotate`, or annotate opens without
  the automatic layer; the leftover `code\landmark_refine\.venv` can be
  deleted (the root `.gitignore` now ignores both).
- **2 Oct, step 6 hand check** (Giulio, `G:\sep_refactor\gui_check\HAND_CHECK.md`):
  passed. Annotate mode with the engine (open a reviewed brain; review a
  proposal with k, u, j, a, U, K, t, p, s), without the engine (plain GUI,
  points only), the cutting angle (refused with points; set and saved on a
  copy), the slice order editor (decisions file unchanged) and the artifact
  annotator (saves; its old close-handler weakness, bug list) all behave as
  before; every difference in the saved files came from the keys pressed.
  The `s` key saves only the points, the close dialog also the affine, as
  in the old GUI.
- **2 Oct, step 7 decisions** (Giulio): `docs/STYLE.md` and the four
  reference examples approved, with the seven open choices as recommended:
  `ruff format`; MATLAB settings stay lower case; `clear all` becomes
  `clear; clc; close all;` in the drivers, QC scripts, tools and test (kind B);
  `sep_setup_paths` once per session; `collect_by_group`'s dead code goes in
  step 8; the structural changes (kind C) are part of step 7, each its own
  commit checked by a rerun; local variables may be renamed (kind B),
  saved names, columns and settings keep theirs. The procedure is in
  `docs/STYLE_PASS.md` (temporary).
- **2-3 Oct, step 7, kinds A and B** (refactor `742d7b8`): nine groups
  restyled in parallel worktrees (comments and layout proven code-identical
  by the identity tools; small code changes in their own commits: `clear all`
  in the drivers, QC scripts and test, local renames, isort, type hints,
  `assert` into `raise`), merged, evened out by a consistency pass, the
  kind A commits in `.git-blame-ignore-revs`. Checks against the reference:
  plasticity (every `.mat`, `.fig` and P7bis output identical), P2/P2bis on
  MG914 (identical), MG914's registration and the GUI drive scripts
  (identical), the Python route (259 of 267 rewritten outputs identical, the
  8 known differences; sheet 08's longer title also moves its panels by a
  fraction of a pixel, proven by redrawing it with the old title).
- **3 Oct (overnight), step 7 kind C and step 8 applied fixes** (refactor
  `47c252f`): Python: one plotting module, the backend set by the run
  scripts, duplicate helpers merged, 63 analysis constants in
  `settings.toml`, `pathlib`, about 25 long functions split, rows read by
  column name, `build_demba_atlas` on argparse; MATLAB: `common/sep_palette.m`,
  duplicate helpers merged, long functions split. 49 fixes that change no
  production output (the `run_order_slices` 'apply' safety fix first; the
  P2/P2bis/P7bis crash fixes; dead code; the "not installed" warning; the
  corrected "not permeabilised" header; tool robustness). Checks against the
  reference all passed (plasticity: every `.mat`, `.fig`, P7bis output
  identical; P2/P2bis identical; MG914 registration identical, the proposal
  in the GPU spread; GUI drive scripts the same; Python 259 of 267 identical,
  the 8 known differences). 21 fixes that change an output wait on
  `step8-pending` for Giulio, listed in `G:\sep_refactor\MORNING_REPORT.md`
  and `G:\sep_refactor\FIXES_STEP8.md`; among them two that change results:
  subref's reference never excluded fiber tracts and ventricles (the
  exclusion names never matched; fixing it moves the young-adult subref
  difference by -0.19 log2 and its q < 0.05 structures from 104 to 86; the
  other readings are unchanged), and the cohort mean counting MG897's 183
  missing sepratio voxels as zero. Open, no commit: P7bis smoothing sets
  voxels outside the tissue to 0, not NaN, so the zeros enter the group
  means of the approved December 2025 comparison.
- **3 Oct, Giulio's answers to the held fixes** (`G:\sep_refactor\MORNING_REPORT.md`):
  yes to 1-9 and 11-21, no to 10 (the automatic annotation stays on the GPU:
  speed matters more than exact repeatability). With them: 1 also takes
  "brain-unassigned" and "unassigned" out of subref's reference; 9 checks the
  whole surprise-bar region list for other assembly errors; 11 also selects
  the behavior mice by name; 15 gets a short title saying what the ratio
  panel's zero means. New fixes, now: 23 P7bis smoothing NaN-aware (changes
  the approved December comparison: measured on the S6 inputs, the S1
  result reported before and after); 24 the video reliability t by exact
  per-mouse folding; 25 NaN-aware region sums ("every time we can we should
  be NaN aware"); 26 P4 register sets its grid explicitly; 27: the sections
  are not cleared (docstrings corrected), test_backward_compat's young count
  not pinned, explore_czi_G's SizeS on adult files, check_demba_to_allen's
  helpers out of `tmp\`, one colour convention for the ratio video, printed
  conclusions that follow the numbers (concise). Later: 22 (P2's slices with
  too few reference pixels: a guard, to discuss), output file names that
  record their settings (ROADMAP: renaming outputs breaks the readers of
  existing data). MG904 is P22 (the age in a raw folder's name is the
  mouse's age); the grant figure only grouped it with the P20 brains, as
  Sami asked.
- **3-4 Oct, step 8 done** (refactor `cfb7cd9`; report
  `G:\sep_refactor\STEP8_REPORT.md`): the accepted held fixes and the new
  fixes 23 to 27, every check passed with each output difference traced to
  its fix (plasticity, P2/P2bis, two Python runs, MG914's registration rerun
  on 4 Oct after C: had filled on 3 Oct with MATLAB and elastix temporary
  folders: stacks, elastix folders and `transform_params.mat` identical, the
  proposal in the GPU spread). Fix 23 on the December 2025 inputs: the RWS
  S1 increase holds (L-R map the same 292 pixels; barrel field first in the
  L-R bars, 13th to 5th in L+R; per mouse +26%, p 0.12, unchanged);
  naive vs behavior changes most, because MG709 has no tissue at slab 565
  but was counted as a fourth mouse of zeros. Open decisions after the
  merge: a minimum number of mice per voxel for P7bis's t (edge voxels with
  2 per group reach |t| 73), the surprise-bar region list (nesting, sums),
  the video mean panel, two more NaN issues (per-brain backgrounds, the young
  brains' warp over zeros). Documents (step 10) on `step10-docs`, synced.
- **4 Oct, temporary folders:** MATLAB (`tp*`, from niftiread and others) and
  elastix (`transformix_*`) leave folders in `%TEMP%` that nothing deletes:
  32,279 folders, 255 GB on C:, deleted with Giulio's OK. They filled C: on
  3 Oct and stopped a registration. For after the merge: the detached runner
  points TMP/TEMP to a scratch folder on the data drive and removes its
  `transformix_*` folders after each registration (bug list).
- **Merge plan (step 11, with Giulio):** a dry run of `step10-docs` into
  `main` conflicts in two files only: `docs/REFACTOR_PLAN.md` (keep `main`'s)
  and `tools/run_matlab_detached.ps1` (keep the branch's, after checking it
  holds `main`'s runner fix). Giulio closes every MATLAB on the code; his 3
  uncommitted settings lines in `P4_register_to_atlas.m` are saved as a patch
  and set again in `registration/run_register_to_atlas.m`; the engine's
  `.venv` moves to `registration/auto_annotation/` and passes its self-test;
  the old `landmark_refine/.venv` and the nested `matlab_elastix-master`
  duplicate go; `sep_test_path` and the Python imports on the production
  code; push.
- **4 Oct, step 11, merged** (`03e1ca0` on `main`, with Giulio): `step10-docs`
  (the refactor and its documents) merged into `main`; conflicts as the dry
  run said (the plan kept from `main`, the runner from the branch). Giulio's
  three uncommitted settings set again in `registration/run_register_to_atlas.m`
  (MG914, annotate, demba_p28; still uncommitted, patch in
  `G:\sep_refactor\merge\`). The engine's `.venv` moved to
  `registration\auto_annotation\` (self-test: models loaded on cuda); the old
  landmark_refine worker stopped by Giulio; the leftovers deleted (the old
  `landmark_refine` and `auto_annotation` folders, the nested
  `matlab_elastix-master`, identical to `third_party/matlab_elastix` apart
  from line endings). On the production code: `sep_test_path` 0 failures,
  `get_cohort('verify')` passes, data root `D:\sep_histology\data`, 32
  sepmap modules import, 26 run scripts answer `--help`. Next: the
  production reruns the fixes call for, the open decisions, A1 to A5.
- **New order** (Giulio, 2 Oct): merge first, new science after. Step 7's
  kind C and the step 8 fixes go in one batch with one set of reruns (a fix
  that changes production numbers still waits for Giulio); then step 10's
  documents and step 11's merge. Step 9 (A1 to A5, then retiring P8 to P10)
  moves after the merge, as normal project work on the merged code; P8 to
  P10 stay in place and keep running until A1 to A5 replace them.
- **For step 10** (Giulio, 2 Oct): the README's cover image is the slice-order
  montage on Giulio's desktop (`slice_order_montage.png`, 1 Oct), copied to
  `assets/` and shown on the README's first lines as the imaging repository
  does (`![...](assets/example.png)`); Giulio can also set it as the GitHub
  social preview (repository settings).
