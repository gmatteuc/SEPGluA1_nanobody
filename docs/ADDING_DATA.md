# Adding new brains

How a new brain goes from its raw `.czi` files to the analyses: copied,
extracted, ordered and corrected (MATLAB, `preprocessing/`), registered to the
atlas of its age (MATLAB, `registration/`), then measured and compared
(Python, `mapping/`). The steps that need someone at the screen are marked
**manual**. A young brain and a new adult for the plasticity comparison take
the same path; where they differ, both are given.

Each MATLAB driver is a script: open it, set its `%% Settings`, run it. Run
`sep_setup_paths` once per MATLAB session first, from the code root. The long
stages (registration, the SEP channel) can run detached with
`tools\run_matlab_detached.ps1` ([`tools/README.md`](../tools/README.md)). The
Python scripts run from the code root in `tools\venv_atlas`, `run_closeup.py`
in `tools\venv_flat`.

## Where things are

- `get_paths` (MATLAB) and `mapping/sepmap/config.py` (Python) take the data
  root as the code folder's sibling `data`: on this machine
  `D:\sep_histology\code` and `D:\sep_histology\data`. The environment
  variable `SEP_DATA_ROOT` moves the whole data tree, inputs and outputs, for
  a run on a copy. A copy of the code is refused the production data.
- The raw files stay on the lab share, which is read only. The pipeline
  copies them and works on the copy.
- Each brain has one folder, `<data>\<group>\<mouse>\`, with the groups
  `naive`, `rws`, `behavior` and `young`. Its `lightsuite\` subfolder holds
  every stage:

```
<data>\young\MG914_SepGluA_P28\
  *.czi                              the raw files, copied from the share
  local_settings.txt                 extraction and atlas settings (young brains)
  lightsuite\
    sliceinfo.mat                    extraction
    volume_centered\                 chan01_DAPI, chan02_Cy5 (nano), chan03_Cy3 (auto), chan04_EGFP (SEP)
    volume_for_ordering.tiff         the colour composite the slice order is curated on
    volume_for_ordering_processing_decisions.txt
                                     the curated order, flips and discards
    slice_order_montage.png          the curated order, one tile per slice
    correction_output\               corrected and scaled volumes, equalized_volume.mat
    volume_centered_processed\       DAPI, NANO, AUTO, DIFF, MASK: what the registration reads
    regopts.mat, volume_aligned\, volume_for_inspection.tiff,
    <mouse>_dim<d>_initial_registration.png                      'align'
    cutting_angle_data.mat           'angle'
    plane_anchors.mat, auto_atlas_planes.mat, auto_proposal_*.mat,
    annotation_provenance.mat        the automatic annotation
    atlas2histology_tform.mat        the control points
    transform_params.mat, elastix_forward\, elastix_reverse\,
    volume_registered\               'register'
    volume_aligned_sep\, volume_registered_sep\                  the SEP channel
```

The atlases sit beside the groups: `atlas\` (the Allen CCF, 10 um, with the
ontology tables) and `atlas_demba_p<age>\` (the DeMBA atlas of one age).

## 1. Name the brain and add it to the cohort table

- A mouse's name is the name of its raw folder on the share, and it is the
  name of its folder under the data root: `MG<n>_SepGluA_P<age>` for a young
  mouse, `<id>_Gria1` for an adult.
- A young mouse's age is the one in that name. MG904_SepGluA_P22 is P22: a
  grant figure grouped it with the P20 brains, as Sami El-Boustani asked, but
  it is registered and analysed at P22.
- The adults' exact ages are not recorded; all are well above P60. Their age
  is left empty in the table, and every adult is registered to the CCF, the
  P56 reference.
- The groups: `naive`; `rws`, an RWS session, then perfusion; `behavior`,
  perfused on the first day of above-chance performance in a single-whisker
  detection task; `young`.

Add one row to `common/cohort.csv`, which `get_cohort` (MATLAB) and
`sepmap/config.py` (Python) both read:

- `name`, `group`, `age_days` (empty for an adult);
- `share_subdir`: the subfolder of the mouse's share folder that holds its
  `.czi` files (`Anatomy\Axioscan` for some brains), empty when they sit at
  its root;
- `mapping_cohort` and `mapping_order`, the Python route's: the cohort the
  brain enters there and its place in the route's stack. Leave them empty
  until step 4, which gives the values.

The first 17 rows, the adults of the plasticity comparison, keep their order:
`run_normalise_groups` and `run_group_differences` select adults by position
(`selected_mice_idx_list`), and `get_cohort('verify')` stops a run when that
order has changed; `run_group_differences` picks the behaviour mice it
analyses by name (`behavior_mice`).
Add new rows at the end. `tests/test_backward_compat` checks the adults'
order, that the young cohort `get_cohort` returns is the table's young rows
in their order (so a new row needs no change to the test), and that each
young brain's age matches its name.

## 2. Preprocessing (MATLAB, `preprocessing/`)

| step | driver | set | writes | look at |
|---|---|---|---|---|
| 1 | `run_copy_raw_data` | `mice_to_process`; `do_copy = false` for a dry run | `<mouse>\*.czi` | the size check printed for every file |
| 2 | `run_extract_and_center` | `mice_to_process`; the mouse's `local_settings.txt` | `volume_centered\`, `sliceinfo.mat`, `volume_for_ordering.tiff` | `volume_for_ordering.tiff` |
| 3 | `run_order_slices`, **manual** | one mouse; `run_mode = 'edit'`, then `'apply'` | the decisions file, `slice_order_montage.png`, `volume_ordered.tiff` | the counts 'apply' prints: flipped, reordered, discarded |
| 4 | `run_residual_correction` | `mice_to_process` | `correction_output\corrected_volume_*`, `scaled_auto_volume_*`, `diagnostic_plots\` | the fit of each slice in `diagnostic_plots\` |
| 5 | `run_nano_equalisation` | `mice_to_process`: the new brains only, never an aligned one | `correction_output\equalized_volume.mat`; `<data>\intensity_diagnostics\` | the profiles and the video of the run |
| 6 | `run_annotate_artifacts`, **manual, optional** | one mouse | `correction_output\artifact_mask_volume_slicewise.mat` | |

**Copy.** `robocopy` copies, never moves or mirrors, and an interrupted copy
resumes. A destination on the share's drive or on a network path is refused.
For a brain from a new batch of acquisitions, `preprocessing/explore_czi`
prints the channels of a `.czi` file: extraction names each channel's file by
its place and dye, and every later step opens the files by those names, so
the four channels must come in the same order as before (traps below).

**Extract.** About ten minutes per brain. The settings come from
`local_settings.txt`, the mouse folder's first, then `lightsuite\`'s; without
either, LightSuite's defaults, which are the adults' values (150 um sections,
`px_process = 5`, `px_register = 20`, `px_atlas = 10`, DAPI for the
registration). A young brain's file says `px_atlas = 20`, with the DeMBA
crop, by 'align' at the latest; `px_register` stays 20 for every brain: see
"The registration grid" below. A brain that fails is reported and the others
go on.

**Order, manual.** `'edit'` opens SliceOrderEditor on
`volume_for_ordering.tiff`: reorder, flip, mark slices for removal, save, close
the window. `make_atlas_reference_sheet` draws the atlas plates one section
apart (`<data>\young\atlas_coronal_reference.png`), to keep open beside it.
Then `'apply'` rebuilds the preview `volume_ordered.tiff`. The order itself is
applied by 'align', which reads the decisions file; steps 4 to 6 work in the
extraction order, so a change of order does not need them again. Settle the
order before annotating (traps below).

**Residual correction.** Slice by slice, the autofluorescence is fitted to the
nano on reference pixels (a robust line), giving the scaled autofluorescence
and the nano minus it. Both variants, `slicewise` and `global`, are written;
every later step reads `slicewise`. A slice whose slope is far from the other
brains' (`slice_data` in `scaled_auto_volume_slicewise.mat`) is one to look
at in `diagnostic_plots\`.

**Nano equalisation.** Each slice's nano is scaled so that its median follows
the moving median of the slice medians over 5 slices. The mice of one run are
stacked and padded to the largest, so a brain's result depends on which
brains run with it. The run rewrites `equalized_volume.mat` of every brain in
it, without a check, and 'align' reads that file (`use_equalized_nano = 1`).
So run new brains that arrive together in one run, as MG909 to MG914 were,
and a single new brain alone; never put a brain that is already aligned in
the run, or its registered NANO no longer matches the file it came from.

**Artifacts, optional.** Outlines artifacts by hand. The masks go only into
the MASK channel, which no current analysis reads. Three RWS brains (MG691,
MG692, MG693) have masks; no other brain has. The annotator has a known
weakness: pressing Escape again while its close dialog is open can break it.

## 3. Registration (MATLAB, `registration/`)

`run_register_to_atlas` runs one mode at a time. Set `mice_to_process` (one
mouse for the modes that open a window), `run_mode` and `atlas_key`:

- a young brain: `'demba_p<age>'`, its own age. The folder
  `<data>\atlas_demba_p<age>\` must exist; `atlas/build_demba_atlas.py <age>`
  builds it (in `tools\venv_atlas`). P16, P20, P22, P28, P32 and P36 are
  built;
- an adult: `'ccf'`.

The driver stops when the atlas's age is not the mouse's, and when the mouse's
`local_settings.txt` gives a `px_atlas` other than the atlas's resolution (20
for DeMBA, 10 for the CCF).

| order | mode | who | what it does | look at |
|---|---|---|---|---|
| 1 | `'align'` | auto | bridges the corrected volumes into LightSuite as five channels, applies the slice order, aligns the slices, fits the atlas rigidly | `<mouse>_dim1..3_initial_registration.png`; `registration/qc/check_registration_error` |
| 2 | `'angle'` | **manual**, optional | the cutting angle by eye: save the plane on 3 to 5 slices spread front to back, close the window to write it | |
| 3 | `'annotate'` | **manual** | with the automatic annotation: the atlas plane on the 4 suggested anchor slices (`j` jumps between them, the wheel moves the plane, `a` fixes it), save with `s` | |
| 4 | `'autoannotate'` | auto | proposes control points on every slice from the anchors, about 30 s on a GPU | |
| 5 | `'annotate'` | **manual** | review: the proposal loads orange, its least confident points marked `?`; `k` accepts a slice (`K` every orange slice), `u` re-proposes it at the plane on screen (`U` every orange slice), the usual tools fix points; save | a `*` beside a pair far off its slice's own affine; an `ORDER CONFLICT` banner |
| 6 | `'register'` | auto | elastix refinement of every slice, the registered volumes of the five channels | the control-point summary and the `CHECK` lines printed before elastix starts |
| 7 | `run_add_sep_channel` | auto, 30 to 45 min | the SEP channel through the same registration | its printed check; `<data>\comparisons_v2\processing_diagnostics\sep_channel\<mouse>.png` |

- **The angle.** Every adult had it set; it must come before the first
  annotation, since it changes the atlas block the control points are
  counted in. 'angle' refuses a mouse that has points.
- **The automatic annotation** needs its engine, installed once per machine
  with `registration\auto_annotation\setup.ps1`
  ([its README](../registration/auto_annotation/README.md)). Without it,
  'annotate' warns and opens the plain GUI: place points on every slice by
  hand (`t` takes the neighbouring slice's points, `p` carries them forward;
  both stay provisional until touched).
- **Saving.** Only accepted or touched slices are saved. A slice still orange
  at save has no control points and registers from the images alone, and so
  does a slice with fewer than 5 points. `s` saves the points; closing the
  window through its dialog saves the affine too.
- **Before a P20 brain's first annotation**, `registration/qc/verify_demba_setup`
  checks the atlas on the path, its label space, its crop and the mouse's
  local settings.
- **'register'** stops on a mouse without control points
  (`allow_image_only_registration = false`; an image-only registration is a
  diagnostic, not a result) and on a slice whose two point lists differ in
  length. It prints `CHECK` for a pair far off its slice's affine (a wrong
  structure, a left-right swap), for points placed on two atlas planes (the
  registration uses the first point's plane for all of them) and for a slice
  off as a whole (usually the wrong plane). Fix them in 'annotate' and
  register again.
- **`run_add_sep_channel`** takes `mice_to_process` or the environment
  variable `P4BIS_MICE`, skips a mouse that is not registered or already has
  its SEP volumes (`overwrite = false`), and re-applies the saved
  registration: nothing is refitted. Every slice should come back at
  r = 1.00000; one below `min_slice_corr` (0.99) is named.

## 4. The Python route (`mapping/`)

The route needs the SEP channel: `run_cohort` stops on a brain without
`volume_registered_sep\chan02_SEP.tiff`. Set the brain's row in the cohort
table:

- `mapping_cohort`: `young_P<age>` for a young brain (`young_P28` for MG914),
  `naive` or `rws` for an adult of those groups; behaviour mice stay empty,
  since the route does not take them;
- `mapping_order`: the next free number (1 to 17 are taken), so the brains
  already in the route keep their order.

A row with `mapping_cohort` set enters the route's brain list
(`config.mapping_mice`): the per-brain steps and `run_diagnostics` take it
when run with no mouse named. Then, from the code root:

```
tools\venv_atlas\Scripts\python.exe mapping\run_per_mouse.py <mouse>
tools\venv_atlas\Scripts\python.exe mapping\run_to_ccf.py <mouse>
tools\venv_atlas\Scripts\python.exe mapping\run_diagnostics.py <mouse>
```

`run_per_mouse` measures the brain on the atlas of its age; `run_to_ccf`
carries a young brain into the adult CCF with CCF Translator (about 2.5 min
per volume, five volumes: nano, autofluorescence, SEP and the tissue mask
twice) and only places an adult. Look at
`processing_diagnostics\01_tissue_<mouse>.png`, `02_levels_<mouse>.png` and,
for a young brain, `04_warp_<mouse>.png`, all under `<data>\comparisons_v2\`;
what to look for is in `processing_diagnostics\README.md`.

Then the cohort steps, in the order the header of every `mapping/run_*.py`
lists: `run_cohort`, `run_compare`, `run_region_plot`, `run_region_groups`,
`run_video`, `run_video_compare`, `run_closeup` (in `tools\venv_flat`), and
`run_diagnostics` with no mouse named. Of the steps after them, those that
read the region table of `run_region_plot` (`region_means_per_mouse.csv`:
`run_ish_compare`, `run_ish_roles`, `run_ish_panel_test`, `run_adult_arms`
and the four `run_beyond_*`) are rerun when it changes, and so are
`run_ish_words` and `run_ish_arms`, which read their outputs (see the zref
trap below).

A brain enters only the cohorts the code defines: `young` (P16, P20 and P22
pooled), `young_P20`, `young_P16`, `young_P22`, `naive`, `rws` and `adult`
(naive and rws). A brain of another age, such as P28, goes through the
per-brain steps, but it enters the cohort volumes and the region tables only
once its cohort is added to `mapping/sepmap/volumes/cohort.py`, to the young
groups of `young_vs_adult/region_plot.py` and `compare.py`, and to
`videos.min_n` in `mapping/settings.toml`. From then on it can move every
brain's zref (traps below). Whether the older ages are pooled with P16 to P22
or form an age series of their own is not decided yet.

## 5. More animals for the plasticity comparison

A new naive, RWS or behavior adult goes through steps 1 to 3 like any brain,
on the CCF (`atlas_key = 'ccf'`). It needs no `local_settings.txt`:
LightSuite's defaults are the adults' values. Then, in `group_comparison/`
([its README](../group_comparison/README.md)):

1. `run_collect_by_group` with the group in `mousetypes_list` and
   `age_filter = []`: the group's stack is rebuilt, mice in table order.
2. `run_normalise_groups`: its mouse lists and the selection of each group
   are written at its top; extend them. `SEP_COHORT_SPECS=rws,naive,behavior`
   normalises the three groups in one run.
3. `run_group_differences`: the same lists, and the behaviour mice to
   analyse by name in `behavior_mice` (among those step 2 saved); `exp_type`
   picks `'rws'` or `'behavior'`.

The comparison approved on 5 December 2025 is in
`<data>\comparisons\naive_vs_rws\` and `naive_vs_behavior\`; today's code
writes `naive_vs_<exp>_nano\` beside it (last on 5 October 2026). A new
naive or RWS adult also changes the Python route's adult cohort once its
`mapping_cohort` is set.

## Traps

- **Channel names change between stages.** Up to `volume_centered\` the files
  carry dye names (`chan02_Cy5` is the nano, `chan03_Cy3` the
  autofluorescence, `chan04_EGFP` the SEP); from `volume_centered_processed\`
  on, role names (`chan01_DAPI`, `chan02_NANO`, `chan03_AUTO`, `chan04_DIFF`,
  `chan05_MASK`), and `chan02_SEP` in `volume_registered_sep\`.
- **The registration grid.** Every registered brain, young and adult, is on
  one grid: 10 um sampling, twice the 20 um registration grid. The Python
  route expects it (it averages 2 x 2 x 2 blocks to 20 um). 'register' sets
  it itself (`registration/pipeline/registered_grid_um.m`, 10 um), whatever
  `px_atlas` extraction copied into `sliceinfo.mat` from
  `local_settings.txt`, and so does `run_add_sep_channel`. Before registering
  any brain, 'register' stops when a brain's `px_register` in `sliceinfo.mat`
  is not 20 (LightSuite places the slices on a grid of half `px_register`),
  and prints a line for a `px_atlas` there other than 10. 'align' takes
  `px_atlas` and the AP crop from the atlas itself, whatever the settings
  file says; the check before every mode stops only when a
  `local_settings.txt` exists and gives another `px_atlas` than the atlas's
  resolution (20 for DeMBA). So a young brain's `local_settings.txt` may say
  20, with the DeMBA crop, from extraction on: 'register' prints its line and
  registers on the 10 um grid. The young brains so far were extracted with 10
  (like the `local_settings_before_demba.txt` kept in the young folders) and
  their file changed to 20 before 'align'.
- **One atlas on the path.** LightSuite finds the atlas with
  `which('average_template_10.nii.gz')`, and every atlas folder holds a file
  of that name; the DeMBA folders hold 20 um data under it. `get_atlas` puts
  the folder a run needs on the path and takes the others off: never
  `addpath` an atlas folder. A brain aligned to DeMBA with `px_atlas = 10`
  would have its AP scale halved without an error, which is why 'align' takes
  the resolution from the atlas and the driver checks the settings file.
- **Control points follow positions, not slices.** `atlas2histology_tform.mat`
  holds one cell per position in the ordered stack. A reorder after
  annotating leaves each point on whatever slice now holds its old position.
  `registration/remap_control_points` carries the points across (a dry run
  unless told to write; points on a slice whose flip changed are dropped).
  It needs the decisions file as it was when the points were placed: copy
  that file aside before editing the order of an annotated brain.
- **Never align an annotated brain again, never extract a registered one
  again.** Aligning rewrites the atlas block that points, anchors and
  proposals are counted in, so 'align' refuses a mouse with any of them and
  names the files. Extraction rewrites `sliceinfo.mat` and the centred
  volumes everything after it was computed from, and it subsamples points at
  random, so a second extraction is not the first.
- **'apply' writes into the mouse's folder.** `run_order_slices` 'apply' sets
  LightSuite's paths to the mouse's folder under the data root. The code at
  the tag `refactor-start` took them from the absolute paths stored in
  `sliceinfo.mat`, and on a copied tree it wrote into the original folder:
  never run that version's 'apply' on a copy.
- **Never equalise an aligned brain again** (step 5 of preprocessing). The
  same brain equalised with other brains gives slightly different volumes,
  and a run rewrites `equalized_volume.mat` of every brain in it without a
  check, so an aligned brain in the run is left with registered volumes made
  from a file that has since changed. A single new brain is run alone.
- **Adding a brain moves zref.** The region tables give every brain's zref
  over the structures that every brain in the tables has, so a new brain can
  change that set and move every brain's zref, the adults' included, and with
  it the numbers of the ISH and beyond-abundance steps. Compare the tables
  before and after.
- **Old names in old files.** Logs written before the reorganisation of
  September and October 2026 name the old scripts
  (`<data>\young\_P4_register_to_atlas.log`); the table at the end gives the
  new ones.
- **LightSuite.** Every brain was processed with the copy in
  `third_party/LightSuite` (upstream commit `2f16206` plus the changes in its
  `PATCHES.md`). The adults were registered in November 2025, before the two
  registration fixes LS1 and LS2, which change nothing for a brain at
  `px_atlas = 10`; a brain on DeMBA needs LS2. A newer upstream changes
  registration results (`PATCHES.md`), so it is not swapped in for some
  brains only.

## Old and new script names

The refactor of September and October 2026 renamed the drivers from
pipeline steps (P0 to P10) to `run_...` names, in folders by pipeline. The
code before it is at the git tag `refactor-start`
(`git show refactor-start:P4_register_to_atlas.m`);
the code behind the grant figures of September 2026 at `grant-2026-09`.
Every file the reorganisation moved, helpers and vendored folders included,
is in [`refactor_name_map.csv`](refactor_name_map.csv), one row per file.
Since then the body of each MATLAB driver has moved into a function in its
folder's `pipeline/`, and `landmark_refine` into `archive/`; after the merge
(4 October), P8, P9 and P10 moved into `adult_matlab/` under `run_...` names,
with new headers and their computations unchanged. `sep_setup_paths` warns
about any `.m` file at the code root other than `get_paths.m` and itself (an
editor tab saved after a move recreates the old file there).

| before | now |
|---|---|
| `P0_copy_raw_data.m` | `preprocessing/run_copy_raw_data.m` |
| `P1_extract_and_center_data.m` | `preprocessing/run_extract_and_center.m` |
| `P1bis_order_slices.m` | `preprocessing/run_order_slices.m` |
| `P2_residual_correction_analysis.m` | `preprocessing/run_residual_correction.m` |
| `P2bis_nano_correction_analysis.m` | `preprocessing/run_nano_equalisation.m` |
| `P3_annotate_artifacts.m` | `preprocessing/run_annotate_artifacts.m` |
| `P4_register_to_atlas.m` | `registration/run_register_to_atlas.m`, the same five modes |
| `P4bis_add_sep_channel.m` | `registration/run_add_sep_channel.m` |
| `P5_collect_data_by_group.m` | `group_comparison/run_collect_by_group.m` |
| `P6bis_analyze_group_averages_and_normalize.m` | `group_comparison/run_normalise_groups.m` |
| `P7bis_analyze_group_differences.m` | `group_comparison/run_group_differences.m` |
| `P6_analyze_group_averages_and_normalize.m`, `P7_analyze_group_differences.m` | `archive/` (replaced by P6bis and P7bis before the refactor) |
| `P8_characterize_merged_distribution.m` | `adult_matlab/run_characterize_distribution.m`, until A1 and A4 replace it |
| `P9_compare_nano_vs_allen_ish.m` | `adult_matlab/run_compare_with_allen_ish.m`, until A1 to A3 replace it |
| `P10_compare_nano_vs_auto.m` | `adult_matlab/run_compare_nano_with_autofluorescence.m`, until A5 replaces it (A1 to A5: `ROADMAP.md`, section 3) |
| `v2_per_mouse.py`, `v2_to_ccf.py`, `v2_cohort.py` | `mapping/run_per_mouse.py`, `run_to_ccf.py`, `run_cohort.py`; the code in `mapping/sepmap/volumes/` |
| `v2_compare.py`, `v2_replot.py`, `v2_region_plot.py`, `v2_region_groups.py`, `v2_video.py`, `v2_video_compare.py` | `mapping/run_compare.py` and so on, each `run_` plus its old name; the code in `mapping/sepmap/young_vs_adult/` |
| `v2_inspect.py` | `mapping/run_closeup.py`; the code in `sepmap/young_vs_adult/closeup.py` |
| `v2_beyond_density.py`, `v2_beyond_controls.py`, `v2_beyond_figures.py`, `v2_beyond_regression.py`, `v2_sep_channel_check.py`, `v2_adult_arms.py` | `mapping/run_beyond_density.py` and so on; the code in `mapping/sepmap/adult/` |
| `v2_panel_build.py`, `v2_panel_fetch.py`, `v2_ish_regions.py`, `v2_ish_reliability.py`, `v2_ish_compare.py`, `v2_ish_words.py`, `v2_ish_roles.py`, `v2_ish_panel_test.py`, `v2_ish_arms.py` | `mapping/run_panel_build.py` and so on; the code in `mapping/sepmap/ish/` |
| `v2_diagnostics.py` | `mapping/run_diagnostics.py`; the code in `sepmap/diagnostics.py` |
| `v2_paths.py` | `mapping/sepmap/config.py` |
| `V2_ISH_PANEL`, `V2_ISH_TABLE` (environment) | `--panel targets` or `--panel ontology` of `run_ish_regions.py` and `run_ish_reliability.py`; the old variables are refused |
| `auto_annotation/`, `setup_auto_annotation.ps1` | `registration/auto_annotation/`, its `setup.ps1` |
| `landmark_refine*` (the GUI's `r` key) | `archive/`, retired |
| `LightSuite-main/`, `matlab_elastix-master/`, `yamlmatlab/`, `BioformatsImage/` | `third_party/LightSuite/`, `matlab_elastix/`, `yamlmatlab/`, `BioformatsImage/` |

Ages (P20), atlas keys (`demba_p20`), cohort tags (`young_P20`) and data file
names (`nano_4d_P20.mat`) kept their P-numbers: there P stands for the
postnatal day. Output folders and file names under the data root did not
change, apart from the logs of the detached runner, named after the new
scripts, and `processing_diagnostics\README.md`, which names them too.
