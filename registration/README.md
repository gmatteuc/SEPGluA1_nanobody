# Registration

Registers each brain's sections to the atlas of its age with LightSuite
(elastix, affine and B-spline), driven by the DAPI channel and by control
points between each section and the atlas. The adults go to the Allen CCF,
the young brains to the DeMBA atlas of their own age. Control points are
placed by hand, or proposed by the automatic annotation and reviewed by hand.

A new brain end to end, with the settings to edit and the checks:
[`../docs/ADDING_DATA.md`](../docs/ADDING_DATA.md).

## Run order

`run_register_to_atlas` runs one mode at a time (`run_mode`), on the mice of
`mice_to_process` and the atlas of `atlas_key`:

| step | mode | what it does |
|---|---|---|
| 1 | `'align'` | bridge the corrected volumes of `../preprocessing/` into LightSuite (DAPI, NANO, AUTO, DIFF, MASK), apply the slice order, align the slices, fit the atlas rigidly |
| 2 | `'angle'` | **manual**, optional: the cutting angle by eye (every adult had it); refused once the mouse has control points |
| 3 | `'annotate'` | **manual**, one mouse: the control-point GUI. With the automatic annotation installed: set the plane on the suggested anchor slices, then (after step 4) review the proposal |
| 4 | `'autoannotate'` | propose control points on every slice from the anchors (the engine in `auto_annotation/`, about 30 s per brain on a GPU) |
| 5 | `'register'` | check the control points, elastix refinement of every slice, the registered volumes of the five channels |

Then `run_add_sep_channel` carries the SEP (green) channel through the same
registration: it re-applies the saved transforms, nothing is refitted (30 to
45 min per brain, run detached). `P4BIS_MICE` in the environment replaces its
mouse list, so the cohort can be split across MATLAB sessions.

By hand, when needed:

| script | what it does |
|---|---|
| `remap_control_points` | carry saved control points across a change of slice order (a dry run unless told to write) |
| `qc/check_registration_error` | the atlas fit error of every aligned brain (`errall` in `regopts.mat`), to read against brains on the same atlas with a similar section count |
| `qc/verify_demba_setup` | before annotating a P20 brain: the atlas on the path, its label space and crop, the brains' local settings, the LightSuite fixes |

## Outputs

Per brain, in `<data>\<group>\<mouse>\lightsuite\`:

- 'align': `volume_centered_processed\`, `regopts.mat`, `volume_aligned\`,
  `volume_for_inspection.tiff`, `<mouse>_dim<d>_initial_registration.png`
- 'angle': `cutting_angle_data.mat`
- 'annotate': `atlas2histology_tform.mat` (the control points); with the
  automatic annotation also `plane_anchors.mat`, `auto_atlas_planes.mat`,
  `annotation_provenance.mat`
- 'autoannotate': `auto_proposal_controlpoints.mat`, `auto_proposal_info.mat`
  (never `atlas2histology_tform.mat`: the review saves what it accepts)
- 'register': `transform_params.mat`, `elastix_forward\`, `elastix_reverse\`,
  `volume_registered\chan01_DAPI.tiff`, `chan02_NANO`, `chan03_AUTO`,
  `chan04_DIFF`, `chan05_MASK`
- `run_add_sep_channel`: `volume_aligned_sep\`, `volume_registered_sep\`
  (`chan01_DAPI`, the check, and `chan02_SEP`); the check figure and text in
  `<data>\comparisons_v2\processing_diagnostics\sep_channel\`

## Where the code is

- `pipeline/`: `register_to_atlas` (the five modes, with the checks before
  each), `registered_grid_um` (the 10 um grid of every registered volume),
  `add_sep_channel`, `auto_annotate` (the MATLAB side of the automatic
  annotation), `recompute_backvalues`
- `annotation_gui/`: the automatic annotation's part of LightSuite's
  control-point GUI, a plugin (`auto_annotation_plugin`), and its settings
  (`annotation_settings`: the suggested anchor count, the outlier rule shared
  by the GUI's `*` mark and the check before 'register')
- `auto_annotation/`: the Python engine, its weights, its setup script; see
  its [README](auto_annotation/README.md)
- `qc/`: the checks above
- the atlases: `../atlas/` (`get_atlas` puts the one a run needs on the path)
- LightSuite, with our changes listed in
  `../third_party/LightSuite/PATCHES.md`

## Notes

- The age check: the driver stops when a mouse's age is not the atlas's (an
  adult counts as P56) and when a `local_settings.txt` of the mouse exists
  and gives a `px_atlas` other than the atlas's resolution: 20 um for DeMBA,
  10 um for the CCF. 'align' itself takes `px_atlas` and the AP crop from the
  atlas.
- The grid: 'register' writes every brain on the same 10 um-sampled grid,
  twice the 20 um registration grid, which the Python route expects. It sets
  the grid itself (`registered_grid_um`), whatever `px_atlas` the brain's
  `sliceinfo.mat` holds, and so does `run_add_sep_channel`. Before
  registering any brain it stops when a `px_register` in `sliceinfo.mat` is
  not 20 (LightSuite places the slices on a grid of half `px_register`), and
  prints a line for a `px_atlas` there other than 10.
- 'align' refuses a mouse with control points, anchors or a proposal: a new
  alignment rewrites the atlas block they are counted in. Control points are
  stored by slice position, so a change of order after annotating needs
  `remap_control_points` and the decisions file as it was before.
- `atlas_extent_slices` sets the atlas margin around the stack for a new
  mouse; a mouse with control points keeps the margin it was annotated with.
- 'register' stops on a mouse without control points unless
  `allow_image_only_registration` is set, which is for diagnostics only. A
  slice with fewer than 5 points, or left orange (provisional) at save,
  registers from the images alone.
- The automatic annotation's matcher is not repeatable on the GPU: two runs
  give points that differ slightly. The proposals are reviewed by hand.
- elastix 5.1.0 must be on the system path (the README of the repository).
