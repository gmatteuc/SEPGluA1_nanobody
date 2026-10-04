# Preprocessing

From the raw `.czi` files of a brain to the corrected volumes the registration
reads. Each section is recorded in four channels: DAPI, the nanobody (Cy5,
called nano), autofluorescence (Cy3, called auto) and the green channel of
SEP (filter EGFP). The steps run one brain or a batch at a time; two are
manual. The mice come from the cohort table (`common/cohort.csv`, through
`get_cohort`).

A new brain end to end, with the settings to edit and the checks:
[`../docs/ADDING_DATA.md`](../docs/ADDING_DATA.md).

## Run order

| step | script | what it does |
|---|---|---|
| 1 | `run_copy_raw_data` | copy the `.czi` files from the lab share (read only) into `<data>\<group>\<mouse>\`, and check their sizes |
| 2 | `run_extract_and_center` | find the sections in the `.czi` files (LightSuite's `getSliceInfo`), write every channel centred, and the colour composite for the ordering |
| 3 | `run_order_slices` | **manual**: `'edit'` opens SliceOrderEditor to reorder, flip and discard slices; `'apply'` rebuilds the preview of the order |
| 4 | `run_residual_correction` | per slice, fit the autofluorescence to the nano on reference pixels; the scaled autofluorescence and the nano minus it |
| 5 | `run_nano_equalisation` | equalise the nano intensity across the slices of each brain, for a batch of brains run together |
| 6 | `run_annotate_artifacts` | **manual, optional**: outline artifacts slice by slice (ArtifactAnnotator) |

Then `../registration/run_register_to_atlas`, mode `'align'`, which reads the
outputs of steps 3 to 6.

Tools, run by hand:

| script | what it does |
|---|---|
| `explore_czi_G` | print the series, scenes, channels and tile positions of `.czi` files (Bio-Formats); only reads |
| `make_atlas_reference_sheet` | coronal plates of the atlas one section apart, to keep beside SliceOrderEditor |
| `make_ordering_volume` | rebuild the ordering composite with other colours, without extracting again; the order is untouched |

## Outputs

Per brain, in `<data>\<group>\<mouse>\lightsuite\`:

- `sliceinfo.mat`, `volume_centered\chan01_DAPI.tiff`, `chan02_Cy5.tiff`
  (nano), `chan03_Cy3.tiff` (auto), `chan04_EGFP.tiff` (SEP),
  `volume_for_ordering.tiff` (step 2)
- `volume_for_ordering_processing_decisions.txt` (columns `OriginalIndex`,
  `FlipState`, `NewOrderOriginalIndex`), `slice_order_montage.png`,
  `volume_ordered.tiff` (step 3)
- `correction_output\corrected_volume_<type>.mat`,
  `scaled_auto_volume_<type>.mat`, `scaled_difference_video_<type>.mp4` and
  `diagnostic_plots\` (step 4), `equalized_volume.mat` (step 5),
  `artifact_mask_volume_<type>.mat` (step 6)

Across brains: `<data>\intensity_diagnostics\`, the statistics, profiles,
heatmaps and videos of step 5, named with the time of the run.

Steps 2 to 6 write only inside the data root; step 1 reads the share and
refuses a destination on the share's drive or on a network path.

## Where the code is

Each driver sets its settings and calls one function in `pipeline/`
(`copy_raw_data`, `extract_and_center`, `order_slices`,
`residual_correction`, `nano_equalisation`, `annotate_artifacts`), with
`ArtifactAnnotator` and `select_reference_pixels` beside them. Shared:
`../common/` (`get_cohort`, `loadVolume`, `select_background_pixels`) and
LightSuite in `../third_party/`.

## Notes

- The extraction settings come from `local_settings.txt`, the mouse folder's
  first, then `lightsuite\`'s; without either, LightSuite's defaults, which
  are the adults' values. Its `px_atlas` ends up in `sliceinfo.mat` but no
  longer sets the grid of the registered volumes: 'register' puts every brain
  on the 10 um grid (`registered_grid_um`) and checks that `px_register` is
  20; see "The registration grid" in `ADDING_DATA.md`.
- The slice order is applied by the registration's 'align' mode, which reads
  the decisions file; steps 4 to 6 work in the extraction order.
- Step 4 writes two variants, `slicewise` and `global`; every later step
  reads `slicewise`.
- Step 5 stacks the brains of one run, padded to the largest, so a brain's
  result depends on the brains run with it. The young brains MG909 to MG914
  were run together. It rewrites `equalized_volume.mat` of every brain in the
  run without a check: never include a brain that is already aligned, and run
  a single new brain alone.
- Step 6 is not part of the current runs. Its masks go only into the MASK
  channel, which no current analysis reads; three RWS brains (MG691, MG692,
  MG693) have masks, from November 2025, and no other brain has.
  ArtifactAnnotator can break when Escape is pressed again while its close
  dialog is open.
- Do not extract a registered brain again: it rewrites `sliceinfo.mat` and
  the centred volumes, and the extraction subsamples points at random.
