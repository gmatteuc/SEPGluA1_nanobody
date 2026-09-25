# sep_histology

Whole-brain mapping of SEP-GluA1 (nanobody-labelled AMPA receptor subunit) surface
expression in coronal mouse brain sections, registered to the Allen CCF.

**This repository contains code only.** No data, at any stage.

## Data locations (not in this repo)

| What | Where | Access |
|------|-------|--------|
| Raw acquisition (`.czi`) | `S:\ElboustaniLab\#SHARE\Data\<mouse>\` | **READ-ONLY — never write or modify** |
| Derived / processed | `D:\sep_histology\data\` | read-write |
| Allen atlas | `D:\sep_histology\data\atlas\` | read-only in practice |

> **Raw data safety.** Nothing in this pipeline may write to `S:`. Raw microscopy data
> is irreplaceable — reacquisition means re-perfusing and re-sectioning animals.
> Note that some *vendored* scripts (`LightSuite-main/compare_mice.m`,
> `LightSuite-main/scripts/protocol_manuscript_generate_plots.m`,
> `BioformatsImage/extractAxioscanImages.m`) contain hard-coded save paths into `S:`.
> They are third-party and **must not be run as-is**.

## Pipeline

Scripts run in order; each stage writes into `D:\sep_histology\data\`.

| Stage | Script | Purpose |
|-------|--------|---------|
| P1 | `P1_extract_and_center_data.m` | Extract slices from raw `.czi`, centre volumes |
| P2 | `P2_residual_correction_analysis.m` | Residual / tiling correction |
| P2bis | `P2bis_nano_correction_analysis.m` | Nano-channel correction |
| P3 | `P3_annotate_artifacts.m` | **Manual** artifact annotation (`ArtifactAnnotator.m`) |
| P4 | `P4_register_to_atlas.m` | Register to Allen CCF (elastix, affine + B-spline) |
| P4bis | `P4bis_add_sep_channel.m` | Carry the SEP (green) channel into registered space by **re-applying** the saved transforms — adds `volume_registered_sep\`, changes nothing that exists |
| P5 | `P5_collect_data_by_group.m` | Assemble per-group 4D volumes across mice |
| P6bis | `P6bis_analyze_group_averages_and_normalize.m` | Per-mouse equalisation + normalisation |
| P7bis | `P7bis_analyze_group_differences.m` | Group-difference analyses |
| P8 | `P8_characterize_merged_distribution.m` | Region-wise distribution over Allen ontology |
| P9 | `P9_compare_nano_vs_allen_ish.m` | Correlate against Allen ISH (100-gene panel) |
| P10 | `P10_compare_nano_vs_auto.m` | Nano vs autofluorescence control (paired) |

`P6bis` / `P7bis` / `P8` / `P9` take a `channel` parameter (`'nano'` | `'auto'`) at the
top of the script; the channel is rolled into output folder names so runs never collide.

## The v2 route (young vs adult)

P5–P8 were built for adults on one atlas, and reused across ages they answer the
wrong question (details in `data\comparisons_v2\README.md`). The cross-age
comparison runs on a separate chain of Python scripts, which reads the registered
volumes directly and writes only under `data\comparisons_v2\`:

| Script | Purpose |
|---|---|
| `v2_per_mouse.py` | per brain, on the atlas of **its own age**: tissue mask, background-subtracted nano, auto and SEP |
| `v2_to_ccf.py` | each young brain carried DeMBA → CCF at its own age; adults are placed, not warped |
| `v2_cohort.py` | per-voxel cohort mean, SD and n, in the adult CCF |
| `v2_compare.py` | young against adult: maps, the per-structure table |
| `v2_region_plot.py` | the statistics, per-mouse region means with **no warping anywhere** |
| `v2_region_groups.py` | the same by system (primary vs higher sensory, frontal…) and by cortical layer |
| `v2_video.py`, `v2_video_compare.py` | plane-by-plane videos, per cohort and young beside adult |
| `v2_inspect.py` | one reading looked at closely: a coronal plane, its video and the cortical flatmaps (whole depth and by layer), all from one set of volumes. Runs in `tools\venv_flat` — see its docstring |
| `v2_diagnostics.py` | the sheets that make each step checkable by eye |

Five readings run through all of it, and none of them replaces another: `ratio`
(nano per unit autofluorescence), `sepratio` (nano per unit SEP, i.e. surface
receptor per unit receptor expressed), `cref` and `subref` (relative to the
brain's own isocortex / subcortex) and `zref` (range-matched).

Run them with the project venv:
`tools\venv_atlas\Scripts\python.exe v2_per_mouse.py`

### Important caveat on what the pipeline measures

`P6bis` equalises **per mouse** and `P8` z-scores **within brain**. Both deliberately
destroy absolute scale. Results describe the **relative spatial distribution** of
GluA1, not absolute expression level. Any claim of the form "GluA1 increases/decreases"
is *not* supported by this pipeline as configured.

## Vendored dependencies

Committed in-tree for reproducibility rather than pinned as submodules:

- `matlab_elastix-master/` — elastix/transformix MATLAB wrapper (registration)
- `LightSuite-main/` — slice handling, Allen CCF helpers
- `yamlmatlab/` — YAML parsing (MIT; bundles snakeyaml, Apache-2.0)
- `BioformatsImage/` — Bio-Formats reader for `.czi`

Each retains its upstream `LICENSE`.

## Requirements

MATLAB (Image Processing + Statistics toolboxes), plus `elastix`/`transformix` binaries
available to the elastix wrapper. Java is required for the Bio-Formats reader.
