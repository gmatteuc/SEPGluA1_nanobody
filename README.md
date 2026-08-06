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
| P5 | `P5_collect_data_by_group.m` | Assemble per-group 4D volumes across mice |
| P6bis | `P6bis_analyze_group_averages_and_normalize.m` | Per-mouse equalisation + normalisation |
| P7bis | `P7bis_analyze_group_differences.m` | Group-difference analyses |
| P8 | `P8_characterize_merged_distribution.m` | Region-wise distribution over Allen ontology |
| P9 | `P9_compare_nano_vs_allen_ish.m` | Correlate against Allen ISH (100-gene panel) |
| P10 | `P10_compare_nano_vs_auto.m` | Nano vs autofluorescence control (paired) |

`P6bis` / `P7bis` / `P8` / `P9` take a `channel` parameter (`'nano'` | `'auto'`) at the
top of the script; the channel is rolled into output folder names so runs never collide.

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
