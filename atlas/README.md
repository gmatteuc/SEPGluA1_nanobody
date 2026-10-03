# Atlases

Where each reference atlas is, which one a brain belongs to, and the checks
that measured how the young and adult atlases relate. The adults are
registered to the Allen CCFv3 (10 um), each young brain to the DeMBA atlas of
its own age (Carey 2025; 20 um). Both annotations are in the CCF's
parcellation index, so a region means the same in either.

## Functions

| file | what it does | used by |
|---|---|---|
| `get_atlas.m` | an atlas by key (`'ccf'`, `'demba_p<age>'`): its folder, files, resolution, age and AP crop; puts that folder on the path and takes every other atlas folder off | the preprocessing and registration drivers, the QC scripts |
| `cohort_atlas_key.m` | the atlas key of a group and age (`'ccf'` for the adults); an error for an age whose atlas is not built | `get_cohort_spec`, `verify_demba_setup` |
| `get_atlas_crop.m` | the annotation on the grid of the registered volumes of a cohort (900 AP planes for the adults, 994 for P20) | `run_collect_by_group`, `run_normalise_groups`, P8 |
| `get_allen_region_mask.m` | the voxels of named regions and all their descendants, through the ontology tables | the plasticity comparison, the QC scripts, P8, P9 |
| `build_demba_atlas.py` | builds `<data>\atlas_demba_p<age>\` from BrainGlobe's DeMBA at that age | by hand, once per age |

To build the atlas of a new age, from the code root:

```
tools\venv_atlas\Scripts\python.exe atlas\build_demba_atlas.py 28
```

It downloads the atlas on first use, remaps the annotation to the
parcellation index (keeping the original as
`annotation_structureids_original.nii.gz`), measures the AP crop that holds
the adults' anatomy and writes it to `aplims.txt`, with `source.txt`.
`get_atlas('demba_p28')` then finds the folder by its key. P16, P20, P22, P28,
P32 and P36 are built.

## QC, run by hand (`qc/`)

| script | what it does |
|---|---|
| `atlas_diagnostics` | measures both atlases (resolution, grid, crop, label space, the AP coverage of each brain) and writes `ATLAS_PARAMETERS.md` and `atlas_crop_and_ap_mapping.png`; after changing an atlas, a crop or the section thickness |
| `compare_atlas_regions` | per region, volume fraction, centroid and label overlap between DeMBA P20 and the CCF; the measurement behind registering the young brains to an age-matched atlas |
| `compare_atlases_montage` | the two atlases side by side, coronal, at matched levels |
| `check_demba_to_allen` | the DeMBA P20 annotation carried to the CCF by CCF Translator, against the CCF annotation |

They write to `<data>\young\registration_qc\`. `draw_overlap_matrix`,
`region_overlap_matrix` and `pad_to_canvas` are their helpers.

## Data

- `<data>\atlas\`: `average_template_10.nii.gz`, `annotation_10.nii.gz`,
  `annotation_boundary_10.nii.gz`, the ontology tables (`parcellation*.csv`)
  and `ccf_translator_fields\`, the copy of the DeMBA-to-CCF deformation
  fields that CCF Translator downloads into its own package folder
  (`tools/requirements_atlas.txt` says how to restore them).
- `<data>\atlas_demba_p<age>\`: `average_template_10.nii.gz`,
  `annotation_10.nii.gz`, `annotation_structureids_original.nii.gz`,
  `aplims.txt`, `source.txt`.

## Notes

- LightSuite finds the atlas with `which('average_template_10.nii.gz')`, so
  every atlas folder holds files of that name, and only one may be on the
  path. Never `addpath` an atlas folder: `get_atlas` does it.
- The DeMBA files hold 20 um data under the `_10` names that LightSuite
  needs. The real resolution is `res_um` in `get_atlas` and `px_atlas` in a
  young brain's `local_settings.txt`; both must say 20, or the AP scale of the
  registration halves without an error.
- BrainGlobe ships DeMBA in Allen structure ids, which collide numerically
  with the parcellation index without meaning the same region; the build
  remaps every id.
- DeMBA is about 11% longer in AP than the CCF for the same anatomy (the
  slope of the region-centroid regression, `get_atlas`'s help), so the adult
  crop does not map to the young atlas by dividing by two; the crop is
  measured.
- `check_demba_to_allen` reads a transformed annotation made by a script that
  is not in the repository (`tmp/demba_to_allen.py`); the numbers in
  `get_atlas`'s help came from two more such scripts, and
  `atlas_diagnostics` measures them again.
