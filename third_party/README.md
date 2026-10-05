# Third-party code

Four packages are vendored here, committed as they came so that the pipeline
runs without downloads. `sep_setup_paths` puts all four on the MATLAB path,
each with its subfolders, below our own folders. Only LightSuite has been
changed, and only in five files (`LightSuite/PATCHES.md`).

| folder | what it is | used for | licence |
|---|---|---|---|
| `LightSuite/` | LightSuite, Dimokratis Karamanlis ([github.com/dimokaramanlis/LightSuite](https://github.com/dimokaramanlis/LightSuite)), upstream commit `2f16206` (`LightSuite/UPSTREAM`) | extraction, slice ordering, alignment, the control-point GUI, registration to the atlas | GPL-3.0 (`LightSuite/LICENSE`) |
| `matlab_elastix/` | MelastiX, MATLAB wrappers for the elastix and transformix programs | called by LightSuite for the B-spline registration and to apply saved transforms | LGPL-3.0 (`matlab_elastix/LICENSE`) |
| `yamlmatlab/` | yamlmatlab, YAML reading and writing (CTU in Prague and Energocentrum PLUS) | matlab_elastix's YAML parameter files (`elastixYAML2struct`) | MIT (`yamlmatlab/MIT-license.txt`); it bundles the Java library SnakeYAML 1.9 (`yamlmatlab/+yaml/external/snakeyaml-1.9.jar`), Apache-2.0 (`+yaml/external/LICENSE-2.0.txt`) |
| `BioformatsImage/` | the `BioformatsImage` class, a reader built on Bio-Formats, with Bio-Formats' own MATLAB toolbox (`bfmatlab/`, with `bioformats_package.jar`) from the Open Microscopy Environment | reading the raw `.czi` files (extraction) | Bio-Formats: GPL-2.0-or-later (the notice in the `bfmatlab/*.m` files); the `BioformatsImage` class has no licence file in this copy |

Inside LightSuite, a few helper files come from elsewhere and keep their own
notices: `helpers/saveastiff.m` (copyright 2012 YoonOh Tak, a BSD-style
licence in the file), `helpers/print2array.m` (copyright Oliver Woodford
2008-2014 and Yair Altman 2015-, from export_fig), and
`helpers/fastSavePNG/` (copyright 2014-2015 Stefan Slonevskiy, MIT,
`helpers/fastSavePNG/LICENSE.md`).

Our changes to LightSuite are covered by its GPL-3.0: each modified file says
so in a comment, and `LightSuite/PATCHES.md` lists every change with its date.
The code of this project has no licence yet (decision L6 of
`docs/history/REFACTOR_PLAN.md`); the GPL obligations come back when the code is
released.

## Do not run

These files are kept as they came, and must not be run from this project:

- `LightSuite/compare_mice.m` and `LightSuite/scripts/protocol_manuscript_generate_plots.m`
  save figures into folders on `S:` (the lab share, whose raw data is
  read-only).
- `LightSuite/ls_analyze_slice_volume.m`, the upstream demo, has paths written
  for one machine (our copy points it at `D:\sep_histology`) and walks a whole
  brain through LightSuite on its own. The pipeline's drivers call LightSuite's
  functions themselves.

Also: `extractAxioscanImages.m`, which sat in `BioformatsImage/` but is not
part of Bio-Formats, writes folders next to the `.czi` files it reads, which
must never happen on raw folders. It has moved to `archive/`, and
`explore_czi_G.m`, the other stray, to `preprocessing/explore_czi.m`.

## Two duplicates to know about

- `matlab_elastix/MelastiX_examples/` holds two example scripts both named
  `RUN_ALL.m`, and both are on the path. Nothing calls either;
  `tests/sep_test_path.m` lists them as a known duplicate.
- A second copy of matlab_elastix nested inside itself
  (`matlab_elastix/matlab_elastix-master/`) was deleted on 4 October 2026;
  `.gitignore` still ignores it, since on the path it would shadow
  `transformix.m` and the rest.
