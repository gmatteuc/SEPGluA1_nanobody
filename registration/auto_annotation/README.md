# auto_annotation

Automatic control points for the LightSuite GUI (`matchControlPointsInSlices`):
the human sets the atlas plane on a few anchor slices, this proposes every
slice, the human reviews the proposal in the GUI. LightSuite's registration is
untouched: the result is an ordinary `atlas2histology_tform.mat`, written by
the GUI from what was accepted.

Validated on the 24 hand-annotated brains with models that never saw the brain
scored, through LightSuite's own registration and by rerunning the paper's
analysis (cohort agreement equal to the manual annotation on 14 of 15 brains;
conclusions unchanged). The method, every experiment and the validation are in
the separate repository `SEPGluA1_autoannotation` (its `LOG.md`).

## Workflow (registration/run_register_to_atlas.m)

1. `run_mode = 'annotate'`: on the suggested anchor slices (`j` jumps between
   them) scroll to the right plane and press `a`; save with `s`. This writes
   `plane_anchors.mat` (the anchors and the plane the GUI interpolates for every
   slice) and `auto_atlas_planes.mat` (the warped atlas exactly as drawn).
2. `run_mode = 'autoannotate'`: writes `auto_proposal_controlpoints.mat` and
   `auto_proposal_info.mat`. About 30 s per brain on a GPU.
3. `run_mode = 'annotate'`: the proposal loads orange (provisional), the least
   confident quarter of each slice's points marked `?`. `k` accepts a slice and
   moves on, `u` re-proposes the slice at the plane on screen, the usual tools
   fix points. Only accepted or touched slices are saved; `annotation_provenance.mat`
   records how many saved pairs are the proposal's, unchanged.
4. `run_mode = 'register'` as always.

No file written here ends in `tform.mat`, which `registerSlicesToAtlas` globs for.

The keys `a`, `j`, `k`, `K`, `u`, `U`, the `?` flags and the files written on
save are a plugin of the GUI (`registration/annotation_gui/`), which
`annotate` passes in only when this engine is installed
(`auto_annotate('check')`); without it the GUI opens without them. The plugin
and the hook it uses are described in its help and in
`third_party/LightSuite/PATCHES.md` (LS7).

## What it does, per slice

1. image-only registration onto the slice's plane: NGF affine, then a smooth
   deformation (12 x 17 grid)
2. a small U-Net reads the atlas plane and picks 32 landmarks where a human
   would click
3. a small ConvNet corrects each landmark's partner point locally (96 px
   patches, averaged over 9 shifted crops); the spread of those 9 answers is the
   confidence

## Layout

```
registration/auto_annotation/
  core.py            the algorithms; importable, no file I/O
  cli.py             python cli.py propose <folder>
                       | section <folder> <slice> <plane> <out.mat>   (the u key)
                       | sections <folder> <request.mat> <out.mat>    (the U key)
  weights/           landmark.pt, matcher.pt, VERSION.txt (which models, which brains they never saw)
  requirements.txt
  setup.ps1          creates the venv (CUDA torch if a GPU), self-test
registration/pipeline/auto_annotate.m   MATLAB wrapper; out.ok / out.message on any failure
registration/annotation_gui/
  auto_annotation_plugin.m   the GUI's automatic-annotation keys, flags and files (opts.plugin)
  annotation_settings.m      suggested anchor count, the outlier rule (the GUI's * and register's check)
```

## Setup (once per machine)

```
cd D:\sep_histology\code
.\registration\auto_annotation\setup.ps1
```
