# landmark_refine

Automatic help for the control-point GUI (`matchControlPointsInSlices`). Two operations,
one piece of machinery, nothing trained:

| operation | what it does | measured (MG903 P20, 19 slices) |
|---|---|---|
| **refine** | *adjusted carry-over*: keep the previous slice's atlas landmarks, nudge them onto the nearest structure boundary of the new plane, and re-derive the histology coordinates through the image match field | **9.7 px** median vs **15.2** for a plain copy; beats it on 16/18 slices; keeps the annotator's mirror-pair structure |
| **suggest** | points from scratch: multi-scale matches scored by confidence + boundary + cornerness, chosen as left/right pairs across lateral bands | 10.2 px at 10 points, 9.3 at 15 unpaired |

The floor — the annotator's own affine residual — is 7.4 px. Every number and the full
experimental record are in `D:\sep_histology\sandbox_landmark_matching\README.md`.

Matching is LoFTR via `kornia` (pretrained `outdoor`), run at four blur levels and pooled.
Runs on CPU in ~2.5 s per slice; a cold Python process adds ~1.5 s.

## Layout

```
landmark_refine/
  core.py            the algorithms; importable, no file I/O
  cli.py             python cli.py request.mat response.mat   (what MATLAB calls)
  requirements.txt
  README.md
../landmark_refine.m          MATLAB wrapper: writes the request, calls cli.py, reads the response
../setup_landmark_refine.ps1  creates the venv and downloads the weights
```

## Setup (once per machine)

```powershell
cd D:\sep_histology\code
.\setup_landmark_refine.ps1
```

This creates `landmark_refine\.venv` (torch CPU, kornia, scipy, numpy — about 300 MB),
downloads the LoFTR weights into `~/.cache/torch/hub/checkpoints/`, and runs a self-test.
The wrapper also accepts `LANDMARK_REFINE_PYTHON` in the environment to point at any other
interpreter that has the requirements installed.

## From MATLAB

```matlab
% images as the GUI shows them: uint8 HxW, atlas plane and histology slice
% labels: the annotation plane, squeeze(gui_data.av(plane,:,:)); pass [] to skip snapping
out = landmark_refine('refine', atlas_im, hist_im, atlas_xy, labels);
%   out.atlas_pts  Nx2 [x y]   nudged onto boundaries
%   out.hist_pts   Nx2 [x y]   transferred through the match field
%   out.n_local    N           matches behind each transfer (0 = global-affine fallback)
%   out.local_rms  N           agreement of those matches, px

out = landmark_refine('suggest', atlas_im, hist_im, [], labels, struct('n', 10, 'mirror', true));
```

`out.ok` is false and `out.message` says why if anything went wrong; the wrapper never
throws on a matching failure, so the GUI can fall back to a plain copy.

## The MATLAB/Python bridge

Deliberately the simplest thing that works: a `.mat` file each way and one `system()`
call. Measured cost of that choice is ~1.5 s per call for process start-up and torch
import, against ~2.5 s of actual matching, so a persistent worker was not worth its
moving parts. If that ever changes, `core.py` is import-safe and stateless apart from the
cached model, so a worker is a small addition.

## Provenance

Built 2–3 September 2026 from the sandbox experiments. Things that were tried and
rejected are recorded there, with numbers: hard contrast gating (collapses to the
midline), snapping to image gradients (no gain — LoFTR already refines at full
resolution), scale-stability filtering (starves the selection), CLAHE, gradient-magnitude
input, cycle consistency, robust fitting, upsampling.
