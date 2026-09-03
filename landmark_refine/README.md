# landmark_refine

Automatic help for the control-point GUI (`matchControlPointsInSlices`). Two operations,
one piece of machinery, nothing trained:

| operation | what it does | measured (MG903 P20, 19 slices) |
|---|---|---|
| **refine** | *adjusted carry-over*: keep the previous slice's atlas landmarks, nudge them onto the nearest structure boundary of the new plane, and re-derive the histology coordinates through the image match field | **9.7 px** median vs **15.2** for a plain copy; beats it on 16/18 slices; keeps the annotator's mirror-pair structure |
| **suggest** | points from scratch: multi-scale matches scored by confidence + boundary + cornerness, chosen as left/right pairs across lateral bands | 10.2 px at 10 points, 9.3 at 15 unpaired |

The floor — the annotator's own affine residual — is 7.4 px. Every number and the full
experimental record are in `D:\sep_histology\sandbox_landmark_matching\README.md`.

Matching is LoFTR via `kornia` (pretrained `outdoor`), run at several blur levels and pooled:
four for `suggest`, two for `refine` (measured as accurate as four on the carry-over, at half
the cost).

## Speed

Measured on slice 6 of MG903, 400 x 570 px:

| | CPU | GPU (RTX A2000) |
|---|---|---|
| `refine`, 2 blur levels | 1.3 s | **0.4 s** |
| `suggest`, 4 blur levels | 2.4 s | **0.7 s** |
| torch import + model load, once per process | 2.2 s | 2.1 s |

So without anything else a GUI press costs several seconds, nearly all of it Python
start-up. Hence the **worker**: `serve.py` imports torch and loads the model once, then
answers requests dropped into a folder under `tempdir`. Measured from MATLAB, through the
worker: **0.47 s per `refine`, 0.76 s per `suggest`** on the GPU, identical results.

```matlab
landmark_refine_worker('start')     % the GUI launcher does this; idempotent, ~3 s
landmark_refine_worker('status')
landmark_refine_worker('stop')      % when done for the day; it idles at no CPU otherwise
```

`landmark_refine()` uses the worker when its heartbeat is fresh and falls back to
`cli.py` (same code, one process per call, ~4 s) when it is not — so nothing ever depends
on the worker being up. The protocol is in `serve.py`'s docstring: requests are written
under a temporary name and renamed into place, so neither side ever reads a half-written
file. The setup script installs the CUDA build of torch when `nvidia-smi` finds a GPU.

## Layout

```
landmark_refine/
  core.py            the algorithms; importable, no file I/O
  handler.py         one request dict -> one response dict; shared by the two below
  cli.py             python cli.py request.mat response.mat   (one process per call)
  serve.py           python serve.py <folder>                  (persistent worker)
  requirements.txt
  README.md
../landmark_refine.m          MATLAB wrapper: worker if alive, else cli.py; same answer
../landmark_refine_worker.m   start / stop / status of the worker from MATLAB
../setup_landmark_refine.ps1  creates the venv (CUDA torch if a GPU) and downloads the weights
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
call. The MATLAB side of that costs nothing measurable (bare `system()` 0.1 s, the `.mat`
round trip 0.03 s); what a call pays for is Python start-up, torch import and model load,
about 2 s. See Speed above for when a persistent worker would be worth having.

## Provenance

Built 2–3 September 2026 from the sandbox experiments. Things that were tried and
rejected are recorded there, with numbers: hard contrast gating (collapses to the
midline), snapping to image gradients (no gain — LoFTR already refines at full
resolution), scale-stability filtering (starves the selection), CLAHE, gradient-magnitude
input, cycle consistency, robust fitting, upsampling.
