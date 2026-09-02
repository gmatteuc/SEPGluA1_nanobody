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

So a carry-over press in the GUI costs about **2.5 s on the GPU**, nearly all of it Python
start-up. The setup script installs the CUDA build of torch when `nvidia-smi` finds a GPU.
If that 2 s ever matters, the fix is a persistent worker process that imports once and
answers requests; `core.py` is stateless apart from the cached model, so it is a small
addition. It was not built because with the CPU numbers it would only have saved 1.5 of ~4.5 s.

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
call. The MATLAB side of that costs nothing measurable (bare `system()` 0.1 s, the `.mat`
round trip 0.03 s); what a call pays for is Python start-up, torch import and model load,
about 2 s. See Speed above for when a persistent worker would be worth having.

## Provenance

Built 2–3 September 2026 from the sandbox experiments. Things that were tried and
rejected are recorded there, with numbers: hard contrast gating (collapses to the
midline), snapping to image gradients (no gain — LoFTR already refines at full
resolution), scale-stability filtering (starves the selection), CLAHE, gradient-magnitude
input, cycle consistency, robust fitting, upsampling.
