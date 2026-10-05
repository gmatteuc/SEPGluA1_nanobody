# Hand check, step 6: the annotation GUI and the two preprocessing GUIs

Kept as the record of the hand check of 2 October 2026, which passed
([REFACTOR_PLAN.md](REFACTOR_PLAN.md), Progress). The paths are those of the
check trees on G: at the time and may have been cleared since; since the
merge, the engine's environment is `registration\auto_annotation\.venv`.

For Giulio, on copies only: nothing here reads or writes `D:\sep_histology`.
About 40 minutes, in a MATLAB desktop (the GUIs need a display). Everything
that can run headless has already been driven, old code against new; see
"Already checked by script" at the end.

The code is `G:\sep_refactor\check\code` (branch `refactor`, which now holds
step 6). The data are:

- `G:\sep_refactor\gui_check\hand\data`: a copy of MG914's lightsuite files
  as reviewed on 30 Sep, and the P28 atlas. The helper `hc` refills it from a
  read-only source (`gui_check\src`) and keeps a snapshot before each part;
- `G:\sep_refactor\pre\gui\data`: MG914's slice-order and artifact inputs,
  for the two preprocessing GUIs.

## Setup (once per MATLAB session)

Close any other GUI of this project, then in a fresh MATLAB:

```matlab
restoredefaultpath
setenv('SEP_DATA_ROOT', 'G:\sep_refactor\gui_check\hand\data')
setenv('AUTO_ANNOTATION_PYTHON', 'D:\sep_histology\code\auto_annotation\.venv\Scripts\python.exe')
cd('G:\sep_refactor\check\code')
sep_setup_paths
addpath('G:\sep_refactor\gui_check\scripts', '-end')   % hc, the copy helper
```

Check before starting:

```matlab
p = get_paths(); p.data             % G:\sep_refactor\gui_check\hand\data
which matchControlPointsInSlices    % G:\sep_refactor\check\code\third_party\LightSuite\control_point_gui\...
which auto_annotation_plugin        % G:\sep_refactor\check\code\registration\annotation_gui\...
auto_annotate('check')              % ok: 1
!git -c safe.directory=G:/sep_refactor/check/code -C G:/sep_refactor/check/code log -1 --oneline
```

## The helper `hc` (backup and compare)

| type | does |
|---|---|
| `hc('restore', 'reviewed')` | empties the copy and refills it with MG914 as reviewed (every file); the snapshot of it goes to `gui_check\hand\before` (the backup) |
| `hc('restore', 'review')` | the same without `atlas2histology_tform.mat` and `annotation_provenance.mat`: the proposal opens orange on every slice |
| `hc('restore', 'angle')` | the same as `'review'`, for the cutting-angle GUI, which refuses a brain with points |
| `hc('annotate')`, `hc('angle')` | `run_register_to_atlas` with `run_mode` `'annotate'` or `'angle'` on MG914 (atlas `'demba_p28'`, `atlas_extent_slices = 15`), typed in `hc` so the driver in the check tree is not edited |
| `hc('compare')` | the copy against its snapshot, file by file (`tools\sep_compare_outputs`), then variable by variable for each `.mat` that differs |
| `hc('list')` | the copy's files, dates and sizes |

The source in `gui_check\src` is never written, so any part can be redone
from `hc('restore', ...)`.

## A. Annotate a reviewed brain, nothing edited (`r`, `s`, close)

```matlab
hc('restore', 'reviewed')
hc('annotate')
```

- Console: `NOTE: control points already exist ...`, `cutting angle: from
  cutting_angle_data.mat`, **no** line `the automatic annotation is not
  installed`, then `Loaded 7 plane anchor(s) from ...\hand\data\...` and
  `Loaded the automatic proposal onto 0 slice(s), provisional (orange).`
- Window title: `MG914_SepGluA_P28 >> REVIEW DONE: s to save, then
  run_register_to_atlas 'register' << [ h ] controls`. Slice 1: `Slice 1/46
  ACCEPTED, Npoints = 31, mse = 5.67`, `anchored at atlas 10.00 h-slice
  widths`, solid markers.
- Controls window: the right column starts with `AUTOMATIC ANNOTATION (two
  GUI sessions)`, then `TAKE FROM THE NEIGHBOUR [ t ]` and `CARRY FORWARD
  [ p ]`. No `r` anywhere.
- Right arrow through a few slices: yellow `?` after some numbers. Slice 4:
  `13*` in cyan and `-- 1 pair(s) * far off this slice's affine` (stars on
  slices 4, 20, 22 to 25, 39 to 46, as the drive found).
- Press `r`, then `x`: nothing happens, no message.
- Back on slice 1, press `s`. Console:
  `Saved 46 annotated slice(s) to ...`, `Saved 7 plane anchor(s) to ...`,
  `Provenance: 1442 of 1453 saved pair(s) are the automatic proposal,
  unchanged.`
- Close the window, `Save?` Yes: the same three lines, the window goes.

```matlab
hc('compare')
```

Expected: 10 files, 8 same. `annotation_provenance.mat` differs only in
`saved` (the time). `atlas2histology_tform.mat`: `atlas_control_points` and
`histology_control_points` same; `atlas2histology_tform` differs. That
variable is the affine of the slice on screen when the window was closed
(registration reads only the points); closed on slice 1 it differs from the
30 Sep file by the same amount in the old code. `plane_anchors.mat`,
`auto_atlas_planes.mat` (not rewritten) and the rest: same.

## B. Review a proposal (`k`, `u`, `j`, `a`, `U`, `K`, `t`, `p`, `s`)

```matlab
hc('restore', 'review')
hc('annotate')
```

- Console: `Loaded 7 plane anchor(s) ...`, `Loaded the automatic proposal
  onto 46 slice(s), provisional (orange).`
- Window title: `>> STEP 3 REVIEW: 46 slice(s) left, k accepts, u
  re-proposes <<`. Slice 1: hollow orange markers, orange numbers, red `?`;
  title `Slice 1/46 PROPOSED (orange, not saved: k accepts), Npoints = 32,
  mse = 6.77`.

Then:

1. `k`: `Slice 1 accepted (32 point(s)); 45 provisional slice(s) left.`,
   slice 2 shown, title bar `45 slice(s) left`. Left arrow: slice 1 now
   `ACCEPTED`, solid markers, its `?` yellow. Right arrow back to slice 2.
2. Slice 2: wheel two notches. The second title line ends `-- plane changed,
   press u to re-propose here`. `u`: `Proposing slice 2 at atlas plane
   ...`, then `Proposed 32 point(s), n marked ?. k to accept.` (a few
   seconds, GPU); points orange at the new plane. `k` accepts it.
3. `j`: slice 16 (the suggested slices are 1, 16, 31, 46). Wheel to another
   plane, `a`: `Slice 16 anchored at atlas plane ... (8 anchor(s) set).`
   `a` again on the same plane: `Anchor removed from slice 16.`
4. `U` (shift+u): `Re-proposing 44 orange slice(s) from ... anchor(s) /
   accepted slice(s)...`, about half a minute, then `Done: 44 slice(s)
   re-proposed, still orange.` Slices 1 and 2 untouched.
5. `K` (shift+k): dialog `Accept all 44 orange slice(s) as proposed?`.
   Cancel: nothing changes. `K` again, `Accept all`: `Accepted 44 slice(s)
   as proposed: [...]. Save with s.`; title bar `REVIEW DONE: s to save,
   then run_register_to_atlas 'register'`.
6. `t` and `p`: on slice 10 press `c` (points gone, `not anchored,
   predicted ...`), then `t`: `Took 32 point(s) from slice 9 as they are, at
   atlas plane .... Provisional: touch one to keep, c to discard.` (orange).
   `e`, click one of its points and let go (solid again), `e`. `p`: `CARRY
   FORWARD on: ...` and `CARRY FORWARD` in the title bar; `p` again: off.
7. `s`, then close, `Save?` Yes: `Saved 46 annotated slice(s) ...`, `Saved 7
   plane anchor(s) ...`, `Provenance: X of Y saved pair(s) are the automatic
   proposal, unchanged.`

```matlab
hc('compare')
```

Expected: `atlas2histology_tform.mat` and `annotation_provenance.mat` only
in the copy (the review state had none); `plane_anchors.mat`:
`anchor_slices` and `anchor_planes` same (once the anchor on slice 16 is
removed again), `planes` may differ (the planes of the slices as now saved);
`auto_atlas_planes.mat`, `auto_proposal_*.mat` and the rest same (the GUI
does not write them). In `annotation_provenance.mat`, `n_auto_unchanged`
equals `n_points` on every slice except slice 10.

## C. No engine installed

```matlab
setenv('AUTO_ANNOTATION_PYTHON', '')
hc('restore', 'reviewed')
hc('annotate')
```

- Console: `the automatic annotation is not installed (no Python
  interpreter found: run registration\auto_annotation\setup.ps1 or set
  AUTO_ANNOTATION_PYTHON):` and `the GUI opens without its keys (a, j, k,
  u).`; no `Loaded ...` lines.
- The GUI as LightSuite with our fixes: window title `MG914_SepGluA_P28`
  (`[ h ] controls` appears after the first arrow key); slice 1 drawn without its
  points until the first arrow key (LightSuite has always opened that way;
  the plugin draws it at start-up); titles without `ACCEPTED`; no `?`; the
  `*` marks as in A; controls window without the `AUTOMATIC ANNOTATION`
  block.
- `a`, `j`, `k`, `u`, `r`: nothing. `s`: only `Saved 46 annotated slice(s)
  to ...` (no anchor, no provenance line). Close, `Save?` Yes.

```matlab
hc('compare')
setenv('AUTO_ANNOTATION_PYTHON', 'D:\sep_histology\code\auto_annotation\.venv\Scripts\python.exe')
```

Expected: as in A, but `annotation_provenance.mat` and `plane_anchors.mat`
same (not written).

## D. Cutting angle (`run_mode = 'angle'`)

```matlab
hc('restore', 'reviewed')
hc('angle')
```

Refused: `run_register_to_atlas: MG914_SepGluA_P28 already has control
points. The cutting angle changes the atlas block ...`. Nothing written.

```matlab
hc('restore', 'angle')
hc('angle')
```

- Console: `NOTE: a cutting angle is already saved and will be overwritten
  on close: ...\hand\data\...\cutting_angle_data.mat`, `opening the
  cutting-angle GUI against atlas 'demba_p28'.`, the keys line.
- The cutting-angle GUI (LightSuite, unchanged). Shift+arrows tilt, the
  wheel moves the plane, return saves the plane on the current slice: do 3
  slices. Close: the warning `Not all slices have a saved atlas position`,
  OK, `Save the alignment data?` Yes: `Data saved to:
  G:\sep_refactor\gui_check\hand\data\young\MG914_SepGluA_P28\lightsuite\cutting_angle_data.mat`.

```matlab
hc('compare')
```

Expected: only `cutting_angle_data.mat` differs (the new angle); every
other file same, `regopts.mat` included.

## E. Slice order, `run_mode = 'edit'` (never `'apply'`)

```matlab
setenv('SEP_DATA_ROOT', 'G:\sep_refactor\pre\gui\data')
p = get_paths(); p.data             % G:\sep_refactor\pre\gui\data
order_slices(struct('mice_to_process', {{'MG914_SepGluA_P28'}}, 'run_mode', 'edit'))
```

- Console: `=== MG914_SepGluA_P28 (group young) ===`, `NOTE: a decisions
  file already exists and will be overwritten on save:`, `opening
  SliceOrderEditor on: G:\sep_refactor\pre\gui\data\young\MG914_SepGluA_P28\lightsuite\volume_for_ordering.tiff`,
  `Loaded saved order from decisions file.`
- The close-up and the montage windows. Arrows move, `m` toggles the
  montage. Change nothing. `s`: `Processing decisions saved to: ...` and
  `Montage snapshot saved to: ...`. Esc, `Save & Close`. MATLAB returns.

```matlab
d = 'G:\sep_refactor\pre\gui\data\young\MG914_SepGluA_P28\lightsuite\volume_for_ordering_processing_decisions.txt';
isequal(fileread(d), fileread('G:\sep_refactor\pre\gui\decisions_before.txt'))   % 1
```

(If 0, `visdiff` the two.) Never `'apply'`: it rebuilds `volume_ordered.tiff`
in the folder `sliceinfo.mat` names, the production one (bug list); this
tree has no `sliceinfo.mat` on purpose, so it would stop.

## F. Artifacts

```matlab
p = get_paths();
annotate_artifacts(struct('paths', p, 'groups_to_process', {{'young'}}, ...
    'mice_to_process', {{'MG914_SepGluA_P28'}}, 'correction_type', 'slicewise'))
```

- Console: `run_annotate_artifacts: 1 mouse/mice selected.` (loads 0.6 GB).
- Window: `Slice 1/46 - Red: Scaled Autofluo, Green: Nano | Keys: ...`.
  Right arrow to some slice z, `d`, click a few points, double-click to close
  the polygon; `s`: `Saved artifact masks for MG914_SepGluA_P28 in
  G:\sep_refactor\pre\gui\data\young\MG914_SepGluA_P28\lightsuite\correction_output\artifact_mask_volume_slicewise.mat`.
  Esc, `Save & Close`. MATLAB returns.
- The same command again: `Loaded existing annotations for
  MG914_SepGluA_P28 from ...`, and the polygon is back on slice z. Esc,
  `Discard & Close`.

```matlab
m = 'G:\sep_refactor\pre\gui\data\young\MG914_SepGluA_P28\lightsuite\correction_output\artifact_mask_volume_slicewise.mat';
S = load(m); find(squeeze(any(any(S.artifact_mask_vol, 1), 2)))'   % z only
delete(m)                                                            % the tree as before
```

## G. Finish

```matlab
setenv('SEP_DATA_ROOT', ''); setenv('AUTO_ANNOTATION_PYTHON', '')
!git -c safe.directory=G:/sep_refactor/check/code -C G:/sep_refactor/check/code status --short
```

The status prints nothing (the GUIs write only into the data copies). Quit
MATLAB.

## Already checked by script (1 Oct)

The sandbox's drive scripts (`drive_gui`, `drive_session1`,
`drive_provenance`, `drive_stars`) and four added ones (`drive_gui_reopen`,
`drive_hand`, `drive_save_unchanged`, and `drive_hand` and `drive_stars`
without the engine), copied to `gui_check\scripts` with their paths adapted
and the GUI opened through each code's own `annotate` mode, run in a fresh
MATLAB each: the code at `586430f` on `gui_check\old\data`, the new code on
`gui_check\new\data`, inputs identical by hash. Logs in `results\logs`, the
folders after each drive in `results\<drive>\<side>\files`. With the engine:
every file identical, apart from the time saved in the provenance, the time
a point was clicked, and the GPU proposals (the `u` on slice 2 moved its
points by a median of 0.19 px, two old runs by 0.13 px; the whole proposal
by a median of 0.09 px, at most 13 px, as between two old runs); the
messages identical apart from the old `annotate` mode's line about the
landmark_refine worker and the `r` hint after a carry-forward. Without the
engine: the same files; the new GUI is the plain one by design (no
suggested-anchor title, no stage in the title bar, the hint about grabbing
a point instead of `k`).
