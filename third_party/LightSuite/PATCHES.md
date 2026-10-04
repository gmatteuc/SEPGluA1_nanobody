# Our changes to LightSuite

LightSuite (Dimokratis Karamanlis, GPL-3.0, see `LICENSE`) is vendored from
upstream commit `2f16206` of 29 October 2025 (see `UPSTREAM`). Of its 153
files, five differ from that commit. Each carries the line

```
% Modified for the SEP-GluA1 project; see third_party/LightSuite/PATCHES.md.
```

at its top: in the four functions, just below the function's own help text
(a comment above the `function` line would replace that text in `help`), and
on the first line of the demo script. Everything else is upstream as
downloaded.

| file | changes | kind |
|---|---|---|
| `ls_analyze_slice_volume.m` | P0a (reverted) | demo script, back to upstream |
| `slice_module/alignSliceVolume.m` | P0b, LS2 | registration |
| `slice_module/registerSlicesToAtlas.m` | LS1 | registration |
| `control_point_gui/matchControlPointsInSlices.m` | LS3, LS4, LS5, LS7 (LS8 removed, AA moved out) | control-point GUI |
| `slice_module/SliceOrderEditor.m` | LS6 | slice-order GUI |

The history of each change is in this repository's git log: before the
refactor the folder was `LightSuite-main/`, so
`git log refactor-start -- LightSuite-main` lists every commit that touched it.
The two edits P0a and P0b predate the repository and were found by comparing
the first commit (`b9f2343`) with upstream `2f16206`.

## Changes that can affect results

Both registration fixes are still missing upstream (checked 30 Sep 2026).

**LS1. `registerSlicesToAtlas`: 0x4 placeholders when there are no control
points** (commit `00138a3`, 2 Sep 2026). When a brain has no control-point
file, the per-slice placeholders were `zeros(0,2)`, but the loop reads columns
`[3 2]` of them, so the first slice threw "Index in position 2 exceeds array
bounds" and the image-only branch could never run. Now `zeros(0,4)`. Effect on
results: none on any result in use. Every registered brain has control points,
so the branch ran only in the image-only trials of 2 Sep, which were redone.
(Two later commits, `4729628` and `e729787`, first changed the atlas-plane
extrapolation in this file and then put it back, so that the 17 adults stay
reproducible; LS1 is all that remains.)

**LS2. `alignSliceVolume`: atlas resolution from `px_atlas`** (`00138a3`,
2 Sep 2026). `regopts.allenres` was hard-coded to 10 um; it is now
`sliceinfo.px_atlas`. `registerSlicesToAtlas` divides by it, so a 20 um atlas
was downsampled to 40 um and every AP coordinate was off by a factor of two,
with nothing in the log. Effect on results: none for the adults
(`px_atlas = 10`); required for any brain registered to a 20 um DeMBA atlas.

**P0b. `alignSliceVolume`: PNG export of the initial registration** (before
the first commit). The three `<mouse>_dim<d>_initial_registration` figures
are written with `exportgraphics` (painters renderer, 300 dpi) instead of
upstream's `savepngFast`. Why: not recorded. Effect on results: none; only the
look of those diagnostic PNGs. Upstream has since changed this export its own
way, so this edit is dropped when the copy is replaced by upstream (Z1 below).

## GUI changes (no effect on the registration code)

**LS3. Control-point GUI fixes** (`cc75e01`, `9fd3140`, `5ae1fc1`, `c5f4936`;
2 to 29 Sep 2026). Clicks on a point reach the image handler (markers, rings
and grid lines are click-through), dragging works (the motion callbacks go on
the figure, not the image), the title of a slice with a partial pair says so,
and saving with no points and no annotation yet writes no empty
`atlas2histology_tform.mat` (an empty file would let `register` run from
images alone instead of stopping).

**LS4. Atlas-plane prediction without backward extrapolation, and the
order-conflict warning** (`5ea8aed`, `4729628`, `e729787`, `cc75e01`;
2 Sep 2026). Upstream extrapolated the plane of unannotated slices from the
last two annotated ones, so one step backwards sent every later slice off the
front of the atlas (in the case that showed it, 13 of 43 slices were clamped
to plane 1). The GUI
now interpolates between anchors and steps beyond them at the section spacing,
and a slice anchored behind its predecessor shows a one-line
`ORDER CONFLICT` banner. Registration keeps upstream's rule (see LS1). Effect:
none when every slice has control points, as every registered brain has.

**LS5. Editing and the review aids** (`5ea8aed`, `4729628`, `151b489`,
`5ae1fc1`, `2dcf116`, `7f42fdf`, `04c0484`; 2 to 30 Sep 2026). `e` edit mode
(grab and drag a point; `d` deletes it with its pair), ctrl+z takes back the
last point, `p` carries a slice's points forward and `t` copies the
neighbour's, both provisional and never saved unless touched; numbered labels
on both panels (since step 6 a plugin can add a mark after them, LS7: the
automatic annotation's `?` flags); `h` reopens the controls window; the slice
title gives the atlas plane and whether it is anchored; a `*` marks pairs far
off the slice's own affine, by the same rule the registration driver applies
before registering (`report_suspect_pairs` in
`registration/pipeline/register_to_atlas.m`). Since step 6 the rule's three
thresholds are `gui_data.outlier_rule` (3 times the slice's median residual
and 30 px for a pair, a 20 px median for the whole slice), which a plugin may
set: ours sets them from `annotation_settings`, the definition
`report_suspect_pairs` reads too. Effect: none on the registration code. One
consequence to know: a slice left orange (provisional) at save time is saved
empty and registers from images alone.

**LS6. SliceOrderEditor montage** (`606d935`, `8ab0873`, `733b4bc`, `8fee90c`;
12 and 13 Aug 2026). A second window with every slice as a tile, click to
jump, drag to move, `m` to toggle it, every shortcut in the title, and
`slice_order_montage.png` saved next to the decisions file. The new code is
one block at the bottom of the file plus seven one-line calls marked
`% montage`; `showMontage = false` restores the original editor. Effect: none.
(Our copy saves the decisions without upstream's later crop columns; see Z1.)

**LS7. A plugin hook in the control-point GUI** (step 6 of the refactor,
1 Oct 2026). `opts.plugin`, empty by default, is a function handle that adds
to the GUI without changing it. The GUI calls it as

```
value = plugin(event, gui_fig, gui_data, value, info)
```

with `gui_fig` the window, `gui_data` its state at that moment, `value` what
the GUI goes on with (the plugin returns it, changed or not; `[]` where
nothing is asked), and `info.gui` three of the GUI's functions:
`update_slice(gui_fig)`, `update_window_title(gui_fig)` and
`gui_data = step_slice(gui_data, step)`. A plugin returns `value` unchanged
for an event it does not know. The calls, eleven in all, each one or two
lines through the helper `call_plugin`:

| event | where | value, info |
|---|---|---|
| `open` | start-up, once the window is built and its state stored, before the first draw | `[]`. The plugin may add fields of its own to `gui_data`, fill slices with points and their provisional flags, set the three fields below, store it with `guidata`, and draw |
| `key` | `keypress`, any key the GUI does not use (`otherwise`) | `[]`; `info.key`, the key event |
| `planes` | `update_slice`, before the plane prediction | the atlas plane of each slice with points, NaN elsewhere; the plugin may fix more |
| `title` | `align_ccf_to_histology`, before the fit's numbers and the `*` warning are added | the two lines of the slice title |
| `labels` | `draw_point_numbers` | `rows` (one flag per point, all false), `text`, `color`: the numbers of flagged points get `text` after them, in `color` |
| `window` | `update_window_title`, after `EDIT` and `CARRY FORWARD` | the modes shown in the window title |
| `edit` | `d` after a pair is removed (`info.row`), `c` and `t` (all of a slice's points replaced, `info.row` empty) | `gui_data`, for a plugin that keeps data per row; `info.slice` |
| `save` | `s`, and a close that saves, after `atlas2histology_tform.mat` | `[]`; `info.histology`, `info.atlas`, the points as saved |

At `open` a plugin may also set three fields the GUI reads: `outlier_rule`
(the `*` mark's thresholds, LS5; `factor` 3, `min_px` 30, `slice_px` 20 by
default), `controls_extra` (lines on top of the right column of the controls
window; none by default) and `unsaved_hint` (the last line of the note
printed when a save leaves orange slices out; by default "Grab any point on
one to accept it as it stands."). Without a plugin every call hands its value
back as given, so the GUI is LS3 to LS5 alone. The change to the file: the
contract in the help text, the four fields set at start-up, the eleven
calls, the helper, `affine_outliers` taking the rule as an argument, and
`label_points` drawing a plugin's mark where it drew `?`. The file names no
code of this project. Effect: none.

Our plugin is `registration/annotation_gui/auto_annotation_plugin.m` (the AA
layer below), with its settings in `registration/annotation_gui/annotation_settings.m`.
`run_register_to_atlas`'s `annotate` mode passes it when
`auto_annotate('check')` finds the engine's interpreter, the test every call
to the engine starts with; without the engine it says so and the GUI opens
without the plugin.

**LS8. The `r` key: landmark proposals** (`7f42fdf`, 3 Sep 2026; removed in
step 6 of the refactor, decision L3). `r` proposed points for a slice by
matching the neighbouring slice's landmarks through an image matcher
(Python, then at the code root). Removed from the GUI: the key, the function
that called the matcher (`auto_refine_points`), the channel it matched on,
the key's block in the controls window (now one for `t`), and the hint
"Press r for a refined proposal" after a carry-forward. `t`, `p`, the
numbered labels and the `?` flags stay (LS5, AA). The matcher (its MATLAB
wrapper and worker, its Python folder, its setup script and the frozen
requirements of its environment) is in `archive/`, and the registration
driver's `annotate` mode no longer starts its worker.

**AA. The automatic-annotation layer** (`f437c8a`, `9fc2700`, `7e90bbd`,
`26ee86d`, `710d6b2`; 29 Sep 2026). About 550 lines: `a` sets a plane anchor
and `j` jumps between the suggested ones; a proposal from
`run_register_to_atlas`'s `autoannotate` mode loads provisional, with the
least confident points marked `?`; `k` / `K` accept, `u` / `U` re-propose
through `auto_annotate`; saving also writes `plane_anchors.mat`,
`auto_atlas_planes.mat` and `annotation_provenance.mat`. A brain with no
anchors and no proposal behaves as before. Effect: none on the registration
code; the annotation is still an ordinary `atlas2histology_tform.mat`, written
from what was accepted.

Moved out of this file in step 6 of the refactor (1 Oct 2026), into
`registration/annotation_gui/auto_annotation_plugin.m`, behind the hook LS7.
What moved, statements unchanged: the loading of the anchors and the proposal
(at `open`), the keys `a`, `j`, `k`, `K`, `u`, `U` (at `key`), saving the three
files (`save_auto_files`, at `save`) and its helpers `remember_proposal`,
`predict_planes` and `write_atlas_planes`, the review stage in the window
title (at `window`), and the start-up draw of a brain with anchors or a
proposal (at `open`). Rewritten as the hook's answers, same results: the
anchors in the plane prediction (at `planes`), `PROPOSED`/`ACCEPTED`,
`PLANE ANCHOR`, `SUGGESTED ANCHOR` and "plane changed" in the slice title (at
`title`), the `?` flags (at `labels`, and kept in step with `d`, `c` and `t`
at `edit`), the keys in the controls window (`controls_extra`) and the note's
hint about `k` (`unsaved_hint`). The suggested anchor count (4) and the
outlier rule are `annotation_settings`. With the plugin the GUI behaves as
before step 6, without `r`.

## Changes that do nothing

**P0a. `ls_analyze_slice_volume.m`, the demo script** (before the first
commit). Its data folder and mouse name were pointed at this project
(`D:\sep_histology`, `CG027`), its `SliceOrderEditor` call commented out, and
a few blank lines moved. The pipeline never runs it: the drivers call
LightSuite's functions themselves. It is on the do-not-run list
(`third_party/README.md`). Returned to upstream's version in the refactor,
so it no longer differs from upstream.

**Step 4 of the refactor** (this commit): the one-line notice at the top of
each file above (a comment, so the code is identical), and in
`matchControlPointsInSlices.m` the name of the registration driver in its
messages and comments (`P4` became `run_register_to_atlas`). Since step 6
those messages are the plugin's, and the GUI file no longer names the driver.

## Upstream contributions (Z2, after the refactor)

LS1 to LS7 are kept small and separate so each can become one pull request,
after rebasing on upstream (its GUI now has its own numbering and binds
backspace to "delete last point"; its SliceOrderEditor was rewritten). Also
worth asking for, from our workarounds: saving the per-slice transforms, an
output-folder option, relative elastix paths, a seed in the align step, the
demo's `settings` bug, and extentfactor, the CPD iteration count and the
B-spline bins and samples as options with upstream's defaults. The mechanics:
a fork holding the patch series, vendored here with `git subtree`.

## Replacing this copy with upstream (Z1, after the refactor)

Upstream was 168 commits ahead on 30 Sep 2026. Changes there that alter
results: the B-spline settings (every registration), the align extent (6 to
10, so plane indices shift by 30 and re-aligning an annotated brain
invalidates its points), a lower CPD iteration count, and an extraction that
no longer crops and that edits the curated decisions file (which would
invalidate the artifact masks and the residual-correction outputs). Also:
`run_order_slices` must then use upstream's SliceOrderEditor (ours drops the
crop columns on save); `run_add_sep_channel` gains the crop step or refuses;
`cpwt` stays 0.2 (upstream's demo uses 0.4); the path starts at
`<LightSuite>/src`; our matlab_elastix is checked against the fork upstream
now asks for. So: never re-run extraction or align on an existing brain;
measure the noise floor of the current code first (register only); register
the reference mice with upstream from their existing align outputs and control
points and compare voxelwise, per region, and on the headline numbers;
re-validate the automatic annotation, whose models were trained on the current
geometry; record which version processed which cohort.
