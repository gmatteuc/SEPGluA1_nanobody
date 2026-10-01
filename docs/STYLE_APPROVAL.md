# Style guide and reference examples: for approval

A note for Giulio, to delete once the guide is approved. It says what the
start of step 7 changed, what to decide, and what was left for later.

## What was written

- `docs/STYLE.md` (447 lines): shared rules (111 lines), the MATLAB half
  (110: the imaging repository's guide adapted, Y6, with a table of where
  the P8/P10 conventions differ and which wins) and the Python half (202: the
  draft finalised, Y1, Y2, Y5, Y7). A rule that holds in both languages is
  stated once, in Shared rules; each half keeps the language's form of it,
  with one short example per rule, quoted from the reference examples or
  the reference projects (Y7). A first version ran to 1,318 lines; the
  condensed one keeps every rule and drops the repetition, the extra
  examples, the line citations into the reference projects and the
  explanations a reader does not need to follow a rule.
  `docs/STYLE_python_draft.md` is removed.
- `docs/STYLE_PASS.md` (292 lines), temporary, deleted when step 7 is done:
  the procedure of the style pass (kinds A, B and C of change, the checks
  for each file), the code of the reference examples not to copy, what
  applies to old code only, the reference projects' habits the guide does
  not follow, the pipeline README template and the environment recipes.
- `ruff.toml` at the root: line length 90, Python 3.12, rules E, W, F,
  double quotes, `third_party` and `archive` excluded.
- The four reference examples, restyled. Only comments, headers,
  docstrings, quotes, spacing, line breaks and section titles changed.
- `.git-blame-ignore-revs`, listing the restyle commit, which holds the
  four files and nothing else.

Checks: `sep_check_code_identity` gives `same code` for both MATLAB files,
`check_code_identity.py` for both Python files; a planted change to one
number (`min_num_contrib = 4`, `SMOOTH = 1.5`) is reported as `CODE CHANGED`
by each tool. `checkcode` gives the same 1 and 8 messages as before (one
unnecessary bracket; eight unused variables), at shifted lines. `ruff format
--check` and `ruff check` pass on both Python files (before: 43 ruff
messages in `compare.py`, 33 long lines and 10 `;`).

## What changed in each reference file, and why

**`group_comparison/run_collect_by_group.m`.** The header moved above
`clear all`, so `help run_collect_by_group` shows it, and took the imaging
layout: name, title banner, the three steps with this one marked, what it
does, the setup. `%% User-defined parameters` became `%% Settings`. Each
setting has a lowercase comment of one or two lines above it instead of a
sentence-case one plus an end-of-line note. The comment of
`correction_type` now says it is not passed to the code.

```matlab
% before
clear all
close all
clc

% /// Plasticity comparison, step 1 of 3: collect the registered volumes by group ///
% For each group, stacks the registered nano volumes of its mice into one 4D

% after
%% run_collect_by_group
% ===== Collect the registered volumes of each group =====
%
% Plasticity comparison, step 1 of 3:
%   1. run_collect_by_group    stack each group's registered volumes  <- this script
```

```matlab
% before
% Cohort selection (mice come from the shared registry get_cohort.m).
mousetypes_list = {'young'};                    % 'rws' | 'naive' | 'behavior' | 'young'

% after
% groups to collect ('rws', 'naive', 'behavior', 'young'); the mice of each come
% from the cohort registry, get_cohort
mousetypes_list = {'young'};
```

**`group_comparison/pipeline/collect_by_group.m`.** Comments lowercase,
one or two lines, one per step, a blank line before each. The `%% Atlas`
section, which held only a comment, is gone: its point (every volume of a
cohort is on one grid, which is what lets `cat(4, ...)` work) now sits above
the `cat`. The 17-line warning about the stale `auto_4d.mat` and the
averages became a paragraph of the help block, with its dates and numbers;
the story of how it was found went (the plan's known defects keep it). The
comment above the commented-out saves says that why they are off is not
recorded. Long `fullfile` lines and the `clear` list continue with `...`.
Two printed messages stay at 99 and 101 characters: a string cannot be
split without changing the code.

**`mapping/run_compare.py`.** The header gives the run order before what
the script does, as the MATLAB headers do, and now lists the four files it
writes (moved from the module's docstring); `main()`'s docstring says what
it does instead of repeating the module's.

**`mapping/sepmap/young_vs_adult/compare.py`.** The module docstring lost
its history ("v2, step 4", "Nothing is transformed here any more") and the
interpreter path, and now gives the method and the reason in 22 lines, and
the run script that calls it. Every function has a docstring; the two
5- and 7-line comments went into docstrings. The end-of-line comments of
the constants moved above them. Statements joined by `;` were split, the
backslash continuation became brackets, quotes became double, and
`ruff format` set the layout. Two stale statements were corrected: the
pooled young group is the P16, P20 and P22 brains, not "P20 + P16" (the
constant's comment) nor "every registered young brain" (the docstring;
MG914, P28, is registered but not in `YOUNG`).

```python
# before
CCF_AP0 = 180                        # the adult registered crop, 10 um planes
...
SMOOTH = 1.0                         # voxels at 20 um, applied to the log2 map only
YOUNG = 'young'                      # pooled P20 + P16

# after
# start of the adult registered crop along AP, in 10 um planes; added to the plane
# number in each slice title
CCF_AP0 = 180
...
# Gaussian sigma in 20 um voxels, applied to the log2 map only
SMOOTH = 1.0
```

```python
# before
def fold_n(n):
    ml = n.shape[2]; h = ml // 2

# after
def fold_n(n):
    """Fold a count map as fold does, keeping the larger count of the two sides."""
    ml = n.shape[2]
    h = ml // 2
```

```python
# before
        vmax = float(np.nanpercentile(adult_v[iso], 99)) if not signed else \
            float(np.nanpercentile(np.abs(adult_v[both]), 98))

# after
        vmax = (
            float(np.nanpercentile(adult_v[iso], 99))
            if not signed
            else float(np.nanpercentile(np.abs(adult_v[both]), 98))
        )
```

## Open choices

1. **A formatter for Python.** The draft said no formatter, since the
   references are formatted by hand (black would change 361 of their 2,918
   lines). The guide now uses `ruff format`: it does Y2's double quotes and
   the 90-character lines mechanically, it never changes the code (the
   identity check confirms it), and a newcomer needs one command instead of
   hand rules. Its one cost here is visible above: a long conditional
   expression is spread over several lines, which the guide's own rule then
   turns into `if`/`else` (kind B). **Recommendation: keep `ruff format`.**
   Without it, delete the `[format]` table of `ruff.toml`, the first bullet
   of the Python half's Layout and naming, and `ruff format` in STYLE_PASS.md
   (kind A, Python check 4); the reference files stay valid.
2. **MATLAB settings in lower case.** The imaging guide writes settings in
   upper case (`SESSIONS`). The guide keeps `snake_case` as a project rule:
   each setting goes to the pipeline function under its own name, and
   `sep_run_driver_copy` commands and `production_settings.md` name them.
   **Recommendation: keep lower case.** Upper case would rename the
   settings of 11 drivers and every command that names one.
3. **`clear all`.** The guide follows the imaging guide,
   `clear; clc; close all;`. `clear` empties the workspace as `clear all`
   does; what it keeps is the state inside functions: `persistent` and
   `global` variables. Our own MATLAB code has none. The third-party code
   has some, and none of it changes a result: LightSuite's `panel.m` keeps
   its property defaults, a debug flag and its figure registry, and locks
   itself with `mlock`, so `clear all` never cleared it anyway; `inifile.m`'s
   global `NL_CHAR` is set again on every call, from the platform, and only
   LightSuite's Jetraw scripts call it; yamlmatlab's `ReadYamlRaw` keeps its
   two call options and a verbosity level, and its one caller here
   (`elastixYAML2struct`, through `ReadYaml`) always gives the defaults. `clear all` is in
   23 live files: the 11 drivers, the QC scripts (4 in `atlas/qc`, 2 in
   `registration/qc`), `make_ordering_volume`, `make_atlas_reference_sheet`,
   `tests/test_backward_compat.m`, and P8 to P10. **Recommendation: change
   the drivers, the QC scripts, the two tools and the test in the style
   pass (kind B, one rerun of each pipeline), and the "clear all" trap in
   `tools/README.md` with them; leave P8 to P10, which retire.**
4. **Path setup.** The imaging drivers call `setup_paths` themselves; ours
   rely on `sep_setup_paths` once per session, because a driver copy run by
   `sep_run_driver_copy` sits in `driver_copies\`, where `mfilename` cannot
   find the code root. **Recommendation: keep it as is** (the guide says
   so).
5. **Dead code in `collect_by_group`.** It loads the auto and mask volumes
   and computes the averages and the relative difference only for saves that
   are commented out (eight of them), and sets four folder names it never
   uses. **Recommendation: remove them in step 8, one commit, outputs
   unchanged; the help block keeps the warning about the old files.**
6. **Who does the structural changes.** The plan leaves `matplotlib.use` and
   the module constants "for the style pass" (Progress, step 5, Python
   route), while the guide's first version put them outside it, so no step
   owned them. STYLE_PASS.md now has a kind C for them: structural changes
   (backend into the run scripts, constants into `settings.toml`, `pathlib`,
   merging the five `save_figure` and the three `fold_n`, a plotting module
   and a MATLAB palette function, splitting the functions still over 60
   lines), each in its own commit after the kind A and B commits of the
   files it touches, each checked by a rerun. The alternative is a named step
   between 7 and 8. **Recommendation: kind C inside step 7, as STYLE_PASS.md
   now says; the plan's step 7 then names it.**

Decided in the guide while reconciling the earlier conventions, to confirm:

- A comment is one or two lines, but a workaround for how MATLAB behaves
  (as P10's on `barh`) or a justified setting may take more.
- Purple, pink and green stay out of new colours; `PuOr_r`, magma and
  plasma stay, as existing figures use them.
- The palette is typed in each script until one function holds it (kind C).
- Message pieces in `error([...])` stay aligned under the first piece, as
  the code already does everywhere; other continuations indent by four.

## Not changed: renames and code the guide would want

Renames of local names (saved variables, struct fields, columns and
settings keep theirs):

| file | now | suggested |
|---|---|---|
| `collect_by_group.m` | `current_mouse_type`, `mousetype_idx` | `group`, `group_idx` |
| | `num_current` | `n_mice` |
| | `file1_name`, `file2_name`, `file4_name` | `nano_file`, `auto_file`, `mask_file` |
| | `nanoVols_type` (and auto, mask) | `nano_vols` |
| | `allenDir` (not a LightSuite name) | `atlas_dir` |
| | `has_reg`, `min_num_contrib` | `is_registered`, `min_n_mice` |
| | `diff_4d_new`, `avg_diff_new` | `rel_diff_4d`, `avg_rel_diff` |
| `compare.py` | `m` (a reading), `z` (the saved maps) | `reading`, `maps` |
| | `both` (the voxels compared) | `compared` |
| | `r`: the comparison map at 332, a row at 403 | `contrast`, `row` |
| | `w`: the weights at 339, the CSV writer at 422 | `weights`, `writer` |
| | `a`, `p` (the adult and young maps, 327-328; `p` reads as a p-value) | `adult_v`, `young_v` |
| | `eps` (the floor; `eps` is also the EPS path in `save_figure`) | `floor` |
| | `ann`, `ann_h`, `iso`, `cov`, `ok`, `labs` | `annotation`, `annotation_left`, `isocortex`, `coverage`, `covered`, `structure_of_voxel` |
| | `nm, ac, dv`, `g`, `tot`, `cnt`, `lim2` | `name, acronym, division`, `structure`, `total`, `count`, `diff_lim` |

Renames of names used outside their module:

| name | used by | suggested |
|---|---|---|
| `compare.fold_n` | imported by `video_compare.py:43`; defined again in `closeup.py` and `video.py` | `fold_count`, with the import changed (better: one copy, kind C) |
| `compare.CCF_AP0` | a module constant | `CROP_START_PLANE` |

Also: the `log2_zref` column of `region_table.csv` holds a difference, not a
log2 ratio (a column name, so it stays unless the table changes anyway).

Code changes the "do not copy" list (STYLE_PASS.md) names: `clear all`
(choice 3); in `compare.py`, imports in isort order with
`LinearSegmentedColormap` at the top, the four multi-line conditional
expressions into `if`/`else` and type hints (kind B), and the backend,
`os.path`, `SMOOTH` and the 100-voxel minimum, the settings table, the
drawing and `save_figure`, and the two long functions (kind C, choice 6).

## Possible bugs noticed (for the plan's bug list)

- `compare.draw_figures` titles each slice `plane {zc * 2 + CCF_AP0}`. `zc`
  indexes the full 660-plane CCF grid at 20 um (annotation planes 1-658 in
  the reference run's `volumes_ccf20.npz`), so `zc * 2` is already the CCF
  plane at 10 um, and the titles read 180 planes too high: 526 to 1074 for
  CCF planes 346 to 894. `closeup.py` converts with `plane // 2`, no offset.
- `run_collect_by_group`'s `correction_type` is set but not passed, so it
  does nothing.
- P8 colours its bars with the first 200 levels of `flipud(gray(256))`, and
  its comment says this avoids pure white at the low end; it does the
  opposite: level 1 is white (t = 0) and the dark end is cut. Bars of
  unreliable structures are white on white. P8 retires, but A4 should not
  copy it.
- The draft said `matplotlib.use` between imports makes every later import
  an E402. ruff 0.16.9 does not flag them (checked); the guide now gives the
  real reason to set the backend in run scripts.

## Left for the plan

- `docs/REFACTOR_PLAN.md` (line 22) links to `STYLE_python_draft.md`, which
  no longer exists; to update where the plan is maintained.
- The plan's rule puts "the style-pass commits" into
  `.git-blame-ignore-revs`. STYLE_PASS.md narrows it to kind A commits: kind
  B and C change code, and blame should show them.
- Step 7 of the plan to name kind C, if choice 6 is taken, and to point to
  `docs/STYLE_PASS.md` for its procedure.
