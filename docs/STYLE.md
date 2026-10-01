# Code style

Reference for all code in this repository, MATLAB and Python. New code
follows it; old code is brought in line in the style pass (step 7 of
[REFACTOR_PLAN.md](REFACTOR_PLAN.md)), whose rules are at the end.

The rules shared by both languages come first, then the MATLAB half, then the
Python half. The MATLAB half is the style guide of the imaging repository
(`D:\dendrites\code\docs\STYLE.md`), adapted to this project (decision Y6).
The Python half follows two reference projects by the same author, the
Genedata and MaxWell exercises, quoted wherever a rule comes from them, so
the style pass works from this file alone (Y7).

Before this guide, the code followed conventions taken from P8 and P10
(script skeleton, comments, naming, figures). Where they differ from the
imaging guide, the imaging guide wins (Y6), unless a rule specific to this
project applies, such as the colormaps. The table at the end of the MATLAB
half lists each difference and which rule wins.

## Reference examples

Files already in final style, to copy from:

| file | shows |
|---|---|
| `group_comparison/run_collect_by_group.m` | a MATLAB driver: header, settings, one call |
| `group_comparison/pipeline/collect_by_group.m` | the function a driver calls: help block, settings unpacked, `%%` steps, comments |
| `mapping/run_compare.py` | a Python run script: header with the run order, `main()`, options |
| `mapping/sepmap/young_vs_adult/compare.py` | a package module: docstring with the method, constants, docstrings, comments |

Their comments, headers, docstrings and layout follow this guide. Some of
their code does not yet, because a comment and layout pass does not change
code; do not copy these:

- `run_collect_by_group.m` starts with `clear all` on three lines, and its
  setting `correction_type` is not passed to the code;
- `collect_by_group.m` sets folder variables and computes averages it never
  uses, sets `min_num_contrib = 3` in the middle of the code, and adds the
  atlas folder to the path though nothing in it reads the atlas; two of its
  printed messages run to 99 and 101 characters, since a string cannot be
  split without changing the code;
- `compare.py` sets the figure backend in the module, imports
  `LinearSegmentedColormap` inside a function, keeps its imports out of isort
  order, builds paths with `os.path`, keeps parameters (`SMOOTH`, the
  100-voxel minimum) in the code, writes conditional expressions over several
  lines, and its public functions have no type hints. It also departs in
  structure: it draws figures although it computes (figures belong in a
  plotting module), has its own `save_figure` (one of five copies in the
  package), reads its settings key by key into constants
  (`MIN_N_YOUNG = SETTINGS["young_vs_adult"]["min_n_young"]`) instead of
  naming the table, and its `draw_figures` and `main` are 142 and 174 lines
  long.

## Shared rules

### Voice

The code reads as if written by a careful scientist for a colleague.

- Plain and short. No capitals for emphasis, no "NOTE:", "IMPORTANT" or
  "!!!", no bold.
- No history: no "previously", "was changed on", "fixed a bug where", no
  story of how a bug was found. History belongs in the git log, and in
  `docs/` when a reader needs it.
- No inflated or stock wording: no "comprehensive", "robust", "elegantly",
  "leverages", "ensures", "Here we".
- No over-engineering: no options, layers or defensive branches nobody needs.
- Keep expressions simple: one step per line, intermediate variables with
  readable names rather than nested calls. A plain loop is better than a
  clever vectorised line that nobody can read.
- Numbers where there are numbers (a threshold and why, a count, a size).
- Text the program prints, writes or draws (messages, figure titles, file
  names, column names) is code, not a comment: changing it is a code change.

### Comments

- Short and lowercase: `% load traces`, `# per-plate quality control`. Verb
  first where it reads naturally, otherwise a noun phrase. A comment that
  starts with a capital letter is wrong, unless it starts with a proper
  name, an acronym or a code identifier (`% Gaussian sigma ...`,
  `% DeMBA grid ...`).
- On their own line, just before the code they describe. One per step: inside
  a loop each sub-step (load, skip, initialise, accumulate) gets its own.
- No comments at the end of code lines. The exceptions are pragmas read by
  tools (`%#ok<NANMEAN>`, `# noqa: E402`) and `settings.toml`, where the end
  of line is the units column.
- Give the reason only when it isn't obvious, in one line:
  `# the atlas in flat grey under the data, so no data reads as grey`.
- For important decisions (a threshold, a window, what gets excluded) say
  briefly why, and what the alternative would get wrong, in one or two lines:

```matlab
% weight by frame count: files differ in length after bad-frame removal,
% and a plain mean of means would favour the short ones
```

- Don't restate the code (`% loop over files` above a loop over files), but
  do put a short comment before any line that is dense or not obvious
  (`% rescale to 0-1 and clip values outside the limits`).
- A comment is one or two lines, wrapped at 90 characters like the code. Two
  kinds may take the lines they need: one above a setting or a constant that
  justifies its value (`mea/trajectory.py:10-14`), and one that explains a
  workaround for how MATLAB or a library behaves. A longer explanation of a
  method goes into the help block (MATLAB) or the docstring (Python) of the
  function or module that implements it.
- Papers are cited by author and year, in one line: "(Hill 2011, as in
  SpikeInterface)" (`mea/units.py:64`). What a paper shows and what it means
  here goes into `docs/SCIENTIFIC_CONTEXT.md` (decision L7), not into a
  header, help block, docstring or comment.
- No commented-out code in new code; git keeps old versions. Existing
  commented-out code stays until the owner decides, with a comment saying why
  it is off, or that the reason is not recorded.

### Messages

- An error says what was found, what was expected, the path, and the fix
  when there is one.
- Progress is printed, one lowercase line per step, sub-items indented by
  two spaces. A case that is skipped while the run goes on is a warning.

### Data safety

- Raw data is read only: nothing writes under `S:\` or to the raw `.czi`
  files in `data\<group>\<mouse>\`.
- Outputs go under the data root, and only there: `get_paths().data`
  (MATLAB) or `config.DATA` (Python). No code writes a `D:\` path of its own.
- `SEP_DATA_ROOT` moves the whole data tree, inputs and outputs, for checks.
  `get_paths.m` and `sepmap/config.py` refuse the production data to a copy of
  the code, and any data root inside the snapshot on G:. A check never writes
  into the production data.
- Code does not delete files. A step overwrites only its own outputs.
- Excluded data, and a value that cannot be determined, is NaN, never zero,
  with the reason in a comment; NaN is left out of averages.
- An expensive result is cached, and a switch recomputes it
  (`force_recompute_*` in MATLAB, an option named for what it redoes in
  Python). Loaded or cached data is checked against what it must hold before
  it is used (the mice it holds against the mice asked for), and the run
  stops if it differs. A cached table that was sorted is put back in the
  order of its list on loading (`ismember`).

### Shared definitions

Each of these is defined in one place; never copy a value into a script.

- Paths: `get_paths.m` (MATLAB) and `mapping/sepmap/config.py` (Python),
  both worked out from where the code sits (the parent of the code folder,
  plus `data`).
- The MATLAB path: `sep_setup_paths`, once per session.
- The cohort: `get_cohort` (MATLAB; `get_cohort_spec` for the mice a stack
  holds) and `COHORTS` in `mapping/sepmap/volumes/cohort.py` (Python). One
  table for both languages comes in step 8 of the plan (Y4). New code takes
  its mice from these, never from a list of its own.
- The atlases: `get_atlas(key)`, `cohort_atlas_key` and `get_atlas_crop`.
  `get_atlas` adds the one atlas folder a run needs to the path; no other code
  should add one. LightSuite finds the atlas with `which`, so a second atlas
  folder on the path is picked up without a warning.
- Parameters: a MATLAB driver's `%% Settings`; for Python,
  `mapping/settings.toml`, read by `sepmap/config.py`. A number someone could
  reasonably choose differently is a parameter, never a number in the middle
  of the code. A number fixed by the method stays in the code, with a comment
  saying what it is, and so does the layout of a figure (sizes, font sizes).
- Colours and colormaps: the values are in Figures, below. Neither language
  has a palette module yet: the MATLAB scripts name the colours in their
  settings block, the Python modules type them where they draw. A Python
  plotting module and a MATLAB function in `common/` are structural changes
  of the style pass (kind C); until then, new code copies the values below
  exactly.

### Figures

Colormaps:

| what | colormap |
|---|---|
| intensity (nano, autofluorescence, a reading) | hot |
| atlas, anatomy, raw images | gray |
| a difference (young minus adult, control against experimental, t maps) | blue-red, symmetric limits, zero in the middle: `get_color2color_colormap([0 0 1], [1 0 0])` in MATLAB, `RdBu_r` in Python |
| a signed position within a brain (zref) | purple-orange, `PuOr_r` (purple low, orange high), so that blue-red keeps meaning a difference |
| counts of brains (n maps) | magma |
| coverage lines | plasma |

- Never parula. jet or turbo only when asked for, into a subfolder of their
  own, never as the default (`run_closeup --cmap`); the difference panel
  stays blue-red, whose zero is a real null.
- No data is flat grey (`#bfbfbf`), a colour no data colormap above produces
  (gray is for the atlas and raw images only). In the young-against-adult
  maps hot stops at 0.82 of its range, so saturation reads as bright yellow
  and never as white.
- Bars of a value per structure across mice: the height is the value, the
  colour its reliability across mice (a t value, clamped): grey, darker for
  more reliable, never pure white, so every bar stays visible on the white
  background. Bars that compare categories (models, groups of genes) take
  colours of the palette.
- Palette: nano orange `[0.95 0.55 0.10]`, autofluorescence yellow
  `[0.95 0.85 0.20]`, their per-mouse dots darker (`[0.65 0.30 0.00]`,
  `[0.70 0.60 0.00]`), lines joining paired mice grey `[0.6 0.6 0.6]`. Groups
  in the Python figures: young `#c0392b`, naive `#555555`, rws `#9a9a9a`. New
  colours extend the orange, yellow, blue, red and grey family: no green,
  pink or purple. The only exceptions are the colormaps of the table
  (`PuOr_r`, magma, plasma), which existing figures use.
- Scatter plots of many structures: 35-point dots, no edge, alpha 0.85.
- Every count shown in a title is computed, never typed.
- Coronal planes are drawn dorsal up: a plane is (DV, ML), never transposed.
- An underscore in a MATLAB title is a TeX subscript: give the text object
  `'Interpreter', 'none'`, or replace the underscores
  (`strrep(tag, '_', ' ')`).
- A file name carries every setting that changes what the file holds
  (`_smooth5` or `_nosmooth`, the channel, `_P20`), so runs with different
  settings sit side by side. Existing output names do not change.
- White figure background, except panels that show images.

## MATLAB

### Files and folders

- A pipeline is one folder: its drivers (`run_*.m`) at the top, the
  functions they call in `pipeline/`, checks in `qc/`. Shared functions are
  in `common/`, atlas functions in `atlas/`.
- A driver is a script: a header, the settings, one call into `pipeline/`.
  Everything else is a function.

### Drivers (run_*.m)

- Header first (layout under Headers), then `clear; clc; close all;`, the one
  line with several statements.
- No path setup in the driver: `sep_setup_paths` runs once per MATLAB session,
  before any driver, and the detached runner runs
  `restoredefaultpath; cd(CodeDir); sep_setup_paths` before each stage. This
  differs from the imaging repository, whose drivers find the root with
  `mfilename` and call `setup_paths` themselves: here a check runs a copy of
  the driver from `<data_root>\..\driver_copies\` (`sep_run_driver_copy`),
  where `mfilename` would point at the copy's folder.
- Then a `%% Settings` section with everything the user is expected to
  change. Each setting has a comment above it, with the allowed values where
  useful:

```matlab
% groups to collect ('rws', 'naive', 'behavior', 'young'); the mice of each come
% from the cohort registry, get_cohort
mousetypes_list = {'young'};
```

- The first setting is `paths = get_paths();`, with no folder written out.
- Settings have lowercase `snake_case` names, like every other variable.
  This differs from the imaging guide, which writes them in upper case: here
  each setting goes to the pipeline function under its own name, as a field
  of `run_settings` and then as a local variable, and `sep_run_driver_copy`
  and `docs/production_settings.md` address settings by those names.
- Then `%% Run`: the settings go into the struct `run_settings` under the
  same names, and one function call does the work:

```matlab
%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.mousetypes_list = mousetypes_list;
run_settings.age_filter = age_filter;
run_settings.skip_missing = skip_missing;
collect_by_group(run_settings);
```

- A driver computes nothing beyond its settings. A value derived from the
  settings is computed in the function.

### Pipeline functions

- The function a driver calls takes `run_settings`, and its help block names
  the driver. The settings are unpacked first, one per line:

```matlab
% settings of run_collect_by_group, under the names the code below uses
paths = run_settings.paths;
mousetypes_list = run_settings.mousetypes_list;
```

- A body longer than about 60 lines is split by `%%` sections into steps
  (`%% Collect each group`, `%% Save`). A long function is a short main
  part of `%%` steps that calls local functions, each with explicit inputs and
  outputs; no nested functions sharing a workspace, no globals, no `evalin`
  or `assignin`.
- `load` always into a struct: `S_nano = load(file, 'nano_4d');`, never
  into the workspace, where a loaded name can shadow a function.
- Options other than a driver's settings are passed as a struct, and a
  missing field gets its default at the top of the function.
- A number in a pipeline function that someone could choose differently is a
  setting of its driver (Shared definitions), not a variable set halfway
  down.
- Local functions come after the main part, after one separator line in the
  style of the header banner (not a `%%` section, since they are not steps to
  run):

```matlab
% ===== Local functions =====
```

  A long file with several groups of local functions uses one per group:
  `% ===== Local functions: loading and alignment =====`.

### Headers

Every file has one, local functions included. Sentence case (capitalised),
unlike the lowercase inline comments. Keep them short.

- Drivers: the big picture first: which pipeline, which step this is, and
  what the step does in general. Only then a few words on the data it is
  currently set up for. The header comes before any code, so
  `help run_collect_by_group` shows it:

```matlab
%% run_collect_by_group
% ===== Collect the registered volumes of each group =====
%
% Plasticity comparison, step 1 of 3:
%   1. run_collect_by_group    stack each group's registered volumes  <- this script
%   2. run_normalise_groups    bring the mice of a group onto one intensity scale
%   3. run_group_differences   left-right differences, control against experimental group
%
% Stacks the registered nano volumes of a group's mice into one 4D array
% (AP x DV x ML x mouse), saved in the group's folder for run_normalise_groups.
% ...
%
% Setup: the young cohort, P20 brains only, for the comparison with the
% adults. For the plasticity comparison, the groups are rws, naive and
% behavior, with no age filter. Run sep_setup_paths first, once per MATLAB
% session; the code is in pipeline\collect_by_group.m.
```

- Functions: the usual MATLAB help block: the name in capitals and a summary
  line, the call syntax, inputs, outputs and options with defaults. No
  essays.

```matlab
function collect_by_group(run_settings)
%COLLECT_BY_GROUP  Stack the registered volumes of each group's mice.
%   COLLECT_BY_GROUP(run_settings) does the work of run_collect_by_group,
%   which sets the fields of run_settings (paths, mousetypes_list,
%   age_filter, skip_missing) and says what each one does.
```

- Local functions: one or two sentence-case lines under the signature, then
  a blank line:

```matlab
function [base_dir, global_diagnostics_dir, data_4d] = load_cohort_stack(paths, ...
    current_mouse_type, channel, S, subset_indices)
% The cohort's folders, and its stack of the channel with the selected mice only.
```
  (`group_comparison/pipeline/normalise_groups.m:151-153`)

### Spacing and sections

- A blank line before each comment block, so each block reads on its own:

```matlab
% get the mouse
mouse_name = current_mice{i};

% set its folders (base_dir is the group's)
base_dir = fullfile(paths.data, current_mouse_type);
```

- A blank line at the start of a loop or `try` body that holds several steps,
  and between a function's help block and its code.
- `%%` sections with a short capitalised title, not numbered, split the code
  into steps: `%% Settings`, `%% Run`, `%% Collect each group`. Every script
  has them, and so does the body of any function longer than about 60 lines.
  Short functions have none.
- One statement per line: `if ~exist(out_dir, 'dir')`, `mkdir(out_dir);` and
  `end` on three lines.
- Spaces around `=` and after commas. Lines under about 90 characters; a long
  call continues with `...` and an indent of four. The pieces of a message
  joined in `[...]` line up under the first piece (the `error` below).
- A string literal that does not fit stays whole on its line: splitting it
  into `['...' '...']` changes the code (kind B). A list of names after a
  command (`clear a b c`) can continue with `...`; the code stays the same.
- Text is in single quotes, as character vectors: `'lightsuite'`. Y2 (double
  quotes) is for Python only: in MATLAB `"..."` is a string, another type, so
  changing the quotes changes the code, and the style pass leaves them.

### Naming

- `snake_case` for variables and functions: `nano_dir`, `smooth_suffix`,
  `out_dir`. `camelCase` only for names that come from LightSuite:
  `nanoVol`, `sliceinfo`.
- `n_*` for counts (`n_mice`), `*_v` for a vector with one value per region
  (`nano_mean_v`), `S_<name>` for a loaded struct (`S_nano`).
- Switches start with what they switch: `produce_*`, `show_*`, `save_*`,
  `do_*`, `use_*`, `compute_*`, `force_recompute_*`.
- The style pass renames nothing; a rename this guide would want goes on a
  list for the owner.

### Errors and messages

- An `error` message starts with the name of the function or driver, then
  the shared rule (what was found, the path, the fix):

```matlab
error(['sep_setup_paths: %s is back at the code root (%s) under its name from ' ...
       'before the refactor, which moved it (docs/refactor_name_map.csv). An ' ...
       'editor tab saved after the move recreates it. Carry any edit over to the ' ...
       'moved file, then delete the old one.'], strjoin(back, ', '), code_dir);
```

- The check of loaded or cached data (Data safety), on three lines:
  `if ~isequal(nano_mice, auto_mice)`, `error(...);`, `end`.
- A string option goes through `switch`, with
  `otherwise error('Unknown agg_method: %s', agg_method)`.
- Progress with `fprintf`, a skipped case with `warning`:
  `fprintf('  %d mice paired across channels.\n', n_mice)`.

### Figures

- In code that runs unattended:
  `figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', 'Position', [0 0 1 1])`.
- Each figure is saved as `.fig` (`saveas`) and `.png` (`exportgraphics`,
  300 dpi); file names as in the shared Figures.
- Horizontal grouped bars: `barh(Y, 'grouped')`, then `EdgeColor = 'none'` on
  each series. A grid of panels has three columns and as many rows as needed:
  `n_rows_grid = ceil(n_panels / n_cols_grid)`.

### Where the earlier conventions differ

| topic | earlier (P8, P10) | imaging guide | here |
|---|---|---|---|
| start of a driver | `clear all`, `close all`, `clc` on three lines | `clear; clc; close all;` | imaging |
| driver header | `% /// Pipeline script #N: purpose ///`, then (1), (2) | `%% run_x`, `% ===== Title =====`, step list, setup | imaging |
| settings section | `%% User-defined parameters` | `%% Settings` | imaging |
| settings names | `snake_case` | `UPPER_CASE` | earlier (project rule, above) |
| comments | sentence case | lowercase | imaging |
| comment length | several lines for a why (P10, 9 lines on `barh`'s `YEndPoints`) | one or two lines | imaging, except a workaround for how MATLAB behaves, as in P10, or a justified setting (Shared rules) |
| a setting's units | end of line, `=` aligned in a column | comment above | imaging |
| numbers in the code | every configurable value at the top, no magic numbers | (no rule) | earlier, as a shared rule: a value someone could choose differently is a setting |
| section titles | numbered, `%% Figure 1: per-region paired bars` | short, unnumbered | imaging |
| local functions | `%% Local function: <purpose>` above each | one `% ===== Local functions =====` | imaging |
| short `if` | `if ~exist(d, 'dir'), mkdir(d); end` | three lines | imaging |
| path setup | `%% Add paths` in each driver | the driver calls `setup_paths` | neither: `sep_setup_paths` once per session (project rule, above) |
| saving figures | `.png` and `.fig` | `save_fig`: `.png`, and `.eps` on request | earlier (project rule: the existing outputs) |
| colours | hot, gray, blue-red; no green, pink or purple | condition colours | earlier (project rule), except the `PuOr_r`, magma and plasma colormaps of existing figures |
| where the palette is | named RGB triplets in each script's parameters block | one function, `condition_color` | imaging, once the function exists (kind C); until then earlier, with the values of Figures |
| bar colour | reliability, in grey | condition | earlier (project rule) |
| replaced files | overwritten by their step | moved to `superseded/` | earlier: there is no `set_aside_file` here |

## Python

### Sources

Two reference projects by the same author, read only on
`S:\ElboustaniLab\#SHARE\Documents\Giulio\`:

- the Genedata exercise, package `pka/`
  (`Genedata_exercise\PKA_exercise_Giulio_Matteucci.zip`);
- the MaxWell exercise, package `mea/`
  (`MaxWell_exercise\Retina_exercise_Giulio_Matteucci.zip`).

A module of a package is cited as `pka/fitting.py:75`, any other file with
its project folder, as `genedata/run_pipeline.py:35` or
`maxwell/tests/test_mea.py:79`. A rule with only citations is what both
references do; where only one does it, the rule names it. "Departure" marks
a rule that departs from them, with its source: the owner's additions (no
comment above `def` repeating the docstring, Y1; double quotes everywhere,
Y2; `raise` instead of `assert`; type hints on public functions; ruff at 90
characters), the MATLAB half, or the refactor plan. The rules apply to all
Python code, the verification tools in `tools/` included (Y5).

### Layout of a pipeline

- A Python pipeline is one folder. The entry points (`run_*.py`),
  `settings.toml`, `README.md` and `tests/` sit at the top, and the code they
  call is one package below them, as the MATLAB drivers sit above
  `pipeline/`:

```
pka/                   analysis package (used by the notebooks and the pipeline)
  config.py            loads settings.toml
  loading.py           raw JSON -> tidy well table + design checks
  preprocess.py        control statistics, QC, outlier flagging, normalization
  fitting.py           constrained Hill fits, activity classification
  selection.py         candidate selection
  plotting.py          shared plot style and figure functions
01_exploration.ipynb   data exploration and QC (pre-executed)
run_pipeline.py        runs the full analysis and writes all outputs
settings.toml          every threshold and fit parameter (single source of truth)
tests/                 a few sanity checks on the core functions
```
  (`genedata/README.md:9-25`, shortened; `maxwell/README.md:11-32` has the
  same shape)

- The entry points sit next to the package, never inside it:
  `python run_x.py` puts only the script's own folder on the import path, so
  a script inside the package folder cannot import the package by name.
- One module per stage, named after the stage: `loading.py`, `preprocess.py`,
  `fitting.py`, `selection.py` (pka); `loading.py`, `sorting.py`, `units.py`,
  `responses.py`, `trajectory.py` (mea). Departure: a route with several
  analyses has one sub-package per analysis, each laid out the same way
  (`mapping/sepmap/volumes/`, `young_vs_adult/`, `adult/`, `ish/`).
- `config.py` loads the settings and holds the paths. `plotting.py` holds the
  palette, the style, the save function and every figure function; the
  modules that compute draw nothing. A sub-package with figures of its own
  has its own `plotting.py`.
- Only the run scripts have an `if __name__ == "__main__":` block
  (`genedata/run_pipeline.py:102`, `maxwell/run_pipeline.py:162`).
- `__init__.py` is empty (`mea/__init__.py`) or holds a one-line docstring.
  `pka/__init__.py:3-16` re-exports every public name (28 unused imports for
  pyflakes); not followed.
- Notebooks are optional, numbered in reading order (`01_exploration.ipynb`),
  and only call the package: "All logic is in the `mea` package; cells only
  call it." (`maxwell/01_exploration.ipynb`, first cell). The first cell is
  markdown: a title, what the notebook shows, and that it uses the package;
  the first code cell imports, calls `set_style()` and loads the data. No
  result is computed only in a notebook.

### Run scripts (run_*.py)

- One `run_<step>.py` per step, named for what it does. Departure: each
  reference has one `run_pipeline.py`; our steps take hours, fetch from the
  network or need another environment (decision L1).
- The header is the module docstring: the summary line, the run order with
  this step marked (as in the MATLAB drivers), what the script does and the
  files it writes, the usage line, then its options. From the reference
  example:

```python
"""Young against adult in the adult CCF: maps and the per-structure table.

Python route, in run order (tools\\venv_atlas; run_closeup in tools\\venv_flat):
     1. run_per_mouse          per brain: tissue mask, backgrounds
     2. run_to_ccf             every brain on the adult CCF grid
     3. run_cohort             cohort mean, SD and n per voxel
     4. run_compare            young against adult: maps, table       <- this script
     ...

Folds the hemispheres of the cohort volumes and compares the young cohort with
the adults reading by reading (the method is in sepmap/young_vs_adult/compare.py).
Writes, in comparisons_v2/young_vs_adult/ under the data root:

    volumes_ccf20.npz     adult_*, young_*, log2_* maps, n maps, annot20 (AP, DV, ML half)
    ...

    python run_compare.py

V2_READINGS, a comma-separated subset of the readings (ratio, sepratio,
cref, subref, zref), leaves the others out of this step; the run prints the
readings in force.
"""
```
  The options paragraph as in Genedata: "--inhibitor picks the compound for
  the detailed curve figure; by default the most potent candidate is shown."
  (`genedata/run_pipeline.py:7-8`).

- `main()` reads as the list of stages. Each stage is a blank line, a
  comment naming it, the calls, and one printed line saying what was done,
  with the counts that show it went right:

```python
    # unit table and quality control
    table = units.unit_table(trains, templates, mapping, duration_s)
    table = units.apply_unit_qc(table)
    table.to_csv(out / "units.csv", index=False)
    n_ok = int(table["accepted"].sum())
    print(f"unit QC: {n_ok} accepted of {len(table)}, table written to units.csv")
```
  (`maxwell/run_pipeline.py:51-56`, without the stage number)

- `main()` holds only glue: the calls, a loop over items, output file names,
  the counts for the printed lines. Anything longer than a line or two of
  computing is a function in the package. Stage comments are not numbered
  (pka; mea numbers them, `# stage 4: ...`, not followed). Every run first
  prints the settings in force with `config.print_settings`.
- Run options (which mice, whether to recompute) are `argparse` options,
  parsed under `if __name__ == "__main__":` and passed to `main()` as
  arguments, so a test or a notebook can call `main()` directly
  (`genedata/run_pipeline.py:35, 102-108`; mea parses them inside `main()`,
  not followed). Analysis parameters are never options: they live in
  `settings.toml`.
- The cache switch (Data safety) is an option named for what it redoes:
  `--resort` reruns the spike sorting, otherwise the cache is read
  (`maxwell/run_pipeline.py:21`, `mea/sorting.py:75-84`). The cache file name
  is a setting (`maxwell/settings.toml:26`).
- Progress (Messages) with `print`; neither reference uses `logging`.
  Departure: `flush=True` on prints inside long loops, because long runs are
  detached and write to a log file.
- The figure backend is set only in run scripts, never in a package module,
  so notebooks still show their figures and importing a module changes no
  global state (pka sets it in its run script, `genedata/run_pipeline.py:14-16`;
  mea never sets it). Departure: it goes on the first lines under
  `if __name__ == "__main__":`, not between the imports as in pka, so the
  imports stay one block. `matplotlib.use` switches the backend even after
  pyplot has been imported (matplotlib 3.11,
  `matplotlib/__init__.py:1216-1280`):

```python
if __name__ == "__main__":
    # figures go to files, no display needed
    matplotlib.use("Agg")

    parser = argparse.ArgumentParser(description="per-mouse volumes")
    parser.add_argument("mice", nargs="*", help="mice to process (default: all)")
    args = parser.parse_args()
    main(args.mice)
```

### Settings and paths

- Every analysis parameter lives in `settings.toml`: one table per stage, one
  key per line, with its unit, and the reason when it isn't obvious, in a
  comment aligned in a column. The file starts by saying what it is:

```toml
# Single source of truth for all pipeline parameters.
# Loaded by pka/config.py; no threshold or fit hyperparameter hardcoded in the code.

[qc]
d_prime_min = 12.0        # minimum assay window, in control-noise units
gradient_flag_pct = 3.0   # top-to-bottom change that flags a plate
gradient_p_max = 0.05     # significance required to flag
```
  (`genedata/settings.toml:1-7`). In this file end-of-line comments are the
  rule: they are the units column, and stay on their key's line even past 90
  characters.

- Parameters that change what a figure shows (colour limits, clip quantiles,
  bin sizes, frame rates) are settings too (`genedata/settings.toml:33-34`,
  `maxwell/settings.toml:14-16, 55-58`). Figure sizes and colours stay in
  `plotting.py`.
- A number fixed by the method (Shared definitions): `0.6745` turns a median
  absolute deviation into a standard deviation (`mea/units.py:82`, which
  lacks the comment). Departure: the references keep a few choices in the
  code (`p < 0.05` in `mea/trajectory.py:68`, the 12 gallery compounds in
  `genedata/run_pipeline.py:89`); here those are settings.
- `config.py` loads the file once, with `tomllib`, into `SETTINGS`
  (`mea/config.py:1-5`). A module names the tables it uses after its
  imports, each in upper case after its table: `QC = SETTINGS["qc"]`,
  `FIT = SETTINGS["fit"]` (`pka/preprocess.py:21-22`, `pka/fitting.py:30`);
  not abbreviations such as `REC`, `SRT`, `UNI` (`mea/units.py:9-12`).
- A setting that a caller (a test, a notebook, a sensitivity check) may need
  to change is a keyword argument that defaults to `None`, read from the
  settings in an `if` block when it is not given:

```python
def select_sta_units(table, evoked=None, n=None, n_tuned=2, min_separation_um=None):
    ...
    if n is None:
        n = UNI["sta_units"]
    if min_separation_um is None:
        min_separation_um = UNI["sta_min_separation_um"]
```
  (`mea/sta.py:43-52`). A test then passes its own value:
  `sta.select_sta_units(table, n=2, min_separation_um=100)`
  (`maxwell/tests/test_mea.py:79`). Never a setting as the default value
  (`flag_pct=QC["gradient_flag_pct"]`, `pka/preprocess.py:50`): it is fixed
  when the module is imported, and a notebook that changes `SETTINGS`
  afterwards does not change it. The `if` block, not the one-line
  `x = ... if x is None else x`.
- Paths are `pathlib.Path` objects, joined with `/`; no `os.path`, no paths
  built by adding strings (neither reference uses `os.path`):
  `sorted(Path(data_dir).glob("Plate_*.json"))` (`pka/loading.py:43`),
  `out.mkdir(parents=True, exist_ok=True)` (`maxwell/run_pipeline.py:25`),
  `out / "units.csv"`.
- Departure: the references keep their data inside the project and take the
  folders as options (`genedata/run_pipeline.py:104-105`). Our data is on
  another drive, so `config.py` is the only place the data root is set
  (Shared definitions), and no other module writes a `D:\` path or lists
  mice.
- Fixed values that are not parameters (column lists, layouts, label maps)
  are module constants in upper case, after the imports and the settings
  names, with a comment above when the value needs one: `RESULT_COLUMNS`
  (`genedata/run_pipeline.py:27-31`), `ASSUMED_ANGLES_DEG`
  (`mea/trajectory.py:10-15`). A constant used only inside its module starts
  with an underscore: `_SETTINGS_PATH` (`mea/config.py:5`).

### Module docstrings

- Every module starts with a docstring, except an empty `__init__.py`. It
  explains the method, and why it is done this way, as a MATLAB header
  explains the pipeline. One line for a simple module, up to about 20 lines
  when the method needs them. Sentence case.

  When one line says it all: `"""Baseline, evoked metrics and direction
  tuning per unit."""` (`mea/responses.py:1`).

  When there is a method and a reason for it:

```python
"""Per-plate control statistics, quality control, and normalization.

All normalization is anchored to each well's own plate, because the three
replicates of every inhibitor sit on different plates:

    percent_activity = 100 * (value - blank_mean) / (neutral_mean - blank_mean)

so 0% = plate background and 100% = uninhibited response on every plate.
Anchoring to both controls removes the background offset and the window gain
(i.e. the assay effect range) at once, the two things that drift from plate
to plate. The per-plate d' (window over pooled control sd) is the
assay-quality index.
"""
```
  (`pka/preprocess.py:1-13`). The reference module
  `sepmap/young_vs_adult/compare.py` does the same for its comparison.

- Layout: the summary on the first line, just after the opening quotes, then
  a blank line and the body, the closing quotes on their own line, then a
  blank line and the imports (pka; mea runs the summary into the body,
  `mea/activity.py:1-3`, not followed).
- Formulas go on their own line, indented four spaces. Papers are cited as
  in Shared rules. A workaround gives its link (the issue link in
  `mea/sorting.py:4-5`).
- Say why the method is the right one, and what the obvious alternative
  would get wrong, in a line or two: "the doses span more than four decades,
  so the log keeps the optimizer on a well-scaled axis"
  (`pka/fitting.py:3-5`). When a method replaced another, say what the other
  gets wrong, not how the code got here.
- No bold, no capitals for emphasis, no markdown other than code in
  backticks.
- A package module has no usage line; its docstring ends with the run script
  that calls it ("Run by run_compare.py."). A run script's usage line is
  indented four spaces and starts with plain `python`
  (`genedata/run_pipeline.py:5`), never the path of one environment's
  interpreter.

### Function docstrings

- Every function and class has a docstring, private and nested ones too; one
  line is enough for a small helper. (pka documents every function; mea
  leaves nine without one, such as `_finish`, `mea/plotting.py:55`.)
- Plain prose, with no `Args:`, `Returns:` or `Parameters` sections (neither
  reference uses the NumPy or the Google format). A summary line in sentence
  case that ends with a period and fits on one line; then, only where the
  names don't make it obvious: what the inputs are (shapes, units), what
  comes back, when a value is NaN, and why the method is the right one.
  Argument names go in backticks (pka). A list of rules is numbered in plain
  text (`pka/preprocess.py:91-100`).

```python
def fit_all(matrix, return_flags=False):
    """Fit every inhibitor in the dose matrix; one results row per inhibitor.

    With `return_flags` also returns a table of the individual points dropped
    by the model pass (inhibitor, plate, concentration).
    """
```
  (`pka/fitting.py:111-116`)

  Shapes and units: `"""Traces for the given channels in microvolts, shape
  (n_channels, n_samples)."""` (`mea/loading.py:106`).

  The reason, when the obvious choice would be wrong:

```python
    """Add inhibition_at_cmax and candidate columns to the fit results.

    Efficacy is the fitted curve evaluated at `c_max`, the highest tested
    concentration: for weak compounds the lower plateau is never observed,
    so 100 - S_inf would claim efficacy that was never measured.
    """
```
  (`pka/selection.py:16-21`)

- A one-line docstring stays on one line:
  `"""Four-parameter Hill equation on log10 concentration."""`
  (`pka/fitting.py:35`).
- Where a function returns NaN, the docstring says when: "With
  `apply_exclusions` (default), excluded wells are NaN; turn it off to fit on
  the raw data." (`pka/preprocess.py:133-134`).
- Departure (Y1): no comment line above `def`. The references put a
  lowercase comment above 109 of their 114 functions, followed by a docstring
  that says the same thing:

```python
# read one plate file into a table with one row per well
def load_plate(path):
    """Read one plate JSON file into a tidy DataFrame with one row per well."""
```
  (`pka/loading.py:15-17`). Here the docstring alone carries it; if the
  comment says something the docstring doesn't, that goes into the
  docstring.

- No blank line between the docstring and the first line of code (none of
  the references' 114 functions has one; MATLAB differs).

### Comments

The shared rules, in the references' form:

- A step, named by a noun phrase or a verb: `# per-plate quality control`
  (`genedata/run_pipeline.py:45`). Don't rewrite a clear comment only to put
  a verb first.
- A reason:

```python
    except RuntimeError:
        # the optimizer failed to converge, so the fit is not usable
        return failed, np.zeros(len(y), dtype=bool)
```
  (`pka/fitting.py:75-77`)

- A decision, and what the alternative would get wrong:

```python
    # double-angle vector: strong for cells tuned to a bar axis, where the two
    # opposite response lobes cancel the ordinary direction vector
    axis_vector = np.sum(rates * np.exp(2j * angles))
```
  (`mea/responses.py:147-149`)

- A dense line: `# design matrix: intercept, soma x, soma y`
  (`mea/trajectory.py:60`). A check states what is expected: `# the DIO
  record must alternate on/off starting with on, by stimulation design`
  (`mea/loading.py:75`).
- Departures: the references have 19 end-of-line comments, such as
  `matplotlib.use("Agg")  # figures go to files, no display needed`
  (`genedata/run_pipeline.py:16`), and comments up to 125 characters; here
  the comment goes on the line above, wrapped at 90.

### Layout

- Departure: `ruff format` sets the layout (line breaks, indentation inside
  brackets, blank lines between functions, spaces, quotes), with the
  settings of `ruff.toml` (line length 90, double quotes). The references are
  formatted by hand: black would change 361 of their 2,918 lines. A formatter
  gives one layout without hand rules, and makes the double quotes of Y2
  mechanical.
- What the formatter leaves to you:
  - blank lines inside a function split its body into blocks, one per step.
    In `main()` every block starts with a comment; elsewhere a block gets one
    when what it does is not clear from the names (`mea/units.py:104, 118`);
  - a long string is split into adjacent literals and a long comment
    rewrapped (ruff check reports both, E501);
  - no backslash continuations: use brackets (departure from
    `pka/plotting.py:150-151` and `maxwell/run_pipeline.py:119-120`);
  - quotes: double, as `ruff format` sets them. Single quotes remain inside
    an f-string's braces (`f"{wells['plate'].nunique()} plates"`) and around
    a text that holds double quotes, where `ruff format` itself writes
    `'He said "x"'`. Formatting with f-strings only; neither reference uses
    `%` or `.format`.
- A function longer than about 60 lines, or with more than five or six
  steps, is split into functions (the references' 114 functions have a
  median of 14.5 lines; five exceed 60).
- Functions come in the order the pipeline uses them, each helper just above
  the first function that calls it: `hill`, `_hill_drop`, `fit_series`,
  `fit_all` (`pka/fitting.py:33-134`). This differs from MATLAB, where local
  functions come last.
- No separator lines inside functions. A module that is not yet split by
  stage may open each group of functions with one banner,
  `# ===== Registration =====` (from the MATLAB half).

### Naming

- `snake_case` for modules, functions and variables, `UPPER_CASE` for
  constants and settings tables, a leading underscore for private helpers
  and constants: `_save`, `_hill_drop` (pka), `_finish`, `_MOSAIC_CELLS`
  (mea).
- A function name says what it returns or does. `load_*` reads files,
  `check_*` validates and fails, `plot_*` draws a figure, `*_table` returns a
  table, and there are `fit_*` and `select_*` functions.
- A quantity with a unit carries the unit as a suffix, in variables, table
  columns and settings keys alike: `duration_s`, `fs_hz`, `min_amplitude_uv`,
  `ac50_max_um`, `refractory_pct`. Counts start with `n_`: `n_bins`.
- Table columns are the names readers see in the CSV: `evoked_hz`,
  `latency_s` (`mea/responses.py:89-94`), `top_to_bottom_pct`
  (`pka/preprocess.py:69`).
- Short names only for a loop variable used over a few lines
  (`for u, tr in trains.items()`, `mea/units.py:23`) or for the standard
  symbol of a formula next to its comment (`A` for the design matrix,
  `mea/trajectory.py:60-61`). Anything used for longer has a readable name
  (departure from `maxwell/run_pipeline.py:30`, which keeps the open file as
  `f` through a 143-line `main()`).

### Functions and data

- A function takes data and returns data. It reads no globals other than
  settings and constants, and writes files only when that is its job (the
  run scripts, a cache function, the `save=` of a figure function).
- Tables are pandas DataFrames with one row per item (well, unit,
  structure) and named columns, built from a list of dicts, and read by
  column name, never by position:

```python
    rows = []
    for plate, p in wells.groupby("plate", sort=True):
        ...
        rows.append({
            "plate": plate,
            "blank_mean": blank.mean(),
            ...
        })
    return pd.DataFrame(rows).set_index("plate")
```
  (`pka/preprocess.py:32-46`)

- A function that adds columns works on a copy and returns it:
  `out = wells.copy()` (`pka/preprocess.py:81`). It changes its input in place
  only when memory requires it, and the docstring says so: "Bandpass then ZCA
  whiten, in place, chunked float32." (`mea/sorting.py:32`).
- Right after loading, a `check_*` function compares the data with what it
  must be (counts, layout, alignment) and fails if it differs; the run script
  prints that the checks passed (`pka/loading.py:50-67`,
  `mea/loading.py:59-65`).
- An excluded item stays in the table, with a boolean column and a reason
  column: `excluded` and `exclude_reason` (`pka/preprocess.py:123-124`),
  `accepted` and `reject_reason` (`mea/units.py:55-57`).
- A value that cannot be determined (the NaN rule of Data safety):

```python
            # a cell silent through every baseline chunk has no noise scale,
            # its effect size is reported as not determined
            dprime = (evoked - base_mean) / base_sd if base_sd > 0 else np.nan
```
  (`mea/responses.py:85-86`, comment wrapped at 90)

- Functions, not classes (neither reference defines a class). A class only
  for state that lives across calls, such as a reader holding an open file
  (`Source`, `mapping/sepmap/volumes/per_mouse.py:98`) or a network
  (`registration/auto_annotation/core.py:188`).
- Simple expressions (Voice), in Python also: no comprehension nested more
  than two levels, no lambda assigned to a name (use `def`), no conditional
  expression spread over several lines (use `if`/`else`).

### Random numbers

- Every draw comes from a local generator,
  `rng = np.random.default_rng(seed)`; never `np.random.seed` or the global
  `np.random` functions.
- An analysis function that draws takes `seed=0` as a keyword argument and
  makes its own generator (`mea/sta.py:30-34`). A figure that jitters points
  uses a generator of its own, separate from the analysis
  (`pka/plotting.py:154, 201`). A test makes its own, with a fixed seed.

### Imports

- At the top, in three groups separated by a blank line: standard library,
  third party, this package (`maxwell/run_pipeline.py:8-14`). Within a group,
  in the order isort gives: plain `import` lines first, then `from` lines,
  each alphabetical (`pka/plotting.py:7-16`).
- Run scripts, tests and notebooks import the package's modules and call
  `module.function`, so every call says where the function lives:

```python
from mea import activity, loading, plotting, responses, sorting, sta, trajectory, units
...
    epochs = loading.load_epochs(f)
```
  (`maxwell/run_pipeline.py:13, 33`; Genedata imports names from the
  package's `__init__.py` instead, not followed)

- Inside the package, a module imports the names it uses from a sibling,
  with absolute imports: `from pka.fitting import hill`
  (`pka/selection.py:11`). mea uses relative ones
  (`from .config import SETTINGS`); not followed.
- Departure: an import inside a function only for a dependency that is not
  installed in every environment that runs the module, with a comment saying
  so (the references import inside functions without such a reason,
  `pka/plotting.py:356`).
- The package never edits `sys.path`. The one exception is a script that
  MATLAB runs by path (`registration/auto_annotation/cli.py:40-42`), with a
  comment saying why.

### Errors

- Check what can go wrong with the data, and fail with a message (Messages).
- Departure: `raise`, not `assert`. The references check their data with
  bare asserts (`pka/loading.py:53-67`, `mea/loading.py:40, 64`), which vanish
  under `python -O` and say nothing when they fail. The form to use is the
  one pka uses for a missing input:

```python
    paths = sorted(Path(data_dir).glob("Plate_*.json"))
    if not paths:
        raise FileNotFoundError(f"No Plate_*.json files found in {data_dir}")
```
  (`pka/loading.py:43-45`). A design check, the comment above the old assert
  turned into the message:

```python
    n_plates = wells["plate"].nunique()
    if n_plates != 15:
        raise ValueError(f"expected 15 plates in 384-well format, found {n_plates}")
```
  (adapted from `pka/loading.py:52-53`)

- Built-in exception types: `FileNotFoundError` for a missing input,
  `ValueError` for data that is not what was expected, `RuntimeError` for a
  step that failed. A custom exception only when a caller catches that case
  by name (`NotReferenceGrid`, raised at `mapping/sepmap/ish/regions.py:113`
  and caught at 196). Package code never calls `sys.exit` or raises
  `SystemExit`; only a run script turns a usage mistake into an exit
  (`parser.error`).
- Catch only what you expect, with the reason in a comment
  (`pka/fitting.py:75-77`, above). A broad `except Exception` only at a
  process boundary that passes the message on, with a comment
  (`registration/auto_annotation/cli.py:149-150`: the GUI shows it).
- Tests keep `assert`, which is how pytest works.

### Type hints

- Departure: public functions of the package have type hints on their
  arguments and return value. Private helpers, nested functions, tests and
  run scripts need none. The docstring still gives the shapes and units,
  which hints cannot.
- Built-in generics and `|`: `list[str]`, `dict[str, float]`, `Path | None`;
  `np.ndarray` for arrays, `pd.DataFrame` for tables:

```python
def load_plates(data_dir: str | Path) -> pd.DataFrame:
    """Load all Plate_*.json files in `data_dir` into one tidy DataFrame."""
```
  (`pka/loading.py:41-42`, with hints added)

- On Python 3.12 hints are evaluated when the `def` runs, so a hint that
  names a type that is not imported breaks the import of the module.

### Figures

- One plotting module holds the palette, the colormaps, `set_style()`, the
  save function and the figure functions. A figure function is
  `plot_<what>(data, ..., save=None)` and returns the figure; one that draws
  a single panel takes `ax=None` and returns the axes
  (`pka/plotting.py:370-408, 412-441`).
- Figures are saved through one function, and a saved figure is closed so a
  long run does not accumulate them:

```python
def _finish(fig, save, tight=True):
    if tight:
        fig.tight_layout()
    if save is not None:
        # pipeline use: write and close; notebooks pass save=None and display
        fig.savefig(save, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return fig
```
  (`mea/plotting.py:55-62`; pka's `_save` does not close the figure, not
  followed). The package has five copies of a `save_figure` instead
  (`diagnostics.py:50`; in `young_vs_adult/`, `closeup.py:117`,
  `region_groups.py:91`, `region_plot.py:104`, `compare.py:66`), none of which
  closes the figure. All but the first write an EPS next to the PNG, with the
  image layers rasterised, and all write `<name>_new.png` when the PNG is
  open in a viewer. They are to be
  merged into one save function of the kind above, keeping those two
  behaviours (kind C).
- New figures use the dpi of the save function. Existing figures keep their
  own, since changing it changes every output.
- The palette and colormaps are named constants at the top of the plotting
  module, with a comment above giving the role when it is not obvious
  (`mea/plotting.py:10-28`; pka puts the role at the end of the line, not
  followed). The colours themselves are in the shared Figures section.
- `set_style()` sets the matplotlib rcParams once
  (`mea/plotting.py:32-51`); the run script calls it, and so does the first
  code cell of every notebook.
- Titles and axis labels in lowercase, names and symbols written as they
  are, the unit in parentheses: `"time from onset (s)"`,
  `"peak amplitude (uV)"` (`mea/plotting.py:301, 156`). A computed count
  (Figures, shared): `f"recorded electrodes (n = {len(mapping)})"`
  (`mea/plotting.py:76`).
- New figure files are numbered in run order: `01_raw_plates.png`
  (`genedata/run_pipeline.py:64`). Existing file names stay.

### Tests

- pytest, in `tests/` next to the package, run from the pipeline folder with
  `python -m pytest tests` (`maxwell/README.md:45`). Plain `pytest` puts
  `tests/`, not the pipeline folder, on the import path, and cannot import
  the package (`genedata/README.md:43` says `pytest`; on a copy of that
  project it stops with a collection error).
- Departure: one file per sub-package or topic, `test_<topic>.py`, as the
  MATLAB tests are organised (each reference has a single file).
- A test is a function named for the expected behaviour
  (`test_flat_curve_is_classified_inactive`, `genedata/tests/test_pka.py:55`),
  with a one-line docstring saying what it checks.
- A test has a known answer, on synthetic data: plant something, get it
  back. A tolerance that is not obvious gets its reason:

```python
def test_snippet_average_recovers_planted_waveform():
    """A planted spike shape is recovered on its channel; a clean channel stays flat."""
    rng = np.random.default_rng(1)
    traces = rng.normal(0, 1, size=(200_000, 8)).astype(np.float32)

    # spike-like dip
    wave = -30 * np.exp(-0.5 * ((np.arange(41) - 20) / 3.0) ** 2)
    times = rng.choice(np.arange(1000, 199_000, 400), size=300, replace=False)
    for s in np.sort(times):
        traces[s - 20 : s + 21, 3] += wave.astype(np.float32)
    out = sta.snippet_average(traces, np.sort(times), half=20)

    # planted channel recovered within noise (SEM about 1/sqrt(300))
    assert np.abs(out[:, 3] - wave).max() < 0.5
    assert np.abs(out[:, 0]).max() < 0.5
```
  (`maxwell/tests/test_mea.py:44-55`, its comment above `def` made the
  docstring and its end-of-line comment moved above)

- Include the case that must come out negative: a flat curve called
  inactive (`genedata/tests/test_pka.py:54-59`), an untouched channel that
  stays flat, a corrupted input refused (`maxwell/tests/test_mea.py:35-41`).
  A test that could only pass checks nothing.
- Errors are tested with `pytest.raises(ValueError)`.
- Departure: a test that needs real data is skipped, not failed, when the
  data isn't connected (as the MATLAB tests do):
  `@pytest.mark.skipif(not Path(config.DATA).is_dir(), reason="data not connected")`
  (`config.DATA` is the data root, a string).
- Tests write only into pytest's `tmp_path`, never into the output folders.
  A `tests/README.md` lists the tests (file, what it checks, data).

### README of a pipeline

Every pipeline folder has a README with the sections of the references
(`maxwell/README.md:1-121`; Genedata has the same, without the limitations):

```
# <Pipeline name>

What it does in two to four lines: the question, the data, what comes out.

## Layout
## Setup and run
## Assumptions and preprocessing
## Analyses
## Outputs
## Known limitations
```

- The bullets give reasons and numbers: "The three replicates of each
  inhibitor are on three different plates, so each well is normalized to its
  own plate's controls" (`genedata/README.md:51-55`).
- Limitations are stated plainly: "Three template-similar unit pairs are
  most likely single cells split in two by the sorter (evidence in notebook
  02). They are documented rather than merged, so the accepted-unit count is
  slightly inflated; no conclusion depends on it."
  (`maxwell/README.md:115-118`).
- Departure: the README ends with its known limitations.

### Environments

- One requirements file per environment, in `tools/`, pinned to the versions
  the results were produced with, with a comment at the top saying which
  Python and how to create the environment, and the packages in commented
  groups (as `genedata/environment.yml:4-14`; ours are `venv`s from the
  Anaconda Python 3.12.7, not conda environments):

```
# tools\venv_atlas: Python 3.12 (Anaconda). Create it from the repository root with
#   python -m venv tools\venv_atlas
#   tools\venv_atlas\Scripts\python -m pip install -r tools\requirements_atlas.txt

# scientific stack (the versions the results were produced with)
numpy==2.5.2
```

- The README of each pipeline says which environment runs which script:
  `tools\venv_atlas` for the Python route, `tools\venv_flat` for the
  flatmaps, the automatic annotation's own environment for its engine.
- pytest and ruff live in `tools\venv_dev`, so the analysis environments
  never change.

### Linting and formatting

- `ruff.toml` at the repository root (departure: neither reference lints):

```toml
line-length = 90
target-version = "py312"
extend-exclude = ["third_party", "archive"]

[lint]
select = ["E", "W", "F"]

[format]
quote-style = "double"
```

- `ruff check` catches long lines (E501), imports below other code (E402;
  ruff 0.16 does not count a `matplotlib.use` call as code), several
  statements on one line (E701-E703), a lambda assigned to a name (E731),
  unused imports and undefined names (F). `ruff format` sets the layout.
  Both run from `tools\venv_dev`:

```
tools\venv_dev\Scripts\python -m ruff format <files>
tools\venv_dev\Scripts\python -m ruff check <files>
```

## The style pass

How existing code is brought in line with this guide. The pass works file
by file, from this guide and its reference examples.

### Three kinds of change

Each kind goes into commits of its own, so a commit is either checked for
code identity or by a rerun, never half of each.

**Kind A: comments and layout.** The code must stay identical, which the
identity check proves without running anything:

- MATLAB: `[T, ok] = sep_check_code_identity(new_dir, ref_dir)` compares the
  parse trees (`mtree`) of every `.m` file, comments and layout removed.
- Python: `python tools\check_code_identity.py NEW_DIR REF_DIR` compares the
  `ast` of every `.py` file, comments and docstrings removed.

Every file with kind A changes only must come out `same code`; whatever the
check reports is the answer. Kind A is: comments (add, rewrite, delete, move
an end-of-line comment above), headers, help blocks and docstrings (fold a
comment above `def` into the docstring), section titles, blank lines,
indentation, spaces, line breaks (`...` in MATLAB, brackets in Python),
statements split onto their own lines (`a = 1; b = 2`); in Python also
quotes (Y2), `ruff format`, `1.` written `1.0`, a long string split into
adjacent literals (run the check for an f-string).

One exception: a docstring used at run time is program output, so editing it
is kind B. Search for `__doc__` first: `registration/auto_annotation/cli.py:157`
prints its whole docstring as its usage, and
`atlas/build_demba_atlas.py:241` the second-to-last line of its own;
`landmark_refine/cli.py:50` and `landmark_refine/serve.py:109` print theirs.

**Kind B: small code changes inside one file.** Only where the assignment
allows it, verified by rerunning and comparing the outputs with
`sep_compare_outputs` or `compare_outputs.py` (`tools/README.md`): contents
identical, and figures identical or `same render` as `tools/README.md`
defines it.

- rename a local variable, every use (never a name saved to a file, a
  struct field or dict key, a column name, a keyword callers pass, or a
  setting `sep_run_driver_copy` addresses);
- split a dense expression into named intermediates, with the same
  arithmetic in the same order, the same reductions, NaN handling and dtypes
  (the `float16` and `float32` casts stay where they are);
- a conditional expression or a lambda assigned to a name becomes
  `if`/`else` or a `def`;
- remove an unused import (check it has no side effect) or an unused
  variable;
- an f-string without placeholders becomes a plain string;
- add type hints to public functions (then import every changed module);
- `assert` becomes `raise`, the comment above it the message; `SystemExit` in
  package code becomes a built-in exception (check first what catches the
  old type: `registration/auto_annotation/cli.py` hands every message to
  MATLAB);
- MATLAB: a driver's `clear all`, `close all` and `clc` become
  `clear; clc; close all;`;
- printed text changes only if nothing reads it.

Figures keep the same drawing calls, in the same order, with the same
properties. Random state: never add, remove or reorder a call that draws
from a generator (`rand`, `randn`, `randperm`, `kmeans`, a shuffle, any
`rng.` call). In `mapping/sepmap/ish/panel_test.py`
one `rng` (line 196) feeds up to five permutation tests in turn (lines 203,
215, 225, 240), so reordering them changes every p-value; the figure has its
own generator (line 257).

When in doubt, add a comment instead. A possible bug is not fixed in the
pass: it goes on the bug list of the plan. A rename this guide would want
goes on a list for the owner.

**Kind C: structural changes.** The plan gives the style pass
these too (`matplotlib.use` and the module constants are left "for the
style pass" in its Progress). They come after the kind A and B commits of
the files they touch, one change per commit, each checked by rerunning
every step that uses the changed code and comparing as for kind B:

- `matplotlib.use` out of the modules and into the run scripts;
- constants into `settings.toml` (`SMOOTH`, the 100-voxel minimum), and a
  module's settings read through its table
  (`YOUNG_VS_ADULT = SETTINGS["young_vs_adult"]`) rather than key by key;
- `os.path` to `pathlib` (a `Path` has no string methods:
  `path.replace(".png", "_new.png")` in `save_figure` would rename a file);
- changing a signature (settings as keyword arguments, a `seed` argument);
- merging helpers defined more than once: the five `save_figure` (Figures,
  above), `per_unit` (`mapping/sepmap/volumes/cohort.py:190` and
  `young_vs_adult/region_plot.py:195`), `fold_n` (`compare.py`, `closeup.py`,
  `video.py`);
- one Python plotting module and one MATLAB palette function in `common/`,
  and moving the drawing out of the modules that compute;
- reordering functions, splitting a function longer than about 60 lines
  (step 6 split only the longest);
- tuples read by position (`r[11]`, `young_vs_adult/region_plot.py:373`)
  into named columns;
- the environment variable `V2_READINGS` into a setting or an option.

**Not part of the style pass**: moving files (step 4); the cohort lists
into one table (step 8, Y4); fixing a possible bug (step 8); figure dpi and
output file names, which change every output (the owner decides).

### Checks for each file

MATLAB:

1. header before any code; a driver has `%% Settings` and `%% Run`, a
   function its help block, every local function its one or two lines;
2. comments lowercase, on their own line, at most two lines (except a
   justified setting or constant, or a workaround), a blank line before each
   block, `%#ok` pragmas left at the end of their line;
3. `%%` sections in scripts and in function bodies over about 60 lines; one
   `% ===== Local functions =====` line;
4. lines under about 90 characters, except a string literal that does not
   fit;
5. `sep_check_code_identity`: `same code`;
6. `checkcode`: no message that the file did not have before.

Python:

1. module docstring: summary on the first line, method and reason, no
   history, no bold or capitals, the run script that calls it (or the usage
   line with plain `python`);
2. every function and class has a docstring, with no comment line above
   `def` and no blank line after the docstring;
3. comments lowercase, on their own line, at most two lines (except a
   justified setting or constant, or a workaround), a blank line before each
   block, one per stage in `main()`;
4. `ruff format --check` passes; no backslash continuation;
5. `ruff check` shows no warning the file did not have before;
6. `check_code_identity.py`: `same code`.

The kind A commits are listed in `.git-blame-ignore-revs`, so `git blame`
skips them; a kind A commit therefore holds nothing else, not even a new
file. Kind B and C commits change code, so blame keeps them. At the end, a
search that no file credits an author other than Giulio.
