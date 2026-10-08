# Code style

Reference for all code in this repository, MATLAB and Python. New code
follows it; the code that existed was brought in line by the style pass,
step 7 of [history/REFACTOR_PLAN.md](history/REFACTOR_PLAN.md).

Shared rules come first. The MATLAB half is the imaging repository's guide
(`D:\dendrites\code\docs\STYLE.md`) adapted (Y6); the Python half quotes two
reference projects wherever a rule comes from them, so this file is enough to
work from (Y7). Where the earlier conventions differ (those of the scripts in
`adult_matlab/`, which keep them until they retire), the imaging guide wins
unless a project rule applies (table at the end of the MATLAB half).

## Reference examples

| file | shows |
|---|---|
| `group_comparison/run_collect_by_group.m` | a MATLAB driver: header, settings, one call |
| `group_comparison/pipeline/collect_by_group.m` | the function it calls: help block, settings unpacked, `%%` steps |
| `mapping/run_compare.py` | a Python run script: header with the run order, `main()` |
| `mapping/sepmap/young_vs_adult/compare.py` | a package module: docstring with the method, constants, docstrings |

They are in final style, with one departure: `compare.py`, like every
analysis module of the Python route, still draws its own figures, which new
code puts in a plotting module (Python, Pipelines and run scripts).

## Shared rules

### Voice and messages

The code reads as if written by a careful scientist for a colleague: plain
and short, no capitals for emphasis, no "NOTE:", "IMPORTANT" or "!!!", no
bold, no stock words ("comprehensive", "robust", "elegantly", "leverages",
"ensures", "Here we"). No history ("previously", "fixed a bug where"): it
belongs in the git log, and in `docs/` when a reader needs it.

- No options, layers or defensive branches nobody needs. One step per line,
  named intermediate variables rather than nested calls, a plain loop rather
  than a clever vectorised line.
- Numbers where there are numbers: a threshold and why, a count, a size.
- Text the program prints, writes or draws is code: changing it is a code change.
- An error says what was found, what was expected, the path, and the fix when
  there is one. Progress is one lowercase line per step, sub-items indented
  by two spaces; a case skipped while the run goes on is a warning.

### Comments

- Short and lowercase, verb first where it reads naturally, else a noun phrase:
  `% load traces`, `# per-plate quality control`. A capital only for a proper
  name, an acronym or an identifier (`# Gaussian sigma in 20 um voxels ...`).
- On their own line, just before the code, one per step (in a loop, one per
  sub-step); never at the end of a line, except tool pragmas (`%#ok<NANMEAN>`,
  `# noqa: E402`) and the units column of `settings.toml`.
- The reason only when it isn't obvious, in one line (`# the atlas in flat
  grey under the data, so no data reads as grey`). For an important decision
  (a threshold, a window, an exclusion), why, and what the alternative would
  get wrong:

```matlab
% weight by frame count: files differ in length after bad-frame removal,
% and a plain mean of means would favour the short ones
```

- Don't restate the code (`% loop over files`), but comment a dense line
  (`% rescale to 0-1 and clip values outside the limits`), and say what a
  check expects (`# the DIO record must alternate on/off starting with on`).
- One or two lines, wrapped at 90 characters; more only to justify a setting
  or constant, or to explain a workaround. A longer account of a method goes
  into the help block or docstring. Papers by author and year, in one line,
  "(Hill 2011, as in SpikeInterface)"; what a paper shows and what it means
  here goes into `docs/SCIENTIFIC_CONTEXT.md` (L7), not into the code.
- No commented-out code in new code. Old commented-out code stays until the
  owner decides, with a comment saying why it is off, or that it is not recorded.

### Data safety

- Raw data is read only: nothing writes under `S:\` or to the raw `.czi`
  files. Outputs go under the data root only (`get_paths().data`,
  `config.DATA`); no code writes a `D:\` path of its own.
- `SEP_DATA_ROOT` moves the whole data tree for checks; `get_paths.m` and
  `sepmap/config.py` refuse the production data to a copy of the code, and a
  data root inside the G: snapshot. A check never writes into production data.
- Code deletes no files; a step overwrites only its own outputs.
- Excluded data, and a value that cannot be determined, is NaN, never zero,
  with the reason in a comment; NaN is left out of averages.
- An expensive result is cached, with a switch that recomputes it. Loaded or
  cached data is checked against what it must hold (its mice against the mice
  asked for) and the run stops if it differs; a cached table that was sorted
  is put back in its list's order on loading (`ismember`).

### Shared definitions

Each is defined in one place; never copy a value into a script.

- Paths: `get_paths.m` and `mapping/sepmap/config.py`, worked out from where
  the code sits. The MATLAB path: `sep_setup_paths`, once per session.
- The cohort: one table for both languages, `common/cohort.csv` (Y4), read by
  `get_cohort` (`get_cohort_spec` for the mice a stack holds) and by
  `sepmap.config.mapping_mice`, which `COHORTS` in `sepmap/volumes/cohort.py`
  groups. New code takes its mice from these, never from a list of its own.
- The atlases: `get_atlas(key)`, `cohort_atlas_key`, `get_atlas_crop`. Only
  `get_atlas` adds an atlas folder to the path, the one a run needs.
- Parameters: a MATLAB driver's `%% Settings`; in Python `mapping/settings.toml`.
  A number someone could reasonably choose differently is a parameter. One fixed
  by the method stays in the code, with a comment saying what it is (`0.6745`
  turns a median absolute deviation into an SD); a figure's layout stays too.
- Colours: `common/sep_palette.m` in MATLAB (`sep_palette('nano')`,
  `sep_palette('difference')`; `help sep_palette` lists the names) and
  `mapping/sepmap/plotting.py` in Python (`GROUP_COLOURS`, `NO_DATA_GREY`,
  `hot_cut()`), with the values under Figures. A new colour goes into the
  palette of its language, never into a script.

### Figures

| what | colormap |
|---|---|
| intensity (nano, autofluorescence, a reading) | hot: `sep_palette('intensity')`, `plotting.hot_cut()` |
| atlas, anatomy, raw images | gray: `sep_palette('anatomy')` |
| a difference (young minus adult, control against experimental, t maps) | blue-red, symmetric limits, zero in the middle: `sep_palette('difference')`, `RdBu_r` |
| a signed position within a brain (zref) | `PuOr_r`, purple low, so blue-red keeps meaning a difference |
| counts of brains (n maps) | magma |
| coverage lines | plasma |

- Never parula; jet or turbo only when asked for, into a subfolder of their
  own (`run_closeup --cmap`), and a difference stays blue-red.
- No data is flat grey `#bfbfbf` (`plotting.NO_DATA_GREY`), which no data
  colormap produces. In the young-against-adult maps hot stops at 0.82 of its
  range, never reaching white (`plotting.hot_cut()`).
- Bars of a value per structure across mice: the height is the value, the
  colour its reliability (t, clamped) in grey, darker for more reliable, never
  pure white (`sep_palette('bars')`). Bars that compare categories take
  palette colours.
- Palette (`sep_palette`; the groups also in `plotting.py`): nano
  `[0.95 0.55 0.10]`, autofluorescence `[0.95 0.85 0.20]`, their per-mouse dots
  `[0.65 0.30 0.00]` and `[0.70 0.60 0.00]`, lines joining paired mice
  `[0.6 0.6 0.6]`; groups young `#c0392b`, naive `#555555`, rws `#9a9a9a`; the
  plasticity comparison's control and experimental groups, MATLAB's default
  blue and orange (`'control'`, `'experimental'`, darker for their means); the
  receptor subunits dark blue `#1f3b73` (`plotting.DARK_BLUE`). New colours
  extend this orange, yellow, blue, red and grey family: no green, pink or
  purple outside the colormaps above.
- In Python the channels are `plotting.NANO`, `AUTO`, `NANO_DOT`, `AUTO_DOT`
  and `PAIR_LINE`. The ISH figures add the groups of divisions of a scatter of
  structures (`DIVISION_GROUP`, `DIVISION_GROUP_COLOURS`: cortex `#e07b00`,
  hippocampal formation `#c0392b`, thalamus `#3a6db5`, other grey matter
  `#555555`), the gene sets (`SET_COLOURS`: subunits dark blue, localisation
  red, other postsynaptic `#555555`, presynaptic `#9a9a9a`, GABAergic markers
  `#7f9cc9`, glia `#c8c8c8`), the 95% band of a null distribution behind the
  data (`NULL_BAND`, `#c9d6ea`), the density step of a variance budget
  (`DENSITY_BLUE`, `#7f9cc9`), and the green (SEP) channel in blue so the three
  channels stay apart (`SEP` `#3a6db5`, its per-mouse dots `SEP_DOT` `#24427f`,
  what is left of it once autofluorescence is regressed out `SEP_REMAINDER`
  `#8fb3e0`). `plotting.bars_grey(t, t_max)` is
  `sep_palette('bars')`: grey 0.78 at t = 0 to black at t_max.
- Scatter plots of many structures: 35-point dots, no edge, alpha 0.85.
- Names beside dots never sit on each other: `ish/plotting.spread_labels` (a
  scatter) and `spread_positions` (a column of names) move them apart and join
  a moved name to its dot. A spatial p is written `p = 0.012` (three decimals),
  or `p ≤ 0.0001` when no surrogate reached it.
- Counts in titles are computed, never typed. Coronal planes are drawn dorsal
  up, (DV, ML), never transposed. White background, except image panels.
- A file name carries every setting that changes what the file holds
  (`_smooth5` or `_nosmooth`, the channel, `_P20`). Existing names stay.

## MATLAB

### Drivers and functions

- A pipeline is one folder: its drivers (`run_*.m`) at the top, the functions
  they call in `pipeline/`, checks in `qc/`; shared functions in `common/`,
  atlas functions in `atlas/`. Tools run by hand (`make_ordering_volume`,
  `remap_control_points`) sit at the top beside the drivers, and a
  self-contained component gets a subfolder of its own
  (`registration/annotation_gui/`, `registration/auto_annotation/`).
- A driver is a script: the header, `clear; clc; close all;` on one line,
  `%% Settings`, `%% Run`; everything else is a function. No path setup:
  `sep_setup_paths` runs once per session, and before each stage of the
  detached runner (a driver copy cannot find the code with `mfilename`).
- `%% Settings` holds everything the user is expected to change, first
  `paths = get_paths();`, each with a comment above giving the allowed values
  where useful (`% groups to collect ('rws', 'naive', 'behavior', 'young')`).
- Settings are lowercase `snake_case`, unlike the imaging guide, because
  `sep_run_driver_copy` and `docs/history/production_settings.md` use their
  names.
- `%% Run` copies the settings into `run_settings` under the same names
  (`run_settings.mousetypes_list = mousetypes_list;`) and makes one call,
  `collect_by_group(run_settings);`. A driver computes nothing else.
- The function a driver calls takes `run_settings`, names the driver in its
  help block, and first unpacks the settings, one per line
  (`mousetypes_list = run_settings.mousetypes_list;`).
- A long function is a short main part of `%%` steps calling local functions
  with explicit inputs and outputs: no nested functions, no globals, no
  `evalin` or `assignin`.
- `load` into a struct, never the workspace: `S_nano = load(file, 'nano_4d');`.
  Other options come as a struct, a missing field defaulted at the top.
- Local functions come last, after `% ===== Local functions =====` (not a
  `%%` section: they are not steps); a long file has one such line per group,
  `% ===== Local functions: loading and alignment =====`.

### Headers

Every file has one, local functions included: sentence case, short.

- A driver's comes before any code, so `help` shows it: the pipeline and its
  steps, what this one does in general, then briefly the data it is set up for:

```matlab
%% run_collect_by_group
% ===== Collect the registered volumes of each group =====
%
% Plasticity comparison, step 1 of 3:
%   1. run_collect_by_group    stack each group's registered volumes  <- this script
%   2. run_normalise_groups    bring the mice of a group onto one intensity scale
%   3. run_group_differences   left-right differences, control against experimental group
%
% Stacks the registered nano volumes of a group's mice into one 4D array ...
%
% Setup: the young cohort, P20 brains only, for the comparison with the ...
```

- A function's is the usual help block, no essays: the name in capitals and a
  summary line, the call syntax, inputs, outputs, options with defaults:
  `%COLLECT_BY_GROUP  Stack the registered volumes of each group's mice.`
- A local function's is one or two lines under the signature, then a blank
  line: `% The cohort's folders, and its stack of the channel with the
  selected mice only.`

### Layout and naming

- A blank line before each comment block, at the start of a loop or `try`
  body with several steps, and after the help block.
- `%%` sections with a short capitalised title, not numbered
  (`%% Collect each group`), in every script and in function bodies over about
  60 lines; short functions have none. One statement per line:
  `if ~exist(out_dir, 'dir')`, `mkdir(out_dir);` and `end` on three lines.
- Spaces around `=` and after commas. Lines under about 90 characters,
  continued with `...` and an indent of four; the pieces of a message joined
  in `[...]` line up under the first. A string literal that does not fit stays
  whole on its line.
- Text in single quotes (`'lightsuite'`); `"..."` is a string, another type.
- `snake_case`; `camelCase` only for LightSuite's names (`nanoVol`). `n_*` for
  counts, `*_v` for one value per region (`nano_mean_v`), `S_<name>` for a
  loaded struct. Switches start with what they switch: `produce_*`, `show_*`,
  `save_*`, `do_*`, `use_*`, `compute_*`, `force_recompute_*`.

### Errors, messages, figures

- An `error` starts with the name of the function or driver, then the shared
  rule: `error(['sep_test_path: sep_setup_paths resolves to %s, not to this
  code folder (%s). Run restoredefaultpath, ...'], ...)`.
- A check on three lines: `if ~isequal(nano_mice, auto_mice)`, `error(...);`,
  `end`. A string option goes through `switch`, with
  `otherwise error('Unknown agg_method: %s', agg_method)`.
- Progress with `fprintf`, a skipped case with `warning`:
  `fprintf('  %d mice paired across channels.\n', n_mice)`.
- In code that runs unattended, `figure('Visible', 'off', 'Color', 'w',
  'Units', 'Normalized', 'Position', [0 0 1 1])`; each figure saved as `.fig`
  (`saveas`) and `.png` (`exportgraphics`, 300 dpi).
- An underscore in a title is a TeX subscript: `'Interpreter', 'none'`, or
  `strrep(tag, '_', ' ')`.
- Grouped horizontal bars: `barh(Y, 'grouped')`, then `EdgeColor = 'none'` on
  each series. A grid of panels has three columns:
  `n_rows_grid = ceil(n_panels / n_cols_grid)`.

### Where the earlier conventions differ

| topic | earlier (`adult_matlab/`) | here |
|---|---|---|
| driver | `clear all`, `close all`, `clc` on three lines; `% /// Pipeline script #N ///` header; `%% User-defined parameters` | imaging: `clear; clc; close all;`, banner header, `%% Settings` |
| settings | `snake_case` (imaging: `UPPER_CASE`), units at the line's end | earlier names (project rule); units in a comment above (imaging) |
| comments | sentence case, several lines for a why | imaging: lowercase, one or two lines, more for a workaround or a setting's value |
| sections | numbered titles; `%% Local function: <purpose>` above each | imaging: unnumbered titles; one separator line |
| short `if` | on one line | imaging: three lines |
| numbers in the code | every configurable value at the top | earlier, now a shared rule |
| path setup | `%% Add paths` in each driver (imaging: `setup_paths`) | neither: `sep_setup_paths` once per session |
| figures | `.png` and `.fig`; hot, gray, blue-red, grey bars; RGB triplets in each script (imaging: `save_fig`, `condition_color`) | earlier (project rules), plus `PuOr_r`, magma, plasma; the colours from one palette function, `common/sep_palette.m` |
| replaced files | overwritten (imaging: moved to `superseded/`) | earlier: no `set_aside_file` here |

## Python

The reference projects, by the same author, are the Genedata exercise
(package `pka/`) and the MaxWell exercise (`mea/`), read only in
`S:\ElboustaniLab\#SHARE\Documents\Giulio\`; most examples below are quoted
from them. Where they differ from this guide, the guide wins. Every Python
file follows it, `tools/` included (Y5).

### Pipelines and run scripts

- One folder: the entry points (`run_*.py`), `settings.toml`, `README.md` and
  `tests/` at the top, next to the package they call, never inside it.
- One module per stage, named after it (`loading.py`, `preprocess.py`,
  `fitting.py`); a route with several analyses has one sub-package per
  analysis, laid out the same way (`sepmap/volumes/`, `young_vs_adult/`).
- `config.py` loads the settings and holds the paths. `plotting.py` holds the
  palette, the style, the save function and every figure function (one per
  sub-package with figures); the modules that compute draw nothing. In the
  Python route, `mapping/sepmap/plotting.py` holds the palette, the
  colormaps, `save_figure` and the drawing several modules share (the coronal
  frame); the figure functions of each analysis are still in its module, and
  new figure code goes into a plotting module.
- Only run scripts have an `if __name__ == "__main__":` block. `__init__.py`
  is empty or a one-line docstring.
- Notebooks are optional, numbered (`01_exploration.ipynb`), and only call the
  package, so no result exists only there: first a markdown cell (title, what
  it shows), then a code cell that imports, calls `set_style()` and loads.
- One `run_<step>.py` per step (L1). Its header is the module docstring: the
  summary line, the run order with this step marked, what it does and writes,
  the usage line, its options (`mapping/run_compare.py`).
- `main()` first prints the settings in force (`config.print_settings`), then
  holds only glue (calls, a loop over items, file names, counts), one stage
  per block: an unnumbered comment, the calls, and one printed line with the
  counts that show it went right:

```python
    # unit table and quality control
    table = units.unit_table(trains, templates, mapping, duration_s)
    table = units.apply_unit_qc(table)
    table.to_csv(out / "units.csv", index=False)
    n_ok = int(table["accepted"].sum())
    print(f"unit QC: {n_ok} accepted of {len(table)}, table written to units.csv")
```

- Run options (which mice, whether to recompute) are `argparse` options,
  parsed under `if __name__ == "__main__":` and passed to `main()`, so a test
  can call it; analysis parameters never are. The cache switch is named for
  what it redoes (`--resort`); the cache file name is a setting.
- Progress with `print`, not `logging`, with `flush=True` inside long loops.
- The figure backend is set only in run scripts (`matplotlib.use("Agg")`, first
  under `if __name__ == "__main__":`), never in a module, so notebooks still
  show their figures and importing a module changes no global state.

### Settings and paths

- Every analysis parameter is in `settings.toml`, which says at its top what
  it is: a table per stage, a key per line, the unit (and reason) in an
  end-of-line comment aligned in a column, even past 90 characters:
  `d_prime_min = 12.0        # minimum assay window, in control-noise units`.
- What changes a figure's content (colour limits, clip quantiles, bin sizes)
  is a setting, and so is a choice such as `p < 0.05`; figure sizes and
  colours stay in `plotting.py`.
- `config.py` loads the file once (`tomllib`) into `SETTINGS`; a module names
  each table it uses after its imports, upper case after the table:
  `QC = SETTINGS["qc"]`, never an abbreviation such as `REC`.
- A setting a caller may change is a keyword argument defaulting to `None`,
  read in an `if` block (`if n is None:`, `n = UNITS["sta_units"]`); never a
  setting as the default value, which is fixed at import.
- Paths are `pathlib.Path` joined with `/` (`out / "units.csv"`), never
  `os.path` or strings added together.

### Docstrings

- Every module but an empty `__init__.py` has one, explaining the method and
  why, as a MATLAB header explains the pipeline, in sentence case: one line
  (`"""Baseline, evoked metrics and direction tuning per unit."""`) up to
  about 20 (`compare.py:1-22`). The summary on the first line, a blank line,
  the body, the closing quotes on their own line.
- It says why the method is right and what the obvious alternative gets
  wrong ("the doses span more than four decades, so the log keeps the
  optimizer on a well-scaled axis"), never how the code got here. Formulas on
  their own line, indented four; a workaround gives its link; no markdown but
  backticks. A package module's ends with the run script that calls it ("Run
  by run_compare.py."); a run script's has a usage line, indented four, with
  plain `python`.
- Every function and class has one, private and nested ones too; a small
  helper's is one line, kept on one line. Plain prose, no `Args:`, `Returns:`
  or `Parameters` sections: a one-line summary in sentence case ending with a
  period, then only what the names don't say (shapes and units, what comes
  back, when a value is NaN, why the method is right), argument names in
  backticks, a list of rules numbered:

```python
def fit_all(matrix, return_flags=False):
    """Fit every inhibitor in the dose matrix; one results row per inhibitor.

    With `return_flags` also returns a table of the individual points dropped
    by the model pass (inhibitor, plate, concentration).
    """
```

- No comment line above `def` (Y1), and no blank line after the docstring.

### Layout and naming

- `ruff format` sets the layout (`ruff.toml`: 90 characters, double quotes;
  single only where it keeps them, as inside an f-string's braces). Left to
  you: blank lines splitting a body into steps, each opened by a comment in
  `main()`, elsewhere only when the names leave it unclear; long strings split
  into adjacent literals; brackets, never a backslash, to continue a line;
  f-strings, never `%` or `.format`.
- A function over about 60 lines, or five or six steps, is split. Functions
  come in the order the pipeline uses them, each helper just above its first
  caller. No separator lines inside functions; a module not yet split by
  stage may open each group of functions with `# ===== Registration =====`.
- `snake_case`; `UPPER_CASE` for settings tables and for module constants
  that are not parameters (column lists, layouts, label maps), after the
  settings tables; a leading underscore for private names (`_hill_drop`).
- A function name says what it does: `load_*` reads files, `check_*`
  validates and fails, `plot_*` draws, `*_table` returns a table.
- A unit is a suffix, in variables, columns and settings keys alike
  (`duration_s`, `min_amplitude_uv`); counts start with `n_`. Table columns
  are the names readers see in the CSV.
- Short names only for a loop variable over a few lines
  (`for u, tr in trains.items()`) or a formula's symbol beside its comment.

### Functions and data

- A function takes data and returns data: no globals but settings and
  constants, no files written unless that is its job (a run script, a cache,
  a figure's `save=`).
- Tables are DataFrames, one row per item, built from a list of dicts and
  read by column name, never by position. A function that adds columns works
  on a copy (`out = wells.copy()`); in place only when memory requires it,
  and its docstring says so.
- Right after loading, a `check_*` function compares the data with what it
  must be (counts, layout, alignment) and fails if it differs; the run script
  prints that the checks passed. An excluded item stays in the table, with a
  boolean and a reason column (`excluded`, `exclude_reason`).
- Functions, not classes, but for state that lives across calls (a reader
  holding an open file, `Source` in `sepmap/volumes/per_mouse.py`). No
  comprehension nested more than two levels, no lambda assigned to a name, no
  conditional expression over several lines.
- Random draws come from a local `rng = np.random.default_rng(seed)`, never
  `np.random.seed` or the global functions. An analysis function that draws
  takes `seed=0`; a figure that jitters points has its own generator, a test
  one with a fixed seed.
- Imports at the top, in three groups (standard library, third party, this
  package), in isort order. Run scripts, tests and notebooks call
  `module.function` (`loading.load_epochs(f)`); within the package, absolute
  imports from a sibling (`from pka.fitting import hill`). An import inside a
  function only for a dependency missing from some environment, with a
  comment; no `sys.path` edit but in a script MATLAB runs by path.

### Errors and type hints

- `raise`, never `assert` (gone under `python -O`), with the shared message:
  `raise FileNotFoundError(f"No Plate_*.json files found in {data_dir}")`;
  `ValueError` for unexpected data, `RuntimeError` for a failed step, a
  custom exception only when a caller catches it by name. Package code never
  exits (`sys.exit`, `SystemExit`); a run script uses `parser.error`.
- Catch only what you expect, with the reason
  (`# the optimizer failed to converge, so the fit is not usable`); a broad
  `except Exception` only at a process boundary that passes the message on.
- Public package functions have type hints, with built-in generics and `|`:
  `def load_plates(data_dir: str | Path) -> pd.DataFrame:`; private helpers,
  tests and run scripts need none. A hinted type must be imported (3.12).

### Figures and tests

- One plotting module holds the palette and colormaps as constants at its
  top, a `set_style()` called by the run script and every notebook, the save
  function and the figure functions: `plot_<what>(data, ..., save=None)`
  returns the figure, or with `ax=None` draws one panel and returns the axes.
  `mapping/sepmap/plotting.py` has the constants (`RED`, `DARK_GREY`,
  `MID_GREY`, `DARK_BLUE`, `NO_DATA_GREY`, `GROUP_COLOURS` and the ISH
  additions above), the colormaps (`hot_cut`, `transparent_bad`),
  `set_style()` (white ground, fonts of 7 to 10 points, no top or right
  spine, Type 42 fonts in an EPS), `bars_grey`, `draw_plane` (one coronal
  plane in a panel: the atlas dark grey under the data, borders on top), the
  save function and `tidy` for the axes. The figure functions of the ISH
  analysis are in `mapping/sepmap/ish/plotting.py`.
- The save function closes the figure after `fig.savefig(save, dpi=150,
  bbox_inches="tight")`; new figures use its dpi, existing ones keep theirs.
  The route's `plotting.save_figure(fig, path, dpi)` writes the PNG and, by
  default, an EPS beside it with the image layers rasterised (PostScript has
  no transparency), and `<name>_new.png` when the PNG is open in a viewer. It
  leaves the figure open, because a video draws into the figure of its still
  again; the caller closes it.
- Titles and labels lowercase, names and symbols as written, the unit in
  parentheses (`"peak amplitude (uV)"`), counts computed
  (`f"recorded electrodes (n = {len(mapping)})"`). New figure files are
  numbered in run order (`01_raw_plates.png`).
- Tests use pytest, in `tests/` next to the package, run from the pipeline
  folder as `python -m pytest tests`; one `test_<topic>.py` per sub-package or
  topic, and a `tests/README.md` listing them (file, what it checks, data).
- A test is named for the expected behaviour
  (`test_flat_curve_is_classified_inactive`), has a one-line docstring, and
  a known answer on synthetic data: plant something, get it back, with the
  reason for a tolerance that is not obvious
  (`# planted channel recovered within noise (SEM about 1/sqrt(300))`).
- Include the case that must come out negative (a flat curve called
  inactive, a corrupted input refused); errors with `pytest.raises`. Tests
  keep `assert`, and write only into pytest's `tmp_path`.
- A test on real data is skipped when the data is not connected:
  `@pytest.mark.skipif(not Path(config.DATA).is_dir(), reason="data not connected")`.

### README, environments, linting

- Every pipeline folder's `README.md`: what it does in two to four lines,
  then Layout, Setup and run (which environment runs which script),
  Assumptions and preprocessing, Analyses, Outputs, and last Known
  limitations, stated plainly, with reasons and numbers.
- One pinned requirements file per environment in `tools/`, its top comment
  saying which Python and how to create it (`tools\venv_atlas` for the Python
  route, `tools\venv_flat` for the flatmaps); pytest and ruff in `tools\venv_dev`.
- `ruff.toml` at the root: line length 90, Python 3.12, rules E, W and F,
  double quotes. Run `tools\venv_dev\Scripts\python -m ruff format <files>`,
  then the same with `ruff check`.
