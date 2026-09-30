## Python

The Python code follows two reference projects by the same author, both read
only on `S:\ElboustaniLab\#SHARE\Documents\Giulio\`:

- the Genedata exercise, package `pka/`
  (`Genedata_exercise\PKA_exercise_Giulio_Matteucci.zip`);
- the MaxWell exercise, package `mea/`
  (`MaxWell_exercise\Retina_exercise_Giulio_Matteucci.zip`).

Citations: a module of a package as `pka/fitting.py:75`, any other file with
its project folder, as `genedata/run_pipeline.py:35` or
`maxwell/tests/test_mea.py:79`.

A rule with only reference citations is what the references do. Where only
one of them does it, the rule names that one. Where this guide departs from
them, the rule says "Departure" and names the source of the rule instead: the
four additions the owner asked for (no comment above `def` repeating the
docstring, explicit errors instead of `assert`, type hints on public
functions, ruff at 90 characters), the MATLAB half of the guide, the
style-pass section of the guide, or the refactor plan
(`docs/REFACTOR_PLAN.md`).

This draft becomes the Python half of `docs/STYLE.md`. The shared voice,
comment and data-safety rules, which open that file, apply unchanged.

Reference examples in this repository: one `run_*.py` and one module, named
at the start of the style pass from the first files restyled and approved.

### Layout of a pipeline

- A Python pipeline is one folder. The entry points (`run_*.py`),
  `settings.toml`, `README.md` and `tests/` sit at the top, and the code they
  call is one package below them. This matches the MATLAB layout (drivers on
  top, `pipeline/` below) and both references:

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

- The entry points sit next to the package, never inside it. `python run_x.py`
  puts only the script's own folder on the import path, so a script inside
  the package folder cannot import the package by name.
- One module per stage, named after the stage: `loading.py`, `preprocess.py`,
  `fitting.py`, `selection.py` (pka); `loading.py`, `sorting.py`, `units.py`,
  `responses.py`, `trajectory.py` (mea). Departure: a route with several
  analyses has one sub-package per analysis, each laid out the same way
  (the plan's target layout; neither reference has sub-packages).
- `config.py` loads the settings and holds the paths. `plotting.py` holds the
  palette, the style, the save function and every figure function; the
  modules that compute draw nothing. A sub-package with figures of its own
  has its own `plotting.py`.
- Only the run scripts have an `if __name__ == "__main__":` block; package
  modules have none (the only two in the references are
  `genedata/run_pipeline.py:102` and `maxwell/run_pipeline.py:162`).
- `__init__.py` is empty (`mea/__init__.py`) or holds a one-line docstring and
  nothing else. `pka/__init__.py:3-16` re-exports every public name, which
  pyflakes reports as 28 unused imports; not followed.
- Notebooks are optional, and numbered in reading order: `01_exploration.ipynb`
  (both references). A notebook only calls the package: "All logic is in the
  `mea` package; cells only call it." (`maxwell/01_exploration.ipynb`, first
  cell). Its first cell is markdown: a title, what the notebook shows, and that
  it uses the package; its first code cell imports, calls `set_style()` and
  loads the data (all six reference notebooks). No result is computed only in
  a notebook.

### Entry points (run_*.py)

- One `run_<step>.py` per step, named for what it does, listed in run order
  in the pipeline README. Departure: each reference has a single
  `run_pipeline.py`, because its whole analysis runs in one go. Ours has steps
  that take hours, fetch from the network, or need another environment
  (the plan, decision L1).
- The header is the module docstring (layout under Headers, below). As in the
  MATLAB drivers, it gives the pipeline and this script's step in it, what
  the script does, then the usage line and its options:

```python
"""Run the full analysis: load -> QC -> normalize -> fit -> select.

Writes the results table and all figures to the output folder.

    python run_pipeline.py [--data-dir data] [--out results] [--inhibitor CP_xxxxxxx]

--inhibitor picks the compound for the detailed curve figure; by default
the most potent candidate is shown.
"""
```
  (`genedata/run_pipeline.py:1-9`). With the step list of the MATLAB header,
  a sketch for this repository:

```python
"""Per-mouse volumes from the registered stacks.

Python route, volumes, step 1 of 3:
    1. run_per_mouse   tissue mask and backgrounds per brain  <- this script
    2. run_to_ccf      every brain on the adult CCF grid
    3. run_cohort      cohort mean, SD and count per voxel

What it does, and what it can be used for (2-3 lines).

    python run_per_mouse.py [mouse ...]
"""
```

- `main()` reads as the list of stages. Each stage is a blank line, a comment
  naming it, the calls, and one printed line saying what was done, with the
  counts that show it went right:

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
  and the counts for the printed lines (`genedata/run_pipeline.py:75-81`,
  `maxwell/run_pipeline.py:69-76`). Anything longer than a line or two of
  computing is a function in the package.
- Stage comments name the stage and are not numbered, as in pka
  (`genedata/run_pipeline.py:40, 45, 51`) and as the MATLAB `%%` titles. mea
  numbers them (`# stage 4: unit table and quality control`,
  `maxwell/run_pipeline.py:51`); not followed.
- Run options (which mice, whether to recompute) are `argparse` options. They
  are parsed under `if __name__ == "__main__":` and passed to `main()` as
  arguments, so a test or a notebook can call `main()` directly
  (`genedata/run_pipeline.py:35, 102-108`; mea parses them inside `main()`,
  `maxwell/run_pipeline.py:18-22`, not followed). Analysis parameters are
  never options: they live in `settings.toml`.
- Expensive steps are cached, with an option named for what it redoes, as
  MATLAB's `force_recompute_*`: `--resort` reruns the spike sorting, otherwise
  the cache is read (`maxwell/run_pipeline.py:21`, `mea/sorting.py:75-84`).
  The cache file name is a setting (`maxwell/settings.toml:26`).
- Progress is reported with `print`: one lowercase line per stage, sub-items
  indented by two spaces (`maxwell/run_pipeline.py:121`, as `fprintf` does in
  the MATLAB code). Neither reference uses the `logging` module. Departure:
  `flush=True` on prints inside long loops, because our long runs are
  detached and write to a log file.
- The figure backend is set only in run scripts, never in a package module,
  so notebooks still show their figures (pka sets it in the run script,
  `genedata/run_pipeline.py:14-16`; mea never sets it). Departure: it goes on
  the first lines under `if __name__ == "__main__":`, not between the imports
  as in pka, where it makes every later import an E402. `matplotlib.use`
  switches the backend even after pyplot has been imported (matplotlib 3.11,
  `matplotlib/__init__.py:1216-1280`), so every import stays at the top:

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
  (`genedata/settings.toml:1-7`)

```toml
baseline_guard_post_s = 0.5    # gap time discarded after each offset (late OFF responses), checked on the offset PSTH
```
  (`maxwell/settings.toml:43`). In this file end-of-line comments are the
  rule, unlike in code: they are the units column. A comment stays on its
  key's line even past 90 characters.

- Parameters that change what a figure shows (colour limits, clip quantiles,
  bin sizes, frame rates) are settings too
  (`genedata/settings.toml:33-34`, `maxwell/settings.toml:14-16, 55-58`).
  Figure sizes and colours stay in `plotting.py` (both references).
- A number goes to the settings when someone could reasonably choose another
  value. A number fixed by the method stays in the code, with a comment
  saying what it is: `0.6745` turns a median absolute deviation into a
  standard deviation (`mea/units.py:82`, which lacks the comment). Departure:
  the references keep a few choices in the code despite their headers
  (`p < 0.05` in `mea/trajectory.py:68`, the top and bottom four rows in
  `pka/preprocess.py:62-63`, the 12 gallery compounds in
  `genedata/run_pipeline.py:89`); here those are settings.
- `config.py` loads the file once, with `tomllib`, into `SETTINGS`:

```python
"""Load settings.toml once and expose it as a plain dict."""
import tomllib
from pathlib import Path

_SETTINGS_PATH = Path(__file__).resolve().parent.parent / "settings.toml"
```
  (`mea/config.py:1-5`; `pka/config.py:1-7` is the same in three lines)

- A module names the tables it uses at the top, after the imports, each
  named after its table in upper case (the last part of a dotted name):
  `QC = SETTINGS["qc"]`, `OUTLIERS = SETTINGS["outliers"]`,
  `FIT = SETTINGS["fit"]` (`pka/preprocess.py:21-22`, `pka/fitting.py:30`).
  Not abbreviations: `REC`, `SRT`, `UNI`, `DGN` in `mea/units.py:9-12`.
- In new code, a setting that a caller (a test, a notebook, a sensitivity
  check) may need to change is a keyword argument that defaults to `None`,
  and is read from the settings in an `if` block when it is not given:

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
  (`flag_pct=QC["gradient_flag_pct"]`, `pka/preprocess.py:50, 87-88`): that
  value is fixed when the module is imported, and a notebook that changes
  `SETTINGS` afterwards does not change it. The `if` block, not the one-line
  `x = ... if x is None else x` (mea uses the block 16 times and the one-liner
  8 times, as in `mea/responses.py:17-19`). The style pass never changes an
  existing signature (the style-pass section).
- Paths are `pathlib.Path` objects, joined with `/`. No `os.path` and no
  paths built by adding strings (neither reference uses `os.path`):

```python
    paths = sorted(Path(data_dir).glob("Plate_*.json"))
```
  (`pka/loading.py:43`); `out.mkdir(parents=True, exist_ok=True)`
  (`maxwell/run_pipeline.py:25`); `out / "units.csv"`.

- Departure: the references keep their data inside the project and take the
  folders as options with relative defaults
  (`genedata/run_pipeline.py:104-105`). Our data is on another drive, so
  `config.py` is the only place the data root is set. As `get_paths.m`
  does for MATLAB, it derives it from where the code sits (the parent of
  the code folder, plus `data`), and the same variable, `SEP_DATA_ROOT`,
  redirects the whole data tree for checks (the plan, Verification
  design). The cohort (each mouse's group, age and atlas) is also defined
  once (the plan, decision Y4). No other module writes a
  `D:\` path or lists mice.
- Code writes only into its own output folders under the data root from
  `config.py`. Nothing writes
  under the raw data on `S:\`.
- Fixed values that are not parameters (column lists, layouts, label maps)
  are module constants in upper case, after the imports and the settings
  names, with a comment above when the value needs one: `RESULT_COLUMNS`
  (`genedata/run_pipeline.py:27-31`), `ASSUMED_ANGLES_DEG`
  (`mea/trajectory.py:10-15`). A constant used only inside its module starts
  with an underscore: `_SETTINGS_PATH` (`mea/config.py:5`), `_MOSAIC_CELLS`
  (`mea/plotting.py:258-260`).

### Headers (module docstrings)

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
  (`pka/preprocess.py:1-13`; `pka/fitting.py:1-21` is the longest one)

- Layout: the summary on the first line, just after the opening quotes, then
  a blank line and the body, with the closing quotes on their own line. Then
  a blank line and the imports. This is the pka form. mea runs the summary
  into the body, closes on the last line of text and puts the imports on the
  next line (`mea/activity.py:1-3`); not followed.
- Formulas go on their own line, indented four spaces
  (`pka/preprocess.py:6`). Papers are cited by author and year: "(Hill 2011,
  as in SpikeInterface)" (`mea/units.py:64`). A workaround gives its link: the
  issue link in `mea/sorting.py:4-5`.
- Say why the method is the right one, and what the obvious alternative
  would get wrong, in a line or two: "the doses span more than four decades,
  so the log keeps the optimizer on a well-scaled axis"
  (`pka/fitting.py:3-5`). When a method replaced another, say in one or two
  lines what the other one gets wrong, not how the code got here. Which
  version did what, and when, belongs in the git log or in `docs/` (the
  history rule of the MATLAB half).
- No bold, no capitals for emphasis, no markdown other than code in
  backticks.
- The usage line is indented four spaces and starts with plain `python`:
  `python run_per_mouse.py [mouse ...]`
  (`genedata/run_pipeline.py:5`). Never the path of one environment's
  interpreter; the README says which environment runs the script.

### Docstrings

- Every function has a docstring, private and nested ones too, and so does
  every class. One line is enough for a small helper. Source:
  the style-pass section. pka documents every function; mea leaves
  nine without one (`open_recording`, `mea/loading.py:12`; `_finish`,
  `mea/plotting.py:55`; the nested `far`, `mea/sta.py:56`); neither documents
  its tests.
- Plain prose, with no `Args:`, `Returns:` or `Parameters` sections (neither
  reference uses the NumPy or the Google format). A summary line in sentence
  case that ends with a period and fits on one line; then, only where the
  names don't make it obvious: what the inputs are (shapes, units), what
  comes back, when a value is NaN, and why the method is the right one.
  Argument names go in backticks (pka; mea writes them bare, as in
  `mea/sta.py:44`). A list of rules is numbered in plain text
  (`pka/preprocess.py:91-100`).

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
- Departure (owner's addition, no comment repeating the docstring; decision
  Y1): no comment line above `def`. The references put a lowercase comment above 109
  of their 114 functions, followed by a docstring that says the same thing:

```python
# read one plate file into a table with one row per well
def load_plate(path):
    """Read one plate JSON file into a tidy DataFrame with one row per well."""
```
  (`pka/loading.py:15-17`). Here the docstring alone carries it. If the
  comment says something the docstring doesn't, that goes into the docstring.
  A test's comment becomes its docstring.

- No blank line between the docstring and the first line of code (none of
  the references' 114 functions has one). This differs from MATLAB, where a
  blank line follows the help block.

### Comments

The shared rules apply: short, lowercase, on their own line, the reason
rather than the action. In Python, from the references:

- A step, named by a noun phrase or a verb, whichever reads naturally:
  `# per-plate quality control` (`genedata/run_pipeline.py:45`). Don't rewrite
  a clear comment only to put a verb first.
- A reason, in one line:

```python
    except RuntimeError:
        # the optimizer failed to converge, so the fit is not usable
        return failed, np.zeros(len(y), dtype=bool)
```
  (`pka/fitting.py:75-77`)

```python
    # free the whitened array before re-reading raw, keeps peak memory low
    del traces, rec
```
  (`mea/sorting.py:89-90`)

- A decision, and what the alternative would get wrong:

```python
    # double-angle vector: strong for cells tuned to a bar axis, where the two
    # opposite response lobes cancel the ordinary direction vector
    axis_vector = np.sum(rates * np.exp(2j * angles))
```
  (`mea/responses.py:147-149`)

- A dense line gets a comment naming what it builds:
  `# design matrix: intercept, soma x, soma y` (`mea/trajectory.py:60`).
- A check gets a comment stating what is expected: `# the DIO record must
  alternate on/off starting with on, by stimulation design`
  (`mea/loading.py:75`).
- No end-of-line comments in `.py` files (shared rule). Departure: the
  references have 19, such as
  `matplotlib.use("Agg")  # figures go to files, no display needed`
  (`genedata/run_pipeline.py:16`) and
  `assert len(chosen) == n  # accepted population is large, ...`
  (`mea/sta.py:85`); here the comment goes on the line above. Tool pragmas
  such as `# noqa: E402` are not comments for the reader and stay where the
  tool needs them.
- Comments are wrapped at 90 characters like the code (owner's addition; the
  references run to 125, `pka/plotting.py:98`). A comment is one or two
  lines, as in MATLAB; a longer explanation of a method goes into the
  docstring of the function or module that implements it. A comment above a
  constant that justifies its value may be longer
  (`mea/trajectory.py:10-14`).

### Spacing and line layout

- Two blank lines between top-level functions and classes, one between
  methods.
- Blank lines split a function body into blocks, one per step. In `main()`
  every block starts with a comment (both references). Elsewhere a block gets
  a comment when what it does is not clear from the names: `unit_diagnostics`
  comments two of its five blocks (`mea/units.py:104, 118`). Inside a long
  loop, each step of several gets its own comment
  (the style-pass section).
- A function longer than about 60 lines, or with more than five or six
  steps, is split into functions (new code and the structural steps of the
  plan; the style pass does not split functions). The references have 114
  functions, median 14.5 lines; five exceed 60: `main` in
  `maxwell/run_pipeline.py` (143), `plot_unit_diagnostics` (80), `main` in
  `genedata/run_pipeline.py` (65), `unit_diagnostics` (64) and `fit_series`
  (62).
- Functions come in the order the pipeline uses them, each helper just above
  the first function that calls it: `hill`, `_hill_drop`, `fit_series`,
  `fit_all` (`pka/fitting.py:33-134`); `unpack_frame_numbers`,
  `frame_number_at`, `check_frame_contiguity` (`mea/loading.py:44-65`). This
  departs from the MATLAB template, which puts private helpers last (the
  MATLAB half follows decision Y6). The style
  pass does not reorder functions.
- No separator lines inside functions. Modules are split by stage, so the
  references need none between functions either. Departure, from the MATLAB
  half: a module that is not yet split may open each group of functions with
  one banner, `# ===== Registration =====`.
- One statement per line: no `;`.
- Lines are at most 90 characters (owner's addition; ruff E501; the
  references have 59 longer lines). Break inside brackets, either with the
  continuation lines aligned on the opening bracket:

```python
    ax1.errorbar(x, controls["neutral_mean"], yerr=controls["neutral_sd"], fmt="o-",
                 color=BLUE, capsize=2, label="neutral (100% activity)")
```
  (`pka/plotting.py:171-172`), or with a hanging indent of four and the
  closing bracket on its own line (`pka/plotting.py:432-435`,
  `mea/plotting.py:92-95`). Both references use both; keep the form a
  statement already has. No backslash continuations: use brackets (departure
  from `pka/plotting.py:150-151, 379-380` and
  `maxwell/run_pipeline.py:119-120`).

- Strings in double quotes, single quotes only inside an f-string's braces
  (existing code included, decision Y2):
  `f"{wells['plate'].nunique()} plates"` (both references). Formatting with
  f-strings only; neither reference uses `%` or `.format`.

### Naming

- `snake_case` for modules, functions and variables, `UPPER_CASE` for
  constants and settings tables, and a leading underscore for private
  helpers and constants: `_save`, `_hill_drop` (pka), `_finish`,
  `_footprint_grid`, `_MOSAIC_CELLS` (mea).
- A function name says what it returns or does. `load_*` reads files,
  `check_*` validates and fails, `plot_*` draws a figure, `*_table` returns a
  table, and there are `fit_*` and `select_*` functions (`pka/loading.py`,
  `mea/units.py`, `mea/responses.py`).
- A quantity with a unit carries the unit as a suffix, in variables, table
  columns and settings keys alike: `duration_s`, `fs_hz`, `min_amplitude_uv`,
  `ac50_max_um`, `refractory_pct`, `detect_channel_radius_um`. Counts start
  with `n_`: `n_bins`, `n_spikes`.
- Table columns are the names readers see in the CSV: `evoked_hz`,
  `latency_s`, `baseline_sd_hz` (`mea/responses.py:89-94`),
  `top_to_bottom_pct` (`pka/preprocess.py:69`).
- Short names only for a loop variable used over a few lines
  (`for u, tr in trains.items()`, `mea/units.py:23`) or for the standard
  symbol of a formula next to its comment (`A` for the design matrix,
  `mea/trajectory.py:60-61`). Anything used for longer has a readable name.
  Departure, from the owner's rule against cryptic abbreviations: mea keeps
  the open file as `f` through its 143-line `main()`
  (`maxwell/run_pipeline.py:30`).

### Functions and data

- A function takes data and returns data. It reads no globals other than
  settings and constants, and writes files only when that is its job
  (`sort_and_cache`, `mea/sorting.py:75-99`; the run scripts; the `save=` of a
  figure function).
- Tables are pandas DataFrames with one row per item (well, unit, structure)
  and named columns, built from a list of dicts:

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
  (`pka/preprocess.py:32-46`; `mea/units.py:22-39`). Columns are read by
  name, never by position.

- A function that adds columns works on a copy and returns it:
  `out = wells.copy()` (`pka/preprocess.py:81, 102`, `pka/selection.py:23`,
  `mea/units.py:55`). It changes its input in place only when memory
  requires it, and the docstring says so: "Bandpass then ZCA whiten, in
  place, chunked float32." (`mea/sorting.py:32`).
- Right after loading, a `check_*` function compares the data with what it
  must be (counts, layout, alignment) and fails if it differs; the run script
  prints that the checks passed (`genedata/run_pipeline.py:41-43` with
  `pka/loading.py:50-67`; `maxwell/run_pipeline.py:30-37` with
  `mea/loading.py:59-65`).
- An excluded item stays in the table, with a boolean column and a reason
  column: `excluded` and `exclude_reason` (`pka/preprocess.py:123-124`),
  `accepted` and `reject_reason` (`mea/units.py:55-57`); "each with its
  reason in the output" (`genedata/README.md:61`).
- Values that cannot be determined are NaN, never zero, and a comment gives
  the reason:

```python
            # a cell silent through every baseline chunk has no noise scale,
            # its effect size is reported as not determined
            dprime = (evoked - base_mean) / base_sd if base_sd > 0 else np.nan
```
  (`mea/responses.py:85-86`, comment wrapped at 90;
  `mea/trajectory.py:71-73`)

- Functions, not classes: neither reference defines a class. A class only
  for state that lives across calls, such as a reader holding an open file
  (`Source`, `v2_per_mouse.py:99`) or a network (`auto_annotation/core.py:188`).
- Keep expressions simple, as in MATLAB: one step per line, intermediate
  variables with readable names rather than nested calls. No comprehension
  nested more than two levels, no lambda assigned to a name (use `def`), no
  conditional expression spread over several lines (departure from
  `maxwell/run_pipeline.py:119-120` and `pka/plotting.py:379-380`).

### Random numbers

- Every draw comes from a local generator,
  `rng = np.random.default_rng(seed)`; never `np.random.seed` or the global
  `np.random` functions (both references).
- An analysis function that draws takes `seed=0` as a keyword argument and
  makes its own generator (`mea/sta.py:30-34`, `mea/trajectory.py:45-54`).
- A figure that jitters points uses a generator of its own, separate from
  the analysis (`pka/plotting.py:154, 201`).
- A test makes its own generator with a fixed seed (under Tests).

### Imports

- At the top, in three groups separated by a blank line: standard library,
  third party, this package (`maxwell/run_pipeline.py:8-14`,
  `pka/loading.py:8-12`). Within a group, in the order isort gives: plain
  `import` lines first, then `from` lines, each alphabetical
  (`pka/plotting.py:7-16`). `mea/sorting.py:13-24` gives the sorter libraries
  a group of their own; not followed.
- Run scripts, tests and notebooks import the package's modules and call
  `module.function`, so every call says where the function lives:

```python
from mea import activity, loading, plotting, responses, sorting, sta, trajectory, units
...
    epochs = loading.load_epochs(f)
```
  (`maxwell/run_pipeline.py:13, 33`; `maxwell/tests/test_mea.py:6`; the four
  MaxWell notebooks). Genedata imports names from the package's
  `__init__.py` instead (`genedata/run_pipeline.py:18-25`); not followed.

- Inside the package, a module imports the names it uses from a sibling:
  `from pka.fitting import hill` (`pka/selection.py:11`,
  `mea/activity.py:7`).
- Absolute imports within the package: `from pka.config import SETTINGS`
  (pka). mea uses relative ones (`from .config import SETTINGS`,
  `mea/loading.py:6`); not followed.
- Departure: an import inside a function only for a dependency that is not
  installed in every environment that runs the module, with a comment saying
  so. The references import inside functions without such a reason
  (`pka/plotting.py:356`, `mea/plotting.py:391, 485, 529`).
- The package never edits `sys.path` (neither reference does). The one
  exception is a script that MATLAB runs by path
  (`auto_annotation/cli.py:40-42`), with a comment saying why.

### Errors

- Check what can go wrong with the data, and fail with a message that says
  what was found, what was expected and, when there is one, the fix (as the
  MATLAB `error()` messages do).
- Departure (owner's addition): `raise`, not `assert`. The references check
  their data with bare asserts (`pka/loading.py:53-67`,
  `mea/loading.py:40, 64, 76-77`, `mea/sta.py:85`, `mea/trajectory.py:93`).
  Those vanish under `python -O` and say nothing when they fail. The form to
  use is the one pka uses for a missing input:

```python
    paths = sorted(Path(data_dir).glob("Plate_*.json"))
    if not paths:
        raise FileNotFoundError(f"No Plate_*.json files found in {data_dir}")
```
  (`pka/loading.py:43-45`). A design check, written this way, the comment
  above the old assert turned into the message:

```python
    n_plates = wells["plate"].nunique()
    if n_plates != 15:
        raise ValueError(f"expected 15 plates in 384-well format, found {n_plates}")
```
  (adapted from `pka/loading.py:52-53`)

- Built-in exception types: `FileNotFoundError` for a missing input,
  `ValueError` for data that is not what was expected, `RuntimeError` for a
  step that failed. A custom exception only when a caller catches that case
  by name (`NotReferenceGrid`, raised at `v2_ish_regions.py:102` and caught at
  184). Package code never calls `sys.exit` or raises `SystemExit`; only a
  run script turns a usage mistake into an exit (`parser.error`).
- Catch only what you expect, with the reason in a comment
  (`pka/fitting.py:75-77`, above). A broad `except Exception` only at a
  process boundary that passes the message on, with a comment:
  `auto_annotation/cli.py:149-150` (the GUI shows it).
- Tests keep `assert`, which is how pytest works.

### Type hints

- Departure (owner's addition): public functions of the package have type
  hints on their arguments and return value. Private helpers, nested
  functions, tests and run scripts need none. Neither reference uses hints;
  the docstring still gives the shapes and units, which hints cannot.
- Built-in generics and `|`: `list[str]`, `dict[str, float]`, `Path | None`;
  `np.ndarray` for arrays, `pd.DataFrame` for tables:

```python
def load_plates(data_dir: str | Path) -> pd.DataFrame:
    """Load all Plate_*.json files in `data_dir` into one tidy DataFrame."""
```
  (`pka/loading.py:41-42`, with hints added)

- On Python 3.12 (our environments) hints are evaluated when the `def` runs,
  so a hint that names a type that is not imported breaks the import of the
  module. Adding hints is a code change (kind B, below).

### Figures

- One plotting module holds the palette, the colormaps, `set_style()`, the
  save function and the figure functions. Figure functions are
  `plot_<what>(data, ..., save=None)` and return the figure. A function that
  draws one panel takes `ax=None` and returns the axes
  (`pka/plotting.py:370-408, 412-441`).

```python
"""Shared figure style and one function per figure. Figure functions take
save=None; the pipeline passes a path, notebooks do not."""
```
  (`mea/plotting.py:1-2`, a docstring in the mea layout)

- Figures are saved through one function, and a saved figure is closed so
  that a long run does not accumulate them:

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
  (`mea/plotting.py:55-62`). pka's `_save` does not close the figure
  (`pka/plotting.py:35-39`), and `genedata/run_pipeline.py:80` saves around
  it; not followed. Ours, `save_figure`, also writes an EPS next to the PNG
  for the figures that are edited for a paper, with the image layers
  rasterised (`v2_region_plot.py:104-130`).

- New figures use the dpi of the save function. Existing figures keep their
  own, since changing it changes every output.
- The palette and colormaps are named constants at the top of the plotting
  module, with a comment above giving the role when it is not obvious
  (`mea/plotting.py:10-28`). pka gives the role at the end of the line
  (`PURPLE = "#9333EA"  # marks excluded wells in all figures`,
  `pka/plotting.py:24`); not followed. The colours themselves: see Figures
  (both languages).
- `set_style()` sets the matplotlib rcParams once (`pka/plotting.py:43-63`,
  `mea/plotting.py:32-51`). The run script calls it
  (`genedata/run_pipeline.py:38`, `maxwell/run_pipeline.py:26`), and so does
  the first code cell of every notebook.
- Titles and axis labels in lowercase, names and symbols written as they
  are, the unit in parentheses: `"time from onset (s)"`,
  `"peak amplitude (uV)"` (`mea/plotting.py:301, 156`), as our figures already
  do. pka's sentence-case titles (`"Fit quality"`, `pka/plotting.py:222`) are
  not followed. A count shown in a title is computed, never typed:
  `f"recorded electrodes (n = {len(mapping)})"` (`mea/plotting.py:76`).
- New figure files are numbered in run order: `01_raw_plates.png`
  (`genedata/run_pipeline.py:64`), `fig01_raw_overview.png`
  (`maxwell/run_pipeline.py:42`). Existing file names stay as they are.

### Tests

- pytest, in `tests/` next to the package. Run them from the pipeline folder:

```
python -m pytest tests
```
  (`maxwell/README.md:45`). Plain `pytest` puts `tests/`, not the pipeline
  folder, on the import path, and cannot import the package:
  `genedata/README.md:43` says `pytest`, and on a copy of that project
  `pytest` stops with a collection error where `python -m pytest` passes all
  six tests.

- Departure: one file per sub-package or topic, `test_<topic>.py`, as the
  MATLAB tests are organised. Each reference has a single file named after
  its package (`genedata/tests/test_pka.py`, `maxwell/tests/test_mea.py`).
- A test is a function named for the expected behaviour:
  `test_planted_outlier_is_flagged_and_parameters_recovered`,
  `test_flat_curve_is_classified_inactive` (`genedata/tests/test_pka.py:43,
  55`). A one-line docstring says what it checks (the references use a
  comment above instead; see Docstrings).
- A test has a known answer, on synthetic data: plant something, get it back.
  When a tolerance is not obvious, a comment gives its reason:

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

- Each test makes its own generator with a fixed seed:
  `rng = np.random.default_rng(4)` (`genedata/tests/test_pka.py:44`).
- Include the case that must come out negative: a flat curve called inactive
  (`genedata/tests/test_pka.py:54-59`), an untouched channel that stays flat
  (`maxwell/tests/test_mea.py:55`), a corrupted input refused
  (`maxwell/tests/test_mea.py:35-41`). A test that could only pass checks
  nothing.
- Errors are tested with `pytest.raises(ValueError)`.
  `maxwell/tests/test_mea.py:40` expects an `AssertionError`, which the
  "raise, not assert" rule replaces.
- Departure (no reference test needs data): a test that needs real data is
  skipped, not failed, when the data isn't connected
  (the plan, step 10, as the MATLAB tests do):
  `@pytest.mark.skipif(not config.DATA_ROOT.exists(), reason="data drive not connected")`.
- Tests write only into pytest's `tmp_path`, never into the output folders.
- A `tests/README.md` lists the tests in a table (file, what it checks,
  data), as the MATLAB tests README does.

### README of a pipeline

Every pipeline folder has a README, with the sections of the references:

```
# <Pipeline name>

What it does in two to four lines: the question, the data, what comes out.

## Layout
(the tree, one line per file)

## Setup and run
(which environment, the commands in run order, how long, what is cached)

## Assumptions and preprocessing
(one bullet per choice: what was done, why, and the numbers behind it)

## Analyses
(one bullet or paragraph per analysis, same content)

## Outputs
(every file written, one line each)

## Known limitations
```
(`maxwell/README.md:1-121`; Genedata has the same sections, without the
limitations)

- The bullets give the reasons and the numbers, as in "The three replicates
  of each inhibitor are on three different plates, so each well is
  normalized to its own plate's controls" (`genedata/README.md:51-55`) and
  "d' ... 15-47 on all plates" (`genedata/README.md:56-58`).
- Limitations are stated plainly: "Three template-similar unit pairs are
  most likely single cells split in two by the sorter (evidence in notebook
  02). They are documented rather than merged, so the accepted-unit count is
  slightly inflated; no conclusion depends on it."
  (`maxwell/README.md:115-118`).
- Departure: the README ends with its known limitations; the closing
  section of the references is not used.

### Environments

- One requirements file per Python environment, pinned to the versions the
  results were produced with, with a comment at the top saying which Python
  and how to create the environment, and the packages in commented groups
  (`genedata/environment.yml:4-14`; `maxwell/environment.yml:5-19` pins minor
  versions only). Departure: the references use conda `environment.yml`; ours
  are `venv`s made from the Anaconda Python 3.12.7
  (`tools/venv_atlas/pyvenv.cfg`), so the same content goes in a
  requirements file, kept in `tools/` with the environments
  (the plan, step 0):

```
# tools\venv_atlas: Python 3.12 (Anaconda). Create it from the repository root with
#   python -m venv tools\venv_atlas
#   tools\venv_atlas\Scripts\python -m pip install -r tools\requirements_atlas.txt

# scientific stack (the versions the results were produced with)
numpy==2.5.2
pandas==3.0.5
scipy==1.18.1
matplotlib==3.11.1

```

- The README of each pipeline says which environment runs which script:
  `tools\venv_atlas` for the Python route, `tools\venv_flat` for the
  flatmaps, the auto-annotation's own environment for its engine.
- pytest and ruff live in `tools\venv_dev`, built from
  `requirements_atlas.txt` plus a `requirements_dev.txt` that adds them, so
  the analysis environments never change (the plan, Environments). Tests
  run with `tools\venv_dev\Scripts\python -m pytest tests` from the
  pipeline folder.

### Linting

- Departure (owner's addition): `ruff check`, configured once at the top of
  the repository:

```toml
# ruff.toml
line-length = 90
target-version = "py312"
extend-exclude = ["third_party", "archive"]

[lint]
select = ["E", "W", "F"]
```
  This catches long lines (E501), imports below other code (E402), several
  statements on one line (E701-E703), a lambda assigned to a name (E731),
  unused imports and undefined names (F). It runs from `tools\venv_dev`.

- No automatic formatter (no `ruff format`, no black). The references are
  formatted by hand: black would change 361 of their 2,918 lines at line
  length 90, and still 345 at line length 200, so its changes are not about
  line length. On our code it would rewrite most lines at once.

### Applying this to existing code

The style pass has two kinds of change (see the style-pass section of the guide).
In Python, `tools/check_code_identity.py` compares the `ast` of each file
with comments and docstrings removed. Whatever it reports is the answer.

Kind A, which should pass the identity check:

- comments: add, rewrite, delete; move an end-of-line comment to its own
  line above;
- docstrings: add, rewrite, rewrap; fold a comment that sits above `def`
  into the docstring and delete the comment (a test's comment becomes its
  docstring); rewrite a module docstring (history out, method and reason in,
  no bold or capitals, usage line with plain `python`);
- blank lines, indentation inside brackets, spaces;
- `a = 1; b = 2` onto two lines;
- long lines broken inside brackets; a backslash continuation replaced by
  brackets; a long plain string split into adjacent literals (for an
  f-string, run the check);
- `1.` written `1.0`, `.14` written `0.14`;
- single quotes changed to double quotes, a whole file at a time (decision
  Y2);
- exception: in a module that uses its docstring at run time, the docstring
  is program output, so editing it is kind B. `auto_annotation/cli.py:157`
  prints the whole docstring as its usage. `build_demba_atlas.py:202` prints
  the second-to-last line of its docstring. `landmark_refine/cli.py:50` and
  `landmark_refine/serve.py:109` print theirs. Search for `__doc__` before
  editing a docstring.

Kind B, which the check reports as `CODE CHANGED` and which is verified by
rerunning:

- rename a local variable (every use; never a name saved to a file, a
  keyword that callers pass, a dict key or a column name);
- split a dense expression into named intermediates, with the same
  arithmetic in the same order and the same dtypes (the `float16` and
  `float32` casts stay where they are);
- a lambda assigned to a name becomes a `def`; a conditional expression
  becomes `if`/`else`, including `x = default if x is None else x`;
- remove an unused import (first check that it has no side effect);
- an f-string without placeholders becomes a plain string;
- add type hints to public functions (then import every changed module);
- `assert` becomes `raise`, with the comment above the assert as the
  message; `SystemExit` in package code becomes a built-in exception. First
  check what catches the old type: `auto_annotation/cli.py` catches
  everything and hands the message to MATLAB;
- change printed text only if nothing reads it.

Not part of the style pass (structural steps of the plan, each checked with
a rerun):

- moving files into the package, and the new import lines;
- `os.path` to `pathlib`. A `Path` has no `str` methods, so
  `path.replace('.png', '_new.png')` (`v2_region_plot.py:114`) and
  `mouse + '.npz'` fail on one, and `Path.replace` renames a file;
- hard-coded roots into `config.py`, module constants into `settings.toml`,
  cohort lists into one place;
- `matplotlib.use("Agg")` out of the modules and into the run scripts;
- changing a signature: settings as keyword arguments defaulting to `None`,
  a `seed` argument;
- reordering functions, splitting a long `main()`, merging helpers defined
  twice (`per_unit` in `v2_cohort.py:188` and `v2_region_plot.py:203`);
- tuples read by position (`r[11]`, `v2_region_plot.py:330`) into named
  columns;
- figure dpi and output file names;
- environment variables that choose what to compute (`V2_READINGS`,
  `v2_cohort.py:96`; `V2_ISH_PANEL`, `V2_ISH_TABLE`,
  `v2_ish_regions.py:71-72`) into settings or options.

Random state: never add, remove or reorder a call that draws from a
generator. In `v2_ish_panel_test.py`, one `rng` (line 196) feeds up to five
permutation tests in turn (lines 203, 215, 225, 240), so reordering them
changes every p-value. The figure has its own generator (line 257).

For each file, check:

1. module docstring: summary on the first line, method and reason, no
   history, no bold or capitals, usage line with plain `python`;
2. every function and class has a docstring, with no comment line above
   `def` and no blank line after the docstring;
3. comments lowercase, on their own line, at most two lines, a blank line
   before each block, one per stage in `main()`;
4. no `;`, no line over 90 characters, no backslash continuation, double
   quotes;
5. `ruff check` (or `pyflakes`) shows no new warning compared with the
   baseline file;
6. `check_code_identity.py`: `same code` for every file with kind A changes
   only.
