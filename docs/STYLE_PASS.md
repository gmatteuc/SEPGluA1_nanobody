# The style pass

Temporary: the procedure of step 7 of [REFACTOR_PLAN.md](REFACTOR_PLAN.md),
which brings existing code in line with [STYLE.md](STYLE.md). It is deleted
when step 7 is done; the rules themselves are in STYLE.md, and only there.
The pass works file by file, from STYLE.md and its reference examples.

## Reference examples: code not to copy

The four reference examples have their comments, headers, docstrings and
layout in final style. Some of their code is not, because a comment and
layout pass does not change code; do not copy these:

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

## Old code

- Existing commented-out code stays until the owner decides, with a comment
  saying why it is off, or that the reason is not recorded.
- Don't rewrite a clear comment only to put a verb first.
- The pass renames nothing unless the assignment allows it (kind B); a rename
  STYLE.md would want goes on a list for the owner.
- MATLAB text in double quotes stays: `"..."` is a string, not a character
  vector, so changing the quotes changes the code. Y2 (double quotes) is for
  Python only.
- A MATLAB string literal that does not fit on a line stays whole: splitting
  it into `['...' '...']` changes the code (kind B). A list of names after a
  command (`clear a b c`) can continue with `...`; the code stays the same.

## The Python reference projects

They are `Genedata_exercise\PKA_exercise_Giulio_Matteucci.zip` and
`MaxWell_exercise\Retina_exercise_Giulio_Matteucci.zip` in
`S:\ElboustaniLab\#SHARE\Documents\Giulio\`; the first version of STYLE.md,
in the git history, cites each of its quotes from them by file and line. Their
habits that STYLE.md does not follow: a lowercase comment above 109 of their
114 functions repeating the docstring (Y1); single quotes (Y2); data checked
with bare `assert` (`pka/loading.py:53-67`, `mea/loading.py:40, 64`); no type
hints; 19 end-of-line comments and comments up to 125 characters; layout by
hand (black would change 361 of their 2,918 lines), with backslash
continuations (`pka/plotting.py:150-151`, `maxwell/run_pipeline.py:119-120`);
no linter; one `run_pipeline.py` each and a single test file each;
`pka/__init__.py:3-16` re-exporting every public name; numbered stage comments
and options parsed inside `main()` (mea); the backend set between the imports
(pka); a docstring summary run into its body (`mea/activity.py:1-3`); nine
functions without a docstring (mea, such as `_finish`); a save function that
does not close the figure (pka's `_save`); a colour's role at the end of its
line (pka); names imported from the package's `__init__.py` (Genedata) or
relative imports (mea); imports inside functions without a reason
(`pka/plotting.py:356`); data inside the project, its folders taken as options
(`genedata/run_pipeline.py:104-105`); choices kept in the code (`p < 0.05` in
`mea/trajectory.py:68`, the 12 gallery compounds in
`genedata/run_pipeline.py:89`); a setting as a default value
(`pka/preprocess.py:50`); table names abbreviated (`REC`, `SRT`, `UNI`,
`mea/units.py:9-12`); `0.6745` without its comment (`mea/units.py:82`); the
open file kept as `f` through a 143-line `main()`
(`maxwell/run_pipeline.py:30`); `pytest` run plainly
(`genedata/README.md:43`), which on a copy of that project stops with a
collection error.

## Three kinds of change

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
- add type hints to public functions (then import every changed module:
  Python 3.12 evaluates hints when the `def` runs, so a hint naming a type
  that is not imported breaks the import of the module);
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
pass: it goes on the bug list of the plan.

**Kind C: structural changes.** The plan gives the style pass
these too (`matplotlib.use` and the module constants are left "for the
style pass" in its Progress). They come after the kind A and B commits of
the files they touch, one change per commit, each checked by rerunning
every step that uses the changed code and comparing as for kind B:

- `matplotlib.use` out of the modules and into the run scripts (it switches
  the backend even after pyplot has been imported: matplotlib 3.11,
  `matplotlib/__init__.py:1216-1280`);
- constants into `settings.toml` (`SMOOTH`, the 100-voxel minimum), and a
  module's settings read through its table
  (`YOUNG_VS_ADULT = SETTINGS["young_vs_adult"]`) rather than key by key;
- `os.path` to `pathlib` (a `Path` has no string methods:
  `path.replace(".png", "_new.png")` in `save_figure` would rename a file);
- changing a signature (settings as keyword arguments, a `seed` argument);
- merging helpers defined more than once: `per_unit`
  (`mapping/sepmap/volumes/cohort.py:190` and
  `young_vs_adult/region_plot.py:195`), `fold_n` (`compare.py`, `closeup.py`,
  `video.py`), and the five copies of `save_figure` (`diagnostics.py:50`; in
  `young_vs_adult/`, `closeup.py:117`, `region_groups.py:91`,
  `region_plot.py:104`, `compare.py:66`), none of which closes the figure.
  All but the first write an EPS next to the PNG, with the image layers
  rasterised, and all write `<name>_new.png` when the PNG is open in a
  viewer; the one save function of STYLE.md keeps those two behaviours;
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

## Checks for each file

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
5. `ruff check` shows no warning the file did not have before (it catches
   long lines, E501; imports below other code, E402, though ruff 0.16 does
   not count a `matplotlib.use` call as code; several statements on one
   line, E701-E703; a lambda assigned to a name, E731; unused imports and
   undefined names, F);
6. `check_code_identity.py`: `same code`.

The kind A commits are listed in `.git-blame-ignore-revs`, so `git blame`
skips them; a kind A commit therefore holds nothing else, not even a new
file. Kind B and C commits change code, so blame keeps them. At the end, a
search that no file credits an author other than Giulio.

## README of a pipeline

Every pipeline folder gets a README with the sections of the references
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

The Layout section as in Genedata:

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

- The bullets give reasons and numbers: "The three replicates of each
  inhibitor are on three different plates, so each well is normalized to its
  own plate's controls" (`genedata/README.md:51-55`).
- Limitations are stated plainly: "Three template-similar unit pairs are
  most likely single cells split in two by the sorter (evidence in notebook
  02). They are documented rather than merged, so the accepted-unit count is
  slightly inflated; no conclusion depends on it."
  (`maxwell/README.md:115-118`).
- Setup and run says which environment runs which script: `tools\venv_atlas`
  for the Python route, `tools\venv_flat` for the flatmaps, the automatic
  annotation's own environment for its engine.
- The README ends with its known limitations.

## Environments

One requirements file per environment, in `tools/`, pinned to the versions
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

pytest and ruff live in `tools\venv_dev`, so the analysis environments never
change:

```
tools\venv_dev\Scripts\python -m ruff format <files>
tools\venv_dev\Scripts\python -m ruff check <files>
```
