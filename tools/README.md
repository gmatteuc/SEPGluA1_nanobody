# Tools

Checks that a change to the code does not change the results, and the runner
for long MATLAB stages. The checks follow the verification design of
[docs/REFACTOR_PLAN.md](../docs/REFACTOR_PLAN.md): old and new code each run
in a fresh session on their own copy of the data, never on the production
data, and their outputs are compared file by file.

## Where a run reads and writes: SEP_DATA_ROOT

`get_paths.m` (MATLAB) and `mapping/sepmap/config.py` (Python) take the data root from the
environment variable `SEP_DATA_ROOT` when it is set, otherwise from the
folder next to the code (`<root>\data`). It moves the whole data tree, inputs
and outputs. Both refuse the production data (`D:\sep_histology\data`) to a
copy of the code, and any data root inside the snapshot on G:.

A check points `SEP_DATA_ROOT` at its own copy, through the runner below or
`sep_run_driver_copy`. `sep_run_driver_copy` refuses the production data;
the runner also serves production runs, so it takes the production data as
a data root but refuses a check's log there. Both refuse the snapshot, before
anything runs.

## Long runs: run_matlab_detached.ps1

```
powershell -File tools\run_matlab_detached.ps1 -Script run_collect_by_group -CodeDir D:\sep_histology\code -DataRoot D:\sep_histology\data -LogDir D:\sep_histology\data\young
```

`-CodeDir`, `-DataRoot` and `-LogDir` have no defaults. The runner passes the
data root to MATLAB as `SEP_DATA_ROOT` and appends everything to
`<LogDir>\_<script>.log`. Each stage runs in a fresh session:
`restoredefaultpath; cd(CodeDir); sep_setup_paths`, then it prints the code
root, which `get_paths` it resolved, the data root and the commit (also
written at the top of the log, with the number of uncommitted changes).

- Several stages, separated by commas, run in order; a failing stage stops
  the chain. `-WaitForPid` holds the start until another process has ended
  (at most `-MaxWaitMinutes`, 480 by default).
- A stage can be a call with arguments; commas inside brackets or quotes do
  not split it. Use single quotes only: double quotes do not reach MATLAB.
- Refused (exit 1, message on the console): a log folder in the production
  data when the data root is a copy, a data root or log folder inside the
  snapshot on G:, a data root that does not exist.

## After a comment or layout-only edit

```matlab
addpath('D:\sep_histology\code\tools')
[T, ok] = sep_check_code_identity('D:\sep_histology\code', 'C:\...\copy_before_the_edit')
```

```
python tools\check_code_identity.py NEW_DIR REF_DIR [--map name_map.csv]
```

Every `.m` (or `.py`) file is parsed and compared without comments and layout
(and, for Python, docstrings). A changed file fails, and so do a file that
does not parse, a new file with no counterpart and a reference file nothing
was compared with. Two folders with no file of the language stop the check,
since an empty comparison would pass. Files that were moved or renamed are
paired through the name map, one table for both languages (the names below
only show the format):

```
old_path,new_path
P5_collect_data_by_group.m,group_comparison/run_collect_by_group.m
,tests/sep_test_path.m
scratch.m,
```

An empty old path is a file added on purpose, an empty new path one removed
on purpose; a path in the map that does not exist stops the check.

## After a code change

1. Keep the old code in its own folder (a worktree at the old commit), next to
   its own copy of the inputs: the check trees of the plan.
2. Run each version in a fresh session on its own copy, with the runner. To
   change a driver's settings without editing it, run a copy of it:

   ```
   powershell -File tools\run_matlab_detached.ps1 -CodeDir G:\sep_refactor\check\code -DataRoot G:\sep_refactor\check\data -LogDir G:\sep_refactor\check\logs -Script "addpath('tools'); sep_run_driver_copy('run_group_differences', 'G:\sep_refactor\check\data', struct('exp_type', '''rws''', 'generate_diff_videos', 'false'))"
   ```

   Each field of the struct replaces the assignment to that setting (the
   value as MATLAB text); other lines change through regular expressions,
   `{pattern, replacement; ...}` as a fourth argument. A setting or pattern
   that matches nothing stops the run. `sep_run_driver_copy` writes the copy
   to `<data_root>\..\driver_copies\`, sets `SEP_DATA_ROOT` for the run and
   restores it after, even on an error, and runs only if `get_paths` then
   resolves the requested data root.
3. Compare the outputs, MATLAB's with MATLAB and the Python route's with
   Python:

   ```matlab
   T = sep_compare_outputs('G:\sep_refactor\ref\data\comparisons', ...
       'G:\sep_refactor\check\data\comparisons', ...
       struct('replace_text', {{'G:\sep_refactor\ref', 'G:\sep_refactor\check'}}, ...
              'newer_than', datetime(2026, 10, 1, 9, 0, 0)));
   ```

   ```
   python tools\compare_outputs.py OLD_DIR NEW_DIR --replace G:\sep_refactor\ref G:\sep_refactor\check --newer-than "2026-10-01 09:00"
   ```

   Each file is compared by type (tables, `.mat` and `.fig` contents, arrays,
   pixels, figures up to rendering, EPS without dates, video frames).
   `replace_text` / `--replace` makes the paths a run saves, which name its own
   tree, comparable. `newer_than` / `--newer-than` reports a file the run did
   not rewrite, so a run that wrote nothing cannot pass: compare output
   folders only with it. Only differences you expected, and can explain,
   should be listed.

A check only counts if it could have failed: the file dates show the run
rewrote its outputs, and for a fix, the old code fails the same test.

## Traps

- **The current folder comes first.** A function in the current folder
  shadows the path, so a session started in one code folder uses that
  folder's `get_paths` whatever else is on the path. The runner always starts
  in `-CodeDir`.
- **clear all.** Drivers start with it. `sep_run_driver_copy` runs the copy in
  a workspace of its own, so the variable it restores afterwards survives.
- **Renderer noise.** Anti-aliasing can put an edge or a glyph a fraction of
  a pixel elsewhere, which changes the pixels along it. `sep_compare_outputs`
  calls such an image `same render` only if every changed pixel, in each
  colour channel, is a mix of the original's colours within 1 pixel (to 8
  levels of 255); anything else is `DIFFERENT`: a colour changed at the same
  brightness, a changed digit, an image of another size. Every `same render`
  is listed with its count of changed pixels. A change confined to the
  anti-aliased rim of ink passes, and JPEG compression error is no mix, so
  compare figures saved as PNG. `compare_outputs.py` tolerates 1 level on
  under 0.1% of the pixels.
- **Dates inside files.** `.mat`, `.fig`, `.eps` and videos hold the time they
  were written, so their bytes always differ; the tools compare their
  contents. `compare_outputs.py` reads none of MATLAB's formats: `.mat`,
  `.fig` and videos go to `sep_compare_outputs`.
- **Random numbers.** Registration is not deterministic (elastix's sampler,
  the point subsampling in extraction and align): compare it on the plan's
  named measures, never file by file.
- **Caches** are reused when present. Delete the computed caches in a check
  tree before each run (the plan lists them).

## Files

| file | what |
|---|---|
| `run_matlab_detached.ps1` | runs MATLAB stages in fresh sessions, with the data root and a log |
| `sep_run_driver_copy.m` | runs a copy of a driver on a copy of the data tree |
| `sep_make_script_copy.m` | copy of a driver script with some settings changed |
| `sep_compare_outputs.m` | file-by-file comparison of two output folders |
| `sep_compare_figure_images.m` | are two figures the same up to rendering |
| `sep_check_code_identity.m` | code identity of two folders of `.m` files |
| `sep_summarise_value.m`, `sep_hash_array.m` | helpers to compare large results |
| `sep_struct_defaults.m` | fills an options struct with defaults |
| `check_code_identity.py` | code identity of two folders of `.py` files |
| `compare_outputs.py` | file-by-file comparison of the Python route's outputs |

The MATLAB tools come from the imaging repository's `tools/`, renamed with
`sep_` so that its copies, if they are on the path, cannot shadow them. They
call only each other. `sep_setup_paths` does not add `tools\`: a check adds
it with `addpath('tools')` from the code root. The Python tools need numpy,
pandas and pillow (`tools\venv_atlas` has them).
