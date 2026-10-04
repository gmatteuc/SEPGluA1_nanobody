# Tests

Checks of the shared definitions that every pipeline relies on, where a
mistake would be silent: the MATLAB path and the cohort table. Each runs in a
fresh MATLAB session, from the code root:

```matlab
restoredefaultpath; cd('D:\sep_histology\code'); sep_setup_paths
run(fullfile('tests', 'sep_test_path.m'))
run(fullfile('tests', 'test_backward_compat.m'))
```

`sep_setup_paths` does not put `tests\` on the path, so each test is run by
its file. Both print what they check and stop with an error when a check
fails.

| file | what it checks | data |
|---|---|---|
| `sep_test_path` | no function name is defined twice on the project path (compared without case, vendored code included, the known vendored duplicates listed), none shadows or is shadowed by a MATLAB or toolbox function, and, when `D:\dendrites\code` is there, none is also defined by the imaging repository, which runs in the same MATLAB | none |
| `test_backward_compat` | the cohort table lists the 17 adults in their fixed order and groups, positions 1 to 17 are still those adults, their folders exist; `get_atlas('ccf')` and `get_paths` give the folders and files the drivers used; the young cohort `get_cohort` returns is the young rows of `common\cohort.csv`, in their order (read from the table, so the test needs no change when a brain is added); each young brain's age matches the age in its name; each young brain's `share_subdir` holds its `.czi` files (skipped when the share is not reachable; the share is only read) | the data root of `get_paths`, the lab share |

Run `sep_test_path` after adding, moving or renaming a function, and
`test_backward_compat` after any change to the cohort table, `get_cohort`,
`get_atlas` or `get_paths`.

Other checks live with their code:

- the automatic annotation's engine: `registration\auto_annotation\setup.ps1`
  ends with a self-test that loads both models in its environment;
- the plasticity comparison, the Python route and the registration are
  checked against their earlier outputs with the tools of
  [`../tools/README.md`](../tools/README.md), on copies of the data, never by
  a test on the production data;
- the Python route has no tests yet. They go in `mapping/tests/`, pytest
  run from `mapping\` as `python -m pytest tests` in `tools\venv_dev`
  (`docs/STYLE.md`).

## Adding a test

One script per topic, `test_<topic>.m`, with the header and the checks of
`test_backward_compat`: each check prints PASS or FAIL with what it
compared, and the script ends with an error when one failed. A good test has
a known answer and fails on the broken code: when fixing a bug, write the
test first and check that it fails. Use synthetic data where possible; a test
that needs data which is not connected says SKIP, not FAIL. A test never
writes into the data root: it uses a temporary folder, or `SEP_DATA_ROOT`
pointed at a copy.
