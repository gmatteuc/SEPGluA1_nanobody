# Tests

Checks of the shared definitions that every pipeline relies on, where a
mistake would be silent: the MATLAB path and the cohort table; and a
known-answer check of the plasticity comparison's permutation test, on
synthetic data. Each runs in a fresh MATLAB session, from the code root:

```matlab
restoredefaultpath; cd('D:\sep_histology\code'); sep_setup_paths
run(fullfile('tests', 'sep_test_path.m'))
run(fullfile('tests', 'test_backward_compat.m'))
run(fullfile('tests', 'test_region_permutation.m'))
```

`sep_setup_paths` does not put `tests\` on the path, so each test is run by
its file. Each prints what it checks and stops with an error when a check
fails.

| file | what it checks | data |
|---|---|---|
| `sep_test_path` | no function name is defined twice on the project path (compared without case, vendored code included, the known vendored duplicates listed), none shadows or is shadowed by a MATLAB or toolbox function, and, when `D:\dendrites\code` is there, none is also defined by the imaging repository, which runs in the same MATLAB | none |
| `test_backward_compat` | the cohort table lists the 17 adults in their fixed order and groups, positions 1 to 17 are still those adults, their folders exist; `get_atlas('ccf')` and `get_paths` give the folders and files the drivers used; the young cohort `get_cohort` returns is the young rows of `common\cohort.csv`, in their order (read from the table, so the test needs no change when a brain is added); each young brain's age matches the age in its name; each young brain's `share_subdir` holds its `.czi` files (skipped when the share is not reachable; the share is only read) | the data root of `get_paths`, the lab share |
| `test_region_permutation` | the permutation test of the plasticity comparison's region measures (`region_permutation_test`) on synthetic stacks of 5 against 5 mice: all 252 splits enumerated with the observed one among them; the observed split gives the bars' rolling median and share to the bit (the arithmetic of `group_differences`, written out in the test); its sum, 99th percentile, top volume and heaviest cluster equal those written out region by region, each region's clusters from `bwconncomp` of its own voxels of one sign, also on a map whose clusters cross a region border and touch the other sign and whose heaviest cluster differs at 6, 18 and 26-connectivity; a split's mirror image equals its own scores with the signs swapped; a strong 400-voxel bump inside a large region and a shift over a small one: cluster and top volume put the large region first at a corrected p < 0.05 (the share and the 99th percentile put the small one first); pure noise under six seeds: no measure reaches a corrected p < 0.05 in more than 3 of 12 maps. Opens a pool of 4 thread workers, about a minute | none |

Run `sep_test_path` after adding, moving or renaming a function,
`test_backward_compat` after any change to the cohort table, `get_cohort`,
`get_atlas` or `get_paths`, and `test_region_permutation` after any change to
`region_permutation_test` or to the t, the surprise or the rolling median of
`group_differences`.

Other checks live with their code:

- the automatic annotation's engine: `registration\auto_annotation\setup.ps1`
  ends with a self-test that loads both models in its environment;
- the plasticity comparison, the Python route and the registration are
  checked against their earlier outputs with the tools of
  [`../tools/README.md`](../tools/README.md), on copies of the data, never by
  a test on the production data;
- the Python route's tests are in `mapping/tests/` (pytest, run from
  `mapping\` as `python -m pytest tests` in `tools\venv_dev`;
  [`../mapping/tests/README.md`](../mapping/tests/README.md)). They cover the
  ISH analysis; the per-brain steps have none yet.

## Adding a test

One script per topic, `test_<topic>.m`, with the header and the checks of
`test_backward_compat`: each check prints PASS or FAIL with what it
compared, and the script ends with an error when one failed. A good test has
a known answer and fails on the broken code: when fixing a bug, write the
test first and check that it fails. Use synthetic data where possible; a test
that needs data which is not connected says SKIP, not FAIL. A test never
writes into the data root: it uses a temporary folder, or `SEP_DATA_ROOT`
pointed at a copy.
