# Tests of the Python route

Known-answer checks of the ISH analysis on synthetic data, plus a few checks of
today's tables, which are skipped when the data root does not hold them. Run
from `mapping\` with the dev environment, which reads the packages of
`tools\venv_atlas`:

```
cd mapping
..\tools\venv_dev\Scripts\python -m pytest tests
```

With `SEP_DATA_ROOT` set to a data tree where the run scripts have run, the
checks of today's tables run too; without it they are skipped.

| file | what it checks | data |
|---|---|---|
| `test_structures.py` | the declared set keeps a grey structure measured in all ten adults and drops one seen in nine, a fibre tract and a "..., unassigned" label, each with its reason; a nano mean below background in one adult takes a structure out; zref has median 0 and p90 - p10 of 1 over its reference in every brain, and a structure added outside the reference moves no other zref; on a symmetric label volume a bilateral structure's centroid lies in its hemisphere with one hemisphere and on the midline with both; erosion leaves the 3 x 3 x 3 core of a 5 x 5 x 5 cube. Today's set: 204 structures, each with a centroid | synthetic; `structure_set.csv`, `centroids.csv` |
