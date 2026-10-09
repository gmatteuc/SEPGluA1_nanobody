# Tests of the Python route

Run from `mapping\` with the development environment:

```
..\tools\venv_dev\Scripts\python.exe -m pytest tests
```

The tests read the data root as every run does (`SEP_DATA_ROOT`, else the
folder beside the code); the one test on real data is skipped when its table
is not there.

| file | what it checks | data |
|---|---|---|
| `test_adult_layers.py` | the adult map by depth (`adult.layers`, `adult.layers_plotting`): layers read from both namings of the CCF labels, layers pooled into bands and L6, the published Harris 2019 scores and the order of the areas, a hierarchy table missing an area refused, the reason a cell is missing, half against half agreement on a shared profile and on noise, mean, SEM and t with `min_mice`, the laminar contrast, the check against the region table, the per-pixel mean and SD of the band maps, the bars' grey; and two close-up helpers it shares (`closeup.smooth_within`, `closeup.band_average`) | synthetic, the hierarchy table beside the code; the last test reads `adult_v2\layers\area_layers_per_mouse.csv` |
