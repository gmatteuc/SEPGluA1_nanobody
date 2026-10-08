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
| `test_section_qc.py` | on a grid with a tenfold gradient along AP: a section planted ten times dimmer is flagged and set missing, the gradient itself is not; a dim section on the exceptions list is kept as absence; a sagittal series is judged along ML (a dim ML plane flagged along ML, nothing along AP); sections not flagged keep their values to the bit and missing voxels stay missing; an end section with too little brain and a section with no data are not judged; the exceptions file's statuses (proposed and accepted keep, rejected does not, anything else refused) | synthetic |
| `test_gene_table.py` | two panels sharing an experiment give one row for it, with both labels, and a repair row; a grid outside the reference box and a missing one are excluded with their reason; a gene measured once has a profile in its experiment's order; GO ancestry follows is_a and part_of (and alt ids) from an obo file; NOT annotations are left out; the gene sets follow their rules (panel roles, post- or presynaptic but not both, marker lists). Today's table: 451 genes, one row per experiment, all of P9's genes, every exclusion with a reason, P9's category on P9's genes only; with no section set missing, 15 genes' reliabilities equal `gene_reliability.csv` of 5 October | synthetic; `gene_table.csv`, `gene_documentation.csv`, the grids, `adult_v2/ish/gene_reliability.csv` |
| `test_spatial_null.py` | on 200 random points in a 10 mm box, maps like the adult one (a smooth Gaussian random field, length 4 mm, plus 30% structure-level noise): over the matched range the surrogates' mean variogram is within 15% of the map's in the median of 8 maps and always at least three times closer than a shuffled map's; resampled surrogates are permutations of the map's ranks; a seed fixes the surrogates; a map with a missing value is refused; a gene that is the map plus noise (rho about 0.5) gets p < 0.05 in at least 9 of 10 maps; of 100 pairs of independent maps, the spatial p is below 0.05 in at most 10 (5 expected) and the ordinary p in at least 25; the p counts the observed value once; the vectorised null equals scipy's Spearman. 500 surrogates per map, about 30 s | synthetic |

The tests take about a minute in all.
