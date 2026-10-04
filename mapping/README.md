# Mapping: the Python route

Measures the nanobody signal of every registered brain on the atlas of its own
age, and asks three questions of it: how it is distributed in the adult brain,
what it measures against the Allen in situ hybridisation maps, and how young
brains differ from adults. It reads the registered volumes that
`../registration/` writes and writes only under `<data>\comparisons_v2\`,
`<data>\adult_v2\` and the ISH grid cache `<data>\atlas_ish\`.

## Layout

```
mapping/
  run_*.py              one entry point per step; the header of each lists the run order
  settings.toml         every analysis parameter, a table per stage
  sepmap/               the package the run scripts call
    config.py           the data root (SEP_DATA_ROOT and its guards), the settings, the cohort table
    plotting.py         the palette, the colormaps, the save function, the coronal drawing
    diagnostics.py      the sheets that audit each step
    volumes/            per_mouse, to_ccf, cohort: per-brain volumes, the adult CCF, cohort means
    young_vs_adult/     compare, replot, region_plot, region_groups, video, video_compare,
                        closeup, hemispheres
    adult/              sep_channel_check, arms, beyond_density, beyond_controls,
                        beyond_figures, beyond_regression
    ish/                regions, compare, words, roles, arms, panel_build, panel_fetch,
                        reliability, panel_test
```

## Setup and run

- `tools\venv_atlas` runs every script but `run_closeup.py`, which needs
  `ccf_streamlines` and runs in `tools\venv_flat`. Their pinned requirements,
  and how to create them, are `tools/requirements_atlas.txt` and
  `requirements_flat.txt`. CCF Translator, in `venv_atlas`, downloads its
  DeMBA-to-CCF deformation fields on first use; the copy the results were
  made with is kept in `<data>\atlas\ccf_translator_fields\`.
- Run from the code root, in the order of the scripts' headers:

  ```
  tools\venv_atlas\Scripts\python.exe mapping\run_per_mouse.py
  ```

- Every run prints the data root, where it came from and the settings in
  force before anything else. `SEP_DATA_ROOT` moves the data root;
  `V2_READINGS`, a comma-separated subset of the readings, limits one run to
  them.
- Run options are command-line options (`--help` on each script): mice for
  the per-brain steps, cohorts for `run_video`, readings and colour ranges for
  the figures, `--panel targets|ontology` for the two ISH passes. Analysis
  parameters are in `settings.toml` only.
- `run_panel_build`, `run_panel_fetch` and `run_ish_words` call the network
  for what their caches lack.

| stage | steps | what |
|---|---|---|
| volumes | 1 `run_per_mouse`, 2 `run_to_ccf`, 3 `run_cohort` | per brain, the tissue mask and the background-subtracted channels; every brain on the adult CCF grid; per voxel, the cohort mean, SD and n of every reading, over the whole brain and with each brain's hemispheres averaged first (for the videos' t) |
| young against adult | 4 `run_compare`, 5 `run_region_plot`, 6 `run_region_groups`, 7 `run_video`, 8 `run_video_compare`, 9 `run_closeup` | maps and the per-structure table; the region statistics, measured on each brain's own atlas, which are the numbers to quote; systems and layers; videos; a coronal plane and the cortical flatmaps (`run_replot` redraws step 4's maps) |
| diagnostics | 10 `run_diagnostics` | one sheet per question, so each step can be checked by eye |
| ISH, 100-gene panel | 11 `run_ish_regions --panel targets`, 12 `run_ish_compare`, 13 `run_ish_words`, 14 `run_ish_roles` | Allen ISH energy per gene and structure; each gene against the adult map; the annotation words and the curated roles of the ranking |
| channel arms | 15 `run_adult_arms`, 16 `run_ish_arms`, 17 `run_sep_channel_check` | SEP/auto, nano/auto and nano/SEP per adult; against the genes; what the green channel reports |
| ISH, 390-gene panel | 18 `run_panel_build`, 19 `run_panel_fetch`, 20 `run_ish_regions --panel ontology`, 21 `run_ish_reliability`, 22 `run_ish_panel_test` | a panel from Gene Ontology terms; its grids; how reliable one ISH map is; localisation genes against controls |
| beyond abundance | 23 `run_beyond_density`, 24 `run_beyond_controls`, 25 `run_beyond_figures`, 26 `run_beyond_regression` | how much of the adult map receptor abundance and synaptic density explain, the controls, the figures |

What fixes the order: step 5's `region_means_per_mouse.csv` is read by steps
12, 14, 15, 22 and 23 to 26; pass 1's gene table by 12, 14, 16 and 17; step
21's merged gene table by 22 to 26; step 25 draws its panel D from step 24's
`controls.csv`. Sheet 07 of step 10 reads the tables of steps 4 and 5, and
sheet 08 the background masks of `../group_comparison/run_normalise_groups`
(naive and P20), so in a full run the plasticity chain comes first.

## Assumptions and preprocessing

- Every brain is read from its registered tiffs (`volume_registered\`
  `chan02_NANO`, `chan03_AUTO`; `volume_registered_sep\chan02_SEP`), on the
  atlas of its age: `ccf` for the adults, `demba_p<age>` for a young brain.
  A brain without its SEP volume stops `run_cohort`.
- Per brain, at 20 um (2 x 2 x 2 block means of the registered grid): the
  tissue mask is the autofluorescence above its off-tissue level (median + 4
  MAD, `tissue.mad_k`), where the nano is non-zero and inside the atlas
  brain. The nano intensity never enters the mask. Each channel has its own
  off-tissue background subtracted.
- Each brain is scaled on its own (the readings below), not fitted onto a
  group as in the plasticity comparison.
- Region statistics are computed in each brain's own atlas, whose labels
  share the CCF's parcellation index, so nothing is warped. Only the maps
  carry a young brain into the CCF, with CCF Translator at its age; an adult
  is placed, never warped.
- The cohorts are those of the cohort table's `mapping_cohort` column: young
  (P16, P20 and P22 pooled), young_P20, young_P16, young_P22, naive, rws,
  and adult (naive and rws). `sepmap/volumes/cohort.py` builds them.

## Analyses

- **Readings**, defined in `sepmap/volumes/cohort.py`: `ratio` (nano per unit
  autofluorescence), `sepratio` (nano per unit SEP), `cref` (over the brain's
  isocortex mean), `subref` (over its subcortex mean), `zref` (log2 of
  `cref`, minus the brain's median structure, over its p90 - p10 spread).
  None replaces another.
- **Young against adult**: per structure and reading, Welch and Mann-Whitney
  tests with their Benjamini-Hochberg q, the P20-only and naive-RWS contrasts
  beside them (`region_stats.csv`); the same by system, division and layer
  (`group_stats.csv`).
- **Adult distribution and the ISH comparison**: the adult map against each
  gene's ISH map, over structures; the reliability of single ISH experiments;
  the localisation genes against expression-matched controls with a
  permutation null; the variance of the map that abundance, synaptic markers,
  postsynaptic-density genes and autofluorescence explain, against the map's
  own reliability, with seven controls.
- **The green channel**: what it tracks across structures, per adult.
- The results and what they mean: `../docs/SCIENTIFIC_CONTEXT.md`.

## Outputs

Under `<data>\comparisons_v2\`:

- `per_mouse\<mouse>.npz` and the cache `<mouse>_scalars.npz`;
  `per_mouse_ccf\<mouse>.npz`
- `ccf\<cohort>\<reading>_{mean,sd,n}.npy`, `<reading>_folded_{mean,sd,n}.npy`
  (each brain's hemispheres averaged first, the left half; `run_video` stops
  without them), `mice.txt`, the cohort videos
- `young_vs_adult\`: `volumes_ccf20.npz`, `region_table.csv`,
  `cortex_table.txt`, `region_means_per_mouse.csv`, `region_stats.csv`,
  `group_stats.csv`, the slice, region, group and laminar figures, the
  side-by-side videos, the close-ups and flatmaps (`detail_*`)
- `processing_diagnostics\`: sheets 01 to 09 and the index `README.md`

Under `<data>\adult_v2\`:

- `ish\`: the gene tables of both passes and their drops, `gene_correlations.csv`,
  `feature_enrichment.csv` (with the cache `annotation\`), `gene_roles.csv`,
  `role_summary.csv`, `gene_reliability.csv`, `gene_region_table_merged.csv`,
  `panel_test.csv`, and their figures
- `panel\`: `panel_v2.csv`, `panel_genes.csv`, `fetch_failures.csv`, the API
  cache `cache\`
- `arms\`: `region_means_arms.csv`, `arm_gene_correlations.csv`,
  `sep_channel_check.csv`, and their figures
- `beyond\`: `structures_used.csv`, `variance_partition.csv`,
  `residual_by_structure.csv`, `controls.csv`, figures; `for_sami\`, the
  figure panels with `numbers_for_the_caption.txt` and
  `regression_table.csv`

Figures are PNG with an EPS beside most of them. Caches are reused when
present: `per_mouse\*_scalars.npz` when its recorded date matches its source,
the ISH tables, the API answers and the grids.

## Known limitations

- `sepratio` is not a surface fraction. The green channel tracks the
  autofluorescence in these sections (rho 0.79 across the ten adults, header
  of `registration/run_add_sep_channel.m`; per mouse in
  `adult_v2\arms\sep_channel_check.csv`), so nano/SEP behaves as a second
  nano over autofluorescence.
- The level of `ratio` depends on exposure. Its zero is a structure as bright
  in nano as in autofluorescence (`region_plot`'s title: "0 = equally
  bright"), but only at the two channels' exposures, so it is no biological
  null; only the order of structures within a brain carries meaning.
- Young and adult brains were imaged in different sessions, and the
  autofluorescence rises with age, so only patterns compare across ages, not
  levels.
- In the region tables, zref's median and spread are taken over the
  structures every brain in the tables has. A new brain can change that set
  and so move every brain's zref, the adults' included. The declared
  structure set of the refactor plan (A1) replaces it.
- The ISH gene ranking is descriptive. Its p and q values are
  anticonservative: genes within a set are co-expressed and brain maps are
  spatially autocorrelated (Fulcher 2021). No spatial null is built yet (A7).
  The ranking includes structures seen in only a few adults until A1 to A3
  restrict it.
- The young cohort holds P16, P20 and P22 brains. A brain of another age
  (`mapping_cohort` `young_P28` and so on) goes through the per-brain steps,
  but no cohort takes it until one is added (see `../docs/ADDING_DATA.md`,
  step 4).
- The route does not yet answer everything P8, P9 and P10 answer (enrichment
  calls per structure, the autofluorescence control per structure: A4, A5),
  so those three stay at the code root until it does.
- There are no tests of this package yet.
