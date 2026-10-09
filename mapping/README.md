# Mapping: the Python route

Measures the nanobody signal of every registered brain on the atlas of its own
age, and asks three questions of it: how it is distributed in the adult brain,
what it measures against the Allen in situ hybridisation maps, and how young
brains differ from adults. It reads the registered volumes that
`../registration/` writes and writes only under `<data>\comparisons_v2\`,
`<data>\adult_v2\`, the ISH grid cache `<data>\atlas_ish\` and the
third-party data it downloads, `<data>\reference\`.

## Layout

```
mapping/
  run_*.py              one entry point per step; the header of each lists the run order
  settings.toml         every analysis parameter, a table per stage
  gene_targets.csv      the hand-written 100-gene panel of P9 (April)
  ish_section_exceptions.csv  dim ISH sections kept as true absence, with reason and status
  tests/                pytest, known answers on synthetic data (tests/README.md)
  sepmap/               the package the run scripts call
    config.py           the data root (SEP_DATA_ROOT and its guards), the settings, the cohort table
    plotting.py         the palette, the style, the colormaps, the save function, the coronal drawing
    structures.py       the declared structure set (A1), its zref, one-hemisphere centroids
    hemispheres.py      the two hemispheres of a volume folded onto one (volumes, young_vs_adult)
    diagnostics.py      the sheets that audit each step
    volumes/            per_mouse, to_ccf, cohort: per-brain volumes, the adult CCF, cohort means
    young_vs_adult/     compare, replot, region_plot, region_groups, video, video_compare,
                        closeup
    adult/              profiles (per adult and structure, the channels and zref),
                        synaptome (the measured synapse density, Zhu 2018),
                        beyond_density, beyond_controls, beyond_calibration,
                        beyond_regression, beyond_figures, sep_channel_check
    ish/                panel_build, panel_fetch, regions, reliability (the inputs);
                        section_qc, gene_table, spatial_null, gene_ranking, robustness,
                        divisions, gene_sets, top_genes, overview (the ISH analysis);
                        plotting (its figures); compare, words, roles, panel_test
                        (the run of 5 October, retiring)
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
- `run_panel_build`, `run_panel_fetch`, `run_ish_section_qc`,
  `run_ish_gene_table` and `run_synaptome` call the network for what their
  caches lack; `--offline` makes the last three stop instead.

| stage | steps | what |
|---|---|---|
| volumes | 1 `run_per_mouse`, 2 `run_to_ccf`, 3 `run_cohort` | per brain, the tissue mask and the background-subtracted channels; every brain on the adult CCF grid; per voxel, the cohort mean, SD and n of every reading, over the whole brain and with each brain's hemispheres averaged first (for the videos' t) |
| young against adult | 4 `run_compare`, 5 `run_region_plot`, 6 `run_region_groups`, 7 `run_video`, 8 `run_video_compare`, 9 `run_closeup` | maps and the per-structure table; the region statistics, measured on each brain's own atlas, which are the numbers to quote; systems and layers; videos; a coronal plane and the cortical flatmaps (`run_replot` redraws step 4's maps) |
| diagnostics | 10 `run_diagnostics` | one sheet per question, so each step can be checked by eye |
| ontology panel | 11 `run_panel_build`, 12 `run_panel_fetch` | a 390-gene panel from Gene Ontology terms, and its ISH grids (network, once; not rerun for an analysis) |
| ISH inputs | 13 `run_structure_set`, 14 `run_ish_section_qc`, 15 `run_ish_gene_table`, 16 `run_ish_spatial_null` | the declared structures and the adult profiles (A1); section QC of every experiment (A2); one gene table, merged profiles, gene sets, documentation (A9); surrogate maps with the map's smoothness (A7) |
| the genes against the map | 17 `run_ish_gene_ranking`, 18 `run_ish_robustness`, 19 `run_ish_divisions`, 20 `run_ish_gene_sets` | each gene against the map and the autofluorescence map, with the null (A8); the ranking under other choices (A3); between or within divisions (A6); gene sets and localisation against matched controls |
| the measured synapse density | 21 `run_synaptome` | fetches the PSD95 and SAP102 punctum densities of Zhu et al. 2018 as Hansen et al. share them, from their repository at a pinned commit (network, once; `--offline` stops instead), into `<data>\reference\synaptome\`; places the samples in the CCF structures; how much of the declared set and of analysis 4's fit the PSD95 density covers, and how it agrees with the mRNA density terms, Gria1 and the map |
| beyond Gria1 expression and synapse density | 22 `run_beyond_density`, 23 `run_beyond_controls`, 24 `run_beyond_calibration`, 25 `run_beyond_regression`, 26 `run_beyond_figures` | how much of the map the main model of 8 October (Gria1, synapse density, autofluorescence) predicts, and the genes against what it leaves; seven controls, and the leftover under other folds, structures and the check rows of `[beyond.variants]` (PSD95 or SAP102 puncta, the four subunits); the same model on maps whose answer is known, and the nano map against them on the same structures; the regression per structure; the figures of part 1 |
| the genes that follow the map | 27 `run_ish_top_genes` | the genes past the map's null, Gria1, Cacng8 and the AMPA receptor complex family described on every axis of the ISH line, each added to the main model against maps of its smoothness; the tests named on 8 October for the leftover (Cacng8, then the family as a group and gene by gene); per-gene sheets |
| the green channel, the overview | 28 `run_sep_channel_check`, 29 `run_ish_overview` | what the green channel reports; April's headline, the numbers for the text, the overview and the index of the figures |

What fixes the order: step 5's `region_means_per_mouse.csv` is read by steps
13 (the stored `cref` and `zref`, to measure what the declared reference
moves), 18 (the stored `zref` and `ratio` rows) and 23 (the readings of
control G); step 13's structure set, profiles and centroids by every later
step; step 14's flags by 15; step 15's gene table and profiles by 17 to 28;
step 16's surrogates by 17, 19 and 20; step 17's `gene_ranking.csv` and
`null_rho.npz` by 18 to 20, 27 and 29; steps 18 and 19's tables by 27; step
21's synapse density by 22 to 27; step 22's tables by 23 to 27, and steps 23 to
25's by 26; step 29 reads every step's numbers, so it comes last. Sheet 07
of step 10 reads the tables of steps 4 and 5, and sheet 08 the background
masks of `../group_comparison/run_normalise_groups` (naive and P20), so in a
full run the plasticity chain comes first. The scripts of the ISH run of 5
October (`run_ish_regions`, `run_ish_compare`, `run_ish_words`,
`run_ish_roles`, `run_ish_reliability`, `run_ish_panel_test`) are out of the run
order: they write the frozen folder `adult_v2\ish\`, and move to `archive/`
next. `run_ish_arms` and `run_adult_arms` are in `../archive/` already (9
October): the green channel is mostly autofluorescence, so the channel ratios
they tested against the genes have no premise.

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
- **The ISH analysis** (steps 13 to 29), on the declared structures (grey
  matter measured in all ten adults) and one gene table: how much of the adult
  map Gria1 expression and synapse density predict, against the map's own
  reliability and a calibration floor, with seven controls (part 1); each
  gene, each gene set fixed in advance and the localisation genes against the
  map, with a spatial null of surrogate maps, between and within divisions,
  and on the autofluorescence map, and the genes that follow the map most
  described one by one, with the tests named on 8 October for what the main
  model leaves (part 2); what the green channel reports.
  Synapse density is also measured, not only read from mRNA: PSD95 puncta per
  structure in one adult mouse (Zhu et al. 2018), placed in the CCF structures
  by Allen id (`run_synaptome`). They cover 77 of the 126 structures of the
  fit, under the 80% the rule of 8 October asks, so the main model keeps the
  mRNA terms and PSD95 is a check row. Method, figures and results:
  `../docs/ISH_ANALYSIS.md`.
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

- `ish_analysis\`: the ISH analysis. `tables\` (every table of steps 13 to
  29, the surrogates, `numbers_<step>.csv` and `numbers_for_the_text.csv`),
  `synaptome\` (the measured synapse density per sample and per structure,
  its coverage and agreement), `beyond\` (part 1's tables and working
  figures), `green_channel\`, `figures\` (the guided figures
  `00_overview.png` to `16_april_headline.png`, each main figure with its
  detailed versions, `03s1_beyond_budget.png` and so on, the index `README.md`,
  `qc\`, `genes\` and `top_genes\`), `cache\` (Allen
  experiment lists, mygene records, `go-basic.obo`). The table of every file:
  `../docs/ISH_ANALYSIS.md`, section 10.
- `panel\`: `panel_v2.csv`, `panel_genes.csv`, `fetch_failures.csv`, the API
  cache `cache\` (steps 11 and 12)
- `ish\`, `arms\`, `beyond\`: the run of 5 October, frozen and not to quote
  (`../docs/FIGURES.md`, Superseded outputs)

Under `<data>\reference\`: third-party data, each folder with a
`fetch_log.txt` (source, version, citation, checksums); `synaptome\` is
written by `run_synaptome` (step 21).

Figures are PNG with an EPS beside most of them. Caches are reused when
present: `per_mouse\*_scalars.npz` when its recorded date matches its source,
the ISH tables, the API answers and the grids.

## Known limitations

- `sepratio` is not a surface fraction. The green channel tracks the
  autofluorescence in these sections (rho 0.72 to 0.83 in every adult,
  `adult_v2\ish_analysis\green_channel\sep_channel_check.csv`), so nano/SEP
  behaves as a second nano over autofluorescence.
- The level of `ratio` depends on exposure. Its zero is a structure as bright
  in nano as in autofluorescence (`region_plot`'s title: "0 = equally
  bright"), but only at the two channels' exposures, so it is no biological
  null; only the order of structures within a brain carries meaning.
- Young and adult brains were imaged in different sessions, and the
  autofluorescence rises with age, so only patterns compare across ages, not
  levels.
- In the young-against-adult region tables and maps, zref's median and
  spread are taken over the structures every brain in the tables has, so a
  new brain can move every brain's zref, the adults' included. The ISH
  analysis takes them over its declared set (A1); the young-against-adult
  tables move to it with A1's second half.
- The spatial null of the ISH analysis rests on structure centroids in one
  hemisphere and a variogram matched at short range; structures of very
  different sizes count as points. It is calibrated on random maps of another
  kind (about 4% false positives at 0.05). A difference between two genes'
  correlations is tested against maps related to both alike, and still needs
  to be large to pass.
- The calibration floor of part 1 errs low: a gene measured by one Allen
  experiment is the same in both halves, so its mismatch is not in it.
- The measured synapse density is one adult mouse, sampled in a few coronal
  sections: it covers 77 of the 126 structures of part 1's fit and 96 of the
  204 declared ones (8 October), and a region the source gives only above
  several structures (PTLp, the midbrain's motor part) is never spread onto
  them. As shared, each of its 37 punctum subtypes is scaled to 0..1, so the
  PSD95 density is a mean of subtype maps, not a count of puncta.
- The true absences of the section QC (`ish_section_exceptions.csv`, one
  section today) and the sections kept at a step in expression are a human
  call, proposed and not yet reviewed.
- The ISH outputs quoted in the documents were made on a full copy of the
  production inputs (`SEP_DATA_ROOT`); the production data root has them once
  steps 13 to 29 run there after the merge.
- The young cohort holds P16, P20 and P22 brains. A brain of another age
  (`mapping_cohort` `young_P28` and so on) goes through the per-brain steps,
  but no cohort takes it until one is added (see `../docs/ADDING_DATA.md`,
  step 4).
- The route does not yet answer everything the MATLAB scripts of
  `../adult_matlab/` answer (enrichment calls per structure, the
  autofluorescence control per structure: A4, A5), so they keep running
  until it does.
- The tests (`tests/`) cover the ISH analysis; the per-brain steps and the
  young-against-adult tables have none yet.
