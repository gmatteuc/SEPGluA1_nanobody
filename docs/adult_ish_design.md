# Adult nano characterisation and the ISH gene comparison — design

**Status: design agreed 25 Sep 2026, implementation starting. Nothing here is built yet.**
This document is the running record of *what we decided and why*, plus a map of the
pipeline. Update it as steps land; if the code and this file disagree, the file is wrong
and should be fixed in the same commit.

## The question

Two halves of one story:

1. **Characterise the adult distribution of surface GluA1** across the whole brain — the
   thing the old `merged_naive_rws` outputs did.
2. **Compare that spatial pattern against a large panel of genes** and ask *which kinds of
   genes match it best*. This is also the strongest validation we have that the nanobody
   reports the SURFACE pool rather than total receptor, because the answer is already
   suggestive: the nano map correlates better with AMPAR trafficking and anchoring
   machinery (Cacng8 = TARP γ-8 at ρ = 0.78, Dlg2 = PSD-93 at 0.76) than with *Gria1*
   itself (ρ = 0.52, rank 22 of 97). See `project_measurement_validation_strategy` in
   memory.

## Where this sits

It reuses the v2 chain rather than rebuilding it. The adult per-mouse volumes, tissue
masks, background subtraction and cohort means are **already computed and audited**:
`comparisons_v2/per_mouse/*.npz` (10 adults, `sig`/`auto`/`sep`/`tissue`),
`comparisons_v2/ccf/{adult,naive,rws}/` (five readings), and the per-mouse region table.
Nothing in `comparisons_v2/` or `comparisons/` is written to. New outputs go to
`data/adult_v2/`.

Naming is `v2_adult_*` / `v2_ish_*` deliberately: this is the same generation of code
answering a different question, and the two should converge into one codebase.

## What the old route got wrong, and must not be inherited

From `project_p8_p9_pipeline` and the audit of the v2 work:

| problem | consequence | fix |
|---|---|---|
| P6bis maps each mouse with an **affine** (slope *and* intercept) | ratios between regions do not survive it | use the v2 per-mouse volumes, which are a pure scale |
| `abs()` in P8 over zero-filled planes | unreached voxels became bright slabs | v2 tissue mask + n maps |
| tissue mask taken from the **nano** channel | circular, and a partial section reads as all-background | v2 mask (auto channel, + 4 MAD) |
| group `.mat` files silently stale | nano and auto from different registrations | read the registered tiffs |
| within-brain z-score as the only reading | destroys absolute scale | the five explicit readings, each with its blindness documented |
| 100 hand-picked genes, ad-hoc categories | no background, no null | see below |
| naive p-values on map-to-map correlation | ~875-fold false-positive inflation in mouse | ensemble null, see below |

## Decisions and why

**D1 — Region-level correlation is primary.** Voxel-level Pearson stays as a secondary
diagnostic. The null model is defined over regions, so that is where inference happens.

**D2 — Keep boundary protection when aggregating ISH to regions.** ISH is 200 µm and
registration is imperfect, so a region mean picks up its neighbours. P9's eroded and
distance-weighted masks were a good idea; keep one, having *tested* which matters rather
than assuming.

**D3 — The null randomises the MAP, not the gene annotations.** Fulcher, Arnatkeviciute &
Fornito 2021 show the conventional random-gene null inflates false positives **875-fold in
mouse**. They recommend comparing against an ensemble of spatially autocorrelated
surrogate *phenotype* maps (their "SBP-spatial"), built with a spatial-lag model; for
mouse they use ρ = 0.8 and a characteristic length d₀ = 1.46 mm. The question then becomes
the one we actually mean: *are genes in this category more correlated with the nano map
than with a random map of the same spatial smoothness?*

*Why not brainSMASH:* same idea (surrogate maps with matched autocorrelation), but built
for human neuroimaging with surface/volumetric distances. The spatial-lag ensemble has
published mouse parameters and a validated false-positive analysis on Allen ISH, which is
our exact setting. Revisit if the lag model fits badly.

**D4 — The null must be shown to be calibrated, not assumed.** A random phenotype should
come out significant about 5% of the time. Producing that curve — and the conventional
null's much worse one beside it — is a required diagnostic, not an optional one.

**D5 — The gene panel becomes a background, not a universe.** Hand-picking 100 genes is
the setup the Fulcher paper warns about. The Allen coronal set is ~4,000 genes at ~1 MB
each (~4 GB), and the existing 100 become *one category among many* that we test like any
other. Sagittal (~20,000 genes, left hemisphere only) is a possible extension because our
nano data is hemisphere-folded; check quality before relying on it.

**D6 — Gene sets are defined before looking.** SynGO (Koopmans 2019, PDF in
`data/ref_papers/`) for curated synaptic terms with a proper background, plus GO. The
existing hand-made categories are kept for comparison so we can say what curation changed.
Defining "trafficking" after seeing the ranking would destroy the argument's value.

**D7 — Three arms, predicted in advance.** Run the identical analysis for `nano`, `sep`
and `nano/sep`:

| quantity | predicted best correlate |
|---|---|
| SEP (total receptor) | *Gria1* expression |
| nano (surface pool) | both, intermediate |
| nano / SEP (surface fraction) | trafficking machinery, **not** *Gria1* |

**D8 — Reproduce before improving.** First milestone is the new code reproducing the old
ranking (Cacng8 top, Gria1 ≈ 22nd). Only then do the improvements land, so every later
difference is an improvement rather than a bug.

**D9 — ISH quality control is a rule applied before correlations are seen.**
`n_bad_sections` is already recorded per gene (28 for Cacng8) and was never used. Decide
thresholds, apply, document how many genes drop and why.

## Pipeline map

Nothing below is implemented yet; steps are ticked as they land.

```
  [ existing, reused, not rewritten ]
  comparisons_v2/per_mouse/*.npz         10 adults: sig, auto, sep, tissue
  comparisons_v2/ccf/{adult,naive,rws}/  cohort mean / sd / n, five readings

  [ new ]
  v2_ish_prepare.py     Allen ISH grids -> a checked, region-aggregated table
                          . read MetaImage (67 x 41 x 58 @ 200 um, RAI)
                          . verify orientation and resolution against the atlas
                          . aggregate to structures with boundary protection
                          . per-gene QC flags
                        -> adult_v2/ish/gene_region_table.parquet

  v2_adult_profile.py   the adult distribution itself
                          . per-structure profile, one dot per mouse, CI
                          . ranked tables, by division and macro group
                          . cortical laminar profile
                          . naive vs rws as the built-in negative control
                        -> adult_v2/profile/

  v2_ish_null.py        the spatial-lag ensemble and its calibration test
                        -> adult_v2/ish/null/

  v2_ish_compare.py     nano / sep / nano-sep against every gene
                          . region Spearman, ensemble p, FDR
                          . category enrichment (SynGO, GO, the old panel)
                        -> adult_v2/ish/

  v2_adult_diagnostics.py   one sheet per step, as in v2
                        -> adult_v2/processing_diagnostics/
```

## References

In `data/ref_papers/` (PDFs) — summaries in `PAPER_SUMMARIES.md` there.

- **Fulcher, Arnatkeviciute & Fornito 2021**, *Overcoming false-positive gene-category
  enrichment in the analysis of spatially resolved transcriptomic brain atlas data*,
  Nat Commun 12:2669. `Fulcher_2021.pdf`. The null-model argument, D3 and D4.
- **Koopmans et al. 2019**, *SynGO: an evidence-based, expert-curated knowledge base for
  the synapse*, Neuron 103:217. `Koopmans_2019.pdf`. The gene sets, D6.
- Allen Mouse Brain Atlas ISH: ~20,000 genes, all with sagittal, a subset coronal;
  200 µm grid, expression *energy*. <https://brain-map.org/support/documentation/api-for-mouse-brain-atlas>
- GCEA false-positives toolbox (the authors' implementation):
  <https://brain-map.org/support/community-tools/gcea>
- Burt et al. 2020, brainSMASH — considered and set aside, see D3.
