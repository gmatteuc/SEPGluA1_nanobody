# Adult nano characterisation and the ISH gene comparison — design

**Status: design agreed 25 Sep 2026, implementation starting.**
This document is the running record of *what we decided and why*, plus a map of the
pipeline. Update it as steps land; if the code and this file disagree, the file is wrong
and should be fixed in the same commit.

## The question

Two halves of one story:

1. **Characterise the adult distribution of surface GluA1** across the whole brain — what
   the old `merged_naive_rws` outputs did.
2. **Compare that spatial pattern against a panel of genes** and ask *which kinds of genes
   match it best*. This is also the strongest validation we have that the nanobody reports
   the SURFACE pool rather than total receptor: the nano map already correlates better
   with AMPAR trafficking and anchoring machinery (Cacng8 = TARP γ-8 at ρ = 0.78,
   Dlg2 = PSD-93 at 0.76) than with *Gria1* itself (ρ = 0.52, rank 22 of 97). See
   `project_measurement_validation_strategy` in memory.

## Shape of the thing

**Simple, like v2.** The v2 young-vs-adult route works because each step is one file, the
intermediate results sit on disk, and every step has a diagnostic. This is the same shape
applied to a different question, not a bigger machine.

It **reuses** the v2 chain rather than rebuilding it. The adult per-mouse volumes, tissue
masks, background subtraction and cohort means are already computed and audited:
`comparisons_v2/per_mouse/*.npz` (10 adults: `sig`, `auto`, `sep`, `tissue`),
`comparisons_v2/ccf/{adult,naive,rws}/` (five readings), and
`young_vs_adult/region_means_per_mouse.csv`, which already contains every adult.

Nothing in `comparisons_v2/` or `comparisons/` is written to. New outputs go to
`data/adult_v2/`. Naming is `v2_adult_*` / `v2_ish_*` deliberately: the same generation of
code answering a different question, converging on one codebase.

## The core — what gets built now

```
  [ existing, reused, not rewritten ]
  comparisons_v2/per_mouse/*.npz                10 adults: sig, auto, sep, tissue
  comparisons_v2/ccf/{adult,naive,rws}/         cohort mean / sd / n, five readings
  young_vs_adult/region_means_per_mouse.csv     per mouse x structure x reading

  [ new, three scripts ]
  v2_adult_profile.py    the adult distribution
                           . per structure: mean over the 10 adults, SEM, one dot per mouse
                           . ranked; grouped by division and macro group
                           . naive vs rws printed beside it as the built-in null
                         -> adult_v2/profile/

  v2_ish_regions.py      the Allen grids -> one region table          [the only new work]
                           . read MetaImage (67 x 41 x 58 @ 200 um, RAI)
                           . VERIFY orientation and resolution against the atlas
                           . aggregate to the same structures as the nano table
                           . per-gene coverage and quality flags
                         -> adult_v2/ish/gene_region_table.csv

  v2_ish_compare.py      the comparison
                           . join on structure; Spearman per gene
                           . run it three times: nano, sep, nano/sep
                           . rank; summarise by the existing category column
                         -> adult_v2/ish/

  diagnostics live beside each step, as in v2 -> adult_v2/processing_diagnostics/
```

The three arms are in the core because they cost nothing — the readings already exist, so
it is the same correlation run three times — and they are the whole validation argument:

| quantity | predicted best correlate |
|---|---|
| SEP (total receptor) | *Gria1* expression |
| nano (surface pool) | both, intermediate |
| nano / SEP (surface fraction) | trafficking machinery, **not** *Gria1* |

## Milestone D3 — does zref reproduce the old ranking? **Yes** (25 Sep 2026)

`v2_ish_regions.py` then `v2_ish_compare.py`, 95 genes x 280 adult structures, ten adults.
Diagnostic: `data/adult_v2/ish/ish_old_vs_new.png`, every gene's old rho against its new
one, machinery genes in red.

| v2 reading | Spearman against the old MATLAB ordering |
|---|---|
| `cref` | **+0.914** |
| `subref` | +0.913 |
| `zref` | **+0.910** |
| `ratio` | +0.869 |
| `sepratio` | +0.840 |

**Cacng8 is first under every nano reading** (zref +0.803, old +0.780), as it was in the
old route. So the affine fit and the within-brain z-score were not doing work that the
ranking depended on: one documented transform gets the same answer, and the method section
can lose a page. `zref` is adopted here on the same terms as in young-vs-adult.

Where the two routes differ, and it matters for what we claim:

| | old route | v2 `zref` |
|---|---|---|
| Gria1 rho | +0.520 | +0.665 |
| Gria1 rank | 22 / 97 | **10 / 95** |
| best machinery − Gria1 | 0.260 | 0.138 |
| machinery genes above Gria1 | 8 / 33 | 4 / 33 |

The dissociation **survives in direction but is roughly half the size**. Cacng8, Grm5,
Dlg2 and Grin2b still beat Gria1, but Gria1 is now inside the top ten rather than
twenty-second. Two readings of that, and we cannot yet choose between them: either the old
normalisation was suppressing Gria1 (its within-brain z-score flattens amplitude, which is
exactly where Gria1's cortex-heavy map lives), or the old gap was partly an artefact and
the honest effect is the smaller one. Either way the headline should be stated from the v2
numbers, not the old ones. Note also that `ratio` — the reading closest in spirit to the
old measurement — gives Gria1 +0.529, almost exactly the old +0.520.

`sepratio` (nano/SEP, the surface fraction) does **not** sharpen the dissociation as
predicted: Gria1 climbs to rank 7 and the gap narrows to 0.102. Against 95 genes with no
null, that is not evidence against the measurement model — it is a reason not to lean on
the surface-fraction arm until the third arm exists (see below).

**The third arm is missing.** All five readings have nano in the numerator, so "SEP alone
should look most like Gria1" cannot be tested yet. It needs one more reading in
`v2_cohort` — `sep / auto`, the same denominator as `ratio` — after which the three-arm
table in this document is a single figure. This is the next thing worth building.

**Panel: 95 of 100 genes.** `gene_drops.csv` lists the five and why. Sst, Chat and Tph2
return HTTP 404 from the Allen grid endpoint (retried 25 Sep, not a local problem — the
experiment ids in `gene_targets.csv` have no downloadable grid). Olig2 and Calb2 download
fine but come back in a box of their own (68 x 40 x 50 and 73 x 41 x 53 against the
reference 67 x 41 x 58) with `Offset = 0 0 0`, so there is no honest way to place them
against the atlas; cropping them to fit would assign every voxel to the wrong structure
without complaining. Both genes have other experiments in the Allen collection, so the
panel can be repaired by choosing different ids — a curation decision, deliberately not
made silently in code.

## Decisions and why

**D1 — Region-level correlation.** Voxel-level stays out of the core; it can be a
diagnostic later. Region level is what the nano table already is, and what any future null
will be defined on.

**D2 — Keep boundary protection when aggregating ISH.** ISH is 200 µm and registration is
imperfect, so a region mean picks up its neighbours. P9 eroded the masks and distance-
weighted them; keep something equivalent, having *tested* whether it matters rather than
assuming it does.

**D3 — Reproduce before improving.** First milestone is reproducing the old ranking
(Cacng8 top, Gria1 ≈ 22nd). If the new code agrees, every later difference is an
improvement rather than a bug.

**D4 — The ranking is descriptive, and says so.** No p-values in the core. See "What the
core cannot claim" below; the figures and tables must carry that caveat in writing, so
nobody quotes significance that was never computed.

**D5 — Use the 100 genes already on disk.** They are local (94 MB), coronal, verified, and
carry a `category` column. Expanding the panel is a later question, not a prerequisite.

**D6 — Existing categories now, curated ones later.** `gene_targets.csv` has
AMPAR_core / auxiliary / trafficking / scaffold / plasticity / excitatory / control_*.
Good enough to look at; not good enough to claim with.

## What the core cannot claim, and why

Correlating two brain maps over regions gives inflated significance, because everything is
high in cortex and hippocampus. Fulcher, Arnatkeviciute & Fornito (Nat Commun 2021)
measured the conventional gene-category null at **875-fold** false-positive inflation in
mouse. Our own ranking already shows the signature: interneuron markers (Htr3a +0.663,
Chrm1 +0.645) outrank genes we expected on biological grounds.

So the core produces a **ranking and a picture**, not a test. That is honest and useful —
it is what the old route produced too, only cleaner — but the claim "trafficking genes
match better than Gria1" needs the work below before it goes in a paper.

## Later, if the simple version holds up

In rough order of value, none of it a prerequisite for the core:

1. **A spatial null.** Randomise the *map*, not the gene annotations: an ensemble of
   spatially autocorrelated surrogate maps (Fulcher's "SBP-spatial"), spatial-lag model,
   mouse parameters ρ = 0.8, d₀ = 1.46 mm. Then show it is calibrated — a random phenotype
   should be significant ~5% of the time — as a figure.
   *brainSMASH* does the same thing in spirit but is built for human neuroimaging; the
   lag ensemble has published mouse parameters. Revisit if the lag model fits badly.
2. **A background panel.** ~4,000 Allen coronal genes (~1 MB each) so the curated 100
   become one category among many rather than the universe.
3. **SynGO** (Koopmans 2019) for curated synaptic terms with a proper background, with the
   sets fixed *before* looking, so "trafficking" is a prediction and not a post-hoc
   reading of a ranked list.
4. **ISH quality control as a rule.** `n_bad_sections` is recorded per gene (28 for
   Cacng8) and has never been used. Decide thresholds, apply them, report what drops.

## What the old route got wrong, and must not be inherited

| problem | consequence | fix |
|---|---|---|
| P6bis maps each mouse with an **affine** (slope *and* intercept) | ratios between regions do not survive it | use the v2 per-mouse volumes, a pure scale |
| `abs()` in P8 over zero-filled planes | unreached voxels became bright slabs | v2 tissue mask + n maps |
| tissue mask taken from the **nano** channel | circular, and a partial section reads as all-background | v2 mask (auto channel, +4 MAD) |
| group `.mat` files silently stale | nano and auto from different registrations | read the registered tiffs |
| within-brain z-score as the only reading | destroys absolute scale | the five explicit readings |

## References

PDFs in `data/ref_papers/`, summaries in `PAPER_SUMMARIES.md` there.

- **Fulcher, Arnatkeviciute & Fornito 2021**, Nat Commun 12:2669. `Fulcher_2021.pdf`.
  Why the core cannot claim significance, and what the later null should be.
- **Koopmans et al. 2019**, Neuron 103:217. `Koopmans_2019.pdf`. SynGO, for later.
- Allen Mouse Brain Atlas ISH, 200 µm grid, expression *energy*:
  <https://brain-map.org/support/documentation/api-for-mouse-brain-atlas>
