# Adult nano characterisation and the ISH gene comparison — design

This file describes the code before step 7 of the refactor (`v2_*.py` scripts, constants in the code, panels chosen by environment variable); the current names are in [refactor_name_map.csv](refactor_name_map.csv), the run order in the headers of `mapping/run_*.py`, and the constants in `mapping/settings.toml`.

**Status: designed on 25 Sep 2026 and built and run on 25 and 26 Sep 2026; a record since,
kept as the specification of A1 to A10 ([ROADMAP.md](ROADMAP.md), sections 3 and 7) until
they are done.**
This document is the record of *what we decided and why*, plus a map of the pipeline as it
was then.

## The question

Two halves of one story:

1. **Characterise the adult distribution of surface GluA1** across the whole brain — what
   the old `merged_naive_rws` outputs did.
2. **Compare that spatial pattern against a panel of genes** and ask *which kinds of genes
   match it best*. This is also the strongest validation we have that the nanobody reports
   the SURFACE pool rather than total receptor: the nano map already correlates better
   with AMPAR trafficking and anchoring machinery (Cacng8 = TARP γ-8 at ρ = 0.78,
   Dlg2 = PSD-93 at 0.76) than with *Gria1* itself (ρ = 0.52, rank 22 of 97). See
   [SCIENTIFIC_CONTEXT.md](SCIENTIFIC_CONTEXT.md), "What the map measures".

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

  v2_ish_words.py        what KIND of gene is at the top
                           . GO terms per gene from mygene.info, cached
                           . features = whole GO terms, and words from them
                           . per feature: Mann-Whitney on rho, bootstrap interval, BH q
                         -> adult_v2/ish/feature_enrichment.csv

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
| `subref` | +0.913 (+0.912 after fix 1 of step 8, which takes fibre tracts, ventricles and unassigned labels out of subref's reference) |
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

## What kind of gene is at the top — annotation, not our own labels (26 Sep 2026)

`v2_ish_words.py`. The hand-written `category` column cannot answer this honestly: we
wrote it while choosing the genes, so explaining the ranking with it is close to circular,
and it is coarse anyway ("trafficking" holds Cacng8 beside the presynaptic Bsn and Syn1).
So the grouping comes from GO terms pulled per gene from mygene.info and cached, and two
feature sets are tested against each gene's rho: whole GO terms, and single words taken
from those terms and from the gene name. Figure: `ish_word_enrichment.png`.

**The predicted contrast, named in advance** (point 3 of the validation strategy — "a
generic synaptic marker, or a SynGO presynaptic set that should NOT win, separates the two
readings"), under `zref`:

| word | genes | gap in median rho | BH q |
|---|---|---|---|
| **postsynaptic** | 68 | **+0.339** | 0.009 |
| glutamatergic | 64 | +0.358 | 0.009 |
| spine | 31 | +0.238 | 0.096 |
| dendritic | 38 | +0.219 | 0.165 |
| vesicle | 47 | +0.208 | 0.361 |
| **presynaptic** | 49 | **+0.042** | 0.723 |
| axon | 47 | +0.009 | 0.983 |
| inhibitory | 16 | −0.055 | 0.852 |
| gabaergic | 9 | −0.065 | 0.870 |

This is the useful shape of the result. Both compartments are well represented — 68
postsynaptic genes against 49 presynaptic — so the split is **not** forced by how the
panel was built, and it is a within-panel contrast, which is the only kind this design
supports. The genes that track the adult nano map are the postsynaptic glutamatergic ones;
the presynaptic and axonal vocabulary is flat. That is what the nanobody is supposed to be
reporting.

Top discovered features tell the same story from the other end: "glutamatergic synapse"
(61 genes, +0.356), "dendritic spine" (24, +0.298), "postsynaptic density", "neuron spine",
"positive regulation of ampa receptor activity".

**How to read the figure, and how not to.** Bars are ordered by effect size, which puts
five-gene groups on top for free — "asymmetric synapse", "temperature", "sleep", "early".
The 95% bootstrap whisker is there so that is visible rather than buried: those intervals
run from about −0.6 to +0.5, while `postsynaptic` sits at [+0.17, +0.54]. And the q values
are anticonservative throughout, by exactly the mechanism in Fulcher 2021 — genes sharing
a GO term are co-expressed, so their rho values are not independent. Read effect size and
interval width; treat q as a tie-breaker between features of similar size.

**D7 — annotation comes from outside, and is cached.** mygene.info, one JSON per gene in
`adult_v2/ish/annotation/`, evidence codes kept but not filtered on: dropping IEA would
strip most of the annotation off the less-studied genes and quietly bias the comparison
towards genes somebody has already studied by hand.

## Membrane pool or total receptor — and why the SEP arm cannot answer it (26 Sep 2026)

Giulio's objection, and it is correct: Gria1, Cacng8, Dlg2 and Grip1 are **all** postsynaptic,
so the postsynaptic-versus-presynaptic result above separates compartments but says nothing
about surface pool versus total receptor. No gene grouping can — SynGO included — because
both hypotheses predict a postsynaptic map. The discriminating axis is between **channels**.

`v2_adult_arms.py` puts three on one footing, per adult per structure, checked against
`v2_region_plot` on the two arms they share (max difference 5e-5, which is half the last
digit that table stores):

| arm | meaning it was supposed to have |
|---|---|
| `sepauto` = SEP / auto | total receptor |
| `ratio` = nano / auto | surface receptor |
| `sepratio` = nano / SEP | the surface fraction |

**It failed, and the failure is the result.** `SEP / auto` correlates with Gria1 expression
at **−0.11**, where plain nano manages +0.62. That sent us to the channels themselves, with
no denominators anywhere — `v2_sep_channel_check.py`, figure `sep_channel_check.png`:

| measured on structure means, 10 adults | |
|---|---|
| dynamic range (p90−p10, log2): nano / autofluo / **SEP** | 1.93 ± 0.26 / 1.07 ± 0.16 / **0.95 ± 0.16** |
| SEP ~ autofluorescence | **+0.791 ± 0.036** |
| SEP ~ nano | +0.487 ± 0.105 |
| nano ~ autofluorescence | +0.257 ± 0.110 |
| against Gria1: nano / autofluo / SEP / SEP−autofluo | +0.604 / +0.232 / +0.299 / +0.120 |

**In this fixed tissue the green channel is mostly autofluorescence.** It varies
less across the brain than the autofluorescence channel does, and it tracks it at 0.79 in
every one of the ten brains while nano tracks it at 0.26. Subtracting the autofluorescence
component linearly leaves a residual that correlates with Gria1 at only +0.12.

Three consequences, and they reach backwards through the project:

1. **`sepratio` is not a surface fraction.** It is nano divided by a second
   autofluorescence-like channel, which is exactly why it tracks `ratio` at ρ 0.89–0.97
   within every mouse and why it never sharpened the Gria1 dissociation. Wherever
   `nano/SEP` is described as "receptor on the membrane per unit receptor expressed" —
   including in `v2_cohort`'s own header and in the young-vs-adult README — that reading
   needs the caveat.
2. **The three-arm validation cannot be run on this data.** It is not a coding problem and
   no reanalysis fixes it; it needs a channel that actually reports total receptor.
3. **What survives is the partial correlation**, which never used the green channel:

**Test 2 — Gria1 partialled out of the nano map.** Even a perfect surface map correlates
with Gria1, because there is no surface receptor where there is no receptor. So: what does
the nano map still explain once Gria1 expression is removed? Over 33 auxiliary /
trafficking / scaffold genes, on ranks:

| arm | machinery median partial ρ | positive |
|---|---|---|
| `ratio` (nano / auto) | **+0.116** | 20 / 33 |
| `sepratio` | +0.138 | 21 / 33 |
| `sepauto` | −0.054 | 12 / 33 |

with Cacng8 **+0.539**, Cnih2 +0.497, Dlg2 +0.412, Grm5 +0.378 in the `ratio` arm. So the
nano map carries spatial structure that AMPAR anchoring and trafficking genes predict and
that GluA1 mRNA does not account for. That is the strongest version of the argument the
current data supports — and it is suggestive, not decisive, because mRNA is not protein:
a regional translation or turnover gradient would look the same.

**What would settle it,** in order of how decisive: a total-GluA1 antibody stain (or
autoradiography) on a subset of the same brains, giving a genuine total-receptor channel
and with it the three-arm figure; or a knockout / no-primary control, which tests specificity
rather than surface-versus-total. Both are wet-lab asks, not analysis.

## Subunit against localisation — the contrast that *does* discriminate (26 Sep 2026)

An earlier version of this document said no gene annotation could separate surface pool
from total receptor. **That was too strong and Giulio was right to push back.** What cannot
discriminate is the *compartment* axis — pre- versus postsynaptic — because both hypotheses
predict a postsynaptic map. The axis that can is **function**: what sets receptor
ABUNDANCE against what sets its LOCALISATION, and both sides of that line are postsynaptic.

| | under "nano = total receptor" | under "nano = surface receptor" |
|---|---|---|
| subunit genes Gria1–4 | carry the map | carry part of it |
| TARPs, PSD scaffolds, NSF | add nothing beyond the subunits | explain variance the subunits do not |

`v2_ish_roles.py`. The existing `category` column cannot express this — it puts the
metabotropic receptors Grm1–5 in "auxiliary" beside the TARPs and fills "trafficking"
mostly with presynaptic vesicle machinery — so the roles are re-curated by protein
function, written to `gene_roles.csv` so they can be argued with. No ρ value was consulted
in drawing them, though they were drawn after the per-gene ranking had been seen.

**Result: the direction is as predicted, and there is no evidence for it.** Both halves
matter.

| | |
|---|---|
| subunit composite alone | R² = 0.327 |
| localisation composite alone | R² = 0.427 |
| both | R² = 0.441 |
| **unique to subunits** | **+0.014** |
| **unique to localisation** | **+0.114** |
| shared | +0.312 (the two composites correlate at ρ +0.735) |

So the subunits add almost nothing once the localisation genes are in, which is what
hypothesis B says. But the two sets differ in size as well as in meaning, and a 15-gene
composite is less noisy than a 4-gene one. The test that separates those is a **within-family
permutation**: split the same 19 genes into 4 and 15 every possible way — 3,876 of them —
and recompute. The functional split lands at **p = 0.4877**, the median of that
distribution. Splitting those genes by function is no better than splitting them at random.

Sensitivity, because the objection is obvious: the localisation set carries passengers —
Cacng5, Cacng7 and Grip2 are close to absent from the forebrain. Filtering the family on
**expression only** (median energy at or above the family median, applied to both sides,
never looking at a ρ) leaves 10 genes, and sharpens the asymmetry — unique to subunits
+0.000, unique to localisation +0.131 — while the permutation still gives **p = 0.3128**.

**Specificity control**: presynaptic vesicle machinery, which is membrane trafficking of a
completely different kind, correlates at median +0.363 against localisation's +0.385. If
the map preferred AMPAR localisation specifically, it should not.

**And note what this does to the earlier headline.** Cacng8 +0.80 against Gria1 +0.66 is a
single-gene comparison that does *not* generalise to its category: the subunit median
(+0.582, with Gria1 +0.66, Gria3 +0.65, Gria2 +0.51, Gria4 −0.02) is *above* the
localisation median (+0.385, spanning Cacng8 +0.80 down to Cacng7 −0.26). Quote Cacng8 as
one gene, never as evidence about trafficking genes in general.

**This is "no evidence", not "evidence of no effect", and the test is weak.** There are
only four AMPAR subunits and one of them is anti-correlated, so the surrogate distribution
is wide: its p95 is +0.44 against an observed +0.10, meaning only an effect roughly four
times larger could have been detected. What would give it power: a much larger panel so
the composites are built from tens of genes rather than four; a finer parcellation or
voxel-level maps for more independent observations; and — still the only thing that settles
the question rather than sharpening it — a real total-receptor channel.

## The powered version of the test — built, run, and negative (26 Sep 2026)

The previous section ended with "no evidence, and the test is weak". This builds the panel
that could carry it and runs it properly. Scripts, in order: `v2_panel_build.py`,
`v2_panel_fetch.py`, `v2_ish_regions.py` (same code, new panel via `V2_ISH_PANEL`),
`v2_ish_reliability.py`, `v2_ish_panel_test.py`.

### The panel comes from the ontology, not from us

Gene Ontology terms queried through mygene.info and cached; every gene in
`panel_genes.csv` carries the term ids that placed it, so the panel is arguable term by
term rather than gene by gene.

| set | terms | genes |
|---|---|---|
| subunit | GO:0004971 ∩ GO:0032281 | 4 |
| localisation | GO:0099072, GO:0099645, GO:0097113, GO:0098970, plus GO:0032281 minus the subunits | 84 |
| control_psd | GO:0014069 postsynaptic density, minus the above | 300 |

**One hand-made decision in the whole panel, and it is written into the code as
`OVERRIDE`:** GO annotates the delta receptors Grid1 and Grid2 to both AMPA terms. They
are a different receptor family and their transcript level does not set how much AMPA
receptor a region has, so they get their own role instead of counting as subunits. They
stay in the panel so the test can be run either way.

**Both planes of section**, where the old panel was coronal only: 390 genes have at least
one experiment, 669 experiments in all, 667 downloaded (two pre-2005 experiments ship a zip
with no energy grid; both genes have others). 479 MB.

### One Allen map is more reliable than we assumed — including Gria1's

218 genes have more than one experiment, so their map can be correlated against itself:

| | |
|---|---|
| reliability, median over 218 genes | **+0.691**, quartiles [+0.519, +0.770] |
| coronal against sagittal (n = 261 pairs) | +0.689 |
| same plane (n = 80 pairs) | +0.567 |
| genes below 0.3 | 22 |

**Gria1's own map scores 0.91**, Gria2 0.90, Gria3 0.81, Gria4 0.57; Cacng8 0.91, Dlg2
0.91, Dlg4 0.76, Nsf 0.66. That closes a caveat open since April: Gria1 ranking below the
trafficking genes is **not** an artefact of a bad Gria1 experiment. Replicates are merged
on ranks (expression energy carries an arbitrary per-experiment scale) into
`gene_region_table_merged.csv`, which every test below uses.

### The statistic had to change too

Comparing a 4-gene composite with a 15-gene one confounds meaning with noise, and no panel
fixes it because the mouse has four AMPA subunit genes. So the question is asked per gene:

> for every non-subunit gene, how much of the map does it explain **once the subunit
> composite is removed** — `partial ρ(map, gene | subunits)` — and do localisation genes
> retain more than other postsynaptic genes?

Both sides are now large, and the subunit composite enters once as a covariate where its
noise costs both sides equally. Controls are matched one-to-one to localisation genes on
median expression energy, because reliability rises steeply with expression and an
unmatched comparison would partly measure which set contains louder genes.

### Result: no difference, and this time that means something

| contrast (zref) | localisation | control | difference | p |
|---|---|---|---|---|
| partial ρ, expression-matched | +0.139 (84) | +0.149 (84) | **−0.010** | 0.74 (0.75 after fix 5 of step 8: each permutation test gets its own generator, a move within Monte Carlo error) |
| partial ρ, all controls | +0.139 (84) | +0.123 (300) | +0.016 | 0.65 |
| plain ρ, before partialling | +0.390 | +0.379 | +0.011 | 0.89 |
| partial ρ, reliability ≥ 0.3 only | 45 genes | 39 genes | −0.004 | 0.91 |

**Sensitivity, so the null is interpretable**: a difference of **±0.062** would have been
detected at p < 0.05; the observed is 0.16 of that. **Positive control**: the same
machinery, on control genes split at their median reliability, finds |ρ| 0.509 against
0.354 — difference **+0.155, p = 0.0007** (0.0006 after fix 5). So the test detects a real effect of that size
with these sample sizes and does not detect this one.

**AMPAR localisation genes, as a class, explain no more of the adult nano map than
expression-matched postsynaptic genes do, once receptor abundance is taken out.**

### What the map actually tracks

The top of the control list answers that, and it is not a list of hidden TARPs: Arpc5
+0.602, Cdk5r1 +0.560, Ptk2b +0.552, Slc8a2 +0.547, Rgs14 +0.513, Cap2 +0.512, Baiap2
+0.501, Grin2a +0.499. General postsynaptic signalling and cytoskeletal genes of forebrain
excitatory neurons. Cacng8 is still first overall at partial +0.695, with Dlg2 +0.546 and
Igsf11 +0.558 — but Arpc5 and Cdk5r1 sit in the same band, and they have nothing to do with
AMPA receptor trafficking.

So the honest reading of the whole ISH arm: **the adult nano map is predicted about equally
well by any well-measured forebrain postsynaptic gene.** Cacng8 topping the ranking is real
and reproducible, and it is not evidence that the map is about trafficking, because the
class it belongs to carries no more signal than the neighbours.

### What is left

Three things, none of them another gene panel:

1. **A real total-receptor channel** — a total-GluA1 antibody stain or autoradiography on a
   subset of the same brains. Still the only thing that settles surface-versus-total rather
   than sharpening it.
2. **A spatial null** (Fulcher SBP-spatial), which would tell us how much of the +0.39
   baseline correlation shared by *every* postsynaptic gene is simply "everything is high in
   cortex and hippocampus".
3. Accepting the descriptive result as descriptive, and not building the paper's
   measurement-validation argument on it.

## Is the map just abundance, or just synaptic density? **Neither** (26 Sep 2026)

The two things a reader will say the nano map is. Ranking genes cannot answer either — and
the flat ranking in the previous section is exactly what you would expect if the map were
only tracking how much synapse a region has. So the question is turned around: **take both
explanations out, and ask whether what is left is signal or noise.**

    nano zref  ~  receptor abundance + synaptic density   ->   residual
    does the residual come out the same in two independent halves of the cohort?

A residual is never zero; what makes one interesting is reproducibility, because noise
cannot replicate across independent animals. Two scripts: `v2_beyond_density.py` does the
analysis in four numbered steps with a figure each, and **`v2_beyond_controls.py` exists to
break it** — seven named ways the claim could be wrong, each with the number that decides.

### Step 0 — which structures, decided by rule before any fitting

Grey matter only (by division), and no `..., unassigned` entries — those are voxels the
atlas could not place, not structures, and what lands in them depends on how each brain
registered. This matters: in the first run "brain, unassigned" (+57 ranks) and "cerebrum
related" (+49) sat near the top of the residual. **125 structures** survive the rule;
`structures_used.csv` lists every one kept and dropped with the reason, and `fig0` draws it.

### Step 1 — the ceiling

Over all 126 five-against-five splits of the ten adults, the half-cohort maps agree at
**ρ = 0.974**; Spearman-Brown gives **0.987** for the full cohort. So **97.4% of this map's
variance is explainable in principle**. Every R² below is quoted against that, not against
1.0.

### Step 2 — what the explanations buy

Scored by cross-validation on held-out structures, because a model with more terms always
fits better on the data it was fitted to:

| model | CV R² | of the ceiling |
|---|---|---|
| abundance (Gria1–4) | 0.253 | 26% |
| synaptic markers (Syp, Syn1, Vamp2, Bsn, Syt1, Dlg4, Homer1, Shank2/3, Nlgn1, Camk2a) | 0.148 | 15% |
| psd_pc1 (first PC of 188 postsynaptic-density genes) | 0.262 | 27% |
| autofluorescence, same brains | −0.043 | 0% |
| all four, straight lines | 0.416 | 43% |
| **all four, allowed to bend** | **0.597** | **61%** |

The last row is the one to quote — it is the covariates' best shot. **It leaves 39% of the
explainable variance unaccounted for.**

### Step 3 — and the leftover replicates

Residualise each half-cohort map on that bending model and compare the two leftovers:
**ρ = 0.934**, Spearman-Brown **0.966**, against the map's own 0.974 / 0.987. The leftover
is very nearly as reproducible as the map itself.

### Step 4 — where it lives

Medial geniculate +48 ranks, subthalamic nucleus +44, ventral LGN +38, lateral habenula
+38, septofimbrial +35, VISal +32, CA3 +31, subiculum +30. Against VPM −49, VPL −47, RSPd
−46, posterior complex −43. No single gene of the 390 accounts for it (best: Cacng8 +0.301).

### The controls — and the one that failed

| | control | result |
|---|---|---|
| A | smooth spatial gradient (an illumination artefact; the sections are not cleared) | quadratic in AP/DV/ML explains R² 0.18; with position as a covariate the leftover still replicates at 0.934 — **pass** |
| B | small or noisy structures | ρ with log volume −0.116 — **pass** |
| C | one or two odd animals | every pair of mice agrees, median 0.780, worst 0.595 — **pass** |
| D | the whisker manipulation | naive leftover vs RWS leftover ρ +0.884 — **pass** |
| E | curvature modelled away as a straight line | **FAILED FIRST TIME** — see below |
| F | our choice of covariates | 20 components of the whole 390-gene space reach CV R² 0.826 (85% of the ceiling), and their leftover *still* replicates at 0.879 — **pass** |
| G | specific to zref | every reading gives the same picture; leftover replicates 0.912–0.936 — **pass** |

**Control E failed, and the analysis changed rather than the wording.** Fitted as straight
lines, the covariates reached CV R² 0.416; allowed squares and cubes, 0.597. A fifth of the
map was being credited to the leftover that the covariates could in fact explain. So the
headline model now bends, the claim moved from "half the explainable variance" to **39%**,
and E now asks the follow-up instead: going further, to fifth powers, buys 0.602 against
0.597 — the model is bent enough.

That is the whole reason for the controls, and it is worth saying plainly: the first
version of this result was too generous to us by about ten points.

### What to claim, and what not to

**Claim:** *about 40% of the explainable variance of the adult surface-GluA1 map is not
accounted for by receptor abundance or synaptic density, even when those are given a
flexible relationship and the whole cohort's own autofluorescence, and what remains
replicates across independent animals at 0.97. A twenty-component model of the entire
390-gene panel — far richer than either explanation — still leaves a remainder that
replicates at 0.88.*

**Do not claim** that the leftover *is* the surface fraction. A residual is only ever "not
explained by what we put in": mRNA is an imperfect proxy for protein, and anything
spatially organised that was left out lands in it — regional differences in translation,
turnover, subunit composition, or nanobody access to tissue. Separating those is what a
real total-receptor channel is for.

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
