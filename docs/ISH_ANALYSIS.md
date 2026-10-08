# The ISH analysis

What the adult nano map is, read against the Allen in situ hybridisation (ISH)
maps of the adult mouse brain: the question, the argument, each analysis with its
method, figure, result and meaning, the limits, and how to rerun it. Built on 8
October 2026 from the decisions of the ISH discussion (all five as recommended,
`G:\sep_refactor\ish_discussion\ISH_DISCUSSION.md`, section 8).

- Code: `mapping/sepmap/structures.py`, `mapping/sepmap/adult/` and
  `mapping/sepmap/ish/`, run by steps 13 to 27 of the Python route
  ([mapping/README.md](../mapping/README.md)).
- Outputs: `<data>\adult_v2\ish_analysis\`. The numbers below are those of the
  run of 8 October 2026 on a full copy of the production inputs
  (`SEP_DATA_ROOT=G:\sep_refactor\check\data`); the production data root holds
  the same once steps 13 to 27 run there. Every number below is in
  `tables\numbers_for_the_text.csv` there (a `.txt` beside it reads more
  easily), written by `run_ish_overview.py` from the numbers each step writes;
  where a number comes from another table, the table is named.
- Figures: `figures\00_overview.png` to `14_april_headline.png`, PNG and EPS.
  `figures\README.md` walks through them in order, each with its question,
  what to look at, what to take from it and what it means, with the numbers of
  the run. The
  figures are not versioned, so this document names each by its file.

## 1. The question and the argument

The adult nano map orders the brain's structures, and two halves of the cohort
order them the same way (rho 0.974 between half-cohort maps, 126 grey-matter
structures). A reader asks first whether that order is simply how much GluA1
mRNA a structure makes, or how many synapses it has. Allen's ISH maps answer
that, because they measure both: the receptor subunits' mRNA, and the mRNA of
synaptic and postsynaptic-density genes.

The ISH line serves one argument, in two connected parts:

1. **The map is not satisfactorily explained by Gria1 expression or by synapse
   density.** Predicted from receptor mRNA, synaptic markers and
   postsynaptic-density genes, a sizeable part of the map is left over; it
   replicates across mice and stands above the floor that Allen-to-Allen
   mismatch alone would leave. This is the core claim (section 4).
2. **The map is therefore something else, and the reading the data support is
   the surface fraction of the receptor**, shaped by trafficking regulation and
   scaffolding. The gene analyses ask whether they corroborate that reading:
   the genes that regulate surface AMPA receptors (Cacng8, a TARP;
   trafficking and scaffolding genes) should track the map better than
   abundance genes or unrelated genes, beyond a null that respects the brain's
   smoothness (section 5).

The limit, once: the surface fraction is the interpretation the data support,
not a measurement. A total-GluA1 stain on the same brains would measure it; the
green SEP channel cannot (section 6.3).

What the comparison cannot do on its own: separate surface from total receptor
(both follow the postsynaptic side of synapses), read protein from mRNA (mRNA
sits in cell bodies, receptor on dendrites), or give a p value without a null
that keeps the brain's smooth gradients (cortex and hippocampus high, thalamus
and hypothalamus low). Each Allen experiment is one P56 mouse on a 200 um grid.

`figures\00_overview.png`: the question, the argument, one row per step with
the numbers of the run and what stays open.

## 2. Words used here

| word | meaning |
|---|---|
| structure | one atlas region (CA1, VPM). Every comparison runs across structures, one value per structure on each side |
| rho | Spearman correlation across structures: do two maps put the structures in the same order |
| division | a coarse part of the brain: Isocortex, OLF, HPF, CTXsp, STR, PAL, TH, HY, MB, P, MY, CB |
| `zref` | per brain, log2 of each structure's nano mean over the brain's isocortex mean, minus its median over the declared set, over its p90 - p10 spread there; the adult map is the mean of the ten adults' `zref` |
| the declared set | the structures every comparison uses: grey matter measured in all ten adults (section 3) |
| spatial null | random maps with the nano map's smoothness, the surrogates; a spatial p is the share of surrogates that correlate with a gene at least as strongly as the map does |
| partial rho | rho once a third map (the subunit composite) has been regressed out of both sides, on ranks |
| ceiling | the share of the map that two halves of the cohort reproduce (Spearman-Brown of their agreement) |
| leftover | the part of the map a fitted model does not predict |
| calibration floor | what the same model leaves of a map made only of receptor mRNA and density: the leftover Allen-to-Allen mismatch produces alone |
| P9 | the MATLAB comparison of April 2026, `adult_matlab/run_compare_with_allen_ish.m` |

## 3. The inputs

Every analysis reads the same inputs, fixed before any correlation.

**The structures** (A1; `run_structure_set.py`, `figures\01_structures.png`).
A structure enters when it is grey matter and every adult measures it: of the
280 structures of the adult table, 235 are measured in all 10 adults and 204 of
those are grey matter (TH 45, HY 43, Isocortex 39, MB 25, STR 14, HPF 11, OLF
9, PAL 9, CTXsp 7, P 2). The structures seen in one to four adults, mostly pons
and medulla, are left out: their values do not replicate across mice.

**The adult map.** The mean of the ten adults' `zref`, with its zero and spread
taken over the declared set, so they depend only on the brain and the list,
not on which other brains are in a run. Against the 17-brain reference used
before, each adult's zero moves by 0.012 to 0.103, and the cohort map's order
does not change (rho 0.99995). The autofluorescence map of the same sections is
read the same way.

**The genes** (A9; `run_ish_gene_table.py`, `figures\02_genes.png`). One table,
`tables\gene_table.csv`, holds P9's 100 genes and the 390 genes of the ontology
panel, 451 genes, with a row per Allen experiment: 757 experiments with the
repair (other experiments of Chat, Tph2, Olig2 and Calb2), 740 of them usable,
each excluded one with its reason. All 451 genes have a profile, P9's 100
included. Each gene carries its section QC, its reliability, its gene sets
(from GO, release 2026-07-26) and P9's category as a label only.
`tables\gene_documentation.csv` is the table for a supplement.

**Section QC** (A2; `run_ish_section_qc.py`, `figures\qc\`). Along each
experiment's own section axis, a section is set missing when its median energy
is below 0.2 of the median of its three neighbours on each side; it is never
filled in. Of 25,687 sections judged, 196 were set missing, in 133 experiments
(12 of them P9's own); 5 dim sections are kept as true absence, from
`mapping/ish_section_exceptions.csv` (anterior Tac1, posterior Glra1),
proposed and not yet reviewed (section 8). `figures\qc\00_flagged.png` lists
every flag; one sheet per experiment shows its profile and orientation.

**Reliability.** For the 222 genes measured more than once, two Allen
experiments of the same gene agree at a median rho of 0.689 (quartiles 0.516 to
0.772; 23 genes below 0.3). Gria1, Cacng8 and Dlg2 are at 0.91. A gene's rho
with the nano map is capped by its own reliability, so a low rho of an
unreliable gene says little.

**The gene sets** (fixed in `mapping/sepmap/ish/gene_sets.py` and committed
before any rho on these inputs): subunits (4 genes), localisation (84: the GO
terms of AMPA receptor transport, anchoring, clustering and the auxiliary
subunits), other postsynaptic (196), presynaptic (16), GABAergic neuron markers
(5) and glia (7, P9's list). The presynaptic set is small because GO annotates
the vesicle genes to both sides of the synapse, and the rule keeps a gene in
one side only (section 8).

## 4. Part 1: not explained by receptor mRNA or synaptic density

`run_beyond_density.py` to `run_beyond_figures.py`; `figures\03_beyond_budget.png`
and `04_beyond_where.png`.

**Method.**

- The structures: the declared set, less those where a subunit or marker gene
  has no ISH value: 126 structures, 10 adults.
- The map and every predictor as ranks across structures. Predictors:
  abundance, Gria1 to Gria4 as four separate terms; density, 11 synaptic marker
  genes and the first principal component of 186 postsynaptic-density genes;
  the autofluorescence of the same brains. Each enters as x, x² and x³.
- Each model is scored on structures it has not seen (cross-validated R²),
  as a share of the ceiling: two halves of the cohort, in all 126 splits,
  agree at 0.974, and Spearman-Brown makes that 0.987 for ten adults (97.9% to
  99.0% over resampled structures). Left = 1 - share.
- Replication: the model fitted to each half-cohort map, and the two leftovers
  correlated. The map's reliability and the fit alone imply a value; the
  observed one is set beside it.
- Calibration: the same model on maps whose answer is known, built from half of
  each gene's Allen experiments and predicted from the other half, with ten
  made-up adults as noisy as ours: a map made only of receptor mRNA and density
  (the floor), and a map that is one Allen Gria1 experiment. 113 structures,
  where both halves measure every subunit and marker.
- Seven controls try to break the result: a spatial gradient, structure size,
  single animals, naive against RWS, curvature, the whole gene table, the
  reading.

**Result** (`numbers_for_the_text.csv`, step `beyond`; `beyond\variance_partition.csv`,
`calibration.csv`, `controls.csv`):

| predictors, each bent, held out | share of the reproducible map |
|---|---|
| Gria1 alone | 47% |
| Gria1 to Gria4, four terms | 64% |
| density alone (markers, postsynaptic-density component) | 58% |
| autofluorescence alone | -6% |
| the model: Gria1 to Gria4, then density (+8%), then autofluorescence (+0.5%) | 73% |
| left | 27% (17% to 38% over resampled structures) |

- The leftover replicates: the two half-cohort leftovers agree at 0.928 (lowest
  split 0.873). The map's reliability and the fit alone imply 0.883; two
  unrelated leftovers fall between -0.17 and +0.17.
- The calibration floor: a map made only of receptor mRNA and density leaves 15%
  (12% to 19% over the draws), its leftover replicating at 0.80. The nano map,
  predicted the same way on the same 113 structures, leaves 29% to 31%. A map
  that is one Allen Gria1 experiment leaves 34% (28% to 41%).
- All 7 controls pass. Control F, 25 components of 311 genes with no
  hypothesis, predicts 84%; its leftover still replicates at 0.895.
- The April model, with the four subunits averaged into one term, leaves 36%
  (30% to 53%): Gria4 runs against the map (rho -0.12, `tables\gene_ranking.csv`)
  and dilutes the average.
- Where it sits (`beyond\regression_table.csv`): above prediction in the medial
  geniculate (+50 ranks), the rostrolateral visual area, the septofimbrial
  nucleus and the subiculum; below it in the substantia innominata, piriform
  cortex, ventral retrosplenial cortex and VAL.

**What it means.** About a quarter of the map's reproducible pattern, 27%,
is not predicted by receptor mRNA or synaptic density, and it is reproducible
across mice. It is about twice the floor, so it is more than what receptor
mRNA and density would leave through Allen-to-Allen mismatch; its range over
resampled structures reaches down to 17%, close to the floor's highest draw
(19%).

**What it does not mean.** It is no larger than what a map of one Allen Gria1
experiment leaves (34%), so a pessimistic benchmark, in which the truth itself
carries one experiment's error, could produce it. Its replication follows from
the map's reliability (0.883 expected), so it is not separate evidence. It is
"not predicted by receptor mRNA or synaptic density", not "beyond gene
expression": the whole gene table leaves 16%. And a leftover says what the
predictors miss, not what it is; a claim about one structure needs its own null.

## 5. Part 2: what the leftover looks like, through the genes

If the map is the surface fraction, genes that set how much receptor reaches
and stays at the membrane should follow it better than abundance genes or
unrelated genes. Each test below asks that, against the spatial null.

### 5.1 One comparison

`figures\05_one_comparison.png`. A gene's map and the nano map are each reduced
to one value per declared structure (for a gene, the mean rank of its usable
experiments), and the two orders correlated. Cacng8 gives +0.807 on 172
structures, Gria1 +0.646 on 164, the astrocyte gene Aqp4 +0.021. Much of a
whole-brain rho is cortex and hippocampus against thalamus, which is why the
null and the test inside divisions follow.

### 5.2 The spatial null

`run_ish_spatial_null.py`, `mapping/sepmap/ish/spatial_null.py`;
`figures\06_spatial_null.png`.

**Why.** Brain maps share smooth gradients, so two unrelated maps correlate by
chance far more than shuffled ones: in the calibration, unrelated smooth maps
correlate with an SD of rho 0.213, a shuffled map with 0.071, and the ordinary
Spearman p falls below 0.05 for 62% of random maps tested against nano.

**How** (Burt et al. 2020, written in numpy and scipy; brainsmash was not
installed, because it imports scikit-learn and joblib, which `venv_atlas` does
not have). Distances between structures are those between their centroids in
one hemisphere, in mm. A surrogate is the map's ranks shuffled, smoothed over
each structure's nearest neighbours with the share of neighbours whose
variogram fits the map's best (up to the 25th percentile of the distances), and
given the map's own values in its new order. 10,000 surrogates per map, nano
and autofluorescence, on the 204 declared structures, cut to each gene's
structures. A spatial p is two-sided; Benjamini-Hochberg runs within P9's genes
and within all genes.

**Checks.** Over the matched range the surrogates' variogram is within a median
6.1% of the map's (23.7% at most, at the shortest distances). The calibration
uses random fields of another kind (Gaussian, exponential covariance with a
4.0 mm range fitted to the map), so the null is not checked against its own
generator: random maps tested against nano get a spatial p below 0.05 in 4.0%
of 2,000 tests (5% expected), pairs of random maps in 3.2%.

**Limits.** Centroids stand for structures of very different sizes; the
variogram is matched at short range; and since smooth maps vary together over
large parts of the brain, a difference between two genes' rho needs to be
large to pass.

### 5.3 Gene by gene (analysis 1)

`run_ish_gene_ranking.py`; `figures\07_gene_ranking.png`.

**Result** (step `gene_ranking`; `tables\gene_ranking.csv`, `gap.csv`):

- 12 of 451 genes pass the spatial null at BH q < 0.05, 6 of P9's 100: Cacng8,
  Arpc5, Igsf11, Htr3a, Grm5, Neurl1a, Mapk1, Add3, Gria1, Cnih2, Grip1, Eps8.
  Four are localisation genes (Cacng8, Igsf11, Cnih2, Eps8), one is a subunit.
- Cacng8 is first, +0.807 (p ≤ 0.0001), and first in every robustness variant
  (section 6.2). Gria1 is +0.646 (p ≤ 0.0001, q 0.0025 within P9's genes), 11th
  of P9's genes and 32nd of all.
- The Cacng8 - Gria1 gap, on the 164 structures both have: +0.167, steady
  across adults (+0.156 to +0.179 over resampled adults), but inside its null:
  spatial p 0.348. It moves with the Allen experiment that stands for each gene,
  from +0.08 to +0.211 over the four pairings.

**What it means.** The map follows a TARP's pattern and Gria1's beyond a map
with the brain's smoothness. Cacng8 leads Gria1 in each of the 10 adults and
in every variant, but by a margin smooth maps reach by chance: these maps cannot tell
which of the two it follows more closely. A description of the map, and
consistent with the surface-fraction reading, not evidence for it.

### 5.4 Kinds of genes (analysis 3)

`run_ish_gene_sets.py`; `figures\08_gene_sets.png` and `09_localisation.png`.

**Method.** Each set's median rho against the median of the same genes' rho
with every surrogate, so co-expressed genes stay together in the null; BH over
the sets with at least 5 genes. Two contrasts named in advance: the three
postsynaptic sets pooled against the presynaptic set, and against glia. Then
the localisation test of 5 October on the new inputs: each gene's partial rho
with the map once the subunit composite is removed from both; the 84
localisation genes against 84 other postsynaptic genes matched on expression,
labels permuted; a positive control (control genes with a reproducible map
against those without) that must come out for the test to say anything; and
the matched difference against the surrogates.

**Result** (step `gene_sets`; `tables\set_tests.csv`, `contrasts.csv`,
`localisation_summary.csv`):

| set | genes | median rho | spatial p | q |
|---|---|---|---|---|
| subunits | 4 | +0.57 | too few to test | |
| localisation | 84 | +0.375 | 0.054 | 0.14 |
| other postsynaptic | 196 | +0.398 | 0.055 | 0.14 |
| presynaptic | 16 | +0.124 | 0.26 | 0.43 |
| GABAergic markers | 5 | +0.022 | 0.87 | 0.87 |
| glia | 7 | +0.063 | 0.57 | 0.72 |

- No set passes at BH q < 0.05.
- Postsynaptic above presynaptic: +0.271, spatial p 0.0174. Postsynaptic above
  glia: +0.333, spatial p 0.173 (7 glial genes).
- Localisation against matched controls, abundance removed: +0.0038, p 0.906
  (spatial p 0.887); a difference of 0.077 would have been detected. But the
  positive control does not come out with the GO control pool (+0.043, p
  0.397), so by its own rule this test says nothing. With the control pool of 5
  October it does (+0.111, p 0.0098), and localisation is still no better than
  its controls (+0.024, p 0.604). The localisation genes with the highest
  partial rho are Cacng8, Igsf11, Dlg2, Dlg3 and Lgi1
  (`localisation_test.csv`).

**What it means.** The map is postsynaptic-like, as any glutamate receptor
label should be (a sanity check that passes). The genes that set surface AMPA
receptors do not stand out from other postsynaptic genes of the same
expression, beyond the null.

### 5.5 Inside divisions (analysis 2)

`run_ish_divisions.py`; `figures\10_between_within.png` and the gene sheets
`figures\genes\<gene>.png` (Cacng8, Gria1, Grm5, Dlg2, Aqp4).

**Method.** Per gene, rho with a map that knows only each structure's division
(its median), and the mean rho inside the divisions with at least 8 structures
(HPF, HY, Isocortex, MB, OLF, PAL, STR, TH), weighted by their structures, so
no contrast between divisions can enter it. Tested against the surrogates of
the whole map; a shuffle of structures inside divisions is shown beside it.

**Result** (step `divisions`; `tables\within_division.csv`):

- The division-only map orders the genes as the real map does (0.967). Inside
  divisions the agreement falls to 0.574, and the median rho from +0.346 to
  +0.103.
- 11 genes pass the within null (7 of P9's): Cacng8 +0.548, Dlg2 +0.513,
  Neurl1a, Htr3a, Gria1 +0.417, Nptx1, Arrb2, Arpc5, Slc8a2, Tfrc, Rap2a.
- The shuffle inside divisions would have passed 32.2% of random smooth maps
  (the surrogates 2.6%), so it is shown and not used.

**What it means.** Cacng8, Dlg2 and Gria1 follow the map at a fine scale, not
only through the cortex-against-thalamus contrast: the stronger claim, and it
holds for both a TARP and the receptor's own mRNA.

### 5.6 The leftover itself

`figures\04_beyond_where.png`, panels C and D; `beyond\leftover_genes.csv`,
`leftover_sets.csv`. Every gene against the leftover of part 1, each against
the leftover's own surrogates: 0 of 451 pass at BH q < 0.05. The closest is
Aqp4 (+0.331, p 0.030 before correction); Cacng8 is +0.150 (p 0.57, 41st). The
localisation set's median is +0.033 (p 0.87); glia +0.188 (p 0.043 before
correction). None of these tests was named before the leftover was seen.

### 5.7 What part 2 says

Consistent with the surface-fraction reading, and not singling it out. The map
follows receptor expression and not the tissue (section 6.1), it is
postsynaptic-like, and Cacng8, Dlg2 and Gria1 lead the ranking inside
divisions too. But Cacng8's lead over Gria1 is inside its null, no gene set
passes its null, the localisation genes predict the map no better than matched
postsynaptic controls, and no gene follows the leftover. The genes neither
contradict nor confirm part 2; only a total-receptor measurement can.

## 6. Controls and limits

### 6.1 The tissue: autofluorescence (A8)

`figures\11_autofluorescence.png`. The autofluorescence map of the same
sections, read as nano is and tested with its own surrogates.

- With Gria1 +0.176 (p 0.447), with Cacng8 +0.296 (p 0.292): neither passes.
  Adult by adult, nano's rho with Gria1 is 0.482 to 0.754, autofluorescence's
  -0.266 to +0.297 (`tables\gene_ranking_per_adult.csv`).
- The tissue has a pattern of its own: 38 genes pass the autofluorescence null
  (7 of P9's), ribosomal proteins, Mapt, Eno2 and GABA-A subunits among them,
  and the nano and autofluorescence gene orders agree at only 0.37.

So the nano ranking is the label's: it follows Gria1 in every adult, and the
autofluorescence ranking is a different one.

### 6.2 The choices made: robustness (A3)

`run_ish_robustness.py`; `figures\12_robustness.png`; `tables\robustness_summary.csv`.
Eleven variants of the primary ranking (Spearman, `zref` with the declared
reference, merged profiles after QC, full means, the declared set): ten change
one choice, the last is the route of 5 October as it ran.

| variant | gene order against the primary (P9's genes) | Gria1's rank | gap |
|---|---|---|---|
| Pearson on log2 | 0.960 | 15 | +0.191 |
| ISH means eroded | 0.995 | 6 | +0.151 |
| both sides eroded | 0.994 | 5 | +0.144 |
| `ratio` reading | 0.966 | 8 | +0.132 |
| P9's single experiment per gene | 0.986 | 17 | +0.211 |
| P9's nine divisions, any number of adults | 0.993 | 15 | +0.180 |
| every structure of the adult table | 0.941 | 5 | +0.126 |
| the route of 5 October | 0.938 | 9 | +0.144 |

The nano erosion, the stored `zref` and no section QC change nothing (0.9996
and above). Cacng8 is first in every variant. Gria1's rank is what moves (5th
to 17th): quote it with its null, not as a place.

### 6.3 The limit: the green channel (analysis 5)

`run_sep_channel_check.py`; `figures\13_green_channel.png`;
`green_channel\sep_channel_check.csv`.

The tag's own green (SEP) fluorescence was meant to show all tagged receptor,
so that nano over SEP would be the surface fraction. In every adult the green
channel follows autofluorescence across the 204 structures (rho 0.723 to
0.827), and varies across them about as little (p90 - p10 of log2: SEP 0.915,
autofluorescence 0.999, nano 1.967). With Gria1: nano +0.619, SEP +0.333,
autofluorescence +0.113. What is left of SEP once its autofluorescence part is
regressed out follows nano at +0.579 and Gria1 at +0.343: tag that survived, or
nano's fluorescence leaking into the green channel; the filter sets decide which.

So the surface fraction cannot be measured in these brains. A total-GluA1
stain (permeabilised) or autoradiography on some of the same brains would
measure it, and would separate it from translation, turnover, subunit
composition and nanobody access. The tests of the channel ratios against genes
(`ish/arms.py`, `arms_vs_genes.png`) are retired: their premise failed.

### 6.4 Limits

- Each Allen experiment is one P56 mouse on a 200 um grid; the adults' age is
  not recorded.
- mRNA sits in cell bodies, receptor on dendrites, which can lie in another
  structure (cortical layer 1).
- The gene panel is 451 genes chosen from GO and by hand, not a genome-wide
  background.
- The spatial null rests on structure centroids and a short-range variogram.
- The exceptions list (true absences) is a human call, still proposed.
- No total-receptor channel, no knockout or no-primary control, no DAPI
  control through the chain.

## 7. What changed from April

`run_ish_overview.py`; `figures\14_april_headline.png`; `tables\april_anova.csv`,
`april_groups.csv`.

P9's headline (22 April) grouped 97 genes by the categories of
`gene_targets.csv`, written while the genes were chosen, with "auxiliary"
split by hand into Aux forebrain (Cacng8, Cacng3, Cnih2, Cnih3, Grm5) and Aux
other, and tested the ten groups with a one-way ANOVA: p 0.032.

- The gene order reproduces: today's rho against April's agrees at 0.96 over 97
  genes.
- The group p rests on the split: 0.20 without it. Today it is 0.011 with the
  split and 0.148 without; 0.002 to 0.033 over structure sets and P9's single
  experiments.
- Every such p treats co-expressed genes as independent draws. Against the
  map's surrogates, which keep them together, the groups' F gives p 0.286. One
  of the nine groups with enough genes passes its own null at BH q < 0.05: Aux
  forebrain, the group written after looking.

So the April result that stands is the gene order, Cacng8 first; the
difference between kinds of genes does not survive the null, and section 5.4
asks it again on sets fixed in advance.

| | P9 (April) | Python, ISH line (October) |
|---|---|---|
| nano map | voxel mean, each adult fitted onto its group | per-adult structure means, `zref` over the declared set, mean of 10 adults |
| ISH onto the atlas | grid stretched; -1 clamped to 0 | atlas sampled onto the grid; -1 missing |
| sections | dim sections interpolated | failed sections set missing, true absences on a reviewed list |
| structures | 206 of nine divisions | 204 grey-matter structures measured in all 10 adults |
| genes | 100 hand-picked, one experiment each | 451 (P9's and the ontology panel), every usable experiment merged |
| statistic | eroded Spearman | Spearman on full means; Pearson and eroded means as robustness rows |
| groups | hand categories and a hand split | GO sets and cited marker lists, committed first |
| tests | one ANOVA over genes | a spatial null for every rho, the gap, the sets and the contrasts |
| controls | none | autofluorescence through the same chain, robustness, positive controls |

## 8. Points to settle

- **The exceptions list.** Tac1 sections 4 and 5 and Glra1 15, 16 and 61 are
  kept as true absence, status "proposed" in
  `mapping/ish_section_exceptions.csv`. To review on `figures\qc\00_flagged.png`
  and the gene's sheets; the numbers are final once reviewed.
- **`adult/arms.py` and `run_adult_arms.py`.** No step of the ISH line reads
  their table now that `ish/arms.py` retires. Proposal: retire them with it.
  Until decided they stay, out of the run order.
- **The presynaptic set.** GO annotates Syn1, Vamp2, Stx1a, Bsn, Syt1, Cplx1,
  Camk2a and Slc17a7 to both sides of the synapse, so the "not both" rule keeps
  16 presynaptic genes, mostly inhibitory and neuromodulatory. The sets stay as
  committed; the pre- and postsynaptic genes are drawn as context in figure 08.
- **The localisation test's positive control** fails with the GO control pool
  and works with that of 5 October. Both are shown; neither turns the
  localisation result positive.
- **Part 1's benchmark.** The leftover is above the floor and inside what one
  Allen Gria1 experiment leaves. Which benchmark a claim uses is to agree with
  Sami El-Boustani before the figure is shown.

## 9. How to rerun

In `tools\venv_atlas`, from the code root, after steps 1 to 10 of the route:

```
tools\venv_atlas\Scripts\python.exe mapping\run_structure_set.py
tools\venv_atlas\Scripts\python.exe mapping\run_ish_section_qc.py
...
tools\venv_atlas\Scripts\python.exe mapping\run_ish_overview.py
```

| step | script | needs | time (check tree, 8 October) |
|---|---|---|---|
| 11, 12 | `run_panel_build`, `run_panel_fetch` | the network, once; not rerun for this analysis | |
| 13 | `run_structure_set` | the per-brain files; `--recompute` redoes the per-adult channel table | 10 s (about 10 minutes with `--recompute`) |
| 14 | `run_ish_section_qc` | the Allen API for the repair's experiment lists, cached in `cache\`; `--offline` stops instead; `--sheets` draws the missing QC sheets | 22 s (about 15 minutes with `--sheets` from none) |
| 15 | `run_ish_gene_table` | mygene.info and the GO ontology, cached in `cache\`; `--offline` stops instead | 4 minutes |
| 16 | `run_ish_spatial_null` | the surrogates are cached; `--recompute` draws them again | 14 minutes, the calibration |
| 17 to 20 | `run_ish_gene_ranking`, `run_ish_robustness`, `run_ish_divisions --sheets`, `run_ish_gene_sets` | the surrogates | 99, 57, 158 and 29 s |
| 21 to 26 | `run_beyond_density` to `run_sep_channel_check` | | 3 to 37 s each |
| 27 | `run_ish_overview` | every step's numbers | 7 s |

A full run of steps 13 to 27 takes about 27 minutes, and a second run gives
the same tables.

Every run prints the data root and the settings in force. `SEP_DATA_ROOT`
moves the data root, for a copy. The settings are `[structures]`, `[ish_qc]`,
`[ish_analysis]`, `[spatial_null]`, `[beyond]`, `[beyond_calibration]`,
`[beyond_figures]` and `[ish_figures]` of `mapping/settings.toml`. The tests:
`cd mapping`, `..\tools\venv_dev\Scripts\python -m pytest tests`
([mapping/tests/README.md](../mapping/tests/README.md)).

The outputs of 5 October in `adult_v2\ish\`, `arms\`, `beyond\` and `panel\`
are frozen; nothing here writes there but steps 11 and 12 (`panel\`). The run
scripts that still write into them (`run_ish_regions`, `run_ish_compare`,
`run_ish_words`, `run_ish_roles`, `run_ish_arms`, `run_ish_reliability`,
`run_ish_panel_test`, `run_adult_arms`) are out of the run order and move to
`archive/` next.

## 10. Files

Under `<data>\adult_v2\ish_analysis\`:

| file | step | what |
|---|---|---|
| `tables\structure_set.csv`, `zref_reference.csv`, `centroids.csv` | 13 | every structure with its division, adults and reason; each adult's zero and spread; one-hemisphere centroids (mm) |
| `tables\adult_per_mouse.csv`, `adult_profile.csv` | 13 | per adult and structure, voxels and log2 means of nano, autofluorescence and SEP, plain and eroded, and `zref`; the cohort profile |
| `tables\experiments.csv`, `section_qc.csv`, `experiment_qc.csv` | 14 | every experiment; every section judged, with its flag; per experiment, the sections set missing and kept |
| `tables\gene_region_table.csv`, `gene_table.csv`, `gene_profiles.csv`, `gene_documentation.csv` | 15 | per experiment and structure, ISH means; per experiment, labels, QC and reliability; per gene and structure, the merged rank profile; the supplement table (UTF-8 with BOM) |
| `tables\surrogates_nano.npy`, `surrogates_auto.npy`, `surrogate_structures.csv`, `variogram.csv`, `null_calibration.csv` | 16 | the surrogates (10,000 x 204), their structures, the variograms, the calibration |
| `tables\gene_ranking.csv`, `gene_ranking_per_adult.csv`, `gap.csv`, `null_rho.npz` | 17 | per map and gene: rho, spatial p, q, null band, spread over adults, ranks; per adult; the gap; every gene's rho with every surrogate |
| `tables\ranking_robustness.csv`, `robustness_summary.csv` | 18 | per variant and gene; per variant |
| `tables\within_division.csv`, `within_division_detail.csv`, `within_calibration.csv` | 19 | per gene, division-only and within rho with both nulls; per gene and division; the two nulls on random maps |
| `tables\gene_sets.csv`, `set_tests.csv`, `contrasts.csv`, `localisation_test.csv`, `localisation_summary.csv` | 20 | set members; set tests; contrasts; partial rho per gene and pool; every test of the localisation design |
| `beyond\` | 21 to 25 | `structures_used.csv`, `variance_partition.csv`, `calibration.csv`, `controls.csv`, `regression_table.csv`, `residual_by_structure.csv`, `replication.csv`, `leftover_genes.csv`, `leftover_sets.csv`, `numbers_for_the_caption.txt`, working figures |
| `green_channel\sep_channel_check.csv` | 26 | per adult, each channel's range and correlations |
| `tables\april_headline.csv`, `april_anova.csv`, `april_groups.csv` | 27 | P9's genes then and now; the ANOVA under each choice; today's groups against the null |
| `tables\numbers_<step>.csv`, `numbers_for_the_text.csv` and `.txt` | 13 to 27 | the numbers of each step, and all of them |
| `cache\` | 14, 15 | Allen experiment lists, mygene records, `go-basic.obo` |

| figure | drawn by |
|---|---|
| `figures\00_overview.png`, `14_april_headline.png`, `README.md` | `run_ish_overview.py` |
| `figures\01_structures.png` | `run_structure_set.py` |
| `figures\02_genes.png` | `run_ish_gene_table.py` |
| `figures\03_beyond_budget.png`, `04_beyond_where.png` | `run_beyond_figures.py` |
| `figures\05_one_comparison.png`, `06_spatial_null.png`, `07_gene_ranking.png`, `11_autofluorescence.png` | `run_ish_gene_ranking.py` |
| `figures\08_gene_sets.png`, `09_localisation.png` | `run_ish_gene_sets.py` |
| `figures\10_between_within.png`, `genes\<gene>.png` | `run_ish_divisions.py` (`--sheets`) |
| `figures\12_robustness.png` | `run_ish_robustness.py` |
| `figures\13_green_channel.png` | `run_sep_channel_check.py` |
| `figures\qc\00_flagged.png`, `qc\<gene>_<experiment>.png` | `run_ish_section_qc.py` (`--sheets`) |

Reference: Burt JB, Helmer M, Shinn M, Anticevic A, Murray JD (2020).
Generative modeling of brain maps with spatial autocorrelation. NeuroImage 220,
117038.
