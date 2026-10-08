# The ISH analysis

What the adult nano map is, read against the Allen in situ hybridisation (ISH)
maps of the adult mouse brain: the question, the argument, each analysis with its
method, figure, result and meaning, the limits, and how to rerun it. Built on 8
October 2026 from the decisions of the ISH discussion (all five as recommended,
[history/ISH_DISCUSSION.md](history/ISH_DISCUSSION.md), section 8), and revised
the same day after three reviews (statistics, figures, facts).

- Code: `mapping/sepmap/structures.py`, `mapping/sepmap/adult/` and
  `mapping/sepmap/ish/`, run by steps 13 to 28 of the Python route
  ([mapping/README.md](../mapping/README.md)), on branch `post-ish`, not merged
  yet.
- Outputs: `<data>\adult_v2\ish_analysis\`. The numbers below are those of the
  run of 8 October 2026 on a full copy of the production inputs (the data root
  set by `SEP_DATA_ROOT`); the production data root holds the same once steps
  13 to 27 run there after the merge. Every number below is in
  `tables\numbers_for_the_text.csv` there (a `.txt` beside it reads more
  easily), written by `run_ish_overview.py` from the numbers each step writes;
  where a number comes from another table, the table is named.
- Figures: `figures\00_overview.png` to `15_april_headline.png`, PNG and EPS.
  `figures\README.md` walks through them in order, each with its question,
  what to look at, what to take from it and what it means, with the numbers of
  the run. The figures are not versioned, so this document names each by its
  file.

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
   postsynaptic-density genes, a sizeable part of the map is left over, and it
   replicates across mice. This is the core claim (section 4).
2. **The map is therefore something else, and the reading the data support is
   the surface fraction of the receptor**, shaped by trafficking regulation and
   scaffolding. The gene analyses ask whether they corroborate that reading:
   the genes that regulate surface AMPA receptors (Cacng8, a TARP; trafficking
   and scaffolding genes) should track the map better than abundance genes or
   unrelated genes, beyond a null that respects the brain's smoothness
   (section 5).

The limit, once: the surface fraction is the interpretation the data support,
not a measurement. A total-GluA1 stain on the same brains would measure it; the
green SEP channel cannot (section 6.3).

Where it stands today:

- **Part 1 holds against the calibration floor.** 27% of the map's
  reproducible pattern is not predicted by receptor mRNA or synaptic density,
  and it replicates across mice. On the same structures it is 14 points above
  what such a map would leave through Allen-to-Allen mismatch alone (the floor,
  15%), with an interval whose lower end is near zero (+0% to +28%).
- **It does not hold against the stricter benchmark.** A map that is one Allen
  Gria1 experiment leaves as much (34%) as the nano map does. Which benchmark a
  claim uses is still to agree (section 8).
- **Part 2 is not borne out by the gene tests.** The map follows Cacng8 and
  Gria1 beyond the null, inside divisions too, but Cacng8's lead over Gria1 is
  inside the null of maps related to both alike, no gene set passes, the
  localisation genes do no better than matched controls, and no single gene
  follows the leftover. The genes are consistent with the surface-fraction
  reading; they do not single it out.

What the comparison cannot do on its own: separate surface from total receptor
(both follow the postsynaptic side of synapses), read protein from mRNA (mRNA
sits in cell bodies, receptor on dendrites), or give a p value without a null
that keeps the brain's smooth gradients (cortex and hippocampus high, thalamus
and hypothalamus low). Each Allen experiment is one P56 mouse on a 200 um grid.

`figures\00_overview.png`: the question, the argument with this run's numbers
and where each part stands, one row per step with what stays open.

## 2. Words used here

| word | meaning |
|---|---|
| structure | one atlas region (CA1, VPM). Every comparison runs across structures, one value per structure on each side |
| rho | Spearman correlation across structures: do two maps put the structures in the same order |
| division | a coarse part of the brain: Isocortex, OLF, HPF, CTXsp, STR, PAL, TH, HY, MB, P, MY, CB |
| `zref` | per brain, log2 of each structure's nano mean over the brain's isocortex mean, minus its median over the declared set, over its p90 - p10 spread there; the adult map is the mean of the ten adults' `zref` |
| the declared set | the structures every comparison uses: grey matter measured in all ten adults (section 3) |
| beyond | short for "beyond abundance and density": how much of the map receptor mRNA and synaptic density predict, and what is left (part 1; the `beyond` folder and figures 03 and 04) |
| ceiling | the share of the map that two halves of the cohort reproduce (Spearman-Brown of their agreement) |
| leftover | the part of the map a fitted model does not predict |
| calibration floor | what the same model leaves of a map made only of receptor mRNA and density: the leftover Allen-to-Allen mismatch produces alone |
| jackknife | here, recomputing a number on subsamples that each leave out a fifth of the structures; the spread of those values, scaled for the subsample size, gives its interval |
| spatial null | random maps with the nano map's smoothness, the surrogates; a spatial p is the share of surrogates that correlate with a gene at least as strongly as the map does |
| surrogate | the map's ranks shuffled, smoothed over near structures and rescaled until its variogram matches the map's (Burt et al. 2020) |
| variogram | how unlike two structures are, on average, as a function of their distance: low for near structures in a smooth map, flat for a shuffled one |
| Freedman-Lane | a null for a residual: each surrogate is put through the same fit before it is compared, so it carries no more of the model's pattern than the residual does |
| partial rho | rho once a third map (the subunit composite) has been regressed out of both sides, on ranks |
| BH, q | Benjamini-Hochberg correction for testing many genes or sets; q is the corrected p, and "past BH" means q < 0.05 |
| positive control | a difference that must exist, run through the same test: if the test cannot find it, a null result means nothing |
| P9 | the MATLAB comparison of April 2026, `adult_matlab/run_compare_with_allen_ish.m` |

## 3. The inputs

Every analysis reads the same inputs, fixed before any correlation.

**The structures** (A1; `run_structure_set.py`, `figures\01_structures.png`).
A structure enters when it is grey matter and every adult measures it: of the
280 structures of the adult table, 235 are measured in all 10 adults and 204 of
those are grey matter (TH 45, HY 43, Isocortex 39, MB 25, STR 14, HPF 11, OLF
9, PAL 9, CTXsp 7, P 2). Of the 45 structures the rule leaves out for their
adults, 22 are seen in one to four adults, mostly pons and medulla, whose values
do not replicate across mice, and 23 in five to nine; all ten is the rule agreed
on 30 September (S1).

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
is below 0.2 of the median of its three neighbours on each side and below 0.2
of the brightest of the three on each side; it is never filled in. A section
dim against one side only sits at a step in expression (Cnih3 going from the
forebrain into the midbrain) and is kept; a section whose neighbours read below
0.1 energy, the noise of the grids, is not judged. Of 23,605 sections judged,
120 were set missing, in 96 experiments (10 of them P9's own); 16 sections in
13 experiments were kept at a step, 2,082 had neighbours at the noise level, and
76 experiments sit at that level throughout (labelled near zero in
`tables\experiment_qc.csv`). One dim section is kept as a true absence, from
`mapping/ish_section_exceptions.csv` (Glra1, section 61), proposed and not yet
reviewed (section 8). `figures\qc\00_flagged.png` lists every flag and step;
one sheet per experiment shows its profile and orientation.

**Reliability.** For the 222 genes measured more than once, two Allen
experiments of the same gene agree at a median rho of 0.689 (quartiles 0.516 to
0.772; 23 genes below 0.3). Gria1, Cacng8 and Dlg2 are at 0.91. A gene's rho
with the nano map is capped by its own reliability, so a low rho of an
unreliable gene says little.

**The gene sets** (fixed in `mapping/sepmap/ish/gene_sets.py` and committed
before any rho on these inputs): subunits (4 genes), localisation (84: the GO
terms of AMPA receptor transport, anchoring, clustering and the auxiliary
subunits), other postsynaptic (196), presynaptic (16), GABAergic neuron markers
(5) and glia (7, P9's list). The presynaptic set is small and mostly cell-type
and peptide markers, because GO annotates the vesicle machinery to both sides
of the synapse and the rule keeps a gene on one side only (section 8).

**The measured synapse density** (step 21; `run_synaptome.py`,
`synaptome\feasibility.png`). The density terms of part 1 are Allen mRNA, which
sits in cell bodies: a presynaptic marker's mRNA marks where the neurons that
make the synapses are, not where their synapses are. Zhu et al. (2018) counted
excitatory synapses where they are, in a knock-in mouse with PSD95 and SAP102
tagged: one adult male (about P80), five coronal sections, every punctum sorted
into 37 subtypes by intensity, size and shape, and the density of each subtype
(puncta per unit area) measured in regions of the Allen Reference Atlas.
Hansen et al. share the table: 37 subtypes by 775 samples, a sample being one
region of one hemisphere in one section. It is downloaded from their
repository at a pinned commit (`0399525`, files as in release v1.0, Zenodo doi
10.5281/zenodo.18201390) into `<data>\reference\synaptome\`, each file checked
against git's hash of it at that commit.

- *What is read.* Subtypes 1 to 11 hold PSD95 alone, 12 to 18 SAP102 alone and
  19 to 37 both (as Hansen et al.'s code indexes them). As shared, each
  subtype's density is scaled to 0..1 across the 775 samples, so absolute
  counts, and their sum, are gone. The PSD95 density is the mean over the 30
  subtypes whose puncta hold PSD95, each subtype's map weighing alike; never a
  punctum's intensity or size, since PSD95 per synapse is scaffolding, the
  surface side. Variants: PSD95 alone (Hansen et al.'s "PSD95 synapses"),
  SAP102 (26 subtypes) and every punctum (37). They order the structures of
  the fit at 0.86, 0.68 and 0.94 with the PSD95 density.
- *Placing the samples.* By the Allen id the source gives each sample (its
  acronym where it gives none), through the CCF 2017 ontology: a sample lies in
  the structure that is its id or the id's nearest ancestor. A unit (a layer,
  or the structure itself) is the mean of its samples over both hemispheres
  and every section; a structure is the mean of its sampled units weighted by
  their voxels in the 20 um CCF annotation, alike when a unit is not drawn
  there (9 structures of the fit, olfactory areas, the lateral entorhinal
  cortex and the subiculum, whose layers the annotation does not draw). A
  region the source gives only above several structures is never spread onto
  them, and a structure with no sample stays missing. Of the 775 samples, 739
  lie in 114 structures; left out are 8 above the structures (the midbrain's
  and the medulla's motor parts, the medial septal complex, the midbrain
  raphe), 12 that the 2017 annotation no longer draws (the layers of PTLp, now
  VISa and VISrl; the cochlear granular lamina), 15 whose names are not CCF
  acronyms (Mop, ZID, ZIV, DCOmo, CENT1) and the right locus coeruleus, with no
  punctum at all.
- *Coverage.* 77 of the 126 structures of the fit (61%) and 96 of the 204
  declared ones. Missing from the fit: 10 thalamic nuclei (AM, AV, CL, MD, PR,
  PT, PoT, RE, SMT, VAL), 9 cortical areas (AIp, FRP, SSp-ll, SSs, VISC, VISa,
  VISl, VISpm, VISrl), 7 hypothalamic, 6 midbrain (all under the midbrain's
  motor part), 5 pallidal, 4 striatal, 3 olfactory, 3 hippocampal and both
  pontine structures. Two are measured in less than half their volume (SCs
  from its zonal layer, BMA from its posterior part).
- *Agreement* (Spearman over the 77): with the marker mRNA composite 0.73, with
  `psd_pc1` 0.80, with Gria1 0.42, with the nano map 0.45, with
  autofluorescence 0.29. The two hemispheres agree at 0.96 (75 structures),
  the one check of reliability a single mouse allows.
- *The rule.* The measured density replaces the mRNA terms in the main model
  only if it covers at least 80% of the fit (101 of 126; `[beyond]
  min_psd95_coverage`, fixed before this was run). It covers 61%, so the main
  model keeps the mRNA terms on the 126 structures, and PSD95 enters as a
  variant on the 77 it covers.

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
- Each model is scored on structures it has not seen (cross-validated R²,
  averaged over 20 shufflings of five folds, since one shuffling alone moves the
  leftover by several points), as a share of the ceiling: two halves of the
  cohort, in all 126 splits, agree at 0.974, and Spearman-Brown makes that
  0.987 for ten adults. Left = 1 - share.
- The interval over structures is a jackknife: 400 subsamples each leave out a
  fifth of the structures, every predictor is rebuilt on the structures kept as
  the production run builds it, and the spread of the subsamples gives the 95%
  interval. A subsample draws no structure twice, so no structure sits in a
  training and a test fold at once (the bootstrap of 26 September let it, and
  flattered the fit).
- Replication: the model fitted to each half-cohort map, and the two leftovers
  correlated. The map's reliability and the fit alone imply a value; the
  observed one is set beside it.
- Calibration: the same model on maps whose answer is known, built from half of
  each gene's Allen experiments and predicted from the other half, with ten
  made-up adults as noisy as ours: a map made only of receptor mRNA and density
  (the floor), and a map that is one Allen Gria1 experiment. 113 structures,
  where both halves measure every subunit and marker. The nano map is read the
  same way on the same structures, and its difference from each known map is
  recomputed on each jackknife subsample, so the interval is that of the
  difference. A gene measured by one experiment (Shank2, Shank3 and Nlgn1 of
  the markers) is the same in both halves, so the floor errs low.
- Seven controls try to break the result: a spatial gradient, structure size,
  single animals, naive against RWS, curvature, the choice of predictors (the
  components of every gene measured in all the structures, judged against its
  own calibration floor), the reading. Beside them, the leftover under other
  folds and on other structures (`beyond\variants.csv`).

**Result** (`numbers_for_the_text.csv`, step `beyond`; `beyond\variance_partition.csv`,
`calibration.csv`, `calibration_jackknife.csv`, `controls.csv`, `variants.csv`):

| predictors, each bent, held out | share of the reproducible map |
|---|---|
| Gria1 alone | 47% |
| Gria1 to Gria4, four terms | 67% |
| density alone (markers, postsynaptic-density component) | 58% |
| autofluorescence alone | -2% |
| the model: Gria1 to Gria4, then density (+7%), then autofluorescence (-1%) | 73% |
| left | 27% (15% to 39% over structures) |

- The leftover replicates: the two half-cohort leftovers agree at 0.928 (lowest
  split 0.873). The map's reliability and the fit alone imply 0.883; two
  unrelated leftovers fall between -0.18 and +0.17.
- The calibration, on 113 structures: a map made only of receptor mRNA and
  density leaves 15% (13% to 17% over the draws), its leftover replicating at
  0.79: the floor. The nano map, read the same way, leaves 29%, 14 points more
  (95% +0% to +28% over resampled structures), about twice the floor. A map
  that is one Allen Gria1 experiment leaves 34% (33% when the truth comes from
  one half of its experiments, 35% from the other); nano minus that is -5%
  (-28% to +17%).
- The leftover under other choices (`beyond\variants.csv`): one shuffling of the
  folds 27% (single shufflings 25% to 33%), ten folds 27%, leave one out 26%,
  and 26% on the 159 structures kept once the markers measured once (Shank2,
  Shank3, Nlgn1) are left out. With folds of spatial blocks, which keep a
  structure's neighbours out of its fit, 42%; on the calibration's structures
  nano then leaves 40%, the floor 21% and the Gria1 map 40%: about twice the
  floor again, and as much as the Gria1 map.
- All 7 controls pass. Control F: the components of the 311 genes measured in
  every structure (most often 21, picked inside each training fold) predict 81%;
  on its own calibration nano leaves 22% to 26%, against 4% to 9% for a map made
  of those genes' expression.
- The composite model of 26 September, the four subunits averaged into one
  term, leaves 37% (24% to 49%): Gria4 runs against the map (rho -0.12,
  `tables\gene_ranking.csv`) and dilutes the average.
- Where it sits (`beyond\regression_table.csv`): above prediction in the medial
  geniculate (+50 ranks), the rostrolateral visual area, the septofimbrial
  nucleus and the subiculum; below it in the substantia innominata, piriform
  cortex, ventral retrosplenial cortex and VAL.

**What it means.** About a quarter of the map's reproducible pattern, 27%,
is not predicted by receptor mRNA or synaptic density, and it is reproducible
across mice. On the same structures it stands above what receptor mRNA and
density would leave through Allen-to-Allen mismatch, by a margin whose interval
reaches down to near zero, under random folds and under spatial blocks alike.

**What it does not mean.** It is no larger than what a map that is one Allen
Gria1 experiment leaves (34%): a map made of nothing but one experiment's
Gria1, with that experiment's own error, leaves as much. Its replication follows
from the map's reliability (0.883 expected), so it is not separate evidence. It
is "not predicted by receptor mRNA or synaptic density", not "beyond gene
expression": the components of many panel genes predict 81% of the map, though
no single gene or set follows the leftover (section 5.6). And a leftover says
what the predictors miss, not what it is; a claim about one structure needs its
own null.

## 5. Part 2: do the genes that set surface receptor follow the map better than abundance genes?

If the map is the surface fraction, genes that set how much receptor reaches and
stays at the membrane should follow it better than abundance genes or unrelated
genes. Three tests carry that question: the Cacng8 - Gria1 gap (section 5.3,
`figures\07_gene_ranking.png` B), the localisation genes against matched
controls once the subunit composite is removed (section 5.4,
`09_localisation.png`), and every gene against the leftover of part 1 (section
5.6, `11_leftover_genes.png`). The others describe the map: which genes look
like it, which kinds, and at what scale.

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
each structure's nearest neighbours (itself left out, as in brainsmash) with
the share of neighbours whose variogram fits the map's best (up to the 25th
percentile of the distances), and given the map's own values in its new order.
10,000 surrogates per map, nano and autofluorescence, on the 204 declared
structures, cut to each gene's structures. A spatial p is two-sided;
Benjamini-Hochberg runs within P9's genes and within all genes.

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

- 12 of 451 genes pass the spatial null after BH (q < 0.05), 6 of P9's 100;
  before correction 128 of 451 and 29 of P9's are past their band. The 12:
  Cacng8, Arpc5, Igsf11, Htr3a, Grm5, Neurl1a, Mapk1, Add3, Gria1, Cnih2, Grip1,
  Eps8. Four are localisation genes, one is a subunit; one of the four, Cnih2,
  has two Allen experiments that disagree (reliability 0.08), so its place is
  weak.
- Cacng8 is first, +0.807 (p ≤ 0.0001), and first of P9's genes in every
  robustness variant (section 6.2; under Pearson on log2 it is 2nd of all
  genes). Gria1 is +0.646 (p ≤ 0.0001, q 0.0025 within P9's genes), 11th of
  P9's genes and 32nd of all.
- The Cacng8 - Gria1 gap, on the 164 structures both have: +0.167, steady
  across adults (+0.156 to +0.179 over resampled adults; Cacng8 leads in 10 of
  10 adults). Tested against maps related to both genes alike, it is inside
  the null at its upper edge (p 0.105; 95% of that null -0.224 to +0.166);
  against maps unrelated to both, a wider null and a conservative bound, p
  0.348. It moves with the Allen experiment that stands for each gene, from
  +0.08 to +0.211 over the four pairings (p 0.064 to 0.385).

**What it means.** The map follows a TARP's pattern and Gria1's beyond a map
with the brain's smoothness. Cacng8 leads Gria1 in each of the 10 adults and in
every variant, but by a margin that maps following both genes alike reach about
one time in ten: these maps cannot tell which of the two the nano map follows
more closely. A description of the map, consistent with the surface-fraction
reading, not evidence for it.

### 5.4 Kinds of genes (analysis 3)

`run_ish_gene_sets.py`; `figures\08_gene_sets.png` and `09_localisation.png`.

**Method.** Each set's median rho against the median of the same genes' rho
with every surrogate, so co-expressed genes stay together in the null; BH over
the sets with at least 5 genes. Two contrasts named in advance, both needed for
the map to count as postsynaptic-like: the three postsynaptic sets pooled
against the presynaptic set, and against glia. Then the localisation test of 5
October on the new inputs: each gene's partial rho with the map once the
subunit composite is removed from both; the 84 localisation genes against 84
other postsynaptic genes matched on expression, labels permuted; the same with
the four subunits removed as separate terms, the abundance term of part 1; a
positive control (control genes with a reproducible map against those
without); the matched difference against maps whose remainder, once the
composite is removed, is a surrogate of the real remainder; and the test's
power: maps that carry the pattern setting the localisation genes apart from
their controls, in growing amounts, and how often the test finds it.

**Result** (step `gene_sets`; `tables\set_tests.csv`, `contrasts.csv`,
`localisation_summary.csv`, `localisation_power.csv`):

| set | genes | median rho | spatial p | q |
|---|---|---|---|---|
| subunits | 4 | +0.568 | too few to test | |
| localisation | 84 | +0.375 | 0.054 | 0.135 |
| other postsynaptic | 196 | +0.407 | 0.050 | 0.135 |
| presynaptic | 16 | +0.124 | 0.256 | 0.427 |
| GABAergic markers | 5 | +0.022 | 0.874 | 0.874 |
| glia | 7 | +0.063 | 0.573 | 0.717 |

- No set passes after BH.
- Postsynaptic above presynaptic: +0.274, spatial p 0.016. Postsynaptic above
  glia: +0.336, spatial p 0.170 (7 glial genes). One of the two contrasts
  named in advance passes, so the criterion (both) is not met. And the
  presynaptic side here is 16 genes, mostly cell-type and peptide markers (Npy,
  Cck, Th, Calb2, Slc17a6, Penk, Gad1, Gad2, Slc32a1), not the vesicle
  machinery, so the first contrast reads "postsynaptic genes against these
  markers".
- Localisation against matched controls, the subunit composite removed:
  +0.0041, p 0.916 (spatial p 0.885); a difference beyond ±0.082 would reach p
  < 0.05. With the four subunits removed as separate terms: -0.0276, p 0.546.
  The power check: maps carrying the localisation pattern give a matched
  difference the test finds in 80% of maps from +0.077; on maps with no such
  pattern it gives p < 0.05 on 0.2% of them, so the label test errs on the
  conservative side.
- The positive control of April's design does not come out with the GO control
  pool (+0.043, p 0.431); with the control pool of 5 October, added after that
  (secondary), it does (+0.111, p 0.0117), and localisation is still no better
  than its controls (+0.0073, p 0.839). The localisation genes with the highest
  partial rho are Cacng8, Igsf11, Dlg2, Dlg3 and Lgi1 (`localisation_test.csv`).

**What it means.** The map is postsynaptic-like only by half the criterion
named in advance. The genes that set surface AMPA receptors do not stand out
from other postsynaptic genes of the same expression: no localisation advantage
larger than about +0.08 in median partial rho, beyond the null.

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
  divisions the agreement falls to 0.576, and the median rho from +0.347 to
  +0.105.
- 11 genes pass the within null (7 of P9's): Cacng8 +0.548, Dlg2 +0.513,
  Neurl1a, Htr3a, Gria1 +0.417, Nptx1, Arrb2, Arpc5, Slc8a2, Tfrc, Rap2a.
- The shuffle inside divisions would have passed 32.2% of random smooth maps
  (the surrogates 2.6%), so it is shown and not used.

**What it means.** Cacng8, Dlg2 and Gria1 follow the map at a fine scale, not
only through the cortex-against-thalamus contrast: the stronger claim, and it
holds for both a TARP and the receptor's own mRNA.

### 5.6 The leftover itself

`figures\11_leftover_genes.png`; `beyond\leftover_genes.csv`, `leftover_sets.csv`.

**Method.** Every gene against the leftover of part 1, each against the
leftover's own surrogates, each surrogate put through the same fit first (the
leftover carries nothing of the model's columns, and a surrogate that did would
give a null far too wide). The gene sets are read the same way. None of these
tests was named before the leftover was seen.

**Result.** 0 of 451 genes pass after BH; 23 have a spatial p below 0.05, about
what chance gives. The closest is Aqp4, an astrocyte gene (+0.332, p 0.0004, q
0.18); Cacng8 is +0.148 (p 0.0024, 41st). The glia set passes (median +0.189, p
0.008, q 0.04); localisation does not (+0.032, p 0.246).

**What it means.** No single gene follows what receptor mRNA and density leave.
The glia set does, which may say something about the leftover (an astrocytic or
glial part of the tissue that the label reaches, or that shapes it), but the
test was not named in advance: a lead to test afresh, not a finding. Cacng8's
small uncorrected p is in the same position.

### 5.7 What part 2 says

Consistent with the surface-fraction reading, and not singling it out. The map
follows receptor expression and not the tissue (section 6.1), and Cacng8, Dlg2
and Gria1 lead the ranking inside divisions too. But Cacng8's lead over Gria1
is inside the null of maps related to both alike, no gene set passes its null,
the postsynaptic criterion named in advance is met by half, the localisation
genes predict the map no better than matched postsynaptic controls (no
advantage larger than about +0.08), and no single gene follows the leftover.
The genes neither contradict nor confirm part 2; only a total-receptor
measurement can.

## 6. Controls and limits

### 6.1 The tissue: autofluorescence (A8)

`figures\12_autofluorescence.png`. The autofluorescence map of the same
sections, read as nano is and tested with its own surrogates.

- With Gria1 +0.176 (p 0.447), with Cacng8 +0.296 (p 0.292): neither passes.
  Adult by adult, nano's rho with Gria1 is 0.482 to 0.754, autofluorescence's
  -0.266 to +0.297, and nano is above autofluorescence for both genes in all 10
  adults (`tables\gene_ranking_per_adult.csv`).
- The tissue has a gene pattern of its own. After BH, 38 genes pass the
  autofluorescence null against 12 for nano (ribosomal proteins, Mapt, Eno2 and
  GABA-A subunits among them); before correction it is the other way round,
  94 past their band against 128. BH counts only the smallest p values, and
  autofluorescence has more very small ones. The two gene orders agree at only
  0.363.

So Gria1 and Cacng8 follow the label, not the tissue, in every adult; the
autofluorescence ranking is a different one, though not an empty one.

### 6.2 The choices made: robustness (A3)

`run_ish_robustness.py`; `figures\13_robustness.png`; `tables\robustness_summary.csv`.
Eleven variants of the primary ranking (Spearman, `zref` with the declared
reference, merged profiles after QC, full means, the declared set): ten change
one choice, the last is the route of 5 October as it ran.

| variant | gene order against the primary (P9's genes) | Gria1's rank | gap |
|---|---|---|---|
| Pearson on log2 | 0.960 | 15 | +0.191 |
| ISH means eroded | 0.995 | 6 | +0.151 |
| both sides eroded | 0.994 | 6 | +0.144 |
| `ratio` reading | 0.966 | 8 | +0.132 |
| P9's single experiment per gene | 0.986 | 17 | +0.211 |
| P9's nine divisions, any number of adults | 0.993 | 15 | +0.180 |
| every structure of the adult table | 0.941 | 5 | +0.126 |
| the route of 5 October | 0.938 | 9 | +0.144 |

The nano erosion, the stored `zref` and no section QC change nothing (0.9996
and above). Cacng8 is first of P9's genes in every variant. Gria1's rank is
what moves (5th to 17th): quote it with its null, not as a place. The gap stays
between +0.126 and +0.211, inside the band of its null (figure 13 C).

### 6.3 The limit: the green channel (analysis 5)

`run_sep_channel_check.py`; `figures\14_green_channel.png`;
`green_channel\sep_channel_check.csv`.

The tag's own green (SEP) fluorescence was meant to show all tagged receptor,
so that nano over SEP would be the surface fraction. In every adult the green
channel follows autofluorescence across the 204 structures (rho 0.723 to
0.827), while nano follows it at only 0.047 to 0.529, and the green channel
varies across them about as little (p90 - p10 of log2: SEP 0.915,
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
- The calibration floor errs low (genes measured once are the same in both
  halves), and part 1's leftover is only as large as what one Allen Gria1
  experiment leaves.
- The exceptions list (true absences) and the sections kept at a step are a
  human call, still proposed.
- No total-receptor channel, no knockout or no-primary control, no DAPI
  control through the chain.

## 7. What changed from April

`run_ish_overview.py`; `figures\15_april_headline.png`; `tables\april_anova.csv`,
`april_groups.csv`.

P9's headline (22 April) grouped 97 genes by the categories of
`gene_targets.csv`, written while the genes were chosen, with "auxiliary"
split by hand into Aux forebrain (Cacng8, Cacng3, Cnih2, Cnih3, Grm5) and Aux
other, and tested the ten groups with a one-way ANOVA: p 0.032.

- The gene order reproduces: today's rho against April's agrees at 0.959 over
  97 genes.
- The group p rests on the split: 0.20 without it. On today's rho it is 0.011
  with the split and 0.149 without; 0.002 to 0.033 over structure sets and
  P9's single experiments.
- Every such p treats co-expressed genes as independent draws. Today's F across
  April's groups, against the F of the map's surrogates, which keep co-expressed
  genes together, gives p 0.289. One of the nine groups with enough genes
  passes its own null after BH: Aux forebrain, the group written after looking.

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
| tests | one ANOVA over genes | a spatial null for every rho, the gap, the sets, the contrasts and the leftover |
| controls | none | autofluorescence through the same chain, robustness, positive controls, a power check |

## 8. Points to settle

- **The rebuild of part 1, its rules fixed first (8 October).** Written down
  and committed before any of its results:
  - the main model is Gria1 + synapse density + autofluorescence, each bent
    and cross-validated as now (`[beyond]` of `mapping/settings.toml`).
    Abundance is Gria1 alone: the stained protein is GluA1, and Gria2 to Gria4
    make partners the nanobody does not see, whose availability sets GluA1's
    assembly and trafficking, the surface side. The four subunits stay as a
    check row;
  - synapse density is the measured PSD95 punctum density (Zhu et al. 2018, as
    shared by Hansen et al. 2026) in place of the mRNA panel (markers and
    `psd_pc1`), if it covers at least 80% of the structures of the fit; below
    that, the panel stays and PSD95 is a variant. Puncta count excitatory
    synapses where they are, while mRNA sits in somata;
  - the variants (`[beyond.variants]`) show how much the choice of density
    measure matters. One measured map in place of the panel's two terms gives
    the model fewer terms, so PSD95 and the panel together is the conservative
    bound. The calibration runs the same main model; PSD95 is one measurement,
    the same in both halves, like the markers measured once, so the floor errs
    low;
  - against the new leftover, Cacng8 is the one gene named in advance, and the
    AMPA receptor complex family (Schwenk et al. 2012 and GO:0032281, with
    Gria2 to Gria4, without Gria1) is tested as a group, then gene by gene
    within it (`mapping/sepmap/ish/gene_sets.py`). Cacng8's p against the
    leftover of the four-subunit model (section 5.6) was seen before it was
    named.
- **The exceptions list.** Glra1 section 61 is kept as true absence, status
  "proposed" in `mapping/ish_section_exceptions.csv`, and 16 sections are kept
  at a step in expression by the rule. To review on `figures\qc\00_flagged.png`
  (red, hatched, dots) and the gene's sheets; the numbers are final once
  reviewed.
- **Part 1's benchmark.** The leftover is above the floor and inside what one
  Allen Gria1 experiment leaves. Which benchmark a claim uses is to agree with
  Sami El-Boustani before the figure is shown.
- **The glia lead.** The glia set follows the leftover (q 0.04), a test not
  named in advance. Whether to test it afresh, and how.
- **`adult/arms.py` and `run_adult_arms.py`.** No step of the ISH line reads
  their table now that `ish/arms.py` retires. Proposal: retire them with it.
  Until decided they stay, out of the run order.
- **The presynaptic set.** GO annotates Syn1, Vamp2, Stx1a, Bsn, Syt1, Cplx1,
  Camk2a and Slc17a7 to both sides of the synapse, so the "not both" rule keeps
  16 presynaptic genes, mostly cell-type and peptide markers. The sets stay as
  committed; the pre- and postsynaptic genes are drawn as context in figure 08.
- **The localisation test's positive control** fails with the GO control pool
  and works with that of 5 October. Both are shown, with the power check;
  neither turns the localisation result positive.
- **The production run.** The outputs quoted here were made on a full copy of
  the production inputs; steps 13 to 28 run on the production data root after
  the merge.

## 9. How to rerun

In `tools\venv_atlas`, from the code root, after steps 1 to 10 of the route:

```
tools\venv_atlas\Scripts\python.exe mapping\run_structure_set.py
tools\venv_atlas\Scripts\python.exe mapping\run_ish_section_qc.py
...
tools\venv_atlas\Scripts\python.exe mapping\run_ish_overview.py
```

| step | script | needs | time (development copy, 8 October) |
|---|---|---|---|
| 11, 12 | `run_panel_build`, `run_panel_fetch` | the network, once; not rerun for this analysis | |
| 13 | `run_structure_set` | the per-brain files; `--recompute` redoes the per-adult channel table | 10 s (about 10 minutes with `--recompute`) |
| 14 | `run_ish_section_qc` | the Allen API for the repair's experiment lists, cached in `cache\`; `--offline` stops instead; `--sheets` draws the missing QC sheets | 20 s (about 15 minutes with `--sheets` from none) |
| 15 | `run_ish_gene_table` | mygene.info and the GO ontology, cached in `cache\`; `--offline` stops instead | 5 minutes |
| 16 | `run_ish_spatial_null` | the surrogates are cached; `--recompute` draws them again | 14 minutes, the calibration |
| 17 to 20 | `run_ish_gene_ranking`, `run_ish_robustness`, `run_ish_divisions --sheets`, `run_ish_gene_sets` | the surrogates | 2, 1, 3 and 1 minutes |
| 21 | `run_synaptome` | the synaptome of Zhu et al. 2018, downloaded once into `<data>\reference\synaptome\`; `--offline` stops instead | 10 s |
| 22 to 27 | `run_beyond_density` to `run_sep_channel_check` | | under 1 minute each; the calibration 2 minutes |
| 28 | `run_ish_overview` | every step's numbers | 10 s |

A full run of steps 13 to 28 takes about half an hour, and a second run gives
the same tables.

Every run prints the data root and the settings in force. `SEP_DATA_ROOT`
moves the data root, for a copy. The settings are `[structures]`, `[ish_qc]`,
`[ish_analysis]`, `[spatial_null]`, `[beyond]`, `[beyond_controls]`,
`[beyond_calibration]`, `[beyond_figures]` and `[ish_figures]` of
`mapping/settings.toml`. The tests: `cd mapping`,
`..\tools\venv_dev\Scripts\python -m pytest tests`
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
| `tables\experiments.csv`, `section_qc.csv`, `experiment_qc.csv` | 14 | every experiment; every section judged, with its flag and both references; per experiment, the sections set missing and kept, its median level and whether it is near zero |
| `tables\gene_region_table.csv`, `gene_table.csv`, `gene_profiles.csv`, `gene_documentation.csv` | 15 | per experiment and structure, ISH means; per experiment, labels, QC and reliability; per gene and structure, the merged rank profile; the supplement table (UTF-8 with BOM) |
| `tables\surrogates_nano.npy`, `surrogates_auto.npy`, `surrogate_structures.csv`, `variogram.csv`, `null_calibration.csv` | 16 | the surrogates (10,000 x 204), their structures, the variograms, the calibration |
| `tables\gene_ranking.csv`, `gene_ranking_per_adult.csv`, `gap.csv`, `null_rho.npz` | 17 | per map and gene: rho, spatial p, q, null band, spread over adults, ranks; per adult; the gap with both nulls; every gene's rho with every surrogate, and the gap of every null map |
| `tables\ranking_robustness.csv`, `robustness_summary.csv` | 18 | per variant and gene; per variant |
| `tables\within_division.csv`, `within_division_detail.csv`, `within_calibration.csv` | 19 | per gene, division-only and within rho with both nulls; per gene and division; the two nulls on random maps |
| `tables\gene_sets.csv`, `set_tests.csv`, `contrasts.csv`, `localisation_test.csv`, `localisation_summary.csv`, `localisation_power.csv` | 20 | set members; set tests; contrasts; partial rho per gene and pool; every test of the localisation design; its power by effect size |
| `synaptome\samples.csv`, `density.csv`, `coverage.csv`, `agreement.csv`, `feasibility.png` | 21 | per sample of the synaptome, its ids, densities and structure or why none; per structure of the adult table, measured or why not, its units, weights and densities; per division, how many declared and fitted structures are measured; each density's Spearman with the mRNA density terms, Gria1 and the maps |
| `beyond\` | 22 to 26 | `structures_used.csv`, `variance_partition.csv`, `calibration.csv`, `calibration_jackknife.csv`, `jackknife.csv`, `controls.csv`, `gene_space.csv`, `gene_space_calibration.csv`, `gene_space_summary.csv`, `variants.csv`, `regression_table.csv`, `residual_by_structure.csv`, `replication.csv`, `leftover_genes.csv`, `leftover_sets.csv`, `numbers_for_the_caption.txt`, working figures |
| `green_channel\sep_channel_check.csv` | 27 | per adult, each channel's range and correlations |
| `tables\april_headline.csv`, `april_anova.csv`, `april_groups.csv` | 28 | P9's genes then and now; the ANOVA under each choice; today's groups against the null |
| `tables\numbers_<step>.csv`, `numbers_for_the_text.csv` and `.txt` | 13 to 28 | the numbers of each step, and all of them |
| `cache\` | 14, 15 | Allen experiment lists, mygene records, `go-basic.obo` |

| figure | drawn by |
|---|---|
| `figures\00_overview.png`, `15_april_headline.png`, `README.md` | `run_ish_overview.py` |
| `figures\01_structures.png` | `run_structure_set.py` |
| `figures\02_genes.png` | `run_ish_gene_table.py` |
| `figures\03_beyond_budget.png`, `04_beyond_where.png`, `11_leftover_genes.png` | `run_beyond_figures.py` |
| `figures\05_one_comparison.png`, `06_spatial_null.png`, `07_gene_ranking.png`, `12_autofluorescence.png` | `run_ish_gene_ranking.py` |
| `figures\08_gene_sets.png`, `09_localisation.png` | `run_ish_gene_sets.py` |
| `figures\10_between_within.png`, `genes\<gene>.png` | `run_ish_divisions.py` (`--sheets`) |
| `figures\13_robustness.png` | `run_ish_robustness.py` |
| `figures\14_green_channel.png` | `run_sep_channel_check.py` |
| `figures\qc\00_flagged.png`, `qc\<gene>_<experiment>.png` | `run_ish_section_qc.py` (`--sheets`) |

Reference data, under `<data>\reference\`, each folder with a `fetch_log.txt`
(source, version, checksums, what was read):

| folder | what | source |
|---|---|---|
| `synaptome\` | PSD95 and SAP102 punctum density, 37 subtypes in 775 samples (a region of one hemisphere in one section), each subtype scaled to 0..1 as shared (Zhu et al. 2018); written by `run_synaptome.py` | github.com/netneurolab/hansen_synaptome, `data/synaptome/mouse_liu2018/` at commit `0399525412b6f50cfdeb5904b96da7fa8e4b507c`; archived as Zenodo doi 10.5281/zenodo.18201390 |
| `go\` | the GO Consortium's mouse annotation, and its lines for GO:0032281 | `current.geneontology.org/annotations/mgi.gaf.gz`, release 2026-08-05 |
| `ampar_complex\` | Schwenk et al. 2012: Figures 1 to 6 and the supplement (Tables S1 to S4) | the publisher's file server, `ars.els-cdn.com` |

References: Burt JB, Helmer M, Shinn M, Anticevic A, Murray JD (2020).
Generative modeling of brain maps with spatial autocorrelation. NeuroImage 220,
117038. Freedman D, Lane D (1983). A nonstochastic interpretation of reported
significance levels. Journal of Business and Economic Statistics 1, 292-298.
Hansen JY, Luppi AI, Qiu Z, Gini S, Fulcher BD, Gozzi A, et al. (2026). Synapse
types are spatially associated with regional hemodynamics in the mouse brain.
PLOS Biology 24, e3003637. Schwenk J, Harmel N, Brechet A, Zolles G, Berkefeld
H, Müller CS, et al. (2012). High-resolution proteomics unravel architecture
and molecular diversity of native AMPA receptor complexes. Neuron 74, 621-633.
Zhu F, Cizeron M, Qiu Z, Benavides-Piccione R, Kopanitsa MV, Skene NG, et al.
(2018). Architecture of the mouse brain synaptome. Neuron 99, 781-799.
