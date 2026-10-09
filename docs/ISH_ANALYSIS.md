# The ISH analysis

What the adult nano map is, read against the Allen in situ hybridisation (ISH)
maps of the adult mouse brain: the question, the argument, each analysis with its
method, figure, result and meaning, the limits, and how to rerun it. Built on 8
October 2026 from the decisions of the ISH discussion (all five as recommended,
[history/ISH_DISCUSSION.md](history/ISH_DISCUSSION.md), section 8), revised the
same day after three reviews (statistics, figures, facts), and run again on 9
October with part 1's main model, fixed on 8 October before its results
(section 8), and with the genes that follow the map described one by one and the
tests named for the leftover run (sections 5.7 and 5.8).

- Code: `mapping/sepmap/structures.py`, `mapping/sepmap/adult/` and
  `mapping/sepmap/ish/`, run by steps 13 to 29 of the Python route
  ([mapping/README.md](../mapping/README.md)), on branch `post-ish`, not merged
  yet.
- Outputs: `<data>\adult_v2\ish_analysis\`. The numbers below are those of the
  run of 9 October 2026 on a full copy of the production inputs (the data root
  set by `SEP_DATA_ROOT`); the production data root holds the same once steps
  13 to 29 run there after the merge. Every number below is in
  `tables\numbers_for_the_text.csv` there (a `.txt` beside it reads more
  easily), written by `run_ish_overview.py` from the numbers each step writes;
  where a number comes from another table, the table is named.
- Figures: `figures\00_overview.png` to `18_april_headline.png`, PNG and EPS.
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
   density.** Predicted from Gria1 mRNA, synapse density and the tissue's
   autofluorescence, a sizeable part of the map is left over, and it
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

- **Part 1 holds against the calibration floor at its point value.** 29% of
  the map's reproducible pattern is not predicted by Gria1 expression or synapse
  density (15% to 42% over structures), and it replicates across mice. On the
  same structures it is 14 points above what such a map would leave through
  Allen-to-Allen mismatch alone (the floor, 12%), about twice it, with an
  interval that reaches just below zero (-2% to +30%). The measured PSD95
  density, where it exists, leaves more, not less.
- **It does not hold against the stricter benchmark.** A map that is one Allen
  Gria1 experiment leaves as much (34%) as the nano map does (26% on the same
  structures). Which benchmark a claim uses is still to agree (section 8).
- **Part 2 is not borne out by the gene tests of the map, and one test of the
  leftover points its way.** The map follows Cacng8 and Gria1 beyond the null,
  inside divisions too, but Cacng8's lead over Gria1 is inside the null of maps
  related to both alike, no gene set passes, and the localisation genes do no
  better than matched controls. The genes that follow the map most are maps
  much like Gria1 and synapse density (section 5.7). Against the leftover of the
  main model, Cacng8, the one gene named in advance, follows it beyond its null
  (p 0.0002); added to the main model it takes 4.6% of the reproducible map,
  more than maps that relate to the model as it does (p 0.009), not more than
  plain surrogates of it (p 0.059). The AMPA receptor complex family named with
  it follows the leftover beyond the surrogates (p 0.047) but no more than
  postsynaptic genes of the same expression (p 0.55), and inside the family
  only Cacng8 passes BH (section 5.8). The genes are consistent with the
  surface-fraction reading; they do not single it out.

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
| beyond | short for "beyond abundance and density": how much of the map Gria1 expression and synapse density predict, and what is left (part 1; the `beyond` folder and figures 03 and 04) |
| main model | part 1's model, fixed on 8 October 2026: Gria1 + synapse density + autofluorescence (section 4) |
| check row | a variant of the main model that changes one thing of it, fixed with it and reported beside it (`[beyond.variants]`) |
| ceiling | the share of the map that two halves of the cohort reproduce (Spearman-Brown of their agreement) |
| leftover | the part of the map a fitted model does not predict |
| tier | a level of the tests against the leftover, named on 8 October: Cacng8 alone (1), the AMPA receptor complex family (2), every other gene (3, exploratory) |
| maps alike to the model | maps that keep a gene's fit on the main model and replace the rest with a surrogate of it, as large: a null for what a gene takes of the leftover |
| calibration floor | what the same model leaves of a map made only of Gria1 and synapse density: the leftover Allen-to-Allen mismatch produces alone |
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

## 4. Part 1: not explained by Gria1 expression or synapse density

`run_beyond_density.py` to `run_beyond_figures.py`; `figures\03_beyond_budget.png`
and `04_beyond_where.png`.

**The main model**, fixed on 8 October 2026 before any of its results
(`[beyond]` of `mapping/settings.toml`; section 8):

```
nano map  ~  Gria1 + synapse density + autofluorescence
```

- *Abundance is Gria1 alone.* The mice carry SEP-GluA1: the stained protein is
  GluA1, which Gria1 alone encodes. Gria2 to Gria4 encode partner subunits the
  nanobody does not see, and their availability sets GluA1's assembly and
  trafficking: the surface side, not abundance. The four subunits as four
  terms, the model of the first run of 8 October, are a check row.
- *Synapse density* is the measured PSD95 punctum density (section 3) if it
  covers at least 80% of the structures of the fit; otherwise the mRNA panel,
  the composite of 11 synaptic marker genes and the first principal component
  of 186 postsynaptic-density genes. PSD95 covers 77 of the 126 structures
  (61%), so the main model keeps the panel, on all 126, and PSD95 is a check
  row on the 77 it covers.

**Method.**

- The structures: the declared set, less those where a subunit or marker gene
  has no ISH value: 126 structures, 10 adults. They are those of the first run,
  so every check row runs on them or on a part of them.
- The map and every predictor as ranks across structures, each entering as x,
  x² and x³: 13 columns with the intercept.
- Each model is scored on structures it has not seen (cross-validated R²,
  averaged over 20 shufflings of five folds, since one shuffling alone moves the
  leftover by several points), as a share of the ceiling: two halves of the
  cohort, in all 126 splits, agree at 0.974, and Spearman-Brown makes that
  0.987 for ten adults. Left = 1 - share. The budget adds the terms in the
  order Gria1, synapse density, autofluorescence.
- The interval over structures is a jackknife: 400 subsamples each leave out a
  fifth of the structures, every predictor is rebuilt on the structures kept as
  the production run builds it, and the spread of the subsamples gives the 95%
  interval. A subsample draws no structure twice, so no structure sits in a
  training and a test fold at once.
- Replication: the model fitted to each half-cohort map, and the two leftovers
  correlated. The map's reliability and the fit alone imply a value; the
  observed one is set beside it.
- Calibration: the main model on maps whose answer is known, built from half of
  each gene's Allen experiments and predicted from the other half, with ten
  made-up adults as noisy as ours: a map made only of Gria1 and synapse density
  (the floor), and a map that is one Allen Gria1 experiment (the benchmark).
  113 structures, where both halves measure every subunit and marker, as the
  fit's own rule asks of the merged profiles. The nano map is read the same way
  on the same structures, and its difference from each known map is recomputed
  on each jackknife subsample, so the interval is that of the difference. A
  gene measured by one experiment (Shank2, Shank3 and Nlgn1 of the markers) is
  the same in both halves, so the floor errs low; the PSD95 density, one mouse
  measured once, would enter the same way if it were in the main model.
- Seven controls try to break the result: a spatial gradient, structure size,
  single animals, naive against RWS, curvature, the choice of predictors (the
  components of every gene measured in all the structures, judged against its
  own calibration floor), the reading.
- The check rows (`[beyond.variants]`, fixed with the main model), each one
  change to it, on its own structures with its ceiling recomputed there: the
  density as PSD95 puncta, SAP102 puncta or every punctum, on the 77 structures
  where they are measured; the mRNA panel on those same 77, so that the two
  kinds of density meet on the same structures; PSD95 and the panel together,
  the conservative bound, since one measured map in place of the panel's two
  terms gives the model fewer terms; the panel without its five presynaptic
  markers; Gria1 to Gria4 in place of Gria1. Beside them, the main model under
  other folds and on more structures (`beyond\variants.csv`).

**Result** (`numbers_for_the_text.csv`, step `beyond`; `beyond\variance_partition.csv`,
`calibration.csv`, `calibration_jackknife.csv`, `controls.csv`, `variants.csv`):

| predictors, each bent, held out | share of the reproducible map |
|---|---|
| Gria1 alone | 47% |
| synapse density alone (markers, postsynaptic-density component) | 58% |
| autofluorescence alone | -2% |
| Gria1 to Gria4 alone, four terms | 67% |
| the main model: Gria1 (47%), + synapse density (+25%), + autofluorescence (-1%) | 71% |
| left | 29% (15% to 42% over structures) |

- The leftover replicates: the two half-cohort leftovers agree at 0.933 (lowest
  split 0.874). The map's reliability and the fit alone imply 0.897; two
  unrelated leftovers fall between -0.17 and +0.17.
- The calibration, on 113 structures: a map made only of Gria1 and synapse
  density leaves 12% (11% to 15% over the draws), its leftover replicating at
  0.78: the floor. The nano map, read the same way, leaves 26%, 14 points more
  (95% -2% to +30% over resampled structures), about twice the floor. A map
  that is one Allen Gria1 experiment leaves 34% (35% when the truth comes from
  one half of its experiments, 34% from the other); nano minus that is -8%
  (-29% to +14%).
- The check rows (`beyond\variants.csv`; density alone and left are shares of
  the reproducible map on the row's own structures):

  | check row | structures | density alone | left |
  |---|---|---|---|
  | the main model | 126 | 58% | 29% |
  | Gria1 to Gria4 in place of Gria1 | 126 | 58% | 27% |
  | the panel without its presynaptic markers | 126 | 58% | 31% |
  | the mRNA panel, on the structures PSD95 covers | 77 | 45% | 39% |
  | PSD95 puncta | 77 | 22% | 46% |
  | PSD95 puncta and the mRNA panel | 77 | 42% | 39% |
  | SAP102 puncta | 77 | 26% | 48% |
  | every punctum | 77 | 24% | 47% |

  On the 77 structures where it is measured, the PSD95 density predicts the map
  less well than the mRNA panel does (22% against 45% alone) and leaves more
  (46% against 39%); the two together leave as much as the panel alone. So the
  measured density does not shrink the leftover: the panel predicts more of the
  map, and the main model, by keeping it, is the harder test of part 1. The 77
  structures leave more than the 126 with the same model (39% against 29%):
  which structures are in moves the leftover as much as the density measure
  does. The partner subunits change little once density is in (27% left
  against 29%).
- The main model under other folds: one shuffling 28% (single shufflings 26%
  to 31%), ten folds 28%, leave one out 28%, and 28% on the 159 structures kept
  once the markers measured once (Shank2, Shank3, Nlgn1) are left out. With
  folds of spatial blocks, which keep a structure's neighbours out of its fit,
  38%; on the calibration's structures nano then leaves 30%, the floor 16% and
  the Gria1 map 37%.
- All 7 controls pass. Control F: the components of the 311 genes measured in
  every structure (most often 21, picked inside each training fold) predict
  81%; on its own calibration nano leaves 22% to 26%, against 4% to 9% for a
  map made of those genes' expression.
- Where it sits (`beyond\regression_table.csv`): above prediction in the medial
  geniculate (+56 ranks), both parts of the lateral geniculate, the lateral
  posterior nucleus, the septofimbrial nucleus, the subthalamic nucleus and the
  subiculum; below it in the nucleus of reuniens, the medial habenula, dorsal
  retrosplenial cortex, the arcuate nucleus, piriform cortex and the submedial
  nucleus of the thalamus.

**What it means.** About 29% of the map's reproducible pattern is not predicted
by Gria1 expression or synapse density, and it is reproducible across mice. On
the same structures it is about twice what Allen-to-Allen mismatch leaves (26%
against the floor's 12%), a margin whose interval reaches just below zero
(-2%); with folds of spatial blocks the margin is the same (30% against 16%).
Measured synapse density, where it exists, leaves more, not less.

**What it does not mean.** It is no larger than what a map that is one Allen
Gria1 experiment leaves (34%): a map made of nothing but one experiment's
Gria1, with that experiment's own error, leaves as much. Its replication follows
from the map's reliability (0.897 expected), so it is not separate evidence. It
is "not predicted by Gria1 expression or synapse density", not "beyond gene
expression": the components of many panel genes predict 81% of the map. And a
leftover says what the predictors miss, not what it is; a claim about one
structure needs its own null.

## 5. Part 2: do the genes that set surface receptor follow the map better than abundance genes?

If the map is the surface fraction, genes that set how much receptor reaches and
stays at the membrane should follow it better than abundance genes or unrelated
genes. Three tests carry that question: the Cacng8 - Gria1 gap (section 5.3,
`figures\07_gene_ranking.png` B), the localisation genes against matched
controls once the subunit composite is removed (section 5.4,
`09_localisation.png`), and every gene against the leftover of part 1 (section
5.6, `11_leftover_genes.png`), where Cacng8 and the AMPA receptor complex family
were named in advance (section 5.8, `14_ampa_family.png`). The others describe
the map: which genes look like it, which kinds, at what scale, and, gene by
gene for those that follow it most, what kind of map each is (section 5.7,
`12_top_genes.png` and `13_cacng8_gria1.png`).

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

**Method.** Every gene against the leftover of part 1's main model, each
against the leftover's own surrogates, each surrogate put through the same fit
first (the leftover carries nothing of the model's columns, and a surrogate that
did would give a null far too wide). The gene sets are read the same way. Three
tiers were named for this leftover on 8 October, before the main model ran
(`mapping/sepmap/ish/gene_sets.py`): Cacng8 alone, its uncorrected spatial p the
test (its p against the leftover of the four-subunit model, 0.0024, was seen
before it was named); the AMPA receptor complex family, as a group and then
gene by gene within it (section 5.8); every other gene exploratory, BH over all
of them.

**Result** (`beyond\leftover_genes.csv`, `leftover_sets.csv`). Cacng8 follows
the leftover: +0.204, spatial p 0.0002 (28th of 451 genes). Over every gene, 0
of 451 pass after BH and 27 have a spatial p below 0.05, about what chance
gives; the closest is Aldh1l1, an astrocyte gene (+0.321, p 0.20, q 0.80). The
glia set has median +0.164 (p 0.043, q 0.19), localisation +0.048 (p 0.116).

**What it means.** The one gene named in advance, a TARP, follows what Gria1 and
synapse density leave, beyond the leftover's own null; no gene does once every
gene is corrected for. The glia set followed the leftover of the first,
four-subunit model (q 0.04), a test not named in advance; against this one it
does not pass, and it stays an unplanned lead, not pursued.

### 5.7 The genes that follow the map

`run_ish_top_genes.py`, `mapping/sepmap/ish/top_genes.py`;
`figures\12_top_genes.png`, `13_cacng8_gria1.png` and one sheet per gene,
`figures\top_genes\<gene>.png`.

**Method.** The 12 genes past the map's null after BH over every gene (section
5.3), with Gria1 and Cacng8 always and every member of the AMPA receptor complex
family (section 5.8), are each described on every axis of the line:

- with the map: rho, spatial p, q and rank; the mean rho inside divisions and
  its null (section 5.5); the range of its rho over the robustness variants,
  and of its rank over the variants that hold every gene (section 6.2);
- its Allen map: how many experiments, and how well they agree;
- whether it is just a Gria1-like or density-like map: its rho with Gria1, the
  two density terms, autofluorescence and the main model's prediction on the
  126 structures of the fit, and with the PSD95 density on those of them where
  it is measured;
- with the leftover: its rho and spatial p, each surrogate through the same fit
  (section 5.6);
- what it takes of the leftover: the main model with the gene's ranks added,
  bent as every term is, scored on held-out structures as the main model is, as
  a share of the reproducible map. A smooth map takes some by chance, so the
  gene is set beside 1,000 maps of its smoothness put in its place, of two
  kinds. Plain surrogates of the gene are unrelated to the model, so more of
  each is new to it than of a gene that shares the model's pattern: a wide null,
  the conservative bound, and the one this step was first given. Maps that
  relate to the model as the gene does keep the gene's fit on the model and
  replace its remainder with a surrogate of it, as large; they were added on 9
  October after the plain null's numbers were seen, and both are reported;
- what it is: its name, the GO terms that put it in the ontology panel and its
  cellular-component terms at the synapse (which carry SynGO's curation where it
  reached GO), its gene sets and P9's category, all from the gene table.

Nothing about a gene's biology is written here from memory but one line:
Cacng8 encodes TARP γ-8, an AMPA receptor auxiliary subunit expressed mainly in
the hippocampus, where it sets AMPA receptor protein levels and their
extrasynaptic surface expression (Rouach et al. 2005).

**Result** (step `top_genes`; `tables\top_genes.csv`). The genes past the map's
null, best first; "past" marks a within rho past its null after BH, "-" a gene
measured once; the share taken is of the reproducible map, with its p against
plain surrogates and against maps alike to the model:

| gene | with the map, rank | inside divisions | rank over the variants | reliability | with Gria1, with the model's prediction | with the leftover (p) | takes of the leftover (p plain, alike) |
|---|---|---|---|---|---|---|---|
| Cacng8 | +0.807, 1 | +0.548 (past) | 1 to 2 | 0.91 | +0.74, +0.80 | +0.204 (0.0002) | 4.6% (0.059, 0.009) |
| Arpc5 | +0.751, 2 | +0.356 (past) | 1 to 7 | - | +0.58, +0.76 | +0.208 (0.0005) | 3.8% (0.129, 0.010) |
| Igsf11 | +0.743, 3 | +0.210 | 2 to 15 | 0.85 | +0.55, +0.72 | +0.122 (0.069) | 2.2% (0.261, 0.107) |
| Htr3a | +0.737, 4 | +0.443 (past) | 3 to 31 | - | +0.57, +0.65 | +0.287 (0.007) | 5.9% (0.011, 0.002) |
| Grm5 | +0.725, 6 | +0.382 | 2 to 14 | - | +0.77, +0.81 | +0.078 (0.175) | -1.2% (0.886, 0.771) |
| Neurl1a (in psd_pc1) | +0.718, 8 | +0.467 (past) | 3 to 34 | 0.81 | +0.64, +0.69 | +0.134 (0.036) | 0.4% (0.526, 0.346) |
| Mapk1 | +0.687, 17 | +0.213 | 8 to 38 | - | +0.58, +0.69 | +0.168 (0.030) | 1.1% (0.383, 0.031) |
| Add3 | +0.681, 20 | +0.292 | 4 to 89 | - | +0.38, +0.61 | +0.260 (0.035) | 4.2% (0.088, 0.029) |
| Gria1 (the abundance term) | +0.646, 32 | +0.417 (past) | 18 to 47 | 0.91 | +1.00, +0.80 | +0.041 (0.036) | 0.0% |
| Cnih2 | +0.637, 36 | +0.336 | 5 to 43 | 0.08 | +0.56, +0.63 | +0.002 (0.978) | -1.1% (0.785, 0.724) |
| Grip1 (in psd_pc1) | +0.615, 45 | +0.207 | 22 to 75 | 0.65 | +0.74, +0.69 | +0.095 (0.133) | -0.1% (0.588, 0.440) |
| Eps8 | +0.605, 54 | +0.286 | 8 to 79 | 0.55 | +0.53, +0.61 | +0.136 (0.057) | 1.6% (0.220, 0.141) |

- They are maps much like Gria1 and synapse density: rho 0.38 to 0.77 with
  Gria1 and 0.61 to 0.81 with the main model's prediction. That is why the main
  model predicts most of the map.
- 6 of them follow the leftover before correction (Cacng8, Arpc5, Htr3a,
  Neurl1a, Mapk1, Add3); only Cacng8's p is a test named in advance. Gria1 is a
  term of the model, so its rho with the leftover is near zero by construction
  and its p says nothing.
- Added to the main model, Htr3a alone takes more than 95% of plain surrogates
  of it (5.9%, p 0.011); against maps alike to the model, Cacng8, Arpc5, Htr3a,
  Mapk1 and Add3 do (p 0.002 to 0.031).
- Cacng8 against Gria1 (figure 13): the map follows Cacng8 at +0.807 and Gria1
  at +0.646, inside divisions +0.548 and +0.417; the two genes agree at +0.74
  over the structures of the fit. The lead, +0.167, is inside the null of maps
  that follow both alike (p 0.105; section 5.3). Against what Gria1 and synapse
  density leave, Cacng8 gives +0.204 (p 0.0002), and added to the main model it
  takes 4.6% of the reproducible map (28.6% left, 24.0% with it).

**What it means.** The genes that follow the map most describe it; they are
not separate evidence, since most of what they share with the map is what Gria1
and synapse density already predict. A few carry part of what the model leaves.
Of these, Cacng8 is the one named in advance, and it holds against both kinds of
null for its rho, and against the narrower one for what it takes; the wider one
(p 0.059) is the bound a sceptic would quote. Htr3a, a serotonin receptor
measured by one Allen experiment, takes more than any other, an exploratory
finding.

### 5.8 The AMPA receptor complex family

`run_ish_top_genes.py`; `figures\14_ampa_family.png`; `tables\named_tests.csv`,
`top_genes.csv`, `family_members.csv`.

**Method.** Named on 8 October before the main model ran (section 8;
`mapping/sepmap/ish/gene_sets.py`), from sources outside this analysis: the
constituents of native AMPA receptor complexes found by proteomics (Schwenk et
al. 2012, Figure 1D and Table S2), the mouse genes annotated to GO:0032281 AMPA
glutamate receptor complex (GO release 2026-08-05), and the partner subunits
Gria2 to Gria4, less Gria1, the abundance term: 37 genes, 31 with a usable map.
Not tested, being in neither panel of the gene table: Gsg1l, Olfm3, Prrt1,
Prrt2, Rap2b and Shisa8. Tier 2, after Cacng8 alone (tier 1):

- the family's median rho with the leftover against the same genes' median
  over the leftover's 10,000 surrogates, each through the same fit;
- the family against 31 genes of the other postsynaptic set outside it, each
  the one closest in median energy (the localisation test's matching), labels
  permuted 20,000 times;
- then gene by gene, BH within the family only.

Every other gene is exploratory (section 5.6). The same is read on the map
itself, to describe the family, not as a test of the leftover. Every gene's rho
with the leftover was in `leftover_genes.csv` (step 22) before this step ran;
the group tests are those committed on 8 October.

**Result:**

| test | on the leftover | on the map (a description) |
|---|---|---|
| Cacng8 alone | +0.204, p 0.0002 | +0.807, p ≤ 0.0001 |
| the family's median against the surrogates | +0.064 (null 95% -0.063 to +0.063), p 0.047 | +0.414, p 0.044 |
| the family against matched controls | +0.064 against +0.046: +0.018, p 0.55 (±0.062 would pass) | +0.414 against +0.365: +0.049, p 0.50 |
| members past BH within the family | Cacng8 (q 0.006) | Cacng8 and Cnih2 (q 0.002) |

The members closest to the leftover after Cacng8 are Lrrtm4 (+0.314, p 0.074),
Vwc2l (+0.252), Shisa6 (+0.210) and Cpt1c (+0.193); Abhd12 (p 0.013) and Dlg3
(p 0.036) have the next smallest p, q 0.21 and 0.33 within the family. The
partner subunits sit near zero or below (Gria2 -0.064, Gria3 +0.062, Gria4
-0.105).

**What it means.** The family follows what Gria1 and synapse density leave
more than maps of the leftover's smoothness do, just past the null; but other
postsynaptic genes of the same expression follow it about as much, so what the
family shares with the leftover is postsynaptic, not particular to the AMPA
receptor complex. Inside the family it is Cacng8 that carries it.

### 5.9 What part 2 says

Consistent with the surface-fraction reading, and not singling it out. The map
follows receptor expression and not the tissue (section 6.1), and Cacng8, Dlg2
and Gria1 lead the ranking inside divisions too. But Cacng8's lead over Gria1
is inside the null of maps related to both alike, no gene set passes its null,
the postsynaptic criterion named in advance is met by half, and the
localisation genes predict the map no better than matched postsynaptic controls
(no advantage larger than about +0.08). The genes that follow the map most are
maps much like Gria1 and synapse density. Against the leftover, the TARP named
in advance follows it (p 0.0002) and takes part of it, the one result here that
points the way the reading does; the AMPA receptor complex family named with it
follows the leftover no more than postsynaptic genes of the same expression.
The genes neither contradict nor confirm part 2; only a total-receptor
measurement can.

## 6. Controls and limits

### 6.1 The tissue: autofluorescence (A8)

`figures\15_autofluorescence.png`. The autofluorescence map of the same
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

`run_ish_robustness.py`; `figures\16_robustness.png`; `tables\robustness_summary.csv`.
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
between +0.126 and +0.211, inside the band of its null (figure 16 C).

### 6.3 The limit: the green channel (analysis 5)

`run_sep_channel_check.py`; `figures\17_green_channel.png`;
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
- The measured synapse density is one mouse, in 77 of the 126 structures of
  the fit, as shared scaled per subtype: it orders the structures, it does not
  count synapses across them.
- The exceptions list (true absences) and the sections kept at a step are a
  human call, still proposed.
- No total-receptor channel, no knockout or no-primary control, no DAPI
  control through the chain.

## 7. What changed from April

`run_ish_overview.py`; `figures\18_april_headline.png`; `tables\april_anova.csv`,
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

  Run on 9 October: PSD95 covers 77 of the 126 structures, so the main model
  keeps the mRNA panel and PSD95 is a check row (section 4). The outputs of the
  first run of 8 October, with the four subunits, are kept on the development
  copy in `adult_v2\ish_analysis_8oct_four_subunits\`. The tests named for the
  leftover ran the same day (sections 5.7 and 5.8): Cacng8 follows it (p
  0.0002); the family passes against the surrogates (p 0.047), not against
  matched controls (p 0.55).
- **What a gene takes of the leftover.** Two nulls (section 5.7): plain
  surrogates of the gene, the null this step was first given and a wide one,
  and maps that relate to the model as the gene does, added after the first's
  numbers were seen. Cacng8 passes the second (p 0.009), not the first (p
  0.059). Which one a claim quotes is to agree; the first is the bound.
- **The exceptions list.** Glra1 section 61 is kept as true absence, status
  "proposed" in `mapping/ish_section_exceptions.csv`, and 16 sections are kept
  at a step in expression by the rule. To review on `figures\qc\00_flagged.png`
  (red, hatched, dots) and the gene's sheets; the numbers are final once
  reviewed.
- **Part 1's benchmark.** The leftover is above the floor at its point value
  (the interval of the difference reaches -2%) and inside what one Allen Gria1
  experiment leaves. Which benchmark a claim uses is to agree with Sami
  El-Boustani before the figure is shown.
- **The glia lead.** The glia set followed the leftover of the four-subunit
  model (q 0.04), a test not named in advance; against the main model's
  leftover it does not pass (q 0.19). An unplanned lead, not pursued now
  (figure 11 says so in one line).
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
  the production inputs; steps 13 to 29 run on the production data root after
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
| 22 to 26 | `run_beyond_density` to `run_beyond_figures` | | under 1 minute each; the calibration 2 minutes |
| 27 | `run_ish_top_genes` | the tables of steps 15 to 22; `--sheets` draws one sheet per gene | 8 minutes, the nulls of the share taken |
| 28 | `run_sep_channel_check` | | under 1 minute |
| 29 | `run_ish_overview` | every step's numbers | 10 s |

A full run of steps 13 to 29 takes about 40 minutes, and a second run gives
the same tables.

Every run prints the data root and the settings in force. `SEP_DATA_ROOT`
moves the data root, for a copy. The settings are `[structures]`, `[ish_qc]`,
`[ish_analysis]`, `[spatial_null]`, `[beyond]`, `[beyond_controls]`,
`[beyond_calibration]`, `[beyond_figures]`, `[top_genes]` and `[ish_figures]` of
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
| `beyond\` | 22 to 26 | `structures_used.csv` (with whether PSD95 is measured), `variance_partition.csv`, `calibration.csv`, `calibration_jackknife.csv`, `jackknife.csv`, `controls.csv`, `gene_space.csv`, `gene_space_calibration.csv`, `gene_space_summary.csv`, `variants.csv` (the main model, its other folds, its check rows and its wider structures: structures, terms, budget, density alone, left), `regression_table.csv`, `residual_by_structure.csv`, `replication.csv`, `leftover_genes.csv`, `leftover_sets.csv`, `numbers_for_the_caption.txt`, working figures |
| `tables\top_genes.csv` | 27 | per gene characterised (past the map's null, Gria1, Cacng8, the family): why it is there and its tier, its GO terms and gene sets, its rho, p, q and rank with the map, inside divisions, over the robustness variants, its reliability, its rho with each term of the main model, its prediction and PSD95, its rho and p with the leftover (q over all genes and within the family), its matched control, and what it takes of the leftover with both nulls |
| `tables\named_tests.csv`, `family_members.csv` | 27 | the tests named for the leftover, tier 1 and the two group tests of tier 2, on the leftover and on the map; every member of the family with its sources, tested or why not |
| `green_channel\sep_channel_check.csv` | 28 | per adult, each channel's range and correlations |
| `tables\april_headline.csv`, `april_anova.csv`, `april_groups.csv` | 29 | P9's genes then and now; the ANOVA under each choice; today's groups against the null |
| `tables\numbers_<step>.csv`, `numbers_for_the_text.csv` and `.txt` | 13 to 29 | the numbers of each step, and all of them |
| `cache\` | 14, 15 | Allen experiment lists, mygene records, `go-basic.obo` |

| figure | drawn by |
|---|---|
| `figures\00_overview.png`, `18_april_headline.png`, `README.md` | `run_ish_overview.py` |
| `figures\01_structures.png` | `run_structure_set.py` |
| `figures\02_genes.png` | `run_ish_gene_table.py` |
| `figures\03_beyond_budget.png`, `04_beyond_where.png`, `11_leftover_genes.png` | `run_beyond_figures.py` |
| `figures\05_one_comparison.png`, `06_spatial_null.png`, `07_gene_ranking.png`, `15_autofluorescence.png` | `run_ish_gene_ranking.py` |
| `figures\08_gene_sets.png`, `09_localisation.png` | `run_ish_gene_sets.py` |
| `figures\10_between_within.png`, `genes\<gene>.png` | `run_ish_divisions.py` (`--sheets`) |
| `figures\12_top_genes.png`, `13_cacng8_gria1.png`, `14_ampa_family.png`, `top_genes\<gene>.png` | `run_ish_top_genes.py` (`--sheets`) |
| `figures\16_robustness.png` | `run_ish_robustness.py` |
| `figures\17_green_channel.png` | `run_sep_channel_check.py` |
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
PLOS Biology 24, e3003637. Rouach N, Byrd K, Petralia RS, Elias GM, Adesnik H,
Tomita S, et al. (2005). TARP γ-8 controls hippocampal AMPA receptor number,
distribution and synaptic plasticity. Nature Neuroscience 8, 1525-1533. Schwenk
J, Harmel N, Brechet A, Zolles G, Berkefeld
H, Müller CS, et al. (2012). High-resolution proteomics unravel architecture
and molecular diversity of native AMPA receptor complexes. Neuron 74, 621-633.
Zhu F, Cizeron M, Qiu Z, Benavides-Piccione R, Kopanitsa MV, Skene NG, et al.
(2018). Architecture of the mouse brain synaptome. Neuron 99, 781-799.
