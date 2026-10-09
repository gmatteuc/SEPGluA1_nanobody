# The ISH analysis

What the adult nano map is, read against the Allen in situ hybridisation (ISH)
maps of the adult mouse brain. The document follows the one argument the ISH
line serves: part 1, the map is not fully explained by Gria1 expression and
synapse density; part 2, it is therefore something else, and the reading the data
support is the surface fraction of the receptor, with the genes that follow the
map as the corroboration. Each step has its method, figure, result and meaning;
then the controls, the limit, April's headline, how to rerun it and the files.

The line was built on 8 October 2026 from the decisions of the ISH discussion
([history/ISH_DISCUSSION.md](history/ISH_DISCUSSION.md), section 8), and run
again on 9 October with part 1's main model and the tests named for its leftover,
both fixed on 8 October before they were run (section 8.1 says what was known by
then).

- Code: `mapping/sepmap/structures.py`, `mapping/sepmap/adult/` and
  `mapping/sepmap/ish/`, run by steps 13 to 29 of the Python route
  ([mapping/README.md](../mapping/README.md)), on branch `post-ish`, not merged
  yet.
- Outputs: `<data>\adult_v2\ish_analysis\`. The numbers below are those of the
  run of 9 October 2026 on a full copy of the production inputs (the data root
  set by `SEP_DATA_ROOT`); the production data root holds the same once steps
  13 to 29 run there after the merge. Every number is in
  `tables\numbers_for_the_text.csv` there (a `.txt` beside it reads more
  easily), written by `run_ish_overview.py` from the numbers each step writes;
  each section names the step and the table its numbers come from.
- Figures: `figures\00_overview.png` to `16_april_headline.png`, PNG and EPS.
  A main figure answers one question in two to four panels, with one line
  under its title saying what to take from it; its detailed version carries the
  same number with an s (`05s_one_comparison_detail.png`; `03s1` and `03s2`
  where a main figure has two) and holds every panel and number of the analysis.
  `figures\README.md` walks through them in order, each with its question, what
  to look at and the takeaway, with this run's numbers. The figures are not
  versioned, so this document names each by its file.

## 1. The question and the argument

The adult nano map orders the brain's structures, and two halves of the cohort
order them the same way (rho 0.974 between half-cohort maps, over 126
grey-matter structures). A reader asks first whether that order is simply how
much GluA1 mRNA a structure makes, or how many synapses it has. Allen's ISH maps
answer that, because they measure both: the mRNA of Gria1, which encodes the
stained protein, and the mRNA of synaptic and postsynaptic-density genes.

The argument, in two connected parts:

1. **The map is not fully explained by Gria1 expression and synapse
   density.** Predicted from Gria1 mRNA, synapse density and the tissue's
   autofluorescence, a reproducible part of the map is left over: above what
   Allen-to-Allen mismatch alone leaves at its point value, with an interval
   that reaches just below zero, and no larger than what one Allen Gria1
   experiment leaves. This is the core claim (section 4).
2. **It is therefore something else, and the reading the data support is the
   surface fraction of the receptor**, set by trafficking and scaffolding. The
   genes that follow the map most, Cacng8 first (TARP γ-8, an AMPA receptor
   auxiliary subunit), are characterised as the corroboration, and the genes
   that follow what Gria1 and synapse density leave are tested in three tiers
   named on 8 October (section 5).

The limit, once: the surface fraction is an interpretation, not a measurement.
A total-GluA1 stain on the same brains would measure it; the green SEP channel
cannot (section 6.4).

Where it stands today:

- **Part 1 holds against the calibration floor at its point value only.**
  Gria1 and synapse density predict 71% of the map's reproducible pattern; 29%
  is left (15% to 42% over structures), and it replicates across mice (0.93),
  as the map's reliability alone implies (0.90). On the 113 structures of the
  calibration nano leaves 26%, 14 points more than what a map made only of
  Gria1 and synapse density leaves through Allen-to-Allen mismatch alone (the
  floor, 12%); the interval of that margin reaches just below zero (-2% to
  +30%), and the floor errs low (section 6.5). A measured synapse density,
  where it exists, predicts the map less well than the mRNA panel and leaves
  about as much.
- **It does not hold against the stricter benchmark.** A map that is one Allen
  Gria1 experiment leaves more (34%) than the nano map does (26% on the same
  structures). Which benchmark a claim uses is still to agree (section 8).
- **Part 2: the genes are consistent with the surface-fraction reading and do
  not single it out.** Cacng8 follows the map most of 451 genes (+0.807),
  inside divisions too, and Gria1 follows it (+0.646); the genes that follow
  the map most are maps much like Gria1 and synapse density. Cacng8's lead over
  Gria1 sits at the edge of the null of maps related to both alike (just past
  its 95% band, p 0.105 counting a lead either way), no gene set passes, and the
  localisation genes do no better than matched controls. Against the leftover,
  Cacng8 follows it (p 0.0002): a re-test, since its p against the
  near-identical leftover of 8 October was seen before it was named. Added to
  the model it takes 4.6 points of the reproducible map, more than maps alike
  to the model take (p 0.009, a null added after seeing the first) but not more
  than its plain surrogates (p 0.059). The AMPA receptor complex family follows
  the leftover beyond its surrogates only at the edge (p 0.047; 0.050 without
  Cacng8) and no more than postsynaptic genes of the same expression (p 0.55);
  no gene passes once every gene is corrected for.

What the comparison cannot do on its own: separate surface from total receptor
(both follow the postsynaptic side of synapses), read protein from mRNA (mRNA
sits in cell bodies, receptor on dendrites), or give a p value without a null
that keeps the brain's smooth gradients (cortex and hippocampus high, thalamus
and hypothalamus low). Each Allen experiment is one P56 mouse on a 200 um grid.

The figures, in the order of the argument (`figures\00_overview.png` puts the
question, the two parts with their key numbers, the limit and this map on one
page):

| part | main figure | detailed versions |
|---|---|---|
| the inputs | 01 which structures enter every comparison; 02 which genes, and how good one Allen map is | 02s |
| part 1 | 03 do Gria1 expression and synapse density explain the map; 04 where the map sits above or below what they predict | 03s1, 03s2 |
| part 2 | 05 what a gene's rho with the map is; 06 how large a rho unrelated smooth maps give; 07 which genes follow the map, and what kind of maps they are; 08 whether the map follows Cacng8 more closely than Gria1; 09 inside divisions or only between them; 10 whether the genes that set surface receptor follow the map better; 11 whether any gene follows what Gria1 and synapse density leave | 05s, 06s, 07s, 09s, 10s1, 10s2, 11s1, 11s2 |
| controls | 12 the label or the tissue and 13 the choices made, of part 2; 14 a measured synapse density, of part 1 | 12s, 13s, 14s |
| the limit | 15 the green channel | 15s |
| April | 16 what is left of April's headline | 16s |

## 2. Words used here

| word | meaning |
|---|---|
| structure | one atlas region (CA1, VPM). Every comparison runs across structures, one value per structure on each side |
| rho | Spearman correlation across structures: do two maps put the structures in the same order |
| division | a coarse part of the brain: Isocortex, OLF, HPF, CTXsp, STR, PAL, TH, HY, MB, P, MY, CB |
| `zref` | per brain, log2 of each structure's nano mean over the brain's isocortex mean, minus its median over the declared set, over its p90 - p10 spread there; the adult map is the mean of the ten adults' `zref` |
| the declared set | the structures every comparison uses: grey matter measured in all ten adults (section 3) |
| beyond | short for "beyond Gria1 expression and synapse density": how much of the map they predict, and what is left (part 1; the `beyond` folder and figures 03 and 04) |
| main model | part 1's model, fixed on 8 October 2026 before it was run: Gria1 + synapse density + autofluorescence (section 4) |
| check row | a variant of the main model that changes one thing of it, fixed with it and reported beside it (`[beyond.variants]` of `mapping/settings.toml`) |
| ceiling | the share of the map that two halves of the cohort reproduce (Spearman-Brown of their agreement) |
| leftover | the part of the map a fitted model does not predict |
| calibration floor | what the same model leaves of a map made only of Gria1 and synapse density: the leftover Allen-to-Allen mismatch produces alone. It does not hold the mismatch of Allen's P56 mice with these brains (strain, age, grid, registration), nor that of the genes measured once, so it errs low |
| benchmark | what the same model leaves of a map that is one Allen Gria1 experiment |
| jackknife | here, recomputing a number on subsamples that each leave out a fifth of the structures; the spread of those values, scaled for the subsample size, gives its interval |
| spatial null | random maps with the nano map's smoothness, the surrogates; a spatial p is the share of surrogates that correlate with a gene at least as strongly as the map does |
| surrogate | the map's ranks shuffled, smoothed over near structures and rescaled until its variogram matches the map's (Burt et al. 2020) |
| variogram | how unlike two structures are, on average, as a function of their distance: low for near structures in a smooth map, flat for a shuffled one |
| Freedman-Lane | a null for a residual: each surrogate is put through the same fit before it is compared, so it carries no more of the model's pattern than the residual does |
| tier | a level of the tests against the leftover, named on 8 October: Cacng8 alone (1), the AMPA receptor complex family (2), every other gene (3, exploratory) |
| re-test | a test of a result already seen: Cacng8's tier 1, whose p against the near-identical leftover of the first run of 8 October was seen before it was named |
| check row of a test | a reading of a named test that changes one thing of it after the test ran, to see what the result rests on (the family without Cacng8, controls outside the model's terms), reported beside it and never in its place |
| maps alike to the model | maps that keep a gene's fit on the main model and replace the rest with a surrogate of it, as large: a null for what a gene takes of the leftover |
| partial rho | rho once a third map (the subunit composite) has been regressed out of both sides, on ranks |
| BH, q | Benjamini-Hochberg correction for testing many genes or sets; q is the corrected p, and "past BH" means q < 0.05 |
| positive control | a difference that must exist, run through the same test: if the test cannot find it, a null result means nothing |
| P9 | the MATLAB comparison of April 2026, `adult_matlab/run_compare_with_allen_ish.m` |
| A1 to A9, S1, S5 | the decisions of the ISH discussion that a step carries out ([history/ISH_DISCUSSION.md](history/ISH_DISCUSSION.md), section 8, and [REFACTOR_COVERAGE.md](REFACTOR_COVERAGE.md)): A1 the declared structures, A2 section QC, A3 the robustness rows, A6 between or within divisions, A7 the spatial null, A8 autofluorescence through the same steps, A9 the gene table and its documentation; S1 all ten adults, S5 the Cacng8 - Gria1 gap named in advance |

## 3. The inputs

Every analysis reads the same inputs, fixed before any correlation. Numbers:
steps `structure_set`, `section_qc` and `gene_table` of
`numbers_for_the_text.csv`.

**The structures** (A1; `run_structure_set.py`, `figures\01_structures.png`).
A structure enters when it is grey matter and every adult measures it: of the
280 structures of the adult table, 235 are measured in all 10 adults and 204 of
those are grey matter (TH 45, HY 43, Isocortex 39, MB 25, STR 14, HPF 11, OLF
9, PAL 9, CTXsp 7, P 2; `tables\structure_set.csv`). Of the 45 structures the
rule leaves out for their adults, 22 are seen in one to four adults, mostly pons
and medulla, whose values do not replicate across mice, and 23 in five to nine;
all ten is the rule agreed on 30 September (S1).

**The adult map.** The mean of the ten adults' `zref`, with its zero and spread
taken over the declared set, so they depend only on the brain and the list,
not on which other brains are in a run. Against the 17-brain reference used
before, each adult's zero moves by 0.012 to 0.103, and the cohort map's order
does not change (rho 0.99995; `tables\zref_reference.csv`). The
autofluorescence map of the same sections is read the same way.

**The genes** (A9; `run_ish_gene_table.py`, `figures\02_genes.png` and
`02s_genes_detail.png`). One table, `tables\gene_table.csv`, holds P9's 100
genes and the 390 genes of the ontology panel, 451 genes, with a row per Allen
experiment: 757 experiments with the repair (other experiments of Chat, Tph2,
Olig2 and Calb2), 740 of them usable, each excluded one with its reason. All
451 genes have a profile, P9's 100 included. Each gene carries its section QC,
its reliability, its gene sets (from GO, release 2026-07-26) and P9's category
as a label only. `tables\gene_documentation.csv` is the table for a supplement.

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

## 4. Part 1: not fully explained by Gria1 expression and synapse density

`run_beyond_density.py` to `run_beyond_figures.py` (steps 22 to 26);
`figures\03_beyond.png`, `03s1_beyond_budget.png`, `03s2_beyond_controls.png`
and `04_beyond_where.png`.
Numbers: step `beyond` of `numbers_for_the_text.csv`, and the tables of
`beyond\` named below.

### 4.1 The main model

Fixed on 8 October 2026 before it was run (`[beyond]` of
`mapping/settings.toml`; section 8). Not before every number about it: the first
run that day had already given each term's share alone (Gria1 47%, synapse
density 58%, autofluorescence -2%) and the four-subunit model (27% left), and
choosing Gria1 alone moved the leftover by about a point, toward the claim. The
reason for the choice is the biology below, and the four-subunit model stays as
a check row.

```
nano map  ~  Gria1 + synapse density + autofluorescence
```

- *Abundance is Gria1 alone.* The mice carry SEP-GluA1: the stained protein is
  GluA1, which Gria1 alone encodes. Gria2 to Gria4 encode partner subunits the
  nanobody does not see, and their availability sets GluA1's assembly and
  trafficking: the surface side, not abundance. The four subunits as four
  terms, the model of the first run of 8 October, are a check row.
- *Synapse density* is the measured PSD95 punctum density (section 4.2) if it
  covers at least 80% of the structures of the fit; otherwise the mRNA panel,
  the composite of 11 synaptic marker genes (Syp, Syn1, Vamp2, Bsn, Syt1, Dlg4,
  Homer1, Shank2, Shank3, Nlgn1, Camk2a) and the first principal component of
  186 postsynaptic-density genes (`psd_pc1`). PSD95 covers 77 of the 126
  structures (61%), so the main model keeps the panel, on all 126, and PSD95 is
  a check row on the 77 it covers.
- *Autofluorescence* of the same sections, so a pattern the tissue makes
  cannot pass for one of the label.

### 4.2 Synapse density, measured

`run_synaptome.py` (step 21), `mapping/sepmap/adult/synaptome.py`;
`figures\14_synaptome.png` and `14s_synaptome_detail.png`. Numbers: step
`synaptome`; `synaptome\density.csv`, `coverage.csv`, `agreement.csv`.

**Why.** The density terms of the panel are Allen mRNA, which sits in cell
bodies: a presynaptic marker's mRNA marks where the neurons that make the
synapses are, not where their synapses are. Zhu et al. (2018) counted
excitatory synapses where they are.

**Source.** In a knock-in mouse with PSD95 and SAP102 tagged (Dlg4-eGFP,
Dlg3-mKO2), Zhu et al. detected every punctum in five coronal sections of one
adult male (about P80), sorted each into 37 subtypes by intensity, size and
shape, and measured the density of each subtype (puncta per unit area) in
regions of the Allen Reference Atlas. Hansen et al. share the table: 37
subtypes by 775 samples, a sample being one region of one hemisphere in one
section. It is downloaded from their repository, netneurolab/hansen_synaptome,
`data/synaptome/mouse_liu2018/` at the pinned commit
`0399525412b6f50cfdeb5904b96da7fa8e4b507c` (the files are the same at release
v1.0, archived as Zenodo doi 10.5281/zenodo.18201390), into
`<data>\reference\synaptome\`, each file checked against git's hash of it at that
commit; `fetch_log.txt` there records the source, the commit, the citation and
the hashes.

**What is read.** Subtypes 1 to 11 hold PSD95 alone, 12 to 18 SAP102 alone and
19 to 37 both (as Hansen et al.'s code indexes them). As shared, each subtype's
density is divided by its largest value over the 775 samples (each runs from 0,
the right locus coeruleus with no punctum at all, to 1), so absolute counts, and
their sum, are gone. The PSD95 density is the mean over the 30 subtypes whose puncta
hold PSD95, each subtype's map weighing alike; never a punctum's intensity or
size, since PSD95 per synapse is scaffolding, the surface side. Variants: PSD95
alone (Hansen et al.'s "PSD95 synapses"), SAP102 (26 subtypes) and every
punctum (37). They order the structures of the fit at 0.86, 0.68 and 0.94 with
the PSD95 density.

**Placing the samples.** By the Allen id the source gives each sample (its
acronym where it gives none), through the CCF 2017 ontology: a sample lies in
the structure that is its id or the id's nearest ancestor. A unit (a layer, or
the structure itself) is the mean of its samples over both hemispheres and
every section; a structure is the mean of its sampled units weighted by their
voxels in the 20 um CCF annotation, a unit that holds another sampled unit
weighing only the voxels outside it (PAG, sampled whole and in its nucleus ND).
Where a unit is not drawn in the annotation the units weigh alike (9 structures
of the fit: olfactory areas, the lateral entorhinal cortex and the subiculum,
whose layers the annotation does not draw); their covered share is then
unknown, and two of them are measured in part only, SUB in its dorsal part and
AON in three of its parts (AON1, AONm, AONpv). In COAp the posterolateral part
is read both whole and by its layers, each reading counting once. A region the
source gives only above several structures is never spread onto them, and a
structure with no sample stays missing. One thalamic sample is written 'lMD',
which Hansen et al. read as IMD, the intermediodorsal nucleus; it could equally
be the lateral part of MD. The source's id is followed, so MD stays missing; read
the other way, the fit would have one structure more measured (78), still far
below the rule. Of the 775 samples, 739 lie in 114 structures;
left out are 8 above the structures (the midbrain's and the medulla's motor
parts, the medial septal complex, the midbrain raphe), 12 that the 2017
annotation no longer draws (the layers of PTLp, now VISa and VISrl; the cochlear
granular lamina), 15 whose names are not CCF acronyms (Mop, ZID, ZIV, DCOmo,
CENT1) and the right locus coeruleus, with no punctum at all.

**The rule, and how it came out.** The measured density replaces the mRNA
panel in the main model only if it covers at least 80% of the structures of the
fit (101 of 126; `[beyond] min_psd95_coverage`). The threshold was Giulio's, set
on 8 October; an exploratory match of the samples to the structures ran nine
minutes before the rules were committed, but no choice of matching reaches 101
(spreading every region given above several structures would reach 87). It
covers 77 (61%) and 96 of the 204 declared structures. Missing from the fit: 10
thalamic nuclei (AM, AV, CL, MD, PR, PT, PoT, RE, SMT, VAL), 9 cortical areas
(AIp, FRP, SSp-ll, SSs, VISC, VISa, VISl, VISpm, VISrl), 7 hypothalamic, 6
midbrain (all under the midbrain's motor part), 5 pallidal, 4 striatal, 3
olfactory, 3 hippocampal and both pontine structures. Two are measured in less
than half their volume (SCs from its zonal layer, BMA from its posterior
part). So the main model keeps the mRNA panel, and PSD95 enters as a check row
on the 77 structures it covers (section 4.4).

**Agreement** (Spearman over the 77): with the marker mRNA composite 0.73, with
`psd_pc1` 0.80, with Gria1 0.42, with the nano map 0.45, with autofluorescence
0.29. The two hemispheres agree at 0.96 (75 structures), the one check of
reliability a single mouse allows; it says nothing about how a structure varies
from section to section, and 32 of the 77 are measured in at most two samples.

### 4.3 Method

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
  on each jackknife subsample, so the interval is that of the difference. The
  floor errs low in two ways: a gene measured by one experiment (Shank2, Shank3
  and Nlgn1 of the markers, and 25 of the 186 genes behind `psd_pc1`) is the
  same in both halves, so its own mismatch is not in it; and it holds only the
  mismatch of one Allen map with another, not that of Allen's P56 mice with these
  brains (strain, age, the 200 um grid against 20 um masks, registration). A
  small effect goes the other way: each half builds `psd_pc1` from the genes it
  measures everywhere (199 from half A, 221 from half B, 235 merged, on the 113
  structures). The PSD95 density, one mouse measured once, would enter like a
  marker measured once if it were in the main model.
- Seven controls try to break the result: a spatial gradient, structure size,
  single animals, naive against RWS, curvature, the choice of predictors (the
  components of every gene measured in all the structures, judged against its
  own calibration floor), the reading.
- The check rows (`[beyond.variants]`, fixed with the main model), each one
  change to it, on its own structures with its ceiling recomputed there: the
  density as PSD95 puncta, SAP102 puncta or every punctum, on the 77 structures
  where they are measured; the mRNA panel on those same 77, so that the two
  kinds of density meet on the same structures; PSD95 and the panel together,
  every density measure at once; the panel without its five presynaptic
  markers; Gria1 to Gria4 in place of Gria1. Beside them, the main model under
  other folds and on more structures. The rules of 8 October called PSD95 and
  the panel together the conservative bound; scored on held-out structures it is
  not one, since a term that adds nothing costs a little through overfitting.
- The rows that change only the density measure meet the same jackknife
  subsamples of the 77 structures (400, each leaving out 15), so their
  differences carry intervals of their own
  (`beyond\jackknife_density_rows.csv`); the interval of each row alone is far
  wider than the differences between them.

### 4.4 Result

`figures\03_beyond.png`: A, the map against the main model's held-out
prediction; B, the budget; C, the leftover beside the floor and the benchmark;
D, its replication. `03s1_beyond_budget.png` adds the map against each
predictor, the four-subunit check row and control F, every calibration draw and
every check row; `03s2_beyond_controls.png` draws the seven controls. Tables:
`beyond\variance_partition.csv`, `calibration.csv`, `calibration_jackknife.csv`,
`replication.csv`, `controls.csv`, `variants.csv`, `jackknife_density_rows.csv`.

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
  0.78: the floor. The nano map, read the same way, leaves 26% (13% to 40% over
  resampled structures), 14 points more than the floor at its point value; the
  interval of the difference runs from -2% to +30%. A map that is one Allen
  Gria1 experiment leaves 34% (35% when the truth comes from one half of its
  experiments, 34% from the other); nano minus that is -8% (-29% to +14%).
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
  less well alone than the mRNA panel does (22% against 45%; PSD95 minus the
  panel -23 points, 95% -41 to -5), and the model leaves about as much with it
  (46% against 39%; +7 points, -12 to +27); the two together leave as much as
  the panel alone (+0 points, -12 to +12) (`figures\14_synaptome.png` C). So the
  measured density does not shrink the leftover, and keeping the panel, the
  stronger density term alone, makes part 1 no easier to pass. The same model
  leaves 39% on the 77 and 29% on all 126; the 77 are not a random part of the
  126 (few hypothalamic, pallidal and midbrain structures), so the check rows on
  them are compared with each other, not with the main value. The partner
  subunits change little once density is in (27% left against 29%).
- The main model under other folds: one shuffling 28% (single shufflings 26%
  to 31%), ten folds 28%, leave one out 28%, and 28% on the 159 structures kept
  once the markers measured once (Shank2, Shank3, Nlgn1) are left out. With
  folds of spatial blocks, which keep a structure's neighbours out of its fit,
  38%; on the calibration's structures nano then leaves 30%, the floor 16% and
  the Gria1 map 37%.
- All 7 controls pass (`figures\03s2_beyond_controls.png`; `beyond\controls.csv`):
  a smooth gradient explains R² 0.215 of the leftover, which still replicates at
  0.916 with position in the model; its rho with structure volume is -0.065;
  single adults' leftovers agree at a median of +0.755 (lowest pair +0.564); the
  naive and RWS groups' at +0.877; bending to fifth powers buys 0.004 of held-out
  R² over cubes (0.705 to 0.709); and every reading of the map gives a leftover
  that replicates at 0.90 to 0.93. Control F: the components of the 311 genes
  measured in every structure (most often 21, picked inside each training fold)
  predict 81%; on its own calibration nano leaves 22% to 26%, against 4% to 9%
  for a map made of those genes' expression.

**Where it sits** (`figures\04_beyond_where.png`; `beyond\regression_table.csv`,
`residual_by_structure.csv`): above prediction in the medial geniculate (+56
ranks), both parts of the lateral geniculate, the lateral posterior nucleus, the
septofimbrial nucleus, the subthalamic nucleus and the subiculum; below it in
the nucleus of reuniens, the medial habenula, dorsal retrosplenial cortex, the
arcuate nucleus, piriform cortex and the submedial nucleus of the thalamus.

### 4.5 What it means

About 29% of the map's reproducible pattern is not predicted by Gria1
expression or synapse density, and it is reproducible across mice. On the same
structures it is above what Allen-to-Allen mismatch leaves at its point value
(26% against the floor's 12%), a margin whose interval reaches just below zero
(-2%), against a floor that errs low; with folds of spatial blocks the margin is
about the same (30% against 16%). Measured synapse density, where it exists,
leaves about as much.

What it does not mean. It is no larger than what a map that is one Allen Gria1
experiment leaves (34%): a map made of nothing but one experiment's Gria1, with
that experiment's own error, leaves as much. Its replication follows from the
map's reliability (0.897 expected), so it is not separate evidence. It is "not
predicted by Gria1 expression or synapse density", not "beyond gene
expression": the components of many panel genes predict 81% of the map. And a
leftover says what the predictors miss, not what it is; a claim about one
structure needs its own null.

## 5. Part 2: what else the map is, through the genes

If the map is the surface fraction, the genes that set how much receptor
reaches and stays at the membrane should follow it, and should follow what
Gria1 and synapse density leave. The section first says what one comparison is
and why it needs a spatial null (5.1, 5.2), then describes the genes that
follow the map most (5.3), Cacng8 against Gria1 (5.4) and at what scale they
follow it (5.5); then come the tests: kinds of genes fixed in advance (5.6) and
the three tiers named on 8 October against the leftover (5.7).

### 5.1 One comparison

`figures\05_one_comparison.png` and `05s_one_comparison_detail.png`. A gene's
map and the nano map are each reduced to one value per declared structure (for
a gene, the mean rank of its usable experiments), and the two orders correlated.
Cacng8 gives +0.807 on 172 structures, Gria1 +0.646 on 164, the astrocyte gene
Aqp4 +0.021 (step `overview`). Much of a whole-brain rho is cortex and
hippocampus against thalamus, which is why the null and the test inside
divisions follow.

### 5.2 The spatial null

`run_ish_spatial_null.py`, `mapping/sepmap/ish/spatial_null.py`;
`figures\06_spatial_null.png` and `06s_spatial_null_detail.png`. Numbers: step
`spatial_null`; `tables\variogram.csv`, `null_calibration.csv`.

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

### 5.3 The genes that follow the map

`run_ish_gene_ranking.py` and `run_ish_top_genes.py`,
`mapping/sepmap/ish/gene_ranking.py` and `top_genes.py`;
`figures\07_top_genes.png`, `07s_gene_ranking.png` and one sheet per gene,
`figures\top_genes\<gene>.png`. Numbers: steps `gene_ranking` and `top_genes`;
`tables\gene_ranking.csv`, `top_genes.csv`.

**Method.** Every gene's rho with the map against its own null (its rho with
each of the 10,000 surrogates), BH over all 451 genes and within P9's 100. The
genes past the map's null after BH over every gene, with Gria1 and Cacng8
always and every member of the AMPA receptor complex family (section 5.7), are
then each described on every axis of the line:

- with the map: rho, spatial p, q and rank; the mean rho inside divisions and
  its null (section 5.5); the range of its rho over the robustness variants,
  and of its rank over the variants that hold every gene (section 6.2);
- its Allen map: how many experiments, and how well they agree;
- whether it is a Gria1-like or density-like map: its rho with Gria1, the two
  density terms, autofluorescence and the main model's prediction on the 126
  structures of the fit, and with the PSD95 density where it is measured;
- with the leftover: its rho and spatial p, each surrogate through the same fit
  (section 5.7);
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

**Result.** 12 of 451 genes pass the spatial null after BH (q < 0.05), 6 of
P9's 100; before correction 128 of 451 and 29 of P9's are past their band. The
12, best first ("past" marks a within rho past its null after BH, "-" a gene
measured once; the share taken is of the reproducible map, with its p against
plain surrogates and against maps alike to the model):

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
| Cnih2 (two Allen experiments that disagree) | +0.637, 36 | +0.336 | 5 to 43 | 0.08 | +0.56, +0.63 | +0.002 (0.978) | -1.1% (0.785, 0.724) |
| Grip1 (in psd_pc1) | +0.615, 45 | +0.207 | 22 to 75 | 0.65 | +0.74, +0.69 | +0.095 (0.133) | -0.1% (0.588, 0.440) |
| Eps8 | +0.605, 54 | +0.286 | 8 to 79 | 0.55 | +0.53, +0.61 | +0.136 (0.057) | 1.6% (0.220, 0.141) |

- Four of the 12 are localisation genes and one is a subunit; one of the four,
  Cnih2, has two Allen experiments that disagree (reliability 0.08), so its
  place is weak.
- They are maps much like Gria1 and synapse density: rho 0.38 to 0.77 with
  Gria1 and 0.61 to 0.81 with the main model's prediction. That is why the main
  model predicts most of the map.
- Cacng8 is first of all genes, and first of P9's genes in every robustness
  variant (section 6.2; under Pearson on log2 it is 2nd of all genes). Gria1 is
  11th of P9's genes (q 0.0025 within them) and 32nd of all.
- 6 of the 12 follow the leftover before correction (Cacng8, Arpc5, Htr3a,
  Neurl1a, Mapk1, Add3); only Cacng8's p is a test, and a re-test (section
  5.7). Gria1 is a
  term of the model, so its rho with the leftover is near zero by construction
  and its p says nothing. Added to the main model, Htr3a alone takes more than
  95% of plain surrogates of it (5.9%, p 0.011); against maps alike to the
  model, Cacng8, Arpc5, Htr3a, Mapk1 and Add3 do (p 0.002 to 0.031).
- `07s_gene_ranking.png` shows P9's 100 genes one by one with their null bands
  and their rho with the autofluorescence map.

**What it means.** The genes that follow the map most describe it; they are
not separate evidence, since most of what they share with the map is what Gria1
and synapse density already predict. A few carry part of what the model leaves,
against maps alike to the model, the null added after the first was seen; against
plain surrogates only Htr3a does. Of these, Cacng8 is the one named for the
leftover (section 5.7). Htr3a, a serotonin receptor measured by one Allen
experiment, takes the most, an exploratory finding.

### 5.4 Cacng8 against Gria1

`figures\08_cacng8_gria1.png`; the gap in detail in `07s_gene_ranking.png` B.
Numbers: steps `gene_ranking` and `top_genes`; `tables\gap.csv`.

The map follows Cacng8 at +0.807 and Gria1 at +0.646, inside divisions +0.548
and +0.417; the two genes agree at +0.74 over the structures of the fit. The
Cacng8 - Gria1 gap, on the 164 structures both have, is +0.167, steady across
adults (+0.156 to +0.179 over resampled adults; Cacng8 leads in 10 of 10
adults). Against maps related to both genes alike it sits at the upper edge of
the null, just past its 95% band (-0.224 to +0.166): 2.3% of those maps give a
Cacng8 lead this large and 8.2% a Gria1 lead this large, so the test fixed in
advance, which counts a lead as large either way, gives p 0.105. Against maps
unrelated to both, a wider null and a conservative bound, p 0.348. It moves
with the Allen experiment that stands for each gene, from +0.08 to +0.211 over
the four pairings (p 0.064 to 0.385).

**What it means.** The map follows a TARP's pattern and Gria1's beyond a map
with the brain's smoothness. Cacng8 leads Gria1 in each of the 10 adults and in
every variant, by a margin at the edge of what maps following both genes alike
give: past the band one-sided, inside the test fixed in advance. These maps do
not settle which of the two the nano map follows more closely; whether the test
should count a lead one way only is a choice to make before any rerun (section
8), not after. Cacng8's own test is against the leftover (section 5.7).

### 5.5 Inside divisions

`run_ish_divisions.py`; `figures\09_between_within.png`,
`09s_between_within_detail.png` and the gene sheets `figures\genes\<gene>.png`
(Cacng8, Gria1, Grm5, Dlg2, Aqp4). Numbers: step `divisions`;
`tables\within_division.csv`.

**Method.** Per gene, rho with a map that knows only each structure's division
(its median), and the mean rho inside the divisions with at least 8 structures
(HPF, HY, Isocortex, MB, OLF, PAL, STR, TH), weighted by their structures, so
no contrast between divisions can enter it. Tested against the surrogates of
the whole map; a shuffle of structures inside divisions is shown beside it.

**Result.**

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

### 5.6 Kinds of genes

`run_ish_gene_sets.py`; `figures\10_gene_kinds.png`, `10s1_gene_sets.png` and
`10s2_localisation.png`. Numbers: step `gene_sets`; `tables\set_tests.csv`,
`contrasts.csv`, `localisation_summary.csv`, `localisation_power.csv`.

**Method.** Each set's median rho against the median of the same genes' rho
with every surrogate, so co-expressed genes stay together in the null; BH over
the sets with at least 5 genes. Two contrasts named in advance, both needed for
the map to count as postsynaptic-like: the three postsynaptic sets pooled
against the presynaptic set, and against glia. Then the localisation test of 5
October on the new inputs: each gene's partial rho with the map once the
subunit composite is removed from both; the 84 localisation genes against 84
other postsynaptic genes matched on expression, labels permuted; the same with
the four subunits removed as separate terms; a positive control (control genes
with a reproducible map against those without); the matched difference against
maps whose remainder, once the composite is removed, is a surrogate of the real
remainder; and the test's power: maps that carry the pattern setting the
localisation genes apart from their controls, in growing amounts, and how often
the test finds it.

**Result.**

| set | genes | median rho | spatial p | q |
|---|---|---|---|---|
| subunits | 4 | +0.568 | too few to test | |
| localisation | 84 | +0.375 | 0.054 | 0.135 |
| other postsynaptic | 196 | +0.407 | 0.049 | 0.135 |
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

### 5.7 Against the leftover: the three tiers named on 8 October

`run_beyond_density.py` (every gene against the leftover) and
`run_ish_top_genes.py` (the tests named for it); `figures\11_leftover.png`,
`11s1_leftover_genes.png` and `11s2_ampa_family.png`. Numbers: steps `beyond`
and `top_genes`; `beyond\leftover_genes.csv`, `leftover_sets.csv`,
`tables\named_tests.csv`, `family_members.csv`, `leftover_null_check.csv`.

**Method.** Every gene against the leftover of part 1's main model, against the
leftover's own surrogates, each surrogate put through the same fit first (the
leftover carries nothing of the model's columns, and a surrogate that did would
give a null far too wide). Three tiers were named for this leftover on 8
October, before the main model ran (`mapping/sepmap/ish/gene_sets.py`):

1. Cacng8 alone, its uncorrected spatial p the test. Its p against the leftover
   of the four-subunit model, 0.0024, was seen before it was named, and that
   leftover and this one agree at 0.90 over the same 126 structures (every
   gene's rho with the two agrees at 0.87; `ish_analysis_8oct_four_subunits\`
   on the development copy). So tier 1 re-tests a result already seen on a
   near-identical map: it shows the result holds under the main model, it does
   not confirm it.
2. The AMPA receptor complex family as a group: its median rho with the
   leftover against the same genes' median over the leftover's 10,000
   surrogates, and against 31 genes of the other postsynaptic set outside it,
   each the one closest in median energy (the localisation test's matching),
   labels permuted 20,000 times; one test each. Then gene by gene, BH within
   the family only. The family comes from sources outside this analysis: the
   constituents of native AMPA receptor complexes found by proteomics (Schwenk
   et al. 2012, Figure 1D and Table S2), the mouse genes annotated to
   GO:0032281 AMPA glutamate receptor complex (GO release 2026-08-05), and the
   partner subunits Gria2 to Gria4, less Gria1, the abundance term: 37 genes,
   31 with a usable map. Not tested: Olfm3, Prrt1, Prrt2 and Shisa8, which have
   no Allen experiment; Gsg1l and Rap2b, which have two each but are in neither
   panel of the gene table, and adding them would change the rules. Every
   gene's rho with the leftover of the four-subunit model was seen on 8 October;
   the family comes from the paper and GO, not from those numbers. Two of its
   members need saying: Dlg4 is also a marker of the model's density term, so
   its rho with the leftover is near zero by construction, as Gria1's would be
   (the rules took out Gria1 alone, so Dlg4 stays and pulls the median toward
   zero); and most of the matched controls are themselves model inputs (30 of
   the 31 are postsynaptic-density genes, 18 of them inside `psd_pc1`, which the
   model projects out of the leftover).
3. Every other gene, exploratory, BH over all of them, the model's own terms
   (Gria1, the markers) included as the rule set it. With 10,000 surrogates the
   smallest p is 0.0001, just under the line for the first gene (0.05 / 451 =
   0.00011): a gene can cross it only if no surrogate is as extreme. The gene
   sets of section 5.6 are read the same way.

Two check rows of tier 2 were added after the tests ran, to say what the group
result rests on, never in its place: the family without Cacng8 against the same
surrogates, and the family against controls matched afresh from the postsynaptic
genes outside the model's terms. And each surrogate of the leftover, once the
model is projected out of it, is rougher at short range than the leftover
itself (each map's Spearman with the mean of every structure's five nearest
neighbours: leftover 0.636, its surrogates 0.451), which narrows the null; the
tests are read again against Gaussian fields of exponential covariance with
ranges of 2, 6 and 10 mm, 5,000 each, projected the same way
(`tables\leftover_null_check.csv`).

The same group tests are read on the map itself, to describe the family, not as
tests of the leftover.

**Result.**

| test | on the leftover | on the map (a description) |
|---|---|---|
| tier 1, Cacng8 alone | +0.204, p 0.0002 (28th of 451 by rho) | +0.807, p ≤ 0.0001 |
| tier 2, the family's median against the surrogates | +0.064 (null 95% -0.063 to +0.063), p 0.047 | +0.414, p 0.044 |
| tier 2, the family against matched controls | +0.064 against +0.046: +0.018, p 0.55 (±0.062 would pass) | +0.414 against +0.365: +0.049, p 0.50 |
| tier 2, members past BH within the family | Cacng8 (q 0.006) | Cacng8 and Cnih2 (q 0.002) |
| tier 2 check, the family without Cacng8 | +0.063, p 0.050 | |
| tier 2 check, against controls outside the model | +0.064 against +0.034: +0.030, p 0.63 | |
| tier 3, genes past BH over all 451 | 0 (27 below p 0.05) | |

| null of the leftover | smoothness | Cacng8's p | the family's p | genes below p 0.05 |
|---|---|---|---|---|
| its surrogates, through the fit | 0.451 | 0.0002 | 0.047 | 27 |
| Gaussian fields, 2 mm | 0.491 | 0.0002 | 0.033 | 31 |
| Gaussian fields, 6 mm | 0.589 | 0.0002 | 0.042 | 19 |
| Gaussian fields, 10 mm | 0.611 | 0.0002 | 0.051 | 13 |

The leftover's own smoothness is 0.636.

- Beside Cacng8, the members of the highest rho with the leftover are Lrrtm4
  (+0.314, p 0.074), Vwc2l (+0.252), Shisa6 (+0.210) and Cpt1c (+0.193); after
  Cacng8, Abhd12 (p 0.013) and Dlg3 (p 0.036) have the smallest p, q 0.21 and
  0.33 within the family. The partner subunits sit near zero or below (Gria2
  -0.064, Gria3 +0.062, Gria4 -0.105).
- Over every gene the highest rho is Aldh1l1's, an astrocyte gene (+0.321, p
  0.20, q 0.80). The glia set has median +0.164 (p 0.043, q 0.19), localisation
  +0.048 (p 0.116). In the family, the two Allen experiments of Cnih2 and of
  Cacng4 disagree (reliability 0.08 and 0.20), so their place on the map and
  against the leftover is weak; Lrrtm4, the member of the highest rho, has a
  single experiment.
- Added to the main model, Cacng8 takes 4.6 points of the reproducible map
  (28.6% left, 24.0% with it): more than 95% of maps alike to the model take (p
  0.009, a null added after the first was seen), not more than 95% of its plain
  surrogates (p 0.059; section 5.3).

**What it means.** The gene named for the leftover, a TARP, follows what Gria1
and synapse density leave, beyond the leftover's own null and under every
smoother null; but this re-tests what was seen on 8 October on a near-identical
leftover, so it shows the result holds, it does not confirm it. It takes part of
the leftover against maps alike to the model, not against its plain surrogates.
The family follows the leftover more than maps of the leftover's smoothness do
only at the edge (p 0.033 to 0.051 over the nulls), and without Cacng8 its p is
0.050: the group result is not separate from tier 1. Other postsynaptic genes of
the same expression follow the leftover about as much, from inside the model's
terms or outside them, so what the family shares with the leftover is
postsynaptic, not particular to the AMPA receptor complex. No gene passes once
every gene is corrected for, and how many fall below p 0.05 depends on how smooth
the null is (13 to 31). The glia set followed the leftover of the first,
four-subunit model (q 0.04), a test not named in advance; against this one it
does not pass, and it stays an unplanned lead, not pursued.

### 5.8 What part 2 says

Consistent with the surface-fraction reading, and not singling it out. The map
follows receptor expression and not the tissue (section 6.1), and Cacng8, Dlg2
and Gria1 lead the ranking inside divisions too. But the genes that follow the
map most are maps much like Gria1 and synapse density, Cacng8's lead over Gria1
sits at the edge of the null of maps related to both alike, no gene set passes
its null, the postsynaptic criterion named in advance is met by half, and the
localisation genes predict the map no better than matched postsynaptic controls
(no advantage larger than about +0.08). Against the leftover, the TARP named for
it follows it (p 0.0002), the one result here that points the way the reading
does, though a re-test of what was seen on 8 October; it takes part of the
leftover only against the null added after seeing. The AMPA receptor complex
family named with it follows the leftover no more than postsynaptic genes of the
same expression. The genes neither contradict nor confirm part 2; only a
total-receptor measurement can.

## 6. Controls and the limit

### 6.1 The tissue: autofluorescence (A8)

`figures\12_autofluorescence.png` and `12s_autofluorescence_detail.png`.
Numbers: step `gene_ranking`; `tables\gene_ranking.csv`,
`gene_ranking_per_adult.csv`. The autofluorescence map of the same sections,
read as nano is and tested with its own surrogates.

- With Gria1 +0.176 (p 0.447), with Cacng8 +0.296 (p 0.292): neither passes.
  Adult by adult, nano's rho with Gria1 is 0.482 to 0.754, autofluorescence's
  -0.266 to +0.297, and nano is above autofluorescence for both genes in all 10
  adults.
- The tissue has a gene pattern of its own. After BH, 38 genes pass the
  autofluorescence null against 12 for nano (ribosomal proteins, Mapt, Eno2 and
  GABA-A subunits among them); before correction it is the other way round,
  94 past their band against 128. BH counts only the smallest p values, and
  autofluorescence has more very small ones. The two gene orders agree at only
  0.363.

So Gria1 and Cacng8 follow the label, not the tissue, in every adult; the
autofluorescence ranking is a different one, though not an empty one.

### 6.2 The choices made: robustness (A3)

`run_ish_robustness.py`; `figures\13_robustness.png` and
`13s_robustness_detail.png`; `tables\robustness_summary.csv`. Eleven variants of
the primary ranking (Spearman, `zref` with the declared reference, merged
profiles after QC, full means, the declared set): ten change one choice, the
last is the route of 5 October as it ran.

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
what moves (5th to 17th of P9's genes; 18th to 47th of all genes over the
variants that hold every gene): quote it with its null, not as a place. The gap
stays between +0.126 and +0.211, around the upper edge of the primary's null
band (figure 13 C draws the primary's band, not each variant's own): the
primary, the stored `zref` and no section QC sit at its edge, Pearson, P9's
single experiment and P9's divisions above it.

### 6.3 A measured synapse density

`figures\14_synaptome.png` and `14s_synaptome_detail.png`; section 4.2. The
measured PSD95 density covers 61% of the fit, below the 80% the rule asks; on
the 77 structures it covers it predicts the map less well alone than the mRNA
panel (22% against 45%; -23 points, 95% -41 to -5), and the model leaves about
as much with it (46% against 39%; +7 points, -12 to +27). It checks part 1 from
outside the ISH data: the panel the main model keeps is the stronger density
term alone, so keeping it makes part 1 no easier to pass.

### 6.4 The limit: the green channel (analysis 5)

`run_sep_channel_check.py`; `figures\15_green_channel.png` and
`15s_green_channel_detail.png`; `green_channel\sep_channel_check.csv`. Numbers:
step `green_channel`.

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
(`ish/arms.py`, `arms_vs_genes.png`) and the ratio table behind them
(`adult/arms.py`) are retired to `archive/`: their premise failed.

### 6.5 Limits

- Each Allen experiment is one P56 mouse on a 200 um grid; the adults' age is
  not recorded.
- mRNA sits in cell bodies, receptor on dendrites, which can lie in another
  structure (cortical layer 1).
- The gene panel is 451 genes chosen from GO and by hand, not a genome-wide
  background.
- The spatial null rests on structure centroids and a short-range variogram.
- The calibration floor errs low: genes measured once (Shank2, Shank3 and
  Nlgn1 of the markers, 25 of the genes behind `psd_pc1`) are the same in both
  halves, and it holds the mismatch of one Allen map with another, not that of
  Allen's P56 mice with these brains (strain, age, the 200 um grid against 20 um
  masks, registration). Part 1's leftover is above it at its point value only,
  and no larger than what one Allen Gria1 experiment leaves.
- The leftover's surrogates, once the model is projected out of them, are
  rougher than the leftover; the named tests hold against smoother nulls, the
  family's stays at the edge, and per-gene counts depend on the null.
- The measured synapse density is one mouse, in 77 of the 126 structures of
  the fit, as shared scaled per subtype: it orders the structures, it does not
  count synapses across them. Its hemispheres agree, but many structures are
  measured in one or two samples, so how a structure varies across sections is
  not known; SCs and BMA are measured in less than half their volume, SUB and AON
  in part.
- The exceptions list (true absences) and the sections kept at a step are a
  human call, still proposed.
- No total-receptor channel, no knockout or no-primary control, no DAPI
  control through the chain.

## 7. What changed from April

`run_ish_overview.py`; `figures\16_april_headline.png` and
`16s_april_headline_detail.png`. Numbers: step `overview`;
`tables\april_anova.csv`, `april_groups.csv`.

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
difference between kinds of genes does not survive the null, and section 5.6
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
| controls | none | autofluorescence through the same chain, robustness, positive controls, a power check, a measured synapse density |

## 8. Fixed in advance, and what is left to settle

### 8.1 Fixed in advance, and what was known then

- **The rebuild of part 1, its rules fixed first (8 October).** Written down
  and committed before it was run (what was known by then: each term's share
  alone, the four-subunit model, and every gene's rho with its leftover, from the
  first run that day):
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
    measure matters. The rules called PSD95 and the panel together the
    conservative bound, since one measured map in place of the panel's two terms
    gives the model fewer terms; scored on held-out structures it is not a bound
    (section 4.3). The calibration runs the same main model; PSD95 is one
    measurement, the same in both halves, like the markers measured once, so the
    floor errs low;
  - against the new leftover, Cacng8 is the one gene named, and the AMPA
    receptor complex family (Schwenk et al. 2012 and GO:0032281, with Gria2 to
    Gria4, without Gria1) is tested as a group, then gene by gene within it
    (`mapping/sepmap/ish/gene_sets.py`). Cacng8's p against the leftover of the
    four-subunit model was seen before it was named, and the two leftovers agree
    at 0.90, so tier 1 is a re-test. `gene_sets.py` is left as committed, so the
    rules can be read as they were; what the tests rest on, found after they ran,
    is in `top_genes.py` and section 5.7.

  Run on 9 October: PSD95 covers 77 of the 126 structures, so the main model
  keeps the mRNA panel and PSD95 is a check row (section 4). The outputs of the
  first run of 8 October, with the four subunits, are kept on the development
  copy in `adult_v2\ish_analysis_8oct_four_subunits\`. The tests named for the
  leftover ran the same day (section 5.7): Cacng8 follows it (p 0.0002); the
  family passes against the surrogates only at the edge (p 0.047, 0.050 without
  Cacng8), not against matched controls (p 0.55).

### 8.2 To settle

- **What a gene takes of the leftover.** Two nulls (section 5.3): plain
  surrogates of the gene, the null this step was first given and a wide one,
  and maps that relate to the model as the gene does, added after the first's
  numbers were seen. Cacng8 passes the second (p 0.009), not the first (p
  0.059). Which one a claim quotes is to agree; the first is the bound.
- **The gap's test, one way or both.** The Cacng8 - Gria1 gap is tested
  against maps that follow both genes alike, counting a lead as large either
  way (p 0.105), as fixed before; that null is skewed, and one way only 2.3% of
  its maps give a Cacng8 lead as large. Whether the test should count one way
  is Giulio's to decide before any rerun, not after; until then the two-sided p
  stands and the one-sided share is described.
- **The exceptions list.** Glra1 section 61 is kept as true absence, status
  "proposed" in `mapping/ish_section_exceptions.csv`, and 16 sections are kept
  at a step in expression by the rule. To review on `figures\qc\00_flagged.png`
  (red, hatched, dots) and the gene's sheets; the numbers are final once
  reviewed.
- **Part 1's benchmark.** The leftover is above the floor at its point value
  only (the interval of the difference reaches -2%, and the floor errs low) and
  below what one Allen Gria1 experiment leaves. Which benchmark a claim uses is
  to agree with Sami El-Boustani before the figure is shown; a one-sided interval
  of the margin would pass, but that choice was not fixed in advance, so it is
  not taken now.
- **The glia lead.** The glia set followed the leftover of the four-subunit
  model (q 0.04), a test not named in advance; against the main model's
  leftover it does not pass (q 0.19). An unplanned lead, not pursued now
  (figure 11s1 says so in one line).
- **The presynaptic set.** GO annotates Syn1, Vamp2, Stx1a, Bsn, Syt1, Cplx1,
  Camk2a and Slc17a7 to both sides of the synapse, so the "not both" rule keeps
  16 presynaptic genes, mostly cell-type and peptide markers. The sets stay as
  committed; the pre- and postsynaptic genes are drawn as context in figure
  10s1.
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

| step | script | needs | figures | time (development copy, 9 October) |
|---|---|---|---|---|
| 11, 12 | `run_panel_build`, `run_panel_fetch` | the network, once; not rerun for this analysis | | |
| 13 | `run_structure_set` | the per-brain files; `--recompute` redoes the per-adult channel table | 01 | 10 s (about 10 minutes with `--recompute`) |
| 14 | `run_ish_section_qc` | the Allen API for the repair's experiment lists, cached in `cache\`; `--offline` stops instead; `--sheets` draws the missing QC sheets | QC sheets | 20 s (about 15 minutes with `--sheets` from none) |
| 15 | `run_ish_gene_table` | mygene.info and the GO ontology, cached in `cache\`; `--offline` stops instead | 02, 02s | 5 minutes |
| 16 | `run_ish_spatial_null` | the surrogates are cached; `--recompute` draws them again | | 14 minutes, the calibration |
| 17 | `run_ish_gene_ranking` | the surrogates | 05, 05s, 06, 06s, 07s, 12, 12s | 2 minutes |
| 18 | `run_ish_robustness` | | 13, 13s | 1 minute |
| 19 | `run_ish_divisions --sheets` | | 09, 09s, gene sheets | 3 minutes |
| 20 | `run_ish_gene_sets` | | 10, 10s1, 10s2 | 1 minute |
| 21 | `run_synaptome` | the synaptome of Zhu et al. 2018, downloaded once into `<data>\reference\synaptome\`; `--offline` stops instead | 14s | 10 s |
| 22 to 25 | `run_beyond_density` to `run_beyond_regression` | | | under 1 minute each; the calibration 2 minutes |
| 26 | `run_beyond_figures` | the tables of steps 21 to 25 | 03, 03s1, 03s2, 04, 11s1, 14 | 1 minute |
| 27 | `run_ish_top_genes --sheets` | the tables of steps 15 to 22 | 07, 08, 11, 11s2, top-gene sheets | 8 minutes, the nulls of the share taken |
| 28 | `run_sep_channel_check` | | 15, 15s | under 1 minute |
| 29 | `run_ish_overview` | every step's numbers | 00, 16, 16s, `figures\README.md` | 10 s |

A full run of steps 13 to 29 takes about 40 minutes, and a second run gives
the same tables.

Every run prints the data root and the settings in force. `SEP_DATA_ROOT`
moves the data root, for a copy. The settings are `[structures]`, `[ish]`,
`[ish_qc]`, `[ish_analysis]`, `[ish_panel_test]` (the label permutations, and the
structures and reliability a gene needs, which `gene_sets` and `top_genes`
read), `[spatial_null]`, `[beyond]`, `[beyond_controls]`, `[beyond_calibration]`,
`[beyond_figures]`, `[beyond_regression]`, `[top_genes]` and `[ish_figures]` of
`mapping/settings.toml`. The tests: `cd mapping`,
`..\tools\venv_dev\Scripts\python -m pytest tests`
([mapping/tests/README.md](../mapping/tests/README.md)).

The outputs of 5 October in `adult_v2\ish\`, `arms\`, `beyond\` and `panel\`
are frozen; nothing here writes there but steps 11 and 12 (`panel\`). The run
scripts that still write into them are out of the run order
([mapping/README.md](../mapping/README.md) lists them).

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
| `synaptome\samples.csv`, `density.csv`, `coverage.csv`, `agreement.csv` | 21 | per sample of the synaptome, its ids, densities and structure or why none; per structure of the adult table, measured or why not, its units, weights and densities; per division, how many declared and fitted structures are measured; each density's Spearman with the mRNA density terms, Gria1 and the maps |
| `beyond\` | 22 to 26 | `structures_used.csv` (with whether PSD95 is measured), `variance_partition.csv`, `calibration.csv`, `calibration_jackknife.csv`, `jackknife.csv`, `jackknife_density_rows.csv` (the density check rows on the same subsamples), `controls.csv`, `gene_space.csv`, `gene_space_calibration.csv`, `gene_space_summary.csv`, `variants.csv` (the main model, its other folds, its check rows and its wider structures: structures, terms, budget, density alone, left), `regression_table.csv`, `residual_by_structure.csv`, `replication.csv`, `leftover_genes.csv`, `leftover_sets.csv`, `numbers_for_the_caption.txt`, working figures |
| `tables\top_genes.csv` | 27 | per gene characterised (past the map's null, Gria1, Cacng8, the family): why it is there and its tier, its GO terms and gene sets, its rho, p, q and rank with the map, inside divisions, over the robustness variants, its reliability, its rho with each term of the main model, its prediction and PSD95, its rho and p with the leftover (q over all genes and within the family), its matched control, and what it takes of the leftover with both nulls |
| `tables\named_tests.csv`, `family_members.csv`, `leftover_null_check.csv` | 27 | the tests named for the leftover, tier 1 and the two group tests of tier 2, on the leftover and on the map, and the two check rows of tier 2; every member of the family with its sources, tested or why not; the named tests and the count of genes below p 0.05 against smoother nulls, with each null's smoothness |
| `green_channel\sep_channel_check.csv` | 28 | per adult, each channel's range and correlations |
| `tables\april_headline.csv`, `april_anova.csv`, `april_groups.csv` | 29 | P9's genes then and now; the ANOVA under each choice; today's groups against the null |
| `tables\numbers_<step>.csv`, `numbers_for_the_text.csv` and `.txt` | 13 to 29 | the numbers of each step, and all of them |
| `cache\` | 14, 15 | Allen experiment lists, mygene records, `go-basic.obo` |

| figure | drawn by |
|---|---|
| `figures\00_overview.png`, `16_april_headline.png`, `16s_april_headline_detail.png`, `README.md` | `run_ish_overview.py` |
| `figures\01_structures.png` | `run_structure_set.py` |
| `figures\02_genes.png`, `02s_genes_detail.png` | `run_ish_gene_table.py` |
| `figures\03_beyond.png`, `03s1_beyond_budget.png`, `03s2_beyond_controls.png`, `04_beyond_where.png`, `11s1_leftover_genes.png`, `14_synaptome.png` | `run_beyond_figures.py` |
| `figures\05_one_comparison.png`, `05s_one_comparison_detail.png`, `06_spatial_null.png`, `06s_spatial_null_detail.png`, `07s_gene_ranking.png`, `12_autofluorescence.png`, `12s_autofluorescence_detail.png` | `run_ish_gene_ranking.py` |
| `figures\07_top_genes.png`, `08_cacng8_gria1.png`, `11_leftover.png`, `11s2_ampa_family.png`, `top_genes\<gene>.png` | `run_ish_top_genes.py` (`--sheets`) |
| `figures\09_between_within.png`, `09s_between_within_detail.png`, `genes\<gene>.png` | `run_ish_divisions.py` (`--sheets`) |
| `figures\10_gene_kinds.png`, `10s1_gene_sets.png`, `10s2_localisation.png` | `run_ish_gene_sets.py` |
| `figures\13_robustness.png`, `13s_robustness_detail.png` | `run_ish_robustness.py` |
| `figures\14s_synaptome_detail.png` | `run_synaptome.py` |
| `figures\15_green_channel.png`, `15s_green_channel_detail.png` | `run_sep_channel_check.py` |
| `figures\qc\00_flagged.png`, `qc\<gene>_<experiment>.png` | `run_ish_section_qc.py` (`--sheets`) |

Reference data, under `<data>\reference\`, each folder with a `fetch_log.txt`
(source, version, checksums, what was read):

| folder | what | source |
|---|---|---|
| `synaptome\` | PSD95 and SAP102 punctum density, 37 subtypes in 775 samples (a region of one hemisphere in one section), each subtype divided by its largest value as shared (Zhu et al. 2018); written by `run_synaptome.py` | github.com/netneurolab/hansen_synaptome, `data/synaptome/mouse_liu2018/` at commit `0399525412b6f50cfdeb5904b96da7fa8e4b507c`; archived as Zenodo doi 10.5281/zenodo.18201390 |
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
J, Harmel N, Brechet A, Zolles G, Berkefeld H, Müller CS, et al. (2012).
High-resolution proteomics unravel architecture and molecular diversity of
native AMPA receptor complexes. Neuron 74, 621-633. Zhu F, Cizeron M, Qiu Z,
Benavides-Piccione R, Kopanitsa MV, Skene NG, et al. (2018). Architecture of the
mouse brain synaptome. Neuron 99, 781-799.
