# The ISH analysis

What the adult nano map is, read against the Allen in situ hybridisation (ISH)
maps of the adult mouse brain. The document follows the one argument the ISH
line serves: part 1, the map is not fully explained by Gria1 expression and
synapse density; part 2, it is therefore something else, and the reading the data
support is the surface fraction of the receptor, with the genes that follow the
map as the corroboration. Each step has its method, figure, result and meaning;
then the controls, the limit, April's headline, how to rerun it and the files.

The line was built on 8 October 2026 from the decisions of the ISH discussion
([history/ISH_DISCUSSION.md](history/ISH_DISCUSSION.md), section 8), run on the
morning of 9 October with part 1's first model and the tests named for its
leftover, both fixed on 8 October before they were run, and run again the same day
with the second version of part 1, which Giulio decided after seeing the first
(section 8.1 says what was known at each step); three reviews of the second
version that evening changed how it is reported, not what it computes.

- Code: `mapping/sepmap/structures.py`, `mapping/sepmap/adult/` and
  `mapping/sepmap/ish/`, run by steps 13 to 30 of the Python route
  ([mapping/README.md](../mapping/README.md)), on branch `post-ish`, not merged
  yet.
- Outputs: `<data>\adult_v2\ish_analysis\`. The numbers below are those of the
  run of 9 October 2026 with part 1's second version, on a full copy of the
  production inputs (the data root set by `SEP_DATA_ROOT`); the production data
  root holds the same once steps 13 to 30 run there after the merge. Every number is in
  `tables\numbers_for_the_text.csv` there (a `.txt` beside it reads more
  easily), written by `run_ish_overview.py` from the numbers each step writes;
  each section names the step and the table its numbers come from.
- Figures: `figures\00_overview.png` to `15_april_headline.png`, PNG and EPS.
  A main figure answers one question in two to four panels, with one line
  under its title saying what to take from it; its detailed version carries the
  same number with an s (`05s_one_comparison_detail.png`; `03s1` to `03s4`
  where a main figure has several) and holds every panel and number of the
  analysis.
  `figures\README.md` walks through them in order, each with its question, what
  to look at and the takeaway, with this run's numbers. The figures are not
  versioned, so this document names each by its file.

## 1. The question and the argument

The adult nano map orders the brain's structures, and two halves of the cohort
order them the same way (rho 0.972 between half-cohort maps, over the 163
grey-matter structures of part 1). A reader asks first whether that order is
simply how much GluA1 mRNA a structure makes, or how many synapses it has.
Allen's ISH maps answer that, because they measure both: the mRNA of Gria1, which
encodes the stained protein, and the mRNA of postsynaptic genes, chosen by a
measured synapse density.

The argument, in two connected parts:

1. **The map is not fully explained by Gria1 expression and synapse
   density.** Predicted from Gria1 mRNA and synapse density, the mRNA of three
   postsynaptic genes chosen without the map by their agreement with a measured
   synapse density, a reproducible part of the map is left over, well above what
   Allen-to-Allen mismatch alone leaves. This is the core claim (section 4).
2. **It is therefore something else, and the reading the data support is the
   surface fraction of the receptor**, set by trafficking and scaffolding. The
   genes that follow the map most, Cacng8 first (TARP γ-8, an AMPA receptor
   auxiliary subunit), are characterised as the corroboration, and the genes
   that follow what Gria1 and synapse density leave are tested in three tiers
   named on 8 October (section 5).

The limit, once: the surface fraction is an interpretation, not a measurement.
A total-GluA1 stain on the same brains would measure it; the green SEP channel
cannot (section 6.4).

Where it stands today (the second version of part 1, 9 October; section 8.1):

- **Part 1 holds against the calibration floor.** On 163 structures, Gria1 and
  synapse density, two straight terms, predict 56% of the map's reproducible
  pattern (Gria1 only 10%, shared 31%, density only 15%); 44% is left (35% to
  54% over structures; 26% to 62% when whole spatial blocks of neighbours are
  left out), and it replicates across mice (0.93), as the map's reliability
  alone implies (0.94). On the 131 structures of the calibration nano leaves
  42%, 27 points more than a map made only of Gria1 and synapse density leaves
  through Allen-to-Allen mismatch (the floor, 15%; 95% +15 to +40, over spatial
  blocks +10 to +45). The floor errs both ways, by less than that margin: it
  lacks the mismatch of Allen's P56 mice with these brains, and it holds the
  disagreement of two Allen halves where nano meets one. Every check row stays
  above its own floor at its point value, 9 of 10 over spatial blocks; the
  lowest margins are those of Gria1 to Gria4 as abundance (28% left, +13 points)
  and of the model curved (38% left, +20). Control E does not pass: part of the
  straight leftover is curvature. Synapse density is three postsynaptic genes
  chosen without the map; the rule that chose them, re-run on random halves,
  follows the measured PSD95 punctum density at 0.79 on the other half, mostly
  between divisions (inside them 0.08). The leftover is in part a contrast
  between divisions (hippocampus, striatum and olfactory areas above prediction,
  isocortex and hypothalamus below). This is part 1's second version, chosen
  after seeing the first; the first's numbers are kept in section 4.6.
- **Part 2: the genes are consistent with the surface-fraction reading and do
  not single it out.** Cacng8 follows the map most of 451 genes (+0.807), inside
  divisions too, and Gria1 follows it (+0.646); the genes that follow the map
  most are maps much like Gria1 and synapse density. The Cacng8 - Gria1 gap
  (+0.167) is inside the two-sided test fixed in advance (maps related to both
  alike, p 0.105), no gene set passes, and the localisation genes do no better
  than matched controls. Against the leftover, Cacng8 follows it (p 0.0066), a
  re-test of what two earlier leftovers showed; added to the model it takes 9.8
  points of the reproducible map, more than its plain surrogates take (p 0.001,
  the null fixed first; maps alike to the model, the secondary line, p 0.003).
  The AMPA receptor complex family as a group does not follow the leftover (p
  0.92), and the three members past BH within it run against it (Gria4, Olfm2,
  Olfm1). 53 genes pass once every gene is corrected for (45 to 59 under
  smoother nulls), 50 of them with a negative rho: nano sits below prediction
  where pan-neuronal, vesicle and ribosomal mRNA is high, and all 53 follow the
  autofluorescence map too, so neuronal mRNA in general and the tissue are both
  readings of it; an exploratory finding.

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
| part 1 | 03 do Gria1 expression and synapse density explain the map; 04 where the map sits above or below what they predict | 03s1 (the check rows), 03s2 (the controls), 03s3 (the choice of the density genes), 03s4 (the measured synapse density) |
| part 2 | 05 what a gene's rho with the map is; 06 how large a rho unrelated smooth maps give; 07 which genes follow the map, and what kind of maps they are; 08 whether the map follows Cacng8 more closely than Gria1; 09 inside divisions or only between them; 10 whether the genes that set surface receptor follow the map better; 11 whether any gene follows what Gria1 and synapse density leave | 05s, 06s, 07s, 09s, 10s1, 10s2, 11s1, 11s2 |
| controls of part 2 | 12 the label or the tissue; 13 the choices made | 12s, 13s |
| the limit | 14 the green channel | 14s |
| April | 15 what is left of April's headline | 15s |

## 2. Words used here

| word | meaning |
|---|---|
| structure | one atlas region (CA1, VPM). Every comparison runs across structures, one value per structure on each side |
| rho | Spearman correlation across structures: do two maps put the structures in the same order |
| division | a coarse part of the brain: Isocortex, OLF, HPF, CTXsp, STR, PAL, TH, HY, MB, P, MY, CB |
| `zref` | per brain, log2 of each structure's nano mean over the brain's isocortex mean, minus its median over the declared set, over its p90 - p10 spread there; the adult map is the mean of the ten adults' `zref` |
| the declared set | the structures every comparison uses: grey matter measured in all ten adults (section 3) |
| beyond | short for "beyond Gria1 expression and synapse density": how much of the map they predict, and what is left (part 1; the `beyond` folder and figures 03 and 04) |
| main model | part 1's model, the second version, decided on 9 October 2026: Gria1 + synapse density, two straight terms on ranks (section 4) |
| density genes | Rock2, Cap2 and Slc8a2: the postsynaptic genes whose mean rank is synapse density, chosen without the map by their agreement with the measured PSD95 punctum density (section 4.2) |
| partition | the reproducible map in four parts from what Gria1 alone, density alone and both predict held out: Gria1 only, shared, density only and left (section 4.3) |
| check row | a variant of the main model that changes one thing of it, reported beside it with a floor of its own (`[beyond.checks]` of `mapping/settings.toml`) |
| ceiling | the share of the map that two halves of the cohort reproduce (Spearman-Brown of their agreement) |
| leftover | the part of the map a fitted model does not predict |
| calibration floor | what the same model leaves of a map made only of Gria1 and synapse density, predicted from other Allen experiments of the same genes: the leftover Allen-to-Allen mismatch produces alone. It errs both ways: it lacks the mismatch of Allen's P56 mice with these brains (strain, age, grid, registration), and it holds two Allen halves disagreeing where nano meets one half (section 4.3) |
| jackknife | here, recomputing a number on subsamples that each leave out a fifth of the structures, at random; the spread of those values, scaled for the subsample size, gives its interval. Over spatial blocks: each subsample leaves out one of 20 blocks of neighbouring structures, so neighbours, which are not independent, leave together; that interval is wider |
| spatial null | random maps with the nano map's smoothness, the surrogates; a spatial p is the share of surrogates that correlate with a gene at least as strongly as the map does |
| surrogate | the map's ranks shuffled, smoothed over near structures and rescaled until its variogram matches the map's (Burt et al. 2020) |
| variogram | how unlike two structures are, on average, as a function of their distance: low for near structures in a smooth map, flat for a shuffled one |
| Freedman-Lane | a null for a residual: each surrogate is put through the same fit before it is compared, so it carries no more of the model's pattern than the residual does |
| tier | a level of the tests against the leftover, named on 8 October: Cacng8 alone (1), the AMPA receptor complex family (2), every other gene (3, exploratory) |
| re-test | a test of a result already seen: Cacng8's tier 1, whose p against two earlier leftovers (the four-subunit model of 8 October, and the first version of part 1) was seen before |
| check row of a test | a reading of a named test that changes one thing of it after the test ran, to see what the result rests on (the family without Cacng8, controls outside the model's terms), reported beside it and never in its place |
| maps alike to the model | maps that keep a gene's fit on the main model and replace the rest with a surrogate of it, as large: a null for what a gene takes of the leftover |
| partial rho | rho once a third map (the subunit composite) has been regressed out of both sides, on ranks |
| BH, q | Benjamini-Hochberg correction for testing many genes or sets; q is the corrected p, and "past BH" means q < 0.05 |
| positive control | a difference that must exist, run through the same test: if the test cannot find it, a null result means nothing |
| P9 | the MATLAB comparison of April 2026, `adult_matlab/run_compare_with_allen_ish.m` |
| A1 to A9, S1, S5 | the decisions of the ISH discussion that a step carries out ([history/ISH_DISCUSSION.md](history/ISH_DISCUSSION.md), section 8, and [REFACTOR_COVERAGE.md](REFACTOR_COVERAGE.md)): A1 the declared structures, A2 section QC, A3 the robustness rows, A6 between or within divisions, A7 the spatial null, A8 autofluorescence through the same steps, A9 the gene table and its documentation; S1 all ten adults, S5 which ISH numbers to quote |

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
autofluorescence map of the same sections is read the same way. For a check row
of part 1, nano is also measured on the Allen 200 um grid, each 200 um voxel the
mean of its tissue and labelled as the ISH grids are, with zref taken the same way
(`tables\adult_allen_grid.csv`).

**The genes** (A9; `run_ish_gene_table.py`, `figures\02_genes.png` and
`02s_genes_detail.png`). One table, `tables\gene_table.csv`, holds P9's 100
genes and the 390 genes of the ontology panel, 451 genes, with a row per Allen
experiment: 757 experiments with the repair (other experiments of Chat, Tph2,
Olig2 and Calb2), 740 of them usable, each excluded one with its reason. All
451 genes have a profile, P9's 100 included. Each gene carries its section QC,
its reliability, its gene sets (GO annotations from the records of mygene.info,
fetched on 5 and 9 October and cached, with GO's ancestry from `go-basic.obo` of
data-version 2026-07-26) and P9's category as a label only. `tables\gene_documentation.csv` is the table for a supplement.

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

`run_density_markers.py` (step 22) and `run_beyond_density.py` to
`run_beyond_figures.py` (steps 23 to 27); `figures\03_beyond.png`,
`03s1_beyond_budget.png` (the check rows), `03s2_beyond_controls.png` (the
controls), `03s3_density_markers.png` (the choice of the density genes),
`03s4_synaptome_detail.png` (the measured synapse density) and
`04_beyond_where.png`. Numbers: steps `density_markers` and `beyond` of
`numbers_for_the_text.csv`, and the tables of `density_markers\` and `beyond\`
named below.

This is the second version of part 1. Giulio decided it on 9 October 2026 after
seeing every number of the first version, which had run that morning, for
reasons about the measurement; the genes of its density term were committed
before any nano number of it was computed. The first version's numbers are kept
in section 4.6 and its outputs in `adult_v2\ish_analysis_9oct_v1\` on the
development copy; section 8.1 says what was known at each step. Three reviews of
the second version, the same evening, changed how it is reported, not what it
computes: the floor's error is stated both ways, every interval is also given
over spatial blocks, the curved model is a check and not a bound, `psd_pc1`'s
floor takes one list of genes for both halves, and the density genes' agreement
with PSD95 is also read inside divisions (section 8.1).

### 4.1 The main model

```
nano map rank  ~  Gria1 rank + synapse-density rank
```

- *Abundance is Gria1 alone.* The mice carry SEP-GluA1: the stained protein is
  GluA1, which Gria1 alone encodes. Gria2 to Gria4 encode partner subunits the
  nanobody does not see, and their availability sets GluA1's assembly and
  trafficking: the surface side, not abundance. The four subunits are a check
  row.
- *Synapse density* is the mean rank of three postsynaptic genes, Rock2, Cap2
  and Slc8a2, chosen without the map by their agreement with the measured PSD95
  punctum density, after every gene that places or regulates AMPA receptors was
  left out as the surface side (section 4.2).
- *Two straight terms*, ranks across the structures, no interaction, scored on
  structures the fit has not seen (20 shufflings of five folds, as before); the
  share is the held-out R² over the ceiling, and the two weights are reported.
- *Why straight.* The question is whether the map follows the order of the two
  terms. A straight term on ranks credits a term with that and nothing more; a
  term bent to x² and x³ can also fit a map that rises and falls along it, and
  so take up parts of the map that do not follow its order at all, making the
  leftover smaller for a reason the question does not ask about. A straight
  model has a clean floor too: the calibration's known map is made of the two
  terms, so it holds no curvature the model could miss. The curve is not
  ignored: the model with each term curved (x, x², x³) is a check row, and
  control E measures what curving buys. That row is a check, not a bound: terms
  bent further take a little more (to fifth powers 36% is left, against 38% for
  cubes).
- *Why no autofluorescence term.* Autofluorescence is a property of the tissue,
  neither GluA1 abundance nor synapse density, the two things part 1 asks about.
  Alone it predicts nothing of the map held out (-1%), and whether the map is
  the label's or the tissue's is asked on its own (section 6.1). It stays a
  check row, which leaves 43% against the main model's 44%.

### 4.2 Synapse density: three genes chosen by a measured density

**The measured density.** `run_synaptome.py` (step 21),
`mapping/sepmap/adult/synaptome.py`; `figures\03s4_synaptome_detail.png`.
Numbers: step `synaptome`; `synaptome\density.csv`, `coverage.csv`,
`agreement.csv`. Allen mRNA sits in cell bodies; Zhu et al. (2018) counted
excitatory synapses where they are. In a knock-in mouse with PSD95 and SAP102
tagged (Dlg4-eGFP, Dlg3-mKO2), they detected every punctum in five coronal
sections of one adult male (about P80), sorted each into 37 subtypes by
intensity, size and shape, and measured the density of each subtype (puncta per
unit area) in regions of the Allen Reference Atlas. Hansen et al. share the
table: 37 subtypes by 775 samples, a sample being one region of one hemisphere in
one section, downloaded from their repository, netneurolab/hansen_synaptome,
`data/synaptome/mouse_liu2018/` at the pinned commit
`0399525412b6f50cfdeb5904b96da7fa8e4b507c` (the same files as release v1.0, Zenodo
doi 10.5281/zenodo.18201390), into `<data>\reference\synaptome\`, each file checked
against git's hash of it; `fetch_log.txt` there records the source.

Subtypes 1 to 11 hold PSD95 alone, 12 to 18 SAP102 alone and 19 to 37 both. As
shared, each subtype's density is divided by its largest value over the 775
samples, so absolute counts, and their sum, are gone. The PSD95 density is the
mean over the 30 subtypes whose puncta hold PSD95, each subtype's map weighing
alike; never a punctum's intensity or size, since PSD95 per synapse is
scaffolding, the surface side. A sample lies in the structure that is its Allen id
or the id's nearest ancestor (CCF 2017); a structure is the mean of its sampled
units weighted by their voxels in the 20 um CCF annotation; a region the source
gives only above several structures is never spread onto them, and a structure
with no sample stays missing. 739 of the 775 samples lie in 114 structures; the
others are listed with their reason in `synaptome\samples.csv`. One mouse,
sampled in a few sections: the measured density covers 96 of the 204 declared
structures and 89 of the 163 of the fit, so it chooses the density genes and is a
check row, never the term itself. Its two hemispheres agree at 0.97 (87
structures).

**The rule** (`run_density_markers.py`, step 22,
`mapping/sepmap/adult/density_markers.py`; `figures\03s3_density_markers.png`,
and its validation in figure 03 D; `density_markers\`). Fixed by Giulio, written,
run and committed with the genes it gave before any nano number of this version
was computed; the run reads only the gene table, the gene profiles, the declared
set and the synaptome's density table, and a test holds it to that.

- *The pool*: genes of the gene table with two or more usable Allen experiments
  after section QC (222 of 451), measured (`ish.min_voxels`) in at least 90% of
  the 164 declared structures where Gria1 is measured, and annotated in mouse GO
  (GO release 2026-08-05; the GAF of MGI and `go-basic.obo`, each checked against
  its SHA-256, in `<data>\reference\go\2026-08-05\`) to the postsynaptic density
  (GO:0014069) or the postsynaptic specialization (GO:0099572) or a term below
  either by is_a or part_of: 193 genes.
- *Excluded* as the surface side rather than density: Gria1 to Gria4, the
  localisation set, the AMPA receptor complex family, and every gene annotated
  (any evidence) to one of ten terms of AMPA receptors and receptor placement or
  below it (`EXCLUSION_TERMS` of the module; every id checked against its name in
  the release): 55 genes, Dlg4 and Camk2a among them; 138 remain. Each counted
  under its first reason (figure 03s3 E; every reason of every gene in
  `excluded.csv`): the subunits, 4; the AMPA receptor complex family, 12 (Cacng3,
  Cacng4, Cacng5, Cacng7, Cacng8, Cnih2, Dlg4, Grid1, Grid2, Olfm1, Olfm2,
  Shisa9); the localisation set, 14 (Agap3, Camk2g, Ctnnd1, Dlg2, Efnb2, Eps8,
  Erbb4, Gphn, Grid2ip, Igsf11, Neto1, Neto2, Nptx2, Pick1); regulation of AMPA
  receptor activity, 2 (Arc, Grip1); ionotropic glutamate receptor complex, 7
  (Grik1, Grik2, Grin1, Grin2a, Grin2b, Grin3a, Ptk2b); ionotropic glutamate
  receptor binding, 4 (Cdk5, Cdk5r1, Map1a, Nsf); regulation of
  neurotransmitter receptor localization to postsynaptic specialization
  membrane, 6 (Camk2a, Exoc4, Kalrn, Prkcz, Rapgef4, Tmem108); regulation of
  postsynaptic membrane neurotransmitter receptor levels, 6 (Clstn1, Grip2,
  Mapk10, Snap47, Vamp2, Vps35). Broad plasticity terms do not exclude; 26 pool
  genes carry one (`candidates.csv`), Slc8a2 among them. Homer1 passes: its
  receptor terms are metabotropic, and the indirect role of Homer1a in
  homeostatic scaling is not annotated.
- *The choice*: the three pool genes highest in Spearman with the PSD95 density
  over the declared structures where both exist: Rock2 0.808, Cap2 0.807 and
  Slc8a2 0.793 (91 structures each). The next are close (`agreement.csv`), and
  three excluded genes sit among them (Kalrn, Ptk2b, Camk2a).
- *Validation*, the choice repeated on 500 random halves of the 96 declared
  structures with PSD95: Cap2 is chosen on 67% of them, Rock2 63%, Slc8a2 37%
  (Fam81a 29%: the third gene is not stable). On the other half, the genes the
  rule chose on the first agree with PSD95 at a median 0.793 (95% 0.709 to
  0.857), the number to quote. It belongs to the rule, which chooses afresh on
  each half, not to the three genes, whose agreement on the full set, 0.825
  (91 structures), is optimistic; and random halves keep the contrast between
  divisions on both sides. Inside divisions (the mean rho within each division
  of 8 or more of the structures), the rule's choice agrees with PSD95 at 0.077
  on the other half (95% -0.339 to +0.409), the three genes at 0.174 on the full
  set: their agreement is mostly between divisions.
- *Beside it* (`comparison.csv`), each on its own structures, held out and
  inside divisions: the first proposal (Dlg4, Homer1, Camk2a) 0.823 (0.739 to
  0.880), inside divisions 0.398 (0.079 to 0.662); the 11 marker genes 0.732
  (0.597 to 0.839), 0.199; `psd_pc1` (151 postsynaptic-density genes measured in
  every structure with PSD95 and Gria1; 186 in the first version, on its own
  structures) 0.783 (0.694 to 0.849), 0.038. On the 77 structures all four have,
  held out: 0.785, 0.803, 0.732 and 0.784, the same order. So between divisions
  the rule's genes stand for synapse density about as well as the first
  proposal, two of whose genes it leaves out as AMPA-linked; inside divisions
  they hardly do, and the first proposal does best. The first proposal is a
  check row (52% left, section 4.4), so part 1 does not rest on the rule's genes
  standing for density inside divisions.

### 4.3 Method

- *The structures*: the declared ones where Gria1 and the three density genes
  are measured (`ish.min_voxels`, ten 200 um voxels, 0.08 mm³): 163 of 204. Left
  out, 40 without Gria1 (HY 17, MB 11, TH 8, STR 2, PAL 1, HPF 1) and 1 without
  Rock2 (IAD, TH) (`beyond\structures_used.csv`).
- *The ceiling*: two halves of the cohort, in all 126 splits, agree at 0.972;
  Spearman-Brown makes that 0.986 for ten adults. Each model is scored on
  structures it has not seen (held-out R², averaged over 20 shufflings of five
  folds), as a share of the ceiling.
- *The partition*: with A the share of Gria1 alone, D that of density alone and
  AD that of both, Gria1 only is AD - D, density only AD - A, shared A + D - AD
  and left 1 - AD; the four add to 1. Held out, a part can come out negative (a
  term that adds nothing to the other still costs a little through overfitting,
  or two terms together predict more than the sum of each alone); it would be
  reported as computed, signed and never clipped, and drawn leftwards from zero,
  hatched (figure 03 B). None is negative here.
- *The weights*: the straight model fitted on every structure with the map and
  both terms z-scored, so each weight is the map's SD per SD of the term, the
  other term held.
- *Intervals over structures*: a delete-d jackknife, 400 subsamples each leaving
  out a fifth of the structures at random, every predictor rebuilt on the
  structures kept; 95% is the value plus or minus 1.96 jackknife SD. Neighbouring
  structures are not independent (folds of spatial blocks leave more than random
  folds), so every interval is also given *over spatial blocks*: each of 20
  subsamples leaves out one of the 20 k-means blocks of the centroids that the
  folds of spatial blocks use, and the SD is that of the grouped jackknife. It is
  wider, and with 20 blocks itself uncertain (another block count or seed moves
  it by a few points); the random one is the method fixed first, the one over
  blocks the one that allows for neighbours.
- *Replication*: the model fitted to each half-cohort map and the two leftovers
  correlated, beside what the map's reliability and the fit alone imply.
- *The floor*: the same straight model on a map whose answer is known. Each
  gene's Allen experiments are split in two halves, alternately by id; the model
  fitted to nano with one half's predictors gives the known map (what Gria1 and
  synapse density predict of nano, and nothing else); ten made-up adults are the
  known map plus animal noise as large as nano's halves imply; the model with the
  other half's predictors reads them; both ways round, 20 draws each. On the 131
  structures where both halves measure Gria1 and the three genes; nano is read the
  same way on the same structures, and the difference is recomputed on each
  jackknife subsample. Every gene of the model has two or more usable experiments,
  so no gene is the same in both halves. Folds of spatial blocks have a floor of
  their own.
- *Which way the floor errs.* Both ways. What it lacks makes it low: the
  mismatch of Allen's P56 mice with these brains (age, strain, the 200 um grid
  against 20 um masks, registration), which nano meets and the known map does
  not. What it holds makes it high: the known map is built from one half of the
  experiments and read with the other, so it meets the errors of two halves,
  where nano, read with one half, meets one. Read with its own half, the known map
  leaves 1% (the animals' noise alone), so about 14 of the floor's 15 points are
  two halves disagreeing, and one half's error is roughly half of that. The floor
  errs low only if the P56 mismatch costs more than that. Nano itself leaves 42%
  with one half's predictors and 41% with both halves merged.
- *The floor is the only test.* The first version also set nano beside a second
  benchmark, a map that is one Allen Gria1 experiment, with made-up animals built
  from it. Giulio removed it when he decided this version, after seeing the
  first, as unfair: its made-up animals share one Allen brain's quirks, which then
  count as reproducible, while the nano map averages ten brains. With two or more
  experiments behind every gene of the model, the floor has no such part.
- *Seven controls* try to break the result (a spatial gradient, structure size,
  single animals, naive against RWS, curvature, the choice of predictors, the
  reading).
- *The check rows* (`[beyond.checks]`), each one change to the main model, on its
  own structures (the declared ones where its predictors are measured, or the
  main model's cut as it says), with its ceiling recomputed there and a floor of
  its own, built the same way with its own model, and nano minus that floor with
  its two intervals (`beyond\check_rows.csv`): the model curved; +
  autofluorescence; the first proposal for density (Dlg4, Homer1, Camk2a); the 11
  marker genes of the first version; + `psd_pc1`; Gria1 to Gria4 in place of
  Gria1; the measured PSD95 density as the density term, on the structures it
  covers, with the main model on the same structures beside it; only the
  structures of at least 0.4 mm³ (50 voxels of the CCF annotation on the Allen
  grid); and nano measured on the Allen 200 um grid itself, each 200 um voxel the
  mean of its tissue and labelled as the ISH grids are (`ish.regions`), with zref
  taken the same way (`tables\adult_allen_grid.csv`, step 13). `psd_pc1` is the
  first component of the panel's postsynaptic-density genes measured on every
  structure (122 on the 163); in the calibration both halves take one list, the 69
  of them with two halves of experiments measured on every structure in both, so
  no gene of it is the same in both halves. A floor misses the mismatch of what is
  measured once, the same in both halves: the PSD95 density (one mouse), and three
  of the 11 markers (Shank2, Shank3, Nlgn1).

### 4.4 Result

`figures\03_beyond.png`: A, the map against the main model's held-out
prediction; B, the reproducible map in four parts, with the interval of what is
left; C, nano and the floor on the same 131 structures, and nano minus the floor
on an axis of its own with its two intervals; D, the density term against the
measured PSD95 punctum density, with the rule's held-out rho whole brain and
inside divisions. `03s1_beyond_budget.png` leads with the check rows and the main
model under other folds, each with what it leaves and nano minus its own floor,
then the map against each term, every model's share, the parts with their
intervals and the weights, every calibration draw and the replication;
`03s2_beyond_controls.png` draws the seven controls, `03s3_density_markers.png`
the choice of the density genes, `03s4_synaptome_detail.png` the measured
density. Tables: `beyond\variance_partition.csv`, `partition.csv`, `weights.csv`,
`jackknife.csv`, `jackknife_blocks.csv`, `calibration.csv`,
`calibration_jackknife.csv`, `calibration_jackknife_blocks.csv`,
`replication.csv`, `controls.csv`, `folds.csv`, `check_rows.csv`,
`residual_by_division.csv`.

| held out, on 163 structures | share of the reproducible map (95% over structures; over spatial blocks) |
|---|---|
| Gria1 alone | 41% (27% to 55%; 12% to 70%) |
| synapse density alone | 46% (35% to 56%; 17% to 75%) |
| both, the main model | 56% (46% to 65%; 38% to 74%) |
| autofluorescence alone | -1% |
| both curved (the check row curved) | 62% |

| part of the reproducible map | share (95% over structures; over spatial blocks) |
|---|---|
| Gria1 only | 10% (2% to 18%; -10% to +31%) |
| shared | 31% (22% to 39%; 16% to 46%) |
| density only | 15% (4% to 26%; -13% to +43%) |
| left | 44% (35% to 54%; 26% to 62%) |

- The weights, map and terms z-scored: Gria1 +0.39 (+0.23 to +0.54; over
  spatial blocks -0.01 to +0.78), synapse density +0.47 (+0.32 to +0.61; +0.01
  to +0.92). The two terms agree at 0.56 across the structures, which is why a
  third of the reproducible map is shared, and why, once whole blocks of
  neighbours are left out, the split between the two terms is loosely fixed, more
  loosely than what the two leave together.
- The leftover replicates: the two half-cohort leftovers agree at 0.933 (lowest
  split 0.865). The map's reliability and the fit alone imply 0.938; two
  unrelated leftovers fall between -0.15 and +0.15.
- The floor, on 131 structures: a map made only of Gria1 and synapse density
  leaves 15% (11% to 17% over 40 draws), its leftover replicating at 0.81; read
  with its own half, 1%. Nano, read the same way, leaves 42% (31% to 53% over
  resampled structures): **27 points above the floor (95% +15 to +40; over
  spatial blocks +10 to +45)**. With folds of spatial blocks nano leaves 48% and
  the floor 15%.
- The check rows, each against its own floor (nano minus the floor in points, 95%
  over structures; over spatial blocks):

  | check row | structures | left | calibration structures | nano minus its floor |
  |---|---|---|---|---|
  | the main model | 163 | 44% | 131 | +27 (+15 to +40; +10 to +45) |
  | curved: x, x², x³ of each term | 163 | 38% | 131 | +20 (+8 to +32; +5 to +35) |
  | + autofluorescence | 163 | 43% | 131 | +23 (+13 to +34; +11 to +36) |
  | density: Dlg4, Homer1, Camk2a | 163 | 52% | 130 | +30 (+17 to +42; +9 to +50) |
  | density: the 11 marker genes (Shank2, Shank3, Nlgn1 measured once) | 126 | 49% | 115 | +25 (+9 to +40; +6 to +43) |
  | + psd_pc1 | 163 | 43% | 131 | +31 (+18 to +44; +12 to +50) |
  | abundance: Gria1 to Gria4 | 163 | 28% | 126 | +13 (+4 to +22; +0.1 to +26) |
  | density: PSD95 puncta (measured once) | 89 | 50% | 81 | +31 (+10 to +51; +5 to +56) |
  | the main model where PSD95 is measured | 89 | 50% | 78 | +27 (+7 to +46; -0.2 to +53) |
  | structures of 0.4 mm³ or more | 118 | 49% | 105 | +24 (+8 to +40; +3 to +46) |
  | nano on the Allen 200 um grid | 152 | 42% | 126 | +29 (+16 to +41; +9 to +49) |

  Every row stays above its own floor at its point value and over its 95% over
  structures; over spatial blocks 9 of 10 do, the main model on the 89 PSD95
  structures reaching zero (-0.2) and the four subunits just above it (+0.1). The
  partner subunits take the most (Gria1 to Gria4 leave 28%), then curvature
  (38%). Where PSD95 is measured, it predicts the map about as well alone as the
  density genes do (19% against 16%) and the model leaves as much with it (50%
  and 50%). The 89 structures with PSD95 are not a random part of the 163 (few
  hypothalamic, pallidal and midbrain structures), so those two rows are compared
  with each other, not with the main value.
- Under other folds the main model leaves 43% with one shuffling (single
  shufflings 43% to 46%), 44% with ten folds and leave one out, and 53% with
  folds of spatial blocks (`beyond\folds.csv`).
- Six of the seven controls pass (`figures\03s2_beyond_controls.png`;
  `beyond\controls.csv`): a smooth gradient explains R² 0.143 of the leftover,
  which still replicates at 0.934 with position in the model; its rho with
  structure volume is +0.038; single adults' leftovers agree at a median of
  +0.762 (lowest pair +0.578); the naive and RWS groups' at +0.895; every reading
  gives a leftover that replicates at 0.93 to 0.94. Control F: the components of
  the 209 genes measured in every structure (most often 20, picked inside each
  training fold) predict 81%; on its own calibration nano leaves 22% to 24%,
  against 3% to 5% for a map made of those genes; its held-out R² is within 0.01
  of its highest from 30 components on, and highest at 40, the most
  `[beyond_controls] max_pcs` allows. **Control E does not pass**: held out, the
  straight model reaches R² 0.551, the model curved 0.615 and fifth powers 0.636
  (36% left), so curving buys more than the 0.05 the control allows. Part of the
  straight model's leftover is curvature; the curved check row shows how much,
  and it stays above its own floor (+20 points, +8 to +32; over spatial blocks +5
  to +35).

**Where it sits** (`figures\04_beyond_where.png`; `beyond\regression_table.csv`,
`residual_by_structure.csv`, `residual_by_division.csv`). In large part it is a
contrast between divisions, which two straight terms do not reproduce: 40% of the
leftover's variance lies between divisions. Above prediction as a whole (95% over
resampled adults above zero): the hippocampal formation (mean +36 ranks), the
striatum (+33), olfactory areas (+25), the cortical subplate (+14), the pallidum
(+13) and the midbrain (+8); below it the hypothalamus (-23), the isocortex (-18)
and the pons (-18, two structures); the thalamus sits on prediction (+1). Of single
structures, above prediction most in the triangular nucleus of the septum (+91
ranks), the septofimbrial nucleus, the medial geniculate, the indusium griseum,
both parts of the lateral geniculate, the lateral posterior nucleus, the lateral
septum, the subiculum and CA2; below it most in the preoptic part of the
periventricular hypothalamic nucleus, the paraventricular hypothalamic nucleus,
the rhomboid nucleus, the descending division of the paraventricular hypothalamic
nucleus, the nucleus of reuniens and the arcuate nucleus. The leftover follows
the autofluorescence map at -0.19.

### 4.5 What it means

About 44% of the map's reproducible pattern is not predicted by Gria1 expression
and synapse density measured by three postsynaptic genes, and it is reproducible
across mice. On the same structures it is well above what Allen-to-Allen mismatch
leaves (27 points, 95% +15 to +40, and +10 to +45 when whole spatial blocks are
left out), and every check row stays above its own floor at its point value; the
lowest margins are those of the four subunits (+13) and of the curved model
(+20).

What it does not mean. The share depends on the model: curved it is 38%, and with
the partner subunits Gria2 to Gria4 as abundance 28%, so part of the straight
model's leftover is curvature and part is the partner subunits, which both
versions count on the surface side. The share itself is loosely fixed once
neighbouring structures are allowed for (26% to 62% over spatial blocks), and
so is the split between Gria1 and density; the margin over the floor holds under
both resamplings. The floor errs both ways (section 4.3), by less than the
margin. The density genes follow PSD95 mostly between divisions, so inside
divisions the density term is weak; the first proposal, which follows PSD95
inside divisions too, leaves more (52%). Its replication follows from the map's
reliability (0.938 expected), so it is not separate evidence. It is "not
predicted by Gria1 expression or synapse density", not "beyond gene expression":
the components of many genes predict 81% of the map. And a leftover says what the
predictors miss, not what it is; a claim about one structure needs its own null.

### 4.6 The first version, kept for the record

Part 1's first version was fixed on 8 October and run on the morning of 9 October
(section 8.1); this second version replaced it the same day, after its numbers
were seen. Its outputs are on the development copy in
`adult_v2\ish_analysis_9oct_v1\` (`beyond\numbers_for_the_caption.txt` and
`tables\numbers_for_the_text.csv` there); the numbers below are from them. That
folder also holds the first run of this version's step 22 (`density_markers\`,
`tables\numbers_density_markers.csv` and `figures\14s2_density_markers.*`),
written there at 10:24 to 10:35 before the folder was set aside; step 22 reads no
nano value, so they belong to the second version.

- *The model*: Gria1 + synapse density + autofluorescence, each term as x, x² and
  x³. PSD95 punctum density covered 77 of the 126 structures of its fit, fewer
  than the 101 its rule asked, so synapse density was the mRNA panel: the 11
  synaptic marker genes and `psd_pc1`, the first component of 186
  postsynaptic-density genes. 126 structures, where every subunit and marker is
  measured.
- *The result*: two halves of the cohort agreed at 0.974 (98.7% reproducible);
  held out, Gria1 predicted 47% of the reproducible map, synapse density added
  25% and autofluorescence -1%; 29% was left (95% 15% to 42%), its two
  half-cohort leftovers agreeing at 0.933 (0.897 implied by the map's
  reliability and the fit). Density alone predicted 58%, Gria1 to Gria4 as four
  terms 67%; the same model with straight terms left 40%.
- *The floor*, on 113 structures: 12% (11% to 15% over 40 draws); nano 26% (13%
  to 40%), 14 points above it (95% -2 to +30), so above the floor at its point
  value only. That floor missed what is measured once: Shank2, Shank3, Nlgn1 and
  25 of the 186 genes behind `psd_pc1`.
- *Its check rows*, on the 77 structures PSD95 covers: PSD95 puncta as the
  density term 46% left, the mRNA panel 39%, both 39%, SAP102 puncta 48%, every
  punctum 47%; the panel without its five presynaptic markers 31% and Gria1 to
  Gria4 27% on all 126. 7 of 7 controls passed.
- *The second benchmark*, a map of one Allen Gria1 experiment, left 34%, and nano
  did not stand above it (nano minus it -8 points, 95% -29 to +14). Section 4.3
  says why the second version drops it.
- *The genes against its leftover*: Cacng8 +0.20 (p 0.0002); the AMPA receptor
  complex family p 0.047 against the surrogates (0.050 without Cacng8) and 0.55
  against matched postsynaptic genes; no gene past BH over all 451 (27 below p
  0.05). Added to that model, Cacng8 took 4.6 points of the reproducible map, p
  0.059 against its plain surrogates and 0.009 against maps alike to the model.

What changed, and what it did: abundance stays Gria1; synapse density becomes
three genes chosen by a rule that reads no nano value, whose two AMPA-linked
competitors (Dlg4, Camk2a) and every other AMPA-linked gene it leaves out; the
terms become straight and autofluorescence leaves the model; the structures grow
from 126 to 163; the floor becomes the only test. The leftover grows from 29% to
44%, mostly from straightening the terms: the first version's own model with
straight terms already left 40%. The other four points or so come from the new
density term, the structures and autofluorescence leaving the model. The partner
subunits are in neither version's main model: as a check row they bring the
first version's leftover from 29% to 27% and the second's from 44% to 28%, so
the first version's mRNA panel already carried most of what they carry. The
margin over the floor grows from +14 (-2 to +30) to +27 (+15 to +40).

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
- whether it is a Gria1-like or density-like map: its rho with Gria1, the
  density term, autofluorescence and the main model's prediction on the 163
  structures of the fit, and with the PSD95 density where it is measured;
- with the leftover: its rho and spatial p, each surrogate through the same fit
  (section 5.7);
- what it takes of the leftover: the main model with the gene's ranks added as
  one more straight term, as each of the model's is, scored on held-out
  structures as the main model is, as a share of the reproducible map. A smooth
  map takes some by chance, so the gene is set beside 1,000 maps of its
  smoothness put in its place, of two kinds. Plain surrogates of the gene are
  unrelated to the model, so more of each is new to it than of a gene that shares
  the model's pattern: a wide null, the conservative bound, the one this step was
  first given and the one quoted first (Giulio, 9 October). Maps that relate to
  the model as the gene does keep the gene's fit on the model and replace its
  remainder with a surrogate of it, as large; they were added on 9 October after
  the plain null's numbers were seen, and are the secondary line;
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
measured once; the share taken is in points of the reproducible map, with its p
against plain surrogates, the null fixed first, and against maps alike to the
model):

| gene | with the map, rank | inside divisions | rank over the variants | reliability | with Gria1, with the model's prediction | with the leftover (p) | takes of the leftover (p plain, alike) |
|---|---|---|---|---|---|---|---|
| Cacng8 | +0.807, 1 | +0.548 (past) | 1 to 2 | 0.91 | +0.75, +0.88 | +0.200 (0.0066) | 9.8 (0.001, 0.003) |
| Arpc5 | +0.751, 2 | +0.356 (past) | 1 to 7 | - | +0.54, +0.76 | +0.308 (0.0003) | 8.6 (0.004, 0.005) |
| Igsf11 | +0.743, 3 | +0.210 | 2 to 15 | 0.85 | +0.57, +0.76 | +0.230 (0.017) | 6.7 (0.006, 0.007) |
| Htr3a | +0.737, 4 | +0.443 (past) | 3 to 31 | - | +0.57, +0.76 | +0.227 (0.072) | 5.7 (0.013, 0.003) |
| Grm5 | +0.725, 6 | +0.382 | 2 to 14 | - | +0.79, +0.90 | +0.056 (0.43) | 0.1 (0.575, 0.375) |
| Neurl1a | +0.718, 8 | +0.467 (past) | 3 to 34 | 0.81 | +0.69, +0.83 | +0.108 (0.27) | 2.2 (0.143, 0.152) |
| Mapk1 | +0.687, 17 | +0.213 | 8 to 38 | - | +0.61, +0.79 | +0.185 (0.11) | 1.0 (0.352, 0.144) |
| Add3 | +0.681, 20 | +0.292 | 4 to 89 | - | +0.40, +0.58 | +0.283 (0.0052) | 11.1 (0.001, 0.001) |
| Gria1 (the abundance term) | +0.646, 32 | +0.417 (past) | 18 to 47 | 0.91 | +1.00, +0.85 | +0.013 (0.50) | in the model |
| Cnih2 (two Allen experiments that disagree) | +0.637, 36 | +0.336 | 5 to 43 | 0.08 | +0.56, +0.72 | +0.104 (0.19) | 1.2 (0.205, 0.198) |
| Grip1 | +0.615, 45 | +0.207 | 22 to 75 | 0.65 | +0.74, +0.82 | -0.025 (0.71) | -0.5 (0.811, 0.706) |
| Eps8 | +0.605, 54 | +0.286 | 8 to 79 | 0.55 | +0.50, +0.62 | +0.204 (0.010) | 3.5 (0.071, 0.051) |

- Four of the 12 are localisation genes and one is a subunit; one of the four,
  Cnih2, has two Allen experiments that disagree (reliability 0.08), so its
  place is weak.
- They are maps much like Gria1 and synapse density: rho 0.40 to 0.79 with
  Gria1 and 0.58 to 0.90 with the main model's prediction. That is why the main
  model predicts much of the map.
- Cacng8 is first of all genes, and first of P9's genes in every robustness
  variant (section 6.2; under Pearson on log2 it is 2nd of all genes). Gria1 is
  11th of P9's genes (q 0.0025 within them) and 32nd of all.
- 5 of the 12 follow the leftover before correction (Cacng8, Arpc5, Igsf11,
  Add3, Eps8); only Cacng8's p is a test, and a re-test (section 5.7). Gria1 is
  a term of the model, so its rho with the leftover is near zero by construction
  and its p says nothing. Added to the main model as one more straight term, five
  take more than 95% of their plain surrogates and more than 95% of maps alike to
  the model: Cacng8, Arpc5, Igsf11, Htr3a and Add3.
- `07s_gene_ranking.png` shows P9's 100 genes one by one with their null bands
  and their rho with the autofluorescence map.

**What it means.** The genes that follow the map most describe it; they are
not separate evidence, since most of what they share with the map is what Gria1
and synapse density already predict. Five carry part of what the straight model
leaves, against both nulls. Of these, Cacng8 is the one named for the leftover
(section 5.7); the others are exploratory findings, Add3, which takes the most,
and Htr3a, a serotonin receptor, each measured by one Allen experiment.

### 5.4 Cacng8 against Gria1

`figures\08_cacng8_gria1.png`; the gap in detail in `07s_gene_ranking.png` B.
Numbers: steps `gene_ranking` and `top_genes`; `tables\gap.csv`.

The map follows Cacng8 at +0.807 and Gria1 at +0.646, inside divisions +0.548
and +0.417; the two genes agree at +0.75 over the structures of the fit. The
Cacng8 - Gria1 gap, on the 164 structures both have, is +0.167, steady across
adults (+0.156 to +0.179 over resampled adults; above zero in 10 of 10 adults).
The test fixed in advance is two-sided, against maps related to both genes
alike: a gap as large either way, p 0.105, so the gap is inside it; it would
pass beyond ±0.195, which figure 08 C shades. That null is skewed, which the
shares each way describe: 2.3% of those maps give a Cacng8 lead this large and
8.2% a Gria1 lead this large (its central 95% runs from -0.224 to +0.166, so a
one-sided reading would put the gap at the edge; the test fixed in advance is not
one-sided). Against maps unrelated to both, a
wider null and a conservative bound, p 0.348. It moves with the Allen experiment
that stands for each gene, from +0.08 to +0.211 over the four pairings (p 0.064
to 0.385).

**What it means.** The map follows a TARP's pattern and Gria1's beyond a map
with the brain's smoothness. Whether it follows Cacng8 more closely than Gria1
these maps do not settle: the gap is above zero in each of the 10 adults and in
every variant, and inside the test fixed in advance, which stays two-sided
(Giulio, 9 October). Figure 08's line under its title says no more than that
test. Cacng8's own test is against the leftover (section 5.7).

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
   of the four-subunit model, 0.0024, was seen before it was named, and its p
   against the first version's leftover (0.0002) before the second version was
   decided (`ish_analysis_8oct_four_subunits\` and `ish_analysis_9oct_v1\` on
   the development copy). So tier 1 re-tests a result already seen: it shows the
   result holds under the second model, it does not confirm it.
2. The AMPA receptor complex family as a group: its median rho with the
   leftover against the same genes' median over the leftover's 10,000
   surrogates, and against 31 genes of the other postsynaptic set outside it,
   each the one closest in median energy (the localisation test's matching),
   labels permuted 20,000 times; one test each. Then gene by gene, BH within
   the family only. The family comes from sources outside this analysis: the
   constituents of native AMPA receptor complexes found by proteomics (Schwenk
   et al. 2012, Figure 1D and Table S2), the mouse genes annotated to
   GO:0032281 AMPA glutamate receptor complex (in MGI's GAF as the GO Consortium
   served it on 8 October, dated 2026-05-21; release 2026-08-05 lists the same
   31 genes), and the partner subunits Gria2 to Gria4, less Gria1, the
   abundance term: 37 genes, 31 with a usable map. Not tested: Olfm3, Prrt1,
   Prrt2 and Shisa8, which have
   no Allen experiment; Gsg1l and Rap2b, which have two each but are in neither
   panel of the gene table, and adding them would change the rules. Every
   gene's rho with the leftover of the four-subunit model was seen on 8 October;
   the family comes from the paper and GO, not from those numbers. In the second
   version no member is a term of the model: Dlg4, a marker of the first
   version's density, is left out of the density genes by their rule, so it is
   tested as any member.
3. Every other gene, exploratory, BH over all of them, the model's own genes
   (Gria1, the density genes) included as the rule set it. With 10,000 surrogates the
   smallest p is 0.0001, just under the line for the first gene (0.05 / 451 =
   0.00011): a gene can cross it only if no surrogate is as extreme. The gene
   sets of section 5.6 are read the same way.

Two check rows of tier 2 were added after the tests ran, to say what the group
result rests on, never in its place: the family without Cacng8 against the same
surrogates, and the family against controls matched afresh from the postsynaptic
genes outside the model's terms (none of the 31 controls is a model term, so this
check repeats the test). And each surrogate of the leftover, once the model is
projected out of it, is rougher at short range than the leftover itself (each
map's Spearman with the mean of every structure's five nearest neighbours:
leftover 0.804, its surrogates 0.724), which narrows the null; the tests are read
again against Gaussian fields of exponential covariance with ranges of 2, 6, 10,
20 and 40 mm, 5,000 each, projected the same way
(`tables\leftover_null_check.csv`). None reaches the leftover's smoothness: the
smoothest, at 40 mm, 0.759.

The same group tests are read on the map itself, to describe the family, not as
tests of the leftover.

**Result.**

| test | on the leftover | on the map (a description) |
|---|---|---|
| tier 1, Cacng8 alone | +0.200, p 0.0066 (18th of 451 by rho) | +0.807, p ≤ 0.0001 |
| tier 2, the family's median against the surrogates | -0.004 (null 95% -0.085 to +0.083), p 0.92 | +0.414, p 0.044 |
| tier 2, the family against matched controls | -0.004 against -0.067: +0.063, p 0.29 (±0.117 would pass) | +0.414 against +0.365: +0.049, p 0.50 |
| tier 2, members past BH within the family | Olfm1, Olfm2, Gria4, all with a negative rho; Cacng8 q 0.051 | Cacng8 and Cnih2 |
| tier 2 check, the family without Cacng8 | -0.018, p 0.71 | |
| tier 2 check, against controls outside the model | the same 31 controls: no control is a model term | |
| tier 3, genes past BH over all 451 | 53, 50 of them with a negative rho (117 below p 0.05) | |

| null of the leftover | smoothness | Cacng8's p | the family's p | genes below p 0.05 | genes past BH |
|---|---|---|---|---|---|
| its surrogates, through the fit | 0.724 | 0.0066 | 0.92 | 117 | 53 |
| Gaussian fields, 2 mm | 0.646 | 0.0020 | 0.90 | 124 | 59 |
| Gaussian fields, 6 mm | 0.733 | 0.0018 | 0.91 | 118 | 46 |
| Gaussian fields, 10 mm | 0.746 | 0.0030 | 0.91 | 116 | 47 |
| Gaussian fields, 20 mm | 0.756 | 0.0020 | 0.92 | 115 | 50 |
| Gaussian fields, 40 mm | 0.759 | 0.0022 | 0.92 | 115 | 45 |

The leftover's own smoothness is 0.804, above every null's: the smoother nulls
narrow the gap, they do not close it. Cacng8's p stays between 0.0018 and 0.0066,
the family's between 0.90 and 0.92, and tier 3 counts 45 to 59 genes past BH.

- Added to the main model as one more straight term, Cacng8 takes 9.8 points of
  the reproducible map (44.1% left, 34.3% with it): more than 95% of its plain
  surrogates take (p 0.001), the null fixed first, and more than maps alike to
  the model (p 0.003), the secondary line.
- The members of the highest rho are Shisa6 (+0.268, p 0.16), Dlg3 (+0.255, p
  0.021), Cacng8 and Vwc2l (+0.178); within the family Cacng8 misses BH by a hair
  (q 0.051). The three members past BH run the other way: Gria4 (-0.524), Olfm2
  (-0.451) and Olfm1 (-0.215), high where nano sits below prediction. The partner
  subunits Gria2 (-0.182) and Gria3 (+0.068) sit near zero. In the family, the
  two Allen experiments of Cnih2 and of Cacng4 disagree (reliability 0.08 and
  0.20), and Shisa6, Dlg3 and Lrrtm4 have one experiment each.
- Dlg4 and Camk2a, markers of the first version's density, are no terms of this
  model: the rule of the density genes left them out as AMPA-receptor-linked. Dlg4
  is tested as any member (-0.115, p 0.017); Camk2a falls among the genes past BH
  of tier 3 (-0.147).
- Tier 3, exploratory: 53 genes pass BH over all 451, against 0 in the first
  version. Three have a positive rho (Arpc5, Add3, Amot); 50 a negative one, most
  of them pan-neuronal, synaptic-vesicle and ribosomal genes (Gria4, Rpl8, Mapt,
  Olfm2, Snap47, Sema4f, Arrb1, Hsp90aa1, Cd200, Stx4a; Syp, Syt1, Snap25, Stx1a
  and Syn1 among the rest): the nano map sits below what Gria1 and the three
  postsynaptic genes predict where neuronal mRNA in general is high. The count
  moves little with the null (115 to 124 below p 0.05, 45 to 59 past BH). There
  is a second reading: these are largely the genes the autofluorescence map
  follows (section 6.1). All 53 genes past BH have a positive rho with the
  autofluorescence map, 15 of them past its own null, and over every gene a
  gene's rho with the leftover and its rho with autofluorescence agree at -0.36
  (figure 11s1 C), while the leftover itself follows the autofluorescence map at
  only -0.19 and the model with autofluorescence leaves 43%. Which of the two
  readings holds, neuronal mRNA in general or the tissue, is not settled here.
  The highest rho over every gene is Sox10's (+0.329, p 0.029, q 0.13), then Caln1
  and Aqp4; the glia set has a median of +0.128 (p 0.31, q 0.45).

**What it means.** The gene named for the leftover, a TARP, follows what Gria1 and
synapse density leave in this version too, under every null, and taken as one
more term it takes more than its own plain surrogates; but it re-tests what two
earlier leftovers showed, so it shows the result holds, it does not confirm it.
The AMPA receptor complex family as a group does not follow this leftover, and the
members that pass within it run against it. What the straight model leaves is
broader than the surface side: much of it follows genes expressed in every neuron,
negatively, which are also the genes the tissue's autofluorescence follows,
together with the curvature and the contrast between divisions of section 4.4. So
the genes neither contradict nor confirm part 2; only a total-receptor measurement
can.

### 5.8 What part 2 says

Consistent with the surface-fraction reading, and not singling it out. The map
follows receptor expression and not the tissue (section 6.1), and Cacng8, Dlg2
and Gria1 lead the ranking inside divisions too. But the genes that follow the
map most are maps much like Gria1 and synapse density, the Cacng8 - Gria1 gap
is inside the two-sided test fixed in advance, no gene set passes its null, the
postsynaptic criterion named in advance is met by half, and the localisation
genes predict the map no better than matched postsynaptic controls
(no advantage larger than about +0.08). Against the leftover, the TARP named for
it follows it (p 0.0066) and takes part of it beyond its plain surrogates, the one
result here that points the way the reading does, though a re-test of what two
earlier leftovers showed. The AMPA receptor complex family named with it does not
follow the leftover, and what the leftover does follow most, negatively, is the
mRNA of genes expressed in every neuron, which the tissue's autofluorescence
follows too. The genes neither contradict nor confirm part 2; only a
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
stays between +0.126 and +0.211. The primary's two-sided test passes a lead only
beyond ±0.195 (figure 13 C shades it; each variant's own null is not drawn): every
variant but P9's single experiment per gene (+0.211) sits inside it.

### 6.3 A measured synapse density

`figures\03s4_synaptome_detail.png`, figure 03 D and the rows of 03s1; section
4.2. The measured PSD95 density chooses the density genes and is a check row. It
covers 89 of the 163 structures of the fit; there the density genes order the
structures as PSD95 does at 0.83 (89 structures; optimistic, since they were
chosen by it; on the 91 declared structures with PSD95 0.825, and the rule
re-run on random halves 0.79 on the other half, mostly between divisions),
PSD95 alone predicts the map about as well as the density genes (19% against
16%), and the model leaves as much with it (50% against 50%), each above its own
floor (+31 points, 95% +10 to +51; +27, +7 to +46; over spatial blocks +5 to +56
and -0.2 to +53). It checks part 1 from outside the ISH data.

### 6.4 The limit: the green channel (analysis 5)

`run_sep_channel_check.py`; `figures\14_green_channel.png` and
`14s_green_channel_detail.png`; `green_channel\sep_channel_check.csv`. Numbers:
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
- The calibration floor errs both ways: it lacks the mismatch of Allen's P56
  mice with these brains (strain, age, the 200 um grid against 20 um masks,
  registration), which makes it low, and it holds the disagreement of two Allen
  halves where nano meets one, which makes it high (about 7 of its 15 points);
  nano on the Allen grid, a check row, answers the grid only. A check row's floor
  also misses what it measures once (PSD95, three of the 11 markers).
- Neighbouring structures are not independent. Over spatial blocks every
  interval of part 1 is wider: what is left 26% to 62%, nano minus the floor +10
  to +45, and the split between Gria1 and density hardly fixed.
- The straight model leaves curvature (control E) and the partner subunits in
  the leftover; their check rows show how much (not a bound: terms bent further
  take a little more), and both stay above their own floors.
- The density genes follow the measured PSD95 density mostly between divisions;
  inside them the rule's choice agrees at 0.08 held out, so the density term
  stands for synapse density at the scale of divisions.
- The leftover's surrogates, once the model is projected out of them, are
  rougher than the leftover, and no null tried reaches its smoothness; the named
  tests hold against smoother nulls, and per-gene counts move little with the
  null (45 to 59 genes past BH).
- The measured synapse density is one mouse, in 89 of the 163 structures of
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

`run_ish_overview.py`; `figures\15_april_headline.png` and
`15s_april_headline_detail.png`. Numbers: step `overview`;
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

- **The first version of part 1, its rules fixed first (8 October).** Written
  down and committed before it was run (what was known by then: each term's share
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
  - the variants (`[beyond.variants]` of that version) show how much the
    choice of density measure matters. The rules called PSD95 and the panel
    together the conservative bound, since one measured map in place of the
    panel's two terms gives the model fewer terms; scored on held-out structures
    it was not a bound. The calibration runs the same main model; PSD95 is one
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

  Run on the morning of 9 October: PSD95 covered 77 of the 126 structures, so
  that model kept the mRNA panel; it left 29% (15% to 42%), 14 points above its
  floor (-2 to +30).
  Cacng8 followed its leftover (p 0.0002); the family passed against the
  surrogates only at the edge (p 0.047, 0.050 without Cacng8), not against
  matched controls (p 0.55). Its outputs are kept on the development copy in
  `adult_v2\ish_analysis_9oct_v1\`, those of the first run of 8 October, with
  the four subunits, in `adult_v2\ish_analysis_8oct_four_subunits\`.

- **The density term of the second version of part 1, chosen without the map
  (9 October).** Giulio fixed the rule after seeing the first version's
  numbers, and it was written, run and committed with the genes it gave before
  any nano number of the second version was computed (step 22,
  `mapping/sepmap/adult/density_markers.py`, `[density_markers]` of
  `mapping/settings.toml`; the run reads only the gene table, the gene
  profiles, the declared set and the synaptome's density table, and a test
  holds it to that):
  - the pool: genes of the gene table with two or more usable Allen
    experiments after section QC, measured (`ish.min_voxels`) in 90% of the
    declared structures where Gria1 is, and annotated in mouse GO (GO release
    2026-08-05) to the postsynaptic density (GO:0014069) or the postsynaptic
    specialization (GO:0099572) or a term below either by is_a or part_of;
  - excluded, as the surface side rather than density: Gria1 to Gria4, the
    localisation set, the AMPA receptor complex family, and every gene
    annotated (any evidence) to one of ten terms of AMPA receptors and receptor
    placement (`EXCLUSION_TERMS`) or below it. Broad plasticity terms do not
    exclude; the pool genes carrying them are listed in `candidates.csv`.
    Checked by hand before the run: Dlg4 and Camk2a fall to this rule, Homer1
    passes (its receptor terms are metabotropic; the indirect role of Homer1a in
    homeostatic scaling is not annotated);
  - the choice: the three pool genes highest in Spearman with the measured PSD95
    punctum density over the declared structures where both exist, their mean
    rank the term; repeated on 500 random halves of those structures, the
    agreement of the chosen genes' mean with PSD95 on the other half is the one
    quoted.

  The rule chose Rock2, Cap2 and Slc8a2 (figure 03s3; numbers: step
  `density_markers`).

- **The second version of part 1 (9 October).** Decided by Giulio after seeing
  every number of the first version, for reasons about the measurement; the
  density genes were committed before any nano number of this version was
  computed, and the rest written into the code before it ran:
  1. the model is nano rank ~ Gria1 rank + synapse-density rank, two straight
     terms, no autofluorescence term, no interaction, cross-validated as before;
     the share is the held-out R² over the ceiling, split into Gria1 only,
     density only, shared and left, with the two weights;
  2. synapse density is the mean rank of the three genes of the rule above;
  3. the structures are the declared ones where Gria1 and the three genes are
     measured;
  4. the test is the floor only, with the same straight model, and the
     one-experiment Gria1 benchmark is removed everywhere: its made-up animals
     share one Allen brain's quirks, which then count as reproducible, while the
     nano map averages ten brains;
  5. beside the main value, check rows, each with its own floor where one can be
     built: curved, + autofluorescence, the first proposal for density, the 11
     markers, + `psd_pc1`, the four subunits, PSD95 as the density term on its
     structures, only structures of 0.4 mm³ or more, and nano on the Allen grid;
  6. the gene tests in the three tiers already committed, against the new
     leftover, Dlg4 and Camk2a no longer terms of the model; the Cacng8 - Gria1
     gap stays two-sided, as fixed, and figure 08's line under its title does not
     say the map follows Cacng8 more closely than Gria1; Cacng8's share is quoted
     against its plain surrogates, the null fixed first, with maps alike to the
     model as the secondary line.

  What was known: every number of the first version (section 4.6), among them
  its `variance_partition.csv`, where the first model with straight terms left
  40% against 29% curved, that nano stood above the floor at its point value
  only and not above the one-experiment Gria1 map, and the tests' p values of
  section 5.7 against its leftover. Run the same day: 44% left (35% to 54%), 27
  points above the floor (95% +15 to +40), every check row above its own floor
  (section 4.4); Cacng8 follows the new leftover (p 0.0066), the family does not
  (p 0.92), and 53 genes pass BH over all, 50 of them with a negative rho
  (section 5.7).

- **Three reviews of the second version (9 October, evening).** A statistical,
  a figure and a factual review changed how part 1 and part 2 are reported and
  checked, not the model, the density genes, the structures or the tests named
  in advance; every number they touch was computed after the second version's
  run had been seen:
  - the floor's error is stated both ways, with the known map also read with its
    own half (1%, the animals' noise), where decision 4 had said it errs low;
  - every interval of part 1 is also given over spatial blocks of neighbouring
    structures (`jackknife_blocks.csv`, `calibration_jackknife_blocks.csv`,
    `check_rows_jackknife_blocks.csv`), beside the jackknife over structures;
  - the curved model is called a check, not the bound decision 5 had called it:
    terms bent to fifth powers leave less (36%);
  - `psd_pc1`'s floor builds psd_pc1 from one list of genes for both halves,
    those with two halves of experiments, so none is measured once (its margin
    +30.5 points became +30.8);
  - the density genes' agreement with PSD95 is also read inside divisions and on
    the structures every composite has (`comparison.csv`); the quoted 0.79
    belongs to the rule re-run on random halves;
  - the leftover by division (`residual_by_division.csv`) and how the genes past
    BH on the leftover stand to the autofluorescence map (figure 11s1 C);
  - the smoother nulls of the leftover gained ranges of 20 and 40 mm and count
    the genes past BH;
  - the Cacng8 - Gria1 gap is drawn against the band its two-sided test does not
    pass (±0.195) rather than the null's central 95%;
  - figure 14 repeated 03 D and two rows of 03s1, so it is retired; its detailed
    version is 03s4, and the green channel and April's headline are 14 and 15.

### 8.2 To settle

- **The exceptions list.** Glra1 section 61 is kept as true absence, status
  "proposed" in `mapping/ish_section_exceptions.csv`, and 16 sections are kept
  at a step in expression by the rule. To review on `figures\qc\00_flagged.png`
  (red, hatched, dots) and the gene's sheets; the numbers are final once
  reviewed.
- **The genes that follow the straight model's leftover.** 53 genes pass BH
  over all 451 (45 to 59 under smoother nulls), 50 of them with a negative rho,
  most expressed in every neuron (vesicle, ribosomal, cytoskeletal genes;
  section 5.7). They are also the genes the autofluorescence map follows: all 53
  have a positive rho with it, and over every gene the two rhos agree at -0.36.
  Two readings, neuronal mRNA in general (mRNA in somata against receptor at
  synapses, which the floor cannot hold, since both Allen halves share it) and
  the tissue. Tier 3 is exploratory; whether either is worth a named test, for
  example a neuronal-density composite against the leftover, is a decision to
  take before looking further, not after.
- **Control F's range of components.** Its held-out R² is within 0.01 of its
  highest from 30 components on and highest at 40, the edge of the range
  `[beyond_controls] max_pcs` allows (in the first version, on 126 structures, it
  peaked at 21). The control passes either way; a wider range is a choice for
  Giulio.
- **The density term inside divisions.** The rule's genes follow PSD95 mostly
  between divisions (0.08 inside them, held out), the first proposal better
  inside them (0.40). The rule stands as committed and the first proposal is a
  check row (52% left); whether a rule that also asks for agreement inside
  divisions is wanted is a decision for Giulio, to be taken before any nano
  number of it.
- **The glia lead.** The glia set followed the leftover of the four-subunit
  model (q 0.04), a test not named in advance; against the first version's
  leftover it did not pass (q 0.19), against the second's neither (q 0.45). An
  unplanned lead, not pursued now (figure 11s1 says so in one line).
- **The presynaptic set.** GO annotates Syn1, Vamp2, Stx1a, Bsn, Syt1, Cplx1,
  Camk2a and Slc17a7 to both sides of the synapse, so the "not both" rule keeps
  16 presynaptic genes, mostly cell-type and peptide markers. The sets stay as
  committed; the pre- and postsynaptic genes are drawn as context in figure
  10s1.
- **The localisation test's positive control** fails with the GO control pool
  and works with that of 5 October. Both are shown, with the power check;
  neither turns the localisation result positive.
- **The production run.** The outputs quoted here were made on a full copy of
  the production inputs; steps 13 to 30 run on the production data root after
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
| 13 | `run_structure_set` | the per-brain files; `--recompute` redoes the per-adult channel table; nano on the Allen grid is measured every run | 01 | 30 s (about 1 minute from an empty folder or with `--recompute`) |
| 14 | `run_ish_section_qc` | the Allen API for the repair's experiment lists, cached in `cache\`; `--offline` stops instead; `--sheets` draws the missing QC sheets | QC sheets | 20 s (about 15 minutes with `--sheets` from none) |
| 15 | `run_ish_gene_table` | mygene.info and the GO ontology, cached in `cache\`; `--offline` stops instead | 02, 02s | 5 minutes |
| 16 | `run_ish_spatial_null` | the surrogates are cached; `--recompute` draws them again | | 14 minutes, the calibration |
| 17 | `run_ish_gene_ranking` | the surrogates | 05, 05s, 06, 06s, 07s, 12, 12s | 2 minutes |
| 18 | `run_ish_robustness` | | 13, 13s | 1 minute |
| 19 | `run_ish_divisions --sheets` | | 09, 09s, gene sheets | 3 minutes |
| 20 | `run_ish_gene_sets` | | 10, 10s1, 10s2 | 1 minute |
| 21 | `run_synaptome` | the synaptome of Zhu et al. 2018, downloaded once into `<data>\reference\synaptome\`; `--offline` stops instead | 03s4 | 10 s |
| 22 | `run_density_markers` | the mouse GAF and `go-basic.obo` of GO release 2026-08-05, downloaded once into `<data>\reference\go\2026-08-05\`; `--offline` stops instead; reads no nano value, and stops if the rule chooses other genes than `[density_markers] chosen` | 03s3 | 2 minutes |
| 23 to 26 | `run_beyond_density` to `run_beyond_regression` | `run_beyond_calibration` runs the floor and every check row with its own floor, each over structures and over spatial blocks | | under 1 minute each; the calibration and the check rows 10 minutes |
| 27 | `run_beyond_figures` | the tables of steps 17 and 21 to 26 (step 22's for figure 03 D) | 03, 03s1, 03s2, 04, 11s1 | 1 minute |
| 28 | `run_ish_top_genes --sheets` | the tables of steps 15 to 23 | 07, 08, 11, 11s2, top-gene sheets | 6 minutes, the nulls of the share taken and the smoother nulls |
| 29 | `run_sep_channel_check` | | 14, 14s | under 1 minute |
| 30 | `run_ish_overview` | every step's numbers | 00, 15, 15s, `figures\README.md` | 10 s |

A full run of steps 13 to 30 takes about 50 minutes from an empty folder (the
download cache copied in), and a second run gives the same tables.

Every run prints the data root and the settings in force. `SEP_DATA_ROOT`
moves the data root, for a copy. The settings are `[structures]`, `[ish]`,
`[ish_qc]`, `[ish_analysis]`, `[ish_panel_test]` (the label permutations, and the
structures and reliability a gene needs, which `gene_sets` and `top_genes`
read), `[spatial_null]`, `[density_markers]`, `[beyond]`, `[beyond_controls]`,
`[beyond_calibration]`, `[beyond_figures]`, `[beyond_regression]`, `[top_genes]`
and `[ish_figures]` of `mapping/settings.toml`. The tests: `cd mapping`,
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
| `tables\adult_per_mouse.csv`, `adult_profile.csv`, `adult_allen_grid.csv` | 13 | per adult and structure, voxels and log2 means of nano, autofluorescence and SEP, plain and eroded, and `zref`; the cohort profile; per adult and structure, nano on the Allen 200 um grid (voxels, mean, cref, zref) |
| `tables\experiments.csv`, `section_qc.csv`, `experiment_qc.csv` | 14 | every experiment; every section judged, with its flag and both references; per experiment, the sections set missing and kept, its median level and whether it is near zero |
| `tables\gene_region_table.csv`, `gene_table.csv`, `gene_profiles.csv`, `gene_documentation.csv` | 15 | per experiment and structure, ISH means; per experiment, labels, QC and reliability; per gene and structure, the merged rank profile; the supplement table (UTF-8 with BOM) |
| `tables\surrogates_nano.npy`, `surrogates_auto.npy`, `surrogate_structures.csv`, `variogram.csv`, `null_calibration.csv` | 16 | the surrogates (10,000 x 204), their structures, the variograms, the calibration |
| `tables\gene_ranking.csv`, `gene_ranking_per_adult.csv`, `gap.csv`, `null_rho.npz` | 17 | per map and gene: rho, spatial p, q, null band, spread over adults, ranks; per adult; the gap with both nulls and the lead either way its two-sided test needs; every gene's rho with every surrogate, and the gap of every null map |
| `tables\ranking_robustness.csv`, `robustness_summary.csv` | 18 | per variant and gene; per variant |
| `tables\within_division.csv`, `within_division_detail.csv`, `within_calibration.csv` | 19 | per gene, division-only and within rho with both nulls; per gene and division; the two nulls on random maps |
| `tables\gene_sets.csv`, `set_tests.csv`, `contrasts.csv`, `localisation_test.csv`, `localisation_summary.csv`, `localisation_power.csv` | 20 | set members; set tests; contrasts; partial rho per gene and pool; every test of the localisation design; its power by effect size |
| `synaptome\samples.csv`, `density.csv`, `coverage.csv`, `agreement.csv` | 21 | per sample of the synaptome, its ids, densities and structure or why none; per structure of the adult table, measured or why not, its units, weights and densities; per division, how many declared and fitted structures are measured; each density's Spearman with the density term, Gria1 and the maps |
| `density_markers\candidates.csv`, `excluded.csv`, `agreement.csv`, `structures.csv`, `halves.csv`, `selection.csv`, `comparison.csv` | 22 | every gene of the gene table against the rule of the density term, in the pool or why not, with its postsynaptic and plasticity terms; every reason a gene is excluded, with the term, the gene's own term below it and the evidence; the eligible genes by Spearman with PSD95 density, the pool's ranks, the genes chosen; the declared structures, measured for the model or why not; per random half, the genes chosen and each composite's agreement on the other half, also on the structures every composite has and inside divisions; how often each gene is chosen; each composite's genes and their agreement with PSD95, held out and on the full set, whole brain, on the shared structures and inside divisions |
| `beyond\` | 23 to 27 | `structures_used.csv` (with whether PSD95 is measured), `variance_partition.csv` (every model's share held out), `partition.csv` (the four parts), `weights.csv`, `jackknife.csv` and `jackknife_blocks.csv` (the shares, parts and weights on each subsample, at random and one spatial block left out), `calibration.csv` (nano, the floor, and the floor's map read with its own half), `calibration_jackknife.csv`, `calibration_jackknife_blocks.csv`, `check_rows.csv` (each check row: structures, terms, shares, left, its floor and nano minus it with both intervals, what is measured once, the genes behind psd_pc1), `check_rows_calibration.csv`, `check_rows_jackknife.csv`, `check_rows_jackknife_blocks.csv`, `controls.csv`, `gene_space.csv`, `gene_space_calibration.csv`, `gene_space_summary.csv`, `folds.csv` (the main model under other folds), `regression_table.csv`, `residual_by_structure.csv`, `residual_by_division.csv` (the mean leftover per division, 95% over resampled adults), `replication.csv`, `leftover_genes.csv`, `leftover_sets.csv`, `numbers_for_the_caption.txt` |
| `tables\top_genes.csv` | 28 | per gene characterised (past the map's null, Gria1, Cacng8, the family): why it is there and its tier, its GO terms and gene sets, its rho, p, q and rank with the map, inside divisions, over the robustness variants, its reliability, its rho with each term of the main model, its prediction and PSD95, its rho and p with the leftover (q over all genes and within the family), its matched control, and what it takes of the leftover with both nulls |
| `tables\named_tests.csv`, `family_members.csv`, `leftover_null_check.csv` | 28 | the tests named for the leftover, tier 1 and the two group tests of tier 2, on the leftover and on the map, and the two check rows of tier 2 (with how many controls the second changes); every member of the family with its sources, tested or why not; the named tests and the counts of genes below p 0.05 and past BH against smoother nulls, with each null's smoothness |
| `green_channel\sep_channel_check.csv` | 29 | per adult, each channel's range and correlations |
| `tables\april_headline.csv`, `april_anova.csv`, `april_groups.csv` | 30 | P9's genes then and now; the ANOVA under each choice; today's groups against the null |
| `tables\numbers_<step>.csv`, `numbers_for_the_text.csv` and `.txt` | 13 to 30 | the numbers of each step, and all of them |
| `cache\` | 14, 15 | Allen experiment lists, mygene records, `go-basic.obo` |

| figure | drawn by |
|---|---|
| `figures\00_overview.png`, `15_april_headline.png`, `15s_april_headline_detail.png`, `README.md` | `run_ish_overview.py` |
| `figures\01_structures.png` | `run_structure_set.py` |
| `figures\02_genes.png`, `02s_genes_detail.png` | `run_ish_gene_table.py` |
| `figures\03_beyond.png`, `03s1_beyond_budget.png`, `03s2_beyond_controls.png`, `04_beyond_where.png`, `11s1_leftover_genes.png` | `run_beyond_figures.py` |
| `figures\05_one_comparison.png`, `05s_one_comparison_detail.png`, `06_spatial_null.png`, `06s_spatial_null_detail.png`, `07s_gene_ranking.png`, `12_autofluorescence.png`, `12s_autofluorescence_detail.png` | `run_ish_gene_ranking.py` |
| `figures\07_top_genes.png`, `08_cacng8_gria1.png`, `11_leftover.png`, `11s2_ampa_family.png`, `top_genes\<gene>.png` | `run_ish_top_genes.py` (`--sheets`) |
| `figures\09_between_within.png`, `09s_between_within_detail.png`, `genes\<gene>.png` | `run_ish_divisions.py` (`--sheets`) |
| `figures\10_gene_kinds.png`, `10s1_gene_sets.png`, `10s2_localisation.png` | `run_ish_gene_sets.py` |
| `figures\13_robustness.png`, `13s_robustness_detail.png` | `run_ish_robustness.py` |
| `figures\03s4_synaptome_detail.png` | `run_synaptome.py` |
| `figures\03s3_density_markers.png` | `run_density_markers.py` |
| `figures\14_green_channel.png`, `14s_green_channel_detail.png` | `run_sep_channel_check.py` |
| `figures\qc\00_flagged.png`, `qc\<gene>_<experiment>.png` | `run_ish_section_qc.py` (`--sheets`) |

Reference data, under `<data>\reference\`, each folder with a `fetch_log.txt`
(source, version, checksums, what was read):

| folder | what | source |
|---|---|---|
| `synaptome\` | PSD95 and SAP102 punctum density, 37 subtypes in 775 samples (a region of one hemisphere in one section), each subtype divided by its largest value as shared (Zhu et al. 2018); written by `run_synaptome.py` | github.com/netneurolab/hansen_synaptome, `data/synaptome/mouse_liu2018/` at commit `0399525412b6f50cfdeb5904b96da7fa8e4b507c`; archived as Zenodo doi 10.5281/zenodo.18201390 |
| `go\` | the GO Consortium's mouse annotation, and its lines for GO:0032281, read for the AMPA receptor complex family | `current.geneontology.org/annotations/mgi.gaf.gz`, fetched on 8 October while the current release was 2026-08-05; the file's header is dated 2026-05-21 (MGI's source of 2026-05-20). Release 2026-08-05's own mouse GAF (`go\2026-08-05\`) lists the same 31 genes for GO:0032281, which a test checks |
| `go\2026-08-05\` | the mouse annotation (MGI's GAF, `MOUSE-mod.gaf.gz`, its header dated 2026-08-04) and `go-basic.obo` (data-version 2026-07-26) of GO release 2026-08-05, read by `run_density_markers.py`, each checked against its SHA-256 | `release.geneontology.org/2026-08-05/` (`annotations/gaf/`, `ontology/`) |
| `ampar_complex\` | Schwenk et al. 2012: Figures 1 to 6 and the supplement (Tables S1 to S4) | the publisher's file server, `ars.els-cdn.com` |

References: Ashburner M, Ball CA, Blake JA, Botstein D, Butler H, Cherry JM, et
al. (2000). Gene Ontology: tool for the unification of biology. Nature Genetics 25,
25-29. Burt JB, Helmer M, Shinn M, Anticevic A, Murray JD (2020).
Generative modeling of brain maps with spatial autocorrelation. NeuroImage 220,
117038. Freedman D, Lane D (1983). A nonstochastic interpretation of reported
significance levels. Journal of Business and Economic Statistics 1, 292-298. The
Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. Genetics
224, iyad031. Hansen JY, Luppi AI, Qiu Z, Gini S, Fulcher BD, Gozzi A, et al. (2026). Synapse
types are spatially associated with regional hemodynamics in the mouse brain.
PLOS Biology 24, e3003637. Lein ES, Hawrylycz MJ, Ao N, Ayres M, Bensinger A,
Bernard A, et al. (2007). Genome-wide atlas of gene expression in the adult mouse
brain. Nature 445, 168-176. Rouach N, Byrd K, Petralia RS, Elias GM, Adesnik H,
Tomita S, et al. (2005). TARP γ-8 controls hippocampal AMPA receptor number,
distribution and synaptic plasticity. Nature Neuroscience 8, 1525-1533. Schwenk
J, Harmel N, Brechet A, Zolles G, Berkefeld H, Müller CS, et al. (2012).
High-resolution proteomics unravel architecture and molecular diversity of
native AMPA receptor complexes. Neuron 74, 621-633. Zhu F, Cizeron M, Qiu Z,
Benavides-Piccione R, Kopanitsa MV, Skene NG, et al. (2018). Architecture of the
mouse brain synaptome. Neuron 99, 781-799.
