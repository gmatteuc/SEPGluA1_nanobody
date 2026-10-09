# Scientific context

What this code is for: the questions, what each line of work found, the grant
it now serves, and the papers it builds on. The papers are in
`data\ref_papers\`, with summaries in `PAPER_SUMMARIES.md` there.

A path under `data\` is relative to the data root: the folder `data` beside
the code folder (`D:\sep_histology\data` beside `D:\sep_histology\code` on
the lab computer), or the folder the environment variable `SEP_DATA_ROOT`
names (`get_paths.m`, `mapping/sepmap/config.py`). Each result below comes
from the output file named beside it, or from a document in this repository;
the numbers of the papers come from the papers.

The results are those of the output files of 5 October 2026, when the
Python route was rerun with the merged code of the refactor. Its outputs
before that rerun, the state the grant figures were made from, are kept in
`data\backup_before_rerun_2026-10-05\` and, as of 30 September, in the
snapshot `G:\sep_histology_snapshot_2026-09-29` ([ROADMAP.md](ROADMAP.md),
section 2). Some results will move when the scientific additions A1 to A5
run on one declared set of structures (S1 to S5 and A1 to A5 in
[history/REFACTOR_PLAN.md](history/REFACTOR_PLAN.md)); those are marked
provisional. Where a fix of step 8 or a decision of 4 October moved a quoted
number, the earlier value is given beside it. Line of work 3, and the
reproducibility of line 2, quote the ISH analysis built on 8 October 2026 on
branch `post-ish` (`adult_v2\ish_analysis\`, [ISH_ANALYSIS.md](ISH_ANALYSIS.md)),
run on a full copy of the production inputs: the production data root does not
hold those outputs yet; it does once the branch is merged and steps 13 to 30
of the Python route run there.

## The question

AMPA receptors carry most fast excitatory transmission in the brain. They are
assembled from combinations of four core subunits, GluA1 to GluA4
(Lopez-Ortega et al. 2024). Regulating their function and their traffic to and
from the membrane is critical for many forms of synaptic plasticity, and tagged
GluA1 is recruited to spines after long-term potentiation (Huganir & Nicoll
2013). A sizeable surface pool of receptors is a recurring requirement of
current models of LTP (Huganir & Nicoll 2013).

The grant this project serves uses surface GluA1 as a molecular readout of
synaptic remodelling potential: a proxy for plasticity potential, not a
definitive marker of critical periods.

The project maps surface GluA1, labelled with a nanobody in tissue that was
not permeabilised (below), across the whole mouse brain, and asks four
questions of the map:

1. Where does an experience change it?
2. How is it distributed in the adult brain, and how reproducible is that
   across mice?
3. What does it measure: more than Gria1 expression or synapse density?
4. How does it differ between young and adult mice, and could the
   differences mark critical periods?

## The starting point: the method

- **Mice.** SEP-GluA1 knock-ins: the GluA1 subunit is fused to SEP, a
  pH-sensitive fluorophore that reports surface-exposed receptors in living
  cells (the grant, citing Graves et al. 2021, eLife). A mouse's age is the
  one in the name of its raw data folder (`MG904_SepGluA_P22` is P22). The
  adults (folders ending in `_Gria1`) were well above P60; their exact ages
  are not recorded yet.
- **Labelling.** The sections are fixed and were not permeabilised: no
  detergent, or very little, so the label should stay on the membrane, with
  little access to intracellular receptor pools. The GFP part of SEP-GluA1 is
  amplified with a GFP-booster nanobody. In preliminary preparations this
  showed surface GluA1 with no detectable intracellular labelling (grant,
  Fig. 3a). The map is therefore surface GluA1 in principle; no specificity
  or total-receptor control has tested it on these brains (line 3). The code
  calls its channel nano. The sections are not cleared.
- **Channels.** Each section is recorded in four channels: DAPI (nuclei), the
  nanobody (Cy5, called nano), autofluorescence (Cy3, called auto) and the
  green channel of SEP itself (filter EGFP).
- **Registration.** Serial sections are registered to an atlas with
  LightSuite (elastix, affine and B-spline), driven by the DAPI channel and by
  control points between each section and the atlas. The points are placed by
  hand, or proposed by the automatic annotation and reviewed by hand. Adults
  are registered to the Allen CCF, young brains to the DeMBA atlas of their
  own age (decision of 2 September 2026).
- **Two analysis routes.** The plasticity comparison stays in MATLAB
  (`group_comparison/`). The adult map, the ISH comparison and young against
  adult are in Python (`mapping/`). They normalise differently: the Python
  route scales each brain on its own; the plasticity comparison fits each
  brain to its group's median and the experimental group onto the control
  group. If analyses of both routes go into one paper, the two must be made
  comparable, or the difference discussed.

## The project: four lines of work

### 1. Plasticity localisation (the original aim; paused)

- **Question.** The project began as a search for where plasticity lands
  after an experience: surface GluA1 in mice after rhythmic whisker
  stimulation (RWS, the protocol of Gambino et al. 2014) or after behaviour,
  against naive mice, in cortex or elsewhere.
- **Mice.** The approved comparison is the run of 5 December 2025: five naive
  mice (CGF027, CGF028, CGF033, CGF034, CGF035), five RWS mice (MG691, MG692,
  MG693, MG736, MG737) and four behaviour mice (MG705, MG709, MG716, MG718) of
  the seven registered.
- **Protocols.** Each RWS mouse received one RWS session and was then
  perfused. Each behaviour mouse was perfused on the first day it performed
  above chance in a single-whisker detection task. The naive mice had
  neither.
- **Method.** For every mouse the hemispheres are folded (left plus mirrored
  right, and the absolute left-right difference). The experimental group's
  intensity profile is aligned onto the control group's. The groups are then
  compared voxel by voxel, with Welch t and surprise (-log10 p) maps.
- **Result.** The headline is an increase of the nanobody signal in S1 after
  RWS, in the hemisphere-sum t map of the slab around plane 565, masked at
  p < 0.01. That threshold is uncorrected and applied voxel by voxel, on maps
  smoothed in 3D (Gaussian, sigma 5 voxels), with five mice per group; the
  figure shows the median over planes 555 to 575. The increase is small. It is
  where plasticity caused by whisker stimulation was expected. The outputs
  are in `data\comparisons\naive_vs_rws\` and `naive_vs_behavior\`.
- **Reproduced.** On 1 October 2026 today's group-difference step was rerun
  on the normalised volumes of 26 and 27 November 2025, the inputs of the
  approved run. The slab t maps and surprise masks correlate with the
  approved ones at 0.997 to 0.998 for RWS (0.988 to 0.999 for behaviour), the
  regional bars at 0.994 to 0.998, and the S1 increase is there
  ([history/REFACTOR_PLAN.md](history/REFACTOR_PLAN.md), Progress,
  1 October). The small residue comes from the background masks, which were
  regenerated since. This checks the computation on the same mice; it is not
  a replication in new animals.
- **A change to the smoothing.** The smoothing worked over the tissue only,
  but then set the voxels outside it to zero, so the zeros entered the group
  means and the n of the SEM. Fix 23 of step 8 leaves them out, which changes
  the approved comparison. Rerun on the same December 2025 inputs, the RWS
  S1 increase holds: the hemisphere-difference map keeps the same 292
  significant pixels in the barrel field, the hemisphere-sum map 287 of its
  437, and per mouse the barrel field is 26% higher after RWS (p = 0.12, five
  against five), as before. The behaviour comparison changes most: MG709,
  with no tissue at that slab, no longer counts as a fourth mouse
  ([ROADMAP.md](ROADMAP.md), section 1).
- **Edge voxels and the bars (4 October).** A voxel now gets a t only where
  each group has at least three mice with tissue, and each regional bar is
  the share of one of 71 atlas regions' voxels at p < 0.01, no voxel in two
  bars. On the same inputs the barrel field keeps its 292 and 287 pixels and
  is the only significant spot of the hemisphere-difference map in the slab;
  the large t values at the tissue edges are gone (ROADMAP, section 1). The
  merged code's run on production's current normalised volumes (5 October)
  is in `data\comparisons\naive_vs_rws_nano\` and `naive_vs_behavior_nano\`.
  Those volumes are not the approved run's inputs, so that run is not the
  reference; its region shares agree with the rerun on the approved inputs
  (r 0.994 to 0.997 for RWS).
- **Status.** The team was not fully confident that the effect holds, and
  more animals would be a large investment, so the line is paused, not
  closed. The code stays runnable, and documented well enough to resume with
  more animals.

### 2. The adult distribution

- **Question.** How surface GluA1 is distributed across the whole adult
  brain, and how reproducible that is across mice.
- **Cohort.** The five naive and five RWS mice, pooled. Pooling was checked:
  the part of the map that Gria1 expression and synapse density do not
  predict (line 3) agrees between the naive and the RWS mice at rho +0.877
  (`adult_v2\ish_analysis\beyond\controls.csv`, control D). The behaviour mice
  enter only the plasticity comparison.
- **Reproducibility.** Over the 126 ways of splitting the ten adults into two
  halves of five, the two half-cohort maps agree at rho 0.974 over 126
  grey-matter structures
  (`adult_v2\ish_analysis\beyond\numbers_for_the_caption.txt`; the structures
  in `beyond\structures_used.csv` there). Spearman-Brown takes this to 0.987
  for the full cohort. The autofluorescence of the same brains predicts none
  of the map (-2% of its reproducible part, `beyond\variance_partition.csv`).
- **The structures** every adult analysis of the ISH line reads: the 204
  grey-matter structures measured in all ten adults (A1, S1;
  `adult_v2\ish_analysis\tables\structure_set.csv`), with `zref` taken over
  them ([ISH_ANALYSIS.md](ISH_ANALYSIS.md), section 3).
- **Being rebuilt.** Which structures stand out, and with what confidence, is
  being redone. The MATLAB route (`run_characterize_distribution`) called a
  structure enriched above a threshold on its own scale, and Sami El-Boustani
  asked (28 April) for a stricter threshold and a permutation null. The plan
  replaces both with a test per structure against the brain's median
  structure (S2), on a declared set of structures (S1, A1, A4). Today 22
  structures of the adult table are seen in fewer than five of the ten adults
  (history/REFACTOR_PLAN.md, S1). Its outputs carry known defects (listed in
  the plan) and are not quoted here.

### 3. What the map measures

A reader will ask whether the map simply follows how much receptor a region
makes, or how many synapses it has. The project tests this against the Allen
in situ hybridisation (ISH) maps of the adult mouse brain (Lein et al. 2007).
The analysis, its figures and every number are in
[ISH_ANALYSIS.md](ISH_ANALYSIS.md) (built on 8 October 2026: A1 to A3 and A6
to A9 for the ISH line; part 1 run again on 9 October with its main model, and
the same day in its second version). The
outputs are in `adult_v2\ish_analysis\`, the numbers in its
`tables\numbers_for_the_text.csv`. Its figures follow the argument,
`figures\00_overview.png` (the line on one page) to `16_april_headline.png`,
each main figure of two to four panels beside its detailed version (the same
number with an s).

**The map is not fully explained by Gria1 expression and synapse density**
(part 1 of the argument; figures 03 and 04). This is the second version of part
1, decided by Giulio on 9 October 2026 after seeing the first (ISH_ANALYSIS.md,
section 8.1):

```
nano map rank  ~  Gria1 rank + synapse-density rank
```

Gria1 mRNA is abundance (the stained protein is GluA1, which Gria1 alone
encodes; Gria2 to Gria4 make partners the nanobody does not see, whose
availability sets GluA1's assembly and trafficking). Synapse density is the mean
rank of three postsynaptic genes, Rock2, Cap2 and Slc8a2, chosen without the map
by their agreement with the measured PSD95 punctum density of Zhu et al. (2018),
after every gene that places or regulates AMPA receptors was left out as the
surface side (held out on random halves of the structures, the three agree with
PSD95 at 0.79). The two terms enter straight, so a map made of them, the
calibration floor, holds no curvature the model could miss. Over the 163
declared structures where Gria1 and the three genes are measured (40 lack
Gria1, mostly hypothalamic and midbrain, one Rock2), the model is scored on
held-out structures (20 shufflings of five folds) against the map's own
reliability: two halves of the cohort agree at 0.972, so 98.6% of the map is
reproducible (Spearman-Brown for ten adults).

| held out | share of the reproducible map (95% over structures) |
|---|---|
| Gria1 alone | 41% (27% to 55%) |
| synapse density alone | 46% (35% to 56%) |
| both: Gria1 only 10%, shared 31%, density only 15% | 56% (46% to 65%) |
| left | 44% (35% to 54%) |

- The weights, map and terms z-scored: Gria1 +0.39 (+0.23 to +0.54), synapse
  density +0.47 (+0.32 to +0.61).
- The leftover replicates across halves of the cohort at 0.933; the map's
  reliability and the fit alone imply 0.938, so the replication is expected once
  the map replicates.
- The calibration, on the 131 structures where both halves of each gene's Allen
  experiments measure Gria1 and the three genes: a map made only of Gria1 and
  synapse density, predicted from the other half, leaves 15% (the floor that
  Allen-to-Allen mismatch produces); the nano map, read the same way on the same
  structures, leaves 42%, 27 points more (95% +15 to +40, the difference
  resampled with the structures). Every gene of the model has two or more Allen
  experiments, so the floor holds the whole mismatch of one Allen map with
  another; it does not hold that of Allen's P56 mice with these brains (age,
  strain, grid, registration), so it errs low. The first version's second
  benchmark, a map of one Allen Gria1 experiment, is gone: its made-up animals
  share one Allen brain's quirks, which then count as reproducible, while the
  nano map averages ten brains.
- Check rows, each against its own floor (nano minus the floor, in points): the
  model curved leaves 38% (+20, 95% +8 to +32), with autofluorescence 43% (+23),
  with the first proposal for density (Dlg4, Homer1, Camk2a) 52% (+30), with the
  11 marker genes of the first version 49% (+25), with `psd_pc1` 43% (+30), with
  Gria1 to Gria4 28% (+13, +4 to +22); with the measured PSD95 density as the
  density term 50% on its 89 structures (+31), against 50% for the density genes
  there (+27); on the 118 structures of 0.4 mm³ or more 49% (+24); with nano
  measured on the Allen 200 um grid 42% (+29). Every row stays above its own
  floor; curvature and the partner subunits take the most.
- Under other folds: spatial blocks leave 53% (on the calibration's structures
  nano 48%, the floor 15%); one shuffling of the folds 43%, ten folds and leave
  one out 44%.
- Seven controls try to break it (a spatial gradient, structure size, single
  animals, naive against RWS, curvature, the choice of predictors, the reading);
  six pass (figure 03s2; `adult_v2\ish_analysis\beyond\controls.csv`). Control E
  does not: curving the two terms raises the held-out R² from 0.551 to 0.615, so
  part of the straight leftover is curvature, and the curved check row is its
  bound. The components of the 209 genes measured in every structure predict
  81%, and the nano map stands above that model's own floor (22% to 24% left
  against 3% to 5%). So the claim is "not predicted by Gria1 expression or
  synapse density", never "beyond gene expression".
- The leftover is highest relative to prediction in the triangular nucleus of
  the septum (+91 ranks), the septofimbrial nucleus, the medial geniculate, the
  indusium griseum and the lateral geniculate, and lowest in the preoptic and
  paraventricular hypothalamic nuclei, the rhomboid nucleus, the nucleus of
  reuniens and the arcuate nucleus (`beyond\regression_table.csv`); on the planes
  the isocortex sits below prediction and the hippocampal formation above. A
  claim about one structure needs its own null.
- **What it means.** About 44% of the map's reproducible pattern is not
  predicted by Gria1 expression or synapse density, reproducibly across mice,
  and on the same structures it is well above what Allen-to-Allen mismatch would
  leave. The share depends on the model (curved 38%, with the partner subunits
  28%), and the floor errs low by an unmeasured amount. The reading the data
  support is the surface fraction of the receptor, shaped by trafficking
  regulation and scaffolding. It is an interpretation, not a measurement:
  translation, turnover, subunit composition or nanobody access would land in the
  same leftover, and only a total-receptor measure on the same brains separates
  them.

**The gene analyses are consistent with that reading, and do not single it
out** (part 2; figures 05 to 11). Each gene's rho with the map is tested
against 10,000 surrogate maps with the map's smoothness (Burt et al. 2020;
A7), which give a spatial p below 0.05 for 4.0% of random maps where the
ordinary Spearman p gives 62%.

- 12 of 451 genes pass after BH. Cacng8 (TARP γ-8) is first, +0.807, and first
  of P9's genes in every robustness variant (statistic, borders, reading,
  inputs, structure set; A3; under Pearson on log2 it is 2nd of all genes).
  Gria1 passes too, +0.646, 32nd of all genes (11th of P9's 100): S5 is met,
  its rank is quoted with its null, and it moves with the choices (5th to 17th
  of P9's genes, 18th to 47th of all).
- The Cacng8 - Gria1 gap is +0.167 on the structures both have, above zero in
  every adult. The test fixed in advance is two-sided, against maps related to
  both genes alike, and the gap is inside it (p 0.105; against unrelated maps, a
  conservative bound, p 0.348); the null is skewed, and 2.3% of those maps give
  a Cacng8 lead as large, 8.2% a Gria1 lead. These maps do not settle which of
  the two the nano map follows more closely.
- The tissue's autofluorescence, read and tested the same way, passes for
  neither Gria1 (+0.176, p 0.447) nor Cacng8; adult by adult, nano's rho with
  Gria1 is 0.482 to 0.754 and autofluorescence's -0.266 to +0.297 (A8). Gria1
  and Cacng8 follow the label, not the tissue; the tissue has a gene pattern
  of its own (the two gene orders agree at 0.363).
- Inside divisions, where no cortex-against-thalamus contrast can enter, 11
  genes pass, Cacng8 (+0.548), Dlg2 (+0.513) and Gria1 (+0.417) among them
  (A6).
- Gene sets fixed in advance from GO and cited marker lists: no set passes its
  null after BH. Of the two contrasts named in advance, postsynaptic above the
  presynaptic set passes (+0.274, spatial p 0.016) and above glia does not
  (+0.336, p 0.170), so the postsynaptic criterion (both) is not met; the
  presynaptic set here is mostly cell-type and peptide markers. Once the
  subunit composite is removed, the 84 AMPA receptor localisation genes
  predict the map no better than expression-matched postsynaptic genes
  (+0.0041, p 0.916); the test finds a planted localisation advantage of
  +0.077 in 80% of maps, so there is none larger than that. Its positive
  control comes out only with the control pool of 5 October (p 0.0117), where
  localisation is no better either (+0.0073, p 0.839).
- Against the leftover of the main model (surrogates put through the same
  fit), Cacng8, the gene named for it, follows it: +0.200, spatial p 0.0066, and
  so under smoother nulls (p 0.0018 to 0.0030). It is a re-test: its p against
  the four-subunit leftover (0.0024) and against the first version's (0.0002)
  were seen before. Added to the main model as one more term, Cacng8 takes 9.8
  points of the reproducible map from the leftover: more than its plain
  surrogates take (p 0.001), the null fixed first, and more than maps that
  relate to the model as it does (p 0.003), the secondary line.
- The AMPA receptor complex family named with Cacng8 (the native complexes of
  Schwenk et al. 2012 and GO:0032281, with Gria2 to Gria4, without Gria1; 31 of
  37 genes with a usable map) does not follow the leftover (median -0.004, p
  0.92 against its surrogates; +0.063 above 31 postsynaptic genes of the same
  expression, p 0.29); the three members past BH within it run against it
  (Gria4 -0.52, Olfm2 -0.45, Olfm1 -0.21), and Cacng8 misses that BH by a hair
  (q 0.051).
- Over every gene, an exploratory reading: 53 of 451 pass BH against the
  leftover, 50 of them with a negative rho, most of them expressed in every
  neuron (vesicle genes such as Syp, Syt1, Snap25; ribosomal proteins; Mapt): the
  map sits below what Gria1 and the three postsynaptic genes predict where
  neuronal mRNA in general is high. The glia set, which passed against the
  four-subunit leftover (q 0.04, not named in advance), does not (q 0.45); it
  stays an unplanned lead, not pursued.
- The 12 genes past the null are maps much like Gria1 and synapse density (rho
  0.58 to 0.90 with the main model's prediction), so they describe the map
  rather than add evidence; five of them follow the leftover before correction,
  Cacng8 the one named for it (ISH_ANALYSIS.md, section 5.3; figure 07).
- April's headline, P9's category violins with ANOVA p 0.032, rested on a
  split written after looking (p 0.20 without it) and on genes treated as
  independent draws: today's F across April's groups, against the F of the
  surrogates, gives p 0.289. The gene order reproduces (0.959; figure 16).

**The green channel is not total receptor** (figure 15).

- **The plan.** The SEP tag fluoresces green, so the green channel was meant
  to report all SEP-GluA1, surface and internal, and nano divided by it a
  surface fraction.
- **Measured** on the 204 declared structures
  (`adult_v2\ish_analysis\green_channel\sep_channel_check.csv`). In every adult
  the green channel follows autofluorescence (rho 0.723 to 0.827), while nano
  follows it at only 0.047 to 0.529, and the green channel varies across
  structures about as little (p90 - p10 of log2: 0.915 against 0.999; nano
  1.967). Against Gria1: nano +0.619, the green channel +0.333.
- **Consequences.** `sepratio` is not a surface fraction, and the three-way
  test of September cannot be run on these data; the tests of the channel
  ratios against genes (`ish/arms.py`) and the ratio table behind them
  (`adult/arms.py`) are retired to `archive/`. What their partial
  correlation hinted at, localisation genes predicting the map once Gria1 is
  removed, is the localisation test above, which is negative.
- **What would settle it.** A measure of total receptor in the same brains: a
  total-GluA1 antibody stain after permeabilisation, or autoradiography, on a
  subset of them. A knockout or no-primary control would test specificity
  instead. There is no knockout or competition control, no independent
  regional measure and no DAPI control yet.

### 4. Young against adult (grant, September 2026)

- **Question.** Young mice (P16 to P36) against adults, looking for
  differences in cortex that could correspond to critical periods: the
  grant's Aim 2.1 (below).
- **Brains.** Seven young brains so far: five at P20 (MG897, MG903, MG909,
  MG910, MG913), one at P16 (MG911) and one at P22 (MG904), against the ten
  adults. MG912 (P20) was dropped for poor slice quality. Six more, MG896 and
  MG914 at P28, MG906 and MG908 at P32, MG895 and MG907 at P36, are in
  processing and go through the refactored code.
- **Atlases.** Each young brain is registered to the DeMBA atlas of its own
  age, the adults to the CCF. Region statistics are computed in each brain's
  own atlas, whose labels share the CCF ontology, so they need no warping.
  Maps carry each young brain to the CCF with the DeMBA deformation of its
  age.
- **Only the pattern can be compared.** Young and adult brains were imaged in
  different sessions. The autofluorescence used as an internal standard rises
  with age itself. Nano per unit autofluorescence is lower in young cortex
  (-0.42 log2 in retrosplenial to -1.22 in frontal cortex, the `ratio` rows
  of `comparisons_v2\young_vs_adult\group_stats.csv`), but that is a bound,
  not a value.
- **The young brain is flatter.** The spread of its structures (p90 - p10 of
  log2 nano relative to cortex) is about half the adult one: 0.99 ± 0.27 log2
  in the seven young brains against 1.93 ± 0.26 in the ten adults (mean ± SD
  over brains of the spread the maps divide by: `z_spread` in the per-brain
  caches `comparisons_v2\per_mouse\<mouse>_scalars.npz`, which `run_cohort`
  rewrites, over each brain's structures; the ISH analysis's
  `green_channel\sep_channel_check.csv` takes the adults' spread over the 204
  declared structures instead, 1.97). The region tables take
  the spread over the structures every brain shares, so their values differ
  slightly. `zref`
  removes this per brain, so a zref difference is a difference of positions
  within each brain's own range, not a fold change (see Terms).
- **By system** (7 young against 10 adults, `cref`, log2 young minus adult,
  group medians from `group_stats.csv`; Mann-Whitney p < 0.01, uncorrected,
  unless marked): primary somatosensory +0.22, retrosplenial +0.24, frontal
  -0.43, striatum -0.79, hippocampus -1.03; primary visual -0.06, not
  significant. The reading matters: hippocampus is -1.03 under `cref` and
  +0.49 under `zref`.
- **RL and AL.** RL and AL together, the visuo-tactile areas of the grant,
  are the only visual group that survives correction under `zref`: +0.24
  (p = 0.0007, q = 0.014), against -0.01 between naive and RWS adults; the
  P20 brains alone give +0.34. This holds under `zref` only: under `cref`,
  which is blind to a shift of the whole cortex, RL+AL is +0.08 (p = 0.36).
  V1 shows no difference (`zref` +0.07, q = 0.33).
- **By layer.** Under `cref`, in every sensory system, primary and higher
  order alike, the young-minus-adult difference is larger in the
  supragranular layers (+0.11 to +0.30) than in the infragranular ones. The
  infragranular layers are lower in young brains in the visual and auditory
  systems (-0.13 to -0.25), but slightly higher in somatosensory cortex
  (SSp +0.07, SSs +0.11). Under `zref`, RL+AL is higher in both the
  supragranular layers (+0.28) and layer 4 (+0.27, q = 0.006 each), while
  layer 4 does not differ in V1 or in the other visual areas (+0.13 and
  +0.15, q ≥ 0.09): RL and AL have a laminar profile closer to somatosensory
  cortex than to V1 (all from `group_stats.csv`).
- **The laminar contrast.** Taken per mouse, the supragranular-minus-
  infragranular contrast would cancel any scale factor of a brain (exposure,
  staining strength, the choice of reference), so `cref`, `subref` and
  `zref` would agree on its sign. No output of the code computes it yet. A
  read-only check on the per-mouse volumes (1 October, not kept as an output)
  found it larger in young brains in most laminar systems under those three
  readings, but not under `ratio` and `sepratio`, whose denominators vary
  between layers. It is not quoted until `region_groups` writes it.
- **A critical period cannot be claimed yet.** Two ages cannot show a window.
  V1, where the classic mouse critical period sits (Levelt & Hübener 2012),
  shows no difference at P16 to P22, while among the sensory systems
  somatosensory cortex shows the largest. And the superficial-layer excess is
  also what the later maturation of the superficial layers alone would give.
  The P28 to P36 brains can separate the two: if V1's superficial excess
  peaks near P28 and falls by P36 while S1 has already declined, that is a
  window moving across modalities; the same decline everywhere at the same
  rate is maturation.
- **Registration caveats.** The P20 and adult atlases differ most next to
  white matter. With the two aligned by their crops alone, 18% of the adult
  fibre tracts fall in the P20 isocortex; the DeMBA deformation that carries a
  young brain to the CCF recovers 94% of them
  (`young\registration_qc\demba_to_allen_transform_qc.png`). So results next
  to white matter need care.
- **Grant figure.** The grant's Figure 3b to d shows this comparison: maps of
  pups (n = 7) and adults (n = 10), a flatmap of the pup-minus-adult
  difference, and Mann-Whitney tests on V1, RL+AL, the other higher visual
  areas, S1, S2 and medial prefrontal cortex. Its caption gives the pups as
  aged P16 to P20: MG904, at P22, was grouped with the P20 brains, as Sami
  El-Boustani asked. The layer flatmaps (`detail_flatmap_layers_zref.png`)
  show three bands. For RL+AL, layers 5 and 6 are the weakest statistically
  (`zref` q = 0.092, against 0.006 for the supragranular layers and layer 4).
  Both groups are high there: the young median is 0.54 and the adult one
  0.35, against 0.38 and 0.10 in the supragranular layers (`group_stats.csv`).
  So the band is bright in both groups' maps, and its young-adult difference
  is the smallest of the three (+0.19, against +0.28 and +0.27).
- **Provisional.** A1 changes the reference of zref and so moves every zref
  value, the differences behind the grant figures included. The code state
  behind the grant figures is tagged `grant-2026-09`. Fix 1 of step 8 takes
  the fibre tracts, the ventricles and the unassigned labels out of subref's
  reference, which moves each brain's subref by one constant (+0.19 to +0.45
  log2 in the adults, +0.02 to +0.32 in the young brains). The young-adult
  subref difference moves by -0.21 log2, and its structures at Welch
  q < 0.05 fall from 104 to 81 with fix 1 alone, and to 83 in the production
  tables of 5 October, with decision 4 as well ([ROADMAP.md](ROADMAP.md),
  section 1); fix 1 leaves the other readings unchanged. Fix 3 of step 8
  stops the cohort mean from counting MG897's 183 voxels without a sepratio
  value as zero, so the young sepratio mean there rises by a factor of 7/6
  and its n falls from 7 to 6; nothing else changes. Decision 4 of 4 October
  leaves unimaged voxels out of each brain's backgrounds and carries the
  young brains to the CCF without darkening their tissue edges. The
  differences quoted in this section move by at most 0.006 (RL+AL `zref`
  +0.232 to +0.235, which now rounds to +0.24) and their q values by at most
  0.02; in the maps the young cortical surface reads brighter than before.

## The grant: SNSF Weave (El-Boustani, Geneva; Gjorgjieva, Munich)

*Dendritic plasticity rules shaping the emergence of multisensory integration
in the developing cortex.*

- **Central hypothesis.** Activity-dependent dendritic plasticity rules
  organise visual and tactile synaptic inputs along dendrites during a
  postnatal window, and so shape the nonlinear integration behind supramodal
  visuo-tactile representations in adult cortex.
- **Three aims.** Aim 1, how dendritic organisation supports multisensory
  integration in adult associative cortex. Aim 2, the developmental sequence
  and plasticity window in which visuo-tactile integration emerges. Aim 3,
  the dendritic plasticity rules that build it. The two-photon imaging
  repository (`D:\dendrites\code`, its `docs/SCIENTIFIC_CONTEXT.md`) holds the
  dendritic and functional imaging side; this repository is Aim 2.1.
- **Aim 2.1.** *Map developmental windows of synaptic plasticity potential
  using surface GluA1.*
  - Rationale: AMPA receptor trafficking is a central mechanism of synaptic
    strengthening, and GluA1 plays an important role in potentiation. Total
    GluA1 immunostaining and Gria1 mRNA do not report the pool of receptors
    at the surface.
  - Design: SEP-GluA1 mice at stages such as P14, P18, P22, P26, P30, P34 and
    adult; non-detergent processing and a GFP-booster nanobody; registration
    with LightSuite, juvenile brains first to age-specific templates; V1,
    whisker S1 and the associative areas RL, AL and LI at area and layer
    resolution; hippocampus, amygdala and striatum as internal controls.
  - Main output: a developmental plasticity-potential index, surface GluA1
    normalised for area, layer, age, section quality, staining batch and
    imaging conditions.
  - Key prediction: associative visuo-tactile areas show a delayed or
    prolonged window of elevated surface GluA1, compared with primary sensory
    cortices.
  - Risks: surface GluA1 may not map one-to-one onto synaptic plasticity, so
    it is read as a proxy for plasticity potential, not a definitive marker of
    critical periods. If fixed sections show too weak a cortical signal, the
    histological atlas would be complemented by in vivo bulk imaging of
    SEP-GluA1 from about P15 to P35.
- **Preliminary data.** Non-detergent labelling with no detectable
  intracellular signal (Fig. 3a), and region-specific pup-adult differences,
  including in RL and AL (Fig. 3b to d, the comparison of line 4).
- **Expected outputs.** Among at least three, a manuscript or resource paper
  on developmental surface GluA1 mapping and plasticity windows. The grant
  states that processed datasets, analysis code and simulation code will be
  released upon publication, when ethical and legal requirements allow.
- **Where the code stands.** The young brains in hand are P16, P20, P22, P28,
  P32 and P36, not the grant's planned series. The code does not compute the
  plasticity-potential index; the readings below are what exists.

## Related literature

The eleven papers in `data\ref_papers\`, grouped by the question they bear
on. For each: what it shows, and what it means here. Papers the grant cites,
but that are not in the folder, are named through the grant.

### What surface GluA1 reports

**Huganir RL, Nicoll RA (2013).** AMPARs and synaptic plasticity: the last 25
years. *Neuron* 80:704-717. `Huganir_2013.pdf`

- Review. Regulating AMPA receptor function and membrane trafficking is
  critical for many forms of synaptic plasticity; GFP-tagged GluA1 is
  recruited to spines after LTP; a sizeable surface pool of receptors is a
  recurring requirement of LTP models. Auxiliary proteins steer the
  receptors: TARPs (type I: γ-2, γ-3, γ-4, γ-8) ensure their delivery to the
  surface and to synapses. Deleting cornichons 2 and 3 causes a selective
  loss of surface GluA1-containing receptors in the hippocampus, a
  selectivity that appears to be mediated by TARP γ-8. Earlier work found
  GluA1 required for LTP; molecular replacement experiments (Granger et al.
  2013) found no specific subunit requirement.
- Here: why surface GluA1 is a plausible plasticity readout, and why not a
  specific one, since LTP may not need GluA1 in particular. It also names the
  proteins behind the top of the ISH ranking: Cacng8 is TARP γ-8, Cnih2 is
  cornichon-2. Since localisation genes as a class explain no more of the map
  than other postsynaptic genes, their place at the top is not evidence about
  trafficking.

**Lopez-Ortega E, Choi JY, Hong I, et al. (2024).** Stimulus-dependent synaptic
plasticity underlies neuronal circuitry refinement in the mouse primary visual
cortex. *Cell Rep* 43:113966. `Lopez-Ortega_2024.pdf`

- In vivo two-photon imaging of V1 layer 2/3 over five days of repeated
  grating stimulation (spines imaged in mice of 2.5 to 3 months). Fewer
  neurons respond (44% to 32% of active neurons per field), and the
  persistent ones respond 30% more. Spine density falls by 15%. Spine GluA1,
  read as SEP-GluA1 fluorescence after in utero electroporation, rises by
  about 20% and stays up for at least 48 h; potentiated spines are clustered.
- Here: the kind of change the plasticity line looked for. Repeated sensory
  stimulation over days raises spine GluA1 in the stimulated primary cortex
  of adult mice. The paper reads SEP-GluA1 fluorescence in living tissue; in
  our fixed sections the green channel is mostly autofluorescence, so that
  reading is not available here.

### The plasticity experiment

**Gambino F, Pagès S, Kehayas V, et al. (2014).** Sensory-evoked LTP driven by
dendritic plateau potentials in vivo. *Nature* 515:116-119. `Gambino_2014.pdf`

- Whole-cell recordings in layer 2/3 of the mouse barrel cortex, under
  urethane. Rhythmic whisker stimulation (RWS: the principal whisker
  deflected at 8 Hz for 1 min) potentiates the short-latency whisker-evoked
  PSP with no somatic spike (n = 11, P = 0.008). The potentiation needs
  NMDA-receptor plateau potentials, and these depend on the posteromedial
  thalamic nucleus (POm).
- Here: the protocol behind the RWS group (one session, then perfusion), and
  the reason S1 was where its effect was expected. The regional surprise
  bars of the plasticity comparison label six regions chosen in advance as
  expected to change, the same list for RWS and for behaviour (the barrel
  field, the supplemental somatosensory area, VPM, the posterior thalamic
  complex, the zona incerta and the rostrolateral visual area); the
  highlighting is not a result. The paper measures potentiation of single
  cells over minutes in anaesthetised mice; it does not measure surface
  GluA1.

### Critical periods

**Levelt CN, Hübener M (2012).** Critical-period plasticity in the visual
cortex. *Annu Rev Neurosci* 35:309-330. `Levelt_2012.pdf`

- Review of ocular-dominance plasticity, mostly in rodent V1. In mice it is
  strongest at the end of the fourth postnatal week; P28 to P32 is generally
  taken as the critical period. The maturation of inhibitory circuits opens
  it; what ends it is less clear. For the closing, it reviews the receptors
  of myelin-associated growth inhibitors (Nogo-66 receptor, PirB), the
  extracellular matrix and perineuronal nets, the epigenetic regulation of
  CREB-mediated transcription, and neuromodulatory inputs. Plasticity exists
  well before and long after the peak, and closure is often gradual.
- Here: the reference timing for V1. The young brains analysed so far are
  P16 to P22, before that peak, and V1 shows no young-adult difference; the
  P28 brains sit at the peak. The review does not discuss AMPA receptors. The
  grant treats surface GluA1 as a proxy for plasticity potential, not a
  definitive marker of critical periods.

**Larsen B, Sydnor VJ, Keller AS, Yeo BTT, Satterthwaite TD (2023).** A critical
period plasticity framework for the sensorimotor-association axis of cortical
neurodevelopment. *Trends Neurosci* 46(10):847-862. `Larsen_2023.pdf`

- Review, mostly of human imaging. Development proceeds along a
  sensorimotor-to-association axis, and association cortex matures last. In
  animal models, critical periods progress along cortical hierarchies:
  opened by maturing parvalbumin interneurons and a falling
  excitation/inhibition ratio, closed by myelin and perineuronal nets.
  Parvalbumin-cell development starts in primary cortex and cascades to
  higher-order sensory areas (shown in primates). Rodent prefrontal cortex is
  affected by deprivation, enrichment and stress at later stages, which
  suggests a later or prolonged critical period.
- Here: the frame of the grant's prediction that associative visuo-tactile
  areas have a delayed or prolonged window. Our frontal areas are lower in
  young brains (`cref` -0.43); with two ages this cannot be read as timing.

### Comparing a brain map with gene expression

**Lein ES, Hawrylycz MJ, Ao N, et al. (2007).** Genome-wide atlas of gene
expression in the adult mouse brain. *Nature* 445:168-176. `Lein_2007.pdf`

- The Allen Mouse Brain Atlas: in situ hybridisation for about 20,000 genes in
  56-day-old male C57BL/6J mice. Each gene was run on sagittal sections 200 µm
  apart, with a coronal replicate for about 3,500 genes. The images are
  registered to a reference atlas and summarised by structure or on a
  100 to 300 µm grid. About 80% of the genes show some expression in the
  brain.
- Here: the source of the gene maps, used as expression energy on the 200 µm
  grid, downloaded per experiment. Each experiment is one mouse, and it
  measures mRNA, not protein. Hence the reliability of each map was measured
  (median 0.69 over 222 genes with repeats), and a leftover cannot be read as
  a measurement of the surface fraction. The ISH mice are younger than our
  adults (P56, against well above P60).

**Fulcher BD, Arnatkeviciute A, Fornito A (2021).** Overcoming false-positive
gene-category enrichment in the analysis of spatially resolved transcriptomic
brain atlas data. *Nat Commun* 12:2669. `Fulcher_2021.pdf`

- Asking which gene categories correlate with a brain map, by shuffling
  genes between categories, is badly anticonservative: genes in a category
  are co-expressed, and brain maps are spatially autocorrelated. With random
  phenotype maps and real mouse expression data, the mean category
  false-positive rate rose 875-fold. The fix is to compare against ensembles
  of surrogate maps, including spatially autocorrelated ones (for the mouse
  brain, rho = 0.8 and d0 = 1.46 mm).
- Here: why the ISH line tests every correlation, gene set and contrast
  against surrogate maps that keep the map's smoothness (A7), and why April's
  category ANOVA, which treated co-expressed genes as independent, does not
  stand (today's F across April's groups, against the F of the surrogates,
  gives p 0.289, where an ANOVA treating the genes as independent gives 0.011).

**Burt JB, Helmer M, Shinn M, Anticevic A, Murray JD (2020).** Generative
modeling of brain maps with spatial autocorrelation. *NeuroImage* 220:117038.

- Surrogate maps that keep a brain map's spatial autocorrelation: the map's
  values are shuffled, smoothed over each region's nearest neighbours, and the
  amount of smoothing chosen so that the surrogate's variogram matches the
  map's; their values are then replaced by the map's own, in the new order. A
  correlation is tested against the same correlation with the surrogates. The
  method of the brainsmash package.
- Here: the spatial null of the ISH line (A7), written in numpy and scipy from
  the paper (brainsmash needs packages the analysis environment does not
  have): 10,000 surrogates of the nano map on the centroids of its 204
  structures, in one hemisphere. Checked on random maps of another kind: 4.0%
  false positives at p < 0.05, against 62% for the ordinary Spearman p.

**Koopmans F, van Nierop P, Andres-Alonso M, et al. (2019).** SynGO: an
evidence-based, expert-curated knowledge base for the synapse. *Neuron*
103:217-234. `Koopmans_2019.pdf`

- 1,112 genes and 2,922 annotations to 87 synaptic locations and 179 synaptic
  processes, each based on published experimental evidence. Synaptic genes
  are compared against brain-expressed control sets.
- Here: not read directly (its release ships as spreadsheets, and the
  analysis environment has no reader for them). Its curated annotations enter
  Gene Ontology, so the GO gene sets of the ISH line, fixed before looking,
  carry those that reached the GO records.

### Atlases for young brains

**Carey H, Kleven H, Øvsthus M, et al. (2025).** DeMBA: a developmental atlas
for navigating the mouse brain in space and time. *Nat Commun* 16:8108.
`Carey_2025.pdf`

- A 4D atlas for every postnatal day from P4 to P56, interpolated between six
  templates (P4, P7, P14, P21, P28 and P56, the last being the Allen CCFv3)
  registered with elastix. It carries the CCFv3 and DevCCF labels at every
  age, and CCF Translator moves coordinates and volumes between ages. On 54
  landmarks its transformations are as accurate as an average expert. The
  templates come from different sources, so the apparent shrinking of the
  brain between P28 and P56 is likely methodological.
- Here: the atlas of the young brains, one per age, and the deformation that
  carries them to the CCF for the maps. Because the labels are the CCF's,
  region statistics need no warping. Absolute size is not compared across
  P28 and P56.

**Kronman FN, Liwang JK, Betty R, et al. (2024).** Developmental mouse brain
common coordinate framework. *Nat Commun* 15:9072. `Kronman_2024.pdf`

- DevCCF: templates at E11.5, E13.5, E15.5, E18.5, P4, P14 and P56, from MRI
  and light-sheet data, with developmental labels following the prosomeric
  model, and the CCFv3 registered into its P56.
- Here: it has no template between P14 and P56, so on its own it offers no
  age-matched template for the P20 to P36 brains. DeMBA's DevCCF labels come
  from it.

**Chon U, Vanselow DJ, Cheng KC, Kim Y (2019).** Enhanced and unified
anatomical labeling for a common mouse brain atlas. *Nat Commun* 10:5067.
`Chon_2019.pdf`

- Franklin-Paxinos labels brought into the CCF, with boundaries adjusted using
  an MRI atlas and cell-type-specific transgenic mice, and the dorsal
  striatum segmented by cortico-striatal connectivity.
- Here: not used in the code. It is the mapping to use when results need
  Franklin-Paxinos names.

## Where the code fits

| question | code | outputs, under `data\` |
|---|---|---|
| raw sections to registered volumes | `preprocessing/run_copy_raw_data.m` to `run_annotate_artifacts.m`, `registration/run_register_to_atlas.m`, `run_add_sep_channel.m`; atlases in `atlas/` | `<group>\<mouse>\lightsuite\volume_registered\`, `volume_registered_sep\` |
| where an experience changes the map (line 1) | `group_comparison/run_collect_by_group.m`, `run_normalise_groups.m`, `run_group_differences.m` | `comparisons\<ctrl>_vs_<exp>_<channel>\` |
| per-brain volumes and the readings (lines 2 to 4) | `mapping/run_per_mouse.py`, `run_to_ccf.py`, `run_cohort.py` (`sepmap/volumes/`) | `comparisons_v2\per_mouse\`, `per_mouse_ccf\`, `ccf\<cohort>\` |
| the adult distribution (line 2) | `mapping/run_region_plot.py` (the per-mouse region table); `adult_matlab/run_characterize_distribution.m` and `run_compare_nano_with_autofluorescence.m` until A4 and A5 replace them | `comparisons_v2\young_vs_adult\region_means_per_mouse.csv` |
| the inputs of the ISH line (line 3) | `mapping/run_panel_build.py`, `run_panel_fetch.py` (the ontology panel), `run_structure_set.py`, `run_ish_section_qc.py`, `run_ish_gene_table.py`, `run_ish_spatial_null.py` (`sepmap/structures.py`, `sepmap/adult/profiles.py`, `sepmap/ish/`) | `adult_v2\panel\`, `adult_v2\ish_analysis\tables\` |
| more than Gria1 expression or synapse density (line 3) | `mapping/run_beyond_density.py`, `run_beyond_controls.py`, `run_beyond_calibration.py`, `run_beyond_regression.py`, `run_beyond_figures.py` (`sepmap/adult/`) | `adult_v2\ish_analysis\beyond\` |
| the genes against the map (line 3) | `mapping/run_ish_gene_ranking.py`, `run_ish_robustness.py`, `run_ish_divisions.py`, `run_ish_gene_sets.py`, `run_ish_overview.py` (`sepmap/ish/`); [ISH_ANALYSIS.md](ISH_ANALYSIS.md) | `adult_v2\ish_analysis\` |
| what the green channel reports (line 3) | `mapping/run_sep_channel_check.py` | `adult_v2\ish_analysis\green_channel\` |
| young against adult (line 4) | `mapping/run_compare.py`, `run_region_plot.py`, `run_region_groups.py`, `run_video.py`, `run_video_compare.py`, `run_closeup.py` (flatmaps, in `tools\venv_flat`), `run_replot.py` (`sepmap/young_vs_adult/`) | `comparisons_v2\young_vs_adult\` |
| checking each step by eye | `mapping/run_diagnostics.py` | `comparisons_v2\processing_diagnostics\` |

The run order of the Python route is in the header of every `mapping/run_*.py`,
its parameters in `mapping/settings.toml`.

## Terms used in the code

- **P-numbers.** In an age, a cohort tag or an atlas key, `P<n>` is
  postnatal day n (P20, `young_P20`, `demba_p20`). The old script names P0 to
  P10 (P4, P6bis, P7bis) were pipeline steps. The refactor renamed every
  driver `run_...`, P8, P9 and P10 included, which are in `adult_matlab/`
  until A1 to A5 replace them. The table of old and new names is
  [refactor_name_map.csv](refactor_name_map.csv).
- **Channels.** At acquisition the files carry dye names (`chan02_Cy5` for
  nano, `chan03_Cy3` for auto); after registration, role names
  (`chan01_DAPI`, `chan02_NANO`, `chan03_AUTO`, `chan04_DIFF`,
  `chan05_MASK`).
  - **nano**: the nanobody label of SEP-GluA1, at the cell surface in
    principle (see the method). The registered NANO is the slice-equalised nano of
    `run_nano_equalisation`.
  - **auto**: autofluorescence. The registered AUTO is fitted onto the nano
    slice by slice (slope and intercept on reference pixels,
    `run_residual_correction`).
  - **DIFF**: nano minus the scaled autofluorescence, negatives set to zero.
    It is registered, but no analysis here reads it.
  - **MASK**: the invalid pixels, annotated artefacts or background.
  - **DAPI**: nuclei; the channel that drives registration
    (`regchan = 'dapi'`).
  - **SEP**: the green channel (`chan04_EGFP`, named after the filter),
    carried into registered space by `run_add_sep_channel` as
    `volume_registered_sep\chan02_SEP.tiff`. Mostly autofluorescence in this
    tissue (line 3).
  - `channel = 'nano'` or `'auto'` in the MATLAB drivers picks the channel a
    run analyses; it goes into the output folder names.
- **Readings** (Python route, defined in `sepmap/volumes/cohort.py`), all of
  the background-subtracted signal:
  - `ratio`: nano per unit autofluorescence, voxel by voxel. Not absolute:
    autofluorescence rises with age, so it understates a young deficit.
  - `sepratio`: nano per unit SEP. Built as a surface fraction; it is not one.
  - `cref`: nano over that brain's isocortex mean. A pure scale, so region
    ratios within a brain survive exactly; blind to a change of the whole
    cortex.
  - `subref`: nano over that brain's subcortical mean: without the
    isocortex, olfactory areas, cortical subplate, hippocampus, striatum and
    cerebellum, and from fix 1 of step 8 also without the fibre tracts, the
    ventricles and the unassigned labels (`NOT_SUBCORTEX` in `cohort.py`),
    which leaves the thalamus, hypothalamus, pallidum, midbrain, pons and
    medulla.
  - `zref`: below.
  - `sepauto` (SEP per unit autofluorescence) appeared only in the channel
    arms of line 3, now in `archive/`.
- **zref.** zref = [log2(nano / cortex mean) - m] / s, with m the median and
  s the p90 - p10 spread of log2(structure mean / cortex mean) across that
  brain's structures. Zero is the brain's median structure, mostly
  subcortical, not a physical null. One unit is that brain's own spread. A
  zref difference is already a difference of logs and cannot be converted to
  a fold change: with the spreads of line 4, +1 is about 2x in a young brain
  and about 3.8x in an adult. Quote numbers from the region tables and use the
  maps for the pattern: the maps are smoothed voxel means, and the maps and
  the young-against-adult tables still take m and s over different
  structures. The ISH line takes them over its declared set (A1), the 204
  grey-matter structures measured in all ten adults.
- **LR-sum and LR-diff** (MATLAB route, `common/compute_lr_stats.m`): left
  hemisphere plus, or minus, the mirrored right hemisphere. The plasticity
  comparison uses L + R and |L - R| per mouse.
- **Surprise.** -log10 p of the group-difference test, drawn as a map.
- **Cohorts.** The MATLAB groups and data folders are `naive`, `rws`,
  `behavior` and `young` (registry: `common/get_cohort.m`, with each mouse's
  age). The Python cohorts are `young` (P16, P20 and P22 pooled),
  `young_P20`, `young_P16`, `young_P22`, `naive`, `rws` and `adult` (naive
  plus RWS). Adult folders end in `_Gria1`, young ones in `_SepGluA_P<age>`,
  the mouse's age. The adults were well above P60; the registry holds NaN for
  their age, which is not recorded yet.
- **CCF.** The Allen Mouse Brain Common Coordinate Framework v3, the P56
  reference, a population average of 1,675 mice. Adults are registered at
  10 µm, in the crop of CCF planes 180 to 1079. The Python route works on the
  CCF at 20 µm (660 x 400 x 570).
- **DeMBA.** The developmental atlas of Carey et al. 2025.
  `atlas/build_demba_atlas.py <age>` builds one per age from BrainGlobe's
  `demba_allen_seg_dev_mouse_p<age>_20um`, with the labels remapped to the
  CCF's parcellation index. The files keep LightSuite's `_10` names, so a
  young brain's `local_settings.txt` says `px_atlas = 20` for the alignment;
  its registered volumes are on the 10 µm grid all the same
  (`registered_grid_um`). CCF Translator
  (`brainglobe-ccf-translator`) carries a young brain to the CCF.
- **RWS.** Rhythmic whisker stimulation, as in Gambino et al. 2014: one
  session, then perfusion.
- **Structure sets.** The declared set of the ISH line (S1, A1): grey-matter
  structures measured in all ten adults, with no "..., unassigned" labels, 204
  of the 280 of the adult table
  (`adult_v2\ish_analysis\tables\structure_set.csv`, each dropped structure
  with its reason). The beyond-abundance analysis uses those of them that the
  ISH maps of the subunit and marker genes cover (126,
  `beyond\structures_used.csv` there).
- **Spatial null, surrogate.** A surrogate is a random map with the nano
  map's smoothness (its variogram); a spatial p counts how many of 10,000
  surrogates correlate with a gene at least as strongly as the map does. For a
  residual (the leftover, a partial correlation), each surrogate goes through
  the same fit first (Freedman-Lane).
- **Calibration.** In the beyond-abundance analysis, the same model run on
  maps whose answer is known; its floor is what Allen-to-Allen mismatch alone
  leaves, compared with the nano map on the same structures.
- **Jackknife.** Recomputing a number on subsamples that each leave out a
  fifth of the structures; the spread gives the number's interval over
  structures.
- **Reliability and Spearman-Brown.** Agreement between two independent
  halves of a cohort; Spearman-Brown extends it to the full cohort.
- **Cross-validated R².** Variance explained on structures held out of the
  fit.
- **Expression energy.** The Allen ISH value per voxel of its 200 µm grid.
- **A1 to A10, S1 to S6.** The scientific additions and decisions of
  [history/REFACTOR_PLAN.md](history/REFACTOR_PLAN.md).
