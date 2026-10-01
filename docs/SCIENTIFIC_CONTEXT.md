# Scientific context

What this code is for: the questions, what each line of work found, the grant
it now serves, and the papers it builds on. The papers are in
`data\ref_papers\`, with summaries in `PAPER_SUMMARIES.md` there.

The results are those recorded up to 1 October 2026. Some of them will move in
step 9 of the refactor, when the analyses are rerun on one declared set of
structures (S1 to S5 and A1 to A5 in [REFACTOR_PLAN.md](REFACTOR_PLAN.md)).
Those are marked provisional.

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

The project maps the nanobody signal, in principle surface GluA1, across the
whole mouse brain and asks four questions of the map:

1. Where does an experience change it?
2. How is it distributed in the adult brain, and how reproducible is that
   across mice?
3. What does it measure: more than receptor abundance or synaptic density?
4. How does it differ between young and adult mice, and could the
   differences mark critical periods?

## The starting point: the method

- **Mice.** SEP-GluA1 knock-ins: the GluA1 subunit is fused to SEP, a
  pH-sensitive fluorophore that reports surface-exposed receptors (the grant,
  citing Graves et al. 2021, eLife).
- **Labelling.** As the grant describes it: tissue processed without
  detergent, to keep the label on the membrane and away from intracellular
  receptor pools, and the GFP part of SEP-GluA1 amplified with a GFP-booster
  nanobody. In preliminary preparations this showed surface GluA1 with no
  detectable intracellular labelling (grant, Fig. 3a). The project's own
  records do not agree on this: the header of `run_add_sep_channel.m` calls
  the tissue fixed and permeabilised, and an April note reads nano in fixed
  tissue as total SEP-GluA1, not strictly surface. Which holds for these
  brains is to be confirmed from the staining protocol; until then the map is
  called the nanobody signal, in principle surface GluA1.
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
- **Method.** For every mouse the hemispheres are folded (left plus mirrored
  right, and the absolute left-right difference). The experimental group's
  intensity profile is aligned onto the control group's. The groups are then
  compared voxel by voxel, with Welch t and surprise (-log10 p) maps.
- **Result.** The headline is an increase of the nanobody signal in S1 after
  RWS, in the hemisphere-sum t map of the slab around plane 565, masked at
  p < 0.01. That threshold is uncorrected and applied voxel by voxel, on maps
  smoothed in 3D (Gaussian, sigma 5 voxels), with five mice per group; the
  figure shows the median over planes 555 to 575. The increase is small. It is
  where plasticity caused by whisker stimulation was expected.
- **Reproduced.** On 1 October 2026 today's group-difference step was rerun
  on the normalised volumes of 26 and 27 November 2025, the inputs of the
  approved run. The slab t maps and surprise masks correlate with the
  approved ones at 0.997 to 0.998 for RWS (0.988 to 0.999 for behaviour), the
  regional bars at 0.994 to 0.998, and the S1 increase is there. The small
  residue comes from the background masks, which were regenerated since.
  This checks the computation on the same mice; it is not a replication in
  new animals.
- **Status.** The team was not fully confident that the effect is robust, and
  more animals would be a large investment, so the line is paused, not
  closed. The code stays runnable, with identical results, and documented
  well enough to resume with more animals.

### 2. The adult distribution

- **Question.** How the nanobody signal is distributed across the whole
  adult brain, and how reproducible that is across mice.
- **Cohort.** The five naive and five RWS mice, pooled. Pooling was checked:
  the part of the map that abundance and density do not explain (line 3)
  agrees between the naive and the RWS mice at rho +0.88. The behaviour mice
  enter only the plasticity comparison.
- **Reproducibility.** Over the 126 ways of splitting the ten adults into two
  halves of five, the two half-cohort maps agree at rho 0.974 over 125
  grey-matter structures (Spearman-Brown 0.987 for the full cohort). The
  autofluorescence of the same brains explains none of the map
  (cross-validated R² -0.043).
- **Being rebuilt.** Which structures stand out, and with what confidence, is
  being redone. P8 called a structure enriched above a threshold on its own
  scale, and Sami El-Boustani asked (28 April) for a stricter threshold and a
  permutation null. The plan replaces both with a test per structure against
  the brain's median structure (S2), on a declared set of structures (S1,
  A1, A4). Today 22 structures of the adult table are seen in fewer than five
  of the ten adults. P8's outputs carry known defects (listed in the plan) and
  are not quoted here.

### 3. What the map measures

A reader will ask whether the nanobody map simply follows how much receptor a
region makes, or how many synapses it has. The project tests this against the
Allen in situ hybridisation (ISH) maps of the adult mouse brain (Lein et al.
2007).

**The map is more than abundance and density.** Over 125 grey-matter
structures, the adult map was predicted from receptor abundance (Gria1 to
Gria4 mRNA), synaptic markers, the first principal component of 188
postsynaptic-density genes, and the cohort's own autofluorescence. Each model
was scored by cross-validation on held-out structures, against the map's own
reliability (the ceiling, 97.4% of the variance):

| predictors | cross-validated R² | share of the explainable variance |
|---|---|---|
| receptor abundance (Gria1 to Gria4) | 0.253 | 26% |
| synaptic markers | 0.148 | 15% |
| postsynaptic density, first component | 0.262 | 27% |
| autofluorescence | -0.043 | 0% |
| all four, straight lines | 0.416 | 43% |
| all four, allowed to bend | 0.597 | 61% |

- About 39% of the explainable variance is left over. The leftover
  replicates across independent halves of the cohort at rho 0.934
  (Spearman-Brown 0.966), almost as well as the map itself.
- Seven controls tried to break this: a spatial gradient, structure size,
  single animals, naive against RWS, curvature, the whole gene space and the
  choice of reading. For single animals, the leftovers of every pair of mice
  agree, at a median rho of 0.780 (worst 0.595). The curvature control failed
  the first time: fitted with straight lines, a fifth of the map was credited
  to the leftover. The model was changed to allow curvature, and the claim
  fell from about half to 39%.
- The richest model tried, 20 components of the whole 390-gene panel, reaches
  a cross-validated R² of 0.826 (85% of the ceiling). Its leftover still
  replicates at 0.879.
- The leftover is high, relative to what the predictors give, in the medial
  geniculate (+48 ranks), subthalamic nucleus (+44), ventral lateral
  geniculate (+38) and lateral habenula (+38). It is low in VPM (-49), VPL
  (-47), dorsal retrosplenial cortex (-46) and the posterior thalamic complex
  (-43). No single gene accounts for it (best, Cacng8, rho +0.301).
- This result uses all ten adults and the grey-matter rule, so the
  structure-set changes of step 9 do not affect it.
- **What it means.** It is consistent with the stain reporting surface
  GluA1, and it does not show it. A residual is only what the predictors did
  not explain: regional differences in translation, turnover, subunit
  composition or nanobody access to the tissue would also land there.

**The ISH gene ranking is descriptive.**

- The first comparison (P9, April 2026) correlated the map with 100
  hand-picked genes and found AMPA receptor trafficking and anchoring genes,
  Cacng8 first, above Gria1 itself.
- The Python route reproduces that ranking (rho 0.91 against the old
  ordering). Cacng8 (TARP γ-8) is first under every reading. Under `zref` it
  stays first in every structure set tried, at rho 0.77 to 0.82.
- Gria1's own rank moves with the structure set. It is not quoted until A1 to
  A3 have run (S5 in the plan).
- A single Allen ISH experiment is more reliable than assumed. For the 218
  genes measured more than once, the median agreement between experiments is
  rho 0.69, and Gria1's map scores 0.91. Gria1 ranking below Cacng8 is not a
  bad Gria1 experiment.
- Cacng8 above Gria1 is a single-gene result. On the 100-gene panel, the
  subunit genes as a group correlate better with the map than the
  localisation genes. Splitting the same 19 genes by function does no better
  than splitting them at random (an exact permutation over every split).
- **The powered test is negative.** A 390-gene panel was built from Gene
  Ontology terms, not by hand. Once the subunit composite is removed, the 84
  AMPA receptor localisation genes explain no more of the map than
  expression-matched postsynaptic genes. As a positive control, the same test
  detects the difference between control genes of high and low map
  reliability, so the negative is informative. Genes such as Arpc5 and
  Cdk5r1, which have nothing to do with AMPA receptor traffic, sit in the same
  band as Cacng8. The map is predicted about equally well by any
  well-measured forebrain postsynaptic gene.
- **Provisional.** These ISH results were computed on the unrestricted
  structure set; A1 to A3 rerun them. Under S5 the numbers of the powered
  test and of the permutation are held back until then; the values computed
  so far are in `docs/adult_ish_design.md`. Their p values are
  anticonservative, because genes within a set are co-expressed (Fulcher et
  al. 2021); a spatial null is planned (A7).

**The green channel is not total receptor.**

- **The plan.** The SEP tag fluoresces green, so the green channel was meant
  to report all SEP-GluA1, surface and internal, and nano divided by it a
  surface fraction. Three predictions were set in advance: SEP should track
  Gria1, nano both, and nano/SEP the trafficking genes.
- **Measured.** In these fixed, cleared sections the green channel is mostly
  autofluorescence. In all ten adults it tracks the autofluorescence channel
  at rho 0.79 ± 0.04 across structures. It varies less across the brain than
  either autofluorescence or nano: p90 - p10 of log2, 0.95 ± 0.16 against
  1.07 ± 0.16 and 1.93 ± 0.26. Against Gria1 ISH, nano correlates at +0.60 and
  the green channel at +0.30.
- **Consequences.** `sepratio` (nano per unit SEP) is not a surface fraction,
  and the three-way test cannot be run on these data.
- **What remains.** With Gria1 mRNA partialled out of the nano map (`ratio`
  reading), anchoring and trafficking genes still predict the remainder
  (Cacng8 +0.54, Cnih2 +0.50; median +0.12 over 33 genes, 20 of them
  positive). This is weak, and not decisive. The median is small, and the
  powered test above shows that other well-measured postsynaptic genes predict
  the map about as well, so it is not specific to anchoring or trafficking.
  And mRNA is not protein: a regional gradient of translation or turnover
  would look the same.
- **What would settle it.** A wet-lab control: a total-GluA1 antibody stain
  or autoradiography on a subset of the same brains. A knockout or
  no-primary control would test specificity instead. There is no knockout or
  competition control, no independent regional measure and no DAPI control
  yet.

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
  (-0.4 in retrosplenial to -1.2 log2 in frontal cortex), but that is a
  bound, not a value.
- **The young brain is flatter.** The spread of its structures (p90 - p10 of
  log2 nano relative to cortex) is about half the adult one: 0.9 to 1.0 log2
  in young brains against 1.8 to 1.9 in adults, depending on which brains
  and which structures enter. `zref` removes this per brain, so a zref
  difference is a difference of positions within each brain's own range, not
  a fold change (see Terms).
- **By system** (7 young against 10 adults, `cref`, log2 young minus adult,
  group medians from `group_stats.csv`; Mann-Whitney p < 0.01, uncorrected,
  unless marked): primary somatosensory +0.21, retrosplenial +0.24, frontal
  -0.42, striatum -0.79, hippocampus -1.03; primary visual -0.06, not
  significant. The reading matters: hippocampus is -1.03 under `cref` and
  +0.49 under `zref`.
- **RL and AL.** RL and AL together, the visuo-tactile areas of the grant,
  are the only visual group that survives correction under `zref`: +0.23
  (p = 0.0007, q = 0.014), against -0.01 between naive and RWS adults; the
  P20 brains alone give +0.34. This holds under `zref` only: under `cref`,
  which is blind to a shift of the whole cortex, RL+AL is +0.08 (p = 0.36).
  V1 shows no difference (`zref` +0.07, q = 0.35).
- **By layer.** Under `cref`, in every sensory system, primary and higher
  order alike, the young-minus-adult difference is larger in the
  supragranular layers (+0.11 to +0.30) than in the infragranular ones. The
  infragranular layers are lower in young brains in the visual and auditory
  systems (-0.13 to -0.25), but slightly higher in somatosensory cortex
  (SSp +0.07, SSs +0.11). Under `zref`, RL+AL is higher in both the
  supragranular layers (+0.28) and layer 4 (+0.27, q = 0.006 each), while
  layer 4 does not differ in V1 or in the other visual areas (+0.13 and
  +0.15, q ≥ 0.11): RL and AL have a laminar profile closer to somatosensory
  cortex than to V1.
- **The laminar contrast.** Taken per mouse, the supragranular-minus-
  infragranular contrast is larger in young brains in every laminar system
  under `cref` (+0.05 in frontal to +0.51 in RL+AL; p < 0.05 in seven of the
  nine systems, uncorrected). It cancels any scale factor of a brain
  (exposure, staining strength, the choice of reference), so `subref` gives
  the same values and `zref` the same signs. Under `ratio` and `sepratio`,
  whose denominators vary between layers, it is not positive everywhere: it
  is negative in frontal cortex and in the higher visual areas under both,
  and about zero in retrosplenial cortex under `ratio`. These values were
  recomputed for this document; no output of the code holds them yet.
- **A critical period cannot be claimed yet.** Two ages cannot show a window.
  V1, where the classic mouse critical period sits (Levelt & Hübener 2012),
  shows no difference at P16 to P22, while among the sensory systems
  somatosensory cortex shows the largest. And the superficial-layer excess is
  also what the later maturation of the superficial layers alone would give.
  The P28 to P36 brains can separate the two: if V1's superficial excess
  peaks near P28 and falls by P36 while S1 has already declined, that is a
  window moving across modalities; the same decline everywhere at the same
  rate is maturation.
- **Registration caveats.** Registration is driven by DAPI, whose cell
  packing falls between P20 and adulthood. In the P20 atlas, 19% of the volume
  that is fibre tract in the adult atlas is labelled isocortex, so results
  next to white matter need care.
- **Grant figure.** The grant's Figure 3b to d shows this comparison: maps of
  pups (n = 7) and adults (n = 10), a flatmap of the pup-minus-adult
  difference, and Mann-Whitney tests on V1, RL+AL, the other higher visual
  areas, S1, S2 and medial prefrontal cortex. The figure set chosen for the
  grant includes the layer 5-6 flatmap. That band is the weakest of the three
  statistically (RL+AL infragranular q = 0.110, against 0.006 for the other
  two), and it looks strong partly because the adult baseline is high there.
- **Provisional.** A1 changes the reference of zref and so moves every zref
  value, the differences behind the grant figures included. The code state
  behind the grant figures is tagged `grant-2026-09`.

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
- Here: the protocol behind the RWS group, and the reason S1 was where its
  effect was expected. The regional surprise bars of the plasticity
  comparison highlight the same six areas for RWS and for behaviour: the
  barrel field, the supplemental somatosensory area, VPM, the posterior
  thalamic complex, the zona incerta and the rostrolateral visual area. The
  paper measures potentiation of single cells over minutes in anaesthetised
  mice; it does not measure surface GluA1.

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
  young brains (`cref` -0.42); with two ages this cannot be read as timing.

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
  (median 0.69 over 218 genes with repeats), and a residual cannot be called
  "surface".

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
- Here: why the ISH ranking is descriptive, why its p and q values are
  anticonservative, and the design of the planned spatial null (A7).

**Koopmans F, van Nierop P, Andres-Alonso M, et al. (2019).** SynGO: an
evidence-based, expert-curated knowledge base for the synapse. *Neuron*
103:217-234. `Koopmans_2019.pdf`

- 1,112 genes and 2,922 annotations to 87 synaptic locations and 179 synaptic
  processes, each based on published experimental evidence. Synaptic genes
  are compared against brain-expressed control sets.
- Here: not used yet. The 390-gene panel comes from Gene Ontology terms, each
  gene with the terms that placed it. SynGO is the curated alternative for
  fixing synaptic gene sets before looking.

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
| the adult distribution (line 2) | `mapping/run_region_plot.py` (the per-mouse region table); `P8_characterize_merged_distribution.m` and `P10_compare_nano_vs_auto.m` until A4 and A5 replace them | `comparisons_v2\young_vs_adult\region_means_per_mouse.csv` |
| more than abundance or density (line 3) | `mapping/run_beyond_density.py`, `run_beyond_controls.py`, `run_beyond_figures.py`, `run_beyond_regression.py` (`sepmap/adult/`) | `adult_v2\beyond\` |
| what the green channel reports (line 3) | `mapping/run_sep_channel_check.py`, `run_adult_arms.py`, `run_ish_arms.py` | `adult_v2\arms\` |
| the ISH gene comparison (line 3) | 100-gene panel: `mapping/run_ish_regions.py --panel targets`, `run_ish_compare.py`, `run_ish_words.py`, `run_ish_roles.py`; 390-gene panel: `run_panel_build.py`, `run_panel_fetch.py`, `run_ish_regions.py --panel ontology`, `run_ish_reliability.py`, `run_ish_panel_test.py` (`sepmap/ish/`); `P9_compare_nano_vs_allen_ish.m` until A1 to A3 replace it | `adult_v2\ish\`, `adult_v2\panel\` |
| young against adult (line 4) | `mapping/run_compare.py`, `run_region_plot.py`, `run_region_groups.py`, `run_video.py`, `run_video_compare.py`, `run_closeup.py` (flatmaps, in `tools\venv_flat`), `run_replot.py` (`sepmap/young_vs_adult/`) | `comparisons_v2\young_vs_adult\` |
| checking each step by eye | `mapping/run_diagnostics.py` | `comparisons_v2\processing_diagnostics\` |

The run order of the Python route is in the header of every `mapping/run_*.py`.

## Terms used in the code

- **P-numbers.** In an age, a cohort tag or an atlas key, `P<n>` is
  postnatal day n (P20, `young_P20`, `demba_p20`). The old script names P0 to
  P10 (P4, P6bis, P7bis) were pipeline steps. The refactor renames the
  drivers `run_...`; P8, P9 and P10 keep their names until A1 to A5 replace
  them. The table of old and new names goes into `docs/ADDING_DATA.md`.
- **Channels.** At acquisition the files carry dye names (`chan02_Cy5` for
  nano, `chan03_Cy3` for auto); after registration, role names
  (`chan01_DAPI`, `chan02_NANO`, `chan03_AUTO`, `chan04_DIFF`,
  `chan05_MASK`).
  - **nano**: the nanobody label of SEP-GluA1, in principle of the surface
    receptors (see the method). The registered NANO is the slice-equalised
    nano of `run_nano_equalisation`.
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
  - `subref`: nano over that brain's subcortex mean, without hippocampus and
    striatum.
  - `zref`: below.
  - `sepauto` (SEP per unit autofluorescence) appears only in the channel
    arms of line 3.
- **zref.** zref = [log2(nano / cortex mean) - m] / s, with m the median and
  s the p90 - p10 spread of log2(structure mean / cortex mean) across that
  brain's structures. Zero is the brain's median structure, mostly
  subcortical, not a physical null. One unit is that brain's own spread. A
  zref difference is already a difference of logs and cannot be converted to
  a fold change: with the spreads of line 4, +1 is about 2x in a young brain
  and 3.5x to 3.8x in an adult.
  Quote numbers from the tables and use the maps for the pattern (over 238
  areas the two agree at r = 0.985, but the maps read lower).
- **LR-sum and LR-diff** (MATLAB route, `common/compute_lr_stats.m`): left
  hemisphere plus, or minus, the mirrored right hemisphere. The plasticity
  comparison uses L + R and |L - R| per mouse.
- **Surprise.** -log10 p of the group-difference test, drawn as a map.
- **Cohorts.** The MATLAB groups and data folders are `naive`, `rws`,
  `behavior` and `young` (registry: `common/get_cohort.m`). The Python
  cohorts are `young` (P16, P20 and P22 pooled), `young_P20`, `young_P16`,
  `young_P22`, `naive`, `rws` and `adult` (naive plus RWS). Adult folders end
  in `_Gria1`, young ones in `_SepGluA_P<age>`. The adults' age is not
  recorded in the registry.
- **CCF.** The Allen Mouse Brain Common Coordinate Framework v3, the P56
  reference, a population average of 1,675 mice. Adults are registered at
  10 µm, in the crop of CCF planes 180 to 1079. The Python route works on the
  CCF at 20 µm (660 x 400 x 570).
- **DeMBA.** The developmental atlas of Carey et al. 2025.
  `atlas/build_demba_atlas.py <age>` builds one per age from BrainGlobe's
  `demba_allen_seg_dev_mouse_p<age>_20um`, with the labels remapped to the
  CCF's parcellation index. The files keep LightSuite's `_10` names, so a
  young brain needs `px_atlas = 20`. CCF Translator
  (`brainglobe-ccf-translator`) carries a young brain to the CCF.
- **RWS.** Rhythmic whisker stimulation, as in Gambino et al. 2014.
- **Structure sets.** The beyond-abundance analysis keeps grey-matter
  structures seen in all ten adults, with no "..., unassigned" labels, that
  the ISH maps of the subunit and marker genes cover and that have
  autofluorescence values in every adult (125 structures). The plan's
  declared set (S1) keeps grey-matter structures seen in all ten adults, 234
  of 280.
- **Reliability and Spearman-Brown.** Agreement between two independent
  halves of a cohort; Spearman-Brown extends it to the full cohort.
- **Cross-validated R².** Variance explained on structures held out of the
  fit.
- **Expression energy.** The Allen ISH value per voxel of its 200 µm grid.
- **A1 to A10, S1 to S6.** The scientific additions and decisions of
  [REFACTOR_PLAN.md](REFACTOR_PLAN.md).
