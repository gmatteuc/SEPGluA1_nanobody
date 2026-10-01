# Sources checklist (temporary)

For the reviewer and the owner, to be deleted once checked. Every factual
sentence of [SCIENTIFIC_CONTEXT.md](SCIENTIFIC_CONTEXT.md), in order, with
its source and a short quotation (the anchor) found at that place.

How to read the sources:

- `ref_papers/<file>.pdf p. N`: `data\ref_papers\<file>.pdf`, page N of
  the PDF file (not the journal page number).
- `ref_papers/PAPER_SUMMARIES.md:L`: the summaries in the same folder, line L.
- `notes/<file>:L`: the project's working notes, the `project_*.md` files
  kept outside this repository, line L.
- `docs/...`, `mapping/...` and other repository paths: this branch, line L,
  as of the commit before this file (586430f).
- `REFACTOR_PLAN.md:L` (main): `D:\sep_histology\code\docs\REFACTOR_PLAN.md` on
  `main` at d0860ff, line L. This branch carries an older copy of the plan,
  without the Progress section and its S6 result.
- `D:/dendrites/code/...`: the two-photon imaging repository, read only.
- search: a claim that something is absent, checked by searching the code.
- editorial: a sentence about the document itself.

Anchors are quoted as the source has them, including PDF extraction
artefacts ("n 5 11" is "n = 11" in the PDF; "g-8" is "γ-8"). Every anchor
was checked automatically at the cited page or lines (within two lines).

## Scientific context

**1.** What this code is for: the questions, what each line of work found, the grant it now serves, and the papers it builds on.

- `REFACTOR_PLAN.md:731-738` (main): "docs/SCIENTIFIC_CONTEXT.md**: the question; the four lines of work and what each found"

**2.** The papers are in `data\ref_papers\`, with summaries in `PAPER_SUMMARIES.md` there.

- `REFACTOR_PLAN.md:735-737` (main): "the PDFs stay in `data\ref_papers\`, whose `PAPER_SUMMARIES.md` is the start"

**3.** The results are those recorded up to 1 October 2026.

- `REFACTOR_PLAN.md:1010` (main): "**1 Oct, step 6, the big files**"

**4.** Some of them will move in step 9 of the refactor, when the analyses are rerun on one declared set of structures (S1 to S5 and A1 to A5 in [REFACTOR_PLAN.md](REFACTOR_PLAN.md)).

- `REFACTOR_PLAN.md:441-443` (main): "These change numbers on purpose (A1 moves every zref"
- `REFACTOR_PLAN.md:705-710` (main): "A1, A2, A3, A4, A5 in that order"

**5.** Those are marked provisional.

- editorial: describes this document

## The question

**6.** AMPA receptors carry most fast excitatory transmission in the brain, and each is built from four core subunits, GluA1 to GluA4 (Lopez-Ortega et al. 2024).

- `ref_papers/Lopez-Ortega_2024.pdf` p. 2: "AMPARs are the main postsynaptic mediators of fast excitatory synaptic transmission in the mammalian brain, and their structure is formed by the combination of 4 core sub- units (GluA1–GluA4)"

**7.** Regulating their function and their traffic to and from the membrane is critical for many forms of synaptic plasticity, and GluA1-containing receptors are recruited to spines after long-term potentiation (Huganir & Nicoll 2013).

- `ref_papers/Huganir_2013.pdf` p. 1: "the modulation of the AMPA receptor function and membrane trafﬁcking is critical for many forms of synaptic plasticity"
- `ref_papers/Huganir_2013.pdf` p. 3: "GFP-GluA1 was recruited to synaptic spines after LTP induction"

**8.** Current models of LTP need a sizeable pool of receptors already on the cell surface (Huganir & Nicoll 2013).

- `ref_papers/Huganir_2013.pdf` p. 5: "The requirement for a signiﬁcant surface pool of receptors for the expression of LTP is a recent recurring theme in current models of LTP"

**9.** The grant this project serves uses surface GluA1 as a molecular readout of synaptic remodelling potential: a proxy for plasticity potential, not a definitive marker of critical periods.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "Surface expression of GluA1 therefore provides a useful molecular readout of synaptic remodeling potential"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 12: "interpret surface GluA1 as a proxy for plasticity potential, not as a definitive marker of critical periods"

**10.** The project maps surface GluA1 across the whole mouse brain and asks four questions of the map:

- `REFACTOR_PLAN.md:28-55` (main): "The project in four lines of work"
- `notes/project_scientific_story.md:11-16`: "Giulio's framing (30 Sep 2026)"

**11.** Where does an experience change it?

- `REFACTOR_PLAN.md:33-37` (main): "to find where plasticity lands in cortex or elsewhere"

**12.** How is it distributed in the adult brain, and how reproducible is that across mice?

- `REFACTOR_PLAN.md:41-44` (main): "is distributed across the whole brain, and how reproducible that is across mice"

**13.** What does it measure: more than receptor abundance or synaptic density?

- `REFACTOR_PLAN.md:49-51` (main): "the map is not explained by receptor abundance or synaptic density alone"

**14.** How does it differ between young and adult mice, and could the differences mark critical periods?

- `REFACTOR_PLAN.md:53-55` (main): "looking for differences in cortex that could correspond to critical periods"

## The starting point: the method

**15.** **Mice.** SEP-GluA1 knock-ins: the GluA1 subunit is fused to SEP, a pH-sensitive fluorophore that reports surface-exposed receptors (the grant, citing Graves et al. 2021, eLife).

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "in which the GluA1 subunit of AMPA receptors is fused to a pH-sensitive fluorophore that reports surface-exposed receptors"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 5: "Non-detergent labeling in SEP-GluA1 knock-in mice"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 21: "Graves, A. R., Roth, R. H., Tan, H. L."

**16.** **Labelling.** As the grant describes it: tissue processed without detergent, to keep the label on the membrane and away from intracellular receptor pools, and the GFP part of SEP-GluA1 amplified with a GFP-booster nanobody.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "Tissue will be processed under non-detergent conditions to preserve membrane labeling and minimize access to intracellular receptor pools. The GFP component of SEP-GluA1 will be amplified using GFP-booster nanobodies"

**17.** In preliminary preparations this showed surface GluA1 with no detectable intracellular labelling (grant, Fig. 3a).

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "reveals surface GluA1 with no detectable intracellular labeling in preliminary preparations (Fig. 3a)"

**18.** **Channels.** Each section is recorded in four channels: DAPI (nuclei), the nanobody (Cy5, called nano), autofluorescence (Cy3, called auto) and the green channel of SEP itself (filter EGFP).

- `notes/project_preprocessing_recipe.md:58`: "uses **dye** names: `chan01_DAPI.tiff`, `chan02_Cy5.tiff` (= NANO), `chan03_Cy3.tiff` (= AUTO)"
- `registration/run_add_sep_channel.m:8`: "The microscope recorded a sixth, on the green filter it names EGFP"
- `preprocessing/make_ordering_volume.m:46`: "blue_channel = 'dapi'; % nuclei"

**19.** **Registration.** Serial sections are registered to an atlas with LightSuite (elastix, affine and B-spline), driven by the DAPI channel and by control points between each section and the atlas.

- `registration/run_add_sep_channel.m:38`: "per-slice elastix B-spline + affine"
- `notes/project_preprocessing_recipe.md:55`: "Registration is driven by the DAPI channel"
- `registration/pipeline/register_to_atlas.m:346`: "have different numbers of histology and atlas points"

**20.** The points are placed by hand, or proposed by the automatic annotation and reviewed by hand.

- `notes/project_atlas_decision_ccf_for_all.md:15`: "the manual control points carry the correspondence"
- `REFACTOR_PLAN.md:8-9` (main): "the automatic annotation is merged"
- `REFACTOR_PLAN.md:862-863` (main): "the proposals are reviewed by hand anyway"

**21.** Adults are registered to the Allen CCF, young brains to the DeMBA atlas of their own age (decision of 2 September 2026).

- `notes/project_atlas_decision_ccf_for_all.md:13`: "FINAL DECISION 2026-09-02 (Sami): young cohort registers to the age-matched **DeMBA P20**; adults stay on **CCF**"
- `notes/project_atlas_decision_ccf_for_all.md:29`: "One atlas per AGE, not one per cohort"

**22.** **Two analysis routes.** The plasticity comparison stays in MATLAB (`group_comparison/`).

- `REFACTOR_PLAN.md:95-96` (main): "MATLAB stays for preprocessing, registration, and the plasticity comparison"

**23.** The adult map, the ISH comparison and young against adult are in Python (`mapping/`).

- `REFACTOR_PLAN.md:94-95` (main): "Python is the main analysis route: adult distribution, ISH, young against adult"

**24.** They normalise differently: the Python route scales each brain on its own; the plasticity comparison fits each brain to its group's median and the experimental group onto the control group.

- `REFACTOR_PLAN.md:183-186` (main): "The Python route normalises each brain on its own"
- `REFACTOR_PLAN.md:185-187` (main): "the plasticity comparison fits each brain to its group's median and the experimental group onto the control group"

**25.** If analyses of both routes go into one paper, the two must be made comparable, or the difference discussed.

- `REFACTOR_PLAN.md:187-189` (main): "must be made comparable, or the difference discussed"

## 1. Plasticity localisation (the original aim; paused)

**26.** **Question.** The project began as a search for where plasticity lands after an experience: surface GluA1 in mice after rhythmic whisker stimulation (RWS, the protocol of Gambino et al. 2014) or after behaviour, against naive mice, in cortex or elsewhere.

- `REFACTOR_PLAN.md:33-36` (main): "Surface GluA1 after an experience (rhythmic whisker stimulation, RWS, as in Gambino et al. 2014; or behavior) against naive mice, to find where plasticity lands in cortex or elsewhere"
- `notes/project_scientific_story.md:13`: "to find the locus of plasticity in cortex or elsewhere"

**27.** **Mice.** The approved comparison is the run of 5 December 2025: five naive mice (CGF027, CGF028, CGF033, CGF034, CGF035), five RWS mice (MG691, MG692, MG693, MG736, MG737) and four behaviour mice (MG705, MG709, MG716, MG718) of the seven registered.

- `REFACTOR_PLAN.md:160-161` (main): "The reference is the output of 5 Dec 2025, approved by Sami"
- `REFACTOR_PLAN.md:164-166` (main): "Its mice: naive CGF027, CGF028, CGF033, CGF034, CGF035; RWS MG691, MG692, MG693, MG736, MG737; behavior MG705, MG709, MG716, MG718 (of the 7 in the registry"

**28.** **Method.** For every mouse the hemispheres are folded (left plus mirrored right, and the absolute left-right difference).

- `group_comparison/run_group_differences.m:15`: "folds the hemispheres: for every mouse |L - R| and L + R"
- `common/compute_lr_stats.m:13`: "flipped_right = flip(right, 3);"

**29.** The experimental group's intensity profile is aligned onto the control group's.

- `group_comparison/run_group_differences.m:12-13`: "aligns the experimental group onto the control one: a line fitted between the two groups' mean profiles"

**30.** The groups are then compared voxel by voxel, with Welch t and surprise (-log10 p) maps.

- `group_comparison/run_group_differences.m:7`: "voxel by voxel"
- `group_comparison/run_group_differences.m:16-17`: "Welch t and surprise (-log10 p) maps"

**31.** **Result.** The headline is an increase of the nanobody signal in S1 after RWS, in the hemisphere-sum t map of the slab around plane 565, masked at p < 0.01.

- `REFACTOR_PLAN.md:162-164` (main): "Its headline is the increase of the nanobody signal in S1 after RWS (hemisphere-sum t map, slab 565, surprise-masked at p < 0.01)"

**32.** The increase is small.

- `notes/project_scientific_story.md:20`: "a little bump"

**33.** It is where plasticity caused by whisker stimulation was expected.

- `notes/project_scientific_story.md:20`: "as expected for plasticity localised in S1 caused by RWS"

**34.** **Reproduced.** On 1 October 2026 today's code was rerun on the preserved inputs of that run.

- `REFACTOR_PLAN.md:935-936` (main): "**1 Oct, S6, informative**: today's P7bis on the inputs of the approved December 2025 figures reproduces them"

**35.** The slab t maps and surprise masks correlate with the approved ones at 0.997 to 0.998 for RWS (0.988 to 0.999 for behaviour), the regional bars at 0.994 to 0.998, and the S1 increase is there.

- `REFACTOR_PLAN.md:936-939` (main): "correlate 0.997 to 0.998 for RWS, 0.988 to 0.999 for behavior; individual maps 1.000000; regional bars 0.994 to 0.998"
- `REFACTOR_PLAN.md:939` (main): "The S1 increase is there"

**36.** This checks the computation on the same mice; it is not a replication in new animals.

- `REFACTOR_PLAN.md:935-936` (main): "on the inputs of the approved December 2025 figures"

**37.** **Status.** The team was not fully confident that the effect is robust, and more animals would be a large investment, so the line is paused, not closed.

- `notes/project_scientific_story.md:13`: "the team was not fully confident in their robustness, and since continuing that track needs a big investment"
- `REFACTOR_PLAN.md:35-37` (main): "so the line is paused, not closed"

**38.** The code stays runnable, with identical results, and documented well enough to resume with more animals.

- `REFACTOR_PLAN.md:37-40` (main): "It stays runnable, with identical results, and documented well enough to be resumed with more animals"

## 2. The adult distribution

**39.** **Question.** How the nanobody signal is distributed across the whole adult brain, and how reproducible that is across mice.

- `REFACTOR_PLAN.md:41-44` (main): "How the nanobody signal (in principle surface AMPA receptors containing GluA1) is distributed across the whole brain"

**40.** **Cohort.** The five naive and five RWS mice, pooled.

- `REFACTOR_PLAN.md:43-44` (main): "The adult cohort is the 10 naive and RWS mice pooled"

**41.** Pooling was checked: the part of the map that abundance and density do not explain (line 3) agrees between the naive and the RWS mice at rho +0.88.

- `REFACTOR_PLAN.md:44-46` (main): "agrees between the five naive and the five RWS mice, rho +0.88"
- `notes/project_map_beyond_abundance_density.md:46`: "D naive vs RWS (+0.884)"

**42.** The behaviour mice enter only the plasticity comparison.

- `REFACTOR_PLAN.md:46-47` (main): "behavior adults are not characterised here, only in the plasticity comparison"

**43.** **Reproducibility.** Over the 126 ways of splitting the ten adults into two halves of five, the two half-cohort maps agree at rho 0.974 over 125 grey-matter structures (Spearman-Brown 0.987 for the full cohort).

- `docs/adult_ish_design.md:449-451`: "Over all 126 five-against-five splits of the ten adults, the half-cohort maps agree at **ρ = 0.974**; Spearman-Brown gives **0.987**"
- `notes/project_map_beyond_abundance_density.md:17`: "125 grey-matter structures"

**44.** Every pair of single mice agrees, at a median rho of 0.780 (worst 0.595).

- `notes/project_map_beyond_abundance_density.md:45`: "C single animals (every pair agrees, median 0.780, worst 0.595)"

**45.** The autofluorescence of the same brains explains none of the map (cross-validated R² -0.043).

- `notes/project_map_beyond_abundance_density.md:25`: "autofluorescence, same brains | -0.043, 0%"

**46.** **Being rebuilt.** Which structures stand out, and with what confidence, is being redone.

- `REFACTOR_PLAN.md:433` (main): "A4 | adult distribution: every structure by division, per-mouse mean and SEM, reliability, the enrichment call (S2)"

**47.** P8 called a structure enriched above a threshold on its own scale, and Sami El-Boustani asked (28 April) for a stricter threshold and a permutation null.

- `REFACTOR_PLAN.md:207-210` (main): "P8 called a region enriched at LR-sum >= 1.5 or z >= 0 (zero at the voxel mean); Sami asked (28 Apr) for a threshold of 2 on the LR-sum and, in the same feedback, for a permutation null"

**48.** The plan replaces both with a test per structure against the brain's median structure (S2), on a declared set of structures (S1, A1, A4).

- `REFACTOR_PLAN.md:210-212` (main): "per structure, an exact signed-rank test of zref against the brain's median structure across mice"
- `REFACTOR_PLAN.md:193-195` (main): "S1. The declared structure set"
- `REFACTOR_PLAN.md:430` (main): "A1 | the declared structure set and zref reference (S1)"

**49.** Today 22 structures of the adult table are seen in fewer than five of the ten adults.

- `REFACTOR_PLAN.md:196-198` (main): "22 structures of the adult table are seen in fewer than 5 of the 10 adults"

**50.** P8's outputs carry known defects (listed in the plan) and are not quoted here.

- `REFACTOR_PLAN.md:449-455` (main): "P8's adult outputs carry an artefact of the `abs()` over slabs"

## 3. What the map measures

**51.** A reader will ask whether the nanobody map simply follows how much receptor a region makes, or how many synapses it has.

- `docs/adult_ish_design.md:426-427`: "The two things a reader will say the nano map is"

**52.** The project tests this against the Allen in situ hybridisation (ISH) maps of the adult mouse brain (Lein et al. 2007).

- `REFACTOR_PLAN.md:49-50` (main): "Against the Allen in situ hybridisation maps"
- `ref_papers/Lein_2007.pdf` p. 1: "Genome-wide atlas of gene expression in the adult mouse brain"

**53.** **The map is more than abundance and density.** Over 125 grey-matter structures, the adult map was predicted from receptor abundance (Gria1 to Gria4 mRNA), synaptic markers, the first principal component of 188 postsynaptic-density genes, and the cohort's own autofluorescence.

- `notes/project_map_beyond_abundance_density.md:17`: "125 grey-matter structures, ranks, cross-validated"
- `docs/adult_ish_design.md:461-464`: "abundance (Gria1–4) | 0.253"
- `docs/adult_ish_design.md:462`: "synaptic markers (Syp, Syn1, Vamp2, Bsn, Syt1, Dlg4, Homer1, Shank2/3, Nlgn1, Camk2a)"
- `docs/adult_ish_design.md:463`: "psd_pc1 (first PC of 188 postsynaptic-density genes)"
- `docs/adult_ish_design.md:464`: "autofluorescence, same brains"

**54.** Each model was scored by cross-validation on held-out structures, against the map's own reliability (the ceiling, 97.4% of the variance):

- `docs/adult_ish_design.md:456-457`: "Scored by cross-validation on held-out structures"
- `docs/adult_ish_design.md:451`: "So **97.4% of this map's variance is explainable in principle**"

**55.** | receptor abundance (Gria1 to Gria4) | 0.253 | 26% |

- `notes/project_map_beyond_abundance_density.md:22`: "abundance (Gria1-4) alone | CV R2 0.253, 26% of ceiling"

**56.** | synaptic markers | 0.148 | 15% |

- `notes/project_map_beyond_abundance_density.md:23`: "synaptic markers alone | 0.148, 15%"

**57.** | postsynaptic density, first component | 0.262 | 27% |

- `notes/project_map_beyond_abundance_density.md:24`: "psd_pc1 (188 PSD genes) alone | 0.262, 27%"

**58.** | autofluorescence | -0.043 | 0% |

- `notes/project_map_beyond_abundance_density.md:25`: "autofluorescence, same brains | -0.043, 0%"

**59.** | all four, straight lines | 0.416 | 43% |

- `notes/project_map_beyond_abundance_density.md:26`: "all four, straight | 0.416, 43%"

**60.** | all four, allowed to bend | 0.597 | 61% |

- `notes/project_map_beyond_abundance_density.md:27`: "all four, allowed to bend** | **0.597, 61%**"

**61.** About 39% of the explainable variance is left over.

- `notes/project_map_beyond_abundance_density.md:27`: "leaves **39%**"

**62.** The leftover replicates across independent halves of the cohort at rho 0.934 (Spearman-Brown 0.966), almost as well as the map itself.

- `notes/project_map_beyond_abundance_density.md:28`: "leftover replicates** | **0.934 half, 0.966 Spearman-Brown"
- `docs/adult_ish_design.md:474-476`: "The leftover is very nearly as reproducible as the map itself"

**63.** Seven controls tried to break this: a spatial gradient, structure size, single animals, naive against RWS, curvature, the whole gene space and the choice of reading.

- `docs/adult_ish_design.md:485-493`: "smooth spatial gradient (clearing/illumination artefact)"
- `notes/project_map_beyond_abundance_density.md:43-48`: "A spatial gradient"

**64.** The curvature control failed the first time: fitted with straight lines, a fifth of the map was credited to the leftover.

- `notes/project_map_beyond_abundance_density.md:30-33`: "so a fifth of the map was being credited to the leftover that the covariates could explain"

**65.** The model was changed to allow curvature, and the claim fell from about half to 39%.

- `notes/project_map_beyond_abundance_density.md:33-34`: "Fixed the MODEL, not the wording: headline moved from "half" to 39%"

**66.** The richest model tried, 20 components of the whole 390-gene panel, reaches a cross-validated R² of 0.826 (85% of the ceiling).

- `notes/project_map_beyond_abundance_density.md:47-48`: "F whole 390-gene space (20 components reach CV 0.826 = 85% of ceiling"

**67.** Its leftover still replicates at 0.879.

- `notes/project_map_beyond_abundance_density.md:47-48`: "their leftover STILL replicates at 0.879"

**68.** The leftover is high, relative to what the predictors give, in the medial geniculate (+48 ranks), subthalamic nucleus (+44), ventral lateral geniculate (+38) and lateral habenula (+38).

- `notes/project_map_beyond_abundance_density.md:50-51`: "medial geniculate +48 ranks, subthalamic +44, ventral LGN +38, lateral habenula +38"

**69.** It is low in VPM (-49), VPL (-47), dorsal retrosplenial cortex (-46) and the posterior thalamic complex (-43).

- `notes/project_map_beyond_abundance_density.md:52`: "VPM -49, VPL -47, RSPd -46, posterior complex -43"

**70.** No single gene accounts for it (best, Cacng8, rho +0.301).

- `notes/project_map_beyond_abundance_density.md:52-53`: "Best single gene correlate only +0.301 (Cacng8)"

**71.** This result uses all ten adults and the grey-matter rule, so the structure-set changes of step 9 do not affect it.

- `REFACTOR_PLAN.md:237-238` (main): "The beyond-abundance result is not affected (it already uses all adults and the grey-matter rule)"

**72.** **What it means.** It is consistent with the stain reporting surface GluA1, and it does not show it.

- `REFACTOR_PLAN.md:50-51` (main): "which points to surface GluA1"
- `notes/project_map_beyond_abundance_density.md:62`: "Do NOT claim** the leftover IS the surface fraction"

**73.** A residual is only what the predictors did not explain: regional differences in translation, turnover, subunit composition or nanobody access to the tissue would also land there.

- `notes/project_map_beyond_abundance_density.md:62-64`: "Translation, turnover, subunit composition or nanobody access would all land there"

**74.** **The ISH gene ranking is descriptive.**

- `docs/adult_ish_design.md:535`: "D4 — The ranking is descriptive, and says so"
- `notes/project_ish_panel_v2_powered.md:51-52`: "treating the ranking as descriptive"

**75.** The first comparison (P9, April 2026) correlated the map with 100 hand-picked genes and found AMPA receptor trafficking and anchoring genes, Cacng8 first, above Gria1 itself.

- `notes/project_p8_p9_pipeline.md:21-22`: "P9 — batch comparison vs Allen ISH (100 genes"
- `ref_papers/PAPER_SUMMARIES.md:181`: "our current 100-gene panel is hand-picked"
- `docs/adult_ish_design.md:14-19`: "the nano map already correlates better with AMPAR trafficking and anchoring machinery (Cacng8 = TARP γ-8 at ρ = 0.78, Dlg2 = PSD-93 at 0.76) than with *Gria1* itself"
- `notes/project_p8_p9_pipeline.md:55`: "Surface GluA1 tracks trafficking machinery better than its own mRNA"

**76.** The Python route reproduces that ranking (rho 0.91 against the old ordering).

- `notes/project_adult_ish_v2.md:17-18`: "Spearman between the old MATLAB ordering and the v2 one is **+0.910** for `zref`"

**77.** Cacng8 (TARP γ-8) is first under every reading and every structure set tried (rho 0.77 to 0.82).

- `REFACTOR_PLAN.md:236-237` (main): "Cacng8 first under every structure set tried (rho 0.77 to 0.82)"
- `notes/project_adult_ish_v2.md:20-21`: "Cacng8 is first under every nano reading"
- `docs/adult_ish_design.md:17`: "Cacng8 = TARP γ-8"

**78.** Gria1's own rank depends on which structures enter: 10th of 95 with all structures, 17th when each structure must be seen in at least five adults.

- `notes/project_measurement_validation_strategy.md:101`: "all structures -> Gria1 rank 10 (rho 0.665, 4/33 machinery above); >=5, >=8 or all 10 adults -> rank 17"

**79.** It is not quoted until the declared structure set is in place.

- `REFACTOR_PLAN.md:230-232` (main): "None of the following until A1 to A3 have run: Gria1's rank"

**80.** One Allen map is more reliable than assumed.

- `docs/adult_ish_design.md:348`: "One Allen map is more reliable than we assumed"

**81.** For the 218 genes measured more than once, the median agreement between experiments is rho 0.69, and Gria1's map scores 0.91.

- `notes/project_ish_panel_v2_powered.md:22-23`: "(218 genes have >1 experiment): median **+0.691**"
- `notes/project_ish_panel_v2_powered.md:24-25`: "**Gria1 = 0.91**"

**82.** Gria1 ranking below Cacng8 is not a bad Gria1 experiment.

- `notes/project_ish_panel_v2_powered.md:26-27`: "Gria1 ranking below trafficking genes is NOT a bad experiment"

**83.** Cacng8 above Gria1 is a single-gene result.

- `notes/project_adult_ish_v2.md:100-101`: "Cacng8 +0.80 > Gria1 +0.66 is a SINGLE-GENE result"

**84.** On the 100-gene panel, the subunit genes (median rho +0.58) correlate better with the map than the localisation genes (+0.385).

- `docs/adult_ish_design.md:305-308`: "the subunit median (+0.582"
- `docs/adult_ish_design.md:308`: "localisation median (+0.385"

**85.** Splitting the same 19 genes by function does no better than splitting them at random (exact permutation, p = 0.49).

- `docs/adult_ish_design.md:290-293`: "split the same 19 genes into 4 and 15 every possible way — 3,876 of them"
- `docs/adult_ish_design.md:292-293`: "The functional split lands at **p = 0.4877**"

**86.** **The powered test is negative.** A 390-gene panel was built from Gene Ontology terms, not by hand.

- `notes/project_ish_panel_v2_powered.md:14-15`: "Panel from Gene Ontology, not from us"
- `notes/project_ish_panel_v2_powered.md:19-20`: "390 genes, 669 experiments"

**87.** Once the subunit composite is removed, the 84 AMPA receptor localisation genes explain no more of the map than expression-matched postsynaptic genes: difference -0.010, p = 0.74, where a difference of ±0.062 would have been detected.

- `notes/project_ish_panel_v2_powered.md:17`: "(84)"
- `notes/project_ish_panel_v2_powered.md:35-37`: "difference **-0.010, p = 0.74**. Detectable difference at p<0.05 was ±0.062"

**88.** As a positive control, the same test detects the difference between control genes of high and low map reliability (+0.155, p = 0.0007), so the negative is informative.

- `notes/project_ish_panel_v2_powered.md:37-39`: "control genes split at median reliability give |rho| 0.509 vs 0.354, +0.155, p = 0.0007"
- `docs/adult_ish_design.md:391-392`: "So the test detects a real effect of that size with these sample sizes and does not detect this one"

**89.** Genes such as Arpc5 and Cdk5r1, which have nothing to do with AMPA receptor traffic, sit in the same band as Cacng8.

- `notes/project_ish_panel_v2_powered.md:41-42`: "Top controls by partial rho: Arpc5 +0.602, Cdk5r1 +0.560"
- `notes/project_ish_panel_v2_powered.md:43-44`: "in the SAME BAND as Cacng8"
- `docs/adult_ish_design.md:403-404`: "they have nothing to do with AMPA receptor trafficking"

**90.** The map is predicted about equally well by any well-measured forebrain postsynaptic gene.

- `notes/project_ish_panel_v2_powered.md:46-47`: "the adult nano map is predicted about equally well by any well-measured forebrain postsynaptic gene"

**91.** **Provisional.** These ISH numbers were computed on the unrestricted structure set; A1 to A3 rerun them.

- `REFACTOR_PLAN.md:233-235` (main): "the powered panel test (p = 0.74) and the roles permutation, which run on the same unrestricted structure set"
- `REFACTOR_PLAN.md:432` (main): "rerun of the ranking, roles, panel test, arms and words"

**92.** Their p values are anticonservative, because genes within a set are co-expressed (Fulcher et al. 2021); a spatial null is planned (A7).

- `docs/adult_ish_design.md:179-181`: "genes sharing a GO term are co-expressed"
- `mapping/sepmap/ish/panel_test.py:39-41`: "this remains anticonservative in the way Fulcher 2021 describes"
- `REFACTOR_PLAN.md:436` (main): "A7 | spatial null"

**93.** **The green channel is not total receptor.**

- `notes/project_sep_channel_not_total_receptor.md:3`: "the green/SEP channel in fixed cleared tissue is mostly autofluorescence"

**94.** **The plan.** The SEP tag fluoresces green, so the green channel was meant to report all SEP-GluA1, surface and internal, and nano divided by it a surface fraction.

- `registration/run_add_sep_channel.m:9-12`: "so the green channel is the tagged receptor itself. Ex vivo it reports the whole GluA1 pool"
- `docs/adult_ish_design.md:204`: "`sepratio` = nano / SEP | the surface fraction"

**95.** Three predictions were set in advance: SEP should track Gria1, nano both, and nano/SEP the trafficking genes.

- `docs/adult_ish_design.md:77-81`: "SEP (total receptor) | *Gria1* expression"
- `notes/project_measurement_validation_strategy.md:64-66`: "three arms, each a prediction stated in advance"

**96.** **Measured.** In these fixed, cleared sections the green channel is mostly autofluorescence.

- `docs/adult_ish_design.md:218`: "In this fixed, cleared tissue the green channel is mostly autofluorescence"

**97.** In all ten adults it tracks the autofluorescence channel at rho 0.79 ± 0.04 across structures.

- `notes/project_sep_channel_not_total_receptor.md:20`: "SEP ~ autofluorescence | **+0.791 +- 0.036** (every mouse)"

**98.** It varies less across the brain than either autofluorescence or nano: p90 - p10 of log2, 0.95 ± 0.16 against 1.07 ± 0.16 and 1.93 ± 0.26.

- `notes/project_sep_channel_not_total_receptor.md:19`: "dynamic range p90-p10 log2: nano / autofluo / **SEP** | 1.93+-0.26 / 1.07+-0.16 / **0.95+-0.16**"

**99.** Against Gria1 ISH, nano correlates at +0.60 and the green channel at +0.30.

- `notes/project_sep_channel_not_total_receptor.md:23`: "vs Gria1 ISH: nano / autofluo / SEP / SEP-minus-autofluo | +0.604 / +0.232 / +0.299 / +0.120"

**100.** **Consequences.** `sepratio` (nano per unit SEP) is not a surface fraction, and the three-way test cannot be run on these data.

- `notes/project_sep_channel_not_total_receptor.md:30-33`: "It is NOT a surface fraction"
- `notes/project_sep_channel_not_total_receptor.md:35-36`: "The three-arm validation **cannot be run on this data**"

**101.** **What remains.** With Gria1 mRNA partialled out of the nano map (`ratio` reading), anchoring and trafficking genes still predict the remainder (Cacng8 +0.54, Cnih2 +0.50; median +0.12 over 33 genes, 20 of them positive).

- `notes/project_adult_ish_v2.md:75-76`: "(`ratio` arm, ranks, 33 machinery genes)"
- `notes/project_adult_ish_v2.md:75-78`: "Cacng8 **+0.539**, Cnih2 +0.497, Dlg2 +0.412, Grm5 +0.378; machinery median +0.116, 20/33 positive"
- `notes/project_measurement_validation_strategy.md:58-61`: "with Gria1 partialled out of the nano map, anchoring and trafficking genes still predict the residual"

**102.** This is suggestive, not decisive: mRNA is not protein, and a regional gradient of translation or turnover would look the same.

- `notes/project_measurement_validation_strategy.md:61-62`: "Suggestive, not decisive -- a regional translation or turnover gradient would look the same"
- `notes/project_adult_ish_v2.md:78-80`: "mRNA is not protein"

**103.** **What would settle it.** A wet-lab control: a total-GluA1 antibody stain or autoradiography on a subset of the same brains.

- `notes/project_sep_channel_not_total_receptor.md:40-42`: "a total-GluA1 antibody stain or autoradiography on a subset of the SAME brains"

**104.** A knockout or no-primary control would test specificity instead.

- `notes/project_sep_channel_not_total_receptor.md:42-43`: "or a knockout / no-primary control, which tests specificity rather than surface-vs-total"

**105.** There is no knockout or competition control, no independent regional measure and no DAPI control yet.

- `notes/project_measurement_validation_strategy.md:12-16`: "There is currently no knockout or competition control, no independent regional measure (autoradiography, proteomics), no same-brain replicate, and the DAPI control"
- `REFACTOR_PLAN.md:447` (main): "a DAPI control (in neither route)"

## 4. Young against adult (grant, September 2026)

**106.** **Question.** Young mice (P16 to P36) against adults, looking for differences in cortex that could correspond to critical periods: the grant's Aim 2.1 (below).

- `REFACTOR_PLAN.md:53-55` (main): "Young mice (P16 to P36) against adults, looking for differences in cortex that could correspond to critical periods"

**107.** **Brains.** Seven young brains so far: five at P20 (MG897, MG903, MG909, MG910, MG913), one at P16 (MG911) and one at P22 (MG904), against the ten adults.

- `notes/project_young_vs_adult_comparison.md:102-104`: "young = 7 (MG897/903/909/910/913 P20, MG911 P16, MG904 P22) vs 10 adults"

**108.** MG912 (P20) was dropped for poor slice quality.

- `notes/project_developmental_dataset.md:106-108`: "MG912 (P20) dropped for now** -- Giulio and Sami judged the slice quality too poor to annotate"

**109.** Six more, MG896 and MG914 at P28, MG906 and MG908 at P32, MG895 and MG907 at P36, are in processing and go through the refactored code.

- `notes/project_atlas_decision_ccf_for_all.md:66`: "for MG896/MG914 (P28), MG906/MG908 (P32), MG895/MG907 (P36)"
- `REFACTOR_PLAN.md:117-118` (main): "the remaining young mice (MG896, MG906, MG895, then MG907 and MG908), which then go through the refactored code"
- `REFACTOR_PLAN.md:771-772` (main): "and MG914's SEP channel and Python route with them"

**110.** **Atlases.** Each young brain is registered to the DeMBA atlas of its own age, the adults to the CCF.

- `notes/project_atlas_decision_ccf_for_all.md:13`: "young cohort registers to the age-matched **DeMBA P20**; adults stay on **CCF**"
- `notes/project_atlas_decision_ccf_for_all.md:29`: "One atlas per AGE, not one per cohort"

**111.** Region statistics are computed in each brain's own atlas, whose labels share the CCF ontology, so they need no warping.

- `mapping/sepmap/young_vs_adult/region_plot.py:1-2`: "computed on each brain's OWN atlas -- no warping anywhere"
- `mapping/sepmap/young_vs_adult/region_plot.py:6-7`: "the DeMBA annotations were remapped to Allen parcellation_index"

**112.** Maps carry each young brain to the CCF with the DeMBA deformation of its age.

- `mapping/sepmap/volumes/to_ccf.py:12-14`: "DeMBA canvas -> brainglobe_ccf_translator from age_PND to P56"

**113.** **Only the pattern can be compared.** Young and adult brains were imaged in different sessions.

- `notes/project_measurement_validation_strategy.md:32-33`: "(young and adult were imaged in different sessions)"

**114.** The autofluorescence used as an internal standard rises with age itself.

- `mapping/sepmap/volumes/cohort.py:15-21`: "the internal standard"
- `mapping/sepmap/volumes/cohort.py:21`: "autofluorescence rises with age"

**115.** Nano per unit autofluorescence is lower in young cortex (-0.4 to -1.0 log2), but that is a bound, not a value.

- `notes/project_young_vs_adult_comparison.md:130`: "`ratio` in cortex is now -0.4 to -1.0"
- `mapping/sepmap/volumes/cohort.py:23`: "Treat it as a bound, not a value"

**116.** **The young brain is flatter.** The spread of its structures (p90 - p10 of log2 nano relative to cortex) is 0.98 log2, against 1.93 in adults.

- `notes/project_zref_difference_is_already_relative.md:31`: "young (n=7) | **0.98** log2"
- `notes/project_zref_difference_is_already_relative.md:32`: "adult (n=10) | **1.93** log2"
- `notes/project_zref_difference_is_already_relative.md:34-35`: "the zref difference is a difference of *standardised* log positions"

**117.** `zref` removes this per brain, so a zref difference is a difference of positions within each brain's own range, not a fold change (see Terms).

- `notes/project_young_vs_adult_comparison.md:180-182`: "can never be converted to a fold change"

**118.** **By system** (7 young against 10 adults, `cref`, log2 young minus adult, group medians; Mann-Whitney p < 0.01, uncorrected, unless marked): somatosensory +0.21, retrosplenial +0.24, frontal -0.42, striatum -0.79, hippocampus -1.03; primary visual -0.06, not significant.

- `notes/project_young_vs_adult_comparison.md:105-106`: "somatosensory +0.21**, retrosplenial +0.24**, frontal -0.42**, striatum -0.79**, hippocampus -1.03**, primary visual -0.06"
- `notes/project_young_vs_adult_comparison.md:88`: "By system (cref)"
- `notes/project_young_vs_adult_comparison.md:126-128`: "the grouped result is unchanged from 23 Sep"
- `notes/project_young_vs_adult_comparison.md:134-135`: "`group_stats` a difference of MEDIANS"
- `mapping/sepmap/young_vs_adult/region_groups.py:379`: "('**' if r['mannwhitney_p'] < 0.01"

**119.** **RL and AL.** RL and AL together, the visuo-tactile areas of the grant, are the only visual group that survives correction: `zref` +0.23 (p = 0.0007, q = 0.014), against -0.01 between naive and RWS adults; the P20 brains alone give +0.34.

- `notes/project_young_vs_adult_comparison.md:146-148`: "**RL+AL is the only visual group that survives correction** -- zref +0.23, p = 0.0007, q = 0.014, naive-vs-rws null -0.01, P20-only +0.34"

**120.** V1 shows no difference (`zref` +0.07, q = 0.35).

- `notes/project_young_vs_adult_comparison.md:148`: "V1 is the null (+0.07, q = 0.35)"

**121.** **By layer.** In every sensory system, primary and higher order alike, the supragranular layers hold a larger share of the signal in young brains (+0.15 to +0.32) and the infragranular layers the same or a smaller one (-0.13 to -0.22).

- `notes/project_young_vs_adult_comparison.md:110-111`: "supragranular up (+0.15 to +0.32**), infragranular flat or down (-0.13 to -0.22)"

**122.** In RL+AL both the supragranular layers (+0.28) and layer 4 (+0.27) are higher (q = 0.006 each), while layer 4 is flat in V1 and in the other visual areas: RL and AL have a laminar profile closer to somatosensory cortex than to V1.

- `notes/project_young_vs_adult_comparison.md:148-152`: "RL+AL supragranular +0.28 (q = 0.006) AND granular L4 +0.27 (q = 0.006), infragranular ns"
- `notes/project_young_vs_adult_comparison.md:150-151`: "RL/AL have a somatosensory-like laminar profile rather than a V1-like one"

**123.** **The most robust number.** The supragranular-minus-infragranular contrast is positive in every system under all five readings (+0.09 to +0.44).

- `notes/project_young_vs_adult_comparison.md:154-156`: "positive in every system under ALL FIVE readings (+0.09 to +0.44)"

**124.** It cancels any scale factor of a brain (exposure, staining strength, the choice of denominator).

- `notes/project_young_vs_adult_comparison.md:156-157`: "A supra-minus-infra contrast cancels any per-brain scalar -- session exposure, staining strength"

**125.** **A critical period cannot be claimed yet.** Two ages cannot show a window.

- `notes/project_young_vs_adult_comparison.md:160-161`: "Two time points cannot show a window"

**126.** V1, where the classic mouse critical period sits (Levelt & Hübener 2012), shows no difference at P16 to P22, while somatosensory cortex shows the largest.

- `notes/project_young_vs_adult_comparison.md:161-162`: "V1, the canonical CP area, is the null while somatosensory"
- `ref_papers/Levelt_2012.pdf` p. 3: "the time between postnatal day 28 and postnatal day 32 is generally con- sidered as the critical period"

**127.** And the superficial-layer excess is also what the later maturation of the superficial layers alone would give.

- `notes/project_young_vs_adult_comparison.md:162-163`: "inside-out maturation (L2/3 is simply younger at P16-P22) predicts the same gradient"

**128.** The P28 to P36 brains can separate the two: if V1's superficial excess peaks near P28 and falls by P36 while S1 has already declined, that is a window moving across modalities; the same decline everywhere at the same rate is maturation.

- `notes/project_young_vs_adult_comparison.md:164-167`: "if V1's supragranular excess peaks there and falls by P36 while S1 has already declined, that is the CP sweeping across modalities"

**129.** **Registration caveats.** Registration is driven by DAPI, whose cell packing falls between P20 and adulthood.

- `notes/project_developmental_dataset.md:63-65`: "DAPI cell-packing density falls between P20 and adult"

**130.** In the P20 atlas, 19% of the volume that is fibre tract in the adult atlas is labelled isocortex, so results next to white matter need care.

- `notes/project_atlas_decision_ccf_for_all.md:24`: "19% of adult fiber-tract volume is labelled Isocortex at P20"
- `notes/project_atlas_decision_ccf_for_all.md:24`: "makes white-matter-adjacent results suspect"

**131.** **Grant figure.** The grant's Figure 3b to d shows this comparison: maps of pups (n = 7) and adults (n = 10), a flatmap of the pup-minus-adult difference, and Mann-Whitney tests on V1, RL+AL, the other higher visual areas, S1, S2 and medial prefrontal cortex.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 13: "Atlas-registered maps of normalized surface GluA1 in pups (n=7 mice aged between P16 and P20) and adults (n=10 mice)"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 13: "Flattened cortical projection of the pup-minus-adult difference"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 13: "across V1, RL+AL, non visuo-tactile higher visual areas (HVA), S1, S2 and medial prefrontal cortex (mPFC). Mann-Whitney test"

**132.** The figure set chosen for the grant includes the layer 5-6 flatmap.

- `notes/project_young_vs_adult_comparison.md:195-196`: "the RL-level slice pair (young | adult | difference), the L5-6 flatmap"

**133.** That band is the weakest of the three statistically (RL+AL infragranular q = 0.110, against 0.006 for the other two), and it looks strong partly because the adult baseline is high there.

- `notes/project_young_vs_adult_comparison.md:196-199`: "RL+AL infragranular q = 0.110 against 0.006 for supra and granular"
- `notes/project_young_vs_adult_comparison.md:199`: "partly because the adult baseline is high there"

**134.** **Provisional.** A1 changes the reference of zref and so moves every zref value, the differences behind the grant figures included.

- `REFACTOR_PLAN.md:441-442` (main): "A1 moves every zref, so the young-against-adult differences behind the grant figures"

**135.** The code state behind the grant figures is tagged `grant-2026-09`.

- `REFACTOR_PLAN.md:636-637` (main): "Tag the code state behind the grant figures as `grant-2026-09`"

## The grant: SNSF Weave (El-Boustani, Geneva; Gjorgjieva, Munich)

**136.** *Dendritic plasticity rules shaping the emergence of multisensory integration in the developing cortex.*

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 1: "Dendritic plasticity rules shaping the emergence of multisensory integration in the developing cortex"

**137.** **Central hypothesis.** Activity-dependent dendritic plasticity rules organise visual and tactile synaptic inputs along dendrites during a postnatal window, and so shape the nonlinear integration behind supramodal visuo-tactile representations in adult cortex.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 1: "activity-dependent dendritic plasticity rules organize visual and tactile synaptic inputs along dendrites during a postnatal developmental window, thereby shaping the nonlinear integration mechanisms that underlie supramodal visuo-tactile representations in adult cortex"

**138.** **Three aims.** Aim 1, how dendritic organisation supports multisensory integration in adult associative cortex.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 4: "Aim 1: Determine how dendritic organization supports multisensory integration in adult associative cortex"

**139.** Aim 2, the developmental sequence and plasticity window in which visuo-tactile integration emerges.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 4: "Aim 2: Define developmental sequence and plasticity window during which visuo-tactile integration emerges"

**140.** Aim 3, the dendritic plasticity rules that build it.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 4: "Aim 3: Identify dendritic plasticity rules that generate branch-specific visuo-tactile organization"

**141.** The two-photon imaging repository (`D:\dendrites\code`, its `docs/SCIENTIFIC_CONTEXT.md`) holds the dendritic and functional imaging side; this repository is Aim 2.1.

- `D:/dendrites/code/docs/SCIENTIFIC_CONTEXT.md:102-104`: "**Aim 2, when.** Map a developmental window of plasticity potential with surface GluA1"
- `D:/dendrites/code/docs/SCIENTIFIC_CONTEXT.md:262`: "development: pup somata (Aim 2) and pup dendrites (Aim 3)"
- `REFACTOR_PLAN.md:732-734` (main): "the imaging repository describes this project as part of Aim 2 of the SNSF Weave grant"

**142.** **Aim 2.1.** *Map developmental windows of synaptic plasticity potential using surface GluA1.*

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 9: "Aim 2.1. Map developmental windows of synaptic plasticity potential using surface GluA1"

**143.** Rationale: AMPA receptor trafficking is a central mechanism of synaptic strengthening, and GluA1 plays an important role in potentiation.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 9: "AMPA receptor trafficking is a central mechanism of synaptic strengthening and experience-dependent plasticity. In particular, the GluA1 subunit plays an important role in synaptic"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "potentiation and long-term plasticity"

**144.** Total GluA1 immunostaining and Gria1 mRNA do not report the pool of receptors at the surface.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "Existing approaches based on total GluA1 immunostaining or Gria1 mRNA expression do not directly report the pool of receptors available at the neuronal surface"

**145.** Design: SEP-GluA1 mice at stages such as P14, P18, P22, P26, P30, P34 and adult; non-detergent processing and a GFP-booster nanobody; registration with LightSuite, juvenile brains first to age-specific templates; V1, whisker S1 and the associative areas RL, AL and LI at area and layer resolution; hippocampus, amygdala and striatum as internal controls.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "such as P14, P18, P22, P26, P30, P34 and adult"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "Juvenile brains will be aligned using age-specific templates"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "across primary visual cortex, whisker primary somatosensory cortex, and associative regions including RL, AL and LI"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "hippocampus, amygdala and striatum will be included as internal controls"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "area- and layer-specific resolution"

**146.** Main output: a developmental plasticity-potential index, surface GluA1 normalised for area, layer, age, section quality, staining batch and imaging conditions.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "developmental plasticity-potential index derived from surface GluA1 signal normalized for area, layer, age, section quality, staining batch and imaging conditions"

**147.** Key prediction: associative visuo-tactile areas show a delayed or prolonged window of elevated surface GluA1, compared with primary sensory cortices.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "associative areas involved in visuo-tactile integration will show a delayed or prolonged window of elevated surface GluA1 compared with primary sensory cortices"

**148.** Risks: surface GluA1 may not map one-to-one onto synaptic plasticity, so it is read as a proxy for plasticity potential, not a definitive marker of critical periods.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 12: "surface GluA1 may not map one-to-one onto synaptic plasticity"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 12: "interpret surface GluA1 as a proxy for plasticity potential, not as a definitive marker of critical periods"

**149.** If fixed sections show too weak a cortical signal, the histological atlas would be complemented by in vivo bulk imaging of SEP-GluA1 from about P15 to P35.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 12: "does not reveal sufficiently strong signals in the cortex, we will complement the histological atlas with longitudinal bulk in vivo imaging of SEP-GluA1 fluorescence from approximately P15 to P35"

**150.** **Preliminary data.** Non-detergent labelling with no detectable intracellular signal (Fig. 3a), and region-specific pup-adult differences, including in RL and AL (Fig. 3b to d, the comparison of line 4).

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "with no detectable intracellular labeling in preliminary preparations (Fig. 3a)"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "Atlas-registered comparisons of pups and adult brains show region-specific differences, including in visuo-tactile areas RL/AL"
- `notes/project_scientific_story.md:21`: "Fig 3b-d = pups n=7 (P16-P20) vs adults n=10"

**151.** **Expected outputs.** Among at least three, a manuscript or resource paper on developmental surface GluA1 mapping and plasticity windows.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 17: "We anticipate at least three main outputs"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 17: "a manuscript or resource paper on developmental surface GluA1 mapping and plasticity windows"

**152.** The grant states that processed datasets, analysis code and simulation code will be released upon publication, when ethical and legal requirements allow.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 17: "Processed datasets, analysis code and simulation code will be released upon publication whenever compatible with ethical and legal requirements"

**153.** **Where the code stands.** The young brains in hand are P16, P20, P22, P28, P32 and P36, not the grant's planned series.

- `common/get_cohort.m:75-91`: "'MG911_SepGluA_P16', 'young', 16"
- `common/get_cohort.m:75-91`: "'MG895_SepGluA_P36', 'young', 36"
- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "such as P14, P18, P22, P26, P30, P34 and adult"

**154.** The code does not compute the plasticity-potential index; the readings below are what exists.

- search of this branch (mapping, group_comparison, preprocessing, registration, atlas, common) for `(?i)potential.index|plasticity.potential`: no plasticity-potential index in the code

## Related literature

**155.** The eleven papers in `data\ref_papers\`, grouped by the question they bear on.

- listing of `data\ref_papers\`: 12 PDFs in data\ref_papers, one of them the grant

**156.** For each: what it shows, and what it means here.

- editorial: describes this document

**157.** Papers the grant cites, but that are not in the folder, are named through the grant.

- editorial: describes this document

## What surface GluA1 reports

**158.** **Huganir RL, Nicoll RA (2013).** AMPARs and synaptic plasticity: the last 25 years. *Neuron* 80:704-717. `Huganir_2013.pdf`

- `ref_papers/Huganir_2013.pdf` p. 1: "AMPARs and Synaptic Plasticity: The Last 25 Years Richard L. Huganir"
- `ref_papers/Huganir_2013.pdf` p. 1: "704 Neuron 80, October 30, 2013"
- `ref_papers/Huganir_2013.pdf` p. 14: "Neuron 80, October 30, 2013 ª2013 Elsevier Inc. 717"

**159.** Review.

- `ref_papers/Huganir_2013.pdf` p. 1: "Here we review the progress over the last two and a half decades"

**160.** Regulating AMPA receptor function and membrane trafficking is critical for many forms of synaptic plasticity; GFP-tagged GluA1 is recruited to spines after LTP; a sizeable surface pool of receptors is a recurring requirement of LTP models.

- `ref_papers/Huganir_2013.pdf` p. 1: "membrane trafﬁcking is critical for many forms of synaptic plasticity"
- `ref_papers/Huganir_2013.pdf` p. 3: "GFP-GluA1 was recruited to synaptic spines after LTP induction"
- `ref_papers/Huganir_2013.pdf` p. 5: "is a recent recurring theme in current models of LTP"

**161.** Auxiliary proteins steer the receptors: TARPs (type I: γ-2, γ-3, γ-4, γ-8) ensure their delivery to the surface and to synapses.

- `ref_papers/Huganir_2013.pdf` p. 6: "type I ( g-2, g-3, g-4, and g-8)"
- `ref_papers/Huganir_2013.pdf` p. 6: "ensure the proper maturation and delivery of AMPARs to the neuron’s surface and synapses"

**162.** Deleting cornichons 2 and 3 causes a selective loss of surface GluA1-containing receptors in the hippocampus, a selectivity that appears to be mediated by TARP γ-8.

- `ref_papers/Huganir_2013.pdf` p. 6: "Genetic dele- tion of CNIH-2 and -3 together causes a profound and selective loss of synaptic and surface AMPARs in the hippocampus"
- `ref_papers/Huganir_2013.pdf` p. 6: "This deﬁcit is due to the selective loss of surface GluA1-containing AMPARs"
- `ref_papers/Huganir_2013.pdf` p. 6: "The remarkably selective effect of CNIHs on the GluA1 subunit appears to be mediated by TARP g-8"

**163.** Earlier work found GluA1 required for LTP; molecular replacement experiments (Granger et al. 2013) found no specific subunit requirement.

- `ref_papers/Huganir_2013.pdf` p. 5: "with GluA1 being required for LTP and GluA2 being required for LTD"
- `ref_papers/Huganir_2013.pdf` p. 5: "have shown no speciﬁc subunit requirement for LTP"

**164.** Here: why surface GluA1 is a plausible plasticity readout, and why not a specific one, since LTP may not need GluA1 in particular.

- `ref_papers/Huganir_2013.pdf` p. 5: "have shown no speciﬁc subunit requirement for LTP"

**165.** It also names the proteins behind the top of the ISH ranking: Cacng8 is TARP γ-8, Cnih2 is cornichon-2.

- `docs/adult_ish_design.md:17`: "Cacng8 = TARP γ-8"
- `ref_papers/Huganir_2013.pdf` p. 6: "These include cornichon-2 and -3 (CNIH-2 and CNIH-3)"

**166.** Since localisation genes as a class explain no more of the map than other postsynaptic genes, their place at the top is not evidence about trafficking.

- `notes/project_ish_panel_v2_powered.md:46-48`: "Cacng8 topping the ranking is real and reproducible and is NOT evidence that the map is about trafficking"

**167.** **Lopez-Ortega E, Choi JY, Hong I, et al. (2024).** Stimulus-dependent synaptic plasticity underlies neuronal circuitry refinement in the mouse primary visual cortex. *Cell Rep* 43:113966. `Lopez-Ortega_2024.pdf`

- `ref_papers/Lopez-Ortega_2024.pdf` p. 1: "Lopez-Ortega et al., 2024, Cell Reports 43, 113966"
- `ref_papers/Lopez-Ortega_2024.pdf` p. 1: "Stimulus-dependent synaptic plasticity underlies neuronal circuitry reﬁnement in the mouse primary visual cortex"

**168.** In vivo two-photon imaging of V1 layer 2/3 over five days of repeated grating stimulation (spines imaged in mice of 2.5 to 3 months).

- `ref_papers/Lopez-Ortega_2024.pdf` p. 2: "we monitored the neuronal responses of V1 layer 2/3 neurons"
- `ref_papers/Lopez-Ortega_2024.pdf` p. 2: "exposed head-ﬁxed mice to visual stimulation for 5 consecutive"
- `ref_papers/Lopez-Ortega_2024.pdf` p. 5: "mice at 2.5–3 months of age"

**169.** Fewer neurons respond (44% to 32% of active neurons per field), and the persistent ones respond 30% more.

- `ref_papers/Lopez-Ortega_2024.pdf` p. 4: "44% of active neurons were responsive in the visually stim- ulated mice on day 1, while only 32% were responsive by day 5"
- `ref_papers/Lopez-Ortega_2024.pdf` p. 4: "persistent neurons presented a 30% increase in their onset response"

**170.** Spine density falls by 15%.

- `ref_papers/Lopez-Ortega_2024.pdf` p. 5: "overall reduction in spine density of 15% compared to baseline"

**171.** Spine GluA1, read as SEP-GluA1 fluorescence after in utero electroporation, rises by about 20% and stays up for at least 48 h; potentiated spines are clustered.

- `ref_papers/Lopez-Ortega_2024.pdf` p. 4: "E15 mouse embryos were in utero electroporated"
- `ref_papers/Lopez-Ortega_2024.pdf` p. 7: "in spine GluA1 levels in the apical dendrites of L2/3 V1 neurons after visual stim- ulation that was maintained for at least 48 h"
- `ref_papers/Lopez-Ortega_2024.pdf` p. 1: "Potentiated spines are found in well-deﬁned clusters within dendrites"

**172.** Here: the kind of change the plasticity line looked for.

- `REFACTOR_PLAN.md:33-35` (main): "Surface GluA1 after an experience"

**173.** Repeated sensory stimulation over days raises spine GluA1 in the stimulated primary cortex of adult mice.

- `ref_papers/Lopez-Ortega_2024.pdf` p. 2: "an overall increase in spine AMPA receptor levels in the same sub- set of neurons"

**174.** The paper reads SEP-GluA1 fluorescence in living tissue; in our fixed sections the green channel is mostly autofluorescence, so that reading is not available here.

- `ref_papers/Lopez-Ortega_2024.pdf` p. 5: "were longitudinally imaged in awake mice"
- `ref_papers/Lopez-Ortega_2024.pdf` p. 7: "we used SEP-GluA1 signal intensity as a proxy for studying changes in AMPAR levels"
- `notes/project_sep_channel_not_total_receptor.md:3`: "the green/SEP channel in fixed cleared tissue is mostly autofluorescence"

## The plasticity experiment

**175.** **Gambino F, Pagès S, Kehayas V, et al. (2014).** Sensory-evoked LTP driven by dendritic plateau potentials in vivo. *Nature* 515:116-119. `Gambino_2014.pdf`

- `ref_papers/Gambino_2014.pdf` p. 1: "Sensory-evoked L TP driven by dendritic plateau potentials in vivo"
- `ref_papers/Gambino_2014.pdf` p. 1: "116 | NATURE | VOL 515"
- `ref_papers/Gambino_2014.pdf` p. 4: "VOL 515 | NATURE | 119"

**176.** Whole-cell recordings in layer 2/3 of the mouse barrel cortex, under urethane.

- `ref_papers/Gambino_2014.pdf` p. 1: "Whole-cell patch recordings were targeted to cells above the C2 barrel of urethane-anaesthetized mice"

**177.** Rhythmic whisker stimulation (RWS: the principal whisker deflected at 8 Hz for 1 min) potentiates the short-latency whisker-evoked PSP with no somatic spike (n = 11, P = 0.008).

- `ref_papers/Gambino_2014.pdf` p. 1: "deﬂected back and forth (100-ms deﬂections) for 1 min at a frequency of 8 Hz"
- `ref_papers/Gambino_2014.pdf` p. 1: "None of the recorded cells displayed somatic action potentials during RWS"
- `ref_papers/Gambino_2014.pdf` p. 1: "short-latency PSP amplitudes (PSPshort"
- `ref_papers/Gambino_2014.pdf` p. 1: "n 5 11, **P 5 0.008"

**178.** The potentiation needs NMDA-receptor plateau potentials, and these depend on the posteromedial thalamic nucleus (POm).

- `ref_papers/Gambino_2014.pdf` p. 1: "The induction of LTP depended on the occurrence of NMDAR"
- `ref_papers/Gambino_2014.pdf` p. 1: "inhibition of POm activity during rhythmic whisker stimulation suppressed the generation of those potentials and prevented whisker-evoked LTP"

**179.** Here: the protocol behind the RWS group, and the reason S1 was where its effect was expected.

- `REFACTOR_PLAN.md:34` (main): "rhythmic whisker stimulation, RWS, as in Gambino et al. 2014"
- `notes/project_scientific_story.md:20`: "as expected for plasticity localised in S1 caused by RWS"

**180.** The plasticity comparison's regional figures highlight the barrel field, VPM and the posterior thalamic complex.

- `group_comparison/pipeline/group_differences.m:1207`: "case 'rws', highlighted_areas = {'Primary somatosensory area, barrel field', 'Ventral posteromedial nucleus of the thalamus', 'Posterior complex of the thalamus'"

**181.** The paper measures potentiation of single cells over minutes in anaesthetised mice; it does not measure surface GluA1.

- `ref_papers/Gambino_2014.pdf` p. 1: "This LTP lasted for as long as the cells could be recorded from (at least 15 min after RWS)"

## Critical periods

**182.** **Levelt CN, Hübener M (2012).** Critical-period plasticity in the visual cortex. *Annu Rev Neurosci* 35:309-330. `Levelt_2012.pdf`

- `ref_papers/Levelt_2012.pdf` p. 1: "Critical-Period Plasticity in the Visual Cortex Christiaan N. Levelt"
- `ref_papers/Levelt_2012.pdf` p. 1: "Annu. Rev. Neurosci. 2012. 35:309–30"

**183.** Review of ocular-dominance plasticity, mostly in rodent V1.

- `ref_papers/Levelt_2012.pdf` p. 1: "Work in the rodent visual cortex has led to important insights"
- `ref_papers/Levelt_2012.pdf` p. 1: "ocular dominance, monocular deprivation"

**184.** In mice it is strongest at the end of the fourth postnatal week; P28 to P32 is generally taken as the critical period.

- `ref_papers/Levelt_2012.pdf` p. 3: "The effect is strongest during a relatively brief phase at the end of the fourth postnatal week, and the time between postnatal day 28 and postnatal day 32 is generally con- sidered as the critical period for OD plasticity in the mouse visual cortex"

**185.** The maturation of inhibitory circuits opens it; what ends it is less clear.

- `ref_papers/Levelt_2012.pdf` p. 1: "the maturation of speciﬁc inhibitory circuits plays a key role in the opening of the critical period in the visual cortex, it is less clear what puts an end to it"

**186.** Plasticity exists well before and long after the peak, and closure is often gradual.

- `ref_papers/Levelt_2012.pdf` p. 1: "plasticity in the visual cortex is present well before, and long after, the peak of the critical period"
- `ref_papers/Levelt_2012.pdf` p. 2: "the closure of a critical period is frequently a gradual process"

**187.** Here: the reference timing for V1.

- editorial: states how the paper is used here

**188.** The young brains analysed so far are P16 to P22, before that peak, and V1 shows no young-adult difference; the P28 brains sit at the peak.

- `notes/project_young_vs_adult_comparison.md:102-104`: "young = 7 (MG897/903/909/910/913 P20, MG911 P16, MG904 P22)"
- `notes/project_young_vs_adult_comparison.md:148`: "V1 is the null"
- `common/get_cohort.m:76-91`: "'MG896_SepGluA_P28', 'young', 28"

**189.** The mechanisms it reviews are inhibitory.

- `ref_papers/Levelt_2012.pdf` p. 1: "the maturation of speciﬁc inhibitory circuits plays a key role"

**190.** The grant treats surface GluA1 as a proxy for plasticity potential, not a definitive marker of critical periods.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 12: "interpret surface GluA1 as a proxy for plasticity potential, not as a definitive marker of critical periods"

**191.** **Larsen B, Sydnor VJ, Keller AS, Yeo BTT, Satterthwaite TD (2023).** A critical period plasticity framework for the sensorimotor-association axis of cortical neurodevelopment. *Trends Neurosci* 46(10):847-862. `Larsen_2023.pdf`

- `ref_papers/Larsen_2023.pdf` p. 1: "A critical period plasticity framework for the sensorimotor–association axis of cortical neurodevelopment Bart Larsen"
- `ref_papers/Larsen_2023.pdf` p. 1: "Trends in Neurosciences, October 2023, Vol. 46, No. 10"
- `ref_papers/Larsen_2023.pdf` p. 16: "862 Trends in Neurosciences, October 2023"

**192.** Review, mostly of human imaging.

- `ref_papers/Larsen_2023.pdf` p. 1: "Human neuroimaging studies have provided evidence for a hierarchical sensorimotor-to-association"

**193.** Development proceeds along a sensorimotor-to-association axis, and association cortex matures last.

- `ref_papers/Larsen_2023.pdf` p. 1: "brain development progresses hierarchically along an S –A axis in which areas of association cortex are the last to mature"

**194.** In animal models, critical periods progress along cortical hierarchies: opened by maturing parvalbumin interneurons and a falling excitation/inhibition ratio, closed by myelin and perineuronal nets.

- `ref_papers/Larsen_2023.pdf` p. 1: "Animal models have been used to identify periods of enhanced experience-dependent plasticity – 'critical periods' – that progress along cortical hierarchies"
- `ref_papers/Larsen_2023.pdf` p. 3: "initial strengthening of parvalbumin (PV) interneuron cell signaling and declines in the excitation/inhibition (E/I) ratio"
- `ref_papers/Larsen_2023.pdf` p. 3: "including the formation of cortical myelin and perineural nets"

**195.** Parvalbumin-cell development starts in primary cortex and cascades to higher-order sensory areas (shown in primates).

- `ref_papers/Larsen_2023.pdf` p. 7: "starting a cascade of PV cell development along the cortical hierarchy to higher-order sensory cortex (e.g., areas V2 to TE to 7a in primates)"

**196.** Rodent prefrontal cortex is affected by deprivation, enrichment and stress at later stages, which suggests a later or prolonged critical period.

- `ref_papers/Larsen_2023.pdf` p. 5: "the development of the rodent prefrontal cortex is affected by environmental deprivation, environmental enrichment, and elevated stress when experienced in later developmental stages"
- `ref_papers/Larsen_2023.pdf` p. 5: "suggesting a later or prolonged critical period"

**197.** Here: the frame of the grant's prediction that associative visuo-tactile areas have a delayed or prolonged window.

- `ref_papers/SNSF_Weave_dendritic_multisensory_final.pdf` p. 10: "will show a delayed or prolonged window of elevated surface GluA1"

**198.** Our frontal areas are lower in young brains (`cref` -0.42); with two ages this cannot be read as timing.

- `notes/project_young_vs_adult_comparison.md:105`: "frontal -0.42**"
- `notes/project_young_vs_adult_comparison.md:160-161`: "Two time points cannot show a window"

## Comparing a brain map with gene expression

**199.** **Lein ES, Hawrylycz MJ, Ao N, et al. (2007).** Genome-wide atlas of gene expression in the adult mouse brain. *Nature* 445:168-176. `Lein_2007.pdf`

- `ref_papers/Lein_2007.pdf` p. 1: "Genome-wide atlas of gene expression in the adult mouse brain Ed S. Lein"
- `ref_papers/Lein_2007.pdf` p. 1: "Vol 445 | 11 January 2007 | doi:10.1038/nature05453 168"
- `ref_papers/Lein_2007.pdf` p. 9: "NATURE | Vol 445 | 11 January 2007 176"

**200.** The Allen Mouse Brain Atlas: in situ hybridisation for about 20,000 genes in 56-day-old male C57BL/6J mice.

- `ref_papers/Lein_2007.pdf` p. 1: "containing the expression patterns of ,20,000 genes in the adult mouse brain"
- `ref_papers/Lein_2007.pdf` p. 8: "Male, 56-day-old C57BL/6J mice"

**201.** Each gene was run on sagittal sections 200 µm apart, with a coronal replicate for about 3,500 genes.

- `ref_papers/Lein_2007.pdf` p. 8: "Each riboprobe was initially processed on 200- mm spaced sagittal sections spanning an entire hemi- sphere, and a coronal replicate was generated for a subset (,3,500) of genes"

**202.** The images are registered to a reference atlas and summarised by structure or on a 100 to 300 µm grid.

- `ref_papers/Lein_2007.pdf` p. 2: "through registration with a reference atlas created for the project"
- `ref_papers/Lein_2007.pdf` p. 2: "searchable either by anatomical structure or using a finer (100 mm to 300 mm) resolution geometric grid"

**203.** About 80% of the genes show some expression in the brain.

- `ref_papers/Lein_2007.pdf` p. 3: "Approximately 80% of total genes assayed display some cellular expression above background in the brain"

**204.** Here: the source of the gene maps, used as expression energy on the 200 µm grid, downloaded per experiment.

- `docs/adult_ish_design.md:593`: "Allen Mouse Brain Atlas ISH, 200 µm grid, expression *energy*"
- `mapping/sepmap/ish/panel_fetch.py:4`: "One zip per experiment from `api.brain-map.org/grid_data/download/<id>`"

**205.** Each experiment is one mouse, and it measures mRNA, not protein.

- `notes/project_p8_p9_pipeline.md:58`: "ISH = single mouse per gene"
- `docs/adult_ish_design.md:515-516`: "mRNA is an imperfect proxy for protein"

**206.** Hence the reliability of each map was measured (median 0.69 over 218 genes with repeats), and a residual cannot be called "surface".

- `notes/project_ish_panel_v2_powered.md:22-23`: "(218 genes have >1 experiment): median **+0.691**"
- `notes/project_map_beyond_abundance_density.md:62`: "Do NOT claim** the leftover IS the surface fraction"

**207.** **Fulcher BD, Arnatkeviciute A, Fornito A (2021).** Overcoming false-positive gene-category enrichment in the analysis of spatially resolved transcriptomic brain atlas data. *Nat Commun* 12:2669. `Fulcher_2021.pdf`

- `ref_papers/Fulcher_2021.pdf` p. 1: "Overcoming false-positive gene-category enrichment in the analysis of spatially resolved transcriptomic brain atlas data Ben D. Fulcher"
- `ref_papers/Fulcher_2021.pdf` p. 1: "NATURE COMMUNICATIONS | (2021) 12:2669"

**208.** Asking which gene categories correlate with a brain map, by shuffling genes between categories, is badly anticonservative: genes in a category are co-expressed, and brain maps are spatially autocorrelated.

- `ref_papers/Fulcher_2021.pdf` p. 1: "We show that within-category gene–gene coexpression and spatial autocorrelation are key drivers of the false-positive bias"
- `ref_papers/PAPER_SUMMARIES.md:158-160`: "the conventional test randomises **which genes belong to which category**"

**209.** With random phenotype maps and real mouse expression data, the mean category false-positive rate rose 875-fold.

- `ref_papers/Fulcher_2021.pdf` p. 4: "the mean CFPR increased 875-fold in mouse"

**210.** The fix is to compare against ensembles of surrogate maps, including spatially autocorrelated ones (for the mouse brain, rho = 0.8 and d0 = 1.46 mm).

- `ref_papers/Fulcher_2021.pdf` p. 1: "introduce ﬂexible ensemble-based null models"
- `ref_papers/Fulcher_2021.pdf` p. 3: "spatially autocorrelated SBPs, labeled “SBP-spatial”"
- `ref_papers/Fulcher_2021.pdf` p. 10: "We set ρ = 0.8"
- `ref_papers/Fulcher_2021.pdf` p. 10: "yielding d0 = 1.46 mm in mouse brain"

**211.** Here: why the ISH ranking is descriptive, why its p and q values are anticonservative, and the design of the planned spatial null (A7).

- `docs/adult_ish_design.md:548-551`: "measured the conventional gene-category null at **875-fold** false-positive inflation"
- `REFACTOR_PLAN.md:436` (main): "A7 | spatial null"

**212.** **Koopmans F, van Nierop P, Andres-Alonso M, et al. (2019).** SynGO: an evidence-based, expert-curated knowledge base for the synapse. *Neuron* 103:217-234. `Koopmans_2019.pdf`

- `ref_papers/Koopmans_2019.pdf` p. 1: "Koopmans et al., 2019, Neuron 103, 217–234"
- `ref_papers/Koopmans_2019.pdf` p. 2: "Frank Koopmans,1,2 Pim van Nierop, 2 Maria Andres-Alonso"

**213.** 1,112 genes and 2,922 annotations to 87 synaptic locations and 179 synaptic processes, each based on published experimental evidence.

- `ref_papers/Koopmans_2019.pdf` p. 2: "87 synaptic locations and 179 synaptic processes. SynGO annotations are exclusively based on published, expert-curated evidence. Using 2,922 annotations for 1,112 genes"

**214.** Synaptic genes are compared against brain-expressed control sets.

- `ref_papers/Koopmans_2019.pdf` p. 7: "To compare SynGO genes with other brain-expressed genes, we deﬁned two control gene sets"

**215.** Here: not used yet.

- `docs/adult_ish_design.md:570-572`: "**SynGO** (Koopmans 2019) for curated synaptic terms"
- `docs/adult_ish_design.md:592`: "SynGO, for later"

**216.** The 390-gene panel comes from Gene Ontology terms, each gene with the terms that placed it.

- `docs/adult_ish_design.md:328-330`: "Gene Ontology terms queried through mygene.info and cached; every gene in `panel_genes.csv` carries the term ids that placed it"

**217.** SynGO is the curated alternative for fixing synaptic gene sets before looking.

- `docs/adult_ish_design.md:570-572`: "with the sets fixed *before* looking"

## Atlases for young brains

**218.** **Carey H, Kleven H, Øvsthus M, et al. (2025).** DeMBA: a developmental atlas for navigating the mouse brain in space and time. *Nat Commun* 16:8108. `Carey_2025.pdf`

- `ref_papers/Carey_2025.pdf` p. 1: "Harry Carey1,3, Heidi Kleven 1,3, Martin Øvsthus"
- `ref_papers/Carey_2025.pdf` p. 1: "Nature Communications| (2025) 16:8108"

**219.** A 4D atlas for every postnatal day from P4 to P56, interpolated between six templates (P4, P7, P14, P21, P28 and P56, the last being the Allen CCFv3) registered with elastix.

- `ref_papers/Carey_2025.pdf` p. 1: "a 4D atlas encompassing every postnatal day from 4 to 56"
- `ref_papers/Carey_2025.pdf` p. 2: "representing ages P4, P7, P14, P21, P28, and P56"
- `ref_papers/Carey_2025.pdf` p. 2: "Using the P56 template (the Allen mouse brain Common Coordinate Framework (Allen CCFv3"
- `ref_papers/Carey_2025.pdf` p. 2: "registration using elastix"

**220.** It carries the CCFv3 and DevCCF labels at every age, and CCF Translator moves coordinates and volumes between ages.

- `ref_papers/Carey_2025.pdf` p. 2: "interpolated the segmentation volumes from the Allen CCFv3 2 and the Developmental Common Coordinate Fra- mework"
- `ref_papers/Carey_2025.pdf` p. 2: "CCF Translator uses deformation matrices to translate coordinates or image v olumes"

**221.** On 54 landmarks its transformations are as accurate as an average expert.

- `ref_papers/Carey_2025.pdf` p. 5: "(n = 54 landmarks)"
- `ref_papers/Carey_2025.pdf` p. 2: "validating the DeMBA transformations matrix to be on average as accurate as a neuroanatomist expert"

**222.** The templates come from different sources, so the apparent shrinking of the brain between P28 and P56 is likely methodological.

- `ref_papers/Carey_2025.pdf` p. 6: "the templates were sourced from multiple independent sources"
- `ref_papers/Carey_2025.pdf` p. 6: "DeMBA appears to show the brain shrinking between P28 and P56, however, it is likely this is due to methodological differences"

**223.** Here: the atlas of the young brains, one per age, and the deformation that carries them to the CCF for the maps.

- `notes/project_atlas_decision_ccf_for_all.md:29`: "One atlas per AGE, not one per cohort"
- `mapping/sepmap/volumes/to_ccf.py:12-14`: "brainglobe_ccf_translator from age_PND to P56"

**224.** Because the labels are the CCF's, region statistics need no warping.

- `mapping/sepmap/young_vs_adult/region_plot.py:5-7`: "the DeMBA annotations were remapped to Allen parcellation_index"
- `notes/project_atlas_decision_ccf_for_all.md:32`: "Region-level comparison still needs no warping"

**225.** Absolute size is not compared across P28 and P56.

- `notes/project_demba_atlas_integration.md:36`: "do not compare absolute size across the P28/P56 boundary"

**226.** **Kronman FN, Liwang JK, Betty R, et al. (2024).** Developmental mouse brain common coordinate framework. *Nat Commun* 15:9072. `Kronman_2024.pdf`

- `ref_papers/Kronman_2024.pdf` p. 1: "Developmental mouse brain common coordinate framework Fae N. Kronman1, Josephine K. Liwang1, Rebecca Betty1"
- `ref_papers/Kronman_2024.pdf` p. 1: "Nature Communications| (2024) 15:9072"

**227.** DevCCF: templates at E11.5, E13.5, E15.5, E18.5, P4, P14 and P56, from MRI and light-sheet data, with developmental labels following the prosomeric model, and the CCFv3 registered into its P56.

- `ref_papers/Kronman_2024.pdf` p. 1: "spanning embryonic day (E)11.5, E13.5, E15.5, E18.5, and postnatal day (P)4, P14, and P56, featuring undistorted morphologically averaged atlas templates created from magnetic resonance imaging and co-registered high-resolution light sheet"
- `ref_papers/Kronman_2024.pdf` p. 2: "developmentally consistent annotations based on the prosomeric model"
- `ref_papers/Kronman_2024.pdf` p. 1: "we map the Allen CCFv3 and spatial transcriptome cell-type data to our stereotaxic P56 atlas"

**228.** Here: it has no template between P14 and P56, so on its own it would have put every young brain on an adult template.

- `ref_papers/PAPER_SUMMARIES.md:73-75`: "so DevCCF alone would still have meant registering P16–P36 brains to an adult template"

**229.** DeMBA's DevCCF labels come from it.

- `ref_papers/Carey_2025.pdf` p. 2: "DevCCFv001) provided by Kronman and colleagues"

**230.** **Chon U, Vanselow DJ, Cheng KC, Kim Y (2019).** Enhanced and unified anatomical labeling for a common mouse brain atlas. *Nat Commun* 10:5067. `Chon_2019.pdf`

- `ref_papers/Chon_2019.pdf` p. 1: "Enhanced and uni ﬁed anatomical labeling for a common mouse brain atlas Uree Chon 1, Daniel J. Vanselow 2, Keith C. Cheng 2 & Yongsoo Kim"
- `ref_papers/Chon_2019.pdf` p. 1: "NATURE COMMUNICATIONS | (2019) 10:5067"

**231.** Franklin-Paxinos labels brought into the CCF, with boundaries adjusted using an MRI atlas and cell-type-specific transgenic mice, and the dorsal striatum segmented by cortico-striatal connectivity.

- `ref_papers/Chon_2019.pdf` p. 1: "we adopt here the FP labels into the CCF"
- `ref_papers/Chon_2019.pdf` p. 1: "We use cell type-speci ﬁc transgenic mice and an MRI atlas to adjust and further segment our labels. Moreover, detailed segmentations are added to the dorsal striatum using cortico-striatal connec"

**232.** Here: not used in the code.

- search of this branch (mapping, group_comparison, preprocessing, registration, atlas, common) for `(?i)paxinos|chon_?2019`: no Franklin-Paxinos labels in the code

**233.** It is the mapping to use when results need Franklin-Paxinos names.

- `ref_papers/PAPER_SUMMARIES.md:93-95`: "if any region definition ever needs reconciling against Franklin-Paxinos conventions"

## Where the code fits

**234.** | raw sections to registered volumes | `preprocessing/run_copy_raw_data.m` to `run_annotate_artifacts.m`, `registration/run_register_to_atlas.m`, `run_add_sep_channel.m`; atlases in `atlas/` | `<group>\<mouse>\lightsuite\volume_registered\`, `volume_registered_sep\` |

- `preprocessing/run_copy_raw_data.m:1-12`: "Preprocessing, step 1 of 6"
- `preprocessing/run_annotate_artifacts.m:1-12`: "Preprocessing, step 6 of 6"
- `registration/run_add_sep_channel.m:25-27`: "volume_registered_sep"
- `group_comparison/run_collect_by_group.m:11-12`: "lightsuite\volume_registered\"

**235.** | where an experience changes the map (line 1) | `group_comparison/run_collect_by_group.m`, `run_normalise_groups.m`, `run_group_differences.m` | `comparisons\<ctrl>_vs_<exp>_<channel>\` |

- `group_comparison/run_collect_by_group.m:5`: "Plasticity comparison, step 1 of 3"
- `group_comparison/run_normalise_groups.m:5`: "Plasticity comparison, step 2 of 3"
- `group_comparison/run_group_differences.m:5`: "Plasticity comparison, step 3 of 3"
- `group_comparison/run_group_differences.m:18`: "data\comparisons\<ctrl>_vs_<exp>_<channel>\"

**236.** | per-brain volumes and the readings (lines 2 to 4) | `mapping/run_per_mouse.py`, `run_to_ccf.py`, `run_cohort.py` (`sepmap/volumes/`) | `comparisons_v2\per_mouse\`, `per_mouse_ccf\`, `ccf\<cohort>\` |

- `mapping/run_per_mouse.py:5-6`: "Writes comparisons_v2/per_mouse/<mouse>.npz"
- `mapping/run_to_ccf.py:4-6`: "comparisons_v2/per_mouse_ccf/<mouse>.npz"
- `mapping/run_cohort.py:3-4`: "writes comparisons_v2/ccf/<cohort>/"

**237.** | the adult distribution (line 2) | `mapping/run_region_plot.py` (the per-mouse region table); `P8_characterize_merged_distribution.m` and `P10_compare_nano_vs_auto.m` until A4 and A5 replace them | `comparisons_v2\young_vs_adult\region_means_per_mouse.csv` |

- `mapping/run_region_plot.py:3-7`: "Writes region_means_per_mouse.csv"
- `REFACTOR_PLAN.md:433` (main): "P8's bar charts, tables, enrichment"
- `REFACTOR_PLAN.md:434` (main): "all of P10, P8's `auto` run"
- `REFACTOR_PLAN.md:392` (main): "once A1 to A5 are in the Python route"

**238.** | more than abundance or density (line 3) | `mapping/run_beyond_density.py`, `run_beyond_controls.py`, `run_beyond_figures.py`, `run_beyond_regression.py` (`sepmap/adult/`) | `adult_v2\beyond\` |

- `mapping/run_beyond_density.py:5-6`: "under adult_v2/beyond/"
- `mapping/run_beyond_controls.py:5`: "under adult_v2/beyond/"
- `mapping/run_beyond_regression.py:3-4`: "regression_table.csv under adult_v2/beyond/for_sami/"

**239.** | what the green channel reports (line 3) | `mapping/run_sep_channel_check.py`, `run_adult_arms.py`, `run_ish_arms.py` | `adult_v2\arms\` |

- `mapping/run_sep_channel_check.py:5-6`: "sep_channel_check.png under adult_v2/arms/"
- `mapping/run_adult_arms.py:4-5`: "under adult_v2/arms/"
- `mapping/run_ish_arms.py:4-5`: "under adult_v2/arms/"

**240.** | the ISH gene comparison (line 3) | 100-gene panel: `mapping/run_ish_regions.py --panel targets`, `run_ish_compare.py`, `run_ish_words.py`, `run_ish_roles.py`; 390-gene panel: `run_panel_build.py`, `run_panel_fetch.py`, `run_ish_regions.py --panel ontology`, `run_ish_reliability.py`, `run_ish_panel_test.py` (`sepmap/ish/`); `P9_compare_nano_vs_allen_ish.m` until A1 to A3 replace it | `adult_v2\ish\`, `adult_v2\panel\` |

- `mapping/settings.toml:23-26`: "run_ish_regions.py and run_ish_reliability.py choose a pass with --panel"
- `mapping/run_ish_regions.py:5-8`: "targets, the hand-written 100-gene panel"
- `mapping/run_panel_build.py:5-6`: "into adult_v2/panel/"
- `mapping/run_ish_panel_test.py:5-6`: "under adult_v2/ish/"
- `REFACTOR_PLAN.md:392` (main): "once A1 to A5 are in the Python route"

**241.** | young against adult (line 4) | `mapping/run_compare.py`, `run_region_plot.py`, `run_region_groups.py`, `run_video.py`, `run_video_compare.py`, `run_closeup.py` (flatmaps, in `tools\venv_flat`), `run_replot.py` (`sepmap/young_vs_adult/`) | `comparisons_v2\young_vs_adult\` |

- `mapping/run_compare.py:4-5`: "writes into comparisons_v2/young_vs_adult/"
- `mapping/run_closeup.py:6`: "so this runs in tools\\venv_flat"
- `mapping/run_replot.py:1`: "Redraw run_compare's slice figures"

**242.** | checking each step by eye | `mapping/run_diagnostics.py` | `comparisons_v2\processing_diagnostics\` |

- `mapping/run_diagnostics.py:1-3`: "Writes comparisons_v2/processing_diagnostics/"

**243.** The run order of the Python route is in the header of every `mapping/run_*.py`.

- `mapping/run_video.py:6`: "Python route, in run order"

## Terms used in the code

**244.** **P-numbers.** In an age, a cohort tag or an atlas key, `P<n>` is postnatal day n (P20, `young_P20`, `demba_p20`).

- `REFACTOR_PLAN.md:70-72` (main): "ages (P20), atlas keys (`demba_p20`), cohort tags (`young_P20`)"

**245.** The old script names P0 to P10 (P4, P6bis, P7bis) were pipeline steps; the refactor renames them `run_...`, and the table of old and new names goes into `docs/ADDING_DATA.md`.

- `REFACTOR_PLAN.md:91-92` (main): "Runnable scripts are named `run_...` with a short header explaining the flow; the P-numbers go"
- `REFACTOR_PLAN.md:292-294` (main): "a table of old and new script names, with the tag `refactor-start`, goes into `docs/ADDING_DATA.md`"

**246.** **Channels.** At acquisition the files carry dye names (`chan02_Cy5` for nano, `chan03_Cy3` for auto); after registration, role names (`chan01_DAPI`, `chan02_NANO`, `chan03_AUTO`, `chan04_DIFF`, `chan05_MASK`).

- `notes/project_preprocessing_recipe.md:58-59`: "`volume_registered\` use **semantic** names, set in P4: `{'DAPI','NANO','AUTO','DIFF','MASK'}` -> `chan01_DAPI`, `chan02_NANO`, `chan03_AUTO`, `chan04_DIFF`, `chan05_MASK`"

**247.** **nano**: the nanobody label of surface SEP-GluA1.

- `registration/run_add_sep_channel.m:13-14`: "the nanobody stain reports the receptors that sat on the membrane"

**248.** The registered NANO is the slice-equalised nano of `run_nano_equalisation`.

- `REFACTOR_PLAN.md:309` (main): "P4 reads its equalized_volume.mat (use_equalized_nano = 1)"
- `notes/project_preprocessing_recipe.md:84`: "nano median equalization across slices"

**249.** **auto**: autofluorescence.

- `preprocessing/run_residual_correction.m:5-7`: "the autofluorescence channel (Cy3)"

**250.** The registered AUTO is fitted onto the nano slice by slice (slope and intercept on reference pixels, `run_residual_correction`).

- `preprocessing/run_residual_correction.m:10-13`: "fits nano against autofluorescence on the reference pixels (robustfit, bisquare): a slope and an intercept"
- `registration/pipeline/register_to_atlas.m:486`: "cat(4, dapiVol, nanoVol, scaledautoVol, correctedVol, invalid_mask)"

**251.** **DIFF**: nano minus the scaled autofluorescence, negatives set to zero.

- `preprocessing/pipeline/residual_correction.m:453-454`: "corrected = J - scaled_I; corrected(corrected < 0) = 0;"

**252.** It is registered, but no analysis here reads it.

- search of this branch (mapping, group_comparison) for `chan04_DIFF|DIFF\.tif`: no analysis reads the registered DIFF channel

**253.** **MASK**: the invalid pixels, annotated artefacts or background.

- `registration/pipeline/register_to_atlas.m:475`: "invalid_mask = single(or(artifact_mask_vol,bg_mask_vol));"

**254.** **DAPI**: nuclei; the channel that drives registration (`regchan = 'dapi'`).

- `preprocessing/make_ordering_volume.m:46`: "% nuclei"
- `notes/project_preprocessing_recipe.md:52`: "**regchan** | 'dapi' | **'dapi'**"

**255.** **SEP**: the green channel (`chan04_EGFP`, named after the filter), carried into registered space by `run_add_sep_channel` as `volume_registered_sep\chan02_SEP.tiff`.

- `registration/run_add_sep_channel.m:21-23`: "chan04_EGFP.tiff"
- `registration/run_add_sep_channel.m:22-23`: "EGFP names the filter, SEP the molecule"
- `registration/run_add_sep_channel.m:25-27`: "volume_registered_sep"
- `registration/run_add_sep_channel.m:27`: "chan02_SEP.tiff what the analysis reads"

**256.** Mostly autofluorescence in this tissue (line 3).

- `notes/project_sep_channel_not_total_receptor.md:3`: "the green/SEP channel in fixed cleared tissue is mostly autofluorescence"

**257.** `channel = 'nano'` or `'auto'` in the MATLAB drivers picks the channel a run analyses; it goes into the output folder names.

- `group_comparison/run_normalise_groups.m:44-45`: "Channel to normalize: 'nano' (default, surface GluA1) or 'auto' (autofluorescence control)"
- `group_comparison/run_group_differences.m:18`: "<ctrl>_vs_<exp>_<channel>"

**258.** **Readings** (Python route, defined in `sepmap/volumes/cohort.py`), all of the background-subtracted signal:

- `mapping/sepmap/volumes/cohort.py:10`: "Five readings of the same background-subtracted signal"

**259.** `ratio`: nano per unit autofluorescence, voxel by voxel.

- `mapping/sepmap/volumes/cohort.py:15`: "ratio sig / auto per voxel"

**260.** Not absolute: autofluorescence rises with age, so it understates a young deficit.

- `mapping/sepmap/volumes/cohort.py:21-22`: "autofluorescence rises with age (lipofuscin, tissue density), so this reading understates a real pup deficit"

**261.** `sepratio`: nano per unit SEP.

- `mapping/sepmap/volumes/cohort.py:24-25`: "sig / SEP per voxel"

**262.** Built as a surface fraction; it is not one.

- `mapping/sepmap/volumes/cohort.py:25-35`: "It does not mean what it was built to mean"
- `mapping/sepmap/volumes/cohort.py:33-35`: "but it is not a surface fraction"

**263.** `cref`: nano over that brain's isocortex mean.

- `mapping/sepmap/volumes/cohort.py:36`: "cref sig / (that mouse's isocortex mean of sig"

**264.** A pure scale, so region ratios within a brain survive exactly; blind to a change of the whole cortex.

- `mapping/sepmap/volumes/cohort.py:37-39`: "a pure scale, so region ratios within a mouse survive exactly. Cortex-relative by construction, hence blind to a change that moves the whole cortex"

**265.** `subref`: nano over that brain's subcortex mean, without hippocampus and striatum.

- `mapping/sepmap/volumes/cohort.py:40`: "subref sig / (that mouse's subcortex mean, excluding HPF and STR)"

**266.** `zref`: below.

- editorial: pointer within this document

**267.** `sepauto` (SEP per unit autofluorescence) appears only in the channel arms of line 3.

- `docs/adult_ish_design.md:202`: "`sepauto` = SEP / auto | total receptor"
- `mapping/run_adult_arms.py:3`: "SEP/auto, nano/auto and nano/SEP per structure for the ten adults"

**268.** **zref.** zref = [log2(nano / cortex mean) - m] / s, with m the median and s the p90 - p10 spread of log2(structure mean / cortex mean) across that brain's structures.

- `notes/project_young_vs_adult_comparison.md:173-174`: "zref = [ log2(nano) - m ] / s, m = median and s = p90-p10 of that brain's distribution ACROSS STRUCTURES of log2(structure mean)"
- `notes/project_zref_difference_is_already_relative.md:18-20`: "zref = `[ log2(nano / that mouse's isocortex mean) - m ] / s`"

**269.** Zero is the brain's median structure, mostly subcortical, not a physical null.

- `notes/project_young_vs_adult_comparison.md:176-177`: "So zero is the brain's MEDIAN STRUCTURE (mostly subcortical), not the cortex and not a physical null"

**270.** One unit is that brain's own spread.

- `notes/project_young_vs_adult_comparison.md:176-177`: "one unit is that brain's own spread"

**271.** A zref difference is already a difference of logs and cannot be converted to a fold change: +1 is about 2x in a young brain and about 3.8x in an adult.

- `notes/project_zref_difference_is_already_relative.md:21-22`: "It is a LOG quantity, so a difference of two zref maps is already a log ratio"
- `notes/project_young_vs_adult_comparison.md:180-182`: "can never be converted to a fold change"
- `notes/project_zref_difference_is_already_relative.md:31-32`: "about 2x"
- `notes/project_zref_difference_is_already_relative.md:31-32`: "about 3.8x"

**272.** Quote numbers from the tables and use the maps for the pattern (over 238 areas the two agree at r = 0.985, but the maps read lower).

- `notes/project_young_vs_adult_comparison.md:184-185`: "Over 238 areas r = 0.985, same sign and order always, but the map reads LOWER"
- `notes/project_young_vs_adult_comparison.md:191`: "Rule: quote numbers from the tables, use the maps for pattern"

**273.** **LR-sum and LR-diff** (MATLAB route, `common/compute_lr_stats.m`): left hemisphere plus, or minus, the mirrored right hemisphere.

- `common/compute_lr_stats.m:20-21`: "lr_diff = left - flipped_right; lr_sum = left + flipped_right;"

**274.** The plasticity comparison uses L + R and |L - R| per mouse.

- `group_comparison/run_group_differences.m:15`: "for every mouse |L - R| and L + R"

**275.** **Surprise.** -log10 p of the group-difference test, drawn as a map.

- `group_comparison/run_group_differences.m:16-17`: "Welch t and surprise (-log10 p) maps"

**276.** **Cohorts.** The MATLAB groups and data folders are `naive`, `rws`, `behavior` and `young` (registry: `common/get_cohort.m`).

- `common/get_cohort.m:11`: "group 'rws' | 'naive' | 'behavior' | 'young'"
- `common/get_cohort.m:102`: "fullfile(base_root, reg{i, 2}, reg{i, 1})"

**277.** The Python cohorts are `young` (P16, P20 and P22 pooled), `young_P20`, `young_P16`, `young_P22`, `naive`, `rws` and `adult` (naive plus RWS).

- `mapping/sepmap/volumes/cohort.py:56-58`: "Cohorts: young (P16 + P20 + P22 pooled"
- `mapping/sepmap/volumes/cohort.py:57-58`: "and adult = naive + rws"

**278.** Adult folders end in `_Gria1`, young ones in `_SepGluA_P<age>`.

- `common/get_cohort.m:58`: "'MG691_Gria1'"
- `common/get_cohort.m:88`: "'MG911_SepGluA_P16'"

**279.** The adults' age is not recorded in the registry.

- `common/get_cohort.m:58`: "'MG691_Gria1', 'rws', NaN, ''"
- `notes/project_developmental_dataset.md:90`: "what age were the adult Gria1 mice?"

**280.** **CCF.** The Allen Mouse Brain Common Coordinate Framework v3, the P56 reference, a population average of 1,675 mice.

- `ref_papers/Carey_2025.pdf` p. 2: "Using the P56 template (the Allen mouse brain Common Coordinate Framework (Allen CCFv3"
- `ref_papers/Carey_2025.pdf` p. 8: "based on 1675 subjects"

**281.** Adults are registered at 10 µm, in the crop of CCF planes 180 to 1079.

- `notes/project_preprocessing_recipe.md:50-51`: "| atlasaplims | [180 1079] | [180 1079] |"
- `mapping/sepmap/volumes/to_ccf.py:16`: "already registered to the CCF crop [180 1079] at 10 um"

**282.** The Python route works on the CCF at 20 µm (660 x 400 x 570).

- `mapping/sepmap/volumes/cohort.py:4`: "660 x 400 x 570 at 20 um"

**283.** **DeMBA.** The developmental atlas of Carey et al. 2025.

- `ref_papers/Carey_2025.pdf` p. 1: "DeMBA: a developmental atlas"

**284.** `atlas/build_demba_atlas.py <age>` builds one per age from BrainGlobe's `demba_allen_seg_dev_mouse_p<age>_20um`, with the labels remapped to the CCF's parcellation index.

- `notes/project_atlas_decision_ccf_for_all.md:52-56`: "`build_demba_atlas.py <age>` builds"
- `notes/project_atlas_decision_ccf_for_all.md:54`: "fetch `demba_allen_seg_dev_mouse_p<age>_20um`"
- `notes/project_atlas_decision_ccf_for_all.md:55`: "remap structure ids -> parcellation_index"

**285.** The files keep LightSuite's `_10` names, so a young brain needs `px_atlas = 20`.

- `REFACTOR_PLAN.md:745-747` (main): "the 20 um DeMBA volumes are renamed `*_10` and need `px_atlas = 20`"
- `atlas/build_demba_atlas.py:9-10`: "so the files must carry the _10 name"

**286.** CCF Translator (`brainglobe-ccf-translator`) carries a young brain to the CCF.

- `notes/project_atlas_decision_ccf_for_all.md:33`: "`brainglobe-ccf-translator` 0.21 is installed"
- `mapping/sepmap/volumes/to_ccf.py:12-14`: "brainglobe_ccf_translator from age_PND to P56"

**287.** **RWS.** Rhythmic whisker stimulation, as in Gambino et al. 2014.

- `REFACTOR_PLAN.md:34` (main): "rhythmic whisker stimulation, RWS, as in Gambino et al. 2014"
- `ref_papers/Gambino_2014.pdf` p. 1: "rhythmic whisker stimulation (RWS) protocol"

**288.** **Structure sets.** The beyond-abundance analysis uses grey-matter structures only and no "..., unassigned" labels (125 structures).

- `notes/project_map_beyond_abundance_density.md:39-41`: "grey matter divisions only, no "..., unassigned""
- `docs/adult_ish_design.md:444`: "**125 structures** survive the rule"

**289.** The plan's declared set (S1) keeps grey-matter structures seen in all ten adults, 234 of 280.

- `REFACTOR_PLAN.md:201-205` (main): "seen in all 10 adults"
- `REFACTOR_PLAN.md:204-205` (main): "The all-10 rule keeps 234 of 280 structures"

**290.** **Reliability and Spearman-Brown.** Agreement between two independent halves of a cohort; Spearman-Brown extends it to the full cohort.

- `docs/adult_ish_design.md:449-451`: "the half-cohort maps agree at **ρ = 0.974**; Spearman-Brown gives **0.987** for the full cohort"

**291.** **Cross-validated R².** Variance explained on structures held out of the fit.

- `docs/adult_ish_design.md:456-457`: "Scored by cross-validation on held-out structures"

**292.** **Expression energy.** The Allen ISH value per voxel of its 200 µm grid.

- `docs/adult_ish_design.md:593`: "200 µm grid, expression *energy*"

**293.** **A1 to A10, S1 to S6.** The scientific additions and decisions of [REFACTOR_PLAN.md](REFACTOR_PLAN.md).

- `REFACTOR_PLAN.md:428-439` (main): "A1 | the declared structure set"
- `REFACTOR_PLAN.md:157-160` (main): "Decisions on the science (30 Sep, Giulio)"

## Open points for the owner

Found while checking; none is stated as fact in SCIENTIFIC_CONTEXT.md.

1. **The RWS and behaviour protocols used here.** The sources name RWS "as in
   Gambino et al. 2014" and "behavior", nothing more: not the sessions given
   to our RWS mice, nor the task of the behaviour group, nor the delay before
   perfusion.
2. **The naive-against-behaviour result.** The sources record only that the
   rerun reproduces it (0.988 to 0.999), not what it shows.
3. **The adults' age.** Not in the registry (`NaN` in `common/get_cohort.m`),
   and open since August (`notes/project_developmental_dataset.md:90`). The
   Allen ISH mice are 56 days old.
4. **The grant's Figure 3b caption** gives the pups as "aged between P16 and
   P20" (n = 7), while the seven young brains analysed include MG904 at P22
   (`notes/project_young_vs_adult_comparison.md:102-104`). Which brains made
   the figure?
5. **Permeabilised or not.** `registration/run_add_sep_channel.m:11-12` says
   the tissue is "fixed, permeabilised"; the grant describes non-detergent
   labelling. The answer matters for what the SEP channel and the nanobody can
   reach.
6. **The reading of the by-system numbers** (somatosensory +0.21 to
   hippocampus -1.03) is taken as `cref` from
   `notes/project_young_vs_adult_comparison.md:88` (n = 6, "By system (cref)")
   and :126-128 (unchanged at n = 7). To confirm in `group_stats.csv`.
7. **Gria1's rank.** S5 holds it back from any claim; the document gives the
   two ranks (10th, 17th) only to show that it moves with the structure set.
   Drop the numbers if even that is too much.
8. **Sources outside the repository.** Most numbers of lines 2 to 4 rest on
   the working notes. Before this checklist is deleted, they should point to
   a file in the repository or under `data\` (the planned `docs/ROADMAP.md`
   and `docs/FIGURES.md`, or the output tables).
9. **The plan on this branch is older** than `main`'s; the link to
   REFACTOR_PLAN.md in SCIENTIFIC_CONTEXT.md is right once the two meet.
