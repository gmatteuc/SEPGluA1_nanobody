# Coverage of the retiring MATLAB analyses by the Python route

Evidence behind the section "Python route: what it must gain before P8, P9 and P10 retire" of [REFACTOR_PLAN.md](REFACTOR_PLAN.md). Every scientific output, test and option of P8, P9 and P10 was listed from the code and matched to the Python code that answers the same question (in spirit, not necessarily the same output). A second, independent check then read both sides of every item again, tried to refute each match, and looked for anything missed. Read-only throughout: nothing was run except small scripts reading existing tables (30 Sep 2026, code at 04c0484).

Status: **covered** (the Python route answers the question at least as well), **partial** (part of the meaning is lost, said how), **missing**, **dropped** (dropped by decision on 29 Sep: P8's per-mouse videos, the voxel-level ISH correlation, P8's distance-weight diagnostic figure), **file handling** (caches, file I/O). The status shown is the one after the second check; where the two disagreed, both are given. Items numbered X were found only by the second check.

| script | covered | partial | missing | dropped | file handling | total |
|---|---|---|---|---|---|---|
| P8 | 14 | 12 | 5 | 2 | 2 | 35 |
| P9 | 11 | 10 | 4 | 1 | 7 | 33 |
| P10 | 3 | 4 | 10 | 1 | 2 | 20 |
| all | 28 | 26 | 19 | 4 | 11 | 88 |

## Specifications of the additions

Numbered in the order they run; where two proposals were merged into one addition, both are given. These are starting specifications, refined when the addition is built. Numbers marked read-only were recomputed from existing tables without writing anything.

Where a specification differs from the plan, the plan applies: the minimum number of mice is the one decided in S1 (recommended: all 10 adults; the adult table has 280 structures, 22 of them seen in fewer than 5 adults, 15 of those in the gene ranking); the eroded nano mean is required (A3, A4); modules live under `mapping/sepmap/`; the `MIN_N` key error of `v2_video` is fixed in the Fixes step, not in A10.

### A1. Declared structure set and zref reference (required)

**One declared zref reference set, shared by the tables and the maps**

- **Replaces:** P8's z-score over planes 100-700 (P8-11, P8-13, and the second check's missed item on AP extremes in the zref normalisation).
- **Question:** The zero of zref ('the brain's median structure') is what every 'above/below' statement rests on. It should depend only on that brain and a fixed structure list, not on which young brains are in the run or on noisy hindbrain sections.
- **Module:** mapping/volumes/cohort.py and mapping/young_vs_adult/region_plot.py, both reading the set from mapping/structures.py
- **Specification:** Today the table route takes each brain's median and spread over structures shared by all 17 brains (v2_region_plot.py:245-251, 183 structures), while the map route uses every structure the mouse has, extremes included (v2_cohort.py:162-171). Define the reference set once, e.g. keep_structure ∩ present in all adults, written to a CSV, and use it in both. Report the change per mouse: median shift 0.02-0.06 and zero-point shift up to 0.11 in the adult table (read-only check); the maps are expected to move more. The ISH Spearman results are unchanged (affine per mouse); the A4 enrichment calls and the young-adult zref difference are what move.

### A2. ISH per-section quality control (required)

**ISH per-section quality control: flag failed sections, set them to NaN, and report the effect on the ranking**

- **Replaces:** P9's bad-section detection and repair (P9-07), the NativeProfile diagnostic figure (P9-08) and the n_bad_sections column (P9-22).
- **Question:** Are the Allen ISH structure means, and therefore the gene ranking, corrupted by failed or dim sections in single-experiment genes? This is Sami's 28 Apr 'verify ISH preprocessing quality for all genes'.
- **Module:** mapping/ish/section_qc.py (called from mapping/ish/regions.py before region_means)
- **Specification:** Per gene, per 200 um AP plane of the grid, inside the annotation_200 brain mask: median energy and valid fraction. Flag a plane that has data but whose median is below f x the median of its ±3 in-brain neighbours (proposed f = 0.2). Flag, never interpolate: flagged planes become NaN, consistent with read_energy's -1 handling (v2_ish_regions.py:107). A hand-editable exceptions CSV keeps genuine regional absence (e.g. anterior Tac1, posterior Glra1). A read-only scan with this rule already flags 14 of 95 genes, several of them machinery or controls that bear on the headline: Nptx1 (planes 21, 26, 37), Cnih3 (9, 47, 61-63), Cbln2 (21), Cacng3 (55), Cacng4 (5), Nrgn (44, 53), Map2 (12), Kcna1 (11), Hcn2 (47), Grm2 (65), Grik4 (8), Slc32a1 (5), Tac1 (4-5), Glra1 (15-16, 61). Gria1 is clean. Outputs: ish/section_qc.csv (gene, experiment, plane, median, local median, flag, reason); ish/section_qc.png (profiles of flagged genes with flagged planes marked); gene_region_table.csv rebuilt with the mask; and the ranking before vs after (Spearman of the two gene rankings, Gria1 and Cacng8 ranks). Apply it to the 390-gene panel table too (V2_ISH_PANEL route).

### A3. Gene ranking on the declared set, and its robustness (required)

**A declared structure set and minimum-mice rule for the gene ranking, then rerun and re-quote the ISH headline**

- **Replaces:** P9's region set and sampling (P9-04) and the headline gene ranking and dissociation claim (P9-23). It also answers P8-13 (AP extremes) for the ISH side.
- **Question:** Which genes best predict the adult nano map, and does Gria1 sit below AMPAR machinery genes? Answered on structures whose nano value is actually measured in the cohort, rather than on hindbrain structures seen in one to four adults.
- **Module:** mapping/structures.py (shared: keep_structure moved from beyond_density, plus adult_profile(reading, min_mice, rule)), used by mapping/ish/compare.py, roles.py, panel_test.py, arms.py and words.py
- **Specification:** Finding: adult_profile in v2_ish_compare.py:75-83, v2_ish_panel_test.py:75-81 and v2_ish_roles.py:143-149 averages whatever adults have a structure, with no minimum. Of the 234 structures shared with Cacng8's ISH, 15 are measured in fewer than 5 adults (8 pons, 3 medulla, 2 MB, 1 CB, 1 OLF). P/MY/CB structures are seen in a median of 3 adults, and their profile does not replicate across mice: leave-one-out rho is about 0.14, against 0.93 for forebrain+MB (computed read-only from region_means_per_mouse.csv). Recomputed read-only with zref: all structures give Gria1 rank 10/95 (rho 0.665, 4/33 machinery above). At least 5 adults gives rank 17 (0.636, 7/33). At least 8 or all 10 adults gives rank 17 (0.618, 7/33). P9's nine divisions give rank 18 (0.579, 7/33). Cacng8 stays first throughout (0.77-0.82). v2_beyond_density.prepare is already protected (all adults plus the grey rule, v2_beyond_density.py:355-358). Implement a single adult_profile with MIN_MICE (proposed: all 10, or 8) and keep_structure. Rerun compare, roles (the 3,876-split permutation), panel_test (p for localisation vs matched controls, currently 0.74), arms and words. Write ish/region_set_sensitivity.csv with the headline numbers (Gria1 rank and rho, Cacng8-Gria1 gap, machinery above Gria1, old-vs-new agreement, panel_test p, roles p) under three sets: v2 as is, the declared default, and P9-equivalent (the 9 divisions, no '-unassigned'). Then correct the quoted numbers in docs/adult_ish_design.md and the project notes; 'Gria1 climbs to rank 10' is mostly a few-mice hindbrain effect, not a zref effect.

**Robustness of the gene ranking to the statistic and to border handling**

- **Replaces:** The 5-metric bar chart (P9-25), the eroded-vs-distance-weighted comparison (P9-24) and the border-contamination guard (P9-11, P8-14 on the ISH side), and Pearson (P9-13).
- **Question:** Is the gene ranking an artefact of the choice of statistic or of neighbours' signal leaking into small structures? The old data show the statistic mattered: Spearman vs Pearson gene rankings agreed at only 0.81, with Gria1 at rank 22 vs 9.
- **Module:** mapping/ish/robustness.py
- **Specification:** On the A3 declared structure set, per gene: Spearman on full ISH means (primary), Spearman on ish_mean_eroded (already written at v2_ish_regions.py:194 and never read), and Pearson on log2 values of both sides. Optional: a nano-side eroded mean (one 20 um voxel peeled per structure in v2_region_plot's bincount step) as a fourth column. Output: ish/ranking_robustness.csv with one row per gene and each variant's rho and rank, plus the pairwise agreement between rankings and the Gria1/Cacng8 ranks under each. One small figure: rho under each variant against the primary. Read-only indication: full vs eroded ISH agree at 0.998 (Gria1 rank 10 -> 7).

### A4. Adult regional distribution (required)

**Adult regional distribution: every structure by division, per-mouse mean +/- SEM, per-structure reliability and an enrichment call (with a cohort argument)**

- **Replaces:** P8 Region_MeanSum/Zscore_BarByMacro (P8-18, P8-19), the _withSEM BarByMacro and BarAcrossDivi figures (P8-23, P8-24), Region_MeanSum_Table (P8-15), per-region reliability shading (P8-17), the enrichment classification with orange labels (P8-21), the division bars (P8-20), the map-vs-table check (P8-25), and the young_P20 characterisation run (P8-05, merged_young_P20_nano).
- **Question:** How is surface GluA1 distributed across all brain structures in the adult brain? Which structures sit reliably above or below the brain's own typical level, and how reproducible is each structure's value across mice? This is Sami's 30 Apr ask ('how reliable is the spatial pattern across mice') and his 28 Apr enrichment call. The same question is asked of the young cohort on its own.
- **Module:** mapping/adult/distribution.py (called with cohort='adult' by default; 'young' and 'young_P20' supported)
- **Specification:** Inputs: data/comparisons_v2/young_vs_adult/region_means_per_mouse.csv (readings zref and cref, from v2_region_plot.py:305-308). Recompute zref from cref, not from the stored column: the stored adult zref takes its median and spread over the 183 structures shared by all 17 brains, young included (v2_region_plot.py:245-251). Its zero moves by 0.03-0.11 per adult if the set is taken from the 234 structures shared by the adults alone, and it will move again when the next young brains are added. Structure set: v2_beyond_density.keep_structure (grey-matter divisions, no '-unassigned'), plus a declared MIN_MICE (proposed 8 of 10 adults), plus a column listing every structure dropped and why. Per structure: n_mice, mean, SD, SEM, t = mean/SEM, and the naive-minus-RWS difference. In cref, also the CV of the linear value across mice. The enrichment call: an exact two-sided Wilcoxon signed-rank test of zref against 0 across mice, with BH q across the structure family, classifying each structure as above / below / not different from the brain's median structure. Report the direction, and note that the minimum exact p at n=10 is 1/512. Division summary: per-mouse voxel-weighted division means (as in v2_region_groups.py:185-190), with Isocortex kept whole plus the 8 other P8 divisions, then mean +/- SEM. Outputs: adult_v2/distribution/structure_distribution_<cohort>.csv, ranked by mean; division_distribution_<cohort>.csv; one figure per reading (zref, cref) with one panel per division, horizontal bars (mean +/- SEM), a dot per mouse (naive and RWS in two greys), bar grey level set by reliability |t| (house rule: reliability = grayscale), a zero line, and the BH-significant labels highlighted; a division-level figure of the same form. Printed check (P8-25): the median relative difference between the mean of per-mouse structure means and the structure mean of the v2_cohort adult map. For comparison, place the figures next to P8's Region_*_withSEM PNGs.

### A5. Autofluorescence as a parallel control, per structure (required)

**Autofluorescence as a parallel control, per structure and per division: nano vs auto with a centred contrast, paired signed-rank tests, BH, and the auto distribution itself**

- **Replaces:** P10 as a whole: the paired and delta-z signed-rank tests and their tables (P10-02..P10-06), the division pooling, tests and table (P10-07..P10-10), the BarByMacro, DeltaZ and Macro figures (P10-11..P10-14). Also the descriptive autofluorescence characterisation of P8's channel='auto' run (P8-06: auto bar charts and region table).
- **Question:** Structure by structure, is the nano pattern its own signal, distinct from the autofluorescence of the same brains? Which structures are higher in nano than in autofluorescence, relative to each brain's own level, consistently across mice? What does the autofluorescence regional pattern look like when treated exactly like nano (Sami 30 Apr: 'all channels must be treated identically')?
- **Module:** mapping/adult/nano_vs_auto.py (reusing mouse_channels from mapping/adult/sep_channel_check.py)
- **Specification:** Inputs: per-mouse log2 structure means of nano (sig) and auto for the 10 adults, from v2_sep_channel_check.mouse_channels (v2_sep_channel_check.py:75-94). Both channels are background-subtracted identically (v2_per_mouse.py:177-190), within the same tissue mask. Do not use the stored 'ratio' reading: it is a mean of voxelwise ratios clipped at RATIO_CLIP, and its level depends on exposure (v2_region_plot.py:196-208). Structure set: the same declared set as A4. Two contrasts per mouse and structure. (a) Centred log-ratio: [log2 nano_s - log2 auto_s] minus its median over the reference structures, the exposure-free analogue of P10's paired test. (b) delta-zref: zref_nano - zref_auto, each channel with the zref formula (own isocortex mean, then minus median and divided by the p90-p10 spread over the same reference set), the analogue of P10's delta-z. Statistics: an exact two-sided Wilcoxon signed-rank test across mice per structure, reported with mean, SEM and t, and BH q within each contrast family. Also report P10's one-sided nano>auto p, so the old lists can be compared directly. Guard (the second check found 6 of the 21 old delta-z hits were low-auto periventricular structures with nano_z < 0): write nano_z and auto_z beside the contrast, plus a class column (nano above median / below). Division level: per-mouse voxel-weighted means for the 9 P10 divisions with Isocortex whole, the same tests, and Bonferroni over 9 (feasible at division level: 0.0056 > 1/1024). Outputs: adult_v2/nano_vs_auto/structures.csv and divisions.csv. Figures: (1) per-division panels with nano zref and auto zref per structure side by side, mean +/- SEM, per-mouse dots joined by a thin line per mouse (this also serves as the auto distribution figure that replaces P8-06); (2) the contrast bars sorted by contrast, with the zero line and significant labels highlighted; (3) the division-level paired figure and its contrast figure. Checks: print the overlap with the old significant lists (P10 CSVs: 28 uncorrected, 23 BH paired; 21 and 11 delta-z). Differences are expected, because P10's auto came from the stale 19 Nov auto_4d.mat.

### A6. Within-division gene agreement, per-gene scatters (recommended)

**Within-division gene agreement, and per-gene structure scatters**

- **Replaces:** The per-gene region scatters (P9-19), the within-division regressions (P9-20), the per-gene paired bars (P9-21) and the per-gene Delta-z spatial view in its tabular form (P9-17).
- **Question:** Does a gene track the nano map inside each division, or only through between-division differences (cortex vs thalamus)? Which structures agree or disagree for a specific gene such as Gria1 or Cacng8?
- **Module:** mapping/ish/gene_detail.py
- **Specification:** For every gene, on the A3 structure set, Spearman within each division that has at least 8 structures, and the mean within-division rho, weighted by number of structures. Output: ish/within_division.csv. Compare the ranking by mean within-division rho with the brain-wide ranking, reporting Gria1's and Cacng8's ranks. It is also a cheap partial answer to the missing spatial null, since it removes the largest gradient. Figure for a list of genes (default Gria1, Cacng8, Grm5, Dlg2, one control): a scatter of ISH rank vs nano zref per structure coloured by division, with acronyms on outliers, within-division fits and rho in the legend. Optionally, a paired bar strip sorted by nano showing the per-structure rank difference.

### A7. Spatial null (recommended)

**Spatial null for the nano-gene correlations and for regional enrichment**

- **Replaces:** Nothing in either route: Sami's 28 Apr permutation request, never built in P8/P9. v2_ish_compare.py:34-38 lists it as the next piece of work.
- **Question:** Does a gene correlate with the nano map more than a random map with the same spatial smoothness would? Does a structure's enrichment exceed what random spatially autocorrelated maps give?
- **Module:** mapping/ish/spatial_null.py (the enrichment side in mapping/adult/distribution.py)
- **Specification:** Structure-level surrogates that preserve spatial autocorrelation, using structure centroids in CCF (from annotation_10): variogram-matched surrogates (BrainSMASH-style, Burt 2020) of the nano zref structure profile, 10,000 per run, plus within-division label permutation as a simpler alternative. Per gene: p_spatial for rho(nano, gene), BH across genes. Report how many of the 95 genes remain significant, and whether Cacng8 vs Gria1 differ beyond the null (a surrogate distribution of rho_Cacng8 - rho_Gria1). For enrichment: the surrogate distribution of each structure's zref, as a p for the A4 call. Output: ish/spatial_null.csv and a figure of the null band behind the ranked rho values.

### A8. Gene panel against the autofluorescence map (recommended)

**The gene panel ranked against the autofluorescence map on its own**

- **Replaces:** P9's channel='auto' control (P9-03; the flag existed but was never run) and Sami's 30 Apr request to run the whole chain, P9 included, on autofluorescence.
- **Question:** Are the nano-gene correlations specific to the nanobody, or would the tissue's autofluorescence produce the same gene ranking?
- **Module:** mapping/ish/compare.py (an extra 'auto' profile beside the nano readings)
- **Specification:** Auto zref per structure (from A5's inputs, adult mean over the A3 structure set), correlated with all 95 genes exactly as nano is. Output: an extra row set in gene_correlations.csv (reading='auto'), the distribution of rho_auto against rho_nano per gene, and the ranks of Gria1, Cacng8 and the machinery under auto. Expected result, from rho_auto_gria 0.04-0.38 in sep_channel_check.csv: auto tracks genes weakly, which makes the specificity claim explicit per gene.

### A9. Gene documentation table (recommended)

**Gene documentation table for the supplement**

- **Replaces:** Not a P9 output, but Sami's 28 Apr request tied to the P9 panel (symbol, full name, function, why included, reference).
- **Question:** Methods and supplement: what each gene in the 95-gene and 390-gene panels is, why it is there, and which Allen experiment(s) were used.
- **Module:** mapping/ish/gene_table.py
- **Specification:** Join gene_targets.csv (description = why included), ish/gene_roles.csv (curated role), panel/panel_genes.csv (GO terms, n_experiments, planes), ish/gene_reliability.csv, ish/annotation/<symbol>.json (mygene 'name'; only 20/95 records have 'summary'), the experiment ids and gene_drops.csv reasons. Output: ish/gene_documentation.csv (and .xlsx for the supplement). A 'reference' column is left for hand curation.

### A10. Display extras (videos) (optional)

**Display extras: auto and nano-minus-auto cohort videos, and a threshold contour on the zref video**

- **Replaces:** P8 threshold videos (P8-30, P8-31), P10 contrast video (P10-15), P8's auto cohort video (P8-06), and the per-gene Delta-z video (P9-17).
- **Question:** Where in the brain, voxel by voxel, is nano enriched relative to its own brain and relative to autofluorescence?
- **Module:** mapping/young_vs_adult/video.py (add 'autozref' and 'deltaz' readings behind a flag) or mapping/adult/videos.py
- **Specification:** Signed maps on the diverging scale already used for zref (SIGNED_READINGS, v2_cohort.py:85). Auto zref per voxel uses the same per-brain scalars as A1. deltaz = zref_nano - zref_auto per mouse, then the cohort mean with a t panel. Optional contour at the A4 call. Also fix MIN_N so it has a 'young_P22' key (v2_video.py:45; `v2_video.py young_P22` currently raises KeyError).

### Points still to settle

- P8-05 (first check 'covered', second check 'partial'): resolved as PARTIAL. MIN_N lacks 'young_P22' (v2_video.py:45; the default run at l.146 draws only young and young_P20), and no young-only distribution figure exists; young is only ever plotted against adults (v2_region_plot.py:336-389, v2_region_groups.py:305-309). Settled by giving A4 a cohort argument; no separate addition needed.
- P9-23 (first check 'covered', second check 'partial'): resolved as PARTIAL, with a sharper cause than the second check's. I reproduced its numbers (all structures: Gria1 rank 10; nine divisions: rank 18). The shift comes mainly from 15 structures measured in fewer than 5 adults (mostly pons and medulla, leave-one-out rho ~0.14), not from white matter and not from zref: requiring at least 5 adults alone gives rank 17 and 7/33 machinery above Gria1. The quoted v2 headline ('Gria1 climbs to rank 10, gap halves'; the September notes on the Python route and on measurement validation) is provisional until A3 is run. Whether panel_test's p=0.74 and the roles permutation change is unknown until rerun.
- P9-24 (first check 'dropped_by_design', second check 'partial'): the owner's list says 'distance-weight diagnostic' (29 Sep). P8 has a figure of exactly that name (Diagnostic_DistWeight, P8:1397-1491). P9's correlation_eroded_vs_distweight (P9:735-749) is a ranking-robustness test, so it is not clearly covered by the drop. The owner should confirm. Either way, A3 answers the underlying question cheaply, and the figure itself need not be rebuilt.
- P9-17: the first check says the per-gene Delta-z video is listed as optional in the plan. It was not: draft 1 listed only per-gene scatters; only the working notes of 29 Sep listed the threshold / delta-z video (now decision S4). The owner should decide optional vs dropped and record it in the plan (A10 carries it as optional).
- Which structure set is the declared default for the ISH ranking and the adult distribution (A3/A4): all 10 adults vs at least 8, and whether hindbrain and cerebellum stay in. This is an owner decision; I propose keep_structure plus all (or 8 of) 10 adults as default, with P9's nine-division set reported alongside.
- The definition of 'enriched' (A4): P8 used raw LR-sum >= 1.5 or voxel-mean z >= 0; Sami asked (28 Apr) for a threshold of 2 on the LR-sum and, in the same feedback, for a permutation null. The proposal is a per-structure signed-rank test of zref against the median structure with BH, with the A7 spatial null as the stronger version. This needs Giulio's (and possibly Sami's) agreement, since zero now means the median structure, not the voxel mean.
- P9-07 details: the second check's specific failed sections (Eno2 s43, Grik3 s39, Lrfn2 s42) came from the old MATLAB section indexing. My local-neighbour scan of the v2 grids did not flag them, which suggests they are -1 in the grid and already NaN in v2. It did flag Nrgn 44 and 13 other genes. Which genes are affected is only settled once A2 runs with a reviewed exceptions list.
- P10 direction: P10 tested one-sided nano > auto on a mid-brain-referenced scale dominated by HPF. A5 proposes two-sided tests on centred contrasts, which will not reproduce the old lists one to one. The owner should accept that the question is kept, not the output.

## P8: P8_characterize_merged_distribution.m (adult distribution)

### P8-01 · covered

- **What:** What the input means: P5 -> P6bis normalised nano. Each mouse gets a robust affine fit (slope and intercept) of its isocortex L1-L5 pixels, with RSP/ACA/PL left out, against its group's median consensus, and is mapped as (raw - intercept)/slope. The background mask comes from nano intensity and raw zeros are set to NaN. Every brain ends up on its group's cortical intensity scale before pooling.
- **Kind:** normalisation (input)
- **Old code:** P8_characterize_merged_distribution.m:141-177 (loads <channel>_4d_normalized<tag>.mat + bkgmask); P6bis_analyze_group_averages_and_normalize.m:134-140 (cortex mask), 182 (select_background_pixels on nano), 242-303 (consensus + robust fitlm), 783-813 (apply, NaN raw zeros)
- **Python:** v2_per_mouse.py:main (l.141; off-tissue background subtraction, tissue mask from auto) + v2_cohort.py:mouse_scalars (l.122-173) / mouse_modes (l.176-210; cref = sig / own isocortex mean, zref)
- **Judgement:** Same aim (put every mouse on a comparable scale before pooling), met in a better way. Python subtracts the off-tissue background explicitly and then scales each mouse by a pure scale (its own isocortex mean), with no affine fit to a group consensus, so region ratios survive where the intercept broke them. The tissue mask comes from auto, not nano, and RSP/ACA/PL are not excluded from the reference. An earlier check of the Python route found that zref reproduces the old gene ranking at rho 0.91, but single numbers move (Gria1 goes from rank 22 to rank 10).
- **Second check:** Old: P6bis fits each mouse robustly against the median of its group's pooled isocortex L1-L5 pixels, with RSP/ACA/PL excluded (P6bis:134-140, 242-303, fitlm RobustOpts). It applies (raw-intercept)/slope and NaNs raw zeros (783-813). New: v2_per_mouse.py:176-190 subtracts the off-tissue median, and the tissue mask is auto > bg + 4 MAD plus nano coverage (l.189). v2_cohort.py:201-204 divides by the mouse's own isocortex mean (cref) or builds zref. The aim, putting mice on a common scale before pooling, is met with a pure scale, so region ratios survive. The cortex reference now includes all layers and RSP/ACA/PL (v2_cohort.py:157-159), which moves only where 1 sits. Rank 22 -> 10 and rho 0.91 match the earlier check of the Python route (September 2026).

### P8-02 · covered

- **What:** Pooling naive + rws into one adult cohort, behavior excluded. The non-reference group is linearly aligned to the reference group by a polyfit of its mean per-slice background-subtracted tissue profile over adult planes 200-700.
- **Kind:** normalisation (cohort pooling)
- **Old code:** P8:25, 179-197 (per-slice profiles), 227-249 (polyfit alignment), 261-266 (merge)
- **Python:** v2_cohort.py COHORTS['adult'] = NAIVE + RWS (l.109-119); no group alignment needed because of per-mouse scaling in mouse_modes; pooling checked in v2_region_plot.py:main (naive_minus_rws_log2 column, l.295/314) and v2_beyond_controls.py:control_d_groups (l.185-200)
- **Judgement:** Same 10 adults and still no behavior group. Python needs no group-level affine alignment because each mouse is scaled on its own. It also tests explicitly whether pooling hides a group difference (control D, naive vs RWS leftover rho +0.88, plus a per-structure naive-rws null), which P8 never did.
- **Second check:** Same 10 adults: I read per_mouse_mean_dw's mouse list from the P8 cache (CGF027/028/033/034/035 + MG691/692/693/736/737), which matches v2_cohort.py:109-118. Behavior is absent from every Python cohort. Per-mouse scaling makes P8's polyfit group alignment (P8:227-249) unnecessary. The pooling check exists: region_stats naive_minus_rws_log2 (v2_region_plot.py:295, 314) and control D (v2_beyond_controls.py:185-200). adult_v2/beyond/controls.csv gives rho +0.884, as the first check claims.

### P8-03 · covered

- **What:** Common unit: the merged cohort is divided by the mean tissue intensity over adult planes 300-500 (common_fact). 'Raw' LR-sum is therefore in units of mid-brain tissue mean per hemisphere, and raw_threshold = 1.5 is expressed in those units.
- **Kind:** normalisation (unit)
- **Old code:** P8:251-263
- **Python:** v2_cohort.py:mouse_scalars (l.157-160) and mouse_modes (cref l.201, subref l.202); same scalars in v2_region_plot.py:main refs (l.239-240)
- **Judgement:** Both routes express intensity against an internal brain reference. P8 used one cohort factor, taken after the affine step from all tissue in a mid-AP slab. Python uses a per-mouse scalar: the isocortex mean (cref) or the subcortex without HPF/STR (subref). The values are not interchangeable (cortex = 1 in cref, while mid-brain tissue = 1 per hemisphere in P8), so the old raw threshold has no direct cref equivalent.
- **Second check:** P8:253-263 takes one cohort factor, the mean tissue intensity over adult planes 300-500. Python uses per-mouse internal references (v2_cohort.py:157-160 cortex/subcortex means; v2_region_plot.py:239-240). The purpose, intensity relative to the brain itself, is the same. As the first check says, the units differ, so raw_threshold=1.5 has no direct equivalent. That matters only for P8-21.

### P8-04 · covered

- **What:** Optional NaN-robust 3D Gaussian smoothing of each mouse (sigma 5 voxels) before any analysis. Off by default (apply_smoothing=false).
- **Kind:** processing option
- **Old code:** P8:32-34, 111-116, 199-225
- **Python:** v2_inspect.py:prepare (l.166-185, SMOOTH=(3,1,1) l.100), display only; the statistics in v2_region_plot/v2_cohort are never smoothed
- **Judgement:** P8's production outputs were unsmoothed (the _nosmooth suffix), and so are the Python statistics. Python offers light, mask-normalised smoothing for display only. What is gone is only the option to smooth the statistics, which was never used.
- **Second check:** apply_smoothing=false (P8:33), and the production outputs carry the _nosmooth suffix (data/comparisons/merged_naive_rws_nano/*_nosmooth*). Python smooths for display only (v2_inspect.py:166-183, SMOOTH l.100). The only statistics-adjacent smoothing is the denominator of ratio/sepratio (v2_cohort.py:196-199) and the displayed log2 difference map in v2_compare.py:197, neither of which is the nano statistic. What is lost is an option that was never used.

### P8-05 · partial (first pass: covered)

- **What:** Characterising a young cohort on its own atlas: SEP_MERGE_SPECS=young_P20, DeMBA P20, AP indices scaled by ap_scale. This produced data/comparisons/merged_young_P20_nano.
- **Kind:** option (cohort)
- **Old code:** P8:21-30, 118-139
- **Python:** v2_cohort.py COHORTS young / young_P20 / young_P16 / young_P22 (l.111-119); v2_video.py:main (videos per young cohort); v2_region_plot.py:main (per-mouse structure means on each brain's own age atlas)
- **Judgement:** Python characterises every young age, not only P20. Maps are in CCF after a per-mouse DeMBA->CCF warp (v2_to_ccf); region numbers are computed on each brain's own atlas without warping. P8 did one age at a time on DeMBA P20.
- **Second check:** Refuted as 'covered'. P8 run with SEP_MERGE_SPECS=young_P20 produced a full young-only characterisation on DeMBA P20: Region_* bar charts by division, raw and z, the _withSEM versions, threshold videos and a diagnostic at plane 685 (data/comparisons/merged_young_P20_nano). Python has the young per-structure numbers (region_means_per_mouse.csv, v2_region_plot.py:305-308) and CCF-warped videos. But young brains are only ever plotted against adults: region_plot.png shows 35 areas (v2_region_plot.py:84-86, 336-389) and group_plot.png shows systems/divisions (v2_region_groups.py:305-309). There is no young-only by-division distribution figure, and addition A4 was first scoped to adults (the plan now gives it a cohort argument). 'Every young age' is overstated for maps: v2_video.py:146 draws only young and young_P20 by default, and MIN_N (l.45) has no 'young_P22' key, so `v2_video.py young_P22` raises KeyError. The young-vs-adult use of this run (compare_young_vs_adult_lrsum.py:49-50) is fully covered by v2. To settle: give addition A4 a cohort argument (now in the plan).

### P8-06 · partial

- **What:** Autofluorescence control run (channel='auto', Sami 30 Apr): the whole characterisation (cohort video, region tables, bar charts, SEM plots, threshold videos) repeated on the auto channel, to show the nano pattern is not autofluorescence or tissue.
- **Kind:** control analysis
- **Old code:** P8:57-61, 107, 149-163 (outputs in data/comparisons/merged_naive_rws_auto)
- **Python:** v2_sep_channel_check.py:main (l.112ff; per-mouse structure means of nano, auto and SEP, rho_nano_auto, dynamic ranges); v2_beyond_density.py:autofluorescence (l.222-249) + step2_covariates (auto as covariate); v2_diagnostics.py:sheet_scaling / sheet_denominators; ratio reading in v2_cohort.mouse_modes
- **Judgement:** Python answers the question more sharply: nano vs auto rho is about 0.26 per mouse, and auto explains about 0% of the zref map ranks. Lost: the descriptive autofluorescence map itself (no auto cohort video, no auto per-structure table or figure by division). The per-structure nano-vs-auto test is the planned P10 addition. P8's auto run read auto_4d.mat, which P5 no longer refreshes (P5:160-169), so the old auto outputs come from an earlier registration.
- **Second check:** Numbers verified. adult_v2/arms/sep_channel_check.csv gives mean rho_nano_auto 0.257 +- 0.11, range nano 1.93 vs auto 1.07. adult_v2/beyond/variance_partition.csv puts 'autofluorescence only' at CV R2 -0.04 and in-sample 0.004. The auto per-structure means are computed only inside v2_beyond_density.autofluorescence (l.222-249) and are never written to a file. v2_sep_channel_check writes per-mouse summaries only (l.130-149). v2_adult_arms writes sepauto/ratio/sepratio, not plain auto. So there is no auto map, auto table or auto figure. P5:155-169 confirms auto_4d.mat is no longer written, so P8's auto outputs are stale.

### P8-07 · covered

- **What:** Hemispheric fold: per-mouse LR sum (left + mirrored right, compute_lr_stats), then abs(), then the cohort mean over mice.
- **Kind:** analysis step
- **Old code:** P8:270-275; compute_lr_stats.m:1-27
- **Python:** v2_video.py:fold (l.50-52) / v2_compare.py:fold (l.78); v2_region_plot.py:main pools both hemispheres (bilateral labels)
- **Judgement:** Python averages the two hemispheres (a mean, not a sum), which also survives one missing hemisphere, and folds the cohort volumes rather than each mouse. The abs() that turned unreached planes into bright slabs is gone. The old adult P8 outputs still contain those slabs, because P6bis was never re-run after the raw-zero fix (comparisons_v2/README.md 'Still open').
- **Second check:** compute_lr_stats.m:21 is left + flipped right (a sum, NaN if either hemisphere is missing), and P8:275 applies abs(). Python folds with nanmean (v2_video.py:50-52, v2_compare.py:78-82) on cohort-mean volumes, not per mouse. Region means pool bilateral labels (v2_region_plot.py:210-223). The slab artefact is still in the old outputs (comparisons_v2/README.md, 'Still open', l.211-215).

### P8-08 · covered

- **What:** Cohort mean map of the pooled adults: mean LR-sum per voxel.
- **Kind:** map
- **Old code:** P8:275 (saved at 323-340)
- **Python:** v2_cohort.py:main (l.213-236) -> comparisons_v2/ccf/adult/<reading>_mean.npy (plus _sd, _n)
- **Judgement:** Per-voxel mean over the mice that have tissue at that voxel, with an n map. It is on the 20 um CCF grid rather than the 10 um registered crop, in five readings (cref is closest to P8's raw map, zref to its z map).
- **Second check:** v2_cohort.py:213-236 writes <reading>_mean/_sd/_n.npy for cohort 'adult' = NAIVE+RWS, averaging per voxel over mice with tissue (n counts tissue, l.224). It is on the 20 um CCF grid. This is the same question as P8:275.

### P8-09 · covered

- **What:** Voxelwise reliability across mice: t = mean / SEM.
- **Kind:** statistic (voxel)
- **Old code:** P8:277-281
- **Python:** v2_video.py:main l.88-94 (t from v2_cohort mean, sd and n)
- **Judgement:** Same statistic. Python counts only mice with tissue at the voxel (P8 counted every non-NaN value, including slab artefacts), and its SD is computed per hemisphere and then averaged, not taken over per-mouse L+R sums.
- **Second check:** v2_video.py:88-94: t = mean/(sd/sqrt(n_h)), with sd = nanmean of the per-hemisphere SDs and n_h = max(nL, nR) (fold_n l.55-57). This is the same statistic. Two caveats: when the hemispheres rest on different mouse subsets, the n used overstates the effective n; and t is only drawn in the video, never saved or summarised.

### P8-10 · covered

- **What:** Cohort tissue mask: a voxel is analysed where any mouse has non-background tissue (nano-derived background masks, hemispheres OR'd).
- **Kind:** mask
- **Old code:** P8:283-294
- **Python:** v2_per_mouse.py:main (auto > off-tissue median + 4 MAD, registered nano != 0, inside the atlas) + v2_cohort n maps + MIN_N in v2_video.py (l.45) / MIN_N_YOUNG/ADULT in v2_compare.py (l.42)
- **Judgement:** Better than P8. The mask no longer comes from nano intensity (P8's route flagged partial sections as background). Each mean carries its n, and voxels below MIN_N are shown grey instead of accepting a single mouse.
- **Second check:** v2_per_mouse.py:189 defines tissue = brain & nz>=0.5 & auto > bg_a + 4*MAD, so nano intensity is not used. The n maps (v2_cohort.py:224, 231) and MIN_N (v2_video.py:45, v2_compare.py:42) replace P8's 'any mouse' OR (P8:293). Checked.

### P8-11 · covered

- **What:** Within-brain z-score of the cohort mean map: each voxel in SD units above the brain-wide voxel mean, with mean and SD taken over brain voxels in adult planes 100-700.
- **Kind:** statistic (normalisation)
- **Old code:** P8:296-313; zscore_in_mask l.1940-1948
- **Python:** v2_cohort.py:mouse_modes zref (l.203-204), with z_median/z_spread from mouse_scalars (l.167-171); comparisons_v2/ccf/adult/zref_mean.npy
- **Judgement:** Same idea: where a voxel sits relative to the brain's own typical level, in units of its own spread. zref differs in four ways: it is log2, robust (median and p90-p10), computed over structures, and computed per mouse before averaging. P8's z-score is linear mean/SD over the voxels of the cohort-mean volume. Zero is the median structure, not the voxel mean, so 'z > 0' selects different voxels.
- **Second check:** zref (v2_cohort.py:203-204, scalars l.161-171) answers the same question: position relative to the brain's own typical level, in units of its own spread. The first check correctly lists the differences: log2, robust, structure-based, per mouse, and zero at the median structure rather than the voxel mean. One more: v2_cohort's median/spread use every structure the mouse has with at least 250 voxels, fibre tracts and AP extremes included (l.167-168), while P8 restricted z to planes 100-700 (P8:301-305).

### P8-12 · covered

- **What:** Region set: every STRU leaf under 9 DIVI macros (Isocortex, OLF, HPF, CTXsp, STR, PAL, TH, HY, MB) from the Allen ontology, about 228 leaves, shared by tables, bars and video labels.
- **Kind:** region set
- **Old code:** P8:80-92, 342-367; build_stru_leaves_by_macro l.1783-1908
- **Python:** v2_region_plot.py:main (all 'structure'-level terms keyed by name, MIN_VOX=250, l.83, 176-224); v2_beyond_density.py:keep_structure (l.167, 177-188)
- **Judgement:** Same STRU level, but over all divisions: the per-mouse table adds pons, medulla, cerebellum, fibre tracts, ventricles and '-unassigned' labels (about 250 structures per adult). The adult 'beyond' analyses keep 125 grey-matter structures and drop catch-all labels by rule. There is no 9-macro restriction.
- **Second check:** Checked on the outputs. All 228 P8 leaves (Region_MeanSum_Table csv) except MEV, which had no P8 value, appear in the adult cref rows of region_means_per_mouse.csv (280 structures, 234 in all 10 adults). Within the 9 macros, the only Python extras are the '-unassigned' labels. Caveat for addition A4: MIN_VOX=250 per mouse (v2_region_plot.py:83, 223) leaves 25 of the shared leaves in fewer than 10 adults (CUN 2, VTN 2, IV 3, AOB 4, AT 4, DR 4, Pa4 4, VISpl 5, ...). P8 gave every leaf a value from the cohort-mean volume.

### P8-13 · partial

- **What:** AP restriction of every quantitative analysis (z-score, region means, per-mouse means) to adult planes 100-700, to exclude noisy anterior and posterior extremes. The videos still show the full range.
- **Kind:** analysis restriction
- **Old code:** P8:99-102, 139, 301-303, 555-558, 1059-1060
- **Python:** No explicit equivalent. Closest: v2_per_mouse.py:main tissue mask, v2_region_plot.py MIN_VOX (l.83), v2_diagnostics.py:sheet_coverage (l.128)
- **Judgement:** Python guards against planes a section never reached, but not against poor-quality extreme sections. It measures olfactory bulb, cerebellum and brainstem structures that P8 cut off, and whether those extremes are noisy is never re-examined.
- **Second check:** Python has no AP restriction anywhere (grep). Coverage is handled per voxel (v2_per_mouse.py:189) and MIN_VOX. One point strengthens 'partial': the AP extremes now also feed the zref normalisation. v2_cohort.mouse_scalars takes median/spread over all structures the brain has (l.162-169), OB, CB and medulla included, so a poor extreme section shifts every structure's zref in that mouse's maps. region_plot is somewhat protected because it uses only structures shared by all brains (l.245-251).

### P8-14 · partial

- **What:** Border-robust region mean: distance-weighted mean (bwdist to the power 4, primary) and eroded mean (sphere r=3), both saved, to guard against hippocampal bleed-through and contamination from neighbouring regions.
- **Kind:** statistic (region aggregation)
- **Old code:** P8:66-71, 94-97, 546-602, 632-688
- **Python:** v2_region_plot.py:main (plain voxel mean over tissue voxels of each structure at 20 um, l.210-223); an eroded/full comparison exists only on the ISH side (v2_ish_regions.py:region_means)
- **Judgement:** The region mean exists, but the nano side has no border guard and no robustness comparison (weighted vs eroded vs plain). The 29 Sep decision dropped the distance-weight diagnostic figure, not the guard itself. Small structures' means are partly their neighbours' under registration error.
- **Second check:** The nano region means are plain bincount means (v2_region_plot.py:210-223), and no v2 file computes an eroded or distance-weighted nano mean. Eroded means exist only for ISH (v2_ish_regions.py:131-172). The only indirect guard is control B, structure size vs residual, rho -0.116 (v2_beyond_controls.py:149-162), which tests the leftover, not border contamination. The 20 um 2x2x2 block mean (v2_per_mouse.py:136-172) adds some edge mixing. An eroded-mean column is cheap to add alongside addition A4.

### P8-15 · partial

- **What:** Region_MeanSum_Table (.csv and .mat): the cohort mean per STRU leaf (weighted and eroded), voxel counts (raw and eroded), macro, sorted descending. It is the ranking of regions by mean nano.
- **Kind:** table
- **Old code:** P8:604-615 (cache reload 525-545)
- **Python:** v2_region_plot.py:main -> comparisons_v2/young_vs_adult/region_means_per_mouse.csv (l.305-308, per mouse, all readings, n_vox20); region_stats.csv adult_mean_log2 (l.309-315); adult_v2/beyond/for_sami/regression_table.csv observed_rank (v2_beyond_regression.py:main)
- **Judgement:** The numbers exist, and per mouse is better. What is missing is an adult-only, all-structure summary table ranked by mean. region_stats.csv drops structures that fewer than 2 young brains cover (l.287-288). region_table.csv (v2_compare) restricts adults to voxels the young also cover. regression_table.csv has only 125 grey-matter structures, as zref ranks.
- **Second check:** All cited facts check out. region_stats needs at least 2 young and 4 adults (v2_region_plot.py:287-288). region_table uses the 'both' mask (v2_compare.py:183, 235). regression_table.csv covers 125 structures as ranks (v2_beyond_regression.py:211-217). residual_by_structure.csv carries nano_rank for the same 125 (v2_beyond_density.py:551-557). No adult-only, all-structure table ranked by mean exists.

### P8-16 · covered

- **What:** Z-scored region means: region mean minus the brain-wide voxel mean, divided by the brain voxel SD.
- **Kind:** statistic
- **Old code:** P8:618-624
- **Python:** zref per mouse per structure in region_means_per_mouse.csv (v2_region_plot.py:main value() l.256-271)
- **Judgement:** Same question: where a region sits in its brain's own distribution. zref is per mouse, log2 and robust, with its median and spread taken over structures shared by all brains (l.245-251), so the numbers differ from P8's linear, voxel-based z.
- **Second check:** zref per mouse per structure is written for all readings (v2_region_plot.py:256-271, 305-308). Its median/spread come from the structures shared by all 17 brains (l.245-251). It answers the same 'where does a region sit in its brain' question as P8:618-624.

### P8-17 · partial

- **What:** Per-region reliability: the distance-weighted mean of the voxel t-score inside each region, voxel-weighted per macro. It sets the grey level of every bar.
- **Kind:** statistic (reliability)
- **Old code:** P8:626-747, 842-866
- **Python:** No per-structure equivalent. Map-level reliability: v2_beyond_density.py:step1_ceiling (l.438-463; 126 half-splits)
- **Judgement:** Python measures how reproducible the whole regional pattern is across animals (half-cohort rho 0.974, Spearman-Brown 0.987), which is stronger at the map level. It does not give a per-structure reliability or SEM, so which structures are unreliable is never reported. addition A4 would cover the per-structure variability.
- **Second check:** No per-structure SD, SEM or t in any v2 file: the std/sem grep hits are only gene z-scoring, residual spread and ISH reliability. Map-level reliability is v2_beyond_density.step1_ceiling (l.438-463), restricted to the 125 grey-matter structures that also have ISH (prepare/main l.352-361, 614-617). So 'which structures are unreliable' is never reported.

### P8-18 · missing

- **What:** Region_MeanSum_BarByMacro figure: one panel per DIVI macro; bars are STRU leaves sorted by mean, each bar shaded by its reliability t, with a dashed raw threshold line and enriched labels in orange.
- **Kind:** figure
- **Old code:** P8:749-832
- **Python:** None. The nearest is v2_region_plot.py:main region_plot.png (l.336-389), but it shows only 36 hand-picked AREAS (l.84-86) as a young-vs-adult dot plot.
- **Judgement:** There is no adult all-structure distribution figure by division. It is addition A4 of the plan, and the data are already in region_means_per_mouse.csv.
- **Second check:** I searched all v2_*.py for barh/bar/errorbar: none draws an all-structure adult distribution by division. Correction to the first check: region_plot has 35 areas, not 36 (AREAS l.84-86 includes '|'). The nearest adult per-structure display is v2_beyond_regression.panel_f (F_maps: observed zref rank of 125 structures painted on 3 coronal planes, l.136-182), also without divisions, values or reliability. The v2_inspect flatmaps show the isocortex pattern spatially only.

### P8-19 · missing

- **What:** Region_Zscore_BarByMacro figure: the same per-division panels on z-scored region means, with a threshold at z = 0 and enriched labels.
- **Kind:** figure
- **Old code:** P8:925-1006
- **Python:** None (zref values in region_means_per_mouse.csv)
- **Judgement:** Same gap as P8-18. The natural Python version is addition A4 drawn in zref, where zero is the median structure.
- **Second check:** No z-scored by-division figure exists. The zref values are in region_means_per_mouse.csv, but nothing draws them per division.

### P8-20 · partial

- **What:** Region_MeanSum / Zscore_BarAcrossDivi figures: one bar per macro (weighted by eroded voxel count over its leaves), raw and z-scored, shaded by reliability, with threshold highlighting.
- **Kind:** figure
- **Old code:** P8:834-923
- **Python:** v2_region_groups.py:main -> group_plot.png (dotplot l.272-303; SYSTEMS l.61-77, DIVISIONS l.86-87, voxel-weighted per brain l.186-190) and group_stats.csv adult_median (l.242-256)
- **Judgement:** Division-level summaries exist, with per-mouse dots, voxel-weighted per brain, in every reading. They are framed as young vs adult with medians, with no adult ranking, no reliability shading and no threshold. Isocortex is split into 11 systems, so there is no whole-isocortex bar (in cref that bar would be 0 by construction).
- **Second check:** v2_region_groups.py builds voxel-weighted per-brain group means (l.186-190). DIVISIONS (l.86-87) lack Isocortex as a whole; SYSTEMS (l.61-77) split it. The dotplot (l.272-303) shows naive, rws and young medians only, and group_stats keeps a group only with at least 3 young and 5 adults (l.240). There is no adult ranking, reliability shading or threshold.

### P8-21 · missing

- **What:** Enrichment classification (Sami, 28 Apr): a region or voxel is called enriched if its mean is >= the raw threshold (1.5 LR-sum units) or its z is >= 0. Enriched labels are orange on every bar chart, the percentage of voxels above the z threshold is printed, and labels turn cyan in the videos.
- **Kind:** statistic (classification)
- **Old code:** P8:42-48, 310-313, 800-808, 904-912, 973-980, 1210-1216, 1321-1326, 1350-1364; highlight_enriched_labels l.1911-1937
- **Python:** None (a grep of v2_*.py finds no threshold or enrichment classification of the nano map)
- **Judgement:** Python never calls a region enriched or not. The sign of zref (above or below the median structure) is the obvious analogue, but it is never turned into a call or a list, and a threshold would have to be defined in zref or cref units. The permutation test Sami preferred to a fixed threshold was never built in either route.
- **Second check:** A grep of v2_*.py for enrich/threshold finds only ISH word enrichment (v2_ish_words.py:44, 274, 303) and mask thresholds. The only permutation test is the leftover-noise null (v2_beyond_figures.py:116-131), not a regional enrichment test. Sami's area-preserving permutation test (Sami's feedback of 28 Apr) exists in neither route.

### P8-22 · covered

- **What:** Per-mouse per-region means (distance-weighted and eroded), cached as a matrix.
- **Kind:** table
- **Old code:** P8:1016-1096
- **Python:** v2_region_plot.py:main -> region_means_per_mouse.csv (l.274-308)
- **Judgement:** Per mouse and per structure, in all five readings, computed on each brain's own atlas at 20 um. They are plain voxel means, not weighted or eroded (see P8-14).
- **Second check:** region_means_per_mouse.csv (v2_region_plot.py:274-308) holds per mouse, per structure, per reading values on each brain's own atlas. These are plain means over the mouse's own tissue voxels. P8 used distance-weighted and eroded means masked by the cohort brainMask_merged and planes 100-700 (P8:1057-1060).

### P8-23 · missing

- **What:** _withSEM BarByMacro figures (raw and z): per-region mean across mice +/- SEM, per-mouse dots, reliability shading and threshold (Sami, 30 Apr: 'how reliable is the spatial pattern across mice').
- **Kind:** figure
- **Old code:** P8:1105-1116, 1132-1239
- **Python:** None
- **Judgement:** This is exactly addition A4 (every structure by division, one dot per mouse, mean and SEM). The per-mouse data are already in region_means_per_mouse.csv.
- **Second check:** No per-structure mean +/- SEM figure with per-mouse dots exists in any v2 file. It matches addition A4.

### P8-24 · partial

- **What:** _withSEM BarAcrossDivi figures (raw and z): per-mouse macro means (weighted by eroded voxel count over leaves), mean +/- SEM across mice, per-mouse dots.
- **Kind:** figure
- **Old code:** P8:1241-1340
- **Python:** v2_region_groups.py:main -> group_plot.png (per-mouse dots and medians by system and division)
- **Judgement:** Per-mouse dots per division exist, with naive and rws drawn separately. There is no pooled-adult mean +/- SEM bar, and the figure is framed as young vs adult.
- **Second check:** group_plot.png draws per-mouse dots by division/system, naive and rws separately, with the median as the bar (v2_region_groups.py:275-286). There is no pooled-adult mean +/- SEM and no whole-isocortex entry.

### P8-25 · partial

- **What:** Consistency check: the mean of per-mouse region means against the region mean of the cohort-mean volume (median and max relative difference, warning above 10%).
- **Kind:** diagnostic
- **Old code:** P8:1118-1130
- **Python:** v2_diagnostics.py:sheet_route_agreement (07, l.234ff); comparisons_v2/README.md 'Maps and tables do not give the same number' (r 0.985 over 238 areas)
- **Judgement:** Python checks map-vs-table agreement, but for the young-adult difference (cref) and for zref. The adult-level agreement of per-mouse means with the cohort map per structure is not reported.
- **Second check:** sheet_route_agreement (v2_diagnostics.py:234-257) compares only the cref young/adult difference, voxelwise vs region-wise. The zref version is the README note ('Maps and tables do not give the same number', comparisons_v2/README.md:168-180, r 0.985 over 238 areas). v2_adult_arms self-checks against region_plot, which is a different comparison. The adult-level check of per-mouse mean vs cohort-map region mean is absent.

### P8-26 · covered

- **What:** Across-mouse reproducibility of the adult regional pattern: the question the SEM plots and reliability shading were built to answer.
- **Kind:** statistic (reproducibility)
- **Old code:** P8:1105-1116 (pm_sem), 626-747 (reliability colour)
- **Python:** v2_beyond_density.py:step1_ceiling (l.438-463); v2_beyond_figures.py:panel_b (l.162); v2_beyond_controls.py:control_c_mice (l.165-182) and control_d_groups (l.185-200)
- **Judgement:** Answered at the map level and more rigorously: split-half agreement over all 126 splits (rho 0.974, whole cohort 0.987) on 125 grey-matter structures, plus pairwise animal and naive-vs-RWS checks. The last two are run on the leftover, not the raw map. The per-structure part is P8-17 and P8-23.
- **Second check:** The map-level reproducibility asked for on 30 Apr ('how reliable is the spatial pattern across mice') is answered by step1_ceiling over 126 splits (v2_beyond_density.py:438-463) with a bootstrap CI (v2_beyond_figures.py:238-252). Caveat: it is computed only on the 125 structures shared with the ISH panel and on zref ranks. Controls C/D (v2_beyond_controls.py:165-200) run on the leftover (C pairwise median +0.78, D +0.884 in controls.csv).

### P8-27 · covered

- **What:** Merged cohort LR-sum video: mean (clim 0-5) beside reliability t (clim 0-10), with atlas outlines and smoothed region acronyms.
- **Kind:** video
- **Old code:** P8:434-473; write_lr_video_labeled l.1562-1695; write_lr_video.m
- **Python:** v2_video.py:main (l.72-142) -> comparisons_v2/ccf/adult/video_<reading>_adult.mp4 (also naive, rws and the young cohorts)
- **Judgement:** Same layout, as the docstring says, now for five readings. The t panel is scaled to its 95th percentile and voxels with n < MIN_N are grey. There is no autofluorescence-channel counterpart (see P8-06).
- **Second check:** v2_video.py:72-142 draws mean and t panels with outlines and acronyms (docstring l.2-5: P8 layout), with the t panel scaled to the 95th percentile (l.94) and MIN_N grey. P8's panels were 0-5 and 0-10 (P8:452, 1642). There is no auto counterpart (P8-06).

### P8-28 · file handling

- **What:** Region-label centroid table for the video overlays: per-slice centroids, smoothed by a moving average, drawn only above 150 px.
- **Kind:** display helper
- **Old code:** P8:73-78, 369-432
- **Python:** v2_video.py:main l.122-129 (per-plane center_of_mass, MIN_LABEL_AREA=150)
- **Judgement:** Labels for display only; Python recomputes them per plane.
- **Second check:** Display only. P8:369-432 vs v2_video.py:122-129 (center_of_mass, MIN_LABEL_AREA=150, l.47).

### P8-29 · dropped

- **What:** Per-mouse LR-sum videos (off by default).
- **Kind:** video
- **Old code:** P8:38, 475-511; write_single_mouse_sum_video(_labeled) l.1495-1559, 1698-1780
- **Python:** None
- **Judgement:** Dropped on 29 Sep (plan, Decisions already taken). Per-mouse volumes remain inspectable in comparisons_v2/per_mouse*.
- **Second check:** The owner's decision is recorded in the plan ('Dropped on 29 Sep, not rebuilt in Python: P8's per-mouse videos, ...'). The option was also off by default (P8:38). Per-mouse volumes exist in data/comparisons_v2/per_mouse and per_mouse_ccf.

### P8-30 · missing

- **What:** Raw-intensity threshold video: mean LR-sum with a p5-p90 clim, a cyan contour around voxels at or above the raw threshold, and enriched region labels in cyan.
- **Kind:** video
- **Old code:** P8:1342-1383; write_volume_threshold_video l.1951-2057
- **Python:** None (v2_video cref and ratio videos have no threshold contour)
- **Judgement:** It depends on the enrichment call (P8-21), which Python does not make. Project notes list a threshold video as optional.
- **Second check:** No v2 video draws a threshold contour (v2_video.py:106-140 has none, and grep finds no threshold).

### P8-31 · partial

- **What:** Z-scored threshold video: the z map with a contour at z >= 0 and enriched labels.
- **Kind:** video
- **Old code:** P8:1371-1374, 1385-1392
- **Python:** v2_video.py:main zref video (PuOr diverging scale centred on 0, SIGNED_READINGS in v2_cohort.py l.85)
- **Judgement:** The zref video shows above and below the brain's median structure directly through its diverging scale, which gives the same picture. There is no explicit contour or enriched-label highlight, and zero is the median structure, not the voxel mean.
- **Second check:** The zref video uses PuOr_r on -2..+2 for signed readings (v2_video.py:81, 95-96; SIGNED_READINGS v2_cohort.py:85), so above and below the median structure is readable. There is no contour and no enriched-label highlighting.

### P8-32 · dropped

- **What:** Distance-weight diagnostic at adult plane 620 (hippocampal bleed-through), weighting powers 1, 2 and 4: raw, weights, weighted, suppressed.
- **Kind:** diagnostic figure
- **Old code:** P8:1397-1491
- **Python:** None
- **Judgement:** Dropped on 29 Sep. Note that the border guard it justified is also absent from the Python nano means (P8-14).
- **Second check:** Listed as dropped in the plan (Decisions already taken). The note that the border guard itself is also absent is correct (see P8-14).

### P8-33 · file handling

- **What:** Caches for downstream scripts: mean_lr_sum/t/z volumes + masks (for P9), sparse ROI masks in P9's format, region table .mat, label centroids, per-mouse matrix.
- **Kind:** cache / I/O
- **Old code:** P8:315-340, 375-431, 525-545, 632-688, 1019-1040
- **Python:** v2_cohort.py:main (.npy volumes), v2_cohort.py:mouse_scalars (*_scalars.npz cache)
- **Judgement:** Plumbing. The only readers are P9, P10 and compare_young_vs_adult_lrsum.py, all of which are retiring; P7bis reads nothing from P8.
- **Second check:** The readers of P8 files are P9, P10 and compare_young_vs_adult_lrsum.py:49-50 (grep for mean_lr_sum_, Region_MeanSum_Table, roi_masks_r). P7bis reads none (only its own Region_Stats/Surprise names, P7bis:970-1697). One addition: v2_ish_compare.py:61-62, 112-118 optionally reads P9's gene_panel_summary.csv, which descends from P8's caches. It is guarded by os.path.exists, so retiring the code breaks nothing, but deleting that data folder would silently drop the old-vs-new ranking comparison.

### P8-X1 · partial

- **What:** Young-cohort-only characterisation outputs from the SEP_MERGE_SPECS=young_P20 run: region table, per-division bar charts raw and z, the _withSEM versions and threshold videos for the young brains on their own atlas. The first check files this under P8-05 as 'covered' via the cohort option, but these outputs have no Python counterpart.
- **Old code:** P8:21-30 (SEP_MERGE_SPECS), 118-139 (DeMBA P20, ap_scale); outputs in data/comparisons/merged_young_P20_nano/Region_*_merged_young_P20_nano_nosmooth*
- **Python:** v2_region_plot.py:305-308 (young rows in region_means_per_mouse.csv); v2_video.py:146 (young and young_P20 videos only, in CCF after warp); young is only ever plotted against adults (v2_region_plot.py:336-389, v2_region_groups.py:305-315)
- **Judgement:** The young per-structure numbers exist, but no figure shows the young regional distribution on its own. Addition A4 takes a cohort argument in the plan, so the same figure can be drawn for 'young'. Separately, v2_video.MIN_N (l.45) lacks 'young_P22'.

### P8-X2 · partial

- **What:** Effect of AP-extreme structures on the zref normalisation itself. P8 computed its z reference (mean, SD) only over planes 100-700. v2_cohort.mouse_scalars takes each brain's zref median and p90-p10 spread over every structure that brain has with at least 250 voxels, OB, CB, medulla and fibre tracts included, so the extremes P8 excluded now shift the zref maps of every structure in that mouse.
- **Old code:** P8:99-102, 301-305 (brainMask_analysis restricts zscore_in_mask to planes 100-700)
- **Python:** v2_cohort.py:162-171 (per-mouse median/spread over all its structures); v2_region_plot.py:245-251 uses only the structures shared by all brains, which is safer
- **Judgement:** This belongs with P8-13 but was not noted there. The table route is mostly protected by the shared-structure rule; the map route is not. A check of zref median/spread with and without OB/CB/P/MY would settle whether it matters.

## P9: P9_compare_nano_vs_allen_ish.m (nano against Allen ISH)

### P9-01 · covered

- **What:** Gene panel: which genes the adult nano map is compared against (gene_targets.csv, 100 coronal Allen ISH experiments, each with a hand-written category: AMPAR_core, auxiliary, scaffold, trafficking, excitatory, plasticity, control_inhib/glia/struct/vasc)
- **Kind:** table / input definition
- **Old code:** P9_compare_nano_vs_allen_ish.m:31-32, 73-77 (readtable), 294-298 (symbol, experiment_id, category, description per gene)
- **Python:** v2_ish_regions.py:71 (PANEL defaults to gene_targets.csv), main 155-211; v2_panel_build.py:assign_roles 142-168 + main 171-232 (390-gene GO panel, both planes, panel_v2.csv / panel_genes.csv); v2_panel_fetch.py:fetch 52-76
- **Judgement:** The same 100-gene panel is read by default (95 usable: Sst/Chat/Tph2 404 and Olig2/Calb2 off-grid are dropped with reasons in gene_drops.csv). v2 also has a 390-gene panel defined by GO terms rather than by hand, with coronal and sagittal experiments.
- **Second check:** v2_ish_regions.py:71 defaults PANEL to data\gene_targets.csv, and main 155-211 reads symbol, experiment_id and category. On disk, gene_drops.csv lists exactly Sst/Chat/Tph2 ('grid not downloaded'; the 404 is documented elsewhere) and Olig2/Calb2 (off-reference box), which leaves 95 genes. v2_panel_build.py:assign_roles 142-168 and main 171-232 build the GO panel as cited. Minor corrections: gene_targets.csv has 100 coronal genes in 9 categories and no 'control_vasc' gene (the category exists only as a colour key at P9:719-731). The old summary does contain values for Olig2 and Calb2, computed from stretched off-reference grids.

### P9-02 · covered

- **What:** The nano side of the comparison: one adult map, the cohort mean of naive + RWS mice (behavior excluded), P6bis affine-normalised, hemispheres folded (L+R) and z-scored over the whole brain
- **Kind:** analysis (normalisation choice)
- **Old code:** P9_compare_nano_vs_allen_ish.m:24-29 (p8_merged_tag merged_naive_rws_nano, _nosmooth), 92-103 (loads P8 mean_lr_sum), 196-206 (data_z_global via zscore_in_mask 968-973)
- **Python:** v2_ish_compare.py:adult_profile 75-83 (mean over the ten naive+rws adults of per-mouse structure values from v2_region_plot.py:main 186-224, 305-308, region_means_per_mouse.csv), READINGS zref/cref/subref/ratio/sepratio; the replacement itself is tested in v2_ish_compare.py:report 144-153 and figure 183-229 (ish_old_vs_new.png)
- **Judgement:** Same animals (naive+rws, no behavior), but v2 averages per-mouse structure means on each brain's own atlas instead of averaging voxels after the P6bis affine and a whole-brain z-score. It gives five readings instead of one. v2 checks the swap directly: the zref gene ranking reproduces the old one at rho 0.91, while Gria1 moves from rank 22 to rank 10.
- **Second check:** Same ten animals: get_cohort.m:57-66 matches v2_cohort.py:109-110. v2_ish_compare.py:adult_profile 75-83 averages the naive+rws per-mouse structure values written at v2_region_plot.py:305-308. Covered. The first check wrongly credits the rank change to the normalisation. Re-running the v2 zref correlations read-only on P9's own region set (9 divisions, no '-unassigned') gives old-vs-new ranking agreement 0.970 and Gria1 at rank 18/95. Only the full v2 region set gives 0.910 and rank 10. So most of 'Gria1 22 -> 10' comes from the region set (P9-04), not from zref replacing affine + z-score. The P9 nano map also inherits P8's rws->naive linear alignment (P8:227-249), which v2 has no equivalent for and does not need.

### P9-03 · partial

- **What:** channel='auto' setting: run the whole ISH comparison on the autofluorescence cohort map as a control, asking whether the gene correlations are specific to the nanobody rather than tissue autofluorescence
- **Kind:** analysis (control, behind flag)
- **Old code:** P9_compare_nano_vs_allen_ish.m:22-29 (channel flag, selects the P8 auto cache), 568-572 and 626-645 (channel in titles/filenames). Never run: data\comparisons has no merged_naive_rws_auto_vs_ish_* folders.
- **Python:** v2_sep_channel_check.py:main 112-165 (rho_auto_gria, line 138: autofluorescence vs Gria1 per adult); v2_ish_compare.py:correlate 97-109 on the 'ratio' reading (nano/auto, so auto is divided out); v2_beyond_density.py:autofluorescence 222 + step2_covariates 466 (autofluorescence explains ~0% of the nano ordering)
- **Judgement:** Specificity against autofluorescence is shown in two ways: autofluorescence vs Gria1 per mouse, and the full panel rerun on nano/auto. No v2 script ranks the gene panel against the autofluorescence map on its own, which is what this flag would have produced (the old route never ran it either).
- **Second check:** Confirmed that the old route never ran it: data\comparisons has merged_naive_rws_auto (the P8 auto cache) but no *_auto_vs_ish_* folder. v2_sep_channel_check.py:138 gives autofluorescence vs Gria1 per adult. sep_channel_check.csv has rho_auto_gria at about 0.04-0.38 against rho_nano_gria at 0.50-0.72. v2_beyond_density step2 (variance_partition.csv) gives 'autofluorescence only' in-sample R2 0.004. Still, no v2 script ranks the 95 genes against the autofluorescence map. The 'ratio' reading (nano/auto) tests whether the correlations survive dividing auto out, not whether auto alone produces them. Partial is right.

### P9-04 · partial

- **What:** The region set and sampling of the comparison: 228 STRU-leaf structures under 9 forebrain/midbrain divisions (Isocortex, OLF, HPF, CTXsp, STR, PAL, TH, HY, MB), AP planes 100-700 only, both modalities sampled only where the nano cohort has tissue
- **Kind:** analysis (design choice)
- **Old code:** P9_compare_nano_vs_allen_ish.m:56, 59, 105-124, 1034-1090 (build_stru_leaves_by_macro_local); 197-206 and 462-468 (masking to brainMask_merged and slice range)
- **Python:** v2_region_plot.py:main 186-224 (every 'structure'-level term with MIN_VOX=250 in each mouse, line 83); v2_ish_compare.py:gene_profiles 86-95 and correlate 97-109 (MIN_ISH_VOXELS=10, MIN_STRUCTURES=50); the grey-matter rule exists only in v2_beyond_density.py:keep_structure 177-188 (GREY line 167)
- **Judgement:** The v2 gene ranking (v2_ish_compare/words/roles/arms) correlates over up to 237 structures. That set includes hindbrain and cerebellum, fibre tracts (cm, lfbs, mfbs, eps, cbf, scwm), ventricles (VL, V3, AQ) and 17 '..., unassigned' catch-alls. The old route and v2_beyond_density both exclude these. The grey-matter restriction is lost in the ranking scripts. ISH means are also no longer restricted to voxels the nano sections reached.
- **Second check:** Confirmed. gene_correlations.csv has n_structures 213-237 (Cacng8 235). The shared set includes P/MY (17), CB (2), fibre tracts (cm, lfbs, mfbs, eps, cbf, scwm), ventricles (VL, V3, AQ) and 17 '-unassigned' labels. The first check has one error: v2_beyond_density.py:167 GREY keeps P, MY and CB, so it excludes only fibre tracts, ventricles and catch-alls, not hindbrain or cerebellum. The difference matters for the headline. Recomputed read-only from the on-disk tables, Gria1 is rank 10/95 on all structures, 11/95 under the keep_structure grey rule, and 18/95 on P9's nine divisions (gap to Cacng8 0.138 / 0.143 / 0.195; machinery above Gria1 4 / 4 / 7 of 33). The hindbrain and cerebellum inclusion, not white matter, moves Gria1. P9 also sampled ISH only inside the nano tissue (P9:462-465) and only atlas planes 279-879 (crop 180 at line 86 plus range 100-700).

### P9-05 · file handling

- **What:** Download and cache each Allen ISH energy grid, parse the MetaImage header, read the raw volume
- **Kind:** plumbing
- **Old code:** P9_compare_nano_vs_allen_ish.m:330-344, 912-966 (load_allen_ish_grid_native, parse_mhd_header)
- **Python:** v2_panel_fetch.py:fetch 52-76 (404s and off-grid boxes to fetch_failures.csv); v2_ish_regions.py:read_energy 84-107
- **Judgement:** Plumbing, and better in v2. The -1 no-data flag becomes NaN (P9 fed it into interpolation and clamped it to 0). Grids in a non-reference box such as Olig2 and Calb2 are rejected; P9 silently resampled them into the atlas frame, which misplaces every voxel.
- **Second check:** Plumbing, as cited. P9:406-407 interpolates the raw grid (-1 included) and then clamps negatives to 0. v2_ish_regions.py:107 turns -1 into NaN, and :101-102 plus v2_panel_fetch.py:67-69 reject non-reference boxes. P9:406 does stretch any box (including Olig2's 68x40x50 and Calb2's 73x41x53) to allen_full_size.

### P9-06 · partial

- **What:** ISH orientation against the atlas: fixed permute/flip settings plus a per-gene orientation diagnostic (atlas slice, ISH slice, ISH with atlas outlines at mid-AP)
- **Kind:** diagnostic figure
- **Old code:** P9_compare_nano_vs_allen_ish.m:34-37 (ish_raw_permute/flipdim), 413-430 (Diagnostic_ISH_Orientation_<gene>.png, behind save_diagnostic_slices)
- **Python:** v2_ish_regions.py:read_energy 106 (transpose(2,1,0)) and module docstring 20-28 (orientation set by Dice over 48 permutations/flips); v2_ish_regions.py:101-102 and v2_panel_fetch.py:68-69 (header-size guard)
- **Judgement:** Orientation was settled once and is guarded by a header-size check, so the question is answered. The mask-matching check described in the docstring is not code in the repository, and no per-gene visual overlay of ISH on the atlas exists.
- **Second check:** v2_ish_regions.py:106 transposes, and the docstring at 20-28 describes the Dice test over 48 permutations. Grep over all *.py finds no such code: 'Dice' appears only in the docstring and in v2_diagnostics.py:294, which is for nano planes. v2 has no per-gene image of an ISH map on the atlas. Partial is fair (it could also be argued as QC plumbing).

### P9-07 · partial

- **What:** Bad ISH section detection and repair: a native 200 um section whose median energy is below 0.3 x the global median (or has no signal) is replaced by linear interpolation from the nearest good sections; n_bad_sections recorded per gene
- **Kind:** analysis (data-quality correction)
- **Old code:** P9_compare_nano_vs_allen_ish.m:346-383 (threshold 355, repair loop 360-376, n_bad_sections 383)
- **Python:** No per-section QC. Related: v2_ish_regions.py:read_energy 107 (missing voxels to NaN, not interpolated); v2_panel_build.py:allen_experiments 127-139 (failed=false filter, 390 panel only); v2_ish_reliability.py:pair_reliability 74-84 and merge 87-103 (test-retest reliability and rank-merging of replicate experiments)
- **Judgement:** v2 moves the section-quality question to whole-map reliability (median 0.69, Gria1 0.91) and averages replicates, which covers noisy genes on the 390 panel. Dim but non-missing sections are never detected or repaired. The 95-gene table behind v2_ish_compare/words/roles/arms (mostly one coronal experiment per gene) has no protection. The old count was inflated: all 100 genes have n_bad > 0 because empty padding sections were counted.
- **Second check:** v2 has no per-section QC. The first check's 'all 100 genes' should be 97 (3 are NaN), and 'padding' is only part of it. Reading the old gene_result_*.mat files, every gene has at least 2 no-data sections, and inside the analysis range there are 62 missing plus 49 dim flagged sections across 35 genes. Some dim flags are real failures that v2 now silently averages in: Eno2 section 43 at 0.01 against a median of 19.3, Grik3 s39 0.00 against 2.12, Lrfn2 s42 0.00, Nrgn s44 0.03 against 2.78. Others are true regional non-expression that P9 wrongly 'repaired' (for example Slc17a6 s15-18). P9 also copied the last good section over all posterior sections (Cacng8 s48-66, at P9:370-371). Gria1 has 4 missing sections inside the range (36, 38, 39, 41), which v2 NaNs correctly. docs/adult_ish_design.md:573-574 lists ISH QC as still to do.

### P9-08 · missing

- **What:** Native-profile diagnostic figure per gene: median energy per native section before repair, flagged bad sections, profile after repair, threshold line
- **Kind:** diagnostic figure
- **Old code:** P9_compare_nano_vs_allen_ish.m:385-402 (Diagnostic_ISH_NativeProfile_<gene>.png, behind save_diagnostic_slices); n_bad_sections in summary 668-709
- **Python:** none per section; the nearest is v2_ish_reliability.py:figure 166-210 (ish_reliability.png, gene-level reliability)
- **Judgement:** No v2 figure or table shows ISH quality section by section. Reliability is measured only per gene, and only for genes with more than one experiment (218 of 390).
- **Second check:** No v2_*.py reads section medians or plots anything per section. v2_diagnostics.py (sheets 01-09) is nano-only. v2_ish_reliability.py:figure 166-210 is per gene. gene_reliability.csv confirms 218 of 390 genes have a reliability value.

### P9-09 · covered

- **What:** Bring ISH onto the nano grid: upsample 200 um to 10 um (linear), crop to the atlas AP range, mask to the brain, fold hemispheres (L+R) before region means
- **Kind:** analysis (resampling)
- **Old code:** P9_compare_nano_vs_allen_ish.m:404-411 (imresize3), 432-441 (compute_lr_stats.m L+R); helper compute_lr_stats.m:1-28
- **Python:** v2_ish_regions.py:annotation_200 110-118 (annotation sampled every 20th voxel onto the native 200 um grid) and region_means 131-152 (both hemispheres pooled per structure name)
- **Judgement:** Same aim, opposite direction: v2 brings the atlas to the ISH grid instead of interpolating ISH up to 10 um, so it invents no values. A structure needs 3 valid 200 um voxels (10 in the comparison scripts), so the smallest structures drop out. Pooling both hemispheres makes the L+R fold unnecessary.
- **Second check:** v2_ish_regions.py:annotation_200 110-118 and region_means 131-152 match the claim, as do MIN_VOXELS=3 (line 76) and MIN_ISH_VOXELS=10 in the comparison scripts. compute_lr_stats.m:1-28 sums the left half with the flipped right half, so a mean over both hemispheres is equivalent. The v2 direction is also more correct: P9:406 stretches the 67x41x58 grid (13.4 x 8.2 x 11.6 mm) onto the 13.2 x 8.0 x 11.4 mm atlas, a scale error of about 1.5-2.5% (up to about 0.1-0.17 mm at the far end). v2 maps grid voxel i to atlas voxel 20i.

### P9-10 · file handling

- **What:** Whole-brain z-scoring of nano and ISH within the analysis mask; zscore_method setting chooses standard (mean/SD) or robust (median/MAD)
- **Kind:** analysis setting
- **Old code:** P9_compare_nano_vs_allen_ish.m:46, 202-206, 443-450, 968-981 (zscore_in_mask, robust_zscore_in_mask)
- **Python:** no z-scoring of ISH; comparisons are on ranks (v2_ish_compare.py:correlate 97-109); nano uses zref (v2_region_plot.py:value 256-271)
- **Judgement:** Region means of an affinely z-scored volume are an affine transform of the raw region means, so every region correlation in P9 is unchanged by this step and by the robust flag. It only sets the common display scale of the Delta-z video and paired bars, which are separate items (P9-17, P9-21).
- **Second check:** Checked. data_vol_masked (P9:197-200) and data_z_masked (P9:202-206) cover the same voxel set (brainMask_merged within the slice range), and the same holds for ISH at 462-468 against 443-449. zscore_in_mask and robust_zscore_in_mask (968-981) are affine, so region z-means are affine in the raw means, and the region Spearman and Pearson are unchanged. Voxel Pearson uses raw values (453-456). Display only.

### P9-11 · partial

- **What:** Protection against border contamination of region means (registration error, 200 um ISH blur, hippocampal bleed-through): distance-weighted means (bwdist^4, primary) and eroded means (sphere r=3), with Spearman and Pearson computed for both (r_spearman_ero, r_pearson_ero)
- **Kind:** analysis (robustness)
- **Old code:** P9_compare_nano_vs_allen_ish.m:57-58, 126-194 (masks and weights), 208-231 (nano means), 470-493 (ISH means), 509-523 (correlations for both methods)
- **Python:** v2_ish_regions.py:region_means 131-152 and main 164-172, 194 (writes ish_mean_eroded, one 200 um voxel peeled); no v2 script reads ish_mean_eroded, and v2_region_plot.py nano means are plain unweighted means
- **Judgement:** v2 computes the eroded ISH mean so a downstream script could check whether erosion changes the ranking (its docstring says so, lines 38-41), but nothing does, and the nano side is never eroded or weighted. The owner dropped the distance-weight comparison itself (P9-24). The underlying question, whether small structures take their neighbours' signal, is unanswered for the gene ranking. v2_beyond_controls.py:control_b_size 149-162 checks structure size only for the residual claim.
- **Second check:** Citations correct. ish_mean_eroded is written at v2_ish_regions.py:194, and grep shows no other v2 script reads it. The design doc (docs/adult_ish_design.md:526-529, D2) itself asks that this be tested. The check is cheap and the ISH half is currently benign: on the v2 table, full vs eroded ISH gives ranking agreement 0.998 (Gria1 rank 10 -> 7), and the old DW vs eroded Spearman ranking agreed at 0.997. The nano side is never eroded or weighted (v2_region_plot.py:210-223), so it stays partial until the check is in code.

### P9-12 · covered

- **What:** Per-gene region-level Spearman correlation between the nano map and the gene's ISH map (distance-weighted, the primary statistic)
- **Kind:** statistic
- **Old code:** P9_compare_nano_vs_allen_ish.m:508-527 (r_spearman_dw is r_spearman)
- **Python:** v2_ish_compare.py:correlate 97-109 (Spearman over shared structures, per gene and per reading, written to gene_correlations.csv); also v2_ish_roles.py:main 259-276 (role_summary.csv) and v2_ish_panel_test.py:main 163-175 (plain and partial rho on the 390 panel)
- **Judgement:** Same statistic on a different nano normalisation (five readings) and a different region set (see P9-04). It reproduces the old ordering at rho 0.91 for zref (Cacng8 +0.80, Gria1 +0.66).
- **Second check:** v2_ish_compare.py:correlate 97-109 computes Spearman per gene and reading, and v2_ish_roles.py:259-276 and v2_ish_panel_test.py:163-175 are as cited. Recomputed: zref old-vs-new 0.910, Cacng8 +0.803, Gria1 +0.665. Same statistic and question. The values depend on the region set (P9-04).

### P9-13 · covered

- **What:** Per-gene region-level Pearson correlation (distance-weighted)
- **Kind:** statistic
- **Old code:** P9_compare_nano_vs_allen_ish.m:519, 527 (r_pearson_dw is r_pearson)
- **Python:** v2_ish_compare.py:correlate 97-109 (Spearman only)
- **Judgement:** v2 never computes Pearson; the rank-only choice is documented in v2_beyond_density.py (docstring 'Ranks everywhere', about 206-210) because Allen energy has an arbitrary per-experiment scale. The question, whether a gene's regional profile predicts the nano map, is answered at least as well by Spearman. Only the outlier-sensitive linear magnitude (old Pearson far below Spearman for e.g. Dlg2, 0.42 vs 0.75) is not reported.
- **Second check:** Covered in spirit by Spearman. Citation error: 'Ranks everywhere' is at v2_beyond_density.py:103-107, not 206-210. Note the old data: Spearman-DW vs Pearson-DW gene rankings agree at only 0.81, Gria1 is rank 9 under Pearson against 22 under Spearman, and Dlg2 is 0.42 against 0.75 (verified). The old headline was sensitive to the choice of statistic (see P9-25).

### P9-14 · dropped

- **What:** Per-gene voxel-level Pearson correlation between nano and ISH over the analysis mask (ISH>0, at least 100 voxels)
- **Kind:** statistic
- **Old code:** P9_compare_nano_vs_allen_ish.m:279-280, 452-459
- **Python:** none (v2_beyond_density.py docstring 'Structures, not voxels', about 202-204)
- **Judgement:** Dropped by the owner on 29 Sep. v2 gives the reason: ISH exists only at 200 um and in a different mouse, so a voxel comparison would mostly measure registration error.
- **Second check:** The 29 Sep drop decision lists 'voxel-level ISH correlation' as dropped on purpose. Citation error: 'Structures, not voxels' is at v2_beyond_density.py:99-101, not 202-204.

### P9-15 · covered

- **What:** Per-gene region table: for each structure, the nano mean, ISH mean, both z-means and voxel count (Region_NanoVsISH_Table_<gene>.csv)
- **Kind:** table
- **Old code:** P9_compare_nano_vs_allen_ish.m:495-506, 546-547
- **Python:** v2_ish_regions.py:main 190-202 (gene_region_table.csv: ISH mean, eroded mean, n_voxels, coverage per gene and structure) with v2_region_plot.py:main 305-308 (region_means_per_mouse.csv, nano per mouse and structure)
- **Judgement:** The same numbers are in two long tables that join on structure name, and v2 adds per-mouse values and ISH coverage. No single per-gene side-by-side file is written.
- **Second check:** gene_region_table.csv (written at v2_ish_regions.py:190-202) holds ish_mean, ish_mean_eroded, n_voxels and coverage per gene and structure. region_means_per_mouse.csv (v2_region_plot.py:305-308) holds nano per mouse and structure with division. The two join on structure name. The z-mean columns are affine and so redundant (P9-10).

### P9-16 · file handling

- **What:** Per-gene saved results: gene_result_<gene>.mat (both methods' r values, native medians, region tables, settings) and ish_lr_sum_<gene>.mat (the folded 10 um ISH volume)
- **Kind:** plumbing / cache
- **Old code:** P9_compare_nano_vs_allen_ish.m:307-308, 533-544
- **Python:** v2_ish_compare.py:main 245-249 (gene_correlations.csv); ISH volumes are not saved
- **Judgement:** Storage and caching only. v2 keeps the scientific content in CSVs and has no 10 um ISH volumes.
- **Second check:** Storage and cache only (P9:533-544). v2_ish_compare.py:245-249 writes gene_correlations.csv.

### P9-17 · partial

- **What:** Per-gene three-panel comparison video, one coronal frame per plane: nano L+R map | ISH L+R map | voxelwise Delta-z (nano z minus ISH z), with region acronyms on the difference panel. It shows where in the brain nano exceeds or falls short of that gene's expression.
- **Kind:** video
- **Old code:** P9_compare_nano_vs_allen_ish.m:40, 48-53, 62, 450 (diff_z), 559-574, 983-1032 (write_3panel_comparison_video); get_color2color_colormap.m
- **Python:** none per gene. Nearest: v2_beyond_regression.py:panel_f 136-182 (F_maps.png: observed, predicted and residual painted on 3 coronal planes, for the combined abundance + density model)
- **Judgement:** Where the map departs from gene-based prediction is shown spatially in v2, but only for the combined multi-gene model, on three planes, at structure resolution. No per-gene spatial comparison exists. The refactor plan lists this Delta-z video as optional, not as dropped.
- **Second check:** Partial agreed: v2_beyond_regression.py:panel_f 136-182 paints observed, predicted and residual of the combined model on 3 planes. No per-gene spatial view exists, not even for Gria1. The first check was wrong about the plan: draft 1 listed only 'per-gene scatter plots (Gria1, Cacng8)' as optional. The 'threshold / delta-z video' was optional only in the working notes of 29 Sep. The plan should record it (now decision S4), or it will be lost silently at retirement.

### P9-18 · file handling

- **What:** Region label centroids (per plane, smoothed) for the acronym overlay in the video
- **Kind:** plumbing
- **Old code:** P9_compare_nano_vs_allen_ish.m:233-277
- **Python:** v2_video.py / v2_video_compare.py draw acronyms their own way (young-vs-adult videos only)
- **Judgement:** Only a drawing aid for the video in P9-17.
- **Second check:** A drawing aid only (P9:233-277). v2_video.py:47,77 and v2_video_compare.py:64 draw acronyms their own way.

### P9-19 · missing

- **What:** Per-gene region scatter: ISH z vs nano z, one labelled dot per structure coloured by division, identity line, title with rS / rP / rVox
- **Kind:** figure
- **Old code:** P9_compare_nano_vs_allen_ish.m:41, 282-288 (division colours), 576-619 (Scatter_NanoVsISH_<gene>.png/.fig)
- **Python:** none per gene. Nearest: v2_beyond_regression.py:panel_e 92-133 (observed vs predicted for the combined model); v2_ish_compare.py:figure 183-229 (old vs new rho per gene, not per structure)
- **Judgement:** No v2 figure shows the structure-by-structure relationship for any single gene, e.g. Gria1 or Cacng8. The refactor plan lists Gria1 and Cacng8 scatters as an optional addition.
- **Second check:** Grep of every scatter( call in v2_*.py finds no structure-by-structure scatter of one gene against nano. The only structure-level scatters are v2_beyond_regression.py:99 (combined model), channel-vs-channel scatters (v2_sep_channel_check.py:188,196; v2_adult_arms.py:210) and control plots. Draft 1 of the plan listed it as optional.

### P9-20 · missing

- **What:** Within-division relationship: a separate regression of nano on ISH for each division, with 95% confidence band and slope in the legend, asking whether a gene tracks the nano map inside each division or only between divisions
- **Kind:** statistic + figure
- **Old code:** P9_compare_nano_vs_allen_ish.m:585-605 (per-division polyfit, t-based CI, slope label)
- **Python:** none. Related only: v2_beyond_controls.py:control_a_space 122-146 (spatial-gradient control for the residual claim, not per gene)
- **Judgement:** All v2 gene correlations are brain-wide across structures, so a correlation driven purely by cortex vs thalamus differences cannot be told apart. The old slopes were descriptive, but they were the only within-division view of gene-map agreement.
- **Second check:** Grep for division/within in v2_ish_*.py and v2_beyond_*.py finds division used only for the keep_structure filter (v2_beyond_density.py:177-199). No per-division gene-vs-map fit exists. control_a_space (v2_beyond_controls.py:122-146) tests spatial position on the residual, not gene tracking within divisions.

### P9-21 · partial

- **What:** Per-gene paired bar charts: nano z and ISH z for every structure side by side, sorted by nano and then by ISH, with moving-average trend lines, showing which structures agree or disagree
- **Kind:** figure
- **Old code:** P9_compare_nano_vs_allen_ish.m:42, 621-650 (PairedBar_<channel>VsISH_sort<nano|ISH>_<gene>.png/.fig)
- **Python:** none per gene. Nearest: v2_beyond_density.py:step4_where 544-608 (residual_by_structure.csv) and v2_beyond_regression.py:main 207-215 (regression_table.csv: observed, predicted, residual per structure for the combined model)
- **Judgement:** Which structures depart from gene-based prediction is answered in v2 for the combined abundance + density model, with replication across animals, which is stronger. The per-gene structure-by-structure agreement view is gone.
- **Second check:** residual_by_structure.csv (v2_beyond_density.py:551-557) and regression_table.csv (v2_beyond_regression.py:211-217) give structure-level agreement and disagreement for the combined abundance+density model only. No per-gene paired view exists. The 'replication across animals' is a whole-residual statistic (step3, 527-541), not per structure.

### P9-22 · covered

- **What:** Cross-gene summary table: per gene, Spearman and Pearson (distance-weighted and eroded), voxel Pearson and n_bad_sections, sorted by distance-weighted Spearman
- **Kind:** table
- **Old code:** P9_compare_nano_vs_allen_ish.m:656-712 (gene_panel_summary.csv/.mat)
- **Python:** v2_ish_compare.py:correlate 97-109 and main 245-249 (gene_correlations.csv: rho and n_structures per gene and reading); v2_ish_reliability.py (gene_reliability.csv); v2_ish_panel_test.py (panel_test.csv)
- **Judgement:** Covered for the Spearman ranking, across five readings; Pearson, eroded and voxel columns and n_bad are not kept. Dependency to keep when retiring: v2_ish_compare.py:61-62 and old_ranking 112-118 read this exact P9 file to draw ish_old_vs_new.png. The file (not the script) must stay under data\comparisons, or the old-vs-new check is silently skipped.
- **Second check:** The Spearman ranking is covered (gene_correlations.csv), and v2_ish_compare.py:61-62 and 112-118 read the old file. An extra dependency the first check missed: that file sits in merged_naive_rws_vs_ish_summary_nosmooth, written by an older P9 whose tag had no channel. Current P9 builds p8_merged_tag = 'merged_naive_rws_nano' (line 28), so a rerun would write to ..._nano_vs_ish_summary_nosmooth, which v2 never reads. data\comparisons has 101 merged_naive_rws_vs_ish_* folders and no _nano_ ones. The old file must be preserved as a frozen artefact, not regenerated.

### P9-23 · partial (first pass: covered)

- **What:** Headline gene ranking: which genes best predict the nano map, and the claim that surface GluA1 tracks AMPAR auxiliary/trafficking machinery (Cacng8, Grm5, Cnih2, Dlg2) better than its own mRNA (Gria1)
- **Kind:** scientific claim / ranking
- **Old code:** P9_compare_nano_vs_allen_ish.m:710 (sort by r_spearman_dw), 751-768, 816-903 (figures the claim was read from)
- **Python:** v2_ish_compare.py:report 129-180 (top 5 per reading, ranks of Cacng8 and Gria1, best-machinery-minus-Gria1 gap, machinery above Gria1, old vs new); v2_ish_roles.py:main 243-325 (subunit vs localisation commonality + exact 3,876-split permutation + expression sensitivity 209-240); v2_ish_panel_test.py:main 145-252 (390-gene powered test, expression-matched controls, positive control); v2_ish_arms.py:partial_spearman 116-124 (rho(map, gene | Gria1))
- **Judgement:** The same question with tests P9 never had, and the answer changed. Cacng8 is still first, but Gria1 rises to rank 10 of 95 and the gap halves. The powered test finds localisation genes no better than matched postsynaptic controls (p=0.74), so trafficking genes as a class do not beat Gria1.
- **Second check:** The question is answered with stronger tests, and the cited functions exist as described (v2_ish_compare.py:report 129-180, v2_ish_roles.py:243-325, v2_ish_panel_test.py:145-252, v2_ish_arms.py:116-124). But the headline answer depends on a region-set choice that v2 has not settled. All ranking scripts correlate over up to 237 structures including P/MY/CB, white matter and ventricles. On P9's own nine-division set, recomputed read-only from the on-disk tables, Gria1 is rank 18/95 (rho 0.579), the Cacng8-Gria1 gap is 0.195 and 7/33 machinery genes sit above Gria1. The 'answer changed' statement (rank 10, gap halves) is therefore mostly a region-set effect. To settle it: rerun v2_ish_compare, v2_ish_roles and v2_ish_panel_test with keep_structure (or a declared region set) and report both.

### P9-24 · partial (first pass: dropped)

- **What:** Eroded vs distance-weighted comparison: grouped bars of Spearman under both aggregation methods for every gene
- **Kind:** diagnostic figure
- **Old code:** P9_compare_nano_vs_allen_ish.m:735-749 (correlation_eroded_vs_distweight.png/.fig)
- **Python:** none
- **Judgement:** This is the distance-weight diagnostic the owner dropped on 29 Sep. The broader border-contamination question it served is tracked separately in P9-11.
- **Second check:** The 29 Sep drop decision says 'distance-weight diagnostic', covering P8, P9 and P10 together. P8 has a figure with exactly that name: Diagnostic_DistWeight_slice%d_pow%d.png, titled 'Distance weighting diagnostic' (P8:1397-1491; the outputs exist in comparisons\merged_naive_rws_auto). P9's correlation_eroded_vs_distweight is a different thing, a robustness test of the gene ranking to border handling. So the drop cannot be confirmed to cover it. The question is folded into P9-11 and is cheap to answer (ranking agreement old 0.997, v2 full vs eroded 0.998). Ask the owner, or add the one-line eroded check.

### P9-25 · partial

- **What:** Five-metric cross-gene bar chart (Spearman DW/eroded, Pearson DW/eroded, voxel Pearson per gene): is the gene ranking robust to the choice of statistic and aggregation?
- **Kind:** figure
- **Old code:** P9_compare_nano_vs_allen_ish.m:751-768 (correlation_barchart_4metrics.png/.fig)
- **Python:** v2_ish_compare.py:report 144-153 and figure 183-229 (ranking robustness across the five readings and against the old route, ish_old_vs_new.png); gene_correlations.csv holds all readings
- **Judgement:** v2 tests robustness to normalisation (five readings, old vs new) instead of robustness to the statistic. Voxel Pearson and the distance-weighted/eroded pair are dropped by design, but Spearman vs Pearson, and full vs eroded ISH means (computed but unused, P9-11), are compared nowhere.
- **Second check:** v2_ish_compare.py:144-153 and 183-229 test robustness to normalisation, not to the statistic. The old data show that the statistic is the part that mattered: DW vs eroded agree at 0.997, while Spearman vs Pearson agree at only 0.813, with Gria1 rank 22 against 9. So Spearman vs Pearson is a substantive robustness gap, not cosmetic. It is still partial because v2's rank-only choice is argued in the v2_beyond_density.py docstring (lines 103-107).

### P9-26 · covered

- **What:** Split of the 'auxiliary' category into forebrain-expressed auxiliaries (Cacng8, Cacng3, Cnih2, Cnih3, Grm5) and the rest (split_auxiliary setting)
- **Kind:** analysis (gene grouping, behind flag)
- **Old code:** P9_compare_nano_vs_allen_ish.m:772-813
- **Python:** v2_ish_roles.py:ROLES 95-132 (curated by protein function: TARPs, cornichons and Shisa9 as localisation, Grm1-5 as other_glut_r) and sensitivity 209-240 (drop low-expression family genes on expression alone, never on rho); v2_ish_panel_test.py:greedy_match 114-129 (expression matching)
- **Judgement:** The old split was a hand-made, post-hoc grouping by where genes are expressed, and it put a metabotropic receptor (Grm5) among the AMPAR auxiliaries. v2 answers the intended question, whether barely-expressed auxiliaries dilute the class, with a rho-blind expression filter and expression-matched controls. The result is still negative (p=0.31 after filtering).
- **Second check:** The old split (P9:772-813) is a display regrouping by forebrain expression. v2_ish_roles.py:ROLES 95-132 moves Grm1-5 to other_glut_r. sensitivity 209-240 drops low-expression family genes by median energy only (rho-blind) and reruns the permutation (earlier result: p = 0.3128). v2_ish_panel_test.py:greedy_match 114-129 adds expression matching. Same question (dilution by barely expressed auxiliaries), better controlled.

### P9-27 · covered

- **What:** Category distribution plots: correlation of each gene with the nano map grouped by category (violin + dots + gene labels, sorted by median), for DW Spearman, eroded Spearman and voxel Pearson
- **Kind:** figure
- **Old code:** P9_compare_nano_vs_allen_ish.m:816-903 (violin_spearman_dw / violin_spearman_ero / violin_voxel_pearson .png/.fig); helper plot_violinplot.m
- **Python:** v2_ish_roles.py:figure 328-383 panel 1 (rho by curated role, dots + median, zref) and main 259-292 (role_summary.csv, all readings); v2_ish_words.py:figure 259-276 (GO term and word groups, bootstrap intervals); v2_ish_panel_test.py:figure 255-302
- **Judgement:** The DW-Spearman view is reproduced with curated roles (plus the original control labels) instead of the circular hand-made categories, and extended with GO groupings. The eroded and voxel versions fall under dropped items. Gene names are not printed on the dots.
- **Second check:** v2_ish_roles.py:figure 328-383 panel 1 shows rho per curated role with dots and median for zref, keeping the control_* labels (role_of 162-167). role_summary.csv holds all readings. v2_ish_words.py:figure 259-276 and v2_ish_panel_test.py:255-302 are as cited. Gene names are not printed. The old hand categories themselves are no longer plotted, only the curated roles.

### P9-28 · covered

- **What:** One-way ANOVA: does a gene's correlation with the nano map depend on its category? (p and F in each violin title)
- **Kind:** statistic
- **Old code:** P9_compare_nano_vs_allen_ish.m:884-894 (anova1 on the category index)
- **Python:** v2_ish_words.py:test_features 172-200 (Mann-Whitney, bootstrap interval on the median gap, BH across GO features); v2_ish_roles.py:permutation 191-206 (exact within-family split null); v2_ish_panel_test.py:two_sample 132-142 (label-permutation test, expression-matched, positive control 230-244)
- **Judgement:** No omnibus test over the hand-made categories. The question is answered with targeted contrasts on categories fixed in advance (GO, curated roles), with nulls that handle co-expression. The ANOVA treated co-expressed genes as independent and so was anticonservative (Fulcher 2021); the old result was p=6.6e-4.
- **Second check:** The cited contrasts exist: v2_ish_words.py:172-200, v2_ish_roles.py:191-206, v2_ish_panel_test.py:132-142 and 230-244. Covered in spirit by targeted, pre-specified contrasts. Minor: p=6.6e-4 is the 50-gene result (P8/P9 notes of April 2026), not necessarily what the 100-gene violin titles show.

### P9-29 · covered

- **What:** Specificity against control genes: inhibitory, glial, structural and vascular control genes should sit near zero while excitatory/AMPAR genes correlate
- **Kind:** scientific claim
- **Old code:** P9_compare_nano_vs_allen_ish.m:719-731, 787-799 (control categories in the violins); key result 'controls cluster near zero' (P8/P9 notes of April 2026)
- **Python:** v2_ish_roles.py:role_of 162-167 (keeps control_inhib/glia/struct) and figure panel 1 328-348; v2_ish_words.py:PREDICTED 83-85 (gabaergic, inhibitory, presynaptic words named in advance); v2_ish_roles.py:main 318-322 (presynaptic machinery as specificity control); v2_ish_panel_test.py (control_psd pool)
- **Judgement:** The control labels are kept and shown per role, and v2 adds sharper specificity contrasts: postsynaptic vs presynaptic words (+0.34 vs +0.04), and presynaptic vesicle machinery vs AMPAR localisation. control_vasc was already absent from the old violins (not in cats_ordered).
- **Second check:** Control labels are kept (v2_ish_roles.py:162-167) and plotted (328-348). The presynaptic and word contrasts are as cited. Checked on disk: in old DW and new zref alike, control_inhib has a low median (+0.144 / +0.122) and control_glia sits near zero, but 'controls cluster near zero' is not true gene by gene. Chrm1 +0.67, Htr3a +0.67 and Drd1 +0.66 in v2, and v2 shows this. control_vasc has no genes in gene_targets.csv, so its absence from the violins is moot.

### P9-30 · file handling

- **What:** Batch bookkeeping: skip genes already processed, cache ROI masks and distance weights to disk, record NaN for genes whose download or shape check fails
- **Kind:** plumbing
- **Old code:** P9_compare_nano_vs_allen_ish.m:135-194 (roi_masks cache), 310-328 (skip if done), 335-344 and 435-441 (failure handling)
- **Python:** v2_panel_fetch.py:already_there 39-42 and fetch 52-76; v2_ish_regions.py:main 178-187, 207-211 (gene_drops.csv with reasons); v2_panel_build.py:cached 100-112
- **Judgement:** Caching and error handling only. v2 records each dropped gene with its reason rather than a bare NaN.
- **Second check:** Caching and failure handling (P9:135-194, 310-344, 435-441). The v2 equivalents are at the cited lines (v2_panel_fetch.py:39-42, 52-76; v2_ish_regions.py:178-187, 207-211; v2_panel_build.py:100-112).

### P9-X1 · missing

- **What:** Region-set sensitivity of the headline gene ranking. The old answer is a forebrain+midbrain answer: 9 divisions, atlas planes 279-879, ISH only where the nano cohort has tissue. Nothing in v2 restricts the ranking to that set or reports how the ranking depends on it. Read-only recomputation from the on-disk tables: v2 as is gives Gria1 rank 10/95 and old-vs-new 0.910. The keep_structure grey rule gives rank 11 and 0.927. P9's nine divisions give rank 18, gap 0.195, 7/33 machinery above Gria1 and old-vs-new 0.970.
- **Old code:** P9_compare_nano_vs_allen_ish.m:56, 59, 86-87, 105-124, 462-468, 1034-1090
- **Python:** none; v2_ish_compare.py:correlate 97-109, v2_ish_roles.py:259-276 and v2_ish_panel_test.py:158-175 use every shared structure; v2_beyond_density.py:keep_structure 177-188 is not applied in the ranking scripts
- **Judgement:** This is a test to add before retiring. The 'Gria1 climbs to rank 10' change comes mostly from including hindbrain and cerebellum, not from zref. v2 should declare its region set for the ranking and report the P9-equivalent set beside it.

### P9-X2 · file handling

- **What:** Provenance of the old reference ranking. The gene_panel_summary.csv that v2 compares against was written by an older P9 whose output tag had no channel. Today's P9 would write to merged_naive_rws_nano_vs_ish_summary_nosmooth, a folder v2 never reads. The file also holds values for Olig2 and Calb2 computed from stretched off-reference grids.
- **Old code:** P9_compare_nano_vs_allen_ish.m:28, 69-70, 406, 711
- **Python:** v2_ish_compare.py:61-62, old_ranking 112-118
- **Judgement:** A retirement dependency: keep data\comparisons\merged_naive_rws_vs_ish_summary_nosmooth\gene_panel_summary.csv as a frozen artefact, and do not expect a P9 rerun to refresh it. When quoting old-vs-new agreement, drop Olig2 and Calb2 (v2 already excludes them).

### P9-X3 · file handling

- **What:** ISH grid-to-atlas placement. P9 stretches the 200 um grid (67x41x58 = 13.4 x 8.2 x 11.6 mm) onto the 10 um atlas box (13.2 x 8.0 x 11.4 mm) with imresize3. That is a systematic 1.5-2.5% scale error that displaces ISH by up to about 0.1-0.17 mm near the far end of the analysis range, and it is one reason the old and new rho values differ.
- **Old code:** P9_compare_nano_vs_allen_ish.m:406-411
- **Python:** v2_ish_regions.py:annotation_200 110-118 (stride-20 sampling, zero offset per the docstring at 20-28)
- **Judgement:** Fixed in v2 by construction. Record it in the handover notes as an old-route defect, so that differences in old-vs-new numbers are not read as biology.

## P10: P10_compare_nano_vs_auto.m (nano against autofluorescence)

### P10-01 · covered

- **What:** Overarching question (Sami, 30 Apr): is the nano spatial pattern its own signal, distinct from the autofluorescence pattern of the same brains, rather than autofluorescence bleed-through or an acquisition artefact?
- **Kind:** question
- **Old code:** P10_compare_nano_vs_auto.m:5-13 (header), answered through the tests at 117-148 and 181-210
- **Python:** v2_sep_channel_check.py:main (rho_nano_auto per adult L136; dynamic ranges L131-133; figure L168-224) + v2_beyond_density.py:autofluorescence L222-251 and step2_covariates L466-524 (autofluorescence-only model L485) + v2_beyond_figures.py:panel_a L136
- **Judgement:** Python answers this across the whole map, per mouse, with autofluorescence read from the current registered tiffs: nano~auto rho is about +0.26 +- 0.11 per adult, nano's range is 1.93 log2 against auto's 1.07, and autofluorescence alone explains about 0% of the zref structure ordering (cross-validated). That is stronger than P10 for the brain-wide claim, but it does not say WHICH structures have nano above auto (P10-02/05).
- **Second check:** The cited code exists and computes what the first check says. v2_sep_channel_check.py L131-136 computes per-adult range_nano, range_auto and rho_nano_auto. The CSV data/adult_v2/arms/sep_channel_check.csv gives rho_nano_auto = 0.21, 0.21, 0.14, 0.49, 0.43, 0.23, 0.25, 0.28, 0.17, 0.16, mean 0.257 (sd 0.11). Ranges are 1.93 (nano) against 1.07 (auto). v2_beyond_density.py L485 runs the 'autofluorescence only' model; variance_partition.csv gives CV R2 -0.043 and in-sample 0.004. The first check did not cite v2_beyond_controls.py control A (header L19-20, smooth spatial gradient = clearing/illumination artefact), which answers the 'acquisition artefact' half of the question. Caveats: (1) the sep_channel_check figure (L168-224) does not draw rho_nano_auto; it is only printed (L160) and written to the CSV. (2) That rho runs over every structure, fibre tracts and ventricles included. (3) The beyond_density test uses cohort-mean ranks on 125 ISH-covered grey-matter structures only. The brain-wide question is answered. Which structures, and in which direction, is not (P10-02/05).

### P10-02 · missing

- **What:** Per-structure paired one-sided Wilcoxon signed-rank test (nano > auto) across the 10 adults (naive+rws), one test per STRU leaf, run where at least 3 mice are valid; prints the count significant
- **Kind:** test
- **Old code:** P10_compare_nano_vs_auto.m:117-148
- **Python:** none. The raw material exists: per-mouse per-structure log2(nano/auto) for all 10 adults, in comparisons_v2/young_vs_adult/region_means_per_mouse.csv (reading 'ratio', v2_region_plot.py:main L203-208, L256-271, written L305-308) and in adult_v2/arms/region_means_arms.csv (v2_adult_arms.py). Only young-vs-adult Welch/Mann-Whitney tests are run on it.
- **Judgement:** In P10 each channel was first scaled to its own mid-brain tissue mean (P6bis affine per mouse, then the P8 common factor over slices 300-500), so 'nano > auto' means the region is enriched relative to brain level more in nano than in auto. The v2 'ratio' reading is raw background-subtracted counts whose zero depends on exposure (adult cortex nano/auto 0.74-1.46), so the planned 'signed-rank on the per-mouse log2 ratio' must centre each mouse first (e.g. minus its median over structures) or it tests exposure, not enrichment. Also, P10's auto came from the stale adult auto_4d.mat (19 Nov 2025, about 2 voxels off the 26-27 Nov nano registration), so the Python numbers will differ.
- **Second check:** Searched all *.py for wilcoxon/signrank/ttest_rel/ttest_1samp/binomtest. The only hit is v2_ish_arms.py L199, a paired test over GENES (sepratio vs sepauto), not structures. v2_region_plot.py L287-295 and v2_region_groups.py L230-247 test only young vs adult (Welch/Mann-Whitney). Recomputed from data/comparisons/nano_vs_auto/NanoVsAuto_paired_signrank_nosmooth_uncorrected.csv: 228 leaves, 206 tested, 28 at p<0.05. Minimum p = 0.0009765625 = 2^-10, which confirms n = 10 and the exact test. Spirit note is right about centring. v2's 'ratio' is a mean of voxelwise sig/smoothed-auto, clipped to +-20 (v2_region_plot.py L203-208). Its level varies by mouse: adult isocortex is about 0.79-1.52 linear from region_means_per_mouse.csv, and the per-mouse median over structures is -0.30 to +0.38 log2. So it must be centred per mouse. A cleaner input for the port is log2(mean nano) - log2(mean auto) per structure, which v2_sep_channel_check.mouse_channels (L75-94) already computes. Also for the port: P10's reference (mid-brain tissue mean) is dominated by HPF (nano 14.5 vs auto 1.9 at division level), so most non-HPF divisions look nano<auto. With a median-centred contrast, both tails should be reported, not only nano>auto.

### P10-03 · missing

- **What:** Multiple-comparison handling for the per-structure tests: optional Bonferroni across the regions tested (apply_bonferroni, default false, so alpha = 0.05 uncorrected); the saved outputs are the uncorrected version
- **Kind:** correction
- **Old code:** P10_compare_nano_vs_auto.m:36-37, 136-148, 203-208
- **Python:** none for this family; a reusable BH is in v2_region_plot.py:bh_fdr L134-152
- **Judgement:** With n = 10 the smallest possible one-sided exact signed-rank p is 1/1024 = 9.8e-4, which is above the Bonferroni threshold 0.05/206 = 2.4e-4, so P10's Bonferroni option can never reject anything. Recomputed from the old CSVs: paired test 28/206 uncorrected, 23 at BH q<0.05, 0 Bonferroni; delta-z 21, 11, 0. The planned BH is the workable choice.
- **Second check:** Recomputed from the old CSVs. Paired test: 28 uncorrected, 23 at BH q<0.05, 0 Bonferroni. Delta-z: 21, 11, 0. The Bonferroni argument holds: 0.05/206 = 2.43e-4 < 1/1024. bh_fdr exists at v2_region_plot.py L134-152, but it is applied only to the young-vs-adult families (v2_region_plot.py L297-303; v2_region_groups.py L248-252). No nano-vs-auto family exists to correct.

### P10-04 · missing

- **What:** Table NanoVsAuto_paired_signrank_*.csv: per structure the acronym, division, p, significant flag, and cohort-mean nano and cohort-mean auto, sorted by p
- **Kind:** table
- **Old code:** P10_compare_nano_vs_auto.m:150-160 (output data/comparisons/nano_vs_auto/NanoVsAuto_paired_signrank_nosmooth_uncorrected.csv)
- **Python:** none. v2_sep_channel_check.py:mouse_channels L75-94 computes per-mouse log2 structure means of nano and auto separately, but only per-mouse summary rows are written (L144-149).
- **Judgement:** No Python table puts per-structure nano and auto (or their centred contrast) next to a nano-vs-auto p value. The old table flags mostly HPF, STR, CTXsp and TR as nano > auto; isocortex does not pass.
- **Second check:** The status holds. v2_sep_channel_check.py writes only per-mouse summary rows (L121-149). No per-structure nano/auto table exists anywhere in v2, and no v2 script writes the per-structure auto means. But the first check is partly wrong. 'Isocortex does not pass' is true only at the division level (p=1 in NanoVsAuto_macro_signrank). At the structure level ILA passes (p=0.00195, BH q about 0.02) and so does PL (p=0.042 uncorrected). The 28 uncorrected hits by division: HPF 9, STR 9, CTXsp 4, OLF 3 (TR, TT, DP), Isocortex 2 (ILA, PL), PAL 1 (TRS).

### P10-05 · missing

- **What:** Per-structure delta-z test: each mouse's region value z-scored against its own channel's brain mean and SD (cohort-mean LR-sum inside brainMask_analysis, slices 100-700); contrast = nano_z - auto_z; one-sided signed-rank that the contrast is > 0. Asks whether a structure sits further above its brain's average in nano than in autofluorescence, whatever the scale of each channel
- **Kind:** test
- **Old code:** P10_compare_nano_vs_auto.m:162-210
- **Python:** none. Only the nano half exists, as zref (each brain's own median and p90-p10 over structures, v2_region_plot.py:main L245-264). No autofluorescence zref exists anywhere, and nothing contrasts the two per structure.
- **Judgement:** This is the scale-free version of P10-02 and arguably the more meaningful of the two. A v2 equivalent would be zref applied to log2(auto) with the same per-brain formula, then zref_nano - zref_auto per mouse, signed-rank against 0. Note P10 applied cohort-level mu/sigma (not per mouse) on a linear, not log, scale.
- **Second check:** zref is built from sig only: v2_cohort.py mouse_scalars L153-171 and mouse_modes L201-204, v2_region_plot.py L245-264, v2_region_groups.py L205-227. The only auto-side spread is range_auto in v2_sep_channel_check.py L132, the scale half of an auto zref, with no per-structure z. The nearest thing is v2_beyond_density.py step4 (L544-557, residual_by_structure.csv): nano rank minus the rank predicted by abundance + markers + psd_pc1 + autofluorescence. It is rank-based, cohort-level, untested per structure, and driven by the gene covariates (auto alone CV R2 about 0), so it does not answer 'further above brain average in nano than in auto'. Formula check: P10 L187-189 applies one cohort mu/sigma from the voxel distribution of the cohort-mean LR-sum (L172-177) to every mouse, on a linear scale, as the first check says.

### P10-06 · missing

- **What:** Table NanoVsAuto_deltaz_signrank_*.csv: per structure the p, significant flag, mean contrast, mean nano_z and mean auto_z
- **Kind:** table
- **Old code:** P10_compare_nano_vs_auto.m:212-223
- **Python:** none
- **Judgement:** Goes with P10-05; nothing in Python writes an autofluorescence-referenced enrichment per structure. Old result: HPF fields, DG, IG and LSc lead (11 at BH q<0.05).
- **Second check:** Nothing in v2 writes a per-structure autofluorescence-referenced contrast. Old results verified: 21 uncorrected, 11 BH. The p=1/1024 leaders are CA3, DG, FC, IG and LSc. Caveat missing from the first check: 6 of the 21 uncorrected hits are PVpo, PVa, PVi, SFO, ME (HY) and DT (MB). In these nano_z is below the brain mean (about -0.41 to -0.48) and the contrast comes entirely from very low auto_z (about -1.7 to -2.35). So delta-z also flags low-autofluorescence periventricular tissue. A port must keep nano_z and auto_z beside the contrast (as P10 L216 does) or require nano_z > 0.

### P10-07 · partial

- **What:** Division-level pooling: per mouse, the 9 DIVI divisions (Isocortex, OLF, HPF, CTXsp, STR, PAL, TH, HY, MB) averaged over their leaves, weighted by each leaf's eroded voxel count (roi_px_ero), for nano, auto and the delta-z contrast
- **Kind:** aggregation
- **Old code:** P10_compare_nano_vs_auto.m:108-115, 225-277
- **Python:** v2_region_groups.py:main L127-212 (per-mouse voxel-weighted means per division TH/STR/PAL/HPF/HY/MB/OLF/CTXsp and per cortical system, reading 'ratio' = nano/auto and zref; the value() helper is at L213)
- **Judgement:** The Python route pools per mouse by division, weighting by voxel count, which is the same idea, and it carries the nano/auto ratio. But Isocortex is split into 11 systems rather than kept as one division, auto on its own is not pooled, and there is no delta-z.
- **Second check:** v2_region_groups.py DIVISIONS L86-87 has 8 divisions (TH, STR, PAL, HPF, HY, MB, OLF, CTXsp). Isocortex is split into the 11 SYSTEMS of L61-77 and never pooled as one group. Per mouse, per group, L185-190 pools the voxel sums (s_rat[ids].sum()/c). That equals leaf means weighted by voxel count, the same idea as P10's roi_px_ero weights (L259-275), but uneroded. Only sig, ratio and sepratio are pooled. Auto alone is not, and there is no delta-z. value() is at L213 as cited.

### P10-08 · missing

- **What:** Division-level paired one-sided signed-rank test, nano > auto per division; optional Bonferroni over the divisions
- **Kind:** test
- **Old code:** P10_compare_nano_vs_auto.m:279-300
- **Python:** none (v2_region_groups.py:main L231-251 runs only young-vs-adult Mann-Whitney, with BH q)
- **Judgement:** The division-level nano-vs-auto question has no Python test. At the division level Bonferroni can reject (0.05/9 = 0.0056 > 9.8e-4). The old answer: HPF, CTXsp and STR nano > auto; isocortex nano 1.38 vs auto 2.39 in normalised units.
- **Second check:** v2_region_groups.py stats L230-247 are young vs adult only. group_stats.csv carries a descriptive adult_median of log2 ratio per division, with no test against a (centred) zero. Old answer verified in NanoVsAuto_macro_signrank_nosmooth_uncorrected.csv: HPF, CTXsp and STR p = 0.00098. Isocortex nano 1.377 vs auto 2.395. Division-level Bonferroni 0.05/9 = 0.0056 can reject.

### P10-09 · missing

- **What:** Division-level delta-z signed-rank test (pooled nano_z - auto_z > 0) per division
- **Kind:** test
- **Old code:** P10_compare_nano_vs_auto.m:279-300
- **Python:** none
- **Judgement:** Same gap as P10-05, at the coarser scale. The old answer: HPF and STR significant.
- **Second check:** No division-level contrast of nano against auto in any v2 file. Old values verified: HPF p_deltaz = 0.00098 and STR 0.042 significant. CTXsp 0.46, TH 0.65, the rest near 1.

### P10-10 · missing

- **What:** Table NanoVsAuto_macro_signrank_*.csv: per division the paired p/sig, delta-z p/sig, mean nano, mean auto and mean contrast
- **Kind:** table
- **Old code:** P10_compare_nano_vs_auto.m:302-314
- **Python:** none (group_stats.csv from v2_region_groups.py L253 holds only young-vs-adult statistics)
- **Judgement:** No division-level nano-vs-auto table exists in Python.
- **Second check:** group_stats.csv fields (v2_region_groups.py L242-247) are the n values, the medians, diff_median, diff_mean, Mann-Whitney and Welch p, diff_P20only, naive_minus_rws and the BH q. There is no nano-vs-auto p and no auto mean.

### P10-11 · partial

- **What:** Figure Region_NanoVsAuto_BarByMacro: one panel per division, grouped horizontal bars of nano and auto per STRU leaf (mean +- SEM), per-mouse dots, thin grey lines joining each mouse's nano and auto values, regions sorted, significant labels in orange bold
- **Kind:** figure
- **Old code:** P10_compare_nano_vs_auto.m:316-455; helper highlight_significant_labels L744-768
- **Python:** v2_region_plot.py:main figure L336-389 (region_plot.png, 'ratio' panel)
- **Judgement:** The Python figure shows per-mouse log2(nano/auto) dots for the 5 naive and 5 RWS adults (plus the young), but only for 35 hand-picked areas (AREAS L84-86), with medians and young-vs-adult stars. The full per-division structure list, the two channels shown separately with paired lines, mean +- SEM and nano-vs-auto significance are all lost.
- **Second check:** A weak partial. region_plot.png's 'ratio' panel (v2_region_plot.py L344-389) shows per-mouse log2 of the v2 ratio, only for the 35 AREAS of L84-86. Its zero line (L372) is exposure-dependent. Stars are young vs adult (L340, L364-371) and the bars are medians, not mean +- SEM. Neither channel is shown alone, with no paired lines, no nano-vs-auto significance and no full per-division leaf list. Within-mouse ordering across the 35 areas is still readable.

### P10-12 · missing

- **What:** Figure Region_NanoVsAuto_DeltaZ_BarByMacro: per STRU leaf, a bar of the mean nano_z - auto_z with SEM, per-mouse dots, a zero line, significant labels highlighted, one panel per division
- **Kind:** figure
- **Old code:** P10_compare_nano_vs_auto.m:457-532
- **Python:** none (region_plot.png's zref panel shows the nano zref only)
- **Judgement:** No Python figure shows an autofluorescence-referenced contrast per structure.
- **Second check:** The only signed per-structure panels in v2 are the nano zref panels (region_plot.png, group_plot.png) and the beyond_density residual bars (fig3_residual.png, v2_beyond_figures panel_c). Neither is referenced to autofluorescence alone.

### P10-13 · partial

- **What:** Figure Macro_NanoVsAuto_Paired: division-level grouped bars of nano and auto (mean +- SEM), per-mouse dots, paired lines, significant divisions highlighted
- **Kind:** figure
- **Old code:** P10_compare_nano_vs_auto.m:534-613
- **Python:** v2_region_groups.py:main dotplot L272-304, called at L306 (group_plot.png, 'ratio' panel)
- **Judgement:** group_plot.png shows per-mouse nano/auto per division and cortical system, with medians and young-vs-adult stars. It has no nano-vs-auto test, no separate channels and no single Isocortex group.
- **Second check:** The dotplot at v2_region_groups.py L272-303 is called for system/division keys at L305-309 (group_plot.png). Its ratio panel shows per-mouse log2 nano/auto per division and cortical system, with medians and young-vs-adult stars. No separate channels, no paired lines, no SEM, no single Isocortex group, no nano-vs-auto test.

### P10-14 · missing

- **What:** Figure Macro_NanoVsAuto_DeltaZ: division-level delta-z bars with SEM, per-mouse dots, zero line, significant divisions highlighted
- **Kind:** figure
- **Old code:** P10_compare_nano_vs_auto.m:615-662
- **Python:** none
- **Judgement:** Follows from P10-09 being missing.
- **Second check:** Follows from P10-09. No division-level contrast figure in v2.

### P10-15 · partial

- **What:** Voxelwise delta-z video Contrast_video_NanoMinusAuto_z: cohort z_nano - z_auto (each channel's P8 zscore_lr_sum), inside the intersection of the two brain masks, a diverging red-blue map with symmetric colour limits (95th percentile of |contrast|, floor 1.5), atlas boundaries, every 10 um coronal slice of the half brain
- **Kind:** video
- **Old code:** P10_compare_nano_vs_auto.m:664-735; redblue_cmap L770-782
- **Python:** v2_video.py:main L72-112 -> comparisons_v2/ccf/adult/video_ratio_adult.mp4 (cohort mean of voxelwise nano/auto + reliability t = mean/SEM, CCF 20 um, hemispheres folded, atlas outlines); also v2_compare.py:draw_figures L97-167 -> slices_ratio.png
- **Judgement:** Python shows where nano is high per unit autofluorescence, voxel by voxel, for the same 10 adults, with a reliability panel P10 lacked. What is lost is a signed contrast with a meaningful zero (equally enriched in both channels): the ratio's level depends on exposure and is drawn on hot 0-2. This is a cohort video, not a per-mouse one, so the 29 Sep drop does not cover it.
- **Second check:** Verified that data/comparisons_v2/ccf/adult/video_ratio_adult.mp4 exists. mice.txt lists exactly the 10 naive+RWS adults. v2_video.py draws ratio on hot_cut 0..2 (L43, L79, L95-96) with a t = mean/SEM panel (L93, L109-112), hemispheres folded (L50-52, L88-89). v2_compare.draw_figures L97-167 gives slices_ratio.png. video_zref_adult.mp4 is the nano-only analogue of P10's z_nano. There is no signed nano-z minus auto-z map with a meaningful zero. P10 L672-690 uses a diverging map, a mask intersection and a |contrast| p95 limit with a 1.5 floor.

### P10-16 · covered

- **What:** Diagnostic: the brain-wide mean and SD of each channel's cohort-mean LR-sum inside the analysis mask, printed and used as the z normalisers
- **Kind:** diagnostic
- **Old code:** P10_compare_nano_vs_auto.m:171-179
- **Python:** v2_sep_channel_check.py:main L131-133, L154-156 and figure panel L171-181 (per-adult dynamic range, p90-p10 of log2 structure means, for nano, auto and SEP)
- **Judgement:** Python reports how much each channel varies across the brain per mouse, robustly and in log units (nano 1.93 vs auto 1.07 log2), which answers the question behind P10's printed mu/sigma better. P10 only used them as normalisers.
- **Second check:** In P10 mu/sigma (L171-179) are printed and then used only as normalisers for P10-05. The question behind them, how much each channel varies across the brain, is answered per mouse in v2_sep_channel_check.py L131-133 (printed L154-156, figure panel L171-181). Verified values: 1.93 vs 1.07 log2.

### P10-17 · dropped

- **What:** Option agg_method: region value taken as P8's distance-weighted mean (power 4, the default) or the eroded mean (radius 3)
- **Kind:** option
- **Old code:** P10_compare_nano_vs_auto.m:32-33, 78-88
- **Python:** v2_region_plot.py:main L210-223 (plain mean over tissue voxels per structure at 20 um, MIN_VOX 250, no weighting or erosion)
- **Judgement:** The distance-weight vs eroded comparison was dropped on 29 Sep. The Python route uses one plain tissue-voxel mean on each brain's own atlas instead.
- **Second check:** The 29 Sep drop decision drops the 'distance-weight diagnostic'. That diagnostic is P8's Diagnostic_DistWeight_slice620_pow{1,2,4}.png (present in data/comparisons/merged_naive_rws_nano and _auto). P10's agg_method (L32-33, L79-88) only chooses which P8 matrix to read (per_mouse_mean_dw or per_mouse_mean_ero). v2's plain tissue-voxel mean (v2_region_plot.py L210-223, MIN_VOX 250) applies to every reading, so nothing P10-specific is lost. Note: v2 has no border guard at all, so spill from HPF (nano about 7x brain level) into neighbouring leaves is unguarded in both channels.

### P10-18 · file handling

- **What:** Display options: sort_regions_by ('nano' | 'auto' | 'diff'), show_per_mouse_dots, show_paired_lines, the colour palette
- **Kind:** option
- **Old code:** P10_compare_nano_vs_auto.m:39-42, 53-59, 353-359
- **Python:** n/a
- **Judgement:** Presentation choices only. If the figure is rebuilt, sorting by the centred nano-minus-auto contrast is the useful order.
- **Second check:** L39-42, L53-59 and L353-359 are presentation only. Sorting by a centred contrast is a sensible default for a rebuilt figure.

### P10-19 · file handling

- **What:** Plumbing: loading P8's per-mouse and cohort caches for both channels, checking that the mouse lists and ROI acronyms match (pairing), loading the sparse ROI mask cache, filename suffixes (_uncorrected/_bonferroni), .fig/.png saves
- **Kind:** plumbing
- **Old code:** P10_compare_nano_vs_auto.m:61-106, 162-169, 225-238, 448-454
- **Python:** v2_per_mouse.py (nano 'sig' and 'auto' saved together in one per-mouse npz, so they are paired by construction)
- **Judgement:** File I/O and consistency checks. In v2 the pairing check is unnecessary because both channels come from the same brain's registered tiffs in one file.
- **Second check:** v2_per_mouse.py L200 saves sig and auto together per brain, so the pairing check (P10 L90-100) is unnecessary. A side finding that argues for porting rather than re-running: P10's current inputs are no longer mutually consistent. per_mouse_region_means_merged_naive_rws_nano_nosmooth.mat is dated 2026-05-02. mean_lr_sum_merged_naive_rws_nano_nosmooth.mat (source of mu_n/sg_n and z_nano) was rewritten 2026-09-04, after naive/nano_4d_normalized.mat was redone that day. The auto caches date from 2026-05-01 and come from auto_4d.mat of 2025-11-19. Re-running P10 today would mix per-mouse values and normalisers from different runs.

### P10-20 · covered

- **What:** Scope: cohort = merged naive+rws adults (10 mice, behavior excluded); region set = 228 STRU leaves under 9 DIVI divisions, AP slices 100-700 only, left+right hemisphere sum of P6bis-normalised channels at 10 um
- **Kind:** scope
- **Old code:** P10_compare_nano_vs_auto.m:21-26, 44; inherited from P8_characterize_merged_distribution.m:83-102, 263-275
- **Python:** v2_cohort.py (NAIVE + RWS, the same 10 adults); v2_region_plot.py:main L186-223 (every structure with >= 250 20-um tissue voxels in that brain, about 259 per adult)
- **Judgement:** Same 10 adults, and behavior is absent in both routes. Python's structure set is broader (it adds P, MY, CB, fibre tracts and ventricles and has no AP restriction), it is background-subtracted rather than affine-normalised with abs(), and it pools both hemispheres at 20 um. A nano-vs-auto test added to Python should restrict itself to grey matter, as v2_beyond_density.keep_structure (L177) does.
- **Second check:** Same 10 adults: get_cohort.m L57-66 and v2_cohort.py L109-110; data/comparisons_v2/ccf/adult/mice.txt lists those 10, and behavior is absent in both routes. Structure set checked: in the 9 P10 divisions v2 keeps, for CGF027, TH 46, HY 45, Isocortex 42, MB 34, HPF 15, STR 15, OLF 11, PAL 10, CTXsp 8. P10 had TH 45, HY 44, Iso 43, MB 40, HPF 15, STR 14, OLF 11, PAL 9, CTXsp 7 (228 rows, 206 tested). So P10's set is essentially contained, give or take a few MB leaves under MIN_VOX. v2 adds P, CB, fibre tracts, ventricles and unassigned labels (237-280 structures per adult). For the port, restrict to v2_beyond_density.keep_structure (L177-188) or to the 9 divisions.
