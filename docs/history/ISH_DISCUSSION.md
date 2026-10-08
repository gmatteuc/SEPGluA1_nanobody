<!--
Kept as written on 8 October 2026, when Giulio agreed to its section 8 ("all as
recommended"): the principle, step 0, analyses 1 to 5, the keep/merge/drop table and
the five decisions that docs/ISH_ANALYSIS.md builds on. Its figures, scripts and
read-only checks were in a working folder outside the repository and are not kept
here; the figure embeds are given as their file names. The numbers are those of the
run of 5 October and of read-only checks; ISH_ANALYSIS.md has the numbers of the
analysis built from it.
-->

# The ISH line, from the beginning

Written for a discussion of the whole ISH line: comparing the adult nano map with Allen in situ hybridisation maps. Everything here was read or recomputed without changing the repository or the data. Numbers come from the outputs of the 5 October rerun unless a section says otherwise.

- Figures are in `figures\`. Copies of the production figures are in `figures\v2\`.
- The scripts that drew figures 1 to 6 are in `scripts\`, and the read-only checks are in `scratch\` and `p9_original\`.
- Sections 1, 2 and 8 carry the discussion. Sections 3 to 7 hold the evidence.

## Words used here

| word | meaning |
|---|---|
| P9 | the MATLAB comparison of April 2026, now `adult_matlab/run_compare_with_allen_ish.m` (section 2). |
| structure | one atlas region, such as CA1 or VPM. Every comparison here runs across structures, with one value per structure on each side. |
| ρ (rho) | the Spearman correlation across structures: do the two maps put the structures in the same order? |
| division | a coarse part of the brain: Isocortex, OLF, HPF, CTXsp, STR, PAL, TH, HY, MB, and in Python also pons, medulla and cerebellum. |
| `zref` | the Python reading of the nano map. For each mouse it takes log2 of each structure over the isocortex mean, centres it on the brain's median structure and scales it by its p90–p10 spread. The cohort map is the mean of the 10 adults' `zref`. |
| partial ρ | ρ after a third map, for example Gria1, has been taken out of both sides, on ranks. |
| arm | the ratio of two imaging channels of the same brain, per structure: SEP / auto, nano / auto or nano / SEP (section 5). The name comes from a planned three-arm design. |
| beyond | short for "beyond abundance and density" (section 6): how much of the map receptor mRNA and synaptic density predict, and what is left. |
| spatial null | random maps with the same smoothness as the brain map. They show how large a ρ comes about by chance. Not built yet (A7). |
| ceiling | how well two halves of the cohort agree. No explanation can reach higher. |
| leftover | the part of the map that a fitted model does not predict. |

## 1. The question

The nano map gives, for each structure, how much surface GluA1 the nanobody finds. An Allen ISH map gives, for each structure, how much of one gene's mRNA there is. Comparing the two asks whether they put the structures in the same order.

**What the comparison can tell:**
- **Whether the map looks like a gene's expression pattern, and which genes it looks like most.** A rank correlation ignores any scale factor per brain or per session. So it works where absolute numbers do not.
- **Whether a chosen set of expression maps accounts for the map,** and how much is left over (section 6).
- **Whether one kind of gene matches better than another.** This needs two things: the kinds are fixed before looking, and the test allows for co-expressed genes not being independent.

**What it cannot tell:**
- **Surface against total receptor.** Surface and total receptor both follow the postsynaptic side of synapses, so both predict a map that looks like the postsynaptic genes. No ranking of genes separates them. Only a total-receptor measurement in the same brains does.
- **Protein, from mRNA.** mRNA sits in cell bodies. The receptor sits on dendrites, which can lie in a different structure; cortical layer 1 is the plain case. Translation and turnover also differ between regions.
- **Significance, on its own.** Brain maps share smooth gradients: cortex and hippocampus are high, thalamus and hypothalamus low. A map that knows nothing but each structure's division still correlates with Cacng8 at 0.78 (Figure 5A). Without a spatial null a ρ has no p value, and an ANOVA over genes is too optimistic (`docs/ROADMAP.md` §4).
- **Anything about the same animals.** Each Allen ISH experiment is one P56 mouse, on a 200 µm grid. Our adults are older, and their exact age is not recorded (`docs/ROADMAP.md` §12).

So the ISH line can describe the map and test how far expression predicts it. It cannot validate the surface measurement.

*(figure not kept: `figures/fig1_one_gene_comparison.png`)*

**Figure 1. What one gene comparison is.** Rows: the nano map, Cacng8, Gria1, all on CCF plane 700. Columns:
- the map as measured (for nano, each brain over its own isocortex mean, for display);
- the same map as one rank per structure (for nano, from `zref`);
- the scatter of ranks, with ρ.

Cacng8 gives +0.80 on 235 structures and is first of 95 genes. Gria1 gives +0.67 on 225 structures and is 9th. Hollow dots are structures that fewer than 10 adults reach. AOB, seen in 4 mice, is the clearest outlier.

## 2. The original analysis (P9) and its headline

- **Script:** `adult_matlab/run_compare_with_allen_ish.m`. At the tag `refactor-start` it was `P9_compare_nano_vs_allen_ish.m`; the computation is unchanged.
- **Run:** 21–22 April 2026.
- **Outputs:** `S:\ElboustaniLab\#SHARE\Processed\SEP-GluA1_Project\comparisons\comparisons\merged_naive_rws_vs_ish_summary_nosmooth\`, plus one folder per gene. `D:\sep_histology\data\comparisons\` holds an identical copy.

### 2.1 What was computed

1. **The nano map.** 10 adults, naive and RWS pooled: CGF027, CGF028, CGF033, CGF034, CGF035, MG691, MG692, MG693, MG736 and MG737.
   - Each mouse was normalised by a line fitted onto its group's median isocortex.
   - RWS was aligned onto naive by a second line.
   - Left plus mirrored right, then the mean over mice. No smoothing.
   - The whole-brain z-score is for display only. It changes no correlation.
2. **The ISH maps.** One coronal Allen experiment per gene: expression energy on the 200 µm grid.
   - Dim sections were replaced by interpolation. A section counts as dim if its median is below 0.3 × the median section.
   - The grid was stretched onto the 10 µm atlas, negatives clamped to 0, then left plus mirrored right.
3. **The structures.** The 228 finest structures of nine divisions (Isocortex, OLF, HPF, CTXsp, STR, PAL, TH, HY, MB).
   - 206 of them lie inside atlas planes 279–879. These 206 are used for every gene.
   - TH, HY and MB make up 115 of the 206.
   - There is no hindbrain, no cerebellum and no fibre tract.
4. **Two means per structure:**
   - eroded by a 30 µm sphere;
   - distance-weighted by (d/dmax)^4, so the border counts for nothing.
5. **Per gene:** Spearman and Pearson on both means, plus a voxel Pearson. Genes are ranked by the distance-weighted Spearman. No gene has a p value.
6. **Groups.** Each gene carries a hand-written category from `gene_targets.csv`.
   - The figure code splits "auxiliary" into **Aux forebrain** (Cacng8, Cacng3, Cnih2, Cnih3, Grm5) and **Aux other**.
   - A one-way ANOVA across the ten groups goes in the title.

### 2.2 The headline Sami saw and approved

*(figure not kept: `figures/p9_violin_spearman_ero.png`)*

**`violin_spearman_ero.png`, 22 April.** The eroded Spearman of 97 genes, one violin per group, sorted by median. ANOVA p = 0.032 (F = 2.2, df 9 and 87).

| group | n | median ρ | highest | lowest |
|---|---|---|---|---|
| Aux forebrain | 5 | 0.63 | Cacng8 0.80 | Cnih3 0.17 |
| Plasticity | 5 | 0.49 | Arc 0.65 | Fos 0.16 |
| AMPAR core | 4 | 0.46 | Gria3 0.58 (Gria1 0.51) | Gria4 −0.19 |
| Scaffold | 9 | 0.42 | Dlg2 0.78 | Grip2 −0.07 |
| Excitatory | 16 | 0.35 | Nrgn 0.69 | Slc17a6 −0.51 |
| Trafficking | 11 | 0.25 | Nptx1 0.70 | Syt3 0.00 |
| Ctrl struct | 10 | 0.21 | Hcn1 0.64 | Calb2 −0.53 |
| Ctrl inhib | 22 | 0.14 | Htr3a 0.70 | Glra1 −0.50 |
| Aux other | 8 | 0.12 | Grm3 0.41 | Cacng5 −0.12 |
| Ctrl glia | 7 | 0.00 | Cx3cr1 0.34 | Olig2 −0.21 |

Gene by gene: Cacng8 is first (0.80), Dlg2 second, and Gria1 23rd of 97 (0.51). The reading at the time:
- the forebrain AMPA receptor auxiliary subunits match the map best;
- glia match it least;
- the categories differ.

### 2.3 What holds

- **Cacng8 is first.** It is first under all five P9 metrics, and on every structure set tried in Python (0.77–0.82). Only the `sepratio` reading puts Cnih2 a hair ahead.
- **The gene order reproduces.** The Python route orders the same genes as P9 at ρ 0.91 (Figure 2B).
- **The category pattern reproduces, with the same split.** On today's numbers the ANOVA gives p 0.006 to 0.030, depending on the structure set (Figure 2C).
- **Glia sit low and postsynaptic genes high.** The Python word test, which uses annotation nobody here wrote, points the same way. Genes tagged postsynaptic sit +0.34 above the rest; presynaptic ones sit +0.04 above (section 3.3).
- **It is easy to read.** One number per gene, and every gene is named.

### 2.4 What does not hold

1. **The p rests on the split.**
   - Without the Aux forebrain / Aux other split, p = 0.20 (Kruskal-Wallis 0.19); on today's numbers, 0.046. The split was written into the figure code by hand, by forebrain expression.
   - Dropping the two stretched grids (Olig2, Calb2) gives 0.053.
   - An earlier run on the first 50 genes gave 6.6e-4. The p moved with the panel (`docs/REFACTOR_COVERAGE.md`, P9-26, P9-28).
2. **The ANOVA treats genes as independent.** Co-expressed genes are not independent draws, so the p is too small even with the split (`docs/ROADMAP.md` §4).
3. **There is no null.** No gene has a p value, and nothing says how large a ρ two unrelated smooth maps would give (A7).
4. **The categories are mixed, and were written while the genes were being chosen** (`docs/adult_ish_design.md`, l.145–150):
   - "auxiliary" holds the metabotropic receptors Grm1–5 beside the TARPs, and Grm5 lands in Aux forebrain;
   - 9 of the 11 "trafficking" genes are presynaptic vesicle machinery; only Nsf and Nptx1 act on AMPA receptors;
   - "Ctrl inhib" holds neuromodulator receptors, and Htr3a, Chrm1 and Drd1 sit above Gria1;
   - Cacng2 (stargazin) is not in the panel.
5. **The ranking is mostly a contrast between divisions.**
   - Give every structure only its division's median nano value. Cacng8 still gets 0.74 (distance-weighted means), and the gene order agrees with the real one at 0.97.
   - Inside divisions the order changes. Dlg2, Cacng8, Gria1 and Grm5 stay near the top. Nrgn, Camk2a and Slc17a7 fall to ranks 37, 43 and 64.
   - Sources: `p9_original\p9_within_division_readonly.csv`; Figure 5 shows the same on today's numbers.
6. **The result depends on the statistic.**
   - Gria1 is 23rd by eroded Spearman, 22nd by distance-weighted Spearman, 9th by Pearson and 6th by voxel Pearson.
   - Pearson follows the nine hippocampal and septal structures that sit far above the rest (z 1.6 to 3.9; the next is 0.5).
   - Erosion against distance weighting changes nothing (0.997). A 30 µm erosion is small next to a 200 µm grid.
7. **The ISH handling has defects** (`docs/ROADMAP.md` §4; REFACTOR_COVERAGE P9-05, P9-07, P9-X3):
   - the stretched grid has a 1.5–2.5% scale error;
   - Allen's missing-data value −1 was interpolated and clamped to 0, so partly missing sections count as no expression;
   - Olig2 and Calb2 sit in a box of their own and were stretched;
   - Slc17a6's true absence in sections 15–18 was "repaired".
8. **The nano map has defects.**
   - `abs()` turned unreached planes into bright slabs.
   - The per-mouse intercept breaks ratios between regions.
   - RWS was aligned onto naive by a line.
   - The tissue mask came from the nano channel itself.
9. **No autofluorescence control.** The `channel='auto'` option was never run (P9-03).
10. **The figure.** The y-axis stops at −0.4 and cuts off Th, Glra1, Slc17a6 and Calb2. Groups of 4 or 5 genes are drawn as smooth densities.

After seeing the figure, Sami asked (28 and 30 April) for:
- a check of ISH section quality (A2);
- a permutation null (A7);
- a gene table (A9);
- the whole chain run on autofluorescence (A8) and on DAPI;
- Gene Ontology groups instead of hand categories.

P9 did none of them.

*(figure not kept: `figures/fig2_p9_headline_today.png`)*

**Figure 2. The headline, recomputed on today's numbers.**
- A: each group's P9 violin (grey) beside today's (orange), with a line joining each gene's two values.
- B: P9 against today, one dot per gene. The order agrees at 0.91.
- C: the ANOVA p under each choice it depends on.

## 3. What the Python route does now

- **Code:** `mapping/sepmap/ish/` and `mapping/sepmap/adult/`. The run order is in `mapping/README.md`, steps 11–26.
- **Outputs:** `D:\sep_histology\data\adult_v2\{ish,panel,arms,beyond}\`.
- **Two gene panels:** the 100 hand-picked genes of P9 (95 usable), and 390 genes built from Gene Ontology.

### 3.1 ISH onto the atlas (`regions.py`)
- **Does:**
  - samples the atlas onto Allen's 200 µm grid, instead of stretching the grid onto the atlas;
  - −1 becomes missing, and grids in a non-reference box are dropped;
  - one mean per gene and structure, full and eroded by one ISH voxel.
- **Figure:** none.
- **Result:** 95 of 100 genes. Sst, Chat and Tph2 return 404, and Olig2 and Calb2 are off-grid.

### 3.2 Each gene against the map (`compare.py`)
- **Does:**
  - Spearman per gene, between the 10-adult mean of each structure and the gene's structure means;
  - on the structures both sides have, 213–237 per gene;
  - five nano readings, shown side by side: `zref`; `cref` (each brain over its own isocortex mean); `subref` (over its own subcortex mean); `ratio` (nano / auto); `sepratio` (nano / SEP).
- **Figure:** Figure 1 (section 1) shows one comparison. The production figure, `ish_old_vs_new.png`, compares every gene's ρ with its P9 value; Figure 2B is the plainer version of it.
- **Result:**
  - The gene order agrees with P9 at 0.91.
  - Cacng8 is first (+0.80). Gria1 is 9th (+0.67) on all structures, 17th on structures seen in all 10 adults, and 18th on P9's nine divisions. The docs still quote 10th, a number from before 5 October.
  - `zref`, `cref` and `subref` keep each mouse's order, so they give nearly the same ranking (Gria1 9th, 11th, 9th). `ratio` and `sepratio` differ more: 0.87 and 0.84 against P9.

### 3.3 What kind of gene is at the top (`words.py`)
- **Does:**
  - tags each gene with its GO terms, and with single words taken from the terms;
  - for every tag carried by at least 5 genes, compares the ρ of genes with the tag against genes without it;
  - statistics: Mann–Whitney, a bootstrap interval, and Benjamini–Hochberg correction;
  - nine words were named in advance.

*(figure not kept: `figures/v2/ish_word_enrichment.png`)*

- **Result:** the contrast named in advance holds.
  - Postsynaptic genes: +0.34 [+0.18, +0.52], q 0.008.
  - Presynaptic genes: +0.04 [−0.10, +0.31], q 0.75.
  - In the first two panels, 13 of the 24 top bars have intervals that cross zero, nearly all of them groups of 5 to 8 genes. The large groups hold: glutamatergic, postsynaptic and density, q 0.008.

### 3.4 Subunits against localisation genes (`roles.py`)
- **Does:**
  - gives each gene a role by protein function, assigned by hand after the ranking had been seen;
  - builds one rank composite of the 4 subunits and one of 15 localisation genes;
  - measures how much of the map each composite explains alone and how much they share;
  - uses as the null every other way of splitting the same 19 genes 4 against 15.

*(figure not kept: `figures/v2/ish_roles.png`)*

- **Result:** localisation adds more than the subunits: unique R² 0.114 against 0.014. But any split of the same 19 genes does as well (p = 0.49).
  - The test could only have detected an effect about four times larger than the one seen.
  - Median ρ: subunits +0.58, localisation +0.39.

### 3.5 A panel nobody here chose (`panel_build.py`, `panel_fetch.py`)
- **Does:**
  - takes genes from GO terms in three roles: subunits; localisation (AMPA receptor transport, anchoring, clustering and auxiliary proteins); and postsynaptic-density controls;
  - keeps every Allen experiment of each gene, coronal and sagittal.
- **Figure:** none.
- **Result:** 390 genes and 669 experiments: 4 subunits, 84 localisation, 300 controls, 2 delta receptors. 667 experiments downloaded.

### 3.6 How reliable one Allen map is (`reliability.py`)
- **Does:**
  - for each gene with two or more experiments, the median ρ between its experiments;
  - merges each gene's experiments into one profile.

*(figure not kept: `figures/v2/ish_reliability.png`)*

- **Result:** 218 genes can be measured.
  - Median 0.69 (quartiles 0.52–0.77); 22 genes are below 0.3.
  - Gria1 0.91, Cacng8 0.91, Dlg2 0.91, Gria4 0.57, Cnih2 0.08.
  - A gene's ρ with the map is capped by its reliability.
  - Reliability rises only weakly with expression (ρ +0.18). The middle panel's title claims more than that.

### 3.7 The powered test (`panel_test.py`)
- **Does:**
  - removes the subunit composite from the map;
  - asks whether the 84 localisation genes predict what is left better than 84 control genes matched on expression;
  - permutes the labels 20,000 times;
  - adds a positive control whose effect is known to exist: control genes split by their reliability.

*(figure not kept: `figures/v2/ish_panel_test.png`)*

- **Result:** no difference: −0.000, p 0.98. A difference of 0.063 would have been detected. The docs still quote p 0.74, from before 5 October.
  - The positive control is found (+0.16, p 0.0004), so the negative is informative.
  - By partial ρ, the three genes after Cacng8 are controls: Arpc5, Cdk5r1 and Ptk2b.

### 3.8 Arms, and 3.9 beyond abundance
These are in sections 5 and 6.

## 4. P9 against Python

| | P9 | Python | why it changed | what got less clear |
|---|---|---|---|---|
| nano map | voxel mean of 10 adults, each normalised by an affine fit | per-mouse structure means in five readings, then the mean of 10 adults | an affine fit per mouse breaks ratios between regions (`adult_ish_design.md`, D3) | five readings side by side look like five results; three of them give nearly the same ranks |
| ISH onto the atlas | the grid stretched; −1 clamped to 0 | the atlas sampled onto the grid; −1 set to missing | correctness (P9-09, P9-X3) | nothing |
| section quality | dim sections interpolated | no check per section, only reliability per gene | the repair filled true absences | no quality view per gene. The 95-gene table is unprotected; A2 is pending, and an unreviewed scan flags 14 of 95 genes (Gria1 is clean) |
| structures | 206 structures of nine divisions, planes 279–879 | every structure in the per-mouse table, with hindbrain, fibre tracts and "unassigned" labels, and no minimum number of mice | inherited, not decided; S1, A1 and A3 fix it | Gria1 moves from 9th to 17th to 18th with the set, for reasons unrelated to the method |
| borders | distance-weighted mean as the primary; eroded mean in the headline; both sides treated | the full ISH mean; the eroded one is written but never read; nano never eroded | full and eroded agree at 0.998 | no Python figure uses the statistic Sami saw. Eroding the ISH side moves Gria1 from 18th to 11th on P9's nine divisions |
| statistic | Spearman, Pearson and voxel Pearson | Spearman only | expression energy has no common scale; the voxel correlation was dropped on 29 Sep | the Spearman/Pearson sensitivity (Gria1 22nd against 9th) is shown nowhere (A3) |
| panel | 100 hand-picked coronal genes | the same 95, plus the 390-gene GO panel | the 4-against-15 test lacked power | two panels and two tables: the headline genes are in one, the decisive test in the other |
| groups | 9 hand categories plus the hand split | GO terms and words; roles by function; GO roles. The "machinery" of compare and arms is still P9's categories | categories written while choosing the genes make the explanation circular; Sami asked for GO groups | four vocabularies for one question, and none matches the figure Sami approved |
| tests | one ANOVA | Mann–Whitney with BH; an exact split permutation; a matched label permutation with a positive control; a paired Wilcoxon | the ANOVA treats co-expressed genes as independent | four tests with four nulls, none of them spatial. A reader cannot tell which one answers "which kinds of genes match" |
| per-gene views | scatter coloured by division, fits within divisions, paired bars, a three-panel video | none (A6, A10) | not yet built | no picture of what 0.80 against 0.67 looks like |

The data did not change the answer. The gene order agrees at 0.91, and the category pattern survives with the split. What changed is how the result is read. The route also split the line into pieces that nobody has put back together. Section 8 proposes how.

## 5. The arms

### What an arm is

Each brain was imaged in three channels:
- **nano:** the nanobody against the SEP tag. The sections were not permeabilised, so it should reach only tag on the membrane.
- **SEP:** the tag's own green fluorescence. It was meant to show all tagged receptor, inside the cell and on the membrane.
- **auto:** a channel with no label, which records the tissue's autofluorescence.

An **arm** is the ratio of two of these channels, per structure and per adult, in log2 (`mapping/sepmap/adult/arms.py`, which writes `adult_v2\arms\region_means_arms.csv`). There are three:

| arm | ratio | meant to stand for | ρ with Gria1 ISH |
|---|---|---|---|
| `sepauto` | SEP / auto | total receptor | −0.09 |
| `ratio` | nano / auto | surface receptor | +0.53 |
| `sepratio` | nano / SEP | the fraction at the surface | +0.61 |

The name comes from a three-arm design. The prediction was written down before the run:
- Gria1 mRNA should track total receptor (`sepauto`) best.
- The genes that bring receptors to the membrane and hold them there should track the surface fraction (`sepratio`) best.

If both held, the surface reading would be validated in one figure, with no new animals.

*(figure not kept: `figures/fig3_arms.png`)*

**Figure 3. The arms.**
- A: the three channels, and the three ratios between them.
- B: the three raw channels of one adult on plane 700, each on its own scale.
- C: which signal each channel tracks, one dot per adult and a line joining each mouse.

### What happened

- **Test 1: each gene against each arm.** Gria1 went the wrong way: −0.09 on total, +0.61 on the fraction. Its swing towards the fraction (+0.70) is larger than that of 30 of the 33 "machinery" genes.
- **The check that explains it** (`adult/sep_channel_check.py`, Figure 3B and C). It uses the raw channels, with no ratios.
  - SEP tracks autofluorescence at 0.79 ± 0.04, and between 0.70 and 0.82 in every adult. Nano tracks it at 0.26.
  - SEP varies across the brain less than autofluorescence does. The p90 − p10 spread of log2 is 0.95 for SEP, 1.08 for autofluorescence and 1.93 for nano.
  - Against Gria1: nano +0.61, SEP +0.30, auto +0.24. SEP with its autofluorescence part removed gives +0.12.
- **So in this fixed, mounted tissue the green channel is mostly autofluorescence.**
  - `sepauto` is close to noise.
  - `sepratio` behaves like `ratio`; they agree at 0.89–0.97 within each mouse.
  - Test 1 measures how much a gene looks like nano. It says nothing about surface against total.

### Why the old figure reads as obscure

`arms_vs_genes.png`:
- Its axes still say "total receptor", "surface receptor" and "surface fraction". The module's own docstring disproves those meanings.
- Two of its three panels show Test 1, which is void, and the title does not say so.
- The third panel (Test 2) shows only the 33 "machinery" genes, with no comparison group. And "machinery" is P9's hand categories.

### What survives

1. **The finding about the green channel itself.** It is a clean methods result: SEP fluorescence does not survive this protocol well enough to report receptor. It belongs in the methods or a supplement, as Figure 3B and C.
2. **Test 2, weakly.** With Gria1 mRNA taken out of the nano map (`ratio` arm), the 33 machinery genes keep a median partial ρ of +0.12. 20 of 33 are positive, Cacng8 +0.54 and Cnih2 +0.50 among them. But:
   - the effect comes from the auxiliary group alone (+0.23, 12 of 13 positive), and that group includes Grm1–5. Trafficking gives −0.08 and scaffold −0.01;
   - the excitatory markers (+0.10) and two inhibitory controls (Htr3a +0.45, Drd1 +0.39) survive just as well;
   - over all 94 genes the median is +0.03;
   - the powered test (section 3.7) asks the same question properly and finds no difference.

   The panel test replaces Test 2.
3. **An open point for the microscopy side.** With its autofluorescence part removed, the green channel still tracks nano at 0.42, but Gria1 at only 0.12. That is either tag that survived or nano fluorescence leaking into the green channel. The filter sets decide which.

### What would answer the question

A total-receptor measurement in the same brains: a total-GluA1 antibody after permeabilisation, or autoradiography, on some of them. With that, the three-arm design can run as planned. No reanalysis can replace it.

## 6. Beyond abundance and density

### The question

Is the map just how much receptor mRNA a structure has? Or just how many synapses it has? Ranking genes cannot tell, because a map of pure synaptic density would also put the postsynaptic genes on top. So the question is turned round:
- predict the map as well as possible from abundance and density;
- ask whether what is left is real.

### What is done

Code: `adult/beyond_density.py`, `beyond_controls.py`, `beyond_figures.py` and `beyond_regression.py`.

1. **The structures are set by rule, before any fit:** grey matter, measured in all 10 adults and by every covariate gene. 126 are kept: Isocortex 33, TH 24, HY 15, STR 12, OLF 9, MB 9, HPF 9, CTXsp 7, PAL 6, P 2.
2. **The map** is the 10 adults' mean `zref` per structure, as ranks.
3. **Four explanations,** each as ranks:
   - abundance: the mean rank of Gria1–4;
   - markers: the mean rank of 11 synaptic marker genes;
   - PSD: the first principal component of 188 postsynaptic-density genes;
   - the autofluorescence of the same brains.

   Each may bend: x, x² and x³.
4. **Cross-validated R².** Each model is fitted on four fifths of the structures and scored on the remaining fifth.
5. **The ceiling.** The 10 adults are split into two fives, in all 126 ways. The two half-maps agree at 0.974 on average. Spearman-Brown turns that into 0.987 for all ten, and the code uses its square, 0.974, as the explainable share (item 5 of "What needs rewording").
6. **Unexplained** = 1 − R² / ceiling.
7. **Replication.** The model is fitted to each half-map separately, and the two leftovers are correlated.
8. **Seven controls, A to G:** spatial gradient, structure size, single animals, naive against RWS, curvature, the whole gene panel, and the reading.

*(figure not kept: `figures/fig4_beyond.png`)*

**Figure 4. Beyond abundance and density.**
- A: the variance budget, as a share of what can be explained.
- B: half against half, for the map and for its leftover.
- C: the same pipeline run on maps whose answer is known.

### The numbers

Sources: `adult_v2\beyond\variance_partition.csv`, `controls.csv` and `for_sami\numbers_for_the_caption.txt`.

| model | cross-validated R² | share of the ceiling |
|---|---|---|
| abundance (Gria1–4) | 0.267 | 27% |
| markers | 0.162 | 17% |
| PSD | 0.268 | 28% |
| autofluorescence | −0.064 | 0% |
| all four, straight | 0.409 | 42% |
| all four, bent | 0.624 | 64% |

The single rows are fitted straight, as in `variance_partition.csv`. Figure 4A lets every step bend, so its abundance bar reads 30%, not 27%.

- **36% of the explainable variance is left over.** A bootstrap over structures puts it between about 27% and 53%.
- **The leftover replicates across animals at 0.936** (0.887–0.951).
- **All seven controls pass.**
- **The largest leftovers:**
  - above prediction: the medial geniculate, the subthalamic nucleus, the lateral habenula, the ventral LGN, the septofimbrial nucleus and CA3;
  - below prediction: VPM, VPL, RSPd and the posterior complex.
- **Production figures:** `A_what_explains.png`, `E_regression.png`, `F_maps.png` (not kept).

### What holds

- **The map is highly reproducible.** Half-cohort maps agree at 0.974 on 126 grey-matter structures.
- **Autofluorescence of the same tissue predicts none of it.** This is a clean specificity point.
- **The leftover is not a confound.** It is not one animal, the whisker manipulation, structure size or a smooth gradient.
- **The direction holds.** Part of the map is not predicted by these covariates.

### What needs rewording

New read-only checks, in `scratch\arms_beyond\`:

1. **"The leftover replicates" is not separate evidence.**
   - The map's reliability and the in-sample fit alone predict a leftover replication of 0.922. The observed value is 0.936.
   - Once the map replicates, its leftover must too. The noise null rejects something the ceiling has already ruled out (Figure 4B).
2. **The abundance term dilutes Gria1.**
   - Gria4 tracks the map at −0.12, so averaging it into the composite weakens the abundance predictor.
   - With Gria1 alone as the abundance term, abundance explains 47% and the full model leaves 27% (Figure 4A, second bar). With the four subunits as separate predictors, 26% is left.
   - VPM, VPL and the posterior complex then leave the bottom of the list. There nano is as low as Gria1 mRNA, but Gria3 and Gria4 are high (`scratch\arms_beyond\abundance_variants.log`).
3. **Part of the leftover is one Allen map disagreeing with another.** Maps with a known answer were run through the same pipeline, with nano's animal noise and covariates from the other Allen experiments (Figure 4C).
   - A map that is exactly abundance plus density leaves 6–8%, and that leftover replicates at 0.66–0.73.
   - A map that is Gria1 mRNA leaves 27–31%, replicating at 0.88–0.89.
   - So replication cannot tell biology from ISH mismatch.
4. **"Beyond abundance and density" is not "beyond gene expression".** 17 components of 253 panel genes reach 84% of the ceiling. Only 16% escapes the whole panel.
5. **The ceiling is squared once too often.** Spearman-Brown gives 0.987, which is already a share of variance. The code squares it, to 0.974. The effect is under one point (36.0% becomes 36.8%), but "97.4% explainable" should read 98.7%.
6. **Some numbers are stale.**
   - `beyond_density.py`'s docstring still says 61% / 39% / 0.97.
   - `beyond_controls.py` says "about 40%".
   - `docs/adult_ish_design.md` has the numbers of 26 September.

**What the leftover is: unknown.** Nothing here says it is the surface fraction. Translation, turnover, subunit composition or nanobody access would all land in the same place (`docs/adult_ish_design.md`, l.517–521). Only a total-receptor channel separates them.

**A statement that holds today:** about a quarter to a third of the map's reproducible pattern is not predicted by receptor mRNA or synaptic density.
- At least 6–8 points of that would appear even for a map made of nothing else, because Allen maps disagree with each other.
- The whole gene panel, with no hypothesis, leaves only 16% of the map unpredicted.

## 7. The open questions

The sources are `docs/REFACTOR_COVERAGE.md` (A1–A10), `docs/ROADMAP.md` §3, §7 and §12, `docs/history/REFACTOR_PLAN.md` (S1–S5) and `docs/adult_ish_design.md`. None of A1–A10 is built.

**Inputs that every ISH analysis uses**
- **A1:** one declared structure set and `zref` reference, read by every table. It is first in line.
- **S1:** that set is the grey matter seen in all 10 adults. Agreed on 30 Sep, not applied yet.
  - On today's table: 280 structures, 235 seen in all 10 adults, 204 of those grey matter. The plan quotes 234, from before the 5 October rerun (`scratch\s1_count.log`).
  - The rule settles the two questions `REFACTOR_COVERAGE.md` left open: all 10 adults, not 8; hindbrain and cerebellum only where all 10 reach them, which today means 2 pons structures.
  - 202 of the 204 are structures P9 used, so S1 is close to P9's own set. On it, Cacng8 is first (+0.81, 171 structures with enough ISH voxels) and Gria1 17th (+0.60, 161 structures). Held under S5 until A2 and A3 have run.
- **A2:** ISH section quality. Failed sections are set missing, never interpolated, with a reviewed list of true absences. It also touches the beyond analysis, whose covariates come from the 390-gene table.
- **A9:** a gene documentation table (Sami, 28 Apr).
- **Panel repair:** Sst, Chat and Tph2 return 404. Olig2 and Calb2 are off-grid, and other experiments exist for both.
- **The adults' age:** not recorded. The Allen mice are P56.

**The gene ranking**
- **A3, the ranking:** on A1's set, through one adult profile with a minimum-mice rule, with a sensitivity table under three structure sets.
- **A3, its robustness:** Spearman against Pearson on log2, and full against eroded means on both sides.
- **A6:** agreement within divisions, and per-gene scatters. Figure 5 is a preview.
- **A7:** a spatial null for each gene's ρ and for ρ(Cacng8) − ρ(Gria1). Sami asked for a permutation null on 28 Apr.
- **A8:** the 95 genes against the autofluorescence map alone (Sami, 30 Apr).
- **S5:** Gria1's rank, the panel-test p and the roles p are held back until A1–A3 have run.

**Kinds of genes**
- Gene sets fixed before looking (SynGO).
- A background panel of a few thousand Allen coronal genes.
- Whether the category violin stays.

**Surface or total (the arms)**
- A total-receptor channel on some of the same brains.
- The green channel's filter sets: does tag survive, or does nano leak in?
- Rerunning the arms under A3. Today they have no minimum-mice rule.
- Retiring `arms_vs_genes.png`.

**Beyond**
- A2 moves its covariates, and A1 moves its `zref` slightly.
- The abundance term: Gria1 alone, or the four subunits as separate predictors.
- A calibration that uses the same model as nano.
- A7, for any claim about single structures.
- What the leftover is: a total-receptor channel, a specificity control, a DAPI control.

**Not ISH, listed so nothing is lost**
- **A4 and S2:** the adult distribution, meaning which structures sit reliably above the median. What "enriched" means is still to agree with Sami.
- **A5 and S3:** nano against autofluorescence, structure by structure.
- **A10:** videos; optional (S4).
- **Measurement validation:** a knockout or no-primary control, and the DAPI control run through the chain, ISH comparison included (Sami, 28 Apr).

## 8. A proposed structure for the ISH line

*(figure not kept: `figures/fig6_logic_map.png`)*

**Figure 6. The line today, on one page.** Each question, how it was answered, what holds, and what is open.

### The principle

- Ask each question once, with one method, one figure, and one test whose null respects the brain's smoothness.
- Groups of genes come from annotation fixed before looking.
- Set the surface question aside as one these data cannot answer, instead of leaving it half answered.

### Step 0. The inputs (not an analysis)

- **One structure set** (A1, S1), **one adult profile** (the mean `zref` of the 10 adults) and **one reading** (`zref`).
- **One gene table:** the union of the 95 and the 390 genes. Each experiment carries its section QC (A2), its reliability, its GO role, and its P9 category as a label only (A9).
- **One QC sheet per gene:** the native section profile and the orientation check, as P9 drew them.

### Analysis 1. Gene by gene: which expression maps look like the nano map?

- **Question:** which genes' mRNA maps order the structures as the nano map does, beyond what any smooth brain map would?
- **Method:**
  - Spearman per gene on step 0's set. A p per gene from spatial surrogate maps (A7), with BH across genes, and the same for ρ(Cacng8) − ρ(Gria1).
  - The autofluorescence map goes through the same steps (A8). A robustness table covers Spearman and Pearson, full and eroded means, and three structure sets (A3).
- **Figure:** every gene's ρ in rank order, with the null's 95% band and the autofluorescence ρ beside it. Below, scatters for Cacng8, Gria1 and one control, as in Figure 1.
- **What would mean what:**
  - Gria1 above the null and autofluorescence not: the map follows receptor expression and not the tissue. That is a specificity result.
  - Cacng8 − Gria1 above the null: the map follows a TARP's pattern more closely than its own subunit's mRNA. This is a description, not evidence of surface: Cacng8 is rich in hippocampus, and so is the map.
  - Nothing above the null: the ranking is descriptive only, and the text says so.

### Analysis 2. Between divisions or within them?

- **Question:** does a gene follow the map inside divisions, or only through the shared contrast of cortex and hippocampus against thalamus and hypothalamus?
- **Method:**
  - Per gene, ρ with a map that knows only each structure's division, and the mean ρ inside divisions with at least 8 structures.
  - The within-division ρ is tested by shuffling structures within divisions. That is a cheap first spatial null.
- **Figure:** Figure 5.
- **What would mean what:**
  - A high ρ within divisions means the gene and the map vary together at a fine scale. That is the stronger claim.
  - A whole-brain ρ that vanishes within divisions means they share only the gross gradient.
- **Today:**
  - The division-only map orders the genes as the real one does (0.98).
  - Inside divisions the order agrees at only 0.48, and the median ρ drops from +0.29 to +0.05.
  - The top genes within divisions are Cacng8 +0.45, Htr3a +0.43, Dlg2 +0.42, Grm5 +0.38 and Gria1 +0.37. Camk2a and Slc17a7 fall to near zero.
  - The astrocyte gene Aldh1l1 comes 7th (+0.27). Values this small need the null as much as the whole-brain ones.

*(figure not kept: `figures/fig5_between_within.png`)*

**Figure 5. What a whole-brain ρ is made of.**
- A: Cacng8 against the real map, and against a map that knows only the division.
- B: every gene, ρ on the real map against ρ on the division-only map.
- C: every gene, whole-brain ρ against its mean ρ inside divisions.

### Analysis 3. Which kinds of genes?

- **Question:** do kinds of gene that were defined before looking match the map better than others?
- **Method:**
  - Sets come from GO or SynGO and are fixed in a file before any ρ is read: subunits, localisation, other postsynaptic genes, presynaptic genes, glia and inhibitory neuron markers.
  - Each set's median ρ is tested against the spatial null. Localisation is tested against expression-matched postsynaptic controls, with the positive control (the `panel_test` design).
- **Figure:** the violin layout Sami knows: one column per set, one named dot per gene, and each set's null band.
- **What would mean what:**
  - Postsynaptic above presynaptic and glia, beyond the null: the map is postsynaptic-like. Any glutamate receptor label should give this, so it is a sanity check.
  - Localisation above matched controls: the localisation machinery shapes the map.
- **Today:** the first holds (+0.34 against +0.04, though not yet against a spatial null). The second does not (p 0.98, with the positive control at p 0.0004).

### Analysis 4. Is the map more than abundance and density?

- **Question:** how much of the reproducible map do receptor mRNA and synaptic density predict, and is the rest more than one Allen map disagreeing with another?
- **Method:**
  - A cross-validated fit on beyond's 126 structures. The predictors are the four subunits as separate terms, the markers, PSD and autofluorescence.
  - The leftover is set beside the same pipeline run with the same model on maps whose answer is known (the calibration).
- **Figure:** Figure 4, panels A and C: the variance budget and the calibration.
- **What would mean what:**
  - A leftover well above what a pure abundance-and-density map leaves (6–8%): part of the map is set by something else.
  - A leftover no larger than that of a map made of one Allen gene's mRNA: it cannot be told apart from ISH mismatch.
  - Either way, the leftover is not identified as surface.

### Analysis 5. Surface or total: a negative, documented

- **Question:** does the label report membrane receptor, or all receptor?
- **Method:** the green-channel check (`sep_channel_check.py`), on the raw channels with no ratios.
- **Figure:** Figure 3, panels B and C.
- **What it means:** the green channel tracks autofluorescence in every brain, so these data cannot answer the question. It needs a total-receptor channel.

### Keep, merge, drop

| piece | fate |
|---|---|
| P9 category violin | layout kept in Analysis 3. The hand groups, the split and the ANOVA are dropped |
| P9 per-gene scatter | kept in Analysis 1 (A6) |
| P9 fits within divisions | merged into Analysis 2 |
| P9 paired bars | dropped; the scatter shows the same |
| P9 Δz video | dropped (S4). A few per-gene videos stay optional (A10) |
| P9 section repair | replaced by A2: flag, never interpolate |
| P9 QC sheets (native profile, orientation) | kept in step 0 |
| P9 voxel Pearson | dropped |
| `regions.py`, `panel_build.py`, `panel_fetch.py`, `reliability.py` | kept in step 0 |
| `compare.py` | kept in Analysis 1, with one reading. `ratio` is one robustness row |
| `words.py` | the contrast named in advance is kept in Analysis 3; the top-12 bars are dropped |
| `roles.py` | dropped; `panel_test.py` supersedes it |
| `panel_test.py` | kept in Analysis 3 |
| `adult/arms.py` | kept, as the input to the green-channel check |
| `ish/arms.py` and `arms_vs_genes.png` | dropped |
| `sep_channel_check.py` | kept in Analysis 5 |
| `beyond_*` | kept in Analysis 4, with the changes of decision 5 |

**Order of work:**
1. A1.
2. A2 and A9.
3. Analysis 1 without the null (A3, A8).
4. A7.
5. Analyses 1–3 with the null.
6. Rerun Analysis 4.
7. Analysis 5's figure.

The ROADMAP places A6–A10 after the MATLAB scripts retire. This order brings A6, A7 and A8 forward for the ISH line (decision 4).

### Decisions

Five decisions for this discussion. S1 (the structure set) and S5 (no Gria1 rank yet) are already agreed; decisions 1 and 3 only confirm them.

1. **The headline for Sami.**
   - Recommendation: do not show p = 0.032 again.
   - Show Figure 1 (what one comparison is) and Analysis 1's ranked genes with the null band.
   - Quote "Cacng8 first" now. Quote Gria1's rank and the gap between Cacng8 and Gria1 only with A7's null (S5).
   - Show Figure 2 once, to explain why the old p does not stand.
2. **How to group genes.**
   - Recommendation: GO or SynGO sets, fixed in a file before the rerun.
   - P9's categories stay only as a label column in the gene table.
   - Keep the violin layout.
3. **The inputs: one structure set, one reading, one gene table.**
   - Structures: S1 as agreed, grey matter seen in all 10 adults (204 today, nearly P9's set). P9's nine divisions and "all structures" become rows of the sensitivity table.
   - Reading and statistic: `zref` and Spearman. Pearson on log2, eroded means and `ratio` go in the robustness table. Drop the voxel Pearson and the five-reading display.
   - Genes: the union of the two panels, with section QC, reliability, GO role, P9 category and experiment ids. Every analysis reads it. The 95 stay the genes named in figures.
4. **When to build A6–A8.**
   - Recommendation: as part of the ISH rebuild, before the MATLAB adult scripts retire. The ROADMAP puts them after.
   - No ISH p goes to Sami before A7. A6's shuffle within divisions is the cheap first null.
5. **What the channels and the leftover may claim.**
   - Arms: retire the tests of `ish/arms.py` and `arms_vs_genes.png`. Keep `sep_channel_check` as a methods figure. Remove "total receptor" and "surface fraction" from every label.
   - Beyond: enter the four subunits as separate predictors, and run the calibration with the same model. Quote the leftover as a range next to the calibration floor. Write "not predicted by receptor mRNA or synaptic density", never "beyond gene expression". Fix the ceiling and the stale docstrings.

**To raise with Sami** (wet lab and records, not analysis):
- a total-GluA1 stain or autoradiography on some of the same brains, the only answer to surface against total;
- a no-primary control;
- the green channel's filter sets;
- recording the adults' ages.

## Files

- Figures: `G:\sep_refactor\ish_discussion\figures\` (`fig1`–`fig6`; `p9_violin_spearman_ero.png`, a byte-identical copy of the figure on S:; `v2\`, copies of the production figures).
- Scripts for the figures: `G:\sep_refactor\ish_discussion\scripts\`.
- Read-only checks: `G:\sep_refactor\ish_discussion\scratch\` (`ish_numbers.py` and `.log`, `s1_count.py` and `.log`, `arms_by_cat.py`, `arms_beyond\`) and `G:\sep_refactor\ish_discussion\p9_original\`.
- P9 as of the tag `refactor-start`: `G:\sep_refactor\ish_discussion\P9_refactor_start.m`.
- Repository documents: `docs/REFACTOR_COVERAGE.md`, `docs/ROADMAP.md`, `docs/adult_ish_design.md`, `docs/SCIENTIFIC_CONTEXT.md`, `docs/history/REFACTOR_PLAN.md`, `mapping/README.md`.
