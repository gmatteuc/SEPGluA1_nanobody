# Step 8 report, 4 Oct

Written for Giulio on 4 October 2026, when the fixes of step 8 were applied
on the branch `refactor`, before the merge. Kept as the record of what each
fix changed and how it was checked. Paths on G: are working folders of the
refactor and may have been cleared since.

## In short

- Almost. Every fix you said yes to is on `refactor`, now at **`cfb7cd9`** in `G:\sep_refactor\check\code`. The tree is clean and nothing is pushed. All 38 commits since `47c252f` are yours (author and committer).
- Plasticity, preprocessing, both Python runs and the S6 rerun for fix 23 passed. Every output difference is traced to a fix, at the size the fix predicted. None is unexplained.
- **One check did not pass: registration (MG914).** It stopped at 15:55 on 3 Oct because C: was full (0 bytes free). It was never rerun. C: now has 162 GB free, so it can be rerun. See "Still open", item 1.
- Fix 23 keeps the RWS S1 result but changes the behavior comparison a lot. It also makes a new edge effect that needs a decision (item 2).

## What was applied

Numbers are the ones in [MORNING_REPORT.md](MORNING_REPORT.md). Each commit is "Step 8 fix: ...".

**Held fixes, your "yes" (1-9, 11-21; 10 left out as you asked)**

| # | commit(s) | what |
|---|---|---|
| 1 | 6bfba0c, 5d68ecd | subref's reference leaves out fiber tracts, ventricles, "brain-unassigned" and "unassigned". It is now exactly TH, HY, PAL, MB, P and MY |
| 2 | bc400fc, 71fce48 | the subref titles name that reference and what it leaves out |
| 3 | 726d2bc | the cohort mean no longer counts a missing ratio or sepratio voxel as zero |
| 4 | e5bf5fa | ish.words tests its features in sorted order, so the result no longer depends on PYTHONHASHSEED |
| 5 | f58c729 | each ISH panel permutation test has its own random generator |
| 6 | 7bc3bef | adult.arms' self-check bound allows for rounding in both tables, so run_adult_arms exits 0 |
| 7 | ab9b9b8 | adult.arms draws arms_consistency.png before it stops |
| 8 | 346bcf4 | nano_equalisation saves the inter-quartile range as a positive number |
| 9 | 4f41a0e, 2e6a002 | parafascicular is listed once; "mediodorsal" now measures MD, not IMD |
| 11 | eaa65f0, f80ea02 | one cohort table (`common/cohort.csv`); run_group_differences picks the behavior mice by name |
| 12 | 1ee77f5 | the young-vs-adult slice titles give the CCF plane (346 to 894) |
| 13 | a68d7f7 | panel F of the beyond regression names its CCF planes |
| 14 | e06ce73 | fig5 counts its genes (253) instead of saying 390 |
| 15 | aaa1ed2, 5a42eff | the region_plot ratio panel title reads "(log2; 0 = equally bright)" |
| 16 | e7dd275 | diagnostics sheet 01 says what the tissue mask needs |
| 17 | 113ddc8 | the ISH drops table gives the real reason (no energy.mhd in the zip) |
| 18 | 0a2b8f5 | add_sep_channel's panel title is drawn without TeX |
| 19 | ce57598 | "Hemisphere" is spelled right in the video titles |
| 20 | 3b5ffa3 | P6bis's background trace takes its y range from every mouse |
| 21 | e8452e9 | panel A's ceiling band in the EPS is a light band under the bars |

**New fixes (23-27)**

| # | commit(s) | what |
|---|---|---|
| 23 | ad6bf91, f9d6b66 | P7bis leaves voxels outside each mouse's tissue out of every mean, SEM and t (NaN-aware smoothing, per-voxel counts). Its t and surprise videos show only voxels that have a t |
| 24 | b42d398 | the cohort videos' reliability t is taken over each brain folded first; run_cohort writes `<reading>_folded_{mean,sd,n}.npy` |
| 25 | bc1e412, 4f0457f | the region tables and the young-adult P20-only maps leave a missing value out instead of counting it as zero |
| 26 | 68b36e1 | 'register' sets the registered grid itself (10 µm, `registered_grid_um.m`) |
| 27 | 9d9b25d, 560218e | the docstrings say "fixed and mounted", not "cleared" |
| 27 | f4a57dc | test_backward_compat takes the young cohort from the cohort table, not a pinned 14 |
| 27 | f9f87c6 | explore_czi_G finds the scenes in both the adult and the young files |
| 27 | 7c998f8 | check_demba_to_allen's helper is now `atlas/qc/demba_to_allen.py`, with nothing left in tmp\ |
| 27 | d2a85a3 | the ratio video uses the palette's difference map (blue below 1, red above) |
| 27 | cfb7cd9 | printed and drawn conclusions follow the numbers, in short wording |

The stage file `G:\sep_refactor\stages\stage_plasticity.m` now sets `behavior_mice` (fix 11). The old version is kept as `stage_plasticity_step5.m`.

## Checks

| check | code | result |
|---|---|---|
| Static | cfb7cd9 | ruff and format pass; 66 `.py` files compile; 32 sepmap modules import; `--help` exits 0 for all 26 run scripts plus build_demba_atlas and demba_to_allen; `sep_test_path` gives 0 failures (37 folders, 248 names, only the known RUN_ALL duplicate); get_cohort returns 31 mice and its verify passes |
| Plasticity (small set, 40 min) | f9d6b66 | 135 files. naive and rws are the same. behavior has 2 different files (fix 20). comparisons has 16 different files (fixes 23 and 9) |
| Preprocessing, MG914 P2/P2bis (7 min) | f9d6b66 | 114 files: 98 same, 9 same render (1 to 4 anti-aliasing pixels, which vary from run to run), 3 different (fix 8), 4 inputs |
| **Registration, MG914** | f9d6b66 | **Not done.** P4 'register' failed on slice 41 of 46 (`mhd_read` reshape error on a truncated transformix file) while C: had 0 bytes free. The fix-26 guard and get_cohort lines printed correctly before it stopped |
| Python route vs reference (2 h) | f9d6b66 | comparisons_v2: 219 files, 142 same, 77 different. adult_v2: 49 files, 36 same, 12 different, 1 new. All traced below |
| Python route vs the f9d6b66 run (2 h 55 min) | cfb7cd9 | comparisons_v2: 194 same, 17 different, 105 new. adult_v2: 41 same, 8 different. All traced below. Against the reference, the difference set is exactly the f9d6b66 set plus these, with nothing extra or missing |
| S6, production inputs, fix 23 (53 min, exit 0) | cfb7cd9 (no `.m` file differs from f9d6b66) | see "S1 before and after fix 23" |

## Output differences and their fixes

**MATLAB**

| output | what differs | fix |
|---|---|---|
| behavior `Background_mask_diagnostics_trace` (.fig/.png) | upper y limit 535.7 → 700.7; the line data are identical | 20 |
| `Slab_Avg_565_*_surpmask` (both comparisons) | t panels: the old t with the out-of-tissue voxels left out, to the bit | 23 |
| `Indiv_Slab_Avg_565_<group>` (4) | 0.24 to 0.48% of the pixels shown; values inside the tissue are identical | 23 |
| `Region_Surprise_Bar_DiffSum_*` | most bars change; more bars | 23, 9 |
| P2bis IQRs (3 `.mat` files) | exactly the negative of the old values (was −333 to −85, now 85 to 333) | 8 |
| everything else (every P5/P6bis `.mat`, the profiles, 36 Individual_Slice figures) | none | |

**Python, f9d6b66 against the reference**

| output | what differs | fix |
|---|---|---|
| `per_mouse/*_scalars.npz` (17) | `subcortex_mean` only (plus the known cache dates) | 1 |
| `ccf/*/subref_{mean,sd}`, `region_*`, `group_stats`, `cortex_table` | subref only; no other reading changes in any cell | 1 |
| `ccf/young(_P20)/sepratio_{mean,n,sd}` | the same 183 voxels: n 7→6 and 5→4, mean ×7/6 and ×5/4 | 3 |
| `region_table.csv`, `volumes_ccf20.npz` | subref; one ACB sepratio row by 0.0001; young sepratio at the 183 voxels, smoothed into `log2(_alt)_sepratio` (at most 0.13) | 1, 3 |
| `slices_cref/ratio/zref/sepratio` | plane titles only (sepratio also the known script-name title) | 12 |
| `slices_subref` | title, colour range, maps | 12, 2, 1 |
| `region_plot`, `group_plot`, `laminar_plot` | ratio panel title, subref title, subref panel | 15, 2, 1 |
| `video_subref_adult.mp4` | mean and t panels | 1 |
| `01_tissue_*` (17) | title line only | 16 |
| `arms_consistency.png` (new) | now drawn; run_adult_arms exits 0 and reports "agrees" | 6, 7 |
| `sep_channel_check.png` | title line only | 27 |
| `fig5_model_space`, `A_what_explains`, `F_maps` | identical to the figure drawn with that fix alone | 14, 21, 13 |
| `gene_region_table_panel_drops.csv` | reason for Gria1 and Negr1 | 17 |
| `ish_panel_test.png` | p: matched 0.7432→0.7481, positive control 0.0007→0.0006; `panel_test.csv` the same | 5 |
| `feature_enrichment.csv`, `ish_word_enrichment.png` | `gap_lo`/`gap_hi` only (median 0.006, at most 0.164); gap, p, q and n identical | 4 |
| `gene_correlations.csv`, `role_summary.csv` | subref rho only, at most 0.005 | 1 |
| `08_mask_vs_p6bis.png`, diagnostics README | the known differences | known |

**Python, cfb7cd9 against f9d6b66**

| output | what differs | fix |
|---|---|---|
| `ccf/<7 cohorts>/<reading>_folded_{mean,sd,n}.npy` (105 new) | new files; the existing files are identical | 24 |
| `ccf/adult/video_*_adult.mp4` (5) | t panel and its colour bar only | 24 |
| `region_means_per_mouse.csv` | 1 cell: MG897 ACB sepratio −1.6149 → −1.6143 | 25 |
| `region_stats.csv` | 3 cells of ACB sepratio, in the 4th decimal | 25 |
| `region_plot.eps`, `group_plot.eps` | one marker coordinate each; the PNGs are identical | 25 |
| `volumes_ccf20.npz` | the 5 `log2_alt_*` maps only | 25 |
| sheets 06, 08, 09; `fig1_ceiling`, `fig3_residual`, `C_where`, `D_controls`, `numbers_for_the_caption.txt`, `ish_roles.png` | titles and wording only (sheet 09's longer title moves its panels by a fraction of a pixel, proven by redrawing it with the old title) | 27 |
| `per_mouse/*_scalars.npz` (4) | `src_mtime` only | known |

## S1 before and after fix 23

Run on the S6 production inputs. Approved = the December 2025 figures. "Before" = the old code on the same inputs.

**Naive vs RWS, barrel field at slab 565**

| | approved | before | after |
|---|---|---|---|
| L−R map, pixels at p < 0.01 | 271 | 292 | 292 (the same pixels, same mean t 4.32) |
| L+R map, pixels at p < 0.01 | 532 | 437 | 287 (a subset of the 437) |
| L−R surprise bar | 14,155 (rank 1 of 36) | 14,562 (1 of 36) | 15,041 (1 of 42) |
| L+R surprise bar | 89,140 (rank 14 of 38) | 84,521 (13 of 38) | 81,773 (5 of 35) |
| per mouse L+R, naive vs RWS (mean ± SEM, n 5 vs 5) | 1.389 ± 0.143 vs 1.766 ± 0.161, p 0.118 | p 0.120 | 1.414 ± 0.135 vs 1.776 ± 0.157, p 0.119 |
| slab t against approved, r (L−R / L+R) | | 0.998 / 0.998 | 0.867 / 0.932 |

- The big L+R bars that were zeros counted as data collapse: midbrain motor 5.22M → 159, superior colliculus motor 1.50M → not drawn, VISp 0.96M → 4,208, HPF 0.91M → 14,493, RSP 588k → 18,180.
- No mouse's own S1 value moves by more than 0.006. The fix changes which mice enter each voxel's t.
- **New edge effect:** 9 edge pixels now have |t| > 10 (up to 73), where only about 2 mice per group remain. New L−R bars at tissue edges (olfactory tubercle, orbital, VISp, HPF, midbrain motor) are probably the same thing, but this was not checked directly.

**Naive vs behavior.** MG709 has no tissue at slab 565, but the old code counted it as a fourth behavior mouse (its zeros went through the alignment line). With the 3 real mice:
- L+R pixels at p < 0.01 go from 7,888 to 76,879 over the slab (about a quarter of it), and from 1 to 2,944 in S1;
- the L−R map goes from 14 significant pixels, all negative, to 286, mostly positive;
- per mouse, S1 L+R is +55% (1.372 ± 0.131 vs 2.119 ± 0.133, n 5 vs 3, p 0.009), the same before and after.

**Summary for Sami**
1. With the smoothing fixed so that missing tissue is left out instead of counted as zero, the RWS S1 (barrel field) increase at slab 565 is still there. The L−R map has the same 292 significant pixels; the L+R map has 287 instead of 437.
2. The barrel field stays the top region in the L−R surprise bars and rises from 13th to 5th in L+R, because the artefactual midbrain, hippocampal and retrosplenial bars disappear. Per mouse, barrel-field L+R is +26% in RWS (1.78 ± 0.16 vs 1.41 ± 0.14, n = 5 vs 5, Welch p = 0.12), the same as before the fix.
3. Behavior changes most: MG709 has no tissue at slab 565 but was counted as a fourth mouse. With the 3 real mice, behavior shows a broad L+R increase that includes S1 (+55%, p = 0.009), so that figure needs re-reading.

## Region list of the surprise bars (decision 9)

55 names. Each matches an atlas name exactly, and none appears twice.
- **Fixed:** the second "Parafascicular nucleus" (4f41a0e). Also "Mediodorsal nucleus of the thalamus", which had always measured the intermediodorsal nucleus (IMD, 88,205 voxels). The atlas spells MD differently, so the substring fallback found IMD first. It now measures MD (691,695 voxels; MATLAB gives the same counts). In the behavior sum panel, IMD's 47,967 becomes MD's 126,834.
- **Questions for you (nothing changed):**
  1. Nested regions count some voxels twice: SUB inside HPF (4.9% of HPF), STN and ZI inside HY (1.3%, 11.8%), SCm inside MBmot (26.5%). Keep both, drop one, or take the smaller out of the larger?
  2. Bars are sums, so big regions win over small nuclei (HPF 21.3M voxels, STN 0.1M). Sum, mean, or fraction of voxels over threshold?
  3. Was the second PF meant to be another region? PVT, SSp-n, SSp-m, SSp-un, VISpm, VISpl and AUDpo are all absent.
  4. Cosmetic: OT, PIR, SUB, CLA and BLA sit among the cortical names. The bars are sorted by value, so this has no effect.

## Subref numbers (decision 1)

Shift of each brain's subref, in log2.

| brain | fiber tracts and ventricles (6bfba0c) | catch-all labels (5d68ecd) | total |
|---|---|---|---|
| CGF027 | +0.355 | +0.039 | +0.394 |
| CGF028 | +0.344 | +0.039 | +0.383 |
| CGF033 | +0.269 | +0.019 | +0.287 |
| CGF034 | +0.267 | +0.039 | +0.306 |
| CGF035 | +0.242 | +0.026 | +0.268 |
| MG691 | +0.330 | +0.042 | +0.372 |
| MG692 | +0.278 | +0.033 | +0.311 |
| MG693 | +0.186 | +0.008 | +0.194 |
| MG736 | +0.356 | +0.042 | +0.398 |
| MG737 | +0.407 | +0.046 | +0.453 |
| MG897 | +0.286 | +0.030 | +0.316 |
| MG903 | +0.081 | +0.006 | +0.087 |
| MG904 | +0.168 | +0.014 | +0.182 |
| MG909 | +0.101 | +0.005 | +0.106 |
| MG910 | +0.031 | +0.005 | +0.036 |
| MG911 | +0.019 | −0.001 | +0.018 |
| MG913 | +0.175 | +0.018 | +0.193 |

- Adults move +0.19 to +0.45 (mean +0.34); young brains +0.02 to +0.32 (mean +0.13).
- Young−adult subref difference: −0.186, then −0.023, total **−0.209** (P20 only: −0.176 − 0.021; naive−rws: −0.016 − 0.0015).
- Structures with Welch q < 0.05: **104 → 86 → 81** (40 lost, 17 gained). Mann-Whitney: 97 → 78 → 75.
- Ratio (194), sepratio (0), cref (90) and zref (10) are unchanged.

## Fixes 24 and 25, measured

**24 (folded t in the cohort videos)**
- |t| rises in 89.3 to 95.7% of shown voxels, by a median factor of 1.10 to 1.14 for cref and zref. For ratio the factor is 1.02 to 1.06, for sepratio 1.00 to 1.04 and for subref 1.09 to 1.11.
- The t colour limit rises, for example adult cref 22.9 → 28.0, young_P20 cref 29.4 → 45.3, adult zref 14.8 → 16.9.
- zref's t changes sign in 0.5 to 1.1% of voxels, all near zero. cref never changes sign.
- The mean panel, the grey mask and the frames are unchanged, so the videos still match run_compare.
- Cost: 105 new files, 25 GB per run. run_cohort takes 55 min instead of 36.

**25 (NaN-aware region sums and young-adult maps)**
- Only MG897's sepratio is affected (117 voxels in ACB): one cell of `region_means_per_mouse`, three of `region_stats`, 0 of `group_stats`, `region_table` and `region_means_arms`.
- The P20-only `log2_alt_*` maps: 6,378 voxels become NaN in every reading. About 534,500 move, by a median of 2e-5 to 5e-5. Between 37,000 and 59,500 move by more than 0.01, and the largest move is 1.4 to 2.2. No figure draws these maps. The pooled maps are identical.

## Still open

1. **Rerun the registration check (MG914).** C: now has 162 GB free. The failed run left MG914's check-tree `lightsuite` folder half rewritten (no `volume_registered`). The earlier outputs are in `G:\sep_refactor\check\scratch\step8b_registration_before\`. The rerun overwrites the folder, then gets compared with `step8b_registration_compare.m`. Launch it with run_matlab_detached.ps1. MATLAB's temporary files go to C: (niftiread's unzipped atlas, transformix), so C: needs to keep room, or TMP/TEMP should point to G:. Your temp folder still holds 999 old `tp*` folders (about 20 GB) and thousands of `transformix_*` folders. Nothing was deleted.
2. **Fix 23, edge voxels (decision).** Voxels left with 2 to 3 mice per group can give very large t values and new edge bars. The options are a minimum number of mice per voxel (for example 3 per group), and/or summing the bars and the slab opacity only over voxels with their own t. At present the ±10-plane rolling median fills voxels that have no t from their neighbours. Either option moves the December comparison again.
3. **The behavior figure needs re-reading** (MG709, above).
4. **Region list:** the four questions above.
5. **Fix 24:** should the video's mean panel also use the per-brain average? That is a one-line change, and the mean would then differ from run_compare by more than 1% in 12 to 25% of voxels. Also, each frame's header "n = ... mice" still gives the larger side's count, not the t's folded n.
6. **Found under 25, not fixed** because each would move values where data exist: `per_mouse.block_means` averages unsectioned 10 µm voxels in as zero, and those means set the backgrounds; `to_ccf` warps the young brains linearly over the zeros outside the tissue, which darkens tissue edges; compare and video_compare apply the brain-count minimum with `cref_n` for every reading (no voxel is affected today).
7. **Fix 26:** the guard now stops on `px_register ≠ 20` and only prints a line for `px_atlas ≠ 10`, because only `px_register` sets the 10 µm grid. Reverting that one function restores the old check.
8. **Production reruns when the merge is done:** delete the `*_scalars.npz` caches first (they are keyed on file date only), then run_cohort, compare, video and video_compare. run_cohort adds 26 GB on D:, which has 114 GB free. P7bis has been rerun only on the S6 copy on G:.
9. **Later:** question 22 (P2's 5 slices with fewer than two reference pixels need a guard, to discuss); output file names that record their settings; whether "same render" counts pixels or colour channels. Also a side note: zref's median and spread count "brain, unassigned" as one structure among several hundred.
10. **Cleanup when you are happy:** `G:\sep_refactor\check\scratch\step8b_python_after_f9d6b66\` (55 GB), `G:\sep_refactor\s8apply_test\` (28 GB), `G:\sep_refactor\s6\before_step8b\`, and the step8-* worktrees and branches.
11. **Then** step 10 (documents) and step 11 (merge), as planned.

## Where things are

- Python comparisons: `G:\sep_refactor\check\scratch\step8b_python_compare\` (f9d6b66) and `...\step8c_python\` (cfb7cd9; open `set_check.log` first).
- Plasticity: `...\scratch\step8b_plasticity_compare.*` and `step8b_probe\`. Preprocessing: `G:\sep_refactor\pre\compare_step8b.*`.
- S6: outputs in `G:\sep_refactor\s6\data\comparisons\`, numbers in `G:\sep_refactor\s6\scratch\step8b\`, log `G:\sep_refactor\s6\logs\step8b\_stage_s6.log`.
- Merge and static checks: `...\scratch\step8_pymerge\`, `step8_merge\`. Reviews of the fixes: `step8_review\`, `step8_yes\`, `step8_newm\`, `step8_newpy\`.
