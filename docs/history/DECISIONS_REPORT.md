# The four decisions of 4 October: report

Written for Giulio on the evening of 4 October 2026, after he decided the four
points left open by step 8 ([ROADMAP.md](../ROADMAP.md), section 1), and completed at 23:42 with
the measurement of decision 4; all four were merged on 4 and 5 October. Kept
as the record of how each was checked. Paths on G: are working folders of the
refactor and may have been cleared since.

## In short

- **Decisions 1, 2 and 3 are done and checked.** P5 and P6bis are unchanged. P7bis changes only where the decisions predict. The video header gives the folded n.
- **The RWS S1 result holds under decisions 1-2.** Slab 565 keeps the same 292 L−R and 287 L+R significant pixels. All large edge t values are gone. By summed L−R surprise the barrel field is first; by share of significant voxels it is 5th of 46.
- **Decision 4 is done, but not yet measured where it counts.** It was checked on synthetic data and on real files, read only. The headline numbers after the change are still pending. The Python route on the check tree was at step 2 of 24 at 20:35, with about 4 h to go. The recommendation is to adopt both 4a and 4b once those numbers are in (see below).
- **Three branches are ready in main's repo.** They are `post-dec12`, `post-dec3` and `post-dec4`. Nothing is pushed. Every commit is yours as author and committer. The branches share no file, so they merge without conflicts.

## Decision 1: a t only with at least 3 mice per group

A voxel now gets a t only when each group has at least 3 mice with tissue there. The bars and the slab opacity use only voxels that have their own t. The ±10-plane rolling median reads only those voxels and is kept only on them, so it no longer fills in a voxel's neighbours. The new setting is `min_mice_per_group = 3` in `run_group_differences.m`. A value under 2 stops the run.

On voxels that keep their t, nothing changes: the t maps are identical, NaN included, and the opacity only drops to 0.

**Voxels that lose their t** (production stacks, the same as the S6 copy):

| | t at 2 mice | t at 3 mice | lose it | slab 565 | barrel field at 565 |
|---|---|---|---|---|---|
| naive vs RWS | 157,939,079 | 151,060,949 | 6,878,130 (4.35%) | 10,212 of 298,454 | 292 of 10,624 |
| naive vs behavior | 139,702,711 | 121,166,338 | 18,536,373 (13.27%) | 34,312 of 291,415 | 832 of 10,200 |

- **Behavior loses most.** Of its losses, 18.0M are voxels where the behavior group has fewer than 3 mice. MG709 has no tissue at slab 565, so every voxel there needs all three other behavior mice. That is where the black holes in the behavior slab figure come from.
- **The edge t values are gone.** Largest |t| on the slab, before (step 8) → now:
  - RWS L−R: 73.4 → 5.96
  - RWS L+R: 59.0 → 8.02
  - behavior L−R: 144 → 3.88
  - behavior L+R: 145 → 14.6

  The behavior L+R pixels still above 10 lie at least 100 µm inside the tissue, with at least 4 naive and 3 behavior mice each. They belong to the broad L+R increase, not to the edge.
- **Check tree.** It has only 2 mice per group, so at 3 it gives no t at all. `stage_plasticity.m` now sets the minimum to 2 (the old stage is kept as `stage_plasticity_step8.m`). At 2, every P5 and P6bis output equals main's. The minimum of 3 itself ran end to end on the S6 copy (5 vs 5 and 5 vs 4 mice).

## Decision 2: the surprise-bar regions

- **The list** has 71 regions, built from the atlas:
  - every structure-level area of the Isocortex (43);
  - a declared, commented set of subcortical regions (28). These are all of the old ones; none was added.
- **Nesting removed.** A region that contains another region of the list gives up that region's voxels. The run stops if two regions still overlap or if a region is left empty.
- **The bar** is the fraction of the region's voxels with a t that are at p < 0.01. It is NaN when the region has no voxel with a t. The summed surprise moves to a new table, `Region_Surprise_DiffSum_<comp>.csv`, which also holds the voxel counts.
- **Group labels** come from the atlas.

| group | n | regions |
|---|---|---|
| Isocortex | 43 | FRP, MOp, MOs, SSp-n, SSp-bfd, SSp-ll, SSp-m, SSp-ul, SSp-tr, SSp-un, SSs, GU, VISC, AUDd, AUDp, AUDpo, AUDv, VISal, VISam, VISl, VISp, VISpl, VISpm, VISli, VISpor, ACAd, ACAv, PL, ILA, ORBl, ORBm, ORBvl, AId, AIp, AIv, RSPagl, RSPd, RSPv, VISa, VISrl, TEa, PERI, ECT |
| olfactory | 1 | PIR |
| hippocampal formation | 2 | HPF (without SUB), SUB |
| cortical subplate | 2 | CLA, BLA |
| striatum | 3 | OT, ACB, CP |
| pallidum | 1 | GPe |
| thalamus | 13 | VPM, VPL, VM, PO, LP, LD, VAL, MD, PF, RE, CL, RT, GENd |
| hypothalamus | 3 | HY (without STN and ZI), STN, ZI |
| midbrain | 3 | MBmot (without SCm), SCm, SCs |

- **New areas:** FRP, SSp-n, SSp-m, SSp-un, AUDpo, VISpl, VISpm, VISli, VISpor.
- **Split areas.** At the atlas's structure level, ACA, ORB, AI and RSP come as their parts. None of their voxels is lost. Keeping one of them whole would mean declaring it by hand.
- **Labels fixed:** OT is under striatum, PIR under olfactory, SUB under hippocampal formation, CLA and BLA under cortical subplate.
- **Checked independently in Python:** all 71 voxel counts match; no voxel is in two bars; the 43 areas cover all 61.6M isocortex voxels. MD, SUB, SCm and SSp-bfd have their old counts. HPF, HY and MBmot have their old counts minus the regions taken out (HPF 21,320,488 − 1,048,697).
- The full list with voxel counts is in `G:\sep_refactor\wt_scratch_d_p7\region_list_real.csv`.

## S1 under decisions 1-2 (December 2025 inputs, S6 copy)

**Slab 565, pixels at p < 0.01.** On this slab, S1 is the barrel field only.

| | approved (Dec) | step 8 | decisions 1-2 |
|---|---|---|---|
| RWS L−R, barrel field | 271 | 292 | 292 (the same pixels, mean t 4.32) |
| RWS L−R, whole slab | 342 | 306 | 292: the barrel field is all that is left |
| RWS L+R, barrel field | 532 | 287 | 287 (the same pixels, all inside the approved 532) |
| behavior L−R, whole slab | 8 | 286 | 0 |
| behavior L+R, barrel field | 1 | 2,944 | 3,011 (all 2,944 kept, mean t 5.40) |
| behavior L+R, whole slab | 33,083 | 76,879 | 76,109 |

**Barrel-field bars.** Share = the fraction of the region's voxels with a t that are at p < 0.01.

| | approved: sum, rank | step 8: sum, rank | now: share, rank | now: sum, rank |
|---|---|---|---|---|
| RWS L−R | 14,155, 1 of 36 | 15,041, 1 of 42 | 0.0020, 5 of 46 | 14,637, 1 |
| RWS L+R | 89,140, 14 of 38 | 81,773, 5 of 35 | 0.0112, 11 of 40 | 81,332, 5 |
| behavior L−R | 8,423, 4 of 29 | 8,685, 9 of 40 | 0.0013, 7 of 38 | 8,690, 6 |
| behavior L+R | 5.40M, 3 of 52 | 5.86M, 4 of 52 | 0.670, 12 of 59 | 5.84M, 4 |

- **Above the barrel field by share:**
  - RWS L−R: CL, AUDv, VISrl, OT;
  - behavior L+R: the whisker thalamus (PO 0.92, VPL 0.85, VPM 0.84), then SSp-n and SCm.
- **Per mouse, unchanged** (the individual-mouse figures are identical):
  - RWS L+R: +26%, p 0.12, 5 vs 5 mice;
  - behavior L+R: +55%, p 0.009, 5 vs 3 mice.

**Summary for Sami**
1. With a t only where at least 3 mice per group have tissue, the RWS barrel-field increase at slab 565 is unchanged: the same 292 L−R and 287 L+R significant pixels. It is now the only significant L−R spot on the slab, and the large t values at the tissue edge are gone.
2. In the new region bars (the share of a region's voxels at p < 0.01), the barrel field has the largest summed L−R surprise and is 5th of 46 by share. Per mouse, barrel-field L+R is +26% in RWS (p = 0.12), unchanged.
3. Behavior has no significant L−R pixel left. The broad L+R increase stays, with 3,011 significant pixels in the barrel field and +55% per mouse (p = 0.009, 5 vs 3 mice). By share, the top regions are the whisker thalamus (PO, VPM, VPL).
4. The share has no chance level yet. Every RWS L−R share is under 1%, and the rolling median changes what chance would give. Shuffling the group labels would give that level.

## Decision 3: the video header

Each frame's header now gives the n behind the t: the largest folded n among the plane's voxels with a t. Before, it gave the larger side's count. The mean panel stays as run_compare's maps. On the check tree the new header matches the folded-n file exactly. It equals the old header in every adult frame and is higher in 5 of 375 young frames.

## Decision 4: NaN-aware backgrounds (4a) and young warp (4b)

**What changed**
- **4a.** Each channel's 20 µm block mean is now taken only over the imaged 10 µm voxels. Before, it averaged in unimaged voxels as zeros, so a block half outside the section read half its value. The off-tissue medians that set the backgrounds leave NaN out. A new check stops the run if SEP is missing inside the tissue; today it passes on all 17 brains.
- **4b.** Each young channel is warped to the CCF as value × mask, divided by the warped mask. Before, the zeros outside the tissue were interpolated into the edge voxels. The adults are placed, not warped, so this darkening was the young brains' alone.

**Measured so far**
- **4a, fully imaged blocks:** bit-identical to the old code (68M blocks in CGF027, 61M in MG911).
- **4a, backgrounds and tissue:** backgrounds move by less than 1 count on all 17 brains, and cortex means by at most 0.1%. The tissue mask gains about 50,000 blocks and loses none (measured on CGF027 and MG911).
- **4b, interior voxels:** bit-identical, 33.6M in MG903.
- **4b, edge voxels:** about 2.5% of the tissue. They rise by a median of 19% (MG897) to 25% (MG903). Before, 0.17% of MG897's tissue voxels were darkened by more than half. The mean over the whole tissue moves by 0.24%.

**Before and after.** The "after" column fills when the route reaches step 6 (`dec4_compare.txt`).

| measure | before (main) | after |
|---|---|---|
| RL+AL zref, young − adult: difference, MW p, q | +0.232, 0.0007, 0.014 | pending |
| RL+AL zref by band (supra / granular / infra), difference (q) | +0.278 (0.006) / +0.270 (0.006) / +0.186 (0.110) | pending |
| RL+AL cref, difference, p | +0.077, 0.36 | pending |
| laminar contrast, young − adult, cref / zref | VISp +0.48 / +0.25; RL+AL +0.51 / +0.31; 6 of 9 systems at q < 0.05 in cref, 5 in zref | pending |
| structures at q < 0.05 (Welch / MW, 238 tested) | ratio 194/195, cref 90/80, subref 81/75, zref 10/16, sepratio 0/0 | pending |
| zref spread, young / adult | 0.98 / 1.93 | pending |
| L1 cref in the CCF, young / adult (means) | 0.94 / 0.80 | pending |
| outer cortical shell cref, young / adult | 0.56 / 0.38 | pending |
| backgrounds (from the step-1 log) | e.g. CGF033 auto MAD 17.8, cortex 784 | under 1 count; auto MAD by at most 1 (17.8 → 17); cortex means by at most 0.1% |

- **Rows 1-6 cannot see 4b.** They are computed from each brain's own atlas. 4a's shifts are too small to move them beyond rounding. They should stay the same, but that is not yet measured.
- **Rows 7-8 should move.** They are read in the CCF. The young surface should get brighter there and the adults should not.

**Recommendation: adopt 4a and 4b together, once rows 1-6 are measured and unchanged at their rounding.**
- **Why.** Both fix a measured error rather than change a choice. Each leaves untouched every voxel the error did not touch.
- **What 4b changes.** It changes the CCF maps at tissue edges, mainly the cortical surface: the cohort maps, compare's `region_table.csv`, the videos, the closeups and the flatmaps. The young outer cortical shell already reads brighter than the adults' (0.56 vs 0.38), and after 4b it will read brighter still. Any map-level statement about the surface or L1, young against adult, should be read again on the new maps.
- **When not to adopt yet.** If a headline number moves beyond its rounding, or a q crosses 0.05, look at it first.

## Branches and what merging means

All three branches are in main's repo (`D:\sep_histology\code`), on main `b651cfe`. Your uncommitted `registration/run_register_to_atlas.m` is touched by none of them.

| branch | commits | files | merging into main means |
|---|---|---|---|
| `post-dec12` | ab4f8f0, e230218 (decision 1); c9831e9, 3253c11, 44b29fc, a2da518 (decision 2) | `group_comparison/` (P7bis, its rolling video, the driver), its README, `docs/production_settings.md`, `docs/FIGURES.md` | the next P7bis run uses the 3-mice rule and the new bars; nothing else changes |
| `post-dec3` | 847aaf8 | `mapping/sepmap/young_vs_adult/video.py` | only the cohort videos' header text changes, at the next run_video |
| `post-dec4` | 847aaf8, then b4e7829 (4a), ff7b5f9 (4b), 3e51d61 (4a backgrounds and the SEP check) | adds `volumes/per_mouse.py` and `volumes/to_ccf.py` | includes decision 3; every Python output changes at the next rerun from run_per_mouse on |

- **Merge order.** The first merge can fast-forward; the second needs a merge commit, without conflicts.
- **Decision 3 on its own.** Merge `post-dec3` now and `post-dec4` after the measurement.
- **Docs to update after merging, in one commit:**
  - the ROADMAP's "Decisions open from step 8", and its lines 140–145 and 463–469;
  - `docs/ADDING_DATA.md` lines 228–229 ("four volumes");
  - the REFACTOR_PLAN progress entry.

## Production reruns that follow

1. **After `post-dec12`: P7bis, both comparisons, on D:.** Run it with `run_matlab_detached.ps1`, with TMP and TEMP on a data drive.
   - P5 and P6bis need no rerun.
   - This also brings fix 23 to production, which so far has only run on the S6 copy.
   - It writes to `comparisons\naive_vs_<group>_nano`. The approved December figures are in `naive_vs_rws` and `naive_vs_behavior`, which it does not touch.
   - The S6 outputs on G: are what it will give, apart from the videos.
2. **After `post-dec4`: the whole Python route on production.**
   - First delete the `*_scalars.npz` caches, which are keyed on file date only.
   - Then run run_per_mouse and run_to_ccf on all 17 brains, then the 22 steps after them in stage_python's order, plus video_compare. On the check tree this takes about 6 h.
   - This includes the step-8 reruns still due (run_cohort's folded files, compare, video, video_compare). run_cohort adds about 26 GB on D:.
3. **If only `post-dec3` is merged:** run_video on the cohorts. It needs run_cohort's folded files on production.

## Still open

- **Decision 4's after numbers.** The route runs in the background and may stop before it ends. If it stops, resume from the last finished step with `python -u G:\sep_refactor\stages\stage_python_dec4.py <step>`. After step 6, run `dec4_measure.py after` and `dec4_compare.py`; after step 7, run `dec4_video_frame.py`. The commands and environment are in `G:\sep_refactor\check\scratch\`.
- **Check tree housekeeping:**
  - Six adult brains' registered tiffs were copied in for decision 4 (MG693, MG736, MG737, CGF033, CGF034, CGF035). Move them aside before the next `stage_plasticity`, or it will stack 5 mice per group.
  - `dec4_baseline\comparisons_v2\per_mouse\` holds post-dec4 files for 11 brains (copied after the route had rewritten them). The before measurement is not affected. For a file-level comparison, 8 of the 11 can be re-copied from D:. CGF027, MG903 and MG911 need main's run_per_mouse rerun (a few minutes each).
- **Science, not decided:**
  - a chance level for the share bars (label shuffling, Sami's April point);
  - a faint negative RWS L+R band along the ventrolateral cortical edge (474 pixels, 324 within 100 µm of the edge, already in the approved figure, with at least 3 mice per group everywhere);
  - the RWS OT and VISp L−R bars, which stay;
  - the broad behavior L+R increase as a possible normalisation effect (alignment slope 2.93);
  - PVT (the possible intended second PF).
- **Small:**
  - `tools/sep_compare_outputs.m` does not compare `AlphaData`, so it calls an opacity change "same";
  - the fine and coarse region analyses (off in production) still use the old 55 names and allow 2 mice per region;
  - VISpm's label is lower case, as in the atlas;
  - since 4a, the diagnostic sheets draw unimaged voxels white rather than the no-data grey;
  - MG911 and MG904 have patches where auto is 0 inside nano-bright tissue, which leaves holes in their tissue masks, before and after 4a.

## Where things are

- **Implementation:**
  - `G:\sep_refactor\wt\d_p7` (`post-dec12`), with the test and logs in `G:\sep_refactor\wt_scratch_d_p7`;
  - `G:\sep_refactor\wt\d_py` (`post-dec4`), with the tests in `G:\sep_refactor\check\scratch\dec34_py`.
- **Review:** `G:\sep_refactor\review_scratch`.
- **Plasticity check:** `G:\sep_refactor\check\scratch\dec12_*` (main's outputs in `dec12_plasticity_before\`); log in `G:\sep_refactor\check\logs\dec12_plasticity\`.
- **S6:**
  - outputs in `G:\sep_refactor\s6\data\comparisons\`;
  - step 8's in `G:\sep_refactor\s6\before_dec12\`;
  - numbers in `G:\sep_refactor\s6\scratch\dec12\`.
- **Decision 4:**
  - route log `G:\sep_refactor\check\logs\_python_route_dec4.log`;
  - before numbers in `G:\sep_refactor\check\scratch\dec4_measure\before\`;
  - baseline in `G:\sep_refactor\check\scratch\dec4_baseline\`.

## Decision 4, measured (4 Oct, 23:42; full route on all 17 brains, every step exit 0)

Tables: `G:/sep_refactor/check/scratch/dec4_compare.txt`.

| measure | before (main) | after (post-dec4) |
|---|---|---|
| RL+AL zref, young - adult: difference, MW p, q | +0.232, 0.0007, 0.014 | +0.235, 0.0007, 0.014 |
| RL+AL zref by band, supra / granular / infra (q) | +0.278 (0.006) / +0.270 (0.006) / +0.186 (0.110) | +0.281 (0.006) / +0.273 (0.006) / +0.185 (0.092) |
| RL+AL cref, difference, p | +0.077, 0.36 | +0.078, 0.36 |
| laminar contrast young - adult, cref: VISp / RL+AL | +0.477 / +0.506 | +0.481 / +0.506 |
| laminar systems at q < 0.05, cref / zref | unchanged (largest q change: AUDp cref 0.0069 to 0.0102, still < 0.05) | |
| zref spread, young / adult | 0.985 / 1.928 | 0.987 / 1.930 |
| young outer cortical shell relative to interior | 0.548 | 0.599 (+0.12 log2 per brain) |
| young L1 relative to interior | 0.921 | 0.937 (+0.02 log2) |
| adult L1 relative to interior | 0.785 | 0.784 |
| structures at q < 0.05, Welch: ratio / cref / subref / zref | 194 / 90 / 81 / 10 | 194 / 91 / 83 / 6 |
| structures at q < 0.05, Mann-Whitney: ratio / cref / subref / zref | 195 / 80 / 75 / 16 | 198 / 79 / 75 / 17 |

- The headline numbers do not move beyond their rounding.
- The 4 zref structures that leave Welch q < 0.05 (VISal, PPT, RPF, BST) sat exactly on the threshold (q 0.049) and move to q 0.057 to 0.071, with p almost unchanged (e.g. VISal 0.0019 to 0.0018): the Benjamini-Hochberg boundary, not a change in their difference. Mann-Whitney gains one zref structure.
- The largest per-structure changes are small ventral hypothalamic nuclei (SUM, PMv: up to 0.46 log2) measured in 4 to 6 young brains: structures at the tissue surface, which the old zero-darkening and zero block means affected most.
- Before/after correlation of every reading's young - adult differences over 238 structures: 0.993 to 0.998.

**Recommendation: adopt 4a and 4b.** The headline results are unchanged; what moves is what the fixes were meant to correct (the young brains' tissue surface, small surface structures). Do not quote the zref count of q < 0.05 structures (10 or 6): it sits on the BH boundary.
