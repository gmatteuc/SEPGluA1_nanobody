# Morning report, 3 Oct

Written for Giulio on the morning of 3 October 2026, after the overnight run
of step 7's structural changes and the step 8 fixes. Kept as the record of
the held fixes and his answers to them, which follow each item after
`--->`, as he typed them. Paths on G: are working folders of the refactor and
may have been cleared since.

## In short

- `refactor` is at **`47c252f`** (in `G:\sep_refactor\check\code`). The tree is clean and nothing is pushed. There are 63 new commits since `742d7b8`: 14 kind C and 49 fixes.
- Every check passed. The only outputs that differ are the known, explained ones.
- 21 fixes change an output, so they wait for you on `step8-pending` (head `426cfdb`). Answer them by number (see "The decisions" below). If you agree with every recommendation, the answer is **"10 no, all others yes"**.
- Full detail is in [FIXES_STEP8.md](FIXES_STEP8.md), with one row per item.

## What was done overnight

### Kind C (the structural part of the style pass)

- **Python (11 commits, `eba99bc` to `a81182a`):**
  - one plotting module, `sepmap/plotting.py`, holds the palette, `save_figure` and the coronal drawing;
  - the run scripts set the figure backend, instead of each module doing it on import;
  - helpers that were defined twice now live in one place;
  - the analysis parameters are in `settings.toml` (63 constants; values and types checked);
  - paths are `pathlib.Path`;
  - about 25 long functions are split into named steps;
  - CSV rows are read by column name;
  - `build_demba_atlas` takes its age through argparse.
  - `a81182a` fixes a name clash that review found (arms.py's imported `partial`).
- **MATLAB (3 commits, `a78b77d`, `307ff92`, `3e834ae`):**
  - one palette function, `common/sep_palette.m`;
  - duplicate helpers merged (overlap matrices, knee searches, video-writer setup, save blocks);
  - functions over about 70 lines split into a short main of steps.
- **How it was checked:** every moved statement was compared, old against new, and functions were rerun on real and synthetic inputs, with identical results.
- **Checked by reading only:** the video writers, `to_ccf`'s warp, diagnostics sheet 08, the GUI plugin's k/K/u/U keys, and `copy_raw_data`, `extract_and_center` and `order_slices` 'apply'.
- **Two behaviour changes off the checked paths:**
  - importing `sepmap` outside a run script no longer forces the Agg backend;
  - `build_demba_atlas` usage errors now exit 2 with argparse's message.

### Fixes applied (49)

All 49 are merged into `refactor`. None changes an output file of the production pipelines as they run today.

**Preprocessing**

- `83c7bf2` order_slices 'apply': writes into the mouse's own folder, not into whatever folder `sliceinfo.mat` names (safety).
- `f94a239` order_slices 'apply': also reports the discarded slices.
- `e9cb839` nano_equalisation no longer stops when `save_results` is off.
- `759785b` residual_correction no longer stops when `doPlotBkg` is off.
- `ddfb802` residual_correction no longer loads the annotation it never used.
- `ed35622` P2's selection figures are closed after saving (8% faster).
- `c829352` annotate_artifacts makes no empty folder for a skipped mouse.
- `a180301` copy_raw_data also refuses a network (UNC) destination.
- `6ab882d` explore_czi_G reads each file's metadata, and every key.

**Registration**

- `424b83a` 'register' stops if a brain's `px_atlas` is not 10. This is only a guard; setting the grid is question 26.
- `81e2187` "automatic annotation is not installed" is now a warning.
- `ed8b5a5` add_sep_channel's header says the tissue was not permeabilised (a comment only).
- `3cffe3d` the message now reads "run_add_sep_channel: done." instead of "P4bis: done.".
- `6cbc859` verify_demba_setup checks the atlas, not P4's source text (12 passed).
- `364820c` verify_demba_setup's closing note is up to date.
- `a173b04` the engine's section modes answer the GUI on a SystemExit.
- `f824af4` the engine's README lists the `sections` mode.

**Plasticity and shared code**

- `35d2855` collect_by_group's dead code is gone. It reads a third as much, and `nano_4d.mat` is identical.
- `68617ab` the unused `correction_type` setting is removed from run_collect_by_group.
- `952afee` P6bis errors on a missing diagnostics folder instead of writing to the current folder.
- `9926552` P6bis prints "Mouse i / n" from the stack's count, not a fixed 5.
- `c6e9517` P6bis warns when it falls back from the robust fit.
- `87ea18a` P7bis's coarse-only run no longer reads `half_width` before it is set.
- `a4f25d9` compute_lr_stats pairs the planes correctly for an odd width. Every caller uses an even width, so nothing changes.
- `3dc5ef9` P7bis's surprise videos run without the t-score videos.
- `3b57c78` welch_surprise sets every small-group df to NaN, not only the first.
- `e441f86` P7bis names an unknown `exp_type`.
- `5867f7e` P7bis's "Recomputing background masks" message is corrected.
- `dd46723` P7bis's profile figure works with more than five mice in a group.
- `6755076` get_color2color_colormap's dead test is removed.
- `a8c2aa9` atlas_diagnostics no longer writes "Right now NaN of ...".

**Python route**

- `5e18ae4` cohort.py's docstring says which structure set its zref spreads are measured on.
- `0990970` run_ish_compare stops with a clear error, before writing anything, if P9's summary is missing.
- `71c7c91` ish.compare handles an old ranking with no machinery gene.
- `ffba644` `videos.min_n` has a `young_P22` key.
- `114ac75` diagnostics sheets 01 and 02 read `mad_k` from the settings.
- `3e39439` five readers close their files (`with` blocks).
- `ad9c59b` sheet 03 drops an unused `atlas_grid` call.
- `415e860` sheet 07 is skipped, with a printed line, when the young_vs_adult tables are missing.
- `6faebef` region_plot: the unused cohort unpack is removed and the figure is closed.
- `366dcf2` four unused constants are removed.
- `e4abc09` beyond_controls.centroids raises a clear error instead of returning None.
- `32b88d1` ish.words.figure is skipped when no feature was tested.
- `5425da0` region_plot's cortex printout names the ages correctly (P20 + P16 + P22).

**Tools**

- `7a950ec` the Python check tools report an unreadable file instead of crashing.
- `df84961` sep_make_script_copy refuses a setting that is assigned twice.
- `6e588ed` sep_compare_outputs works on files under a drive root.
- `47c252f` compare_outputs.py compares palette images by their colours.
- `781efb8` tools/README names `clear; clc; close all;`.

## Checks on `refactor` at 47c252f

All passed.

- **Static:**
  - `sep_test_path`: 0 failures;
  - checkcode: no new message;
  - ruff, format, isort and compile all pass;
  - all 34 sepmap modules import;
  - all 26 `run_*.py --help` and `build_demba_atlas --help` exit 0.
- **Plasticity (31 min):** 135 files, with **0 different**.
  - 125 are identical, 9 PNGs differ by 1 to 3 anti-aliasing pixels, and 1 is an input, which is not rewritten.
  - Every `.mat` (12), every `.fig` (61) and every P7bis output (20) is identical.
- **Preprocessing, P2 and P2bis on MG914 (7 min):** 114 files, with **0 different**.
  - 99 are identical, 11 PNGs differ by 1 to 3 anti-aliasing pixels, and 4 are inputs.
- **Registration of MG914 (48 min):**
  - every registered volume is identical page for page, and `transform_params.mat` is the same;
  - the GPU proposal moved by a median of 0.06 px (6 px at most) against step 7. That is inside the spread between two runs of the old code (13 px at most).
- **GUI drive scripts:** all 9 give the same files, apart from the click and save times and the usual GPU spread.
- **Python route (112 min):** 267 outputs were rewritten. 259 are identical, and the other 8 are the known ones: the cache dates, and the script names in titles and the README.
  - `run_adult_arms` still fails its self-check exactly as the reference does. Decision 6 fixes that.

## The decisions (fixes held on `step8-pending`)

Each fix is one commit. Each "yes" is cherry-picked onto `refactor` and its check rerun.

### Fixes that change numbers

1. **`95f2354`: subref's reference leaves out fiber tracts and ventricles.**
   - The bug: `NOT_SUBCORTEX` names "fiber tracts" and "VS", which never match, so white matter and ventricles sit in the reference today.
   - Each brain's subref moves by one constant:
     - adults by +0.19 to +0.41 log2 (mean +0.30);
     - young brains by +0.02 to +0.29 (mean +0.12).
   - The young-adult subref difference moves by -0.19 log2, and the number of structures with q < 0.05 drops from 104 to 86.
   - ratio, sepratio, cref and zref are identical.
   - **Recommend yes**, then redo the subref numbers: delete the `*_scalars.npz` caches, then rerun cohort, compare, video and video_compare.
   - Also to answer: should "brain-unassigned" and "unassigned" leave the reference as well? They are 1.5 to 2.3% of it. ---> do it, we should correct
2. **`ccc5985`: the subref titles say "excluding HPF and STR", which is wrong.** Only the titles change. **Yes, together with 1.** ---> labels must be always consitent with what the analyisis is actualyl doing so i would say yes
3. **`426cfdb`: the cohort mean counts a missing ratio or sepratio voxel as zero.**
   - Only MG897's sepratio has such voxels: 183 of them, 0.001%.
   - Fixing it takes the young sepratio n from 7 to 6 there and raises its mean by 14% (young_P20 by 20%).
   - Nothing else changes. **Recommend yes.** ---> yes fix it
4. **`0442402`: the ish.words bootstrap no longer depends on the order of Python sets.**
   - `gap_lo` and `gap_hi` move in about 2215 of 2400 rows (median 0.005, at most 0.16).
   - Every gap, p and q is identical.
   - From then on the output is the same whatever PYTHONHASHSEED. **Recommend yes.** ---> yes fix it
5. **`7259797`: each permutation test of the ISH panel test gets its own random generator.**
   - The p values move within Monte Carlo error: matched 0.7432 to 0.7481, positive control 0.0007 to 0.0006.
   - `panel_test.csv` is identical.
   - **Recommend yes.** The quoted 0.74 becomes 0.75. ---> ok do it
6. **`0fcd1c3`: the adult.arms self-check bound goes from 5.0e-05 to 5.05e-05.**
   - Floating point puts 12 of 5182 values just above 5e-05.
   - `run_adult_arms` then exits 0, and the tables are unchanged. **Recommend yes.** ---> fine do it
7. **`087496e`: adult.arms draws `arms_consistency.png` before it stops on a failed check.**
   - The dashed line sits at the real bound, and the title gives the verdict. **Recommend yes.** ---> ok do it makes sense
8. **`ea70ed7`: P2bis saves the inter-quartile range negative.**
   - The fix flips the sign exactly (MG914: -333..-85 becomes 85..333). Nothing reads these values. **Recommend yes.** ---> sounds wierd, fix it
9. **`3e12dd0`: 'Parafascicular nucleus' appears twice in the surprise bars' region list.**
   - The sum panel goes from 26 bars to 25, with the values unchanged. **Recommend yes.**
   - Also to answer: was the second entry meant to be another nucleus? ---> no ide why there was a second entry, do a doubel check if the panel of regions makes sense, maybe was a manual error aseembling the panel, anyway we do not want double counts so fix it
10. **`848eb57`: the automatic annotation runs on the CPU, so reruns give the same points.**
    - The GPU cannot be made deterministic.
    - CPU against GPU: median 0.11 px, and 2 of 46 sections swap a landmark (about 150 px).
    - It takes 42 to 51 s per brain instead of 25 s.
    - **Recommend no**, unless exactly repeatable proposals matter to you. ---> indeed we do not care so much about exact determinism, speed is imporatnt
11. **`5f8d1dd`: one cohort table, `common\cohort.csv`, for MATLAB and Python (Y4).**
    - Each route keeps its own order, so no output changes. **Recommend yes.**
    - The alternative is one shared order, which brings last-bit changes in the cohort means.
    - Choosing the behavior mice by name, instead of by index, would follow on from this. --- > do as recommended

### Fixes that change only labels or titles

12. **`0c6b75d`: the young-vs-adult slice titles give planes 180 too high.**
    - They now read CCF 346 to 894 (they read 526 to 1074). The maps are identical.
    - **Recommend yes:** the grant's slice figures carry the wrong numbers. ---> correct!!
13. **`ac39887`: panel F of the beyond regression.**
    - The rows are labelled CCF planes 610, 710 and 810 (they were labelled 4.3, 5.3 and 6.3 mm); the panels shift a few pixels. **Yes.** --> fix
14. **`3d55e1e`: `fig5_model_space` says "all 390 genes", but 253 are used.** Only the x label changes. **Yes.** ---> correct
15. **`619e136`: region_plot's ratio panel title says that its zero carries no meaning.** **Yes, or another wording.** ---> no give  a concise and intuitive descriptions of its meaning but since it's a title avoid being verbose
16. **`4945b48`: diagnostics sheet 01's title says what the tissue mask needs.** **Yes.**
    - This is the one cherry-pick conflict: keep `{TISSUE['mad_k']:g}` in the title. ---> ok
17. **`eb144fc`: the ISH drops table gives the real reason for Gria1 and Negr1:** no `energy.mhd` in the downloaded zip. **Yes.** ---> ok
18. **`08418f3`: add_sep_channel's panel title is drawn without TeX**, so "run_register_to_atlas" is no longer garbled. **Yes.** ---> ok
19. **`fba48bc`: "Hemishpere" becomes "Hemisphere"** in the individual videos' titles. **Yes.** ---> sure
20. **`8c9912a`: P6bis's background trace takes its y range from every mouse**, not from the last one only. On behavior, 57 of 1800 points were off the panel. **Yes.** ---> yes
21. **`ab7913f`: in `beyond_figures`' panel A EPS, the ceiling band was printed solid over the bars.** It is now a light grey band under them. **Yes.**
---> sure fix

### Questions with no commit yet (these can wait until after the merge)

22. **P2, 5 of 1301 slices with fewer than two reference pixels.** They draw their overlay with the previous slice's slope. How should such a slice be drawn? Only those 5 PNGs would change. ---> not clear to me maybe we can have a guard special case that resque it a bit we can discuss later
23. **P7bis smoothing sets the voxels outside the tissue to 0, not NaN.** The zeros count in the group means and in the SEM's n. Fixing it changes the approved December 2025 comparison. Fix it, or leave it? ---> this should be definitaly fixed
24. **run_video's reliability t uses n = max(nL, nR).** Keep it, switch to min (conservative), or use exact per-mouse folding (which recomputes every cohort)?
---> i would say use exactly
25. **The region tables fill a missing ratio or sepratio voxel with 0** (the same 183 MG897 voxels as in 3). Leave it as it is, or make the sums NaN-aware?
---> every time we can we should eb nan aware!
26. **P4 'register':** should the registered grid be set explicitly? That decides what a re-extracted brain gets. ---> yes
27. **Smaller ones, any time:**
    - Are the sections cleared? Three docstrings say "fixed, cleared and mounted". ---> no they are not cleraed
    - test_backward_compat's `expected_young = 14`: keep it pinned? ---> no
    - explore_czi_G's SizeS key on an adult file.  ---> ok
    - `check_demba_to_allen` needs scripts that sit in `tmp\`. ---> boh no ideally they shoudl sit together and nothing bnededed in temp
    - Should "same render" count pixels or colour channels? --> no idea
    - write_ratio_video draws red to blue, the reverse of the palette.---> we should jkeep one convention
    - Output file names that do not record their settings. ---> ideally they should provided they do not become overéy verbose and long
    - The wording of conclusions printed whatever the numbers say. ---> sure but again avoid verbosity and overcomplication

## What remains before the merge

- **Your answers above.** Each "yes" is cherry-picked and rechecked one at a time. Three held commit messages have small wording slips (`0fcd1c3`, `3e12dd0`, `ab7913f`; see FIXES_STEP8.md). They get corrected as they are cherry-picked.
- **The MG904 question: is MG904 P20 or P22?**
  - The registry and the cohort table say P22, so the Python route puts it in `young_P22` on `demba_p22`.
  - If it is P20, that row changes, the brain needs registering on `demba_p20`, and the `young_P22` video key (`ffba644`) is revisited.
- **Step 10, documents:**
  - README (cover image `slice_order_montage.png`), SCIENTIFIC_CONTEXT, ADDING_DATA, ROADMAP, FIGURES, the per-folder READMEs, and the tests.
  - Stale names to fix on the way: STYLE.md should name `sep_palette`; `production_settings.md` still lists `correction_type` and the old constants, and so does `REFACTOR_COVERAGE.md`.
- **Step 11, the merge:**
  - close every MATLAB session and the landmark_refine worker;
  - stash, then re-apply, the uncommitted run settings in `P4_register_to_atlas.m` on main;
  - move the engine's `.venv` to `registration\auto_annotation\.venv` and run its self-test;
  - delete the nested duplicate `matlab_elastix-master`.
- **Not yet on the bug list:** three small MATLAB findings from kind C, none of which changes an output:
  - load_allen_atlas and group_differences clear variables that do not exist;
  - an `exist` check in plot_individual_slabs is always true;
  - group_t_and_surprise computes t maps that nothing reads.
- **Not covered by any check:** the engine's `sections` mode (the U key) is driven by no script.
- **Cleanup when done:**
  - `G:\sep_refactor\s8apply_test\` (28 GB);
  - the worktrees `G:\sep_refactor\wt\s7cpy`, `s7cm`, `s8apply` and `s8hold` (keep `s8hold` until the cherry-picks are done).
