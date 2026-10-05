# Step 8 fixes

Kept as the record of every item of the step 8 list and what became of it;
the answers to the held fixes are in [MORNING_REPORT.md](MORNING_REPORT.md),
what was applied in the end in [STEP8_REPORT.md](STEP8_REPORT.md). Paths on
G: are working folders of the refactor and may have been cleared since.

One row per item of the step 8 list: the fixes applied (APPLY, now in `refactor`),
the ones held for Giulio (HOLD, on `step8-pending`), the ones skipped (SKIP, with
the reason) and the ones already done. Written overnight on 2-3 Oct; nothing was
pushed. All commits are authored and committed by Giulio Matteucci.

## APPLY: merged into `refactor`

Branch `step8-apply` (worktree `G:\sep_refactor\wt\s8apply`): 49 commits on
`3e834ae`, each one fix that changes no output file of the production pipelines as
they run today; each message says what it fixes, why production is unaffected, and
how it was checked. The rebase onto `refactor` was a no-op (`refactor` was still at
`3e834ae`), so the commits keep their hashes, and `refactor` was fast-forwarded to
them: **`refactor` head `47c252f`** (in `G:\sep_refactor\check\code`).

Checked per commit: the identity tools on the changed files (a header or docstring
fix comes out `same code`), checkcode (no new message), ruff check and format, and
a direct test of the fixed branch on scratch copies, old code (`3e834ae`) against
new, never on production and never in a check tree's data (the check trees were
only read). Scratch: `G:\sep_refactor\s8apply_test\` (28 GB; the old code is
exported in `base_code\`; delete the folder when done with it). On the branch's
final tree: checkcode over all 77 `.m` files outside archive and third_party gives
no new message against `3e834ae` (one "might be unused" fewer in
residual_correction); the identity tools against `3e834ae` flag exactly the files
the fixes changed code in. Nothing under archive, third_party, docs, or P8/P9/P10
changed.

Checks on `refactor` after the merge (`47c252f`; logs in
`G:\sep_refactor\s8merge_check\`):

- `sep_test_path` in a fresh session (restoredefaultpath, `sep_setup_paths` of the
  check tree): 0 failures (37 folders, 247 names, the one known vendored
  duplicate `RUN_ALL`, nothing shared with `D:\dendrites\code`);
- `ruff check`: all checks passed; `ruff format --check`: 74 files already
  formatted (the 65 `.py` files and the Python blocks of 9 Markdown files);
- import order (`ruff check --select I`, `sepmap` first party): passes (with
  ruff's default classification the 20 `run_*.py` files report I001, as on
  `3e834ae`);
- all 65 tracked `.py` files outside archive and third_party compile;
- the 34 sepmap modules import in venv_atlas; config, closeup, plotting and
  hemispheres import in venv_flat;
- all 26 `mapping\run_*.py --help` and `atlas\build_demba_atlas.py --help` exit 0
  in venv_atlas.

## HOLD: on `step8-pending`, not merged

Branch `step8-pending` (worktree `G:\sep_refactor\wt\s8hold`): 21 commits on
`3e834ae` (head `426cfdb`), never to be merged as a whole: each commit changes an
output, says in its message what and by how much, and ends with a recommendation,
so Giulio can cherry-pick one at a time onto `refactor`.

Against today's `refactor` (`47c252f`), `git merge-tree` finds one textual
conflict: `4945b48` (sheet 01's title) and `114ac75` (sheets 01 and 02 take the
mask's MAD from the settings) change the same line of
`mapping\sepmap\diagnostics.py`. Keep both: the new wording with
`{TISSUE['mad_k']:g}` in place of the typed 4. The other files both branches touch
merge cleanly; each cherry-picked commit still needs its own check.

Measured on copies, never in place: everything under
`G:\sep_refactor\check\scratch\step8_hold\` (`data_py\` holds copies of the
check tree's Python inputs, `plast\` of its normalised plasticity volumes, `gpu\`
of MG914's annotation inputs; `results\` and `logs\` hold the before and after
outputs; `mprobe\` the MATLAB harnesses). The check tree's data was only read;
nothing under `G:\sep_refactor\check\data` was written. Checks on the branch's
final tree: ruff check, ruff format --check and isort pass; every file compiles;
the 32 sepmap modules import in venv_atlas, config and closeup in venv_flat; all 26
`run_*.py --help` exit 0; the engine's core and cli import in its own environment;
checkcode shows no new message in the 6 changed `.m` files; `sep_test_path` passes.

Three commit messages on `step8-pending` carry a wording slip, left as they are
(only the last commit could be amended without a rebase); the rows below have the
right wording:

- `0fcd1c3` gives the float difference as 5.0000000000033e-05; some of the 12
  values are 5.00000000001e-05.
- `3e12dd0` says the other P7bis outputs are identical to the check tree's; they
  were compared with the check tree's own run at the image level only (4 PNGs
  identical); the `.fig` files differ byte for byte by their stored dates.
- `ab7913f` says nothing else in the EPS differs; the diff left out the date and
  title lines.

## The table

| item | sorted as | branch and commit | what it changes | recommendation |
|---|---|---|---|---|
| `run_order_slices` 'apply' reads and writes the folder `sliceinfo.mat` names, not the cohort's (plan, safety, first) | APPLY | refactor `83c7bf2` | sets `sliceinfo.procpath` and `.volorder` to the mouse's folder before `generateReordedVolume`, says so when they differed; help and header updated. Shown on scratch: the old code wrote into the decoy "original" folder (the check tree's MG914 `sliceinfo.mat` names `D:\sep_histology\data`); the new one writes in place, byte-identical tiff to an in-place old run | take it |
| `run_order_slices` 'apply' does not report discarded slices (g2) | APPLY | refactor `f94a239` | the decisions line adds ", N discarded" | take it |
| `run_nano_equalisation` stops with `save_results` off (plan) | APPLY | refactor `e9cb839` | the run's time stamp is set before the saving; when not saving, the statistics figures are closed in its place (with them open the videos came out black after frame 1, found in the test). Save on: 11 outputs identical; save off: finishes, videos identical to a saved run | take it |
| `run_residual_correction` stops with `doPlotBkg` off (plan) | APPLY | refactor `759785b` | `select_reference_pixels` computes `used_clim` (saved in `slice_data`, read by ArtifactAnnotator) outside the figure, same arithmetic; the figure is saved only when drawn. On: all data identical; off: finishes, data identical | take it |
| the annotate mode's "automatic annotation is not installed" line is easy to miss (plan, Giulio) | APPLY | refactor `81e2187` | a `warning` (id `run_register_to_atlas:noAutoAnnotation`) instead of an indented fprintf line | take it |
| `run_add_sep_channel`'s header calls the tissue permeabilised (Giulio, 2 Oct) | APPLY | refactor `ed8b5a5` | header: stained without permeabilising (no detergent, or very little), so the nanobody labels the surface receptors; what the green channel was meant for and what it measured kept. Comment only (same code) | take it |
| `cohort.py`'s docstring quotes zref spreads 0.89 / 1.79 without the structure set | APPLY | refactor `5e18ae4` | names the set (each brain's structures of at least `region_tables.min_vox20` voxels) and gives the context draft's range, 0.9 to 1.0 against 1.8 to 1.9 log2 (production caches: 0.99 +- 0.30, n 7; 1.93 +- 0.28, n 10). Docstring only | take it |
| `explore_czi_G` reads `globalMeta` before setting it; stops at the first metadata key (plan, g2) | APPLY | refactor `6ab882d` | global metadata read per file; keys in a `containers.Map`; a missing SizeS reported as not found. Tested on a copy of an MG914 raw file: old stops, new prints all 20 series | take it |
| `collect_by_group`'s dead code (approved 2 Oct) | APPLY | refactor `35d2855` | auto and mask volumes no longer read, averages, `min_n_mice`, folder variables, the nine commented-out saves and the atlas `addpath` gone; help keeps the stale-`auto_4d.mat` warning. `nano_4d.mat` and `collected_mice.mat` identical; a third of the reads and memory | take it |
| `run_collect_by_group`'s setting `correction_type` is never passed (do-not-copy list) | APPLY | refactor `68617ab` | the setting goes; no stage assigns it (`docs/production_settings.md` still lists it, step 10) | take it |
| `run_residual_correction` loads the atlas annotation and never uses it (plan) | APPLY | refactor `ddfb802` | the `niftiread` goes; outputs identical | take it |
| P6bis writes diagnostics to the current folder when its folder is missing (plan) | APPLY | refactor `952afee` | `checks_folder` errors on an empty folder instead (the branch was dead since step 5, the folder is always passed) | take it |
| `run_ish_compare` skips the old-against-new check silently without P9's summary (plan) | APPLY | refactor `0990970` | `FileNotFoundError` naming the file, read first so nothing is written; outputs identical when it is there | take it; note that after step 9 the frozen summary must stay in the data |
| `ish.compare`'s old-ranking line indexes `mach[0]` unguarded (g8) | APPLY | refactor `71c7c91` | prints that the old ranking has no machinery gene instead of an IndexError | take it |
| `v2_video.MIN_N` has no `young_P22` key (plan) | APPLY | refactor `ffba644` | `videos.min_n.young_P22 = 1` (one brain, as `young_P16`); default cohorts do not include it | take it; revisit if MG904 turns out P20 |
| diagnostics sheets 01/02 type the mask's 4 MAD (g5) | APPLY | refactor `114ac75` | read `tissue.mad_k` (4.0); sheets identical | take it |
| files opened and never closed (g5, g6, g8) | APPLY | refactor `3e39439` | sheet 07, sheet 09, `atlas_grid`, `ish.regions.main`, `closeup.layer_thicknesses` read in `with` blocks; outputs identical | take it |
| diagnostics sheet 03's unused `atlas_grid` call (g5) | APPLY | refactor `ad9c59b` | removed; sheet identical | take it |
| diagnostics sheet 07 crashes without the young_vs_adult tables (g5) | APPLY | refactor `415e860` | skipped with a printed line; sheet identical when they exist | take it |
| `region_plot.structure_means` unpacks an unused cohort; `plot_regions` leaves its figure open (g6) | APPLY | refactor `6faebef` | both fixed; the four outputs identical | take it |
| P7bis coarse analysis reads `half_width` before it is set (plan) | APPLY | refactor `87ea18a` | the coarse block works the width out when the fine one has not run; ran coarse-only on scratch (finishes; ten other figures identical to the check tree's) | take it |
| unused `N_FOLDS`, `RNG` (beyond_controls), `LABEL` (adult.arms), `NICE` (sep_channel_check) (g7) | APPLY | refactor `366dcf2` | removed; syntax trees differ by those four assignments only | take it |
| `beyond_controls.centroids` returns None, then a TypeError (g7) | APPLY | refactor `e4abc09` | ValueError naming the count instead; same centroids otherwise | take it |
| `ish.words.figure` crashes with no tested feature (g8) | APPLY | refactor `32b88d1` | skipped with a printed line; outputs identical | take it |
| "P4bis: done." names the old script (g3) | APPLY | refactor `3cffe3d` | "run_add_sep_channel: done." | take it |
| `verify_demba_setup` checks P4's source text and fails today (plan) | APPLY | refactor `6cbc859` | checks `cohort_atlas_key('young', 20)` against the atlas instead; 12 passed, 0 failed on the check tree | take it |
| `verify_demba_setup`'s closing note is stale (g3) | APPLY | refactor `364820c` | says only `run_group_differences` still hardcodes the CCF and its crop | take it |
| the engine's section modes end without a response on a SystemExit (g3) | APPLY | refactor `a173b04` | `cli.py` also catches SystemExit there, so the GUI shows the message (tested in the engine's environment) | take it |
| the engine's README omits the `sections` mode (g3) | APPLY | refactor `f824af4` | Layout lists the three modes | take it |
| P4 `register` takes the registered grid from `sliceinfo.mat` with no check (plan) | APPLY (guard only) | refactor `424b83a` | 'register' stops before registering when a brain's `sliceinfo.px_atlas` is not 10 (all 31 production brains have 10, read only) | take it; the plan's other half, setting the grid explicitly, is Giulio's call (it decides what a re-extracted brain gets) |
| the Python check tools crash on an unreadable file (g9) | APPLY | refactor `7a950ec` | `check_code_identity.py`: SYNTAX ERROR for a non-UTF-8 file; `compare_outputs.py`: "compare failed" for a truncated .npz/.npy | take it |
| `sep_make_script_copy` replaces only the first of two assignments (g9) | APPLY | refactor `df84961` | refuses a setting assigned more than once (`cohort_specs`, `run_add_sep_channel`'s `mice_to_process`); every stage's settings still copy | take it |
| P6bis prints "Mouse %d / 5" (g4) | APPLY | refactor `9926552` | the stack's count; P6bis on scratch: outputs identical to the check tree's | take it |
| P6bis falls back from the robust fit silently (g4) | APPLY | refactor `c6e9517` | a warning naming the mouse; fits identical (tested with a failing `fitlm` stand-in) | take it |
| `compute_lr_stats` mispairs planes for an odd width (g1) | APPLY | refactor `a4f25d9` | each plane meets its mirror, the midline left out; identical for every caller (all 1140 wide) | take it |
| P7bis surprise videos need the t-score videos on (g4) | APPLY | refactor `3dc5ef9` | `t_lim` set before both branches; surprise-only run finishes, other figures identical | take it |
| `welch_surprise` sets only the first df to NaN for a small group (g4) | APPLY | refactor `3b57c78` | applies to every voxel; results identical (the division by zero already gave NaN p) | take it |
| P7bis `load_groups` falls through on an unknown `exp_type` (g4) | APPLY | refactor `e441f86` | names the bad value | take it |
| P7bis prints "Recomputing background masks" without recomputing (g4) | APPLY | refactor `5867f7e` | "Plane means of the tissue, on the saved background masks"; outputs identical | take it |
| `copy_raw_data`'s share guard compares drive letters only (g2) | APPLY | refactor `a180301` | a network (UNC) destination is refused too | take it |
| `tools/README.md` says drivers start with `clear all` (g9) | APPLY | refactor `781efb8` | names `clear; clc; close all;` | take it |
| `region_plot`'s cortex table says "(P20 + P16)" and "without the P16 brain" (g6) | APPLY | refactor `5425da0` | ages from the cohort lists ("P20 + P16 + P22"), "P20only = the P20 brains alone"; printout only | take it |
| `get_color2color_colormap`'s dead `size(gray, 1)` test (g1) | APPLY | refactor `6755076` | `m = 256`; identical maps | take it |
| `annotate_artifacts` makes an empty folder for a mouse it skips (g2) | APPLY | refactor `c829352` | no mkdir | take it |
| P2's selection figures never closed (g2) | APPLY | refactor `ed35622` | closed after saving; outputs identical, 8% faster slice loop | take it |
| P7bis profile figure crashes with more than five mice in a group (g4) | APPLY | refactor `dd46723` | the five shades repeat; figures identical for five or fewer | take it |
| `sep_compare_outputs.relative_files` drops a character under a drive root (g9) | APPLY | refactor `6e588ed` | trailing separators stripped; same verdicts elsewhere | take it |
| `atlas_diagnostics` writes "Right now NaN of ..." without the regression (g1) | APPLY | refactor `a8c2aa9` | a sentence saying the labels were not counted; note unchanged with the regression | take it |
| `compare_outputs.py` compares palette images by index only (g9) | APPLY | refactor `47c252f` | compared by their colours; no verdict on today's outputs changes | take it |
| P2: a slice with fewer than two reference pixels draws its overlay with the previous slice's slope, or crashes on slice 1 (g2) | HOLD | none (not in step8-pending either) | measured in production (read only): 5 of 1301 slices have no fit (MG904 #40, MG911 #50, MG912 #7, 13, 14); a fix changes those 5 `reference_pix_analysis` PNGs (overlay limits), no data | Giulio decides how such a slice should be drawn |
| `ish.words` bootstrap intervals follow the order of Python sets (plan) | HOLD | step8-pending `0442402` | `feature_enrichment.csv`: `gap_lo` and `gap_hi` move in about 2215 of 2400 rows (median 0.005, at most 0.16), 538 tied rows change place; every gap, p, q identical; the figure's whiskers; from then on identical whatever PYTHONHASHSEED (checked with 0, 1, unset) | take it |
| `adult.arms` stops before drawing `arms_consistency.png` (plan) | HOLD | step8-pending `087496e` | the figure is now written when the check fails; its dashed line sits at the check's bound (was 1e-9), its title gives the verdict; `region_means_arms.csv` identical; the run still exits 1 | take it |
| `adult.arms` self-check bound has no margin (plan) | HOLD | step8-pending `0fcd1c3` | bound 5.0e-05 to 5.05e-05 (half the last digit of both tables, plus 1e-12): 12 of 5182 shared values came out just above 5e-05 in floating point (exactly 0.000050 in decimal), none further; `run_adult_arms` exits 0; the table and `run_ish_arms` unchanged | take it |
| young-against-adult slice titles 180 planes too high (plan) | HOLD | step8-pending `0c6b75d` | `slices_<reading>.png/.eps`: only the six panel-title bands (1.1 to 1.3% of the pixels); planes now read CCF 346 to 894 (were 526 to 1074); maps identical; `run_replot` alike | take it (the grant's slice figures carry the wrong numbers) |
| TeX subscripts in `add_sep_channel`'s panel title (plan) | HOLD | step8-pending `08418f3` | `sep_channel\<mouse>.png`: that title reads `run_register_to_atlas`; the figure is 3 pixel rows shorter; nothing else | take it |
| subref's reference holds the white matter and the ventricles: `NOT_SUBCORTEX` names the categories "fiber tracts" and "VS", which never match a division (found tonight) | HOLD | step8-pending `95f2354` | each brain's subref moves by one constant, +0.19 to +0.41 log2 in adults (mean +0.30), +0.02 to +0.29 in young brains (mean +0.12); young-adult subref difference -0.19 log2; Welch q < 0.05 structures 104 to 86; ratio, sepratio, cref, zref identical. Not rerun: `run_cohort` (V2_READINGS=subref), `run_compare`, `run_video`, `run_video_compare`; delete the `*_scalars.npz` caches first | take it and redo the subref numbers; whether brain-unassigned and unassigned (1.5 to 2.3% of the reference) leave too is Giulio's call |
| P2bis saves the inter-quartile range negative (plan) | HOLD | step8-pending `ea70ed7` | the saved IQRs (`Intensity_Stats_*.mat`, `equalized_volume.mat` `stats_intensity_iqr_*`) change sign, exactly; MG914: -333 to -85 becomes 85 to 333; nothing reads them | take it |
| "Hemishpere" in the individual videos' titles (g4) | HOLD | step8-pending `fba48bc` | every frame of the individual videos: the title rows only | take it |
| P6bis background trace: y limit from the last mouse (g4) | HOLD | step8-pending `8c9912a` | `Background_mask_diagnostics_trace` figure; on the check tree's behavior figure 57 of 1800 points were above the panel, limit 535.7 to 700.7; naive and rws identical; no number | take it |
| subref titles say "excluding HPF and STR" (g5) | HOLD | step8-pending `ccc5985` | `slices_subref` and `region_plot`'s subref panel: the title only (0.18% and 0.09% of the pixels); tables identical | take it with `95f2354` |
| "all 390 genes" where 253 are used (g7) | HOLD | step8-pending `3d55e1e` | `fig5_model_space.png`: panel F's x label (0.4% of the pixels); two log lines give 253; every table identical | take it |
| beyond regression panel F labelled in mm from the cropped grid (g7) | HOLD | step8-pending `ac39887` | `for_sami\F_maps.png/.eps`: rows labelled CCF plane 610, 710, 810 / 10 um (were 4.3, 5.3, 6.3 mm); panels shift a few pixels (1.9%); other 23 outputs identical | take it |
| 'Parafascicular nucleus' twice in the region list (g4) | HOLD | step8-pending `3e12dd0` | `Region_Surprise_Bar_DiffSum_*`: the sum panel 26 to 25 bars (the duplicate gone, values unchanged), measured on naive against rws; the other 8 P7bis outputs identical | take it; Giulio says whether the second entry meant another nucleus |
| diagnostics sheet 01 says the mask needs "a section reached here" (g5) | HOLD | step8-pending `4945b48` | every `01_tissue_<mouse>.png`: the title line only (1.2% of the pixels) | take it; keep `{TISSUE['mad_k']:g}` from `114ac75` in the title when cherry-picking (the one conflict) |
| ISH drops table says "grid not downloaded" for grids `panel_fetch` refused (g8) | HOLD | step8-pending `eb144fc` | `gene_region_table_panel_drops.csv`: Gria1 (526) and Negr1 (695) read "grid not downloaded (no energy.mhd in the zip)"; region tables identical | take it |
| ISH panel test: one generator feeds its five permutation tests (g8) | HOLD | step8-pending `7259797` | every permutation p moves within Monte Carlo error (matched 0.7432 to 0.7481, all controls 0.6542 to 0.6513, positive control 0.0007 to 0.0006); each p no longer depends on which tests ran; `panel_test.csv` identical | take it (the quoted 0.74 becomes 0.75) |
| the ratio reading's zero carries no meaning (plan, "a caption fix at least") | HOLD | step8-pending `619e136` | `region_plot.png`: the ratio panel's title says so (0.08% of the pixels); the video and slices titles not changed | take it or reword |
| automatic annotation not reproducible on the GPU (plan) | HOLD | step8-pending `848eb57` | deterministic algorithms are refused on the GPU (grid_sample's backward has no deterministic CUDA kernel); the engine now runs on the CPU: reruns identical; CPU against GPU median 0.11 px, 2 of 46 sections swap a landmark (about 150 px); 42 to 51 s a brain against 25 s | leave it out unless exactly repeatable proposals matter |
| Y4: one cohort table for MATLAB and Python (plan) | HOLD | step8-pending `5f8d1dd` | `common\cohort.csv`, read by `get_cohort` and `sepmap.config.mapping_mice`; each route keeps its order (`mapping_order`), so `get_cohort`, `get_cohort_spec` and every Python list are identical and no output changes | take it, or drop `mapping_order` for one order and accept last-bit changes in the cohort means |
| `beyond_figures` EPS: transparent band printed solid (g7) | HOLD | step8-pending `ab7913f` | `A_what_explains.eps`: the ceiling interval is a 0.9 grey under the bars, not a solid 0.333 over them and the ceiling line; PNG 1.2% of pixels by one grey level | take it |
| cohort means count a missing ratio or sepratio voxel as zero (g5) | HOLD | step8-pending `426cfdb` | only MG897's sepratio has such voxels, 183 (0.001%); there the young sepratio n goes 7 to 6 and its mean +14% (young_P20 +20%); nothing else | take it |
| P7bis takes the behavior mice by index, not by name (plan) | SKIP | none | no output changes, but the setting `behavior_subset` is what the plasticity stage rewrites (`[1,2]`) and `production_settings.md` documents; by-name selection fits with Y4's cohort table (step8-pending `5f8d1dd`) | do it with Y4 |
| `get_allen_region_mask` takes the first name containing a root (g1) | SKIP | none | not silent (it prints the name taken), and production relies on it: P6bis's 'Retrosplenial' resolves to "Retrosplenial area"; a warning would fire on every run | none, or name the root exactly |
| `test_backward_compat`'s `expected_young = 14` (g1) | SKIP | none | a pinned snapshot of Sami's list; whether it should follow the registry is Giulio's call | Giulio |
| `explore_czi_G` asks for 'Global Information\|Image\|SizeS #1'; MG914's files name it 'Information\|Image\|SizeS' | SKIP | none | so it reports no scene information on those files; the adult files may use the other name | Giulio checks on an adult file |
| `register_to_atlas` align: unused `[H, W, Z]`, an empty `else`, `alignedvol`, loads without outputs (g3) | SKIP | none | the align mode is run by no check (an hour per brain), and the loads belong with kind C's load-into-structs | with that kind C work |
| `add_sep_channel`'s `fprintf(msg)` uses data as a format (g3) | SKIP | none | the text holds only digits and fixed words, so the output is the same; a style point | none |
| `compare_atlas_regions`' unused variables (g1) | SKIP | none | a hand-run QC script; harmless | none |
| `select_reference_pixels`' `idx_max_bis > 90` branch reads past the end (g2) | SKIP | none | unreachable: every caller keeps `p_max` at 75 or below | none |
| `compare_outputs.py`'s "same render" counts colour channels, its docstring pixels (g9) | SKIP | none | changing the code changes the check verdicts the refactor reports rest on; which is meant is Giulio's call | Giulio |
| `run_matlab_detached.ps1` Split-Stages: a transpose `'` flips the quote state (g9) | SKIP | none | no stage uses a transpose; a parser change | note in tools/README |
| `normalise_groups` clears variables not in its workspace and keeps `data_4d` (g4) | SKIP | none | memory only; a clear that misses a later use would stop P6bis, and checking it needs every cohort | with the next P6bis rework |
| "fixed, cleared and mounted tissue" in `adult.arms`, `sep_channel_check`, `cohort.py` docstrings | SKIP | none | only the permeabilisation was corrected (as asked); whether the sections are cleared is for Giulio to confirm | Giulio |
| `addpath(atlas_dir)` after `get_atlas` in three preprocessing functions; `normalise_groups` adds the CCF by hand (consistency) | SKIP | none | path order, not output; the plan's rule "only get_atlas adds an atlas folder" wants a pass of its own | kind C or Giulio |
| stale names in `docs/` (`correction_type` in production_settings, the old constants in REFACTOR_COVERAGE, `sep_palette` missing from STYLE.md, STYLE_PASS's panel_test line numbers and its "do not copy" list for `collect_by_group`) | SKIP | none | documents, step 10 | step 10 |
| P8's reliability bars start at white (plan) | SKIP | none | P8 retires; A4 must not copy the scale | none |
| `ArtifactAnnotator`'s close handler on repeated Esc (plan) | SKIP | none | not used in production; Giulio would replace or retire it (ROADMAP) | none |
| MG904's age | SKIP | none | waiting for Giulio: is MG904 P20 or P22? `common\cohort.csv` (Y4) and the registry say P22, so the Python route puts it in `young_P22` on `demba_p22`; if it is P20, that row changes and the brain needs registering on `demba_p20` | Giulio answers |
| P7bis smoothing sets every non-tissue voxel to 0, not NaN (g4) | SKIP | none | needs Giulio: it changes the approved December 2025 comparison. The zeros count in the group means and in the per-voxel SEM's n; the experimental group's zeros become intercept/scale after the alignment line; a side without tissue gives L - 0. A fix needs NaN off tissue and the per-voxel n carried into the SEM and the Welch df (which uses the largest n). Measurable on the small set in about 20 min per comparison (`plast\run_p7.m`, the outputs of the code before the fix are in `plast\results\p7_rws_base`) | Giulio decides |
| `run_video`'s reliability t uses n = max(nL, nR) (plan) | SKIP | none | needs a method choice: the two hemispheres' counts differ in 22% (adult), 22% (young_P20) and 30% (young) of the voxels with n >= 2, by a median factor 1.1 to 1.3. min(nL, nR) is conservative; per-mouse folding in `run_cohort` is exact but recomputes every cohort | Giulio chooses |
| `write_lr_video` boundaries from the ML gradient only (g4) | SKIP | none | since kind C the three left-right writers share `lr_atlas_boundaries`; only the individual writers use the full gradient, and their overlay is off; a design choice | Giulio decides |
| conclusions printed or titled whatever the numbers say (g7: "so this is real", panel E, fig4 to fig6 titles) | SKIP | none | the wording of scientific claims | Giulio rewords |
| output file names miss settings (`video_compare --dlim`, `closeup --vmax --dlim --smooth`) (g6) | SKIP | none | file names change every output; the owner decides (STYLE_PASS) | Giulio decides |
| `select_background_pixels`' first knee search window differs from the fallbacks' (g1) | SKIP | none | may be intended | Giulio says |
| `check_demba_to_allen` needs a file made by `tmp\demba_to_allen.py`, not in the repository (also `tmp\remap_demba.py`, `tmp\check_crop.py`) (g1) | SKIP | none | Giulio knows where those scripts are | add them or drop the check |
| region tables fill a missing ratio or sepratio voxel with 0 (`per_unit` fill=0; region_plot, region_groups, arms) (found tonight) | SKIP | none | by design per `per_unit`'s docstring; today the same 183 voxels of MG897's sepratio; a NaN fill needs NaN-aware sums in three modules and the arms check | Giulio decides, with `426cfdb` |
| `adult_profile` written three times, no minimum-mice rule (plan) | SKIP | none | replaced by A3 (step 9) | step 9 |
| LightSuite scripts that write to S:, the machine-specific demo (plan) | SKIP | none | marked do-not-run in `third_party\README.md` (step 4); nothing to fix | none |
| nested duplicate `matlab_elastix-master` on disk (plan) | SKIP | none | Giulio deletes it before the merge | Giulio |
| `video_compare` types the 0.02 floor (g6) | done | kind C | the value comes from `readings.log2_floor` since kind C; only its docstring still says 0.02, the same value | none |
| `collect_by_group`'s commented-out saves name renamed variables (check stage) | done | refactor `35d2855` | the commented-out saves are gone | none |
