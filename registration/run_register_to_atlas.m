clear all
close all
clc

% /// Registration, step 1 of 2: register each mouse's sections to its atlas ///
%
% Bridges what preprocessing made of each mouse (the corrected volumes of
% run_residual_correction, the equalised nano of run_nano_equalisation, the
% artifact masks of run_annotate_artifacts) into LightSuite, and registers the
% sections to the atlas of the mouse's age (atlas_key). It runs in five modes,
% one at a time, because the manual steps (the cutting angle, the control
% points) sit between the automatic ones:
%
%   run_mode = 'align'     (auto)   bridge the corrected volumes into LightSuite,
%                                   align the slices and fit the atlas rigidly.
%                                   Writes regopts.mat and volume_for_inspection.tiff.
%   run_mode = 'angle'     (MANUAL) optional, between align and annotate: set the
%                                   cutting angle by eye (determineCuttingAngleGUI).
%                                   Every adult had this done; without it the
%                                   rotation is the automatic rigid fit's, which
%                                   came out at 17-21 deg for MG897 and MG913.
%                                   One angle serves the whole brain: the planes
%                                   saved on individual slices are averaged into a
%                                   single normal by applyAngleToTransform, so save
%                                   it on three to five slices spread front to
%                                   back and let the average cancel the misses.
%                                   Controls: plain left/right change the slice;
%                                   SHIFT + arrows tilt the atlas (0.3 deg per
%                                   press, hold the key); the wheel moves the
%                                   plane along its normal; return saves the plane
%                                   for the current slice, c clears it; space
%                                   toggles the region outlines; 1/2/3 show one
%                                   channel, 0 all. Closing the window writes
%                                   cutting_angle_data.mat. Must come BEFORE
%                                   annotate: it changes the atlas block the
%                                   control points are counted in, so this
%                                   script refuses it once a mouse has points.
%   run_mode = 'annotate'  (MANUAL) open the control-point GUI on one mouse.
%                                   Besides placing points by hand, r proposes
%                                   points for the slice at the atlas plane on
%                                   screen (the neighbouring slice's landmarks,
%                                   refined through the image matcher; a ? marks
%                                   the ones it was unsure of) and t drops in a
%                                   plain copy. Both land provisional and are
%                                   never saved unless touched. Needs the Python
%                                   side once per machine: setup_landmark_refine.ps1.
%                                   Writes atlas2histology_tform.mat.
%
%                                   AUTOMATIC ALTERNATIVE, in three steps:
%                                   1 'annotate': only set the atlas plane on the
%                                     suggested anchor slices (j jumps between
%                                     them, a fixes the plane on screen), save (s).
%                                   2 'autoannotate' (below) proposes every slice.
%                                   3 'annotate' again: the proposal loads orange,
%                                     least confident points marked ?; k accepts a
%                                     slice, u re-proposes it at the plane on
%                                     screen, the usual tools fix points. Only
%                                     accepted or touched slices are saved.
%   run_mode = 'autoannotate' (auto) the automatic control points, from the anchor
%                                   planes: auto_proposal_controlpoints.mat, never
%                                   atlas2histology_tform.mat itself. Needs the
%                                   Python side once per machine:
%                                   registration\auto_annotation\setup.ps1 (a GPU
%                                   makes it minutes).
%   run_mode = 'register'  (auto)   elastix refinement and the registered volumes.
%                                   Picks up the control points if they exist.
%
% The adults were done this way, one mouse at a time, with the GUI lines
% uncommented by hand. Every one of them has control points on every slice, so
% a young brain registered without them is not being treated the same way --
% see the note on 'annotate' below.
%
% Each mouse then goes through run_add_sep_channel, which carries its SEP
% channel through the same registration.
%
% Run sep_setup_paths first, once per MATLAB session. The settings are below,
% the code is in pipeline\register_to_atlas.m.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Cohort selection (mice come from the shared registry get_cohort.m).
% Set mice_to_process to {} to process every mouse in groups_to_process.
groups_to_process = {'young'};                  % 'rws' | 'naive' | 'behavior' | 'young'
mice_to_process   = {'MG904_SepGluA_P22'};   % 'annotate' takes one mouse at a time

% Which half of the script to run. 'annotate' takes one mouse at a time.
% After 'align' (and optionally 'angle'), the control points can be made in
% two ways; both end in the same atlas2histology_tform.mat and 'register':
%   manual     'annotate' (click the points on every slice) -> 'register'
%   automatic  'annotate'     set the plane on 4 suggested slices (j, wheel, a), s
%              'autoannotate' proposes every slice (~30 s on a GPU)
%              'annotate'     review: k / K accept, u / U re-propose, fix points, s
%              'register'
% Details in the header above and in registration/auto_annotation/README.md.
run_mode = 'register';                             % 'align' | 'angle' | 'annotate' | 'autoannotate' | 'register'

% How far the atlas shown in the GUI (and used by the registration) extends
% beyond the slice stack, in slices, on each side. The atlas on screen is
% the rigidly pre-aligned atlas resampled onto the stack's AP range plus this
% margin, so if the automatic rigid fit lands the stack too far back, the
% true plane of the first slices sits outside the margin and the wheel
% cannot reach it (MG912: slice 1 needed ~7 slices beyond LightSuite's 6).
% Widening it costs nothing but memory. It must not change once a mouse has
% control points, because the saved atlas planes are counted from the start
% of this range -- this script leaves such a mouse at the margin it was
% annotated with. The three P20 mice registered before this existed keep
% LightSuite's 6.
atlas_extent_slices = 15;

% Reference atlas.
%
% DECIDED 2026-09-02: the young cohort
% registers to the age-matched DeMBA P20 template. The reasoning is that the
% manual control points carry the correspondence, so the adult template's
% better contrast -- 1.4x the global CV, 1.6x the local smoothed gradient --
% matters less than having a target with P20 proportions. Template contrast
% feeds the image-similarity terms; landmarks do not care about it.
%
% The adults stay on 'ccf' and are NOT re-registered.
%
% Consequence to remember: the two cohorts then live on different grids.
% Registered volumes come out at twice the registration grid, so adults land
% on [900 800 1140] and the young on [994 800 1140]. run_collect_by_group
% onward still assume the adult atlas and crop everywhere, so they must be made
% atlas-aware per cohort before any young data reaches them. Region-level
% comparison across the two is fine once that is done -- both annotations are
% in the same parcellation_index space -- but voxelwise cross-group work would
% need CCF Translator.
atlas_key = 'demba_p22';                        % 'ccf' | 'demba_p20' | 'demba_p16' | any age built

% Choose correction type
correction_type = 'slicewise';

% Set if to use equalized nano volumes
use_equalized_nano = 1;

% Register without manual control points. Off, and it should stay off for
% anything that ends up in a figure: LightSuite relies on the control points to
% register well, so an image-only run is a diagnostic, not a result.
allow_image_only_registration = false;

%% Run

% The settings above go to the code under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.run_mode = run_mode;
run_settings.atlas_extent_slices = atlas_extent_slices;
run_settings.atlas_key = atlas_key;
run_settings.correction_type = correction_type;
run_settings.use_equalized_nano = use_equalized_nano;
run_settings.allow_image_only_registration = allow_image_only_registration;
register_to_atlas(run_settings);
