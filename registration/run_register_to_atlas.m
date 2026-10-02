%% run_register_to_atlas
% ===== Register each mouse's sections to the atlas of its age =====
%
% Registration, step 1 of 2:
%   1. run_register_to_atlas   register the sections to the atlas  <- this script
%   2. run_add_sep_channel     carry the SEP channel into registered space
%
% Bridges what preprocessing made of each mouse (the corrected volumes of
% run_residual_correction, the equalised nano of run_nano_equalisation, the
% artifact masks of run_annotate_artifacts) into LightSuite, and registers the
% sections to the atlas of the mouse's age (atlas_key). It runs in five modes,
% one at a time, because the manual steps (the cutting angle, the control
% points) sit between the automatic ones:
%
%   'align'         (auto) bridge the corrected volumes into LightSuite, align the
%                   slices and fit the atlas rigidly; writes regopts.mat and
%                   volume_for_inspection.tiff
%   'angle'         (manual, optional, between align and annotate) set the cutting
%                   angle by eye (determineCuttingAngleGUI); closing the window
%                   writes cutting_angle_data.mat
%   'annotate'      (manual) open the control-point GUI on one mouse; writes
%                   atlas2histology_tform.mat
%   'autoannotate'  (auto) the automatic control points, from the anchor planes:
%                   writes auto_proposal_controlpoints.mat, never
%                   atlas2histology_tform.mat itself
%   'register'      (auto) elastix refinement and the registered volumes, with the
%                   control points when they exist
%
% The cutting angle: every adult had it set by eye; without it the rotation is
% the rigid fit's, which came out at 17-21 deg for MG897 and MG913. One angle
% serves the whole brain: applyAngleToTransform averages the planes saved on
% single slices into one normal, so save it on three to five slices spread front
% to back and let the average cancel the misses. Keys: left and right change the
% slice; shift and the arrows tilt the atlas (0.3 deg a press, hold the key); the
% wheel moves the plane along its normal; return saves the plane for the current
% slice, c clears it; space shows the region outlines; 1, 2 or 3 show one
% channel, 0 all. It must come before 'annotate', since it changes the atlas
% block the control points are counted in: this script refuses it once a mouse
% has points.
%
% Annotating by hand: besides placing points, t takes the neighbouring slice's
% points as they are, at the atlas plane on screen, and p carries them forward
% as you step. Both land provisional and are saved only when touched.
%
% Annotating automatically, where the engine is installed (once per machine,
% registration\auto_annotation\setup.ps1; a GPU makes it minutes):
%   1. 'annotate'      set the atlas plane on the suggested anchor slices only
%                      (j jumps between them, a fixes the plane on screen), save (s)
%   2. 'autoannotate'  propose every slice
%   3. 'annotate'      review: the proposal loads orange, the least confident
%                      points marked ?; k accepts a slice, u re-proposes it at the
%                      plane on screen, the usual tools fix points; only accepted
%                      or touched slices are saved
% Without the engine the GUI has no such keys.
%
% Every adult went through these modes one mouse at a time and has control
% points on every slice, so a young brain registered without them is not treated
% the same way (allow_image_only_registration below). Each mouse then goes
% through run_add_sep_channel, which carries its SEP channel through the same
% registration.
%
% Setup: one young mouse at a time, on the DeMBA atlas of its age; the adults
% stay on the CCF. Run sep_setup_paths first, once per MATLAB session; the code
% is in pipeline\register_to_atlas.m.

clear all
close all
clc

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% groups to process ('rws', 'naive', 'behavior', 'young'); the mice of each come
% from the cohort registry, get_cohort
groups_to_process = {'young'};

% mice to process, by name ({} = every mouse of groups_to_process); 'angle' and
% 'annotate' take one mouse at a time
mice_to_process = {'MG904_SepGluA_P22'};

% mode ('align', 'angle', 'annotate', 'autoannotate', 'register'), one at a time.
% After 'align' (and 'angle') the control points are made one of two ways, both
% ending in the same atlas2histology_tform.mat and 'register':
%   by hand     'annotate' (click the points on every slice), then 'register'
%   automatic   'annotate'      set the plane on 4 suggested slices (j, wheel, a), s
%               'autoannotate'  propose every slice (about 30 s on a GPU)
%               'annotate'      review: k or K accept, u or U re-propose, fix points, s
%               'register'
% (more in the header and in registration\auto_annotation\README.md)
run_mode = 'register';

% margin of the atlas beyond the slice stack, in slices on each side, for the GUI
% and the registration. The atlas on screen is the rigidly fitted atlas over the
% stack's AP range plus this margin, so when the rigid fit lands the stack too far
% back, the true plane of the first slices is out of the wheel's reach (MG912's
% slice 1 needed about 7 slices beyond LightSuite's 6); a wider margin costs only
% memory. A mouse with control points keeps the margin it was annotated with,
% since the saved atlas planes are counted from the start of this range; the
% three P20 mice registered before this setting existed keep LightSuite's 6.
atlas_extent_slices = 15;

% atlas ('ccf' for the adults; 'demba_p20', 'demba_p16' or any other age built,
% for the young). The young register to the DeMBA template of their age (decided
% 2 Sep 2026): the manual control points carry the correspondence, so the adult
% template's better contrast (1.4x the global CV, 1.6x the local smoothed
% gradient), which feeds only the image-similarity terms, matters less than a
% target with the brain's own proportions. The adults stay on 'ccf' and are not
% registered again. The registered volumes come out at twice the registration
% grid, so the adults land on [900 800 1140] and the young on [994 800 1140]: code
% that reads both must take each cohort's grid from its atlas, and voxelwise work
% across them needs CCF Translator. Both annotations are in the same
% parcellation_index space, so regions compare directly.
atlas_key = 'demba_p22';

% preprocessing correction whose outputs 'align' reads ('slicewise'); the
% correction files hold their own correction_type, which replaces this one when
% they are loaded
correction_type = 'slicewise';

% use the equalised nano of run_nano_equalisation (1) or the corrected one (0);
% read by 'align' only
use_equalized_nano = 1;

% register a mouse without control points; keep it false for anything that ends
% up in a figure: LightSuite relies on the points to register well, so an
% image-only run is a diagnostic, not a result
allow_image_only_registration = false;

%% Run

% pass the settings to the code, under the same names
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
