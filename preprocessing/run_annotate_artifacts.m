%% run_annotate_artifacts
% ===== Outline the artifacts of each mouse by hand =====
%
% Preprocessing, step 6 of 6:
%   1. run_copy_raw_data        copy the .czi files from the lab share
%   2. run_extract_and_center   extract the sections and centre them
%   3. run_order_slices         reorder, flip and discard slices by hand
%   4. run_residual_correction  scale the autofluorescence to the nano
%   5. run_nano_equalisation    equalise nano intensity across slices
%   6. run_annotate_artifacts   outline the artifacts by hand  <- this script
%
% For each selected mouse, opens ArtifactAnnotator on what
% run_residual_correction saved (lightsuite\correction_output\
% scaled_auto_volume_<correction_type>.mat): the scaled autofluorescence in
% red, the nano in green, the background greyed. Draw a polygon around each
% artifact (d), move between slices with the left and right arrows, delete the
% last polygon with delete or the one under the cursor with a right click, save
% with s; closing the window (or escape) asks whether to save. The masks and
% polygons go to correction_output\artifact_mask_volume_<correction_type>.mat,
% which run_register_to_atlas reads and the window reloads next time, to carry
% on. The window holds MATLAB until it is closed; then the next mouse opens. A
% mouse without the correction output is skipped with a note.
%
% Setup: MG903. Run sep_setup_paths first, once per MATLAB session; the code is
% in pipeline\annotate_artifacts.m.

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% mice to annotate, from the cohort registry (get_cohort): the groups ('rws',
% 'naive', 'behavior', 'young'), or the mice named, which take precedence ({} =
% the groups); their windows open one after the other
groups_to_process = {'young'};
mice_to_process = {'MG903_SepGluA_P20'};

% the correction whose output is annotated ('slicewise' or 'global')
correction_type = 'slicewise';

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.correction_type = correction_type;
annotate_artifacts(run_settings);
