clear all
close all
clc

% /// Preprocessing, step 6 of 6: MANUAL annotation of artifacts ///
% For each selected mouse, opens ArtifactAnnotator on what
% run_residual_correction saved (lightsuite\correction_output\
% scaled_auto_volume_<correction_type>.mat): the scaled autofluorescence in
% red, the nano in green, the background greyed. Draw a polygon around each
% artifact (d), move between slices with the left and right arrows, delete
% the last polygon with delete or the one under the cursor with a right
% click, save with s; closing the window (or escape) asks whether to save.
% The masks and polygons go to
% correction_output\artifact_mask_volume_<correction_type>.mat, which
% run_register_to_atlas reads and the window reloads next time, to carry
% on. The window holds MATLAB until it is closed; then the next mouse opens.
% A mouse without the correction output is skipped with a note.
%
% Run sep_setup_paths first, once per MATLAB session. The settings are below,
% the code is in pipeline\annotate_artifacts.m.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Cohort selection (mice come from the shared registry get_cohort.m).
% Set mice_to_process to {} to process every mouse in groups_to_process.
groups_to_process = {'young'};                  % 'rws' | 'naive' | 'behavior' | 'young'
mice_to_process   = {'MG903_SepGluA_P20'};      % {} = all mice in groups_to_process

% Choose correction type
correction_type = 'slicewise';

%% Run

% The settings above go to the code under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.correction_type = correction_type;
annotate_artifacts(run_settings);
