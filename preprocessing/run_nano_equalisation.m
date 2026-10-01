clear all
close all
clc

% /// Preprocessing, step 5 of 6: equalise the nano intensity across each brain's slices ///
% The selected mice are processed together:
%   (1) loads every mouse's centered nano volume
%       (lightsuite\volume_centered\chan02_Cy5.tiff) into one array, padded
%       to the largest height, width and number of slices
%   (2) per slice, the median of the tissue pixels (the background removed
%       with select_background_pixels)
%   (3) scales each slice so that its median becomes the moving median of
%       the slice medians over 5 slices
%   (4) measures the medians again on the equalised volumes
%   (5) saves each mouse's equalised volume in
%       lightsuite\correction_output\equalized_volume.mat, which
%       run_register_to_atlas reads (use_equalized_nano = 1)
% The statistics, profiles and heatmaps before and after equalisation, and a
% video of each mouse's slices before and after, go to base_output_dir
% (intensity_diagnostics\), named with the time of the run (the videos are
% written whatever save_results says, and with it off the run stops at the
% first video, whose name takes the time the saving sets). The padding is
% filled with each slice's mode before the background is selected, so a
% brain's result can depend on which mice run with it.
%
% Run sep_setup_paths first, once per MATLAB session. The settings are below,
% the code is in pipeline\nano_equalisation.m.

%% 1. User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Cohort selection (mice come from the shared registry get_cohort.m).
% Set mice_to_process to {} to process every mouse in groups_to_process.
% NOTE: this used to be a list of numeric INDICES into a hardcoded mouse
% list (mice_to_process = 1:17); it is now a list of mouse NAMES, so the
% selection no longer depends on the order of the registry.
groups_to_process = {'young'};                  % 'rws' | 'naive' | 'behavior' | 'young'
mice_to_process   = {'MG909_SepGluA_P20', 'MG910_SepGluA_P20', 'MG911_SepGluA_P16', ...
                     'MG912_SepGluA_P20', 'MG913_SepGluA_P20', 'MG914_SepGluA_P28'};                         % {} = all mice in groups_to_process

% Reference atlas
atlas_key = 'ccf';

% Output settings
save_results = true;
base_output_dir = fullfile(paths.data, 'intensity_diagnostics');

%% Run

% The settings above go to the code under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.atlas_key = atlas_key;
run_settings.save_results = save_results;
run_settings.base_output_dir = base_output_dir;
nano_equalisation(run_settings);
