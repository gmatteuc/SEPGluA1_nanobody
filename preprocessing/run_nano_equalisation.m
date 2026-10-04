%% run_nano_equalisation
% ===== Equalise the nano intensity across each brain's slices =====
%
% Preprocessing, step 5 of 6:
%   1. run_copy_raw_data        copy the .czi files from the lab share
%   2. run_extract_and_center   extract the sections and centre them
%   3. run_order_slices         reorder, flip and discard slices by hand
%   4. run_residual_correction  scale the autofluorescence to the nano
%   5. run_nano_equalisation    equalise nano intensity across slices  <- this script
%   6. run_annotate_artifacts   outline the artifacts by hand
%
% The selected mice are processed together:
%   1. loads every mouse's centred nano volume
%      (lightsuite\volume_centered\chan02_Cy5.tiff) into one array, padded to
%      the largest height, width and number of slices
%   2. per slice, the median of the tissue pixels (the background removed with
%      select_background_pixels)
%   3. scales each slice so that its median becomes the moving median of the
%      slice medians over 5 slices
%   4. measures the medians again on the equalised volumes
%   5. saves each mouse's equalised volume in
%      lightsuite\correction_output\equalized_volume.mat, which
%      run_register_to_atlas reads (use_equalized_nano = 1)
% The statistics, profiles and heatmaps before and after equalisation, and a
% video of each mouse's slices before and after, go to base_output_dir
% (intensity_diagnostics\), named with the time of the run; save_results switches
% the statistics and figures, not the videos. The padding is filled with each
% slice's mode before the background is selected, so a brain's result can depend
% on which mice run with it.
%
% Setup: the young brains MG909 to MG914, together. Run sep_setup_paths first,
% once per MATLAB session; the code is in pipeline\nano_equalisation.m.

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% mice to equalise, from the cohort registry (get_cohort): the groups ('rws',
% 'naive', 'behavior', 'young'), or the mice named, which take precedence ({} =
% the groups)
groups_to_process = {'young'};
mice_to_process = {'MG909_SepGluA_P20', 'MG910_SepGluA_P20', 'MG911_SepGluA_P16', ...
    'MG912_SepGluA_P20', 'MG913_SepGluA_P20', 'MG914_SepGluA_P28'};

% atlas ('ccf'); only put on the path
atlas_key = 'ccf';

% save the statistics and their figures (the videos and the equalised volumes are
% written either way)
save_results = true;

% folder of the statistics, figures and videos
base_output_dir = fullfile(paths.data, 'intensity_diagnostics');

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.atlas_key = atlas_key;
run_settings.save_results = save_results;
run_settings.base_output_dir = base_output_dir;
nano_equalisation(run_settings);
