%% run_residual_correction
% ===== Scale the autofluorescence onto the nano channel =====
%
% Preprocessing, step 4 of 6:
%   1. run_copy_raw_data        copy the .czi files from the lab share
%   2. run_extract_and_center   extract the sections and centre them
%   3. run_order_slices         reorder, flip and discard slices by hand
%   4. run_residual_correction  scale the autofluorescence to the nano  <- this script
%   5. run_nano_equalisation    equalise nano intensity across slices
%   6. run_annotate_artifacts   outline the artifacts by hand
%
% For each selected mouse, slice by slice, the autofluorescence channel (Cy3)
% is fitted to the nano channel (Cy5) on reference pixels, and the scaled
% autofluorescence is subtracted from the nano:
%   1. loads the centred nano and autofluorescence volumes
%      (lightsuite\volume_centered\chan02_Cy5.tiff and chan03_Cy3.tiff)
%   2. per slice, picks reference pixels and the background on the nano image
%      (select_reference_pixels) and fits nano against autofluorescence on the
%      reference pixels (robustfit, bisquare): a slope and an intercept
%   3. applies the fit two ways, 'slicewise' (each slice's own) and 'global'
%      (the mean over the slices): the scaled autofluorescence, and the nano
%      minus it, with negative values set to 0
%   4. saves, for each way, in lightsuite\correction_output\:
%      corrected_volume_<type>.mat (read by run_register_to_atlas),
%      scaled_auto_volume_<type>.mat (read by run_annotate_artifacts and
%      run_register_to_atlas), and a video of the relative difference,
%      (nano - scaled auto) / scaled auto
% The figures of each slice's fit go to correction_output\diagnostic_plots\,
% with the figures of the reference-pixel selection when savePlotBkg is on. A
% video of the nano / autofluorescence ratio is written when saveRatioMap is on.
%
% Setup: the young brains MG909 to MG914. Run sep_setup_paths first, once per
% MATLAB session; the code is in pipeline\residual_correction.m.

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% mice to correct, from the cohort registry (get_cohort): the groups ('rws',
% 'naive', 'behavior', 'young'), or the mice named, which take precedence ({} =
% the groups); each mouse is corrected on its own
groups_to_process = {'young'};
mice_to_process = {'MG909_SepGluA_P20', 'MG910_SepGluA_P20', 'MG911_SepGluA_P16', ...
    'MG912_SepGluA_P20', 'MG913_SepGluA_P20', 'MG914_SepGluA_P28'};

% atlas ('ccf'); only put on the path
atlas_key = 'ccf';

% draw the figures of the reference-pixel selection (one per slice)
doPlotBkg = true;

% save those figures, when drawn
savePlotBkg = true;

% write a video of the nano / autofluorescence ratio of each slice
saveRatioMap = false;

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.atlas_key = atlas_key;
run_settings.doPlotBkg = doPlotBkg;
run_settings.savePlotBkg = savePlotBkg;
run_settings.saveRatioMap = saveRatioMap;
residual_correction(run_settings);
