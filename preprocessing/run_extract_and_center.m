%% run_extract_and_center
% ===== Extract the sections from the raw .czi files and centre them =====
%
% Preprocessing, step 2 of 6:
%   1. run_copy_raw_data        copy the .czi files from the lab share
%   2. run_extract_and_center   extract the sections and centre them  <- this script
%   3. run_order_slices         reorder, flip and discard slices by hand
%   4. run_residual_correction  scale the autofluorescence to the nano
%   5. run_nano_equalisation    equalise nano intensity across slices
%   6. run_annotate_artifacts   outline the artifacts by hand
%
% For each selected mouse, reads the extraction settings from local_settings.txt
% (the mouse folder first, then lightsuite\; without either, LightSuite's
% defaults, which the adults ran on), finds the sections in the .czi files
% (getSliceInfo), and writes every channel, centred, at px_process resolution
% to lightsuite\volume_centered\chanXX_<dye>.tiff, with sliceinfo.mat and
% volume_for_ordering.tiff, the colour composite that run_order_slices shows.
% Then it writes volume_ordered.tiff, from the decisions file of
% run_order_slices when there is one, otherwise in the extracted order, so a
% rerun after the curation applies it too. A mouse that fails is reported and
% the others go on. The extraction is automatic; the order is curated by hand
% next, in run_order_slices.
%
% local_settings.txt also sets px_atlas, which sliceinfo.mat keeps but which
% does not decide the grid of the registered volume: run_register_to_atlas puts
% every brain on the 10 um grid (registered_grid_um), provided px_register is 20.
%
% Setup: the young brains MG909 to MG914. Run sep_setup_paths first, once per
% MATLAB session; the code is in pipeline\extract_and_center.m.

clear; clc; close all;

%% Settings

% mice to extract, from the cohort registry (get_cohort): the groups ('rws',
% 'naive', 'behavior', 'young'), or the mice named, which take precedence ({} =
% the groups); the first eight young brains are extracted already
groups_to_process = {'young'};
mice_to_process = {'MG909_SepGluA_P20', 'MG910_SepGluA_P20', 'MG911_SepGluA_P16', ...
    'MG912_SepGluA_P20', 'MG913_SepGluA_P20', 'MG914_SepGluA_P28'};

% atlas ('ccf'); only put on the path, the extraction does not use it
atlas_key = 'ccf';

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.atlas_key = atlas_key;
extract_and_center(run_settings);
