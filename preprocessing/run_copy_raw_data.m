%% run_copy_raw_data
% ===== Copy the raw .czi files from the lab share =====
%
% Preprocessing, step 1 of 6:
%   1. run_copy_raw_data        copy the .czi files from the lab share  <- this script
%   2. run_extract_and_center   extract the sections and centre them
%   3. run_order_slices         reorder, flip and discard slices by hand
%   4. run_residual_correction  scale the autofluorescence to the nano
%   5. run_nano_equalisation    equalise nano intensity across slices
%   6. run_annotate_artifacts   outline the artifacts by hand
%
% Copies the .czi files of each selected mouse from the lab share into
% <data>\<group>\<mouse>\, so every later step finds them in one layout, at the
% mouse root, and checks that each file arrived with the same size in bytes.
% The folder on the share comes from the cohort registry (get_cohort): some
% brains keep their .czi files under Anatomy\Axioscan, others at the mouse root.
%
% The copy comes first because the share is read only and getSliceInfo makes its
% lightsuite working folder next to the .czi files it is given: pointed at the
% share, run_extract_and_center would write there. robocopy copies (restartable,
% so an interrupted copy resumes), never with /MIR or /MOV, so the source cannot
% change, and a destination on the share's drive is refused.
%
% Setup: the young brains MG909 to MG914. Run sep_setup_paths first, once per
% MATLAB session; the code is in pipeline\copy_raw_data.m.

clear; clc; close all;

%% Settings

% mice to copy, from the cohort registry (get_cohort): the groups ('rws', 'naive',
% 'behavior', 'young'), or the mice named, which take precedence ({} = the groups)
groups_to_process = {'young'};
mice_to_process = {'MG909_SepGluA_P20', 'MG910_SepGluA_P20', 'MG911_SepGluA_P16', ...
    'MG912_SepGluA_P20', 'MG913_SepGluA_P20', 'MG914_SepGluA_P28'};

% root of the raw data on the lab share, read only: nothing is written there
share_root = 'S:\ElboustaniLab\#SHARE\Data';

% false for a dry run, which reports what would be copied
do_copy = true;

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.share_root = share_root;
run_settings.do_copy = do_copy;
copy_raw_data(run_settings);
