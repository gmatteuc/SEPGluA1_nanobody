close all
clear all
clc

% /// Preprocessing, step 1 of 6: copy the raw .czi from the lab share to local storage ///
% For each selected mouse:
%   (1) Resolves the source dir on the share, honouring the per-mouse
%       share_subdir in the registry (some brains keep their .czi under
%       Anatomy\Axioscan, others at the mouse root)
%   (2) Copies *.czi into <base_root>\<group>\<name>\ so every downstream
%       script sees ONE layout, with the .czi at the mouse root
%   (3) Verifies each file arrived with a byte-identical size
%
% WHY THIS EXISTS: the raw data on the share is READ-ONLY, and getSliceInfo
% creates its 'lightsuite' working folder NEXT TO the .czi it is given. Point
% run_extract_and_center at the share and it would write there. Copying first
% is mandatory, not stylistic.
%
% Transfer uses robocopy (restartable, resumes rather than restarts). No /MIR
% and no /MOV are ever passed, so the source cannot be modified.
%
% Run sep_setup_paths first, once per MATLAB session. The settings are below,
% the code is in pipeline\copy_raw_data.m.

%% User-defined parameters

% Cohort selection (mice come from the shared registry get_cohort.m).
% Set mice_to_process to {} to copy every mouse in groups_to_process.
groups_to_process = {'young'};
mice_to_process   = {'MG909_SepGluA_P20', 'MG910_SepGluA_P20', 'MG911_SepGluA_P16', ...
                     'MG912_SepGluA_P20', 'MG913_SepGluA_P20', 'MG914_SepGluA_P28'};

% Root of the raw data on the lab share (READ-ONLY - never written to)
share_root = 'S:\ElboustaniLab\#SHARE\Data';

% Set false for a dry run that reports what would be copied
do_copy = true;

%% Run

% The settings above go to the code under the same names
run_settings = struct();
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.share_root = share_root;
run_settings.do_copy = do_copy;
copy_raw_data(run_settings);
