close all
clear all
clc

% /// Preprocessing, step 2 of 6: extract the sections from the raw .czi and centre them ///
% For each selected mouse:
%   (1) Reads local_settings.txt (falls back to LightSuite internal defaults)
%   (2) Scans the .czi files and detects valid scenes (getSliceInfo)
%   (3) Extracts and centers all channels at px_process resolution, writing
%       volume_centered\chanXX_*.tiff plus volume_for_ordering.tiff
%   (4) Writes volume_ordered.tiff from the slice-ordering decisions file if
%       one exists, otherwise identity ordering
%
% Steps (1)-(3) are fully automatic. The MANUAL reorder/flip/discard step
% comes after, in run_order_slices (SliceOrderEditor), whose 'apply' mode
% writes volume_ordered.tiff from the saved decisions; a re-run of this
% script would apply them too.
%
% Mice come from the shared registry get_cohort.m rather than a hardcoded
% list, so every cohort (rws / naive / behavior / young) runs through the
% identical code path.
%
% Run sep_setup_paths first, once per MATLAB session. The settings are below,
% the code is in pipeline\extract_and_center.m.

%% User-defined parameters

% Cohort selection. Set mice_to_process to {} to process every mouse in
% groups_to_process; give explicit names to process just those.
groups_to_process = {'young'};                  % 'rws' | 'naive' | 'behavior' | 'young'
mice_to_process   = {'MG909_SepGluA_P20', 'MG910_SepGluA_P20', 'MG911_SepGluA_P16', ...
                     'MG912_SepGluA_P20', 'MG913_SepGluA_P20', 'MG914_SepGluA_P28'};
                                                % {} = all mice in groups_to_process
                                                % the first eight are already extracted

% Reference atlas (not used for extraction itself, only added to the path)
atlas_key = 'ccf';

%% Run

% The settings above go to the code under the same names
run_settings = struct();
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.atlas_key = atlas_key;
extract_and_center(run_settings);
