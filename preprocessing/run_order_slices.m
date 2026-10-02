%% run_order_slices
% ===== Reorder, flip and discard the slices of a mouse by hand =====
%
% Preprocessing, step 3 of 6:
%   1. run_copy_raw_data        copy the .czi files from the lab share
%   2. run_extract_and_center   extract the sections and centre them
%   3. run_order_slices         reorder, flip and discard slices by hand  <- this script
%   4. run_residual_correction  scale the autofluorescence to the nano
%   5. run_nano_equalisation    equalise nano intensity across slices
%   6. run_annotate_artifacts   outline the artifacts by hand
%
% The one manual step between run_extract_and_center and run_residual_correction,
% in two modes:
%   'edit'   opens SliceOrderEditor on the mouse's volume_for_ordering.tiff:
%            reorder, flip and mark slices for removal, then save and close the
%            window. It writes volume_for_ordering_processing_decisions.txt next
%            to the tiff, with the columns OriginalIndex, FlipState and
%            NewOrderOriginalIndex.
%   'apply'  rebuilds volume_ordered.tiff from the decisions file, in the folder
%            that sliceinfo.mat names: never on a copied mouse folder (see
%            pipeline\order_slices.m).
% Run with 'edit', curate, close the window, set 'apply' and run again; then
% run_residual_correction. A step of its own, so the curation does not need the
% ten-minute extraction of run_extract_and_center again.
%
% Setup: the young cohort, one mouse at a time. Run sep_setup_paths first, once
% per MATLAB session; the code is in pipeline\order_slices.m.

clear; clc; close all;

%% Settings

% the mouse to curate, by name (the window shows one mouse, so one name for
% 'edit'); below, the whole young cohort with the number of sections of each:
% work down the list, uncommenting the mouse being curated
%
% the P20 brains come first: the difference from the adults is expected to be
% largest at the youngest ages, so they are the first to need registering; the
% older ages follow, off the critical path

% P20, and the P16 next to it, first
% mice_to_process = {'MG897_SepGluA_P20'};   % 30 sections, done
% mice_to_process = {'MG903_SepGluA_P20'};   % 44 sections, done
% mice_to_process = {'MG909_SepGluA_P20'};   % 46 sections
% mice_to_process = {'MG910_SepGluA_P20'};   % 46 sections
% mice_to_process = {'MG912_SepGluA_P20'};   % 30 sections
mice_to_process = {'MG913_SepGluA_P20'};   % 42 sections
% mice_to_process = {'MG911_SepGluA_P16'};   % 50 sections

% the older ages, after the P20 brains
% mice_to_process = {'MG904_SepGluA_P22'};   % 46 sections, done
% mice_to_process = {'MG896_SepGluA_P28'};   % 43 sections, done
% mice_to_process = {'MG906_SepGluA_P32'};   % 30 sections, done
% mice_to_process = {'MG895_SepGluA_P36'};   % 40 sections, done
% mice_to_process = {'MG914_SepGluA_P28'};   % 46 sections
% mice_to_process = {'MG908_SepGluA_P32'};   % 37 sections
% mice_to_process = {'MG907_SepGluA_P36'};   % 33 sections

% 'edit' opens the window, 'apply' rebuilds volume_ordered.tiff from the decisions
run_mode = 'edit';

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.mice_to_process = mice_to_process;
run_settings.run_mode = run_mode;
order_slices(run_settings);
