close all
clear all
clc

% /// Preprocessing, step 3 of 6: MANUAL slice reorder / flip / discard ///
% Entry point for the one manual step between run_extract_and_center and
% run_residual_correction. Two modes:
%
%   run_mode = 'edit'   opens SliceOrderEditor on the selected mouse's
%                       volume_for_ordering.tiff. Reorder, flip and mark
%                       slices for removal, then save and close the GUI.
%                       It writes, next to that tiff:
%                         volume_for_ordering_processing_decisions.txt
%                       with columns OriginalIndex / FlipState / NewOrderOriginalIndex
%
%   run_mode = 'apply'  reloads the saved sliceinfo and rebuilds
%                       volume_ordered.tiff from the decisions file, in
%                       the folder that sliceinfo.mat names: never use it
%                       on a copied mouse folder (see pipeline\order_slices.m)
%
% Typical use: run with 'edit', curate, close the GUI, switch to 'apply',
% run again. Then continue with run_residual_correction.
%
% This exists so the manual step does not require re-running
% run_extract_and_center's ~10 min extraction just to reach the (previously
% commented-out) GUI call.
%
% Run sep_setup_paths first, once per MATLAB session. The settings are below,
% the code is in pipeline\order_slices.m.

%% User-defined parameters

% Cohort selection (mice come from the shared registry get_cohort.m).
% The GUI is per-mouse, so give exactly one name when run_mode = 'edit'.
%
% Whole young cohort below, with the number of sections in each. Work down the
% list: uncomment the one you are on, comment the previous.
%
% The P20 brains come first now: after the discussion with Sami the youngest
% ages are where the difference from adults is expected to be largest, so they
% are the ones that need to reach the registration stage first. The rest of the
% cohort follows, still worth curating but not on the critical path.

% --- P20 (and the P16 next to it) — priority ---
%  mice_to_process = {'MG897_SepGluA_P20'};   % 30 sections DONE
%  mice_to_process = {'MG903_SepGluA_P20'};   % 44 sections DONE
%  mice_to_process = {'MG909_SepGluA_P20'};   % 46 sections
%  mice_to_process = {'MG910_SepGluA_P20'};   % 46 sections
% mice_to_process = {'MG912_SepGluA_P20'};   % 30 sections
mice_to_process = {'MG913_SepGluA_P20'};   % 42 sections
% mice_to_process = {'MG911_SepGluA_P16'};   % 50 sections

% --- older ages — after the P20s ---
% mice_to_process = {'MG904_SepGluA_P22'};   % 46 sections DONE
% mice_to_process = {'MG896_SepGluA_P28'};   % 43 sections DONE
% mice_to_process = {'MG906_SepGluA_P32'};   % 30 sections DONE
% mice_to_process = {'MG895_SepGluA_P36'};   % 40 sections DONE
% mice_to_process = {'MG914_SepGluA_P28'};   % 46 sections
% mice_to_process = {'MG908_SepGluA_P32'};   % 37 sections
% mice_to_process = {'MG907_SepGluA_P36'};   % 33 sections

% 'edit' = open the GUI, 'apply' = rebuild volume_ordered.tiff from decisions
run_mode = 'edit';

%% Run

% The settings above go to the code under the same names
run_settings = struct();
run_settings.mice_to_process = mice_to_process;
run_settings.run_mode = run_mode;
order_slices(run_settings);
