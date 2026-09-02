close all
clear all
clc

% /// Pipeline script #1bis: MANUAL slice reorder / flip / discard ///
% Entry point for the one manual step between P1 and P2. Two modes:
%
%   run_mode = 'edit'   opens SliceOrderEditor on the selected mouse's
%                       volume_for_ordering.tiff. Reorder, flip and mark
%                       slices for removal, then save and close the GUI.
%                       It writes, next to that tiff:
%                         volume_for_ordering_processing_decisions.txt
%                       with columns OriginalIndex / FlipState / NewOrderOriginalIndex
%
%   run_mode = 'apply'  reloads the saved sliceinfo and rebuilds
%                       volume_ordered.tiff from the decisions file
%
% Typical use: run with 'edit', curate, close the GUI, switch to 'apply',
% run again. Then continue with P2.
%
% This exists so the manual step does not require re-running P1's ~10 min
% extraction just to reach the (previously commented-out) GUI call.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

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
% mice_to_process = {'MG897_SepGluA_P20'};   % 30 sections DONE
% mice_to_process = {'MG903_SepGluA_P20'};   % 44 sections DONE
mice_to_process = {'MG909_SepGluA_P20'};   % 46 sections
% mice_to_process = {'MG910_SepGluA_P20'};   % 46 sections
% mice_to_process = {'MG912_SepGluA_P20'};   % 30 sections
% mice_to_process = {'MG913_SepGluA_P20'};   % 42 sections
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

%% Add paths

lightsuiteDir = paths.lightsuite;
yamlDir = paths.yaml;
addpath(genpath(lightsuiteDir))
addpath(genpath(yamlDir))

%% Resolve cohort

get_cohort('verify');
cohort = get_cohort('names', mice_to_process);

%% Run the selected mode

for mouse_idx = 1:numel(cohort)

    mousename = cohort(mouse_idx).name;
    procpath  = fullfile(cohort(mouse_idx).base_dir, 'lightsuite');
    volorder  = fullfile(procpath, 'volume_for_ordering.tiff');
    decisions = fullfile(procpath, 'volume_for_ordering_processing_decisions.txt');

    fprintf('\n=== %s (group %s) ===\n', mousename, cohort(mouse_idx).group);

    if ~exist(volorder, 'file')
        error('Ordering volume not found:\n  %s\nRun P1 for this mouse first.', volorder);
    end

    switch lower(run_mode)

        case 'edit'
            if numel(cohort) > 1
                error('run_mode ''edit'' opens one GUI at a time; select a single mouse.');
            end
            if exist(decisions, 'file')
                fprintf('  NOTE: a decisions file already exists and will be overwritten on save:\n    %s\n', decisions);
            end
            fprintf('  opening SliceOrderEditor on:\n    %s\n', volorder);
            fprintf('  curate, then SAVE and CLOSE the GUI, set run_mode = ''apply'' and re-run.\n');
            SliceOrderEditor(volorder);

        case 'apply'
            if ~exist(decisions, 'file')
                error(['No decisions file found:\n  %s\n' ...
                       'Run this script with run_mode = ''edit'' first.'], decisions);
            end
            sliceinfo_name = fullfile(procpath, 'sliceinfo.mat');
            if ~exist(sliceinfo_name, 'file')
                error('sliceinfo.mat not found:\n  %s\nRun P1 for this mouse first.', sliceinfo_name);
            end
            S = load(sliceinfo_name);
            sliceinfo = S.sliceinfo;

            T = readtable(decisions);
            fprintf('  decisions: %d slices, %d flipped, %d reordered\n', ...
                height(T), sum(T.FlipState == 1), sum(T.NewOrderOriginalIndex(:)' ~= 1:height(T)));

            generateReordedVolume(sliceinfo);
            fprintf('  rebuilt: %s\n', fullfile(procpath, 'volume_ordered.tiff'));

        otherwise
            error('Unknown run_mode: %s (use ''edit'' or ''apply'').', run_mode);
    end

end
