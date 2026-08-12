close all
clear all
clc

% /// Pipeline script #1: extracts data from raw .czi files and save centered volumes for further processing  ///
% For each selected mouse:
%   (1) Reads local_settings.txt (falls back to LightSuite internal defaults)
%   (2) Scans the .czi files and detects valid scenes (getSliceInfo)
%   (3) Extracts and centers all channels at px_process resolution, writing
%       volume_centered\chanXX_*.tiff plus volume_for_ordering.tiff
%   (4) Writes volume_ordered.tiff from the slice-ordering decisions file if
%       one exists, otherwise identity ordering
%
% Steps (1)-(3) are fully automatic. The MANUAL reorder/flip/discard step
% (SliceOrderEditor, commented at the bottom) comes after; re-running this
% script then applies the saved decisions.
%
% Mice come from the shared registry get_cohort.m rather than a hardcoded
% list, so every cohort (rws / naive / behavior / young) runs through the
% identical code path.

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

%% Add paths

% Define paths
atlas = get_atlas(atlas_key);
allenDir = atlas.dir;
lightsuiteDir = 'D:\sep_histology\code\LightSuite-main';
yamlDir = 'D:\sep_histology\code\yamlmatlab';
elastixDir = 'D:\sep_histology\code\matlab_elastix-master';
% Reader for the raw .czi (BioformatsImage class + bundled bfmatlab).
% Needed by getSliceInfo/generateSliceVolume; no P script added it before,
% so P1 only ever ran interactively with the path already set.
bioformatsDir = 'D:\sep_histology\code\BioformatsImage';
% Add defined paths
addpath(allenDir)
addpath(genpath(lightsuiteDir))
addpath(genpath(yamlDir))
addpath(genpath(elastixDir))
addpath(genpath(bioformatsDir))

%% Resolve cohort

get_cohort('verify');
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end
fprintf('P1: %d mouse/mice selected.\n', numel(cohort));

%% Loop over mice

% If one mouse fails we keep going, otherwise a whole unattended batch can be
% lost to a single bad brain. Any failures are listed again at the end.
failed_mice = {};

for mouse_idx = 1:numel(cohort)

    % Get current mouse metadata
    mousename = cohort(mouse_idx).name;
    dp        = cohort(mouse_idx).base_dir;
    fprintf('\n=== %s (group %s) ===\n%s\n', mousename, cohort(mouse_idx).group, dp);

    try

        %% Initialize extraction via Lightsuite

        if ~exist(dp, 'dir')
            error('Mouse dir not found: %s\nCopy the raw .czi from the lab share first.', dp);
        end

        % Read settings (mouse root first, then the lightsuite subdir where the
        % adult cohort keeps it; parseSettingsFile falls back to its internal
        % defaults if neither exists, which is what the adults actually ran on)
        settings_path = resolve_settings_path(dp);
        sliceinfo     = parseSettingsFile(settings_path);
        fprintf('  settings: %s\n', settings_path);
        fprintf('  slicethickness=%g px_process=%g px_register=%g px_atlas=%g regchan=%s\n', ...
            sliceinfo.slicethickness, sliceinfo.px_process, sliceinfo.px_register, ...
            sliceinfo.px_atlas, sliceinfo.regchan);

        sliceinfo.mousename   = mousename;
        filelistcheck         = dir(fullfile(dp, '*.czi'));
        if isempty(filelistcheck)
            error('No .czi found in %s\nCopy the raw files from the lab share first.', dp);
        end
        filepaths             = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');
        sliceinfo.filepaths   = filepaths;
        fprintf('  %d .czi file(s)\n', numel(filepaths));
        sliceinfo             = getSliceInfo(sliceinfo);

        %% Generate the slice volume (auto)

        % Generate centered volume from raw data
        slicevol = generateSliceVolume(sliceinfo, sliceinfo.regchan); %#ok<NASGU>

        %% Apply slice ordering if it has already been curated (auto)

        % Rebuilds volume_ordered.tiff from volume_for_ordering_processing_decisions.txt
        % when that file exists, otherwise writes identity ordering. Safe to run
        % before curation; re-run afterwards to apply the decisions.
        %
        % The MANUAL reorder/flip/discard step is NOT here: it lives in
        % P1bis_order_slices.m. It used to be a commented-out SliceOrderEditor call
        % at this spot, which was misleading -- the GUI is non-blocking, so inside
        % this loop it would open one window per mouse at once, and it sat *after*
        % generateReordedVolume so its output was never consumed in the same pass.
        generateReordedVolume(sliceinfo);

        fprintf('  done: %s\n', mousename);

    catch err
        fprintf('  FAILED (%s): %s\n', mousename, err.message);
        failed_mice{end+1} = mousename; %#ok<SAGROW>
    end

end

%% Report

if isempty(failed_mice)
    fprintf('\nP1 finished: all %d mouse/mice processed.\n', numel(cohort));
else
    fprintf('\nP1 finished with %d failure(s): %s\n', ...
        numel(failed_mice), strjoin(failed_mice, ', '));
end

%% Local function: locate local_settings.txt for a mouse

function settings_path = resolve_settings_path(mouse_dir)
% Returns the first existing local_settings.txt, checking the mouse root then
% the lightsuite subdir. Returns the mouse-root path when neither exists, so
% parseSettingsFile warns and applies its internal defaults (slicethickness=150,
% px_process=5, px_register=20, px_atlas=10, regchan='dapi') -- the values the
% adult cohort was processed with.

candidates = { ...
    fullfile(mouse_dir, 'local_settings.txt'), ...
    fullfile(mouse_dir, 'lightsuite', 'local_settings.txt')};

settings_path = candidates{1};
for k = 1:numel(candidates)
    if exist(candidates{k}, 'file')
        settings_path = candidates{k};
        return
    end
end
end
