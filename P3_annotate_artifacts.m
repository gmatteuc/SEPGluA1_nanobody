clear all
close all
clc

% /// Pipeline script #3: display outputs of residual correction anlaysis in a GUI where user can annotate artifact for future removal /// 

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Cohort selection (mice come from the shared registry get_cohort.m).
% Set mice_to_process to {} to process every mouse in groups_to_process.
groups_to_process = {'young'};                  % 'rws' | 'naive' | 'behavior' | 'young'
mice_to_process   = {'MG903_SepGluA_P20'};      % {} = all mice in groups_to_process

% Choose correction type
correction_type = 'slicewise';

%% Resolve cohort

get_cohort('verify');
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end
fprintf('P3: %d mouse/mice selected.\n', numel(cohort));

%% Loop over mice

for mouse_idx = 1:numel(cohort)

    % Get current mouse name and type
    mouse_name = cohort(mouse_idx).name;
    mouse_type = cohort(mouse_idx).group;

    % Get dirs
    % NOTE: the 'lightsuite' level was missing here before, so output_dir did
    % not match where P2 writes and P4 reads. As written, P3 could never find
    % scaled_auto_volume_*.mat and silently skipped every mouse.
    base_dir = fullfile(paths.data, mouse_type);
    output_dir = fullfile(base_dir, mouse_name, 'lightsuite', 'correction_output');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Load the saved matfile
    matfile_name = fullfile(output_dir, sprintf('scaled_auto_volume_%s.mat', correction_type));
    if ~exist(matfile_name, 'file')
        fprintf('File not found for %s: %s\n', mouse_name, matfile_name);
        continue;
    end
    load(matfile_name);  % Loads scaledautoVol, nanoVol, bg_mask_vol, slice_data, average_slope, average_intercept, correction_type
    [H, W, Z] = size(nanoVol);

    % Invoke the artifact annotation GUI
    ArtifactAnnotator(scaledautoVol, nanoVol, bg_mask_vol, slice_data, mouse_name, output_dir, correction_type);

    % Clear per-mouse variables to save memory
    clear scaledautoVol nanoVol bg_mask_vol slice_data average_slope average_intercept

end