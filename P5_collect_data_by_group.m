clear all
close all
clc

% /// Pipeline script #5: load output of registartion for all mice and assemble 4D volumes stording data from all mice of each group /// 

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Cohort selection (mice come from the shared registry get_cohort.m).
mousetypes_list = {'young'};                    % 'rws' | 'naive' | 'behavior' | 'young'

% Optional age filter, applied within each group above. Leave empty to take the
% whole group. The young cohort spans P16-P36, but the comparison Sami wants
% first is the youngest ages against the adults, so a P20-only aggregate is
% assembled separately rather than diluting it with the P32/P36 brains.
age_filter = [20];                              % [] = whole group, e.g. [20] or [16 20 22]

% Only mice that actually reached the end of P4 can be collected here; the rest
% are skipped with a warning rather than killing the run, so the aggregate can
% be rebuilt as more brains finish registering.
skip_missing = true;

% Choose correction type
correction_type = 'slicewise';

%% Add paths

% Define paths
allenDir = paths.atlas;
lightsuiteDir = paths.lightsuite;
yamlDir = paths.yaml;
elastixDir = paths.elastix;
% Add defined paths
addpath(allenDir)
addpath(genpath(lightsuiteDir))
addpath(genpath(yamlDir))
addpath(genpath(elastixDir))

%% Get atlas

% Load Allen atlas annotations
AllenFile = [allenDir,'\annotation_10.nii.gz'];
AllenVol = niftiread(AllenFile);
limits = [180 1079];
AllenCrop = AllenVol(limits(1):limits(2),:,:);
brainMask = AllenCrop > 0;

%% Process per mouse type

get_cohort('verify');

for mousetype_idx = 1:numel(mousetypes_list)

    % Get current mouse type and select mice
    current_mouse_type = mousetypes_list{mousetype_idx};
    group_cohort = get_cohort('groups', {current_mouse_type});

    % Narrow to the requested ages, if any. Registry order is preserved so the
    % 4th dimension of the saved volumes stays in the same order as the cohort.
    if isempty(age_filter)
        current_mice = {group_cohort.name};
        subset_tag = '';
    else
        keep = ismember([group_cohort.age_days], age_filter);
        current_mice = {group_cohort(keep).name};
        subset_tag = ['_P' strjoin(arrayfun(@(a) num2str(a), sort(age_filter), ...
                      'UniformOutput', false), 'P')];
    end

    % Drop mice that have not been registered yet, so a partially processed
    % cohort still yields an aggregate of whatever is ready.
    if skip_missing
        has_reg = cellfun(@(m) exist(fullfile(paths.data, current_mouse_type, m, ...
                  'lightsuite', 'volume_registered', 'chan02_NANO.tiff'), 'file') == 2, ...
                  current_mice);
        if any(~has_reg)
            fprintf('P5: %s — skipping %d not-yet-registered mouse/mice: %s\n', ...
                current_mouse_type, nnz(~has_reg), strjoin(current_mice(~has_reg), ', '));
        end
        current_mice = current_mice(has_reg);
    end

    num_current = numel(current_mice);
    if num_current == 0
        warning('P5: no registered mice for group %s%s, nothing to collect.', ...
            current_mouse_type, subset_tag);
        continue
    end
    fprintf('P5: %s%s — collecting %d mouse/mice: %s\n', ...
        current_mouse_type, subset_tag, num_current, strjoin(current_mice, ', '));

    % Initialize cells for volumes of this type
    nanoVols_type = cell(1, num_current);
    autoVols_type = cell(1, num_current);
    maskVols_type = cell(1, num_current);

    for i = 1:num_current

        % Get current mouse
        mouse_name = current_mice{i};

        % Get dirs (base_dir is per type)
        base_dir = fullfile(paths.data, current_mouse_type);
        mouse_dir = fullfile(base_dir, mouse_name, 'lightsuite');
        correction_dir = fullfile(base_dir, mouse_name, 'lightsuite', 'correction_output');
        before_correction_dir = fullfile(base_dir, mouse_name, 'lightsuite', 'volume_centered');
        processed_dir = fullfile(base_dir, mouse_name, 'lightsuite', 'volume_centered_processed');
        aligned_dir = fullfile(mouse_dir, 'volume_aligned');
        registered_dir = fullfile(mouse_dir, 'volume_registered');

        % Load needed data for current mouse
        file1_name = fullfile(registered_dir, sprintf('chan02_NANO.tiff'));
        file2_name = fullfile(registered_dir, sprintf('chan03_AUTO.tiff'));
        file4_name = fullfile(registered_dir, sprintf('chan05_MASK.tiff'));
        nanoVol = single(loadVolume({file1_name}, 1));
        autoVol = single(loadVolume({file2_name}, 1));
        maskVol = single(loadVolume({file4_name}, 1));

        % Store in cells
        nanoVols_type{i} = nanoVol;
        autoVols_type{i} = autoVol;
        maskVols_type{i} = maskVol;

    end

    % Concatenate volumes along the 4th dimension
    nano_4d = cat(4, nanoVols_type{:});
    auto_4d = cat(4, autoVols_type{:});
    mask_4d = cat(4, maskVols_type{:});

    % Clear temporary cells to free memory
    clear nanoVols_type autoVols_type maskVols_type

    % Compute averages along the 4th dimension
    avg_nano = nanmean(nano_4d, 4); %#ok<NANMEAN>
    avg_auto = nanmean(auto_4d, 4); %#ok<NANMEAN>
    sum_mask_4d = nansum(mask_4d, 4); %#ok<NANSUM>

    % Form derived mask for display
    min_num_contrib=3;
    avg_mask=sum_mask_4d<=min_num_contrib;

    % Compute derived diff volumes
    diff_4d_new = (nano_4d-auto_4d)./auto_4d;
    avg_diff_new = nanmean(diff_4d_new, 4); %#ok<NANMEAN>

    % Save 4D volumes and averages in base_dir (use -v7.3 for large matrices).
    % subset_tag is empty for a whole group, so the adult cohorts keep writing
    % plain nano_4d.mat exactly as before; an age-filtered run lands beside it
    % as e.g. nano_4d_P20.mat instead of overwriting it.
    base_dir = fullfile(paths.data, current_mouse_type);
    save(fullfile(base_dir, ['nano_4d' subset_tag '.mat']), 'nano_4d', '-v7.3');
    % save(fullfile(base_dir, ['auto_4d' subset_tag '.mat']), 'auto_4d', '-v7.3');
    % save(fullfile(base_dir, ['mask_4d' subset_tag '.mat']), 'mask_4d', '-v7.3');
    % save(fullfile(base_dir, ['diff_4d_new' subset_tag '.mat']), 'diff_4d_new', '-v7.3');
    % save(fullfile(base_dir, ['avg_nano' subset_tag '.mat']), 'avg_nano', '-v7.3');
    % save(fullfile(base_dir, ['avg_auto' subset_tag '.mat']), 'avg_auto', '-v7.3');
    % save(fullfile(base_dir, ['avg_mask' subset_tag '.mat']), 'avg_mask', '-v7.3');
    % save(fullfile(base_dir, ['sum_mask_4d' subset_tag '.mat']), 'sum_mask_4d', '-v7.3');
    % save(fullfile(base_dir, ['avg_diff_new' subset_tag '.mat']), 'avg_diff_new', '-v7.3');

    % Which mice ended up in the 4th dimension, and in what order. Without this
    % an aggregate built from a partially registered cohort is unreadable later,
    % and the per-mouse work in P8 has no way to label its slices.
    collected_mice = current_mice;
    save(fullfile(base_dir, ['collected_mice' subset_tag '.mat']), 'collected_mice');

    % Clear variables to free memory before next type
    clear nano_4d auto_4d mask_4d diff_4d_new avg_nano avg_auto avg_mask sum_mask_4d avg_diff_new

end