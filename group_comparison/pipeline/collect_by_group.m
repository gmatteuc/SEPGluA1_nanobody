function collect_by_group(run_settings)
%COLLECT_BY_GROUP  Stack the registered volumes of each group's mice.
%   COLLECT_BY_GROUP(run_settings) does the work of run_collect_by_group,
%   which sets the fields of run_settings (paths, mousetypes_list,
%   age_filter, skip_missing) and says what each one does.

% The settings of run_collect_by_group, under the names the code below uses
paths = run_settings.paths;
mousetypes_list = run_settings.mousetypes_list;
age_filter = run_settings.age_filter;
skip_missing = run_settings.skip_missing;

%% Add paths

% The toolboxes are on the path from sep_setup_paths; only the atlas folder
% is added here.
allenDir = paths.atlas;
addpath(allenDir)

%% Atlas

% Nothing below uses the atlas -- the block that used to load and crop the adult
% annotation here fed nothing. What matters is that every registered volume of a
% cohort sits on the same grid, which is what makes cat(4, ...) legal: the
% adults on the CCF (900 x 800 x 1140), the P20 brains on DeMBA (994 x 800 x
% 1140). get_atlas_crop is where that grid lives, and run_normalise_groups and
% P8 take it from there.

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
            fprintf('run_collect_by_group: %s — skipping %d not-yet-registered mouse/mice: %s\n', ...
                current_mouse_type, nnz(~has_reg), strjoin(current_mice(~has_reg), ', '));
        end
        current_mice = current_mice(has_reg);
    end

    num_current = numel(current_mice);
    if num_current == 0
        warning('run_collect_by_group: no registered mice for group %s%s, nothing to collect.', ...
            current_mouse_type, subset_tag);
        continue
    end
    fprintf('run_collect_by_group: %s%s — collecting %d mouse/mice: %s\n', ...
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
    %
    % CAREFUL, and the reason the v2 route stopped reading these files at all.
    % Only nano_4d is still written; every line below it is commented out. So a
    % re-run refreshes nano and leaves whatever auto_4d.mat, mask_4d.mat and the
    % averages happen to be on disk, however old. That is what happened to the
    % adults: the brains were re-registered and run_collect_by_group re-run on
    % 2025-11-27, which rewrote nano_4d.mat, while auto_4d.mat stayed at its
    % 2025-11-19 version, from the registration before. Nano and auto on disk
    % then belonged to two different registrations of the same brain, about two
    % 10 um voxels apart (r 0.94-0.99 plane by plane, every adult), and anything
    % dividing one by the other was quietly mixing them. Uncomment the saves
    % before trusting a file here, or read the registered tiffs directly, as
    % v2_per_mouse does.
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

end
