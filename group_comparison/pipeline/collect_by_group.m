function collect_by_group(run_settings)
%COLLECT_BY_GROUP  Stack the registered volumes of each group's mice.
%   COLLECT_BY_GROUP(run_settings) does the work of run_collect_by_group,
%   which sets the fields of run_settings (paths, mousetypes_list,
%   age_filter, skip_missing) and says what each one does.
%
%   Saves, in each group's folder, nano_4d<tag>.mat (AP x DV x ML x mouse) and
%   collected_mice<tag>.mat (the mice along the fourth dimension, in order);
%   <tag> is empty for a whole group, _P20 for age_filter = [20].
%
%   Only nano_4d is saved: the saves of auto_4d, mask_4d and the averages are
%   commented out, so a rerun leaves those files as they were. The adults'
%   auto_4d.mat (19 Nov 2025) comes from the registration before the one behind
%   their nano_4d.mat (27 Nov), about two 10 um voxels away (r 0.94-0.99 plane
%   by plane), so a ratio of the two would mix two registrations. Uncomment a
%   save before reading its file, or read the registered tiffs, as
%   mapping/sepmap/volumes/per_mouse.py does.

% settings of run_collect_by_group, under the names the code below uses
paths = run_settings.paths;
mousetypes_list = run_settings.mousetypes_list;
age_filter = run_settings.age_filter;
skip_missing = run_settings.skip_missing;

%% Add paths

% the toolboxes are on the path from sep_setup_paths; only the atlas folder is
% added here, although nothing below reads the atlas
allenDir = paths.atlas;
addpath(allenDir)

%% Collect each group

% check that the registry still lists the adults in their legacy order
get_cohort('verify');

for mousetype_idx = 1:numel(mousetypes_list)

    % get the group's mice from the registry
    current_mouse_type = mousetypes_list{mousetype_idx};
    group_cohort = get_cohort('groups', {current_mouse_type});

    % keep the requested ages in registry order, so the 4th dimension follows the
    % cohort; the file names get _P20, or _P16P20P22 for several ages
    if isempty(age_filter)
        current_mice = {group_cohort.name};
        subset_tag = '';
    else
        keep = ismember([group_cohort.age_days], age_filter);
        current_mice = {group_cohort(keep).name};
        subset_tag = ['_P' strjoin(arrayfun(@(a) num2str(a), sort(age_filter), ...
            'UniformOutput', false), 'P')];
    end

    % drop the mice not registered yet, so a partly processed cohort still gives
    % a stack of whatever is ready
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

    % go to the next group if no mouse is ready
    num_current = numel(current_mice);
    if num_current == 0
        warning('run_collect_by_group: no registered mice for group %s%s, nothing to collect.', ...
            current_mouse_type, subset_tag);
        continue
    end
    fprintf('run_collect_by_group: %s%s — collecting %d mouse/mice: %s\n', ...
        current_mouse_type, subset_tag, num_current, strjoin(current_mice, ', '));

    % one cell per mouse, for each channel
    nanoVols_type = cell(1, num_current);
    autoVols_type = cell(1, num_current);
    maskVols_type = cell(1, num_current);

    for i = 1:num_current

        % get the mouse
        mouse_name = current_mice{i};

        % set its folders (base_dir is the group's)
        base_dir = fullfile(paths.data, current_mouse_type);
        mouse_dir = fullfile(base_dir, mouse_name, 'lightsuite');
        correction_dir = fullfile(base_dir, mouse_name, 'lightsuite', ...
            'correction_output');
        before_correction_dir = fullfile(base_dir, mouse_name, 'lightsuite', ...
            'volume_centered');
        processed_dir = fullfile(base_dir, mouse_name, 'lightsuite', ...
            'volume_centered_processed');
        aligned_dir = fullfile(mouse_dir, 'volume_aligned');
        registered_dir = fullfile(mouse_dir, 'volume_registered');

        % load its registered nano, autofluorescence and mask volumes
        file1_name = fullfile(registered_dir, sprintf('chan02_NANO.tiff'));
        file2_name = fullfile(registered_dir, sprintf('chan03_AUTO.tiff'));
        file4_name = fullfile(registered_dir, sprintf('chan05_MASK.tiff'));
        nanoVol = single(loadVolume({file1_name}, 1));
        autoVol = single(loadVolume({file2_name}, 1));
        maskVol = single(loadVolume({file4_name}, 1));

        % store them
        nanoVols_type{i} = nanoVol;
        autoVols_type{i} = autoVol;
        maskVols_type{i} = maskVol;

    end

    % stack the mice along the 4th dimension; a cohort's volumes share one grid,
    % set by get_atlas_crop (CCF 900 x 800 x 1140, DeMBA P20 994 x 800 x 1140)
    nano_4d = cat(4, nanoVols_type{:});
    auto_4d = cat(4, autoVols_type{:});
    mask_4d = cat(4, maskVols_type{:});

    % free the per-mouse copies
    clear nanoVols_type autoVols_type maskVols_type

    % average over mice
    avg_nano = nanmean(nano_4d, 4); %#ok<NANMEAN>
    avg_auto = nanmean(auto_4d, 4); %#ok<NANMEAN>
    sum_mask_4d = nansum(mask_4d, 4); %#ok<NANSUM>

    % mask of the voxels that at most min_num_contrib mice cover, for display
    min_num_contrib = 3;
    avg_mask = sum_mask_4d <= min_num_contrib;

    % nano relative to autofluorescence, (nano - auto) / auto, per mouse and averaged
    diff_4d_new = (nano_4d - auto_4d) ./ auto_4d;
    avg_diff_new = nanmean(diff_4d_new, 4); %#ok<NANMEAN>

    % save the nano stack (-v7.3 for large arrays), as nano_4d_P20.mat beside
    % nano_4d.mat for an age-filtered run; the other saves are off (reason not recorded)
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

    % save the mice along the 4th dimension, in order: without it a stack of a
    % partly registered cohort cannot be labelled later (get_cohort_spec reads it)
    collected_mice = current_mice;
    save(fullfile(base_dir, ['collected_mice' subset_tag '.mat']), 'collected_mice');

    % free memory before the next group
    clear nano_4d auto_4d mask_4d diff_4d_new avg_nano avg_auto avg_mask sum_mask_4d ...
        avg_diff_new

end

end
