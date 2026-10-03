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
%   Only the nano volumes are read and stacked. The auto_4d.mat, mask_4d.mat
%   and average files beside nano_4d.mat in the adults' group folders are not
%   written here, and are not to be read: the adults' auto_4d.mat (19 Nov
%   2025) comes from the registration before the one behind their nano_4d.mat
%   (27 Nov), about two 10 um voxels away (r 0.94-0.99 plane by plane), so a
%   ratio of the two would mix two registrations. For the autofluorescence or
%   the mask, read the registered tiffs, as mapping/sepmap/volumes/per_mouse.py
%   does.

% settings of run_collect_by_group, under the names the code below uses
paths = run_settings.paths;
mousetypes_list = run_settings.mousetypes_list;
age_filter = run_settings.age_filter;
skip_missing = run_settings.skip_missing;

%% Collect each group

% check that the registry still lists the adults in their legacy order
get_cohort('verify');

for group_idx = 1:numel(mousetypes_list)

    % get the group's mice from the registry
    group = mousetypes_list{group_idx};
    group_cohort = get_cohort('groups', {group});

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
        is_registered = cellfun(@(m) exist(fullfile(paths.data, group, m, ...
            'lightsuite', 'volume_registered', 'chan02_NANO.tiff'), 'file') == 2, ...
            current_mice);
        if any(~is_registered)
            fprintf('run_collect_by_group: %s — skipping %d not-yet-registered mouse/mice: %s\n', ...
                group, nnz(~is_registered), strjoin(current_mice(~is_registered), ', '));
        end
        current_mice = current_mice(is_registered);
    end

    % go to the next group if no mouse is ready
    n_mice = numel(current_mice);
    if n_mice == 0
        warning('run_collect_by_group: no registered mice for group %s%s, nothing to collect.', ...
            group, subset_tag);
        continue
    end
    fprintf('run_collect_by_group: %s%s — collecting %d mouse/mice: %s\n', ...
        group, subset_tag, n_mice, strjoin(current_mice, ', '));

    % one cell per mouse
    nano_vols = cell(1, n_mice);

    for i = 1:n_mice

        % load the mouse's registered nano volume
        mouse_name = current_mice{i};
        registered_dir = fullfile(paths.data, group, mouse_name, 'lightsuite', ...
            'volume_registered');
        nano_file = fullfile(registered_dir, sprintf('chan02_NANO.tiff'));
        nanoVol = single(loadVolume({nano_file}, 1));
        nano_vols{i} = nanoVol;

    end

    % stack the mice along the 4th dimension; a cohort's volumes share one grid,
    % set by get_atlas_crop (CCF 900 x 800 x 1140, DeMBA P20 994 x 800 x 1140)
    nano_4d = cat(4, nano_vols{:});

    % free the per-mouse copies
    clear nano_vols

    % save the nano stack (-v7.3 for large arrays), as nano_4d_P20.mat beside
    % nano_4d.mat for an age-filtered run
    base_dir = fullfile(paths.data, group);
    save(fullfile(base_dir, ['nano_4d' subset_tag '.mat']), 'nano_4d', '-v7.3');

    % save the mice along the 4th dimension, in order: without it a stack of a
    % partly registered cohort cannot be labelled later (get_cohort_spec reads it)
    collected_mice = current_mice;
    save(fullfile(base_dir, ['collected_mice' subset_tag '.mat']), 'collected_mice');

    % free memory before the next group
    clear nano_4d

end

end
