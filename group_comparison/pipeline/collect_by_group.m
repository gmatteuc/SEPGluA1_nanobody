function collect_by_group(run_settings)
%COLLECT_BY_GROUP  Stack the registered volumes of each group's mice.
%   COLLECT_BY_GROUP(run_settings) does the work of run_collect_by_group,
%   which sets the fields of run_settings (paths, mousetypes_list,
%   age_filter, skip_missing, channels) and says what each one does.
%
%   Saves, in each group's folder, <channel>_4d<tag>.mat for each channel
%   (AP x DV x ML x mouse, with collected_mice, the mice along the fourth
%   dimension, and source_files, the registered volumes they were read from)
%   and, for the nano channel, collected_mice<tag>.mat (the nano stack's
%   mice, in order); <tag> is empty for a whole group, _P20 for
%   age_filter = [20].
%
%   Each channel is read from the mice's registered volumes,
%   lightsuite\volume_registered\chan02_NANO.tiff for nano and
%   chan03_AUTO.tiff for auto (the autofluorescence), written by one
%   registration, so the two stacks are on one grid voxel for voxel. The
%   auto_4d.mat, mask_4d.mat and average files of 19 November 2025 in the
%   adults' production group folders were not written here, and are not to be
%   read: their auto_4d.mat comes from the registration before the one behind
%   their nano_4d.mat (27 Nov), about two 10 um voxels away (r 0.94-0.99 plane
%   by plane), so a ratio of the two would mix two registrations. A stack
%   written here records its source_files, which that old file lacks, and
%   per_mouse_region_values refuses an auto_4d.mat without them.

% settings of run_collect_by_group, under the names the code below uses
paths = run_settings.paths;
mousetypes_list = run_settings.mousetypes_list;
age_filter = run_settings.age_filter;
skip_missing = run_settings.skip_missing;
channels = run_settings.channels;

% each channel's registered volume, by the channel's name
channel_files = struct('nano', 'chan02_NANO.tiff', 'auto', 'chan03_AUTO.tiff');
is_known = isfield(channel_files, channels);
if ~all(is_known)
    error('run_collect_by_group: unknown channel %s (use ''nano'' or ''auto'').', ...
        strjoin(channels(~is_known), ', '));
end

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
        group_mice = {group_cohort.name};
        subset_tag = '';
    else
        keep = ismember([group_cohort.age_days], age_filter);
        group_mice = {group_cohort(keep).name};
        subset_tag = ['_P' strjoin(arrayfun(@(a) num2str(a), sort(age_filter), ...
            'UniformOutput', false), 'P')];
    end
    base_dir = fullfile(paths.data, group);

    for channel_idx = 1:numel(channels)
        channel = channels{channel_idx};
        channel_file = channel_files.(channel);

        % drop the mice whose volume of the channel is not registered yet, so a
        % partly processed cohort still gives a stack of whatever is ready
        current_mice = group_mice;
        if skip_missing
            is_registered = cellfun(@(m) exist(fullfile(base_dir, m, 'lightsuite', ...
                'volume_registered', channel_file), 'file') == 2, current_mice);
            if any(~is_registered)
                fprintf(['run_collect_by_group: %s %s — skipping %d not-yet-registered ' ...
                         'mouse/mice: %s\n'], group, channel, nnz(~is_registered), ...
                        strjoin(current_mice(~is_registered), ', '));
            end
            current_mice = current_mice(is_registered);
        end

        % go to the next channel if no mouse is ready
        n_mice = numel(current_mice);
        if n_mice == 0
            warning(['run_collect_by_group: no registered %s volumes for group %s%s, ' ...
                     'nothing to collect.'], channel, group, subset_tag);
            continue
        end
        fprintf('run_collect_by_group: %s%s %s — collecting %d mouse/mice: %s\n', ...
            group, subset_tag, channel, n_mice, strjoin(current_mice, ', '));

        % one cell per mouse, and the file each comes from
        vols = cell(1, n_mice);
        source_files = cell(1, n_mice);

        for i = 1:n_mice

            % load the mouse's registered volume of the channel
            mouse_name = current_mice{i};
            source_files{i} = fullfile(base_dir, mouse_name, 'lightsuite', ...
                'volume_registered', channel_file);
            vols{i} = single(loadVolume(source_files(i), 1));

        end

        % stack the mice along the 4th dimension; a cohort's volumes share one grid,
        % set by get_atlas_crop (CCF 900 x 800 x 1140, DeMBA P20 994 x 800 x 1140)
        stack_name = [channel '_4d'];
        S_stack = struct();
        S_stack.(stack_name) = cat(4, vols{:});

        % free the per-mouse copies
        clear vols

        % save the stack (-v7.3 for large arrays), as nano_4d_P20.mat beside
        % nano_4d.mat for an age-filtered run, with its mice and their files
        S_stack.collected_mice = current_mice;
        S_stack.source_files = source_files;
        save(fullfile(base_dir, [stack_name subset_tag '.mat']), '-struct', 'S_stack', ...
            '-v7.3');
        clear S_stack

        % save the nano stack's mice, in order: without it a stack of a partly
        % registered cohort cannot be labelled later (get_cohort_spec reads it)
        if strcmp(channel, 'nano')
            collected_mice = current_mice;
            save(fullfile(base_dir, ['collected_mice' subset_tag '.mat']), ...
                'collected_mice');
        end
    end

end

end
