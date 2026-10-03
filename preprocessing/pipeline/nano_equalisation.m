function nano_equalisation(run_settings)
%NANO_EQUALISATION  Equalise the nano intensity across the slices of the selected mice.
%   NANO_EQUALISATION(run_settings) does the work of run_nano_equalisation,
%   which sets the fields of run_settings (paths, groups_to_process,
%   mice_to_process, atlas_key, save_results, base_output_dir) and says
%   what each one does. The selected mice are processed together.

% settings of run_nano_equalisation, under the names the code below uses
paths = run_settings.paths;
groups_to_process = run_settings.groups_to_process;
mice_to_process = run_settings.mice_to_process;
atlas_key = run_settings.atlas_key;
save_results = run_settings.save_results;
base_output_dir = run_settings.base_output_dir;

%% Add paths

% the toolboxes are on the path from sep_setup_paths; only the atlas folder is
% added here, although nothing below reads the atlas
atlas = get_atlas(atlas_key);
atlas_dir = atlas.dir;
addpath(atlas_dir);

%% Load the volumes

if ~exist(base_output_dir, 'dir')
    mkdir(base_output_dir);
end

% check that the registry still lists the adults in their legacy order, then take
% the mice named, or else every mouse of the groups
get_cohort('verify');
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end

num_mice = numel(cohort);
processed_mouse_names = {cohort.name};
processed_mouse_groups = {cohort.group};
fprintf('run_nano_equalisation: %d mouse/mice selected.\n', num_mice);

% the size of each mouse's volume, and the largest of each dimension
fprintf('--- Phase 1: Scanning Dimensions ---\n');

[dim_store, file_paths] = scan_dimensions(cohort, paths, num_mice);

MAX_H = max(dim_store(:, 1));
MAX_W = max(dim_store(:, 2));
MAX_Z = max(dim_store(:, 3));

% every volume in one array, padded to the largest size
fprintf('--- Phase 2: Loading Volumes into Unified 4D Matrix ---\n');
fprintf('  Max Dimensions: [%d x %d x %d]\n', MAX_H, MAX_W, MAX_Z);
fprintf('  Allocating memory...\n');

nano_4d = load_volumes(file_paths, dim_store, num_mice, MAX_H, MAX_W, MAX_Z);

%% Slice statistics

fprintf('--- Phase 3: Calculating Statistics ---\n');

[intensity_medians, intensity_iqrs] = slice_statistics_raw(nano_4d, dim_store, ...
    num_mice, MAX_Z, processed_mouse_names);

%% Figures before equalisation

fprintf('--- Phase 4: Generating Diagnostic Plots ---\n');

% each slice's median relative to its mouse's mean median
mouse_consensus = repmat(nanmean(intensity_medians, 1), [size(intensity_medians, 1), 1]);
rel_diff_map = (intensity_medians - mouse_consensus) ./ mouse_consensus;

[fig_traces, fig_heatmap_abs, fig_heatmap_rel] = plot_intensity_raw(intensity_medians, ...
    rel_diff_map, processed_mouse_names, num_mice, MAX_Z);

% save them with the statistics, named with the time of the run
if save_results
    timestamp = save_statistics_raw(base_output_dir, intensity_medians, ...
        intensity_iqrs, rel_diff_map, processed_mouse_names, processed_mouse_groups, ...
        fig_traces, fig_heatmap_abs, fig_heatmap_rel);
end

%% Videos before equalisation

fprintf('--- Phase 5: Generating Individual Videos ---\n');

write_videos_raw(nano_4d, dim_store, intensity_medians, processed_mouse_names, ...
    num_mice, base_output_dir, timestamp);

%% Equalise the slices

fprintf('--- Phase 6: Performing Slice Equalization (Window: 5 slices) ---\n');

nano_4d = equalise_slices(nano_4d, intensity_medians, dim_store, num_mice, ...
    processed_mouse_names);

%% Slice statistics after equalisation

fprintf('--- Phase 7: Recalculating Statistics (Equalized) ---\n');

[intensity_medians_eq, intensity_iqrs_eq] = slice_statistics_equalised(nano_4d, ...
    dim_store, num_mice, MAX_Z, processed_mouse_names);

%% Figures after equalisation

fprintf('--- Phase 8: Generating Diagnostic Plots (Equalized) ---\n');

% each slice's median relative to its mouse's mean median
mouse_consensus_eq = repmat(nanmean(intensity_medians_eq, 1), ...
    [size(intensity_medians_eq, 1), 1]); %#ok<*NANMEAN>
rel_diff_map_eq = (intensity_medians_eq - mouse_consensus_eq) ./ mouse_consensus_eq;

[fig_heatmap_abs_eq, fig_heatmap_rel_eq, fig_traces_eq] = plot_intensity_equalised( ...
    intensity_medians_eq, rel_diff_map_eq, processed_mouse_names, num_mice, MAX_Z);

% save them with the statistics, under the time stamp of the first save
if save_results
    save_statistics_equalised(base_output_dir, timestamp, intensity_medians_eq, ...
        intensity_iqrs_eq, rel_diff_map_eq, processed_mouse_names, ...
        processed_mouse_groups, fig_heatmap_abs_eq, fig_heatmap_rel_eq, fig_traces_eq);
end

%% Videos after equalisation

fprintf('--- Phase 9: Generating Individual Videos (Equalized) ---\n');

write_videos_equalised(nano_4d, dim_store, intensity_medians_eq, ...
    processed_mouse_names, num_mice, base_output_dir, timestamp);

%% Save the equalised volumes

fprintf('--- Phase 10: Saving Equalized Volumes to Individual Folders ---\n');

save_equalised_volumes(nano_4d, cohort, paths, dim_store, intensity_medians, ...
    intensity_iqrs, intensity_medians_eq, intensity_iqrs_eq, num_mice);

end

% ===== Local functions: loading =====

function [dim_store, file_paths] = scan_dimensions(cohort, paths, num_mice)
% Height, width and number of slices of each mouse's centred nano volume, one row
% per mouse, and the paths of the volumes; stops when a volume is missing.

% one row per mouse: [H, W, Z]
dim_store = zeros(num_mice, 3);
file_paths = cell(num_mice, 1);

for i = 1:num_mice
    mouse_name = cohort(i).name;
    mouse_type = cohort(i).group;

    nanoPath = fullfile(paths.data, mouse_type, mouse_name, ...
        'lightsuite', 'volume_centered', 'chan02_Cy5.tiff');
    file_paths{i} = nanoPath;

    if ~exist(nanoPath, 'file')
        error('File not found: %s', nanoPath);
    end

    info = imfinfo(nanoPath);
    dim_store(i, 1) = info(1).Height;
    dim_store(i, 2) = info(1).Width;
    dim_store(i, 3) = numel(info);

    fprintf('  Mouse %s: [%d x %d x %d]\n', mouse_name, dim_store(i, 1), ...
        dim_store(i, 2), dim_store(i, 3));
end

end

function nano_4d = load_volumes(file_paths, dim_store, num_mice, MAX_H, MAX_W, MAX_Z)
% All the mice's nano volumes in one 4D array padded with NaN, each volume in the
% top-left corner of its MAX_H x MAX_W x MAX_Z block.

nano_4d = NaN(MAX_H, MAX_W, MAX_Z, num_mice);

for i = 1:num_mice
    fprintf('  Loading [%d/%d] into matrix...\n', i, num_mice);

    % the mouse's own size
    cur_h = dim_store(i, 1);
    cur_w = dim_store(i, 2);
    cur_z = dim_store(i, 3);

    % read its slices into the top-left corner
    for z = 1:cur_z
        nano_4d(1:cur_h, 1:cur_w, z, i) = imread(file_paths{i}, z);
    end
end
fprintf('All volumes loaded.\n');

end

% ===== Local functions: before equalisation =====

function [pmax_val, pmin_val] = background_window(z)
% The percentile window of the background search in slice z: the 15th to the
% 75th percentile in the first nine slices, the 15th to the 50th after.

if z < 10
    pmax_val = 75;
    pmin_val = 15;
else
    pmax_val = 50;
    pmin_val = 15;
end

end

function [intensity_medians, intensity_iqrs] = slice_statistics_raw(nano_4d, ...
    dim_store, num_mice, MAX_Z, processed_mouse_names)
% Median and inter-quartile range of each slice's tissue pixels, before
% equalisation; the padding is set to the slice's mode before the background is
% selected.

intensity_medians = nan(MAX_Z, num_mice);
intensity_iqrs = nan(MAX_Z, num_mice);

for i = 1:num_mice
    fprintf('  Calculating Stats for Mouse %d/%d (%s)...\n', i, num_mice, ...
        processed_mouse_names{i});

    mouse_vol = nano_4d(:, :, :, i);

    slice_medians = nan(MAX_Z, 1);
    slice_iqrs = nan(MAX_Z, 1);

    actual_z = dim_store(i, 3);

    parfor z = 1:MAX_Z

        % padded and empty slices stay NaN
        if z > actual_z
            continue;
        end

        img = im2single(mouse_vol(:, :, z));
        if max(img(:)) == 0
            continue;
        end

        % background between the 15th and the 75th percentile in the first nine
        % slices, the 50th after, with the padding set to the slice's mode
        [pmax_val, pmin_val] = background_window(z);
        img(isnan(img)) = mode(img(:));
        bg_mask = select_background_pixels(img, pmin_val, pmax_val);

        % statistics of the tissue pixels
        fg_pixels = img(~bg_mask);

        if ~isempty(fg_pixels)
            slice_medians(z) = nanmedian(fg_pixels);
            slice_iqrs(z) = quantile(fg_pixels, 0.75)-quantile(fg_pixels, 0.25);
        end
    end

    intensity_medians(:, i) = slice_medians;
    intensity_iqrs(:, i) = slice_iqrs;
end

end

function [fig_traces, fig_heatmap_abs, fig_heatmap_rel] = plot_intensity_raw( ...
    intensity_medians, rel_diff_map, processed_mouse_names, num_mice, MAX_Z)
% Heatmaps and profiles of the slice medians, before equalisation; returns the
% three figures.

% absolute heatmap
fig_heatmap_abs = figure('Name', 'Intensity Heatmap (Absolute)', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.1 0.1 0.5 0.8]);
imagesc(intensity_medians);
colormap(fig_heatmap_abs, sep_palette('intensity'));
c = colorbar;
c.Label.String = 'Median Intensity (Raw)';
xlabel('Mouse');
ylabel('Slice Number');
title('Absolute Intensity Heatmap');
xticks(1:num_mice);
xticklabels(strrep(processed_mouse_names, '_', ' '));
xtickangle(45);
clim([0, max(intensity_medians(:))]);

% relative heatmap
fig_heatmap_rel = figure('Name', 'Intensity Heatmap (Relative)', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.6 0.1 0.5 0.8]);
imagesc(rel_diff_map);
colormap(fig_heatmap_rel, sep_palette('intensity'));
c = colorbar;
c.Label.String = 'Relative Deviation (from Mouse Mean)';
xlabel('Mouse');
ylabel('Slice Number');
title('Relative Deviation Heatmap');
xticks(1:num_mice);
xticklabels(strrep(processed_mouse_names, '_', ' '));
xtickangle(45);
clim([-1, 1]);

% profiles, absolute above and relative below, one line per mouse
fig_traces = figure('Name', 'Intensity Profiles', 'Color', 'w', 'Units', 'normalized', ...
    'Position', [0.1 0.1 0.8 0.8]);
colors = linspace(0.25, 0.75, num_mice)' * [1, 0, 1];

% absolute profiles, with the median over mice
subplot(2, 1, 1);
hold on;
for i = 1:num_mice
    plot(intensity_medians(:, i), 'Color', [colors(i, :) 0.6], 'LineWidth', 2, ...
        'DisplayName', processed_mouse_names{i});
end
group_avg = nanmedian(intensity_medians, 2);
plot(group_avg, 'k-', 'LineWidth', 2, 'DisplayName', 'Group Median');
ylabel('Median Tissue Intensity');
title('Absolute Slice Intensity Profiles');
grid on;
xlim([1 MAX_Z]);

% relative profiles, with the median over mice
subplot(2, 1, 2);
hold on;
for i = 1:num_mice
    plot(rel_diff_map(:, i), 'Color', [colors(i, :) 0.6], 'LineWidth', 2, ...
        'DisplayName', processed_mouse_names{i});
end
group_rel_avg = nanmedian(rel_diff_map, 2); %#ok<*NANMEDIAN>
plot(group_rel_avg, 'k-', 'LineWidth', 2, 'DisplayName', 'Group Median');
yline(0, 'k--', 'LineWidth', 2, 'DisplayName', 'Mouse Mean (Zero Dev)');
xlabel('Slice Number');
ylabel('Relative Deviation');
title('Relative Deviation Profiles');
grid on;
xlim([1 MAX_Z]);
ylim([-0.75, 0.75]);

end

function timestamp = save_statistics_raw(base_output_dir, intensity_medians, ...
    intensity_iqrs, rel_diff_map, processed_mouse_names, processed_mouse_groups, ...
    fig_traces, fig_heatmap_abs, fig_heatmap_rel)
% Save the statistics and the three figures from before equalisation; returns the
% time stamp that names them, which every later file of the run takes.

timestamp = datestr(now, 'yyyymmdd_HHMM'); %#ok<DATST,TNOW1>
savePathData = fullfile(base_output_dir, ['Intensity_Stats_' timestamp '.mat']);

% intensity_iqrs is an input here, so the first branch is the one that runs
if exist('intensity_iqrs', 'var')
    save(savePathData, 'intensity_medians', 'intensity_iqrs', 'rel_diff_map', ...
        'processed_mouse_names', 'processed_mouse_groups', '-v7.3');
else
    save(savePathData, 'intensity_medians', 'rel_diff_map', 'processed_mouse_names', ...
        'processed_mouse_groups', '-v7.3');
end
fprintf('Data saved to: %s\n', savePathData);

exportgraphics(fig_traces, fullfile(base_output_dir, ['Plot_Traces_' timestamp '.png']), ...
    'Resolution', 300);
exportgraphics(fig_heatmap_abs, fullfile(base_output_dir, ...
    ['Plot_Heatmap_Abs_' timestamp '.png']), 'Resolution', 300);
exportgraphics(fig_heatmap_rel, fullfile(base_output_dir, ...
    ['Plot_Heatmap_Rel_' timestamp '.png']), 'Resolution', 300);

close all

end

function write_videos_raw(nano_4d, dim_store, intensity_medians, ...
    processed_mouse_names, num_mice, base_output_dir, timestamp)
% One video per mouse, the background masked, before equalisation; writes
% Video_<mouse>_<timestamp>.mp4 in base_output_dir.

% one hidden figure, reused for every frame
h_fig = figure('visible', 'off', 'units', 'pixels', 'position', [100 100 800 600], ...
    'Color', 'k');
set(h_fig, 'InvertHardcopy', 'off');

for i = 1:num_mice
    mouse_name = processed_mouse_names{i};
    fprintf('  Processing Video for Mouse %d/%d: %s...\n', i, num_mice, mouse_name);

    video_filename = fullfile(base_output_dir, ...
        ['Video_' mouse_name '_' timestamp '.mp4']);
    vidObj = VideoWriter(video_filename, 'MPEG-4');
    vidObj.FrameRate = 5;
    vidObj.Quality = 95;
    open(vidObj);

    actual_z = dim_store(i, 3);

    for z = 1:actual_z

        % the slice, cropped to the mouse's own size
        img_raw = nano_4d(:, :, z, i);
        cur_h = dim_store(i, 1);
        cur_w = dim_store(i, 2);
        img_crop = img_raw(1:cur_h, 1:cur_w);

        % an empty slice gets a black frame, so frame n stays slice n
        img_single = single(img_crop);
        if max(img_single(:)) == 0
            clf(h_fig);
            set(gca, 'Color', 'k');
            axis off;
            frame = getframe(h_fig);
            writeVideo(vidObj, frame);
            continue;
        end

        % the frame: the slice, its background masked
        draw_raw_frame(z, img_single, h_fig, img_crop, mouse_name, intensity_medians, ...
            i);

        % add the frame
        frame = getframe(h_fig);
        writeVideo(vidObj, frame);

        if mod(z, 100) == 0
            fprintf('    Frame %d / %d\n', z, actual_z);
        end
    end

    close(vidObj);
    fprintf('    Video saved: %s\n', video_filename);
end

close(h_fig);
fprintf('Script finished.\n');

end

function draw_raw_frame(z, img_single, h_fig, img_crop, mouse_name, intensity_medians, ...
    i)
% One frame: the slice in grey, its background transparent on black, and its median.

% background mask, as for the statistics
[pmax_val, pmin_val] = background_window(z);
img_single(isnan(img_single)) = mode(img_single(:));
bg_mask = select_background_pixels(img_single, pmin_val, pmax_val);

% draw the slice in grey, the background transparent on black
clf(h_fig);
h_im = imagesc(img_crop);
colormap(sep_palette('anatomy'));
clim([0 5000]);
set(h_im, 'AlphaData', ~bg_mask);
axis image;
axis off;
set(gca, 'Color', 'k');

title([sprintf('%s - Slice %d', strrep(mouse_name, '_', ' '), z), ...
    ' - median = ', num2str(round(intensity_medians(z, i), 2))], ...
    'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');

end

% ===== Local functions: equalisation =====

function nano_4d = equalise_slices(nano_4d, intensity_medians, dim_store, num_mice, ...
    processed_mouse_names)
% Scale each slice so that its median becomes the moving median over 5 slices; in
% place (the same name in and out), padded or empty slices keep a factor of 1.

% 5 slices: two before, the slice, two after
window_size = 5;

for i = 1:num_mice
    mouse_name = processed_mouse_names{i};
    fprintf('  Equalizing Mouse %d/%d: %s...\n', i, num_mice, mouse_name);

    % the slice medians before equalisation, and their moving median
    raw_medians = intensity_medians(:, i);
    target_medians = movmedian(raw_medians, window_size, 'omitnan');

    % scaling factors; 1 where the slice was empty or padded (division by zero, NaN)
    scaling_factors = target_medians ./ raw_medians;
    scaling_factors(isnan(scaling_factors) | isinf(scaling_factors)) = 1;

    % scale the slices in memory
    actual_z = dim_store(i, 3);

    for z = 1:actual_z
        current_factor = scaling_factors(z);

        % a factor of 1 changes nothing
        if current_factor == 1
            continue;
        end

        nano_4d(:, :, z, i) = nano_4d(:, :, z, i) * current_factor;
    end
end
fprintf('Equalization applied to nano_4d in memory.\n');

end

function [intensity_medians_eq, intensity_iqrs_eq] = slice_statistics_equalised( ...
    nano_4d, dim_store, num_mice, MAX_Z, processed_mouse_names)
% Median and inter-quartile range of each slice's tissue pixels, after
% equalisation; the mask comes from a copy with the padding set to its mode, the
% statistics from the slice itself.

intensity_medians_eq = nan(MAX_Z, num_mice);
intensity_iqrs_eq = nan(MAX_Z, num_mice);

for i = 1:num_mice
    fprintf('  Stats (Eq) for Mouse %d/%d (%s)...\n', i, num_mice, ...
        processed_mouse_names{i});

    mouse_vol = nano_4d(:, :, :, i);
    actual_z = dim_store(i, 3);

    slice_medians = nan(MAX_Z, 1);
    slice_iqrs = nan(MAX_Z, 1);

    parfor z = 1:MAX_Z

        % padded and empty slices stay NaN
        if z > actual_z
            continue;
        end

        img = im2single(mouse_vol(:, :, z));
        if max(img(:)) == 0
            continue;
        end

        % the background window, as before equalisation
        [pmax_val, pmin_val] = background_window(z);

        % background mask, from a copy with the padding set to the slice's mode
        img_temp = img;
        img_temp(isnan(img_temp)) = mode(img_temp(:));
        bg_mask = select_background_pixels(img_temp, pmin_val, pmax_val);

        % statistics of the tissue pixels
        fg_pixels = img(~bg_mask);

        if ~isempty(fg_pixels)
            slice_medians(z) = nanmedian(fg_pixels);
            slice_iqrs(z) = quantile(fg_pixels, 0.75)-quantile(fg_pixels, 0.25);
        end
    end

    intensity_medians_eq(:, i) = slice_medians;
    intensity_iqrs_eq(:, i) = slice_iqrs;
end

end

function [fig_heatmap_abs_eq, fig_heatmap_rel_eq, fig_traces_eq] = ...
    plot_intensity_equalised(intensity_medians_eq, rel_diff_map_eq, ...
    processed_mouse_names, num_mice, MAX_Z)
% Heatmaps and profiles of the slice medians, after equalisation; returns the
% three figures.

% absolute heatmap
fig_heatmap_abs_eq = figure('Name', 'Intensity Heatmap (Absolute - Equalized)', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.5 0.8]);
imagesc(intensity_medians_eq);
colormap(fig_heatmap_abs_eq, sep_palette('intensity'));
c = colorbar;
c.Label.String = 'Median Intensity (Equalized)';
xlabel('Mouse');
ylabel('Slice Number');
title('Absolute Intensity Heatmap (Equalized)');
xticks(1:num_mice);
xticklabels(strrep(processed_mouse_names, '_', ' '));
xtickangle(45);
clim([0, max(intensity_medians_eq(:))]);

% relative heatmap
fig_heatmap_rel_eq = figure('Name', 'Intensity Heatmap (Relative - Equalized)', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.6 0.1 0.5 0.8]);
imagesc(rel_diff_map_eq);
colormap(fig_heatmap_rel_eq, sep_palette('intensity'));
c = colorbar;
c.Label.String = 'Relative Deviation';
xlabel('Mouse');
ylabel('Slice Number');
title('Relative Deviation Heatmap (Equalized)');
xticks(1:num_mice);
xticklabels(strrep(processed_mouse_names, '_', ' '));
xtickangle(45);
clim([-1, 1]);

% profiles, absolute above and relative below, one line per mouse
fig_traces_eq = figure('Name', 'Intensity Profiles (Equalized)', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
colors = linspace(0.25, 0.75, num_mice)' * [1, 0, 1];

% absolute profiles, with the median over mice
subplot(2, 1, 1);
hold on;
for i = 1:num_mice
    plot(intensity_medians_eq(:, i), 'Color', [colors(i, :) 0.6], 'LineWidth', 2, ...
        'DisplayName', processed_mouse_names{i});
end
group_avg_eq = nanmedian(intensity_medians_eq, 2);
plot(group_avg_eq, 'k-', 'LineWidth', 2, 'DisplayName', 'Group Median');
ylabel('Median Tissue Intensity');
title('Absolute Profiles (Equalized)');
grid on;
xlim([1 MAX_Z]);

% relative profiles
subplot(2, 1, 2);
hold on;
for i = 1:num_mice
    plot(rel_diff_map_eq(:, i), 'Color', [colors(i, :) 0.6], 'LineWidth', 2, ...
        'DisplayName', processed_mouse_names{i});
end
yline(0, 'k--', 'LineWidth', 2);
xlabel('Slice Number');
ylabel('Relative Deviation');
title('Relative Deviation Profiles (Equalized)');
grid on;
xlim([1 MAX_Z]);
ylim([-0.75, 0.75]);

end

function save_statistics_equalised(base_output_dir, timestamp, intensity_medians_eq, ...
    intensity_iqrs_eq, rel_diff_map_eq, processed_mouse_names, ...
    processed_mouse_groups, fig_heatmap_abs_eq, fig_heatmap_rel_eq, fig_traces_eq)
% Save the statistics and the three figures from after equalisation, named with
% the time stamp of the statistics saved before it.

savePathData = fullfile(base_output_dir, ['Intensity_Stats_Equalized_' timestamp '.mat']);
save(savePathData, 'intensity_medians_eq', 'intensity_iqrs_eq', 'rel_diff_map_eq', ...
    'processed_mouse_names', 'processed_mouse_groups', '-v7.3');
fprintf('Equalized Data saved to: %s\n', savePathData);

exportgraphics(fig_traces_eq, fullfile(base_output_dir, ...
    ['Plot_Traces_Equalized_' timestamp '.png']), 'Resolution', 300);
exportgraphics(fig_heatmap_abs_eq, fullfile(base_output_dir, ...
    ['Plot_Heatmap_Abs_Equalized_' timestamp '.png']), 'Resolution', 300);
exportgraphics(fig_heatmap_rel_eq, fullfile(base_output_dir, ...
    ['Plot_Heatmap_Rel_Equalized_' timestamp '.png']), 'Resolution', 300);

close all

end

function write_videos_equalised(nano_4d, dim_store, intensity_medians_eq, ...
    processed_mouse_names, num_mice, base_output_dir, timestamp)
% One video per mouse, the background masked, after equalisation; writes
% Video_<mouse>_Equalized_<timestamp>.mp4 in base_output_dir.

% one hidden figure, reused for every frame
h_fig = figure('visible', 'off', 'units', 'pixels', 'position', [100 100 800 600], ...
    'Color', 'k');
set(h_fig, 'InvertHardcopy', 'off');

for i = 1:num_mice
    mouse_name = processed_mouse_names{i};
    fprintf('  Processing Video for Mouse %d/%d: %s...\n', i, num_mice, mouse_name);

    video_filename = fullfile(base_output_dir, ['Video_' mouse_name '_Equalized_' ...
        timestamp '.mp4']);
    vidObj = VideoWriter(video_filename, 'MPEG-4');
    vidObj.FrameRate = 5;
    vidObj.Quality = 95;
    open(vidObj);

    actual_z = dim_store(i, 3);

    for z = 1:actual_z

        % the equalised slice, cropped to the mouse's own size
        img_eq = nano_4d(:, :, z, i);
        cur_h = dim_store(i, 1);
        cur_w = dim_store(i, 2);
        img_crop = img_eq(1:cur_h, 1:cur_w);

        % an empty slice gets a black frame, so frame n stays slice n
        img_single = single(img_crop);
        if max(img_single(:)) == 0
            clf(h_fig);
            set(gca, 'Color', 'k');
            axis off;
            frame = getframe(h_fig);
            writeVideo(vidObj, frame);
            continue;
        end

        % the frame: the slice, its background masked
        draw_equalised_frame(z, img_single, h_fig, img_crop, mouse_name, ...
            intensity_medians_eq, i);

        % add the frame
        frame = getframe(h_fig);
        writeVideo(vidObj, frame);

        if mod(z, 100) == 0
            fprintf('    Frame %d / %d\n', z, actual_z);
        end
    end
    close(vidObj);
end

close(h_fig);
fprintf('Full pipeline finished.\n');

end

function draw_equalised_frame(z, img_single, h_fig, img_crop, mouse_name, ...
    intensity_medians_eq, i)
% One frame: the equalised slice in grey, its background transparent on black.

% background mask, on the equalised slice
[pmax_val, pmin_val] = background_window(z);
img_single(isnan(img_single)) = mode(img_single(:));
bg_mask = select_background_pixels(img_single, pmin_val, pmax_val);

% draw the slice in grey, the background transparent on black
clf(h_fig);
h_im = imagesc(img_crop);
colormap(sep_palette('anatomy'));
clim([0 5000]);
set(h_im, 'AlphaData', ~bg_mask);
axis image;
axis off;
set(gca, 'Color', 'k');

title([sprintf('%s (Eq) - Slice %d', strrep(mouse_name, '_', ' '), z), ...
    ' - med = ', num2str(round(intensity_medians_eq(z, i), 2))], ...
    'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');

end

function save_equalised_volumes(nano_4d, cohort, paths, dim_store, intensity_medians, ...
    intensity_iqrs, intensity_medians_eq, intensity_iqrs_eq, num_mice)
% Save each mouse's equalised volume and slice statistics, cropped to its own
% size, in lightsuite\correction_output\equalized_volume.mat.

for i = 1:num_mice

    % the mouse and its correction folder
    current_mouse = cohort(i).name;
    current_type = cohort(i).group;
    base_dir = fullfile(paths.data, current_type);
    correction_dir = fullfile(base_dir, current_mouse, 'lightsuite', ...
        'correction_output');

    if ~exist(correction_dir, 'dir')
        mkdir(correction_dir);
        fprintf('  Created directory: %s\n', correction_dir);
    end

    fprintf('  Saving data for %s (%d/%d)...\n', current_mouse, i, num_mice);

    % the volume, cropped to the mouse's own size
    cur_h = dim_store(i, 1);
    cur_w = dim_store(i, 2);
    cur_z = dim_store(i, 3);
    equalized_volume = nano_4d(1:cur_h, 1:cur_w, 1:cur_z, i);

    % its slice statistics, before and after, without the padded slices
    stats_intensity_median_raw = intensity_medians(1:cur_z, i);
    stats_intensity_iqr_raw = intensity_iqrs(1:cur_z, i);
    stats_intensity_median_eq = intensity_medians_eq(1:cur_z, i);
    stats_intensity_iqr_eq = intensity_iqrs_eq(1:cur_z, i);

    % save
    save_filename = fullfile(correction_dir, 'equalized_volume.mat');

    save(save_filename, ...
        'equalized_volume', ...
        'stats_intensity_median_raw', ...
        'stats_intensity_iqr_raw', ...
        'stats_intensity_median_eq', ...
        'stats_intensity_iqr_eq', ...
        '-v7.3');
end
fprintf('All volumes saved successfully.\n');

end
