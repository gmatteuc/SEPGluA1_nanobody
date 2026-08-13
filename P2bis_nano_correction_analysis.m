clear all
close all
clc

% /// Pipeline script #2bis (alternate): read centered volumes and perform nano channel median equalization across slices  /// 

%% 1. User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Cohort selection (mice come from the shared registry get_cohort.m).
% Set mice_to_process to {} to process every mouse in groups_to_process.
% NOTE: this used to be a list of numeric INDICES into a hardcoded mouse
% list (mice_to_process = 1:17); it is now a list of mouse NAMES, so the
% selection no longer depends on the order of the registry.
groups_to_process = {'young'};                  % 'rws' | 'naive' | 'behavior' | 'young'
mice_to_process   = {'MG909_SepGluA_P20', 'MG910_SepGluA_P20', 'MG911_SepGluA_P16', ...
                     'MG912_SepGluA_P20', 'MG913_SepGluA_P20', 'MG914_SepGluA_P28'};                         % {} = all mice in groups_to_process

% Reference atlas
atlas_key = 'ccf';

% Output settings
save_results = true;
base_output_dir = fullfile(paths.data, 'intensity_diagnostics');

%% 2. Add paths

% Get Allen data path
atlas = get_atlas(atlas_key);
allenDir = atlas.dir;
% Add Allen data path
addpath(allenDir);

%% 3. Scan Dimensions and Load to RAM

if ~exist(base_output_dir, 'dir'), mkdir(base_output_dir); end

% Resolve cohort
get_cohort('verify');
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end

num_mice = numel(cohort);
processed_mouse_names  = {cohort.name};
processed_mouse_groups = {cohort.group};
fprintf('P2bis: %d mouse/mice selected.\n', num_mice);

fprintf('--- Phase 1: Scanning Dimensions ---\n');

dim_store = zeros(num_mice, 3); % [H, W, Z]
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

    fprintf('  Mouse %s: [%d x %d x %d]\n', mouse_name, dim_store(i,1), dim_store(i,2), dim_store(i,3));
end

% Determine Global Max Dimensions
MAX_H = max(dim_store(:,1));
MAX_W = max(dim_store(:,2));
MAX_Z = max(dim_store(:,3));

fprintf('--- Phase 2: Loading Volumes into Unified 4D Matrix ---\n');
fprintf('  Max Dimensions: [%d x %d x %d]\n', MAX_H, MAX_W, MAX_Z);
fprintf('  Allocating memory...\n');

% Allocate 4D Matrix
nano_4d = NaN(MAX_H, MAX_W, MAX_Z, num_mice);

for i = 1:num_mice
    fprintf('  Loading [%d/%d] into matrix...\n', i, num_mice);

    % Get individual dimensions
    cur_h = dim_store(i, 1);
    cur_w = dim_store(i, 2);
    cur_z = dim_store(i, 3);

    % Load slices and place them into the top-left corner
    for z = 1:cur_z
        nano_4d(1:cur_h, 1:cur_w, z, i) = imread(file_paths{i}, z);
    end
end
fprintf('All volumes loaded.\n');

%% 4. Calculate Slice Statistics from Memory

fprintf('--- Phase 3: Calculating Statistics ---\n');

intensity_medians = nan(MAX_Z, num_mice);
intensity_iqrs  = nan(MAX_Z, num_mice);

for i = 1:num_mice
    fprintf('  Calculating Stats for Mouse %d/%d (%s)...\n', i, num_mice, processed_mouse_names{i});

    mouse_vol = nano_4d(:,:,:,i);

    slice_medians = nan(MAX_Z, 1);
    slice_iqrs  = nan(MAX_Z, 1);

    actual_z = dim_store(i, 3);

    parfor z = 1:MAX_Z
        if z > actual_z, continue; end

        img = im2single(mouse_vol(:,:,z));
        if max(img(:)) == 0, continue; end

        if z < 10, pmax_val = 75; pmin_val = 15;
        else, pmax_val = 50; pmin_val = 15;
        end
        img(isnan(img))=mode(img(:));
        bg_mask = select_background_pixels(img, pmin_val, pmax_val);

        fg_pixels = img(~bg_mask);

        if ~isempty(fg_pixels)
            slice_medians(z) = nanmedian(fg_pixels);
            slice_iqrs(z)  = quantile(fg_pixels,0.25)-quantile(fg_pixels,0.75);
        end
    end

    intensity_medians(:, i) = slice_medians;
    intensity_iqrs(:, i)  = slice_iqrs;
end

%% 5. Visualization & Analysis

fprintf('--- Phase 4: Generating Diagnostic Plots ---\n');

% --- Calculation for Relative Metrics ---
mouse_consensus = repmat(nanmean(intensity_medians,1), [size(intensity_medians,1), 1]);
rel_diff_map = (intensity_medians - mouse_consensus) ./ mouse_consensus;

% --- Plot 1: Absolute Heatmap ---
f2 = figure('Name', 'Intensity Heatmap (Absolute)', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.5 0.8]);
imagesc(intensity_medians);
colormap(f2, hot);
c = colorbar; c.Label.String = 'Median Intensity (Raw)';
xlabel('Mouse'); ylabel('Slice Number'); title('Absolute Intensity Heatmap');
xticks(1:num_mice); xticklabels(strrep(processed_mouse_names, '_', ' ')); xtickangle(45);
clim([0, max(intensity_medians(:))]);

% --- Plot 2: Relative Deviation Heatmap ---
f3 = figure('Name', 'Intensity Heatmap (Relative)', 'Color', 'w', 'Units', 'normalized', 'Position', [0.6 0.1 0.5 0.8]);
imagesc(rel_diff_map);
colormap(f3, hot);
c = colorbar; c.Label.String = 'Relative Deviation (from Mouse Mean)';
xlabel('Mouse'); ylabel('Slice Number'); title('Relative Deviation Heatmap');
xticks(1:num_mice); xticklabels(strrep(processed_mouse_names, '_', ' ')); xtickangle(45);
clim([-1, 1]);

% --- Plot 3: Profiles (Absolute & Relative Subplots) ---
f1 = figure('Name', 'Intensity Profiles', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
colors = linspace(0.25, 0.75, num_mice)' * [1, 0, 1];

% Subplot 1: Absolute Profiles
subplot(2,1,1);
hold on;
for i = 1:num_mice
    plot(intensity_medians(:, i), 'Color', [colors(i,:) 0.6], 'LineWidth', 2, 'DisplayName', processed_mouse_names{i});
end
% Group Median line
group_avg = nanmedian(intensity_medians, 2);
plot(group_avg, 'k-', 'LineWidth', 2, 'DisplayName', 'Group Median');
ylabel('Median Tissue Intensity'); title('Absolute Slice Intensity Profiles');
grid on; xlim([1 MAX_Z]);

% Subplot 2: Relative Profiles (Matching f3)
subplot(2,1,2);
hold on;
for i = 1:num_mice
    plot(rel_diff_map(:, i), 'Color', [colors(i,:) 0.6], 'LineWidth', 2, 'DisplayName', processed_mouse_names{i});
end
% Group Median line
group_rel_avg = nanmedian(rel_diff_map, 2); %#ok<*NANMEDIAN>
plot(group_rel_avg, 'k-', 'LineWidth', 2, 'DisplayName', 'Group Median');
yline(0, 'k--', 'LineWidth', 2, 'DisplayName', 'Mouse Mean (Zero Dev)');
xlabel('Slice Number'); ylabel('Relative Deviation'); title('Relative Deviation Profiles');
grid on; xlim([1 MAX_Z]); ylim([-0.75, 0.75]); % Matching heatmap limits

% --- Saving ---
if save_results
    timestamp = datestr(now, 'yyyymmdd_HHMM'); %#ok<DATST,TNOW1>
    savePathData = fullfile(base_output_dir, ['Intensity_Stats_' timestamp '.mat']);

    % Robust save (check if variables exist)
    if exist('intensity_iqrs', 'var')
        save(savePathData, 'intensity_medians', 'intensity_iqrs', 'rel_diff_map', 'processed_mouse_names', 'processed_mouse_groups', '-v7.3');
    else
        save(savePathData, 'intensity_medians', 'rel_diff_map', 'processed_mouse_names', 'processed_mouse_groups', '-v7.3');
    end
    fprintf('Data saved to: %s\n', savePathData);

    exportgraphics(f1, fullfile(base_output_dir, ['Plot_Traces_' timestamp '.png']), 'Resolution', 300);
    exportgraphics(f2, fullfile(base_output_dir, ['Plot_Heatmap_Abs_' timestamp '.png']), 'Resolution', 300);
    exportgraphics(f3, fullfile(base_output_dir, ['Plot_Heatmap_Rel_' timestamp '.png']), 'Resolution', 300);

    close all
end

%% 7. Generate Individual Videos (Masked + Gray)

fprintf('--- Phase 5: Generating Individual Videos ---\n');

% Setup invisible figure once to reuse
h_fig = figure('visible', 'off', 'units', 'pixels', 'position', [100 100 800 600], 'Color', 'k');
set(h_fig, 'InvertHardcopy', 'off');

for i = 1:num_mice
    mouse_name = processed_mouse_names{i};
    fprintf('  Processing Video for Mouse %d/%d: %s...\n', i, num_mice, mouse_name);

    video_filename = fullfile(base_output_dir, ['Video_' mouse_name '_' timestamp '.mp4']);
    vidObj = VideoWriter(video_filename, 'MPEG-4');
    vidObj.FrameRate = 5;
    vidObj.Quality = 95;
    open(vidObj);

    actual_z = dim_store(i, 3);

    for z = 1:actual_z
        % 1. Access Image from RAM
        img_uint16 = nano_4d(:,:,z,i);

        % Cropping to original size
        cur_h = dim_store(i, 1);
        cur_w = dim_store(i, 2);
        img_crop = img_uint16(1:cur_h, 1:cur_w);

        % 2. Compute Mask (Same logic as Stats phase)
        img_single = single(img_crop);
        if max(img_single(:)) == 0
            % Skip empty frames to keep video smooth? Or write black frame?
            % Writing black frame to maintain Z-index alignment
            clf(h_fig); set(gca, 'Color', 'k'); axis off;
            frame = getframe(h_fig); writeVideo(vidObj, frame);
            continue;
        end

        if z < 10, pmax_val = 75; pmin_val = 15; else, pmax_val = 50; pmin_val = 15; end
        img_single(isnan(img_single))=mode(img_single(:));
        bg_mask = select_background_pixels(img_single, pmin_val, pmax_val);

        % 3. Plot with Alpha Mask
        clf(h_fig);

        % Display image
        h_im = imagesc(img_crop);
        colormap(gray);
        clim([0 5000]);

        % Apply Alpha Data
        set(h_im, 'AlphaData', ~bg_mask);

        % Formatting
        axis image; axis off;
        set(gca, 'Color', 'k');

        title([sprintf('%s - Slice %d', strrep(mouse_name, '_', ' '), z),' - median = ',num2str(round(intensity_medians(z,i),2))], ...
            'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');

        % 4. Write Frame
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

%% 8. Perform Slice-wise Equalization

fprintf('--- Phase 6: Performing Slice Equalization (Window: 5 slices) ---\n');

% Set window size (2 before + 2 after)
window_size = 5;

for i = 1:num_mice
    mouse_name = processed_mouse_names{i};
    fprintf('  Equalizing Mouse %d/%d: %s...\n', i, num_mice, mouse_name);

    % Get the raw medians calculated in Phase 3
    raw_medians = intensity_medians(:, i);

    % Calculate the Target Median (Moving median of the neighborhood)
    target_medians = movmedian(raw_medians, window_size, 'omitnan');

    % Calculate Scaling Factors (Target / Raw)
    scaling_factors = target_medians ./ raw_medians;

    % Handle Division by Zero or NaNs (if raw slice was empty/padded)
    scaling_factors(isnan(scaling_factors) | isinf(scaling_factors)) = 1;

    % Apply Scaling to the Volume in RAM
    actual_z = dim_store(i, 3);

    for z = 1:actual_z
        current_factor = scaling_factors(z);

        % Skip if factor is 1 (optimization)
        if current_factor == 1, continue; end

        % Apply multiplicative scaling
        nano_4d(:,:,z,i) = nano_4d(:,:,z,i) * current_factor;
    end
end
fprintf('Equalization applied to nano_4d in memory.\n');

%% 9. Recalculate Statistics on Equalized Data

fprintf('--- Phase 7: Recalculating Statistics (Equalized) ---\n');

% Reset stats matrices
intensity_medians_eq = nan(MAX_Z, num_mice);
intensity_iqrs_eq  = nan(MAX_Z, num_mice);

for i = 1:num_mice
    fprintf('  Stats (Eq) for Mouse %d/%d (%s)...\n', i, num_mice, processed_mouse_names{i});

    mouse_vol = nano_4d(:,:,:,i);
    actual_z = dim_store(i, 3);

    slice_medians = nan(MAX_Z, 1);
    slice_iqrs  = nan(MAX_Z, 1);

    parfor z = 1:MAX_Z
        if z > actual_z, continue; end

        img = im2single(mouse_vol(:,:,z));
        if max(img(:)) == 0, continue; end

        if z < 10, pmax_val = 75; pmin_val = 15; else, pmax_val = 50; pmin_val = 15; end

        % Masking
        img_temp = img;
        img_temp(isnan(img_temp)) = mode(img_temp(:));
        bg_mask = select_background_pixels(img_temp, pmin_val, pmax_val);

        fg_pixels = img(~bg_mask);

        if ~isempty(fg_pixels)
            slice_medians(z) = nanmedian(fg_pixels);
            slice_iqrs(z)  = quantile(fg_pixels,0.25)-quantile(fg_pixels,0.75);
        end
    end

    intensity_medians_eq(:, i) = slice_medians;
    intensity_iqrs_eq(:, i)  = slice_iqrs;
end

%% 10. Visualization (Equalized)

fprintf('--- Phase 8: Generating Diagnostic Plots (Equalized) ---\n');

% --- Calculation for Relative Metrics ---
mouse_consensus_eq = repmat(nanmean(intensity_medians_eq,1), [size(intensity_medians_eq,1), 1]); %#ok<*NANMEAN>
rel_diff_map_eq = (intensity_medians_eq - mouse_consensus_eq) ./ mouse_consensus_eq;

% --- Plot 1: Absolute Heatmap (Eq) ---
f4 = figure('Name', 'Intensity Heatmap (Absolute - Equalized)', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.5 0.8]);
imagesc(intensity_medians_eq);
colormap(f4, hot);
c = colorbar; c.Label.String = 'Median Intensity (Equalized)';
xlabel('Mouse'); ylabel('Slice Number'); title('Absolute Intensity Heatmap (Equalized)');
xticks(1:num_mice); xticklabels(strrep(processed_mouse_names, '_', ' ')); xtickangle(45);
clim([0, max(intensity_medians_eq(:))]);

% --- Plot 2: Relative Deviation Heatmap (Eq) ---
f5 = figure('Name', 'Intensity Heatmap (Relative - Equalized)', 'Color', 'w', 'Units', 'normalized', 'Position', [0.6 0.1 0.5 0.8]);
imagesc(rel_diff_map_eq);
colormap(f5, hot);
c = colorbar; c.Label.String = 'Relative Deviation';
xlabel('Mouse'); ylabel('Slice Number'); title('Relative Deviation Heatmap (Equalized)');
xticks(1:num_mice); xticklabels(strrep(processed_mouse_names, '_', ' ')); xtickangle(45);
clim([-1, 1]);

% --- Plot 3: Profiles (Eq) ---
f6 = figure('Name', 'Intensity Profiles (Equalized)', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
colors = linspace(0.25, 0.75, num_mice)' * [1, 0, 1];

subplot(2,1,1); hold on;
for i = 1:num_mice
    plot(intensity_medians_eq(:, i), 'Color', [colors(i,:) 0.6], 'LineWidth', 2, 'DisplayName', processed_mouse_names{i});
end
group_avg_eq = nanmedian(intensity_medians_eq, 2);
plot(group_avg_eq, 'k-', 'LineWidth', 2, 'DisplayName', 'Group Median');
ylabel('Median Tissue Intensity'); title('Absolute Profiles (Equalized)');
grid on; xlim([1 MAX_Z]);

subplot(2,1,2); hold on;
for i = 1:num_mice
    plot(rel_diff_map_eq(:, i), 'Color', [colors(i,:) 0.6], 'LineWidth', 2, 'DisplayName', processed_mouse_names{i});
end
yline(0, 'k--', 'LineWidth', 2);
xlabel('Slice Number'); ylabel('Relative Deviation'); title('Relative Deviation Profiles (Equalized)');
grid on; xlim([1 MAX_Z]); ylim([-0.75, 0.75]);

% --- Saving ---
if save_results
    savePathData = fullfile(base_output_dir, ['Intensity_Stats_Equalized_' timestamp '.mat']);
    save(savePathData, 'intensity_medians_eq', 'intensity_iqrs_eq', 'rel_diff_map_eq', 'processed_mouse_names', 'processed_mouse_groups', '-v7.3');
    fprintf('Equalized Data saved to: %s\n', savePathData);

    exportgraphics(f6, fullfile(base_output_dir, ['Plot_Traces_Equalized_' timestamp '.png']), 'Resolution', 300);
    exportgraphics(f4, fullfile(base_output_dir, ['Plot_Heatmap_Abs_Equalized_' timestamp '.png']), 'Resolution', 300);
    exportgraphics(f5, fullfile(base_output_dir, ['Plot_Heatmap_Rel_Equalized_' timestamp '.png']), 'Resolution', 300);

    close all
end

%% 11. Generate Individual Videos (Equalized)

fprintf('--- Phase 9: Generating Individual Videos (Equalized) ---\n');

h_fig = figure('visible', 'off', 'units', 'pixels', 'position', [100 100 800 600], 'Color', 'k');
set(h_fig, 'InvertHardcopy', 'off');

for i = 1:num_mice
    mouse_name = processed_mouse_names{i};
    fprintf('  Processing Video for Mouse %d/%d: %s...\n', i, num_mice, mouse_name);

    video_filename = fullfile(base_output_dir, ['Video_' mouse_name '_Equalized_' timestamp '.mp4']);
    vidObj = VideoWriter(video_filename, 'MPEG-4');
    vidObj.FrameRate = 5;
    vidObj.Quality = 95;
    open(vidObj);

    actual_z = dim_store(i, 3);

    for z = 1:actual_z
        % Access Equalized Image from RAM
        img_eq = nano_4d(:,:,z,i);

        cur_h = dim_store(i, 1);
        cur_w = dim_store(i, 2);
        img_crop = img_eq(1:cur_h, 1:cur_w);

        % Compute Mask (on equalized data)
        img_single = single(img_crop);
        if max(img_single(:)) == 0
            clf(h_fig); set(gca, 'Color', 'k'); axis off;
            frame = getframe(h_fig); writeVideo(vidObj, frame);
            continue;
        end

        if z < 10, pmax_val = 75; pmin_val = 15; else, pmax_val = 50; pmin_val = 15; end
        img_single(isnan(img_single))=mode(img_single(:));
        bg_mask = select_background_pixels(img_single, pmin_val, pmax_val);

        clf(h_fig);
        h_im = imagesc(img_crop);
        colormap(gray);
        clim([0 5000]);
        set(h_im, 'AlphaData', ~bg_mask);
        axis image; axis off; set(gca, 'Color', 'k');

        title([sprintf('%s (Eq) - Slice %d', strrep(mouse_name, '_', ' '), z),' - med = ',num2str(round(intensity_medians_eq(z,i),2))], ...
            'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');

        frame = getframe(h_fig);
        writeVideo(vidObj, frame);

        if mod(z, 100) == 0, fprintf('    Frame %d / %d\n', z, actual_z); end
    end
    close(vidObj);
end

close(h_fig);
fprintf('Full pipeline finished.\n');

%% 12. Save Equalized Volumes to Individual Mouse Directories

fprintf('--- Phase 10: Saving Equalized Volumes to Individual Folders ---\n');

for i = 1:num_mice
    % 1. Get Mouse Identity
    current_mouse = cohort(i).name;
    current_type  = cohort(i).group;

    % 2. Define Output Directory
    base_dir = fullfile(paths.data, current_type);
    group_output_dir = fullfile(base_dir, current_mouse, 'lightsuite', 'correction_output');

    if ~exist(group_output_dir, 'dir')
        mkdir(group_output_dir);
        fprintf('  Created directory: %s\n', group_output_dir);
    end

    fprintf('  Saving data for %s (%d/%d)...\n', current_mouse, i, num_mice);

    % 3. Extract and Crop Volume
    cur_h = dim_store(i, 1);
    cur_w = dim_store(i, 2);
    cur_z = dim_store(i, 3);

    % Extract strictly the valid data region
    equalized_volume = nano_4d(1:cur_h, 1:cur_w, 1:cur_z, i);

    % 4. Extract Statistics for this specific mouse (Trimmed to valid Z)
    stats_intensity_median_raw = intensity_medians(1:cur_z, i);
    stats_intensity_iqr_raw    = intensity_iqrs(1:cur_z, i);
    stats_intensity_median_eq  = intensity_medians_eq(1:cur_z, i);
    stats_intensity_iqr_eq     = intensity_iqrs_eq(1:cur_z, i);

    % 5. Save to .mat file
    save_filename = fullfile(group_output_dir, 'equalized_volume.mat');

    save(save_filename, ...
        'equalized_volume', ...
        'stats_intensity_median_raw', ...
        'stats_intensity_iqr_raw', ...
        'stats_intensity_median_eq', ...
        'stats_intensity_iqr_eq', ...
        '-v7.3');
end
fprintf('All volumes saved successfully.\n');