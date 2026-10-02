function normalise_groups(run_settings)
%NORMALISE_GROUPS  Normalise the mice of each group onto a common scale.
%   NORMALISE_GROUPS(run_settings) does the work of run_normalise_groups,
%   which sets the fields of run_settings (paths, mice, mousetypes,
%   mousetypes_list, selected_mice_idx_list, plot_verification_video,
%   channel, cohort_specs) and says what each one does.
%
%   Saves, in each cohort's folder, <channel>_4d_normalized<tag>.mat: the
%   normalised stack (AP x DV x ML x mouse, single) under the name
%   <channel>_4d_normalized, with norm_params (slope and intercept of each
%   mouse), current_mice, selected_mice_idx_list, channel, cohort_spec and
%   atlas_key; and <channel>_4d_normalized_bkgmask<tag>.mat: recomputed_bkg_mask_4d
%   (true for background), current_mice and selected_mice_idx_list. <tag> is the
%   cohort's, empty for an adult group and _P20 for young_P20.
%
%   A raw value of exactly zero is a voxel no section reached: the registered
%   volume is zero-filled beyond the first and last section, and where a section
%   does not cover the plane. The line would map it to -intercept / slope, a
%   large constant that becomes a bright slab once an absolute value is taken
%   wherever a cohort's coverage is partial. It is set to NaN instead, so a mean
%   over mice leaves that mouse out of that voxel; voxels with data are untouched.

% settings of run_normalise_groups, under the names the code below uses
paths = run_settings.paths;
mice = run_settings.mice;
mousetypes = run_settings.mousetypes;
mousetypes_list = run_settings.mousetypes_list;
selected_mice_idx_list = run_settings.selected_mice_idx_list;
plot_verification_video = run_settings.plot_verification_video;
channel = run_settings.channel;
cohort_specs = run_settings.cohort_specs;

%% Add paths

% the toolboxes are on the path from sep_setup_paths; only the atlas folder is
% added here
allenDir = paths.atlas;
addpath(allenDir)

%% Normalise each cohort

for ci = 1:numel(cohort_specs)
    tic

    % the cohort's mice, and its atlas on the grid of its registered volumes
    % (get_atlas_crop: the CCF for the adults, DeMBA for the P20 brains)
    [S, current_mouse_type, current_mice, subset_indices, A, AllenCrop, brainMask, ...
        allenDir] = select_cohort_mice(cohort_specs, ci, mousetypes_list, mousetypes, ...
        mice, selected_mice_idx_list);

    % the cohort's stack of the channel, selected mice only
    [base_dir, global_diagnostics_dir, data_4d] = load_cohort_stack(paths, ...
        current_mouse_type, channel, S, subset_indices);

    % planes pooled for the fit (every fifth), written in adult AP planes (900 in
    % the CCF crop) and scaled to the cohort's AP length: unchanged for the
    % adults, 994 planes for P20
    slices_range_for_norm = unique(round((1:5:900) * A.ap_scale));

    % planes drawn in the diagnostic figures
    slices_to_visualize_list = round([450, 500, 550] * A.ap_scale);

    % colour limit and histogram bins of the diagnostic figures
    plot_limit = 5000;
    hist_num_bins = 150;

    % the cortical voxels the fit uses
    cortex_mask_3d_all = cortex_reference_mask(allenDir, AllenCrop, brainMask);

    % a background mask for every plane of every mouse, and its trace
    [recomputed_bkg_mask_4d, median_vecs, area_vecs] = ...
        recompute_background_masks(data_4d, A);
    plot_background_trace(median_vecs, area_vecs, current_mouse_type, S, channel, ...
        base_dir);

    % pool the cortical tissue of the fit planes, and fit each mouse onto the
    % median mouse
    [cortex_samples_pooled, num_mice_subset] = pool_cortex_samples( ...
        slices_range_for_norm, cortex_mask_3d_all, subset_indices, data_4d, ...
        recomputed_bkg_mask_4d);
    [consensus_pixels_pooled, norm_params] = fit_to_consensus(cortex_samples_pooled, ...
        num_mice_subset);

    % the fit on all the pooled planes
    plot_global_diagnostic(cortex_samples_pooled, consensus_pixels_pooled, ...
        norm_params, num_mice_subset, current_mice, plot_limit, slices_range_for_norm, ...
        global_diagnostics_dir, channel);

    % six figures for each plane drawn
    plot_slice_diagnostics(slices_to_visualize_list, cortex_mask_3d_all, ...
        num_mice_subset, current_mice, recomputed_bkg_mask_4d, data_4d, plot_limit, ...
        hist_num_bins, norm_params, global_diagnostics_dir, channel);

    % a video of every plane, normalised against raw
    if plot_verification_video

        write_verification_video(current_mouse_type, channel, global_diagnostics_dir, ...
            data_4d, cortex_mask_3d_all, num_mice_subset, current_mice, ...
            recomputed_bkg_mask_4d, norm_params, plot_limit);

    end

    % normalise the whole volume and save it
    apply_and_save_normalisation(data_4d, norm_params, num_mice_subset, channel, ...
        base_dir, S, current_mice, selected_mice_idx_list, recomputed_bkg_mask_4d);

    % free memory before the next cohort
    clear data_4d_normalized recomputed_bkg_mask_4d save_struct

    toc

end

end

% ===== Local functions: cohort and masks =====

function [S, current_mouse_type, current_mice, subset_indices, A, AllenCrop, ...
    brainMask, allenDir] = select_cohort_mice(cohort_specs, ci, mousetypes_list, ...
    mousetypes, mice, selected_mice_idx_list)
% The cohort's spec, the mice to normalise, and the cohort's atlas on the grid
% of its registered volumes.

S = get_cohort_spec(cohort_specs{ci});
current_mouse_type = S.group;
mousetype_idx = find(strcmp(mousetypes_list, S.group));

if ~isempty(mousetype_idx)
    % an adult group: the mouse lists and selections of the driver
    all_mice_of_type_idx = find(strcmp(mousetypes, current_mouse_type));
    all_mice_of_type = mice(all_mice_of_type_idx);
    subset_indices = selected_mice_idx_list{mousetype_idx};
    current_mice = all_mice_of_type(subset_indices);
else
    % a young cohort: all the mice run_collect_by_group stacked, in that order
    current_mice = S.mice;
    subset_indices = 1:numel(current_mice);
end

% the cohort's atlas on the grid of its registered volumes
A = get_atlas_crop(S.atlas_key);
AllenCrop = A.annot;
brainMask = A.brainMask;
allenDir = A.csv_dir;

num_current = numel(current_mice);
num_mice_subset = num_current;

fprintf('Processing Group: %s | Selected %d mice: %s\n', ...
    current_mouse_type, num_current, strjoin(current_mice, ', '));
end

function [base_dir, global_diagnostics_dir, data_4d] = load_cohort_stack(paths, ...
    current_mouse_type, channel, S, subset_indices)
% The cohort's folders, and its stack of the channel with the selected mice only.

% the cohort's folder, and its diagnostics folder
base_dir = fullfile(paths.data, current_mouse_type);
global_diagnostics_dir = fullfile(paths.data, current_mouse_type, 'global_diagnostics');
if not(exist(global_diagnostics_dir, 'dir'))
    mkdir(global_diagnostics_dir)
end

fprintf('Loading data for %s...\n', current_mouse_type);

% load the stack (nano_4d or auto_4d) and keep the selected mice
input_var_name = [channel '_4d'];
input_file = fullfile(base_dir, [channel '_4d' S.tag '.mat']);
fprintf('Loading channel ''%s'' from %s ...\n', channel, input_file);
temp_data = load(input_file);
data_4d = temp_data.(input_var_name)(:, :, :, subset_indices);

% free the full stack
clear temp_data
end

function cortex_mask_3d_all = cortex_reference_mask(allenDir, AllenCrop, brainMask)
% Isocortex layers 1 to 5 without the retrosplenial, anterior cingulate and
% prelimbic areas: the voxels the normalisation is fitted on.

% isocortex, layers 1 to 5
target_regions = {'Isocortex'};
name_filter = {'Layer 1', 'Layer 2/3', 'Layer 4', 'Layer 5'};
cortex_mask_3d_inclusion_all = get_allen_region_mask(allenDir, AllenCrop, ...
    target_regions, brainMask, name_filter);

% the same layers of the three midline areas, left out
target_regions = {'Retrosplenial', 'Anterior cingulate area', 'Prelimbic area'};
name_filter = {'Layer 1', 'Layer 2/3', 'Layer 4', 'Layer 5'};
cortex_mask_3d_exclusion_all = get_allen_region_mask(allenDir, AllenCrop, ...
    target_regions, brainMask, name_filter);
cortex_mask_3d_all = and(cortex_mask_3d_inclusion_all, not(cortex_mask_3d_exclusion_all));
end

function [recomputed_bkg_mask_4d, median_vecs, area_vecs] = ...
    recompute_background_masks(data_4d, A)
% A background mask for every plane of every mouse, with percentile bounds that
% change along the AP axis, and the median and area of each mask.

fprintf('Recomputing background masks (per slice/mouse)...\n');
recomputed_bkg_mask_4d = false(size(data_4d));
total_slices = size(data_4d, 1);
median_vecs = NaN(total_slices, size(data_4d, 4));
area_vecs = NaN(total_slices, size(data_4d, 4));
for iii = 1:size(data_4d, 4)
    fprintf('  Processing Mouse %d / 5 ...\n', iii);
    for slice_idx_loop = 1:total_slices
        if mod(slice_idx_loop, 100) == 0
            fprintf('    -> Slice %d / %d\n', slice_idx_loop, total_slices);
        end

        % percentile bounds of the background by AP position; the breakpoints are
        % adult planes, so a longer young volume follows them at the same relative depth
        s_adult = slice_idx_loop / A.ap_scale;
        if s_adult < 100
            pmax_val = 95;
            pmin_val = 0;
        elseif s_adult >= 100 && s_adult < 250
            pmax_val = 75;
            pmin_val = 5;
        elseif s_adult >= 250 && s_adult < 500
            pmax_val = 55;
            pmin_val = 15;
        elseif s_adult >= 500 && s_adult < 700
            pmax_val = 45;
            pmin_val = 15;
        elseif s_adult >= 700 && s_adult < 725
            pmax_val = 55;
            pmin_val = 15;
        else
            pmax_val = 70;
            pmin_val = 20;
        end

        % the plane's background mask
        img_data = squeeze(data_4d(slice_idx_loop, :, :, iii));
        bool_diag_plot = false;
        bg_mask = select_background_pixels(img_data, pmin_val, pmax_val, bool_diag_plot);
        recomputed_bkg_mask_4d(slice_idx_loop, :, :, iii) = bg_mask;

        % its median intensity and area, for the trace figure
        median_vecs(slice_idx_loop, iii) = nanmedian(img_data(bg_mask)); %#ok<NANMEDIAN>
        area_vecs(slice_idx_loop, iii) = sum(bg_mask(:));
    end
end
end

function plot_background_trace(median_vecs, area_vecs, current_mouse_type, S, channel, ...
    base_dir)
% The figure of the background median and area per plane, saved in the cohort
% folder.

h_fig = figure('Visible', 'off', 'Name', ['Background_mask_diagnostics_trace_', ...
    current_mouse_type, S.tag, '_', channel], 'Color', 'w', ...
    'Position', [100 100 1400 600]);
nMice = size(median_vecs, 2);
nSlices = size(median_vecs, 1);
x_vals = 1:nSlices;

% one shade per mouse, from bright to dark magenta
c_bright = [1.0, 0.1, 1.0];
c_dark = [0.3, 0.0, 0.3];
colors = [linspace(c_bright(1), c_dark(1), nMice)', ...
    linspace(c_bright(2), c_dark(2), nMice)', ...
    linspace(c_bright(3), c_dark(3), nMice)'];

% left: the background median of every mouse, labelled M1, M2, ...
subplot(1, 2, 1);
hold on;
grid on;
box on;
title('Background median intensity per slice', 'FontSize', 12);
xlabel('Slice index');
ylabel('Median intensity');
for m = 1:nMice
    y_data = median_vecs(:, m);
    if all(isnan(y_data))
        continue;
    end
    plot(x_vals, y_data, 'Color', colors(m, :), 'LineWidth', 1.5);
    last_idx = find(~isnan(y_data), 1, 'last');
    if ~isempty(last_idx)
        text(x_vals(last_idx), 0.1*m*max(y_data), sprintf('  M%d', m), ...
            'Color', colors(m, :), 'FontSize', 9, ...
            'VerticalAlignment', 'middle');
    end
end
xlim([1 nSlices*1.1]);
ylim([0 max(y_data)*1.1]);

% right: the background area of every mouse
subplot(1, 2, 2);
hold on;
grid on;
box on;
title('Background mask area per slice', 'FontSize', 12);
xlabel('Slice index');
ylabel('Pixel count (sum)');
for m = 1:nMice
    y_data = area_vecs(:, m);
    if all(isnan(y_data))
        continue;
    end
    plot(x_vals, y_data, 'Color', colors(m, :), 'LineWidth', 1.5);
    last_idx = find(~isnan(y_data), 1, 'last');
    if ~isempty(last_idx)
        text(x_vals(last_idx), 0.1*m*max(y_data), sprintf('  M%d', m), ...
            'Color', colors(m, :), 'FontSize', 9, ...
            'VerticalAlignment', 'middle');
    end
end
xlim([1 nSlices*1.1]);
ylim([0 max(y_data)*1.1]);
sgtitle(strrep(['Background_mask_diagnostics_trace:_', current_mouse_type, '_(', ...
    channel, ')'], '_', ' '))

% save it, keeping the figure's background colour
set(h_fig, 'InvertHardcopy', 'off');
clean_fig_name = regexprep(h_fig.Name, '[^a-zA-Z0-9]', '_');
save_path_base = fullfile(base_dir, clean_fig_name);
fprintf('Saving diagnostic plot to: %s\n', save_path_base);
saveas(h_fig, [save_path_base '.fig']);
exportgraphics(h_fig, [save_path_base '.png'], 'Resolution', 300, ...
    'BackgroundColor', 'current');
end

% ===== Local functions: normalisation =====

function [cortex_samples_pooled, num_mice_subset] = pool_cortex_samples( ...
    slices_range_for_norm, cortex_mask_3d_all, subset_indices, data_4d, ...
    recomputed_bkg_mask_4d)
% The cortical tissue values of every mouse on the planes used for the fit, one
% column per mouse (NaN outside that mouse's tissue).

fprintf('Pooling pixels from %d slices for global normalization...\n', ...
    length(slices_range_for_norm));

% one row per pooled voxel, one column per mouse
cortex_samples_pooled = [];
num_mice_subset = numel(subset_indices);

for s_idx = slices_range_for_norm

    % the plane's cortical voxels
    mask_2d_slice = squeeze(cortex_mask_3d_all(s_idx, :, :));
    valid_pixels_indices = find(mask_2d_slice == 1);

    if isempty(valid_pixels_indices)
        continue;
    end

    % each mouse's values there, NaN on its background
    temp_samples = plane_cortex_samples(data_4d, recomputed_bkg_mask_4d, s_idx, ...
        valid_pixels_indices, num_mice_subset);

    % append them to the pool
    cortex_samples_pooled = [cortex_samples_pooled; temp_samples]; %#ok<AGROW>
end
end

function samples = plane_cortex_samples(data_4d, recomputed_bkg_mask_4d, plane, ...
    valid_pixels_indices, num_mice_subset)
% Each mouse's values on the given voxels of one plane, NaN on its background;
% one column per mouse.

samples = NaN(length(valid_pixels_indices), num_mice_subset);
for iii = 1:num_mice_subset
    current_slice = squeeze(data_4d(plane, :, :, iii));
    current_indiv_mask = ~squeeze(recomputed_bkg_mask_4d(plane, :, :, iii));
    vals = current_slice(valid_pixels_indices);
    is_valid_tissue = current_indiv_mask(valid_pixels_indices);
    vals(~is_valid_tissue) = NaN;
    samples(:, iii) = vals;
end
end

function [consensus_pixels_pooled, norm_params] = fit_to_consensus( ...
    cortex_samples_pooled, num_mice_subset)
% The median across mice as the reference, and each mouse's robust line against
% it (slope, intercept).

% the reference: the median mouse of each pooled voxel
fprintf('Calculating Global Median Consensus...\n');
consensus_pixels_pooled = nanmedian(cortex_samples_pooled, 2);

% each mouse's line against it
fprintf('Calculating Global Normalization Parameters...\n');
norm_params = zeros(num_mice_subset, 2);

for i = 1:num_mice_subset
    y_raw = cortex_samples_pooled(:, i);
    x_ref = consensus_pixels_pooled;

    % a robust line through the voxels both have, if there are more than 100
    % (least squares if the robust fit fails); slope 1, intercept 0 if fewer
    valid_idx = ~isnan(x_ref) & ~isnan(y_raw);
    if sum(valid_idx) > 100
        try
            % fitlm's second coefficient is the slope, the first the intercept
            mdl = fitlm(x_ref(valid_idx), y_raw(valid_idx), 'RobustOpts', 'on');
            p(1) = mdl.Coefficients.Estimate(2);
            p(2) = mdl.Coefficients.Estimate(1);
        catch
            p = polyfit(x_ref(valid_idx), y_raw(valid_idx), 1);
        end
    else
        p = [1, 0];
        warning('Not enough valid pixels to fit Mouse %d globally', i);
    end
    norm_params(i, :) = p;
    fprintf('  Mouse %d: Slope=%.2f, Int=%.2f\n', i, p(1), p(2));
end
end

function apply_and_save_normalisation(data_4d, norm_params, num_mice_subset, channel, ...
    base_dir, S, current_mice, selected_mice_idx_list, recomputed_bkg_mask_4d)
% Applies (raw - intercept) / slope to every mouse, sets the voxels no section
% reached to NaN, and saves the volume and the background masks.

fprintf('Applying normalization to full 4D volume and saving...\n');

% the normalised stack, in single
data_4d_normalized = zeros(size(data_4d), 'single');

for i = 1:num_mice_subset
    slope = norm_params(i, 1);
    intercept = norm_params(i, 2);

    % (raw - intercept) / slope
    data_4d_normalized(:, :, :, i) = (data_4d(:, :, :, i) - intercept) / slope;

    % a raw zero is a voxel no section reached, not tissue: NaN, so a mean over
    % mice leaves the mouse out there (the help says why)
    tmp = data_4d_normalized(:, :, :, i);
    tmp(data_4d(:, :, :, i) == 0) = NaN;
    data_4d_normalized(:, :, :, i) = tmp;
    clear tmp

    fprintf('  Applied normalization to Mouse %d/%d (Slope=%.2f, Int=%.2f)\n', ...
        i, num_mice_subset, slope, intercept);
end

% the stack is saved under its channel's name (nano_4d_normalized or
% auto_4d_normalized), which the readers ask for; the masks under one name for both
norm_var_name = [channel '_4d_normalized'];
save_filename = fullfile(base_dir, [channel '_4d_normalized' S.tag '.mat']);
save_filename_bis = fullfile(base_dir, [channel '_4d_normalized_bkgmask' S.tag '.mat']);

save_struct = struct();
save_struct.(norm_var_name) = data_4d_normalized;
save_struct.norm_params = norm_params;
save_struct.current_mice = current_mice;
save_struct.selected_mice_idx_list = selected_mice_idx_list;
save_struct.channel = channel;
save_struct.cohort_spec = S.label;
save_struct.atlas_key = S.atlas_key;
save(save_filename, '-struct', 'save_struct', '-v7.3');

save(save_filename_bis, 'recomputed_bkg_mask_4d', 'current_mice', ...
    'selected_mice_idx_list', '-v7.3');

fprintf('Successfully saved normalized volume to:\n  %s\n', save_filename);
end

% ===== Local functions: diagnostic figures =====

function plot_global_diagnostic(cortex_samples_pooled, consensus_pixels_pooled, ...
    norm_params, num_mice_subset, current_mice, plot_limit, slices_range_for_norm, ...
    global_diagnostics_dir, channel)
% Each mouse against the reference on the pooled planes, before and after the
% normalisation.

fprintf('Generating Global Diagnostic Plot (Pooled Data)...\n');

% at most 50000 points, drawn at random with a fixed seed: millions of
% transparent dots draw very slowly
max_plot_points = 50000;
total_pooled = size(cortex_samples_pooled, 1);
if total_pooled > max_plot_points
    rng(42);
    idx_sub = randperm(total_pooled, max_plot_points);
else
    idx_sub = 1:total_pooled;
end

figure('Visible', 'off', 'Name', 'Global Normalization Diagnostic (Pooled)', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [-0.05 -0.05 0.95 0.95]);

for i = 1:num_mice_subset
    mouse_name = strrep(current_mice{i}, '_', ' ');

    % the drawn points
    y_raw = cortex_samples_pooled(idx_sub, i);
    x_ref = consensus_pixels_pooled(idx_sub);

    % the mouse's slope and intercept
    p = norm_params(i, :);

    % the points normalised
    y_norm = (y_raw - p(2)) / p(1);

    % top row: raw against the median mouse, with the fitted line and identity
    subplot(2, num_mice_subset, i);
    scatter(x_ref, y_raw, 2, 'k', 'filled', 'MarkerFaceAlpha', 0.1);
    hold on;
    axis equal

    valid_idx = ~isnan(x_ref) & ~isnan(y_raw);
    if any(valid_idx)
        min_x = min(x_ref(valid_idx));
        max_x = max(x_ref(valid_idx));
        x_grid = linspace(min_x, max_x, 100);
        y_fit = p(1)*x_grid + p(2);
        plot(x_grid, y_fit, 'b-', 'LineWidth', 2);
    end
    plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1.5);
    text(0.05*plot_limit, 0.85*plot_limit, ...
        sprintf('Slope: %.2f\nInt: %.0f', p(1), p(2)), ...
        'Color', 'b', 'FontSize', 9, 'FontWeight', 'bold');

    title(mouse_name, 'FontSize', 11, 'FontWeight', 'bold');
    xlim([0 plot_limit]);
    ylim([0 plot_limit]);
    grid on;
    set(gca, 'XTickLabel', []);
    if i == 1
        ylabel({'Raw Intensity'; '(Pooled Subset)'}, 'FontSize', 10);
    else
        set(gca, 'YTickLabel', []);
    end

    % bottom row: normalised against the median mouse, with identity
    subplot(2, num_mice_subset, i + num_mice_subset);
    scatter(x_ref, y_norm, 2, 'filled', 'MarkerFaceColor', 'b', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
    hold on;
    axis equal
    plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1.5);

    R = corrcoef(x_ref, y_norm, 'Rows', 'complete');
    if numel(R) > 1
        r_val = R(1, 2);
    else
        r_val = NaN;
    end
    text(0.05*plot_limit, 0.9*plot_limit, sprintf('R = %.2f', r_val), ...
        'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

    title(mouse_name, 'FontSize', 11, 'FontWeight', 'bold');
    xlim([0 plot_limit]);
    ylim([0 plot_limit]);
    grid on;
    xlabel('Global Median Intensity');
    if i == 1
        ylabel({'Normalized'; '(Corrected)'}, 'FontSize', 10);
    else
        set(gca, 'YTickLabel', []);
    end
end

sgtitle(['Global Normalization Diagnostic (Pooled Slices: ' ...
    num2str(min(slices_range_for_norm)) '-' num2str(max(slices_range_for_norm)) ')'], ...
    'FontSize', 14);

% save it in normalization_checks_<channel>\, named after the figure (no plane
% number: it pools them all)
save_check_figure(global_diagnostics_dir, channel, []);
end

function plot_slice_diagnostics(slices_to_visualize_list, cortex_mask_3d_all, ...
    num_mice_subset, current_mice, recomputed_bkg_mask_4d, data_4d, plot_limit, ...
    hist_num_bins, norm_params, global_diagnostics_dir, channel)
% Six diagnostic figures for each of the planes to visualise.

fprintf('Starting visualization loop for %d slices...\n', ...
    length(slices_to_visualize_list));

for viz_idx = 1:length(slices_to_visualize_list)

    slice_to_plot = slices_to_visualize_list(viz_idx);
    fprintf('  Visualizing Slice %d...\n', slice_to_plot);

    % close the figures of the previous plane
    close all

    % the plane's cortical voxels
    cortex_slice_mask = squeeze(cortex_mask_3d_all(slice_to_plot, :, :));
    mask_2d_slice = squeeze(cortex_mask_3d_all(slice_to_plot, :, :));
    valid_pixels_indices = find(mask_2d_slice == 1);

    % plot 1: every mouse's plane, dimmed outside its cortical tissue
    plot_individual_slices(num_mice_subset, current_mice, recomputed_bkg_mask_4d, ...
        slice_to_plot, cortex_slice_mask, data_4d, plot_limit, global_diagnostics_dir, ...
        channel);

    % each mouse's values on the plane's cortical voxels, NaN on its background,
    % for the scatter plots
    cortex_samples = plane_cortex_samples(data_4d, recomputed_bkg_mask_4d, ...
        slice_to_plot, valid_pixels_indices, num_mice_subset);

    % plot 2: every pair of mice against each other
    plot_pairwise_mosaic(cortex_samples, num_mice_subset, current_mice, hist_num_bins, ...
        plot_limit, slice_to_plot, global_diagnostics_dir, channel);

    % plot 3: the median mouse of the plane; and the plane's own reference, the
    % median mouse of each cortical voxel
    slice_stack = squeeze(data_4d(slice_to_plot, :, :, 1:num_mice_subset));
    consensus_slice = median(slice_stack, 3);
    consensus_pixels = nanmedian(cortex_samples, 2); %#ok<NANMEDIAN>

    plot_median_consensus(consensus_slice, recomputed_bkg_mask_4d, slice_to_plot, ...
        num_mice_subset, mask_2d_slice, plot_limit, global_diagnostics_dir, channel);

    % plot 4: every mouse against the median mouse
    plot_individual_vs_median(consensus_pixels, cortex_samples, num_mice_subset, ...
        current_mice, plot_limit, slice_to_plot, global_diagnostics_dir, channel);

    % plot 5: the global fit on this plane, before and after the normalisation
    plot_slice_norm_diagnostic(cortex_samples, consensus_pixels, norm_params, ...
        num_mice_subset, current_mice, plot_limit, slice_to_plot, ...
        global_diagnostics_dir, channel);

    % plot 6: every mouse's plane, normalised and raw
    plot_norm_vs_raw(data_4d, recomputed_bkg_mask_4d, mask_2d_slice, norm_params, ...
        num_mice_subset, current_mice, plot_limit, slice_to_plot, ...
        global_diagnostics_dir, channel);

    fprintf('  Done with Slice %d\n', slice_to_plot);
end
end

function plot_individual_slices(num_mice_subset, current_mice, recomputed_bkg_mask_4d, ...
    slice_to_plot, cortex_slice_mask, data_4d, plot_limit, global_diagnostics_dir, ...
    channel)
% Plot 1: every mouse's plane, dimmed outside its cortical tissue.

for iii = 1:num_mice_subset
    figure('Visible', 'off', 'Name', ['Individual_Slice_' current_mice{iii}], ...
        'Color', 'k');

    % black at 0.75 opacity outside the cortex or the mouse's tissue
    mask_data = not(squeeze(recomputed_bkg_mask_4d(slice_to_plot, :, :, iii)));
    overlay_alpha = zeros(size(cortex_slice_mask));
    overlay_alpha(cortex_slice_mask .* mask_data == 0) = 0.75;

    img_data = squeeze(data_4d(slice_to_plot, :, :, iii));
    imagesc(img_data);
    colormap(sep_palette('intensity'));
    clim([0, plot_limit]);
    hold on;

    black_overlay = zeros(size(img_data));
    h_ov = imagesc(black_overlay);
    set(h_ov, 'AlphaData', overlay_alpha);

    axis image;
    axis off;
    set(gca, 'Color', 'k');
    mouse_name = strrep(current_mice{iii}, '_', ' ');
    t = title(mouse_name);
    set(t, 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');

    cb = colorbar;
    cb.Label.String = 'Intensity (a.u.)';
    cb.Color = 'w';
    cb.Label.Color = 'w';
    hold off;

    % save it in normalization_checks_<channel>\, named after the figure and plane
    save_check_figure(global_diagnostics_dir, channel, slice_to_plot);
end
end

function plot_pairwise_mosaic(cortex_samples, num_mice_subset, current_mice, ...
    hist_num_bins, plot_limit, slice_to_plot, global_diagnostics_dir, channel)
% Plot 2: every pair of mice against each other on the plane's cortex.

figure('Visible', 'off', 'Name', 'Pairwise Cortex Intensity Comparison', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

for row = 1:num_mice_subset
    for col = 1:num_mice_subset
        idx = (row - 1) * num_mice_subset + col;
        subplot(num_mice_subset, num_mice_subset, idx);
        name_row = strrep(current_mice{row}, '_', ' ');
        name_col = strrep(current_mice{col}, '_', ' ');

        % the mouse's histogram on the diagonal, a pair of mice elsewhere
        if row == col
            histogram(cortex_samples(:, row), hist_num_bins, 'EdgeColor', 'none', ...
                'FaceColor', 'r', 'BinLimits', [0, plot_limit]);
            title(name_row, 'FontWeight', 'bold', 'FontSize', 8);
            grid on;
            xlim([0 plot_limit]);
            yticklabels([]);
        else
            x_data = cortex_samples(:, col);
            y_data = cortex_samples(:, row);
            scatter(x_data, y_data, 2, 'filled', 'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
            hold on;
            plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1);
            R = corrcoef(x_data, y_data, 'Rows', 'complete');
            if numel(R) > 1
                r_val = R(1, 2);
            else
                r_val = NaN;
            end
            text(0.05 * plot_limit, 0.9 * plot_limit, sprintf('R = %.2f', r_val), ...
                'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold');
            xlim([0 plot_limit]);
            ylim([0 plot_limit]);
            grid on;
            if col == 1
                ylabel(name_row, 'FontSize', 8, 'FontWeight', 'bold');
            else
                set(gca, 'YTickLabel', []);
            end
            if row == num_mice_subset
                xlabel(name_col, 'FontSize', 8, 'FontWeight', 'bold');
            else
                set(gca, 'XTickLabel', []);
            end
        end
    end
end
sgtitle(['Cortex Pixel Intensity Comparison (Slice ' num2str(slice_to_plot) ...
    ') - Range [0, ' num2str(plot_limit) ']']);

% save it in normalization_checks_<channel>\, named after the figure and plane
save_check_figure(global_diagnostics_dir, channel, slice_to_plot);
end

function plot_median_consensus(consensus_slice, recomputed_bkg_mask_4d, slice_to_plot, ...
    num_mice_subset, mask_2d_slice, plot_limit, global_diagnostics_dir, channel)
% Plot 3: the median mouse of the plane.

figure('Visible', 'off', 'Name', 'Median Consensus Slice');
imagesc(consensus_slice);
colormap(sep_palette('intensity'));
clim([0, plot_limit]);
hold on;

% black at 0.5 opacity outside the cortex, and where no mouse has tissue
slice_bkg_stack = squeeze(recomputed_bkg_mask_4d(slice_to_plot, :, :, 1:num_mice_subset));
at_least_one_tissue = any(~slice_bkg_stack, 3);
median_overlay_alpha = zeros(size(consensus_slice));
is_valid_region = (mask_2d_slice == 1) & at_least_one_tissue;
median_overlay_alpha(~is_valid_region) = 0.5;

black_overlay = zeros(size(consensus_slice));
h_ov = imagesc(black_overlay);
set(h_ov, 'AlphaData', median_overlay_alpha);
axis image;
axis off;
set(gca, 'Color', 'k');
t = title(['Median Mouse (Slice ' num2str(slice_to_plot) ')']);
set(t, 'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
cb = colorbar;
cb.Label.String = 'Intensity (a.u.)';
cb.Color = 'k';
cb.Label.Color = 'k';
hold off;

% save it in normalization_checks_<channel>\, named after the figure and plane
save_check_figure(global_diagnostics_dir, channel, slice_to_plot);
end

function plot_individual_vs_median(consensus_pixels, cortex_samples, num_mice_subset, ...
    current_mice, plot_limit, slice_to_plot, global_diagnostics_dir, channel)
% Plot 4: every mouse against the median mouse on the plane's cortex.

figure('Visible', 'off', 'Name', 'Individual vs Median Comparison', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.1 0.2 0.8 0.4]);

for i = 1:num_mice_subset
    subplot(1, num_mice_subset, i);
    mouse_name = strrep(current_mice{i}, '_', ' ');

    x_data = consensus_pixels;
    y_data = cortex_samples(:, i);

    scatter(x_data, y_data, 2, 'filled', 'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
    hold on;
    plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1);
    axis equal
    R = corrcoef(x_data, y_data, 'Rows', 'complete');
    if numel(R) > 1
        r_val = R(1, 2);
    else
        r_val = NaN;
    end
    text(0.05 * plot_limit, 0.9 * plot_limit, sprintf('R = %.2f', r_val), ...
        'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
    xlim([0 plot_limit]);
    ylim([0 plot_limit]);
    grid on;
    title(mouse_name, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Median Intensity');
    if i == 1
        ylabel('Individual Intensity');
    else
        set(gca, 'YTickLabel', []);
    end
end
sgtitle(['Individual Mice vs. Group Median (Slice ' num2str(slice_to_plot) ')']);

% save it in normalization_checks_<channel>\, named after the figure and plane
save_check_figure(global_diagnostics_dir, channel, slice_to_plot);
end

function plot_slice_norm_diagnostic(cortex_samples, consensus_pixels, norm_params, ...
    num_mice_subset, current_mice, plot_limit, slice_to_plot, global_diagnostics_dir, ...
    channel)
% Plot 5: the global fit of every mouse on the plane, before and after the
% normalisation.

cortex_samples_norm = zeros(size(cortex_samples));

figure('Visible', 'off', 'Name', 'Normalization Diagnostic and Verification', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [-0.05 -0.05 0.95 0.95]);

for i = 1:num_mice_subset
    mouse_name = strrep(current_mice{i}, '_', ' ');
    y_raw = cortex_samples(:, i);
    x_ref = consensus_pixels;

    % the mouse's global slope and intercept, applied to the plane
    p = norm_params(i, :);
    cortex_samples_norm(:, i) = (y_raw - p(2)) / p(1);

    % top row: raw against the plane's median mouse, with the global line
    subplot(2, num_mice_subset, i);
    scatter(x_ref, y_raw, 2, 'k', 'filled', 'MarkerFaceAlpha', 0.1);
    hold on;
    axis equal

    valid_idx = ~isnan(x_ref) & ~isnan(y_raw);
    if any(valid_idx)
        min_x = min(x_ref(valid_idx));
        max_x = max(x_ref(valid_idx));
        x_grid = linspace(min_x, max_x, 100);
        y_fit = p(1)*x_grid + p(2);
        plot(x_grid, y_fit, 'b-', 'LineWidth', 2);
    end
    plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1.5);
    text(0.05*plot_limit, 0.85*plot_limit, ...
        sprintf('Slope: %.2f\nInt: %.0f', p(1), p(2)), ...
        'Color', 'b', 'FontSize', 9, 'FontWeight', 'bold');

    title(mouse_name, 'FontSize', 11, 'FontWeight', 'bold');
    xlim([0 plot_limit]);
    ylim([0 plot_limit]);
    grid on;
    set(gca, 'XTickLabel', []);
    if i == 1
        ylabel({'Raw Intensity'; '(Individual)'}, 'FontSize', 10);
    else
        set(gca, 'YTickLabel', []);
    end

    % bottom row: normalised against the plane's median mouse
    subplot(2, num_mice_subset, i + num_mice_subset);
    scatter(x_ref, cortex_samples_norm(:, i), 2, 'filled', 'MarkerFaceColor', 'b', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
    hold on;
    axis equal
    plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1.5);

    R = corrcoef(x_ref, cortex_samples_norm(:, i), 'Rows', 'complete');
    if numel(R) > 1
        r_val = R(1, 2);
    else
        r_val = NaN;
    end
    text(0.05*plot_limit, 0.9*plot_limit, sprintf('R = %.2f', r_val), 'Color', 'r', ...
        'FontSize', 10, 'FontWeight', 'bold');

    title(mouse_name, 'FontSize', 11, 'FontWeight', 'bold');
    xlim([0 plot_limit]);
    ylim([0 plot_limit]);
    grid on;
    xlabel('Median Intensity (Target)');
    if i == 1
        ylabel({'Normalized'; '(Corrected)'}, 'FontSize', 10);
    else
        set(gca, 'YTickLabel', []);
    end
end
sgtitle(['Normalization diagnostic (Slice ' num2str(slice_to_plot) ')'], 'FontSize', 14);

% save it in normalization_checks_<channel>\, named after the figure and plane
save_check_figure(global_diagnostics_dir, channel, slice_to_plot);
close all
end

function plot_norm_vs_raw(data_4d, recomputed_bkg_mask_4d, mask_2d_slice, norm_params, ...
    num_mice_subset, current_mice, plot_limit, slice_to_plot, global_diagnostics_dir, ...
    channel)
% Plot 6: every mouse's plane normalised (top) and raw (bottom).

figure('Visible', 'off', 'Name', 'Visual Verification Normalized vs Raw', ...
    'Color', 'k', 'Units', 'normalized', 'Position', [0.1 0.05 0.8 0.8]);

for i = 1:num_mice_subset
    mouse_name = strrep(current_mice{i}, '_', ' ');

    % black at 0.5 opacity outside the cortex or the mouse's tissue
    img_raw = squeeze(data_4d(slice_to_plot, :, :, i));
    current_bg_mask = squeeze(recomputed_bkg_mask_4d(slice_to_plot, :, :, i));
    is_valid_tissue = (mask_2d_slice == 1) & (current_bg_mask == 0);

    overlay_alpha = zeros(size(img_raw));
    overlay_alpha(~is_valid_tissue) = 0.5;

    % the plane normalised
    slope = norm_params(i, 1);
    intercept = norm_params(i, 2);
    img_norm = (img_raw - intercept) / slope;

    % top row: normalised
    subplot(2, num_mice_subset, i);
    imagesc(img_norm);
    colormap(sep_palette('intensity'));
    clim([0, plot_limit]);
    axis image;
    axis off;
    set(gca, 'Color', 'k');
    hold on;
    h_ov1 = imagesc(zeros(size(img_norm)));
    set(h_ov1, 'AlphaData', overlay_alpha);
    t = title(['Norm: ' mouse_name]);
    set(t, 'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold');
    if i == num_mice_subset
        cb = colorbar;
        cb.Label.String = 'Normalized Intensity';
        cb.Color = 'w';
        cb.Label.Color = 'w';
        cb.Position = [0.92 0.55 0.01 0.35];
    end
    hold off;

    % bottom row: raw
    subplot(2, num_mice_subset, i + num_mice_subset);
    imagesc(img_raw);
    colormap(sep_palette('intensity'));
    clim([0, plot_limit]);
    axis image;
    axis off;
    set(gca, 'Color', 'k');
    hold on;
    h_ov2 = imagesc(zeros(size(img_raw)));
    set(h_ov2, 'AlphaData', overlay_alpha);
    t = title(['Raw: ' mouse_name]);
    set(t, 'Color', 'w', 'FontSize', 11, 'FontWeight', 'normal');
    if i == num_mice_subset
        cb = colorbar;
        cb.Label.String = 'Raw Intensity';
        cb.Color = 'w';
        cb.Label.Color = 'w';
        cb.Position = [0.92 0.1 0.01 0.35];
    end
    hold off;
end
sgtitle(['Comparison: normalized (top) vs. raw (bottom) - Fixed scale [0, ' ...
    num2str(plot_limit) ']'], 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');

% save it in normalization_checks_<channel>\, named after the figure and plane
save_check_figure(global_diagnostics_dir, channel, slice_to_plot);
close all
end

function write_verification_video(current_mouse_type, channel, ...
    global_diagnostics_dir, data_4d, cortex_mask_3d_all, num_mice_subset, ...
    current_mice, recomputed_bkg_mask_4d, norm_params, plot_limit)
% Video of every plane with cortex: every mouse normalised (top) and raw
% (bottom).

fprintf('Generating Normalization Verification Video (this may take a while)...\n');

% the video file, in normalization_checks_<channel>\
video_filename = ['Normalization_Verification_' current_mouse_type '_' channel '.mp4'];

save_video_path = checks_folder(global_diagnostics_dir, channel);

full_video_path = fullfile(save_video_path, video_filename);

% 15 frames per second, quality 95
vidObj = VideoWriter(full_video_path, 'MPEG-4');
vidObj.FrameRate = 15;
vidObj.Quality = 95;
open(vidObj);

% every plane of the volume
slices_to_video = 1:size(data_4d, 1);

for s_idx = slices_to_video

    % skip the planes with no cortex in the atlas
    mask_2d_slice = squeeze(cortex_mask_3d_all(s_idx, :, :));
    if sum(mask_2d_slice(:)) == 0
        if mod(s_idx, 100) == 0
            fprintf('  Skipping empty slice %d\n', s_idx);
        end
        continue;
    end

    % an invisible figure on black, kept black in the frame
    fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], ...
        'Color', 'k');
    set(fh, 'InvertHardcopy', 'off');

    for i = 1:num_mice_subset
        mouse_name = strrep(current_mice{i}, '_', ' ');

        % the mouse's plane and background
        img_raw = squeeze(data_4d(s_idx, :, :, i));
        current_bg_mask = squeeze(recomputed_bkg_mask_4d(s_idx, :, :, i));

        % dimmed unless in the cortex and in the mouse's tissue
        is_valid_tissue = (mask_2d_slice == 1) & (current_bg_mask == 0);
        overlay_alpha = zeros(size(img_raw));
        overlay_alpha(~is_valid_tissue) = 0.5;

        % the plane normalised
        slope = norm_params(i, 1);
        intercept = norm_params(i, 2);
        img_norm = (img_raw - intercept) / slope;

        % top row: normalised
        subplot(2, num_mice_subset, i);
        imagesc(img_norm);
        colormap(sep_palette('intensity'));
        clim([0, plot_limit]);
        axis image;
        axis off;
        set(gca, 'Color', 'k');
        hold on;

        h_ov1 = imagesc(zeros(size(img_norm)));
        set(h_ov1, 'AlphaData', overlay_alpha);

        t = title(['Norm: ' mouse_name]);
        set(t, 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');

        if i == num_mice_subset
            cb = colorbar;
            cb.Label.String = 'Norm Int';
            cb.Color = 'w';
            cb.Label.Color = 'w';
            cb.Position = [0.92 0.55 0.01 0.35];
        end
        hold off;

        % bottom row: raw
        subplot(2, num_mice_subset, i + num_mice_subset);
        imagesc(img_raw);
        colormap(sep_palette('intensity'));
        clim([0, plot_limit]);
        axis image;
        axis off;
        set(gca, 'Color', 'k');
        hold on;

        h_ov2 = imagesc(zeros(size(img_raw)));
        set(h_ov2, 'AlphaData', overlay_alpha);

        t = title(['Raw: ' mouse_name]);
        set(t, 'Color', 'w', 'FontSize', 10, 'FontWeight', 'normal');

        if i == num_mice_subset
            cb = colorbar;
            cb.Label.String = 'Raw Int';
            cb.Color = 'w';
            cb.Label.Color = 'w';
            cb.Position = [0.92 0.1 0.01 0.35];
        end
        hold off;
    end

    sgtitle(['Slice ' num2str(s_idx) ' - Norm vs Raw - Scale [0 ' num2str(plot_limit) ...
        ']'], 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');

    % write the frame
    frame = getframe(fh);
    writeVideo(vidObj, frame);
    close(fh);

    if mod(s_idx, 20) == 0
        fprintf('  Video Frame: Slice %d written...\n', s_idx);
    end
end

close(vidObj);
fprintf('Video saved successfully: %s\n', full_video_path);
end

function save_check_figure(global_diagnostics_dir, channel, slice_to_plot)
% Save the current figure as .fig and .png in normalization_checks_<channel>\,
% named after the figure, and after the plane when slice_to_plot is not empty.

save_output_dir = checks_folder(global_diagnostics_dir, channel);
fig_handle = gcf;
set(fig_handle, 'InvertHardcopy', 'off');
clean_fig_name = regexprep(fig_handle.Name, '[^a-zA-Z0-9]', '_');
if isempty(clean_fig_name)
    clean_fig_name = 'Untitled_Figure';
end
if isempty(slice_to_plot)
    filename_base = clean_fig_name;
else
    filename_base = sprintf('%s_Slice%d', clean_fig_name, slice_to_plot);
end
saveas(fig_handle, fullfile(save_output_dir, [filename_base '.fig']));
exportgraphics(fig_handle, fullfile(save_output_dir, [filename_base '.png']), ...
    'Resolution', 300, 'BackgroundColor', 'current');
end

function save_output_dir = checks_folder(global_diagnostics_dir, channel)
% The folder normalization_checks_<channel>\ in the diagnostics folder, made if it
% is missing.

if exist('global_diagnostics_dir', 'var')
    save_output_dir = fullfile(global_diagnostics_dir, ['normalization_checks_' channel]);
else
    save_output_dir = fullfile(pwd, ['normalization_checks_' channel]);
end
if ~exist(save_output_dir, 'dir')
    mkdir(save_output_dir);
end
end
