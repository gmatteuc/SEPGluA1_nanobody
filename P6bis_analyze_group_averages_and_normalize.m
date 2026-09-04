clear all
close all
clc

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();


% /// Pipeline script #6: perform final normalization across mice in each group and plot averages ///

%%  Set user-defined parameters

% Set mice and mousetypes
mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', 'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1','MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', 'MG725_Gria1', 'MG727_Gria1'};
mousetypes = {'rws','rws','rws','rws','rws','naive','naive','naive','naive','naive','behavior','behavior','behavior','behavior','behavior','behavior','behavior'};
mousetypes_list = {'rws','naive','behavior'};

% Set list of selected mice
selected_mice_idx_list{1} = 1:5; % rws
selected_mice_idx_list{2} = 1:5; % naive
selected_mice_idx_list{3} = [1,3,4,5]; % behavior (excluding 3 bad mice)

% Set plotting and saving
plot_verification_video = false;

% Channel to normalize: 'nano' (default, surface GluA1) or 'auto' (autofluorescence control).
% Loads <channel>_4d.mat and saves <channel>_4d_normalized.mat in each cohort folder.
% Pipeline scripts downstream (P7bis, P8, P9) read whichever channel they're configured for.
channel = 'nano'; % 'auto'

% Cohorts to normalise, as specs: an adult group by name ('rws', 'naive',
% 'behavior'), or the young group with an age ('young_P20'). An adult group
% keeps the mouse lists and selections above, exactly as before. A young
% cohort takes the mice P5 stacked (collected_mice<tag>.mat), all of them, and
% its own atlas. See get_cohort_spec.
cohort_specs = {'rws', 'naive'};
% Batch runs can pick the cohorts without editing this file:
%   set SEP_COHORT_SPECS=young_P20   (comma-separated for several)
if ~isempty(getenv('SEP_COHORT_SPECS'))
    cohort_specs = strtrim(strsplit(getenv('SEP_COHORT_SPECS'), ','));
end

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

%% Atlas

% Loaded per cohort inside the loop (get_atlas_crop): the adults are on the
% CCF, the P20 brains on DeMBA, and each comes out on the grid of its own
% registered volumes. For an adult group the arrays are the same, byte for
% byte, as the inline crop this block used to do.

%% Plot videos per mouse type

for ci = 1:numel(cohort_specs)
    tic

    % --- Fetch single mouse data for current cohort ---

    S = get_cohort_spec(cohort_specs{ci});
    current_mouse_type = S.group;
    mousetype_idx = find(strcmp(mousetypes_list, S.group));

    if ~isempty(mousetype_idx)
        % An adult group: the lists and selections above, as they always were
        all_mice_of_type_idx = find(strcmp(mousetypes, current_mouse_type));
        all_mice_of_type = mice(all_mice_of_type_idx);
        subset_indices = selected_mice_idx_list{mousetype_idx};
        current_mice = all_mice_of_type(subset_indices);
    else
        % A young cohort: whatever P5 stacked, in that order, all of it
        current_mice   = S.mice;
        subset_indices = 1:numel(current_mice);
    end

    % The cohort's atlas on the grid of its registered volumes
    A = get_atlas_crop(S.atlas_key);
    AllenCrop = A.annot;
    brainMask = A.brainMask;
    allenDir  = A.csv_dir;

    num_current = numel(current_mice);
    num_mice_subset = num_current; % Update subset count for later loops

    fprintf('Processing Group: %s | Selected %d mice: %s\n', ...
        current_mouse_type, num_current, strjoin(current_mice, ', '));

    % Get base_dirs
    base_dir = fullfile(paths.data, current_mouse_type);
    global_diagnostics_dir = fullfile(paths.data, current_mouse_type, 'global_diagnostics');
    if not(exist(global_diagnostics_dir, 'dir'))
        mkdir(global_diagnostics_dir)
    end

    % Load needed 4D data
    fprintf('Loading data for %s...\n', current_mouse_type);

    % Filter the 4D matrices to keep only selected mice
    input_var_name = [channel '_4d'];                   % e.g. 'nano_4d' or 'auto_4d'
    input_file     = fullfile(base_dir, [channel '_4d' S.tag '.mat']);
    fprintf('Loading channel ''%s'' from %s ...\n', channel, input_file);
    temp_data = load(input_file);
    data_4d = temp_data.(input_var_name)(:, :, :, subset_indices);

    % Clear temp variables to free memory
    clear temp_data

    % ----------- Set user-defined parameters for the plot -----------

    % Range of slices to use for CALCULATING the normalization (pooling pixels)
    % Written in adult AP planes (900 in the CCF crop) and scaled to this
    % cohort's AP length: unchanged for the adults, 994 planes for P20.
    slices_range_for_norm = unique(round((1:5:900) * A.ap_scale));

    % List of slices to VISUALIZE and Save (Diagnostics)
    slices_to_visualize_list = round([450,500,550] * A.ap_scale);

    plot_limit = 5000;
    hist_num_bins = 150;

    % ----------- Get atlas mask -----------
    target_regions = {'Isocortex'};
    name_filter = {'Layer 1', 'Layer 2/3', 'Layer 4', 'Layer 5'};
    cortex_mask_3d_inclusion_all = get_allen_region_mask(allenDir, AllenCrop, target_regions, brainMask, name_filter);
    target_regions = {'Retrosplenial', 'Anterior cingulate area', 'Prelimbic area'};
    name_filter = {'Layer 1', 'Layer 2/3', 'Layer 4', 'Layer 5'};
    cortex_mask_3d_exclusion_all = get_allen_region_mask(allenDir, AllenCrop, target_regions, brainMask, name_filter);
    cortex_mask_3d_all=and(cortex_mask_3d_inclusion_all,not(cortex_mask_3d_exclusion_all));

    % ----------- Recompute 4D masks for entire registered volume -----------

    fprintf('Recomputing background masks (per slice/mouse)...\n');
    recomputed_bkg_mask_4d = false(size(data_4d));
    total_slices = size(data_4d, 1);
    median_vecs = NaN(total_slices,size(data_4d,4));
    area_vecs = NaN(total_slices,size(data_4d,4));
    for iii = 1:size(data_4d,4)
        fprintf('  Processing Mouse %d / 5 ...\n', iii);
        for slice_idx_loop = 1:total_slices
            if mod(slice_idx_loop, 100) == 0
                fprintf('    -> Slice %d / %d\n', slice_idx_loop, total_slices);
            end
            % Set settings. The breakpoints are adult AP planes; the index
            % is brought back to that scale so the same schedule applies
            % at the same relative depth in a longer young volume.
            s_adult = slice_idx_loop / A.ap_scale;
            if s_adult<100
                pmax_val=95;
                pmin_val=0;
            elseif  s_adult>=100 && s_adult<250
                pmax_val=75;
                pmin_val=5;
            elseif s_adult>=250 && s_adult<500
                pmax_val=55;
                pmin_val=15;
            elseif s_adult>=500 && s_adult<700
                pmax_val=45;
                pmin_val=15;
            elseif s_adult>=700 && s_adult<725
                pmax_val=55;
                pmin_val=15;
            else
                pmax_val=70;
                pmin_val=20;
            end
            % Get data
            img_data = squeeze(data_4d(slice_idx_loop,:,:,iii));
            % Compute background mask
            bool_diag_plot = false;
            bg_mask = select_background_pixels(img_data, pmin_val, pmax_val, bool_diag_plot);
            recomputed_bkg_mask_4d(slice_idx_loop,:,:,iii) = bg_mask;
            % Store diagnostics
            median_vecs(slice_idx_loop,iii) = nanmedian(img_data(bg_mask)); %#ok<NANMEDIAN>
            area_vecs(slice_idx_loop,iii) = sum(bg_mask(:));
        end
    end

    % Plot diagnostics
    h_fig = figure('Visible', 'off', 'Name', ['Background_mask_diagnostics_trace_',current_mouse_type,S.tag,'_',channel], 'Color', 'w', 'Position', [100 100 1400 600]);
    nMice = size(median_vecs, 2);
    nSlices = size(median_vecs, 1);
    x_vals = 1:nSlices;
    c_bright = [1.0, 0.1, 1.0];
    c_dark   = [0.3, 0.0, 0.3];
    colors = [linspace(c_bright(1), c_dark(1), nMice)', ...
        linspace(c_bright(2), c_dark(2), nMice)', ...
        linspace(c_bright(3), c_dark(3), nMice)'];
    subplot(1, 2, 1);
    hold on; grid on; box on;
    title('Background median intensity per slice', 'FontSize', 12);
    xlabel('Slice index'); ylabel('Median intensity');
    for m = 1:nMice
        y_data = median_vecs(:, m);
        if all(isnan(y_data)), continue; end
        plot(x_vals, y_data, 'Color', colors(m,:), 'LineWidth', 1.5);
        last_idx = find(~isnan(y_data), 1, 'last');
        if ~isempty(last_idx)
            text(x_vals(last_idx), 0.1*m*max(y_data), sprintf('  M%d', m), ...
                'Color', colors(m,:), 'FontSize', 9, ...
                'VerticalAlignment', 'middle');
        end
    end
    xlim([1 nSlices*1.1]);
    ylim([0 max(y_data)*1.1]);
    subplot(1, 2, 2);
    hold on; grid on; box on;
    title('Background mask area per slice', 'FontSize', 12);
    xlabel('Slice index'); ylabel('Pixel count (sum)');
    for m = 1:nMice
        y_data = area_vecs(:, m);
        if all(isnan(y_data)), continue; end
        plot(x_vals, y_data, 'Color', colors(m,:), 'LineWidth', 1.5);
        last_idx = find(~isnan(y_data), 1, 'last');
        if ~isempty(last_idx)
            text(x_vals(last_idx), 0.1*m*max(y_data), sprintf('  M%d', m), ...
                'Color', colors(m,:), 'FontSize', 9, ...
                'VerticalAlignment', 'middle');
        end
    end
    xlim([1 nSlices*1.1]);
    ylim([0 max(y_data)*1.1]);
    sgtitle(strrep(['Background_mask_diagnostics_trace:_',current_mouse_type,'_(',channel,')'],'_',' '))
    set(h_fig, 'InvertHardcopy', 'off'); % Preserve background color
    clean_fig_name = regexprep(h_fig.Name, '[^a-zA-Z0-9]', '_');
    save_path_base = fullfile(base_dir, clean_fig_name);
    fprintf('Saving diagnostic plot to: %s\n', save_path_base);
    saveas(h_fig, [save_path_base '.fig']);
    exportgraphics(h_fig, [save_path_base '.png'], 'Resolution', 300, 'BackgroundColor', 'current');

    %% --- Step 1: GLOBAL POOLING & NORMALIZATION CALCULATION ---

    fprintf('Pooling pixels from %d slices for global normalization...\n', length(slices_range_for_norm));

    cortex_samples_pooled = []; % Will become [TotalPixels x NumMice]
    num_mice_subset = numel(subset_indices);

    for s_idx = slices_range_for_norm
        % 1. Get Atlas Mask for this slice
        mask_2d_slice = squeeze(cortex_mask_3d_all(s_idx, :, :));
        valid_pixels_indices = find(mask_2d_slice == 1);

        if isempty(valid_pixels_indices)
            continue;
        end

        temp_samples = NaN(length(valid_pixels_indices), num_mice_subset);

        for iii = 1:num_mice_subset
            current_slice = squeeze(data_4d(s_idx, :, :, iii));
            current_indiv_mask = ~squeeze(recomputed_bkg_mask_4d(s_idx, :, :, iii)); % Foreground

            vals = current_slice(valid_pixels_indices);
            is_valid_tissue = current_indiv_mask(valid_pixels_indices);

            vals(~is_valid_tissue) = NaN;
            temp_samples(:, iii) = vals;
        end

        % Append to pooled matrix
        cortex_samples_pooled = [cortex_samples_pooled; temp_samples]; %#ok<AGROW>
    end

    % --- Calculate Consensus of POOLED Data ---
    fprintf('Calculating Global Median Consensus...\n');
    consensus_pixels_pooled = nanmedian(cortex_samples_pooled, 2);

    % --- Calculate Normalization Parameters (Slope/Intercept) on POOLED Data ---
    fprintf('Calculating Global Normalization Parameters...\n');
    norm_params = zeros(num_mice_subset, 2);

    for i = 1:num_mice_subset
        y_raw = cortex_samples_pooled(:, i);
        x_ref = consensus_pixels_pooled;

        % Robust Linear Regression
        valid_idx = ~isnan(x_ref) & ~isnan(y_raw);
        if sum(valid_idx) > 100 % Ensure enough points
            try
                mdl = fitlm(x_ref(valid_idx), y_raw(valid_idx), 'RobustOpts', 'on');
                p(1) = mdl.Coefficients.Estimate(2); % Slope
                p(2) = mdl.Coefficients.Estimate(1); % Intercept
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

    %% --- Step 1b: GLOBAL DIAGNOSTIC PLOT (All Pooled Slices) ---

    fprintf('Generating Global Diagnostic Plot (Pooled Data)...\n');

    % Performance Optimization: Subsample if too many points for plotting
    % (Visualizing millions of transparent dots is extremely slow)
    max_plot_points = 50000;
    total_pooled = size(cortex_samples_pooled, 1);
    if total_pooled > max_plot_points
        rng(42); % Constant seed for reproducibility
        idx_sub = randperm(total_pooled, max_plot_points);
    else
        idx_sub = 1:total_pooled;
    end

    figure('Visible', 'off', 'Name', 'Global Normalization Diagnostic (Pooled)', 'Color', 'w', 'Units', 'normalized', 'Position', [-0.05 -0.05 0.95 0.95]);

    for i = 1:num_mice_subset
        mouse_name = strrep(current_mice{i}, '_', ' ');

        % Get Subsampled Data
        y_raw = cortex_samples_pooled(idx_sub, i);
        x_ref = consensus_pixels_pooled(idx_sub);

        % Retrieve Global Slope/Intercept calculated previously
        p = norm_params(i, :);

        % Apply Normalization to the subset
        y_norm = (y_raw - p(2)) / p(1);

        % --- Row 1: Diagnostic (Raw vs Global Median) ---
        subplot(2, num_mice_subset, i);
        scatter(x_ref, y_raw, 2, 'k', 'filled', 'MarkerFaceAlpha', 0.1); hold on;
        axis equal

        % Plot Global Fit Line
        valid_idx = ~isnan(x_ref) & ~isnan(y_raw);
        if any(valid_idx)
            min_x = min(x_ref(valid_idx)); max_x = max(x_ref(valid_idx));
            x_grid = linspace(min_x, max_x, 100);
            y_fit = p(1)*x_grid + p(2);
            plot(x_grid, y_fit, 'b-', 'LineWidth', 2);
        end
        plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1.5);
        text(0.05*plot_limit, 0.85*plot_limit, sprintf('Slope: %.2f\nInt: %.0f', p(1), p(2)), ...
            'Color', 'b', 'FontSize', 9, 'FontWeight', 'bold');

        title(mouse_name, 'FontSize', 11, 'FontWeight', 'bold');
        xlim([0 plot_limit]); ylim([0 plot_limit]); grid on;
        set(gca, 'XTickLabel', []);
        if i == 1, ylabel({'Raw Intensity';'(Pooled Subset)'}, 'FontSize', 10); else, set(gca, 'YTickLabel', []); end

        % --- Row 2: Verification (Normalized vs Global Median) ---
        subplot(2, num_mice_subset, i + num_mice_subset);
        scatter(x_ref, y_norm, 2, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1); hold on;
        axis equal
        plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1.5);

        R = corrcoef(x_ref, y_norm, 'Rows', 'complete');
        if numel(R) > 1, r_val = R(1,2); else, r_val = NaN; end
        text(0.05*plot_limit, 0.9*plot_limit, sprintf('R = %.2f', r_val), ...
            'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

        title(mouse_name, 'FontSize', 11, 'FontWeight', 'bold');
        xlim([0 plot_limit]); ylim([0 plot_limit]); grid on;
        xlabel('Global Median Intensity');
        if i == 1, ylabel({'Normalized';'(Corrected)'}, 'FontSize', 10); else, set(gca, 'YTickLabel', []); end
    end

    sgtitle(['Global Normalization Diagnostic (Pooled Slices: ' num2str(min(slices_range_for_norm)) '-' num2str(max(slices_range_for_norm)) ')'], 'FontSize', 14);

    % Save
    if exist('global_diagnostics_dir', 'var'), save_output_dir = fullfile(global_diagnostics_dir, ['normalization_checks_' channel]); else, save_output_dir = fullfile(pwd, ['normalization_checks_' channel]); end
    if ~exist(save_output_dir, 'dir'), mkdir(save_output_dir); end
    fig_handle = gcf; set(fig_handle, 'InvertHardcopy', 'off');
    clean_fig_name = regexprep(fig_handle.Name, '[^a-zA-Z0-9]', '_');
    if isempty(clean_fig_name), clean_fig_name = 'Untitled_Figure'; end
    filename_base = clean_fig_name; % No slice number needed as this is global
    saveas(fig_handle, fullfile(save_output_dir, [filename_base '.fig']));
    exportgraphics(fig_handle, fullfile(save_output_dir, [filename_base '.png']), 'Resolution', 300, 'BackgroundColor', 'current');


    %% --- Step 2: VISUALIZATION LOOP (Per user-selected slice) ---

    fprintf('Starting visualization loop for %d slices...\n', length(slices_to_visualize_list));

    for viz_idx = 1:length(slices_to_visualize_list)

        slice_to_plot = slices_to_visualize_list(viz_idx);
        fprintf('  Visualizing Slice %d...\n', slice_to_plot);

        close all % Close figures from previous slice iteration

        % Get masks for this specific slice
        cortex_slice_mask = squeeze(cortex_mask_3d_all(slice_to_plot,:,:));
        mask_2d_slice = squeeze(cortex_mask_3d_all(slice_to_plot, :, :));
        valid_pixels_indices = find(mask_2d_slice == 1);

        %% --- Plot 1: Individual Slices with Cortex Highlight ---
        for iii = 1:num_mice_subset
            figure('Visible', 'off', 'Name', ['Individual_Slice_' current_mice{iii}], 'Color', 'k');

            mask_data = not(squeeze(recomputed_bkg_mask_4d(slice_to_plot, :, :, iii)));
            overlay_alpha = zeros(size(cortex_slice_mask));
            overlay_alpha(cortex_slice_mask .* mask_data == 0) = 0.75;

            img_data = squeeze(data_4d(slice_to_plot,:,:,iii));
            imagesc(img_data); colormap(hot); clim([0, plot_limit]); hold on;

            black_overlay = zeros(size(img_data));
            h_ov = imagesc(black_overlay); set(h_ov, 'AlphaData', overlay_alpha);

            axis image; axis off; set(gca, 'Color', 'k');
            mouse_name = strrep(current_mice{iii}, '_', ' ');
            t = title(mouse_name); set(t, 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');

            cb = colorbar; cb.Label.String = 'Intensity (a.u.)'; cb.Color = 'w'; cb.Label.Color = 'w';
            hold off;

            % Save
            if exist('global_diagnostics_dir', 'var'), save_output_dir = fullfile(global_diagnostics_dir, ['normalization_checks_' channel]); else, save_output_dir = fullfile(pwd, ['normalization_checks_' channel]); end
            if ~exist(save_output_dir, 'dir'), mkdir(save_output_dir); end
            fig_handle = gcf; set(fig_handle, 'InvertHardcopy', 'off');
            clean_fig_name = regexprep(fig_handle.Name, '[^a-zA-Z0-9]', '_');
            if isempty(clean_fig_name), clean_fig_name = 'Untitled_Figure'; end
            filename_base = sprintf('%s_Slice%d', clean_fig_name, slice_to_plot);
            saveas(fig_handle, fullfile(save_output_dir, [filename_base '.fig']));
            exportgraphics(fig_handle, fullfile(save_output_dir, [filename_base '.png']), 'Resolution', 300, 'BackgroundColor', 'current');
        end

        %% --- Extraction for Local Slice (For Scatter Plots) ---
        cortex_samples = NaN(length(valid_pixels_indices), num_mice_subset);

        for iii = 1:num_mice_subset
            current_slice = squeeze(data_4d(slice_to_plot, :, :, iii));
            current_indiv_mask = ~squeeze(recomputed_bkg_mask_4d(slice_to_plot, :, :, iii));
            vals = current_slice(valid_pixels_indices);
            is_valid_tissue = current_indiv_mask(valid_pixels_indices);
            vals(~is_valid_tissue) = NaN;
            cortex_samples(:, iii) = vals;
        end

        %% --- Plot 2: Pairwise Cortex Intensity Comparison (Mosaic) ---
        figure('Visible', 'off', 'Name', 'Pairwise Cortex Intensity Comparison', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

        for row = 1:num_mice_subset
            for col = 1:num_mice_subset
                idx = (row - 1) * num_mice_subset + col;
                subplot(num_mice_subset, num_mice_subset, idx);
                name_row = strrep(current_mice{row}, '_', ' ');
                name_col = strrep(current_mice{col}, '_', ' ');

                if row == col
                    histogram(cortex_samples(:, row), hist_num_bins, 'EdgeColor', 'none', 'FaceColor', 'r', 'BinLimits', [0, plot_limit]);
                    title(name_row, 'FontWeight', 'bold', 'FontSize', 8); grid on; xlim([0 plot_limit]); yticklabels([]);
                else
                    x_data = cortex_samples(:, col); y_data = cortex_samples(:, row);
                    scatter(x_data, y_data, 2, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1); hold on;
                    plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1);
                    R = corrcoef(x_data, y_data, 'Rows', 'complete');
                    if numel(R) > 1, r_val = R(1,2); else, r_val = NaN; end
                    text(0.05 * plot_limit, 0.9 * plot_limit, sprintf('R = %.2f', r_val), 'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold');
                    xlim([0 plot_limit]); ylim([0 plot_limit]); grid on;
                    if col == 1, ylabel(name_row, 'FontSize', 8, 'FontWeight', 'bold'); else, set(gca, 'YTickLabel', []); end
                    if row == num_mice_subset, xlabel(name_col, 'FontSize', 8, 'FontWeight', 'bold'); else, set(gca, 'XTickLabel', []); end
                end
            end
        end
        sgtitle(['Cortex Pixel Intensity Comparison (Slice ' num2str(slice_to_plot) ') - Range [0, ' num2str(plot_limit) ']']);

        % Save
        if exist('global_diagnostics_dir', 'var'), save_output_dir = fullfile(global_diagnostics_dir, ['normalization_checks_' channel]); else, save_output_dir = fullfile(pwd, ['normalization_checks_' channel]); end
        if ~exist(save_output_dir, 'dir'), mkdir(save_output_dir); end
        fig_handle = gcf; set(fig_handle, 'InvertHardcopy', 'off');
        clean_fig_name = regexprep(fig_handle.Name, '[^a-zA-Z0-9]', '_');
        if isempty(clean_fig_name), clean_fig_name = 'Untitled_Figure'; end
        filename_base = sprintf('%s_Slice%d', clean_fig_name, slice_to_plot);
        saveas(fig_handle, fullfile(save_output_dir, [filename_base '.fig']));
        exportgraphics(fig_handle, fullfile(save_output_dir, [filename_base '.png']), 'Resolution', 300, 'BackgroundColor', 'current');

        %% --- Plot 3: Median Consensus Slice ---
        % Local consensus for this slice
        slice_stack = squeeze(data_4d(slice_to_plot, :, :, 1:num_mice_subset));
        consensus_slice = median(slice_stack, 3);
        consensus_pixels = nanmedian(cortex_samples, 2); %#ok<NANMEDIAN>

        figure('Visible', 'off', 'Name', 'Median Consensus Slice');
        imagesc(consensus_slice); colormap(hot); clim([0, plot_limit]); hold on;

        slice_bkg_stack = squeeze(recomputed_bkg_mask_4d(slice_to_plot, :, :, 1:num_mice_subset));
        at_least_one_tissue = any(~slice_bkg_stack, 3);
        median_overlay_alpha = zeros(size(consensus_slice));
        is_valid_region = (mask_2d_slice == 1) & at_least_one_tissue;
        median_overlay_alpha(~is_valid_region) = 0.5;

        black_overlay = zeros(size(consensus_slice));
        h_ov = imagesc(black_overlay); set(h_ov, 'AlphaData', median_overlay_alpha);
        axis image; axis off; set(gca, 'Color', 'k');
        t = title(['Median Mouse (Slice ' num2str(slice_to_plot) ')']); set(t, 'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
        cb = colorbar; cb.Label.String = 'Intensity (a.u.)'; cb.Color = 'k'; cb.Label.Color = 'k';
        hold off;

        % Save
        if exist('global_diagnostics_dir', 'var'), save_output_dir = fullfile(global_diagnostics_dir, ['normalization_checks_' channel]); else, save_output_dir = fullfile(pwd, ['normalization_checks_' channel]); end
        if ~exist(save_output_dir, 'dir'), mkdir(save_output_dir); end
        fig_handle = gcf; set(fig_handle, 'InvertHardcopy', 'off');
        clean_fig_name = regexprep(fig_handle.Name, '[^a-zA-Z0-9]', '_');
        if isempty(clean_fig_name), clean_fig_name = 'Untitled_Figure'; end
        filename_base = sprintf('%s_Slice%d', clean_fig_name, slice_to_plot);
        saveas(fig_handle, fullfile(save_output_dir, [filename_base '.fig']));
        exportgraphics(fig_handle, fullfile(save_output_dir, [filename_base '.png']), 'Resolution', 300, 'BackgroundColor', 'current');

        %% --- Plot 4: Individual vs Median Scatter ---
        figure('Visible', 'off', 'Name', 'Individual vs Median Comparison', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.2 0.8 0.4]);

        for i = 1:num_mice_subset
            subplot(1, num_mice_subset, i);
            mouse_name = strrep(current_mice{i}, '_', ' ');

            x_data = consensus_pixels;
            y_data = cortex_samples(:, i);

            scatter(x_data, y_data, 2, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1); hold on;
            plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1);
            axis equal
            R = corrcoef(x_data, y_data, 'Rows', 'complete');
            if numel(R) > 1, r_val = R(1,2); else, r_val = NaN; end
            text(0.05 * plot_limit, 0.9 * plot_limit, sprintf('R = %.2f', r_val), 'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
            xlim([0 plot_limit]); ylim([0 plot_limit]); grid on;
            title(mouse_name, 'FontSize', 10, 'FontWeight', 'bold');
            xlabel('Median Intensity');
            if i == 1, ylabel('Individual Intensity'); else, set(gca, 'YTickLabel', []); end
        end
        sgtitle(['Individual Mice vs. Group Median (Slice ' num2str(slice_to_plot) ')']);

        % Save
        if exist('global_diagnostics_dir', 'var'), save_output_dir = fullfile(global_diagnostics_dir, ['normalization_checks_' channel]); else, save_output_dir = fullfile(pwd, ['normalization_checks_' channel]); end
        if ~exist(save_output_dir, 'dir'), mkdir(save_output_dir); end
        fig_handle = gcf; set(fig_handle, 'InvertHardcopy', 'off');
        clean_fig_name = regexprep(fig_handle.Name, '[^a-zA-Z0-9]', '_');
        if isempty(clean_fig_name), clean_fig_name = 'Untitled_Figure'; end
        filename_base = sprintf('%s_Slice%d', clean_fig_name, slice_to_plot);
        saveas(fig_handle, fullfile(save_output_dir, [filename_base '.fig']));
        exportgraphics(fig_handle, fullfile(save_output_dir, [filename_base '.png']), 'Resolution', 300, 'BackgroundColor', 'current');

        %% --- Plot 5: Normalization Diagnostic (Using Global Norm Params on Local Slice) ---

        cortex_samples_norm = zeros(size(cortex_samples));

        figure('Visible', 'off', 'Name', 'Normalization Diagnostic and Verification', 'Color', 'w', 'Units', 'normalized', 'Position', [-0.05 -0.05 0.95 0.95]);

        for i = 1:num_mice_subset
            mouse_name = strrep(current_mice{i}, '_', ' ');
            y_raw = cortex_samples(:, i);
            x_ref = consensus_pixels;

            % Retrieve Global Slope/Intercept
            p = norm_params(i, :);

            % Apply Normalization
            cortex_samples_norm(:, i) = (y_raw - p(2)) / p(1);

            % Plot Row 1: Diagnostic
            subplot(2, num_mice_subset, i);
            scatter(x_ref, y_raw, 2, 'k', 'filled', 'MarkerFaceAlpha', 0.1); hold on;
            axis equal

            valid_idx = ~isnan(x_ref) & ~isnan(y_raw);
            if any(valid_idx)
                min_x = min(x_ref(valid_idx)); max_x = max(x_ref(valid_idx));
                x_grid = linspace(min_x, max_x, 100);
                y_fit = p(1)*x_grid + p(2);
                plot(x_grid, y_fit, 'b-', 'LineWidth', 2);
            end
            plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1.5);
            text(0.05*plot_limit, 0.85*plot_limit, sprintf('Slope: %.2f\nInt: %.0f', p(1), p(2)), 'Color', 'b', 'FontSize', 9, 'FontWeight', 'bold');

            title(mouse_name, 'FontSize', 11, 'FontWeight', 'bold');
            xlim([0 plot_limit]); ylim([0 plot_limit]); grid on;
            set(gca, 'XTickLabel', []);
            if i == 1, ylabel({'Raw Intensity';'(Individual)'}, 'FontSize', 10); else, set(gca, 'YTickLabel', []); end

            % Plot Row 2: Verification
            subplot(2, num_mice_subset, i + num_mice_subset);
            scatter(x_ref, cortex_samples_norm(:, i), 2, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1); hold on;
            axis equal
            plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1.5);

            R = corrcoef(x_ref, cortex_samples_norm(:, i), 'Rows', 'complete');
            if numel(R) > 1, r_val = R(1,2); else, r_val = NaN; end
            text(0.05*plot_limit, 0.9*plot_limit, sprintf('R = %.2f', r_val), 'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

            title(mouse_name, 'FontSize', 11, 'FontWeight', 'bold');
            xlim([0 plot_limit]); ylim([0 plot_limit]); grid on;
            xlabel('Median Intensity (Target)');
            if i == 1, ylabel({'Normalized';'(Corrected)'}, 'FontSize', 10); else, set(gca, 'YTickLabel', []); end
        end
        sgtitle(['Normalization diagnostic (Slice ' num2str(slice_to_plot) ')'], 'FontSize', 14);

        % Save
        if exist('global_diagnostics_dir', 'var'), save_output_dir = fullfile(global_diagnostics_dir, ['normalization_checks_' channel]); else, save_output_dir = fullfile(pwd, ['normalization_checks_' channel]); end
        if ~exist(save_output_dir, 'dir'), mkdir(save_output_dir); end
        fig_handle = gcf; set(fig_handle, 'InvertHardcopy', 'off');
        clean_fig_name = regexprep(fig_handle.Name, '[^a-zA-Z0-9]', '_');
        if isempty(clean_fig_name), clean_fig_name = 'Untitled_Figure'; end
        filename_base = sprintf('%s_Slice%d', clean_fig_name, slice_to_plot);
        saveas(fig_handle, fullfile(save_output_dir, [filename_base '.fig']));
        exportgraphics(fig_handle, fullfile(save_output_dir, [filename_base '.png']), 'Resolution', 300, 'BackgroundColor', 'current');
        close all

        %% --- Plot 6: Visual Verification (Normalized vs Original Stacked) ---
        figure('Visible', 'off', 'Name', 'Visual Verification Normalized vs Raw', 'Color', 'k', 'Units', 'normalized', 'Position', [0.1 0.05 0.8 0.8]);

        for i = 1:num_mice_subset
            mouse_name = strrep(current_mice{i}, '_', ' ');

            img_raw = squeeze(data_4d(slice_to_plot,:,:,i));
            current_bg_mask = squeeze(recomputed_bkg_mask_4d(slice_to_plot,:,:,i));
            is_valid_tissue = (mask_2d_slice == 1) & (current_bg_mask == 0);

            overlay_alpha = zeros(size(img_raw));
            overlay_alpha(~is_valid_tissue) = 0.5;

            slope = norm_params(i, 1);
            intercept = norm_params(i, 2);
            img_norm = (img_raw - intercept) / slope;

            subplot(2, num_mice_subset, i);
            imagesc(img_norm); colormap(hot); clim([0, plot_limit]);
            axis image; axis off; set(gca, 'Color', 'k'); hold on;
            h_ov1 = imagesc(zeros(size(img_norm))); set(h_ov1, 'AlphaData', overlay_alpha);
            t = title(['Norm: ' mouse_name]); set(t, 'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold');
            if i == num_mice_subset
                cb = colorbar; cb.Label.String = 'Normalized Intensity'; cb.Color = 'w'; cb.Label.Color = 'w';
                cb.Position = [0.92 0.55 0.01 0.35];
            end
            hold off;

            subplot(2, num_mice_subset, i + num_mice_subset);
            imagesc(img_raw); colormap(hot); clim([0, plot_limit]);
            axis image; axis off; set(gca, 'Color', 'k'); hold on;
            h_ov2 = imagesc(zeros(size(img_raw))); set(h_ov2, 'AlphaData', overlay_alpha);
            t = title(['Raw: ' mouse_name]); set(t, 'Color', 'w', 'FontSize', 11, 'FontWeight', 'normal');
            if i == num_mice_subset
                cb = colorbar; cb.Label.String = 'Raw Intensity'; cb.Color = 'w'; cb.Label.Color = 'w';
                cb.Position = [0.92 0.1 0.01 0.35];
            end
            hold off;
        end
        sgtitle(['Comparison: normalized (top) vs. raw (bottom) - Fixed scale [0, ' num2str(plot_limit) ']'], 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');

        % Save
        if exist('global_diagnostics_dir', 'var'), save_output_dir = fullfile(global_diagnostics_dir, ['normalization_checks_' channel]); else, save_output_dir = fullfile(pwd, ['normalization_checks_' channel]); end
        if ~exist(save_output_dir, 'dir'), mkdir(save_output_dir); end
        fig_handle = gcf; set(fig_handle, 'InvertHardcopy', 'off');
        clean_fig_name = regexprep(fig_handle.Name, '[^a-zA-Z0-9]', '_');
        if isempty(clean_fig_name), clean_fig_name = 'Untitled_Figure'; end
        filename_base = sprintf('%s_Slice%d', clean_fig_name, slice_to_plot);
        saveas(fig_handle, fullfile(save_output_dir, [filename_base '.fig']));
        exportgraphics(fig_handle, fullfile(save_output_dir, [filename_base '.png']), 'Resolution', 300, 'BackgroundColor', 'current');
        close all

        fprintf('  Done with Slice %d\n', slice_to_plot);
    end

    %% --- Step 3: Generate Verification Video (Full Volume) ---

    if plot_verification_video

        fprintf('Generating Normalization Verification Video (this may take a while)...\n');

        % 1. Setup Video Writer
        video_filename = ['Normalization_Verification_' current_mouse_type '_' channel '.mp4'];

        if exist('global_diagnostics_dir', 'var')
            save_video_path = fullfile(global_diagnostics_dir, ['normalization_checks_' channel]);
        else
            save_video_path = fullfile(pwd, ['normalization_checks_' channel]);
        end
        if ~exist(save_video_path, 'dir'), mkdir(save_video_path); end

        full_video_path = fullfile(save_video_path, video_filename);

        vidObj = VideoWriter(full_video_path, 'MPEG-4');
        vidObj.FrameRate = 15; % Adjust playback speed (frames per second)
        vidObj.Quality = 95;   % High quality
        open(vidObj);

        % 2. Define range: Iterate over ALL slices available in the volume
        slices_to_video = 1:size(data_4d, 1);

        for s_idx = slices_to_video

            % Optimization: Skip slices where the Atlas Mask is empty
            % (No point plotting empty black screens)
            mask_2d_slice = squeeze(cortex_mask_3d_all(s_idx, :, :));
            if sum(mask_2d_slice(:)) == 0
                if mod(s_idx, 100) == 0, fprintf('  Skipping empty slice %d\n', s_idx); end
                continue;
            end

            % 3. Create invisible figure with black background
            % 'visible', 'off' is crucial for speed so it doesn't pop up on screen
            fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'k');
            set(fh, 'InvertHardcopy', 'off'); % Ensure background stays black

            for i = 1:num_mice_subset
                mouse_name = strrep(current_mice{i}, '_', ' ');

                % --- Prepare Data ---
                img_raw = squeeze(data_4d(s_idx,:,:,i));
                current_bg_mask = squeeze(recomputed_bkg_mask_4d(s_idx,:,:,i));

                % --- Overlay Logic ---
                % Dim if: NOT (In Atlas AND In Foreground)
                is_valid_tissue = (mask_2d_slice == 1) & (current_bg_mask == 0);
                overlay_alpha = zeros(size(img_raw));
                overlay_alpha(~is_valid_tissue) = 0.5;

                % --- Apply Normalization ---
                slope = norm_params(i, 1);
                intercept = norm_params(i, 2);
                img_norm = (img_raw - intercept) / slope;

                % --- Plot Top Row (Normalized) ---
                subplot(2, num_mice_subset, i);
                imagesc(img_norm); colormap(hot); clim([0, plot_limit]);
                axis image; axis off; set(gca, 'Color', 'k'); hold on;

                h_ov1 = imagesc(zeros(size(img_norm)));
                set(h_ov1, 'AlphaData', overlay_alpha);

                t = title(['Norm: ' mouse_name]);
                set(t, 'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');

                if i == num_mice_subset
                    cb = colorbar; cb.Label.String = 'Norm Int'; cb.Color = 'w'; cb.Label.Color = 'w';
                    cb.Position = [0.92 0.55 0.01 0.35];
                end
                hold off;

                % --- Plot Bottom Row (Raw) ---
                subplot(2, num_mice_subset, i + num_mice_subset);
                imagesc(img_raw); colormap(hot); clim([0, plot_limit]);
                axis image; axis off; set(gca, 'Color', 'k'); hold on;

                h_ov2 = imagesc(zeros(size(img_raw)));
                set(h_ov2, 'AlphaData', overlay_alpha);

                t = title(['Raw: ' mouse_name]);
                set(t, 'Color', 'w', 'FontSize', 10, 'FontWeight', 'normal');

                if i == num_mice_subset
                    cb = colorbar; cb.Label.String = 'Raw Int'; cb.Color = 'w'; cb.Label.Color = 'w';
                    cb.Position = [0.92 0.1 0.01 0.35];
                end
                hold off;
            end

            % Main Title
            sgtitle(['Slice ' num2str(s_idx) ' - Norm vs Raw - Scale [0 ' num2str(plot_limit) ']'], ...
                'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');

            % 4. Write Frame
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

    %% --- Step 4: Apply Normalization to Full Volume and Save ---

    fprintf('Applying normalization to full 4D volume and saving...\n');

    % 1. Initialize output variable (same size/type as input to save memory, or single for precision)
    data_4d_normalized = zeros(size(data_4d), 'single');

    % 2. Apply normalization per mouse
    for i = 1:num_mice_subset
        slope = norm_params(i, 1);
        intercept = norm_params(i, 2);

        % Apply formula: (Raw - Intercept) / Slope
        data_4d_normalized(:,:,:,i) = (data_4d(:,:,:,i) - intercept) / slope;

        fprintf('  Applied normalization to Mouse %d/%d (Slope=%.2f, Int=%.2f)\n', ...
            i, num_mice_subset, slope, intercept);
    end

    % 3. Save with channel-aware filenames + variable names (e.g. nano_4d_normalized
    % or auto_4d_normalized) so downstream scripts can pull the correct field via
    % dynamic access. The bkgmask file keeps a generic field name (channel-agnostic
    % across consumers).
    norm_var_name     = [channel '_4d_normalized'];
    save_filename     = fullfile(base_dir, [channel '_4d_normalized' S.tag '.mat']);
    save_filename_bis = fullfile(base_dir, [channel '_4d_normalized_bkgmask' S.tag '.mat']);

    save_struct                        = struct();
    save_struct.(norm_var_name)        = data_4d_normalized;
    save_struct.norm_params            = norm_params;
    save_struct.current_mice           = current_mice;
    save_struct.selected_mice_idx_list = selected_mice_idx_list;
    save_struct.channel                = channel;
    save_struct.cohort_spec            = S.label;
    save_struct.atlas_key              = S.atlas_key;
    save(save_filename, '-struct', 'save_struct', '-v7.3');

    save(save_filename_bis, 'recomputed_bkg_mask_4d', 'current_mice', 'selected_mice_idx_list', '-v7.3');

    fprintf('Successfully saved normalized volume to:\n  %s\n', save_filename);

    % Clear large variable to free memory before next group iteration
    clear data_4d_normalized recomputed_bkg_mask_4d save_struct

    toc

end