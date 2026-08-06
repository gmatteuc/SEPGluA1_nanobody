
% --- Parameters ---
smooth_sigma = 8;       % Controls smoothness (Higher = smoother but loses detail)
expand_threshold = 0.04; % Controls expansion (Lower = mask gets bigger. 0.5 = same size)

% Create a copy to store result (or overwrite if RAM is tight)
smoothed_mask_4d = false(size(mask_4d)); 

[rows, cols, slices, num_mice] = size(mask_4d);

fprintf('Processing mask smoothing...\n');

% Loop over mice
for m = 1:num_mice
    fprintf('  Processing Mouse %d / %d ...\n', m, num_mice);
    
    % Loop over slices (Parallel loop 'parfor' can be used here if available)
    for z = 1:slices
        
        % 1. Extract Slice
        current_mask = mask_4d(:, :, z, m);
        
        % Skip empty slices to save time
        if ~any(current_mask(:))
            continue; 
        end
        
        % 2. Convert to double for filtering
        mask_dbl = double(current_mask);
        
        % 3. Apply Gaussian Blur (2D is very fast)
        % This smooths the jagged edges
        blurred_mask = imgaussfilt(mask_dbl, smooth_sigma);
        
        % 4. Threshold to re-binarize
        % Threshold < 0.5 expands the mask outward
        % Threshold > 0.5 shrinks the mask inward
        new_mask = blurred_mask > expand_threshold;
        
        % 5. Store result
        smoothed_mask_4d(:, :, z, m) = new_mask;
    end
end
fprintf('Done.\n');


% ----------- Visualization smoothnss -----------

slice_check = 500; % Choose your slice
mouse_check = 1;

% Original
mask_orig = squeeze(mask_4d(slice_check, :, :, mouse_check));
img_data = squeeze(nano_4d(slice_check, :, :, mouse_check));

% Smoothed
mask_new = squeeze(smoothed_mask_4d(slice_check, :, :, mouse_check));

figure('Name', 'Mask Smoothing Comparison', 'Color', 'w', 'Position', [100 100 1200 500]);

% Plot 1: Original jagged mask
subplot(1,2,1);
B1 = imoverlay(img_data./quantile(img_data(:),0.9), mask_orig, [1 0 1]); % Pink
imshow(B1);
title('Original (Jagged)');
xlim([100 400]); ylim([100 400]); % Zoom in to a corner to see pixels

% Plot 2: Smoothed & Expanded mask
subplot(1,2,2);
B2 = imoverlay(img_data./quantile(img_data(:),0.9), mask_new, [0 1 0]); % Green
imshow(B2);
title(sprintf('Smoothed (Sigma: %d, Thresh: %.1f)', smooth_sigma, expand_threshold));
xlim([100 400]); ylim([100 400]); % Zoom in to match



















    %% --- Consensus (Median) Analysis ---

    % 1. Compute the Median Slice and Extract Pixels
    fprintf('Computing Median Slice across %d mice...\n', num_mice_subset);

    % Extract the stack for the current slice: [Height x Width x NumMice]
    slice_stack = squeeze(nano_4d(slice_to_plot, :, :, 1:num_mice_subset));

    % Compute Median (change to 'mean' if preferred)
    consensus_slice = median(slice_stack, 3);

    % Extract the specific masked pixels for the scatter plot
    consensus_pixels = consensus_slice(valid_pixels_indices);

    % 2. Visualize the Consensus Slice (with Dimming Mask)
    figure('Name', 'Median Consensus Slice');

    imagesc(consensus_slice);
    colormap(hot);
    clim([0, plot_limit]);
    hold on;

    % Apply the dimming overlay (reusing overlay_alpha from previous section)
    black_overlay = zeros(size(consensus_slice));
    h_ov = imagesc(black_overlay);
    set(h_ov, 'AlphaData', overlay_alpha);

    axis image; axis off;
    set(gca, 'Color', 'k');

    t = title(['Median Mouse (Slice ' num2str(slice_to_plot) ')']);
    set(t, 'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');

    cb = colorbar;
    cb.Label.String = 'Intensity (a.u.)';
    cb.Color = 'k';
    cb.Label.Color = 'k';
    hold off;

    % 3. Plot Scatter: Individual Mouse vs Median

    figure('Name', 'Individual vs Median Comparison', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.2 0.8 0.4]);

    for i = 1:num_mice_subset
        subplot(1, num_mice_subset, i);

        mouse_name = strrep(current_mice{i}, '_', ' ');

        % Data for plotting
        x_data = consensus_pixels;       % Reference (Median)
        y_data = cortex_samples(:, i);   % Individual

        % Scatter plot
        scatter(x_data, y_data, 2, 'filled', 'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);

        hold on;

        % Identity Line
        plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1);
        axis equal

        % Correlation
        R = corrcoef(x_data, y_data);
        r_val = R(1,2);

        text(0.05 * plot_limit, 0.9 * plot_limit, ...
            sprintf('R = %.2f', r_val), 'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

        % Formatting
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



%% --- Normalization & Verification (Combined Figure) ---

% Initialize storage
cortex_samples_norm = zeros(size(cortex_samples));
norm_params = zeros(num_mice_subset, 2); 

fprintf('Computing normalization and generating combined plot...\n');

% Create single figure for both steps (2 Rows x N Columns)
figure('Name', 'Normalization: Diagnostic & Verification', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

for i = 1:num_mice_subset
    
    mouse_name = strrep(current_mice{i}, '_', ' ');
    y_raw = cortex_samples(:, i);      % The individual mouse
    x_ref = consensus_pixels;          % The target (Median)
    
    % --- 1. Calculate Fit (Robust Linear Regression) ---
    try
        % Fit line: y = Slope*x + Intercept
        mdl = fitlm(x_ref, y_raw, 'RobustOpts', 'on'); 
        p(1) = mdl.Coefficients.Estimate(2); % Slope
        p(2) = mdl.Coefficients.Estimate(1); % Intercept
    catch
        p = polyfit(x_ref, y_raw, 1); % Fallback
    end
    norm_params(i, :) = p;
    
    % --- 2. Apply Normalization ---
    % Formula: Normalized = (Raw - Intercept) / Slope
    cortex_samples_norm(:, i) = (y_raw - p(2)) / p(1);
    
    % --- 3. Plot Row 1: Diagnostic (Raw vs Median) ---
    subplot(2, num_mice_subset, i);
    
    % Scatter Raw Data
    scatter(x_ref, y_raw, 2, 'k', 'filled', 'MarkerFaceAlpha', 0.1); hold on;
    axis equal

    % Plot Fitted Line (Blue) - Limited to data range
    x_grid = linspace(min(x_ref), max(x_ref), 100);
    y_fit = p(1)*x_grid + p(2);
    plot(x_grid, y_fit, 'b-', 'LineWidth', 2);
    
    % Plot Target Line (Red Dashed) - Shows goal across full range
    plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1.5);
    
    % Stats Text
    text(0.05*plot_limit, 0.85*plot_limit, sprintf('Slope: %.2f\nInt: %.0f', p(1), p(2)), ...
        'Color', 'b', 'FontSize', 9, 'FontWeight', 'bold');
    
    % Formatting
    title(mouse_name, 'FontSize', 11, 'FontWeight', 'bold');
    xlim([0 plot_limit]); ylim([0 plot_limit]); grid on;
    set(gca, 'XTickLabel', []); % Remove X labels for top row
    
    if i == 1
        ylabel({'Raw Intensity';'(Individual)'}, 'FontSize', 10);
    else
        set(gca, 'YTickLabel', []);
    end
    
    % --- 4. Plot Row 2: Verification (Normalized vs Median) ---
    subplot(2, num_mice_subset, i + num_mice_subset);
    
    % Scatter Normalized Data (Blue)
    scatter(x_ref, cortex_samples_norm(:, i), 2, 'filled', 'MarkerFaceColor', 'b', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1); hold on;
    axis equal

    % Plot Identity Line
    plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1.5);
    
    % Calculate Correlation on Normalized Data
    R = corrcoef(x_ref, cortex_samples_norm(:, i));
    
    % Stats Text
    text(0.05*plot_limit, 0.9*plot_limit, sprintf('R = %.2f', R(1,2)), ...
        'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
        
    % Formatting
    title(mouse_name, 'FontSize', 11, 'FontWeight', 'bold');
    xlim([0 plot_limit]); ylim([0 plot_limit]); grid on;
    xlabel('Median Intensity (Target)');
    
    if i == 1
        ylabel({'Normalized';'(Corrected)'}, 'FontSize', 10);
    else
        set(gca, 'YTickLabel', []);
    end
end

sgtitle(['Normalization diagnostic (Slice ' num2str(slice_to_plot) ')'], 'FontSize', 14);


%% --- Visual Verification: Normalized vs Original (Stacked) ---

figure('Name', 'Visual Verification: Normalized vs Raw', 'Color', 'k', 'Units', 'normalized', 'Position', [0.1 0.05 0.8 0.8]);

% Re-create mask overlay variables (reusing from previous steps)
overlay_alpha = zeros(size(cortex_slice_mask));
overlay_alpha(cortex_slice_mask == 0) = 0.5; % Dimming factor

for i = 1:num_mice_subset

    mouse_name = strrep(current_mice{i}, '_', ' ');

    % 1. Prepare Data
    img_raw = squeeze(nano_4d(slice_to_plot,:,:,i));

    % Retrieve calc parameters
    slope = norm_params(i, 1);
    intercept = norm_params(i, 2);

    % Calculate Normalized Image
    img_norm = (img_raw - intercept) / slope;

    % --- ROW 1: NORMALIZED (The Result) ---
    subplot(2, num_mice_subset, i);

    imagesc(img_norm);
    colormap(hot);
    clim([0, plot_limit]); % Fixed scale
    axis image; axis off;
    set(gca, 'Color', 'k');
    hold on;

    % Dimming Overlay
    h_ov1 = imagesc(zeros(size(img_norm)));
    set(h_ov1, 'AlphaData', overlay_alpha);

    % Titles
    t = title(['Norm: ' mouse_name]);
    set(t, 'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold');

    % Colorbar (Rightmost only)
    if i == num_mice_subset
        cb = colorbar;
        cb.Label.String = 'Normalized Intensity';
        cb.Color = 'w'; cb.Label.Color = 'w';
        cb.Position = [0.92 0.55 0.01 0.35]; % Manual positioning for top row
    end
    hold off;

    % --- ROW 2: ORIGINAL (The Reference) ---
    subplot(2, num_mice_subset, i + num_mice_subset);

    imagesc(img_raw);
    colormap(hot);
    clim([0, plot_limit]); % SAME Fixed scale for comparison
    axis image; axis off;
    set(gca, 'Color', 'k');
    hold on;

    % Dimming Overlay
    h_ov2 = imagesc(zeros(size(img_raw)));
    set(h_ov2, 'AlphaData', overlay_alpha);

    % Titles
    t = title(['Raw: ' mouse_name]);
    set(t, 'Color', 'w', 'FontSize', 11, 'FontWeight', 'normal'); % Lighter font for raw

    % Colorbar (Rightmost only)
    if i == num_mice_subset
        cb = colorbar;
        cb.Label.String = 'Raw Intensity';
        cb.Color = 'w'; cb.Label.Color = 'w';
        cb.Position = [0.92 0.1 0.01 0.35]; % Manual positioning for bottom row
    end
    hold off;

end

sgtitle(['Comparison: normalized (top) vs. raw (bottom) - Fixed scale [0, ' num2str(plot_limit) ']'], ...
    'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');





%


% 3. Plot Scatter: Individual Mouse vs Median

figure('Name', 'Individual vs Median Comparison', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.2 0.8 0.4]);

for i = 1:num_mice_subset
    subplot(1, num_mice_subset, i);

    mouse_name = strrep(current_mice{i}, '_', ' ');

    % Data for plotting
    x_data = consensus_pixels;       % Reference (Median of valid pixels)
    y_data = cortex_samples(:, i);   % Individual (contains NaNs)

    % Scatter plot (Scatter ignores NaNs automatically)
    scatter(x_data, y_data, 2, 'filled', 'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);

    hold on;

    % Identity Line
    plot([0 plot_limit], [0 plot_limit], 'r--', 'LineWidth', 1);
    axis equal

    % Correlation (Handle NaNs)
    R = corrcoef(x_data, y_data, 'Rows', 'complete');
    if numel(R) > 1
        r_val = R(1,2);
    else
        r_val = NaN;
    end

    text(0.05 * plot_limit, 0.9 * plot_limit, ...
        sprintf('R = %.2f', r_val), 'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

    % Formatting
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





% 1. Define Target Regions
roi_list = { ...
    'Somatomotor areas', ...
    'Somatosensory areas', ...
    'Gustatory areas', ...
    'Visceral area', ...
    'Auditory areas', ...
    'Visual areas', ...
    'Anterior cingulate area', ...
    'Prelimbic area', ...
    'Infralimbic area', ...
    'Orbital area', ...
    'Agranular insular area', ...
    'Retrosplenial area', ...
    'Posterior parietal association areas', ...
    'Temporal association areas', ...
    'Perirhinal area', ...
    'Ectorhinal area', ...
    'Claustrum', ...
    'Endopiriform nucleus', ...
    'Lateral amygdalar nucleus', ...
    'Basolateral amygdalar nucleus', ...
    'Basomedial amygdalar nucleus', ...
    'Posterior amygdalar nucleus', ...
    'Olfactory areas', ...
    'Hippocampal formation', ...
    'Striatum dorsal region', ...
    'Striatum ventral region', ...
    'Lateral septal complex', ...
    'Striatum-like amygdalar nuclei', ...
    'Pallidum, dorsal region', ...
    'Pallidum, ventral region', ...
    'Pallidum, medial region', ...
    'Pallidum, caudal region', ...
    'Ventral anterior-lateral complex of the thalamus', ...
    'Ventral medial nucleus of the thalamus', ...
    'Ventral posterolateral nucleus of the thalamus', ...
    'Ventral posterolateral nucleus of the thalamus, parvicellular part', ...
    'Ventral posteromedial nucleus of the thalamus', ...
    'Ventral posteromedial nucleus of the thalamus, parvicellular part', ...
    'Subparafascicular nucleus', ...
    'Subparafascicular area', ...
    'Peripeduncular nucleus', ...
    'Geniculate group, dorsal thalamus', ...
    'Lateral posterior nucleus of the thalamus', ...
    'Posterior complex of the thalamus', ...
    'Posterior limiting nucleus of the thalamus', ...
    'Suprageniculate nucleus', ...
    'Anterior group of the dorsal thalamus', ...
    'Medial group of the dorsal thalamus', ...
    'Midline group of the dorsal thalamus', ...
    'Intralaminar nuclei of the dorsal thalamus', ...
    'Reticular nucleus of the thalamus', ...
    'Geniculate group, ventral thalamus', ...
    'Epithalamus', ...
    'Periventricular zone', ...
    'Periventricular region', ...
    'Hypothalamic medial zone', ...
    'Hypothalamic lateral zone', ...
    'Median eminence', ...
    'Midbrain, sensory related', ...
    'Midbrain, motor related', ...
    'Midbrain, behavioral state related', ...
    'Pons', ...
    'Medulla', ...
    'Cerebellar cortex', ...
    'Cerebellar nuclei' ...
    };



% %% Whole-Brain Ontology Analysis (Fast Vectorized Approach)
% 
% fprintf('Starting Whole-Brain Leaf Region Analysis...\n');
% 
% % 1. Setup Atlas & ID Mapping
% % We work on the Left Hemisphere to match lr_diff data
% atlas_left = AllenCrop(:, :, 1:half_width);
% 
% % Find valid pixels (Atlas > 0)
% valid_mask = atlas_left > 0;
% all_ids_vec = atlas_left(valid_mask); % Vector of IDs for every valid pixel
% 
% % Map large Allen IDs (e.g. 32948) to indices 1..N
% [unique_ids, ~, id_indices] = unique(all_ids_vec);
% n_unique_regions = length(unique_ids);
% 
% fprintf('Found %d unique leaf regions in the atlas volume.\n', n_unique_regions);
% 
% % 2. Pre-slice Data & Masks (Left Hemi)
% target_vol_ctrl_left = abs(lr_diff_ctrl); % Already left hemi? Check dimensions.
% % Note: lr_diff_ctrl is typically [X, Y, Z, Mice]. 
% % If lr_diff was already computed as Left-Right, its Z is half_width.
% % If it matches atlas_left, we are good.
% 
% bg_mask_ctrl_left = logical(recomputed_bkg_mask_4d_ctrl(:, :, 1:half_width, :));
% bg_mask_exp_left  = logical(recomputed_bkg_mask_4d_exp(:, :, 1:half_width, :));
% 
% % 3. Accumulate Means per Region (Vectorized)
% % Preallocate [Regions x Mice]
% leaf_means_ctrl = nan(n_unique_regions, n_ctrl);
% leaf_means_exp  = nan(n_unique_regions, n_exp);
% 
% fprintf('Computing means for all regions simultaneously...\n');
% 
% % --- Process Control Mice ---
% for m = 1:n_ctrl
%     % Get Data & Mask for Mouse m
%     vol_data = target_vol_ctrl_left(:,:,:,m);
%     vol_bg   = bg_mask_ctrl_left(:,:,:,m);
% 
%     % Flatten to match atlas vector
%     data_vec = vol_data(valid_mask);
%     bg_vec   = vol_bg(valid_mask);
% 
%     % Set Background to NaN
%     data_vec(bg_vec) = NaN;
% 
%     % Accumulate: Compute Mean per ID index
%     % accumarray(indices, values, [output_size], function)
%     sums   = accumarray(id_indices, data_vec, [n_unique_regions 1], @nansum);
%     counts = accumarray(id_indices, ~isnan(data_vec), [n_unique_regions 1], @sum);
% 
%     % Avoid div by zero
%     leaf_means_ctrl(:, m) = sums ./ counts;
% end
% 
% % --- Process Exp Mice ---
% for m = 1:n_exp
%     vol_data = target_vol_exp(:,:,:,m);
%     vol_bg   = bg_mask_exp_left(:,:,:,m);
% 
%     data_vec = vol_data(valid_mask);
%     bg_vec   = vol_bg(valid_mask);
% 
%     data_vec(bg_vec) = NaN;
% 
%     sums   = accumarray(id_indices, data_vec, [n_unique_regions 1], @nansum);
%     counts = accumarray(id_indices, ~isnan(data_vec), [n_unique_regions 1], @sum);
% 
%     leaf_means_exp(:, m) = sums ./ counts;
% end
% 
% % 4. Compute T-Scores per Region
% fprintf('Computing T-scores...\n');
% 
% % Means across mice
% mu_ctrl = nanmean(leaf_means_ctrl, 2);
% mu_exp  = nanmean(leaf_means_exp, 2);
% diff_mu = mu_exp - mu_ctrl;
% 
% % SEM
% sem_c = nanstd(leaf_means_ctrl, [], 2) ./ sqrt(n_ctrl);
% sem_e = nanstd(leaf_means_exp, [], 2)  ./ sqrt(n_exp);
% pooled_sem = sqrt(sem_c.^2 + sem_e.^2);
% 
% % T-Score Vector
% t_scores_vec = diff_mu ./ pooled_sem;
% 
% % Filter: Set regions with NaN results (no data points) to 0
% t_scores_vec(isnan(t_scores_vec)) = 0;
% 
% % 5. Reconstruct 3D T-Score Map
% fprintf('Reconstructing 3D T-score map...\n');
% 
% % Initialize empty volume
% t_score_vol = zeros(size(atlas_left), 'single');
% 
% % Map T-scores back to pixels
% % We create a lookup table: 1..N -> T-score
% % Then map the id_indices image through this table
% t_score_valid_pixels = t_scores_vec(id_indices);
% 
% % Fill the volume
% t_score_vol(valid_mask) = t_score_valid_pixels;
% 
% % Visualization: Plot T-Score Map (Montage)
% 
% fprintf('Plotting Whole-Brain T-Map...\n');
% 
% figure('Name', ['WholeBrain_TMap_' comp_tag], 'Color', 'k', 'Position', [50 50 1200 900]);
% 
% % Parameters for visualization
% slices_to_show = 150:50:size(t_score_vol,1)-150; 
% n_cols = 5;
% n_rows = ceil(length(slices_to_show)/n_cols);
% 
% % Global CLim for T-scores (Symmetric around 0)
% t_lims = [-4 4]; 
% 
% % Create Custom Colormap (Blue -> White -> Red)
% n_steps = 256;
% c_blue = [0 0 1];
% c_white = [1 1 1];
% c_red = [1 0 0];
% cmap_neg = [linspace(c_blue(1), c_white(1), n_steps/2)', ...
%             linspace(c_blue(2), c_white(2), n_steps/2)', ...
%             linspace(c_blue(3), c_white(3), n_steps/2)'];
% cmap_pos = [linspace(c_white(1), c_red(1), n_steps/2)', ...
%             linspace(c_white(2), c_red(2), n_steps/2)', ...
%             linspace(c_white(3), c_red(3), n_steps/2)'];
% custom_cmap = [cmap_neg; cmap_pos];
% 
% for k = 1:length(slices_to_show)
%     s_idx = slices_to_show(k);
% 
%     subplot(n_rows, n_cols, k);
% 
%     % Extract slice data
%     im_slice = squeeze(t_score_vol(s_idx, :, :));
% 
%     % Extract Mask (Atlas boundaries)
%     mask_slice = squeeze(atlas_left(s_idx, :, :));
% 
%     % Create Alpha Data: 1 for Brain, 0 for Background
%     alpha_data = double(mask_slice > 0);
% 
%     % Plot
%     imagesc(im_slice);
%     axis image; axis off;
% 
%     % Apply Transparency to Background
%     % This makes the 0 values outside the brain invisible (showing black figure bg)
%     set(gca, 'Color', 'k'); % Set axes background to black
%     h_im = findobj(gca, 'Type', 'image');
%     set(h_im, 'AlphaData', alpha_data);
% 
%     % Apply Colormap & Limits
%     colormap(gca, custom_cmap); 
%     clim(t_lims);
% 
%     % Add Title per slice
%     title(['Slice ' num2str(s_idx)], 'Color', 'w', 'FontSize', 8);
% end
% 
% % Add Colorbar
% c = colorbar;
% c.Position = [0.92 0.1 0.02 0.8];
% c.Color = 'w';
% c.Label.String = ['T-score (' exp_type ' - ' ctrl_type ')'];
% colormap(c, custom_cmap); % Ensure colorbar uses the custom map
% 
% sgtitle(['Whole-Brain Region-based T-Scores (Leaf Nodes)'], 'Color', 'w', 'FontSize', 14);
% 
% % Save Results
% set(gcf, 'InvertHardcopy', 'off'); 
% saveas(gcf, fullfile(comp_out_dir, ['WholeBrain_TMap_Montage_' comp_tag '.fig']));
% % 'BackgroundColor', 'current' ensures the transparency reveals the black figure color in PNG
% exportgraphics(gcf, fullfile(comp_out_dir, ['WholeBrain_TMap_Montage_' comp_tag '.png']), ...
%     'Resolution', 300, 'BackgroundColor', 'current');
% 
% % Also Save the 3D Volume for external viewing
% vol_filename = fullfile(comp_out_dir, ['T_Score_Map_' comp_tag '.mat']);
% save(vol_filename, 't_score_vol', 'unique_ids', 't_scores_vec');
% fprintf('Whole-brain analysis complete. Volume saved to %s\n', vol_filename);


%% 


%% Area-based Analysis: Vectorized Target Region Quantification (Diff & Sum)

fprintf('Starting Area-based quantification (Vectorized)...\n');

% 0. Define metric to plot
choosen_metric = 'P999';

% 1. Define Target Regions
roi_list = { ...
    'Primary somatosensory area, barrel field', ...
    'Primary somatosensory area, trunk', ...
    'Primary somatosensory area, upper limb', ...
    'Primary somatosensory area, lower limb', ...
    'Supplemental somatosensory area', ...
    'Primary motor area', ...
    'Secondary motor area', ...
    'Primary visual area', ...
    'Lateral visual area', ...
    'Anterolateral visual area', ...
    'Anteromedial visual area', ...
    'Primary auditory area', ...
    'Dorsal auditory area', ...
    'Ventral auditory area', ...
    'Anterior cingulate area', ...
    'Olfactory tubercle', ...
    'Prelimbic area', ...
    'Infralimbic area', ...
    'Visceral area',...
    'Gustatory areas',...
    'Piriform area',...
    'Subiculum',...
    'Orbital area', ...
    'Claustrum', ...
    'Agranular insular area', ...
    'Anterior area', ...
    'Rostrolateral visual area', ...
    'Temporal association areas', ...
    'Perirhinal area', ...
    'Ectorhinal area', ...
    'Retrosplenial area', ...
    'Nucleus accumbens', ...
    'Caudoputamen', ...
    'Globus pallidus, external segment', ...
    'Subthalamic nucleus', ...
    'Ventral posteromedial nucleus of the thalamus', ...
    'Ventral posterolateral nucleus of the thalamus', ...
    'Ventral medial nucleus of the thalamus',...
    'Zona incerta',...
    'Posterior complex of the thalamus', ...
    'Lateral posterior nucleus of the thalamus', ...
    'Lateral dorsal nucleus of thalamus', ...
    'Ventral anterior-lateral complex of the thalamus', ...
    'Mediodorsal nucleus of the thalamus', ...
    'Parafascicular nucleus', ...
    'Nucleus of reuniens', ...
    'Central lateral nucleus of the thalamus', ...
    'Reticular nucleus of the thalamus', ...
    'Geniculate group, dorsal thalamus', ...
    'Midbrain, motor related', ...
    'Superior colliculus, motor related', ...
    'Superior colliculus, sensory related', ...
    'Hippocampal formation', ...
    'Basolateral amygdalar nucleus', ...
    'Parafascicular nucleus', ...
    'Hypothalamus' ...
    };
n_rois = length(roi_list);

% 2. Setup Atlas Vectors
fprintf('Mapping Atlas Volume...\n');
atlas_left = AllenCrop(:, :, 1:half_width);
valid_pixels = atlas_left > 0;
pixel_ids = atlas_left(valid_pixels);

% 3. Pre-calculate Masks
fprintf('Building masks for %d regions (Handling overlaps)...\n', n_rois);

roi_masks = false(length(pixel_ids), n_rois);
roi_pixel_counts = zeros(n_rois, 1);

for r = 1:n_rois
    region_name = roi_list{r};
    try
        % Get hierarchical IDs for this region name
        mask_temp = get_allen_region_mask(allenDir, AllenCrop, {region_name}, brainMask, '');
        if size(mask_temp, 3) >= half_width
            mask_temp = mask_temp(:, :, 1:half_width);
        end
        % Get unique IDs belonging to this region
        ids_in_region = unique(atlas_left(mask_temp));
        ids_in_region(ids_in_region == 0) = [];
        % Create boolean mask for the flattened pixel vector
        roi_masks(:, r) = ismember(pixel_ids, ids_in_region);
        roi_pixel_counts(r) = sum(roi_masks(:, r));
    catch
        warning('Could not map region: %s', region_name);
    end
end

% Warning for empty regions
empty_rois = find(roi_pixel_counts == 0);
if ~isempty(empty_rois)
    fprintf('Warning: The following regions have 0 pixels:\n');
    disp(roi_list(empty_rois)');
end

% 4. Initialize Stats Storage
metric_names = {'Mean', 'Median', 'Q1', 'Q3', 'P99', 'P999'};
n_metrics = length(metric_names);
roi_stats_ctrl = nan(n_rois, n_ctrl, n_metrics);
roi_stats_exp  = nan(n_rois, n_exp, n_metrics);

% Pre-slice background masks
bg_mask_ctrl = mask_bg_ctrl;
bg_mask_exp  = mask_bg_exp;

% 5. Analysis Loop
analysis_types = {'Diff', 'Sum'};

for a_idx = 1:length(analysis_types)
    res_type = analysis_types{a_idx};
    fprintf('--- Processing %s Data for ROIs ---\n', res_type);

    if strcmp(res_type, 'Diff')
        input_vol_ctrl = abs(lr_diff_ctrl);
        input_vol_exp  = abs(lr_diff_exp);
    else
        input_vol_ctrl = abs(lr_sum_ctrl);
        input_vol_exp  = abs(lr_sum_exp);
    end

    % --- Control Group ---
    for m = 1:n_ctrl
        tic
        % Flatten Data for this mouse
        vol_data = input_vol_ctrl(:,:,:,m); 
        vol_bg   = bg_mask_ctrl(:,:,:,m);
        
        data_vec = vol_data(valid_pixels);
        bg_vec   = vol_bg(valid_pixels);
        data_vec(bg_vec) = NaN;
        
        % Loop over ROIs
        for r = 1:n_rois
            if roi_pixel_counts(r) == 0, continue; end
            % Extract pixels for this ROI
            roi_vals = data_vec(roi_masks(:, r));
            % Compute Stats
            roi_stats_ctrl(r, m, 1) = mean(roi_vals, 'omitnan'); % Mean
            % Compute other metrics (only if needed, slightly slower)
            roi_vals_clean = roi_vals(~isnan(roi_vals));
            if ~isempty(roi_vals_clean)
                roi_stats_ctrl(r, m, 2) = median(roi_vals_clean);
                roi_stats_ctrl(r, m, 3) = quantile(roi_vals_clean, 0.25);
                roi_stats_ctrl(r, m, 4) = quantile(roi_vals_clean, 0.75);
                roi_stats_ctrl(r, m, 5) = quantile(roi_vals_clean, 0.99);
                roi_stats_ctrl(r, m, 6) = quantile(roi_vals_clean, 0.999);
            end
        end
        toc
    end

    % --- Experimental Group ---
    for m = 1:n_exp
        tic
        vol_data = input_vol_exp(:,:,:,m); 
        vol_bg   = bg_mask_exp(:,:,:,m);
        data_vec = vol_data(valid_pixels);
        bg_vec   = vol_bg(valid_pixels);
        data_vec(bg_vec) = NaN;
        for r = 1:n_rois
            if roi_pixel_counts(r) == 0, continue; end
            roi_vals = data_vec(roi_masks(:, r));
            roi_stats_exp(r, m, 1) = mean(roi_vals, 'omitnan');
            roi_vals_clean = roi_vals(~isnan(roi_vals));
            if ~isempty(roi_vals_clean)
                roi_stats_exp(r, m, 2) = median(roi_vals_clean);
                roi_stats_exp(r, m, 3) = quantile(roi_vals_clean, 0.25);
                roi_stats_exp(r, m, 4) = quantile(roi_vals_clean, 0.75);
                roi_stats_exp(r, m, 5) = quantile(roi_vals_clean, 0.99);
            end
        end
        toc
    end

    % --- Plotting ---
    target_metric_idx = find(strcmp(metric_names, choosen_metric));
    if isempty(target_metric_idx), target_metric_idx = 1; end
    
    roi_means_ctrl = roi_stats_ctrl(:, :, target_metric_idx);
    roi_means_exp  = roi_stats_exp(:, :, target_metric_idx);

    mean_ctrl_roi = nanmean(roi_means_ctrl, 2); 
    mean_exp_roi  = nanmean(roi_means_exp, 2); 
    diff_means    = mean_exp_roi - mean_ctrl_roi;

    sem_c = nanstd(roi_means_ctrl,[],2) ./ sqrt(n_ctrl);
    sem_e = nanstd(roi_means_exp, [],2) ./ sqrt(n_exp);
    pooled_sem = sqrt(sem_c.^2 + sem_e.^2);
    t_score_roi = diff_means ./ pooled_sem;

    valid_rows = ~isnan(t_score_roi);
    t_score_roi = t_score_roi(valid_rows);
    current_rois = roi_list(valid_rows);

    [sorted_t, sort_idx] = sort(t_score_roi, 'ascend'); 
    sorted_rois = current_rois(sort_idx);

    fig_bars = figure('Name', ['Region_Analysis_BarChart_' res_type '_' comp_tag], ...
        'Color', 'w', 'Units', 'Normalized', 'Position', [0 0 0.9 0.9]);

    b = barh(sorted_t); b.FaceColor = 'flat';
    for k = 1:length(sorted_t)
        if sorted_t(k) > 0, b.CData(k,:) = [0.85 0.33 0.10]; else, b.CData(k,:) = [0 0.45 0.74]; end
    end

    yticks(1:length(sorted_rois)); yticklabels(sorted_rois);
    xlabel(['t-score (' exp_type ' - ' ctrl_type ')']);
    title(['Regional LR - ' res_type ' differences - ', strrep(comp_tag,'_',' '), ' - ', choosen_metric]);
    grid on; set(gca, 'FontSize', 10); ylim([0 length(sorted_rois)+1]);

    % Highlighting
    switch exp_type
        case 'rws', highlighted_areas = {'Primary somatosensory area, barrel field', 'Ventral posteromedial nucleus of the thalamus', 'Posterior complex of the thalamus', 'Supplemental somatosensory area', 'Zona incerta'};
        case 'behavior', highlighted_areas = {'Primary somatosensory area, barrel field', 'Ventral posteromedial nucleus of the thalamus', 'Posterior complex of the thalamus', 'Supplemental somatosensory area', 'Zona incerta'};
        otherwise, highlighted_areas = {};
    end
    
    ax = gca; 
    ytl = ax.YTickLabel; colored_labels = repmat({''}, size(ytl));
    for i = 1:length(ytl)
        if ismember(ytl{i}, highlighted_areas), colored_labels{i} = ['\color{magenta} \bf ' strrep(ytl{i}, '_', ' ') ]; else, colored_labels{i} = ['\color{black} ' strrep(ytl{i}, '_', ' ') ]; end
    end
    ax = gca; ax.YTickLabel = colored_labels;

    % Significance lines
    n_regions_tested = length(sorted_t);
    df = n_ctrl + n_exp - 2;
    t_crit_bonf = tinv(1 - 0.05/(2*n_regions_tested), df);
    hold on; xline(t_crit_bonf, 'k:', 'LineWidth', 2); xline(-t_crit_bonf, 'k:', 'LineWidth', 2); hold off;

    set(fig_bars, 'InvertHardcopy', 'off');
    saveas(fig_bars, fullfile(comp_out_dir, ['Region_Stats_Bar_' res_type '_' comp_tag '.fig']));
    exportgraphics(fig_bars, fullfile(comp_out_dir, ['Region_Stats_Bar_' res_type '_' comp_tag '.png']), 'Resolution', 300);
end

clear data_vec roi_masks atlas_left input_vol_ctrl input_vol_exp vol vol_smoothed bg_vec roi_vals roi_vals_clean
fprintf(' Area-based analysis saved to: %s\n', comp_out_dir);




%% 


%% Produce individual LR-asymmetry videos

fprintf('Generating Individual LR-Asymmetry Video (Diff & Sum)...\n');

% 1. Setup Groups & Names
group_list = {ctrl_type, exp_type};
group_mousename_list = {ctrl_mousenames, exp_mousenames};

% Use left-hemi atlas for boundaries
atlas_left = double(AllenCrop(:, :, 1:size(lr_diff_ctrl,3)));

% Visualization Params
clim_diff = [0 5];
clim_sum  = [0 10];

% --- Loop Over Groups (Create 2 separate videos) ---
for g_idx = 1:2

    curr_group_name = group_list{g_idx};
    curr_mouse_names = group_mousename_list{g_idx};

    fprintf('Generating Video for Group: %s ...\n', curr_group_name);

    % Select Data based on group
    if g_idx == 1
        % Control
        curr_diff_data = abs(lr_diff_ctrl);
        curr_sum_data  = lr_sum_ctrl;
        curr_mask_bg   = mask_bg_ctrl;
    else
        % Experimental
        curr_diff_data = abs(lr_diff_exp);
        curr_sum_data  = lr_sum_exp;
        curr_mask_bg   = mask_bg_exp;
    end

    n_mice_in_group = size(curr_diff_data, 4);
    n_slices = size(curr_diff_data, 1);

    % Setup Video Writer
    video_filename = ['Individual_LR_Diff_Sum_' curr_group_name '.mp4'];
    full_video_path = fullfile(comp_out_dir, video_filename);

    vidObj = VideoWriter(full_video_path, 'MPEG-4');
    vidObj.FrameRate = 15;
    vidObj.Quality = 95;
    open(vidObj);

    % Create Figure (Black Background)
    fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'k');
    set(fh, 'InvertHardcopy', 'off');

    % Iterate Slices
    for s_idx = 1:n_slices

        % Skip empty atlas slices
        atlas_slice = squeeze(atlas_left(s_idx, :, :));
        if sum(atlas_slice(:) > 0) == 0
            if mod(s_idx, 50) == 0, fprintf('  Skipping slice %d\n', s_idx); end
            continue;
        end

        % Compute Atlas Boundaries
        [gy, gx] = gradient(atlas_slice);
        boundaries = (abs(gx) + abs(gy)) > 0 & (atlas_slice > 0);
        [b_row, b_col] = find(boundaries);

        % --- Loop over Mice in Current Group ---
        for k = 1:n_mice_in_group

            % Get specific mouse name
            mouse_name = strrep(curr_mouse_names{k}, '_', ' ');

            % Get data
            data_diff = squeeze(curr_diff_data(s_idx, :, :, k));
            data_sum  = squeeze(curr_sum_data(s_idx, :, :, k));
            mask_bg   = squeeze(curr_mask_bg(s_idx, :, :, k));

            % Combine Atlas Mask (must be in brain) AND Tissue Mask (must not be background)
            valid_pixels = (atlas_slice > 0) & (~mask_bg);
            alpha_data = double(valid_pixels);

            % --- Row 1: LR Difference ---
            subplot(2, n_mice_in_group, k);
            imagesc(data_diff);
            clim(clim_diff);
            colormap(gca, hot);
            set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
            axis image; axis off; set(gca, 'Color', 'k');
            hold on;
            title(mouse_name, 'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');
            if k==1
                ylabel('LR Diff', 'Color','w','FontSize',12,'FontWeight','bold');
            end
            cb = colorbar;
            cb.Label.String = '|L - R|';
            cb.Color = 'w'; cb.Label.Color = 'w';
            cb.Location = 'eastoutside';

            % --- Row 2: LR Sum ---
            subplot(2, n_mice_in_group, k + n_mice_in_group);
            imagesc(data_sum);
            clim(clim_sum);
            colormap(gca, hot);
            set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
            axis image; axis off; set(gca, 'Color', 'k');
            hold on;

            if k==1
                ylabel('LR Sum', 'Color','w','FontSize',12,'FontWeight','bold');
            end
            cb = colorbar;
            cb.Label.String = 'L + R';
            cb.Color = 'w'; cb.Label.Color = 'w';
            cb.Location = 'eastoutside';

        end
        % Super Title
        sgtitle(['Slice # ' num2str(s_idx) ' - Left Hemisphere Asymmetry (' curr_group_name ')'], 'Color', 'w', 'FontSize', 14);

        frame = getframe(fh);
        writeVideo(vidObj, frame);
        clf(fh);
        if mod(s_idx, 50) == 0
            fprintf('  Frame %d written...\n', s_idx);
        end
    end
    close(vidObj);
    close(fh);
    clear curr_diff_data curr_sum_data curr_mask_bg
    fprintf('Video saved: %s\n', full_video_path);
end
