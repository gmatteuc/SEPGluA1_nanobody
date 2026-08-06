clear all
close all
clc

% /// Pipeline script #2: read centered volumes and perform residual correction anlayisis to bring autofluorecence and nano channel on the same scale  /// 

%% User-defined parameters

mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', 'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1','MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', 'MG725_Gria1', 'MG727_Gria1'};
mousetypes = {'rws','rws','rws','rws','rws','naive','naive','naive','naive','naive','behavior','behavior','behavior','behavior','behavior','behavior','behavior'};
doPlotBkg = true;
savePlotBkg = true;
saveRatioMap = false;

%% Add paths 

allenDir = 'D:\sep_histology\data\atlas';
lightsuiteDir = 'D:\sep_histology\code\LightSuite-main';
yamlDir = 'D:\sep_histology\code\yamlmatlab'; 
elastixDir = 'D:\sep_histology\code\matlab_elastix-master';
addpath(allenDir)
addpath(genpath(lightsuiteDir))
addpath(genpath(yamlDir))
addpath(genpath(elastixDir))

%% Load atlas

AllenFile = fullfile(allenDir, 'annotation_10.nii.gz');
AllenVol = niftiread(AllenFile);

%% Loop over mice to perform residual analysis (signal correction)

for mouse_idx = [12,13,15,17]; %11%[6,9,10,14,16] %1:numel(mice)  

    % Get curren mouse name and type
    mouse_name = mice{mouse_idx};
    mouse_type = mousetypes{mouse_idx};

    % Get dirs
    base_dir = ['D:\sep_histology\data\', mouse_type, '\'];
    output_dir = fullfile(base_dir, mouse_name, 'lightsuite', 'correction_output');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    plotDir = fullfile(output_dir, 'diagnostic_plots');
    if ~exist(plotDir, 'dir')
        mkdir(plotDir);
    end
    autoPath_centered = fullfile(base_dir, mouse_name, 'lightsuite', 'volume_centered', 'chan03_Cy3.tiff'); % Autofluorescence centered tiff
    nanoPath_centered = fullfile(base_dir, mouse_name, 'lightsuite', 'volume_centered', 'chan02_Cy5.tiff'); % Nanobody centered tiff

    %% Load centered volumes

    tic
    % Load autofluorescence data
    autoVol_centered = loadVolume({autoPath_centered}, 1);
    % Load nanobody data
    nanoVol_centered = loadVolume({nanoPath_centered}, 1);
    toc

    %% Analyze relationship between nano and auto channel

    % Select leading channel
    selectedVol = autoVol_centered; % Autofluorescence (base, I)
    selectedVolSig = nanoVol_centered; % Nanobody (signal, J)

    % Preallocate outputs
    [H, W, Z] = size(selectedVol);
    bg_mask_vol = false(H, W, Z);
    fprintf('Starting within-slice correction analysis on %d slices for %s...\n', Z, mouse_name);
    t0 = tic;

    for z = 1:Z

        % 0) Select input slices and cast to single
        I = selectedVol(:,:,z);
        I = im2single(I);
        J = selectedVolSig(:,:,z);
        J = im2single(J);

        % 1) Select reference pixels for current slice
        range_frac = 0.20;
        if z<10
            rangewinmax = 75;
        else
            rangewinmax = 50;
        end
        [ref_pix_mask_J, ~, bg_mask, ~, used_clim, h_diag_J] = select_reference_pixels(J, 15, rangewinmax, 60, range_frac, doPlotBkg);
        if savePlotBkg
            saveas(h_diag_J, fullfile(plotDir, sprintf('reference_pix_selection_slice_J_%03d.png', z)));
        end
        ref_pix_mask = ref_pix_mask_J;
        bg_mask_vol(:,:,z) = bg_mask;

        % 2) Get values of current reference pixels
        basepix = I(ref_pix_mask);
        sigpix = J(ref_pix_mask);
        base_sig_ratio = nanmean(sigpix ./ basepix); 

        % Visualize result
        figure('Name', sprintf('Analysis Slice %03d', z), 'units', 'normalized', 'outerposition', [0 0 1 1]);
        subplot(1,3,1);
        scatter(basepix, sigpix, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.1);
        hold on;

        % Compute shaded error region with robust regression
        if numel(basepix) > 1

            [b, stats] = robustfit(basepix, sigpix, 'bisquare');
            slope = b(2);
            intercept = b(1);
            y_fit = slope * basepix + intercept;

            % Approximate 95% confidence interval
            resid = sigpix - y_fit;
            resid_std = sqrt(sum(resid.^2) / (length(basepix) - 2));
            df = length(basepix) - 2;
            t_crit = tinv(0.975, df);
            mean_x = mean(basepix);
            delta = t_crit * resid_std * sqrt(1 + 1/length(basepix) + (basepix - mean_x).^2 / sum((basepix - mean_x).^2));

            % Sort for smooth patch
            [basepix_sorted, idx] = sort(basepix);
            y_fit_sorted = y_fit(idx);
            delta_sorted = delta(idx);

            % Patch for CI
            x_patch = [basepix_sorted; flipud(basepix_sorted); basepix_sorted(1)];
            y_patch = [y_fit_sorted + delta_sorted; flipud(y_fit_sorted - delta_sorted); y_fit_sorted(1) + delta_sorted(1)];
            patch(x_patch, y_patch, 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', '95% CI');

            % Fit line
            plot(basepix_sorted, y_fit_sorted, 'r-', 'LineWidth', 2, 'DisplayName', 'Fit');

            % Equation
            eq_str = sprintf('y = %.2f x + %.2f', round(slope, 2), round(intercept, 2));
            rangepixf = abs(0.1 * (max(basepix) - min(basepix)));
            xlim(gca, [min(basepix) - rangepixf, max(basepix) + rangepixf]);
            rangepixf = abs(0.1 * (max(sigpix) - min(sigpix)));
            ylim(gca, [min(sigpix) - rangepixf, max(sigpix) + rangepixf]);
            xlim_vals = get(gca, 'XLim');
            ylim_vals = get(gca, 'YLim');
            text(mean(xlim_vals), ylim_vals(1) + 0.05 * diff(ylim_vals), eq_str, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'Color', 'r', 'FontSize', 10);

        end
        xlim_vals = get(gca, 'XLim');
        ylim_vals = get(gca, 'YLim');
        diag_x = [min(xlim_vals), max(xlim_vals)];
        diag_y = diag_x;
        plot(diag_x, diag_y, '--', 'LineWidth', 1, 'Color', [1,0,1]);
        axis square
        xlabel('Base intensity (I)');
        ylabel('Signal intensity (J)');
        title(sprintf('Scatter plot - Slice %03d - Reference pix # %03d - Ratio %0.3f', z, numel(sigpix), base_sig_ratio));
        grid on;
        axis square;

        subplot(1,3,2);
        imagesc(I); axis image off; colormap gray; title('Base image (I)');
        clim(used_clim*1/slope);
        hold on;
        [rows, cols] = find(ref_pix_mask);
        if ~isempty(rows)
            x = [cols-0.5, cols+0.5, cols+0.5, cols-0.5]';
            y = [rows-0.5, rows-0.5, rows+0.5, rows+0.5]';
            faces = reshape(1:numel(cols)*4, 4, [])';
            patch('Faces', faces, 'Vertices', [x(:), y(:)], ...
                'FaceColor', 'r', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        end
        hold off;
        colorbar;

        subplot(1,3,3);
        imagesc(J); axis image off; colormap gray; title('Signal image (J)');
        clim(used_clim);
        hold on;
        [rows, cols] = find(ref_pix_mask);
        if ~isempty(rows)
            x = [cols-0.5, cols+0.5, cols+0.5, cols-0.5]';
            y = [rows-0.5, rows-0.5, rows+0.5, rows+0.5]';
            faces = reshape(1:numel(cols)*4, 4, [])';
            patch('Faces', faces, 'Vertices', [x(:), y(:)], ...
                'FaceColor', 'r', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        end
        hold off;
        colorbar;

        sgtitle('Reference pixels regression')
        saveas(gcf, fullfile(plotDir, sprintf('reference_pix_analysis_slice_%03d.png', z)));
        close(gcf);

        % Store pixel distributions and regression results
        if ~exist('slice_data', 'var')
            slice_data = struct();
        end
        slice_data(z).basepix = basepix;
        slice_data(z).sigpix = sigpix;
        if numel(basepix) > 1
            slice_data(z).p = b; % [intercept, slope]
            slice_data(z).delta = delta;
            slice_data(z).ratio = base_sig_ratio;
            slice_data(z).slope = slope;
            slice_data(z).intercept = intercept;
            slice_data(z).used_clim = used_clim;
        else
            slice_data(z).p = NaN(2, 1);
            slice_data(z).delta = NaN(size(basepix));
            slice_data(z).ratio = NaN;
            slice_data(z).slope = NaN;
            slice_data(z).intercept = NaN;
            slice_data(z).used_clim = [NaN, NaN];
        end
    end

    fprintf('Analysis done in %.1f s\n', toc(t0));

    %% Visualize aggregate results with jittered scatter plots and curves

    figure('Name', 'Jittered Scatter Plots and Curves of Regression Metrics', 'units', 'normalized', 'outerposition', [0 0 1 1]);

    % Subplot 1: Jittered scatter plots
    subplot(1,2,1);
    categories = {'Ratio', 'Slope', 'Intercept'};
    values = nan(Z, 3);
    for z = 1:Z
        if isfield(slice_data(z), 'ratio')
            values(z, 1) = slice_data(z).ratio;
            values(z, 2) = slice_data(z).slope;
            values(z, 3) = slice_data(z).intercept;
        end
    end
    valid_idx = ~isnan(values(:, 1));
    values = values(valid_idx, :);
    Z_valid = sum(valid_idx);
    slice_indices = find(valid_idx);

    % Add vertical jitter
    jitter = 0.1 * randn(Z_valid, 1);
    x_jittered = repmat(1:3, Z_valid, 1) + [jitter, jitter, jitter];

    % Colormap from black to dark gray
    colors = zeros(Z_valid, 3);
    for i = 1:3
        colors(:,i) = linspace(0, 0.5, Z_valid)';
    end
    colors_expanded = repmat(colors, 3, 1);

    % Scatter
    scatter(x_jittered(:), values(:), 50, colors_expanded, 'filled', 'MarkerFaceAlpha', 1);
    hold on;

    % Median and IQR
    medians = median(values, 1, 'omitnan');
    iqr_vals = iqr(values, 1);
    for i = 1:3
        if i == 3
            yyaxis right
            ylabel('Intercept');
            ylim([0, 10])
            ax = gca;
            ax.YColor = [0 0 0];
        else
            ylabel('Slope / Ratio');
            ylim([0, 2])
            ax = gca;
            ax.YColor = [0 0 0];
        end
        y = medians(i);
        y_min = y - iqr_vals(i) / 2;
        y_max = y + iqr_vals(i) / 2;

        if i == 3
            colline = [0, 0, 1];
        elseif i == 2
            colline = [1, 0, 0];
        elseif i == 1
            colline = [1, 0, 1];
        end
        plot([i-0.2 i+0.2], [y y], '-', 'LineWidth', 2, 'Color', colline);
        if ~isnan(y_min) && ~isnan(y_max) && isfinite(y_min) && isfinite(y_max)
            patch([i-0.2 i-0.2 i+0.2 i+0.2], [y_min y_max y_max y_min], colline, ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end

    set(gca, 'XTick', 1:3, 'XTickLabel', categories);
    xlabel('Channel relationship metric');
    title('Jittered Scatter Plots of Ratios, Slopes, and Intercepts');
    grid on;
    axis square;
    set(gca, 'fontsize', 12);

    % Subplot 2: Curves over slices
    subplot(1,2,2);
    hold on;
    yyaxis left
    plot(slice_indices, values(:, 1), '-', 'LineWidth', 1.5, 'DisplayName', 'Ratio', 'Color', [1, 0, 1]);
    plot(slice_indices, values(:, 2), '-', 'LineWidth', 1.5, 'DisplayName', 'Slope', 'Color', [1, 0, 0]);
    ylabel('Slope / Ratio');
    ylim([0, 2])
    ax = gca;
    ax.YColor = [0 0 0];
    yyaxis right
    plot(slice_indices, values(:, 3), '-', 'LineWidth', 1.5, 'DisplayName', 'Intercept', 'Color', [0, 0, 1]);
    ylabel('Intercept');
    ylim([0, 10])
    ax = gca;
    ax.YColor = [0 0 0];
    xlabel('Slice Index');
    title('Regression Metrics Across Slices');
    grid on;
    axis square;
    set(gca, 'fontsize', 12);

    % Get average slope and intercept
    average_slope = nanmean([slice_data.slope]); %#ok<*NANMEAN>
    average_intercept = nanmean([slice_data.intercept]);

    % Save
    saveas(gcf, fullfile(plotDir, 'jittered_regression_metrics_with_curves.png'));
    close(gcf);

    %% Apply corrections for both types

    correction_types = {'slicewise', 'global'};
    for ct = 1:numel(correction_types)

        correction_type = correction_types{ct};
        use_per_slice = strcmp(correction_type, 'slicewise');

        % Preallocate corrected volume
        correctedVol = zeros(H, W, Z, 'single');
        scaledautoVol = zeros(H, W, Z, 'single');
        nanoVol = zeros(H, W, Z, 'single');

        fprintf('Applying %s correction for %s...\n', correction_type, mouse_name);
        t0 = tic;
        for z = 1:Z

            I = im2single(selectedVol(:,:,z));
            J = im2single(selectedVolSig(:,:,z));

            if use_per_slice
                slope = slice_data(z).slope;
                intercept = slice_data(z).intercept;
            else
                slope = average_slope;
                intercept = average_intercept;
            end

            % Scale auto and subtract
            scaled_I = slope * I + intercept;
            corrected = J - scaled_I;
            corrected(corrected < 0) = 0;
            correctedVol(:,:,z) = corrected;
            scaledautoVol(:,:,z) = scaled_I;
            nanoVol(:,:,z) = J;

        end
        fprintf('Correction done in %.1f s\n', toc(t0));

        % Save matfile
        matfile_name = fullfile(output_dir, sprintf('corrected_volume_%s.mat', correction_type));
        save(matfile_name, 'correctedVol', 'bg_mask_vol', 'slice_data', 'average_slope', 'average_intercept', 'correction_type', '-v7.3');
        % Save matfile
        matfile_name_bis = fullfile(output_dir, sprintf('scaled_auto_volume_%s.mat', correction_type));
        save(matfile_name_bis, 'scaledautoVol', 'nanoVol', 'bg_mask_vol', 'slice_data', 'average_slope', 'average_intercept', 'correction_type', '-v7.3');

        % Write scaled difference video ( relative difference: (J - scaled_I) / scaled_I )
        videoFile_diff = fullfile(output_dir, sprintf('scaled_difference_video_%s.mp4', correction_type));
        vidObj_diff = VideoWriter(videoFile_diff, 'MPEG-4');
        vidObj_diff.FrameRate = 1;
        vidObj_diff.Quality = 100;
        open(vidObj_diff);

        % Precompute global color limits
        clim_min_diff = -1;
        clim_max_diff = 1; % For relative

        for z = 1:Z
            I = im2single(selectedVol(:,:,z));
            J = im2single(selectedVolSig(:,:,z));

            if use_per_slice
                slope = slice_data(z).slope;
                intercept = slice_data(z).intercept;
            else
                slope = average_slope;
                intercept = average_intercept;
            end

            scaled_I = slope * I + intercept;
            diff_map = (J - scaled_I) ./ scaled_I; % NB: relative scaled difference

            % Background mask
            diff_map(bg_mask_vol(:,:,z)) = NaN;

            figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
            h = imagesc(diff_map);
            clim(gca, [clim_min_diff, clim_max_diff]);
            colorbar;
            colormap(get_color2color_colormap([0,0,1],[1,0,0]));

            ax = gca;
            ax.Color = 'k';
            alpha_mask = ~isnan(diff_map);
            set(h, 'AlphaData', alpha_mask);
            axis equal;
            grid on;
            ylim([0, size(I, 1)]);
            xlim([0, size(I, 2)]);
            if use_per_slice
                title(sprintf('Scaled difference map (relative) - Slice # %d (slicewise %0.2f)', z, slope));
            else
                title(sprintf('Scaled difference map (relative) - Slice # %d (global %0.2f)', z, slope)); %#ok<*UNRCH>
            end

            frame = getframe(gcf);
            writeVideo(vidObj_diff, frame);
            close(gcf);
        end
        close(vidObj_diff);

    end

    %% Write ratio map video (J / I)

    if saveRatioMap

        videoFile_ratio = fullfile(output_dir, 'ratio_map_video.mp4');
        vidObj_ratio = VideoWriter(videoFile_ratio, 'MPEG-4');
        vidObj_ratio.FrameRate = 1;
        vidObj_ratio.Quality = 100;
        open(vidObj_ratio);

        % Precompute global color limits
        clim_min_ratio = 0;
        clim_max_ratio = 2;

        for z = 1:Z
            I = im2single(selectedVol(:,:,z));
            J = im2single(selectedVolSig(:,:,z));
            % Background mask
            I_mask = I;
            I_mask(bg_mask_vol(:,:,z)) = NaN;
            J_mask = J;
            J_mask(bg_mask_vol(:,:,z)) = NaN;
            ratio_map = J_mask ./ I_mask;
            figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
            h = imagesc(ratio_map);
            clim(gca, [clim_min_ratio, clim_max_ratio]);
            colorbar;
            colormap(get_color2color_colormap([1,0,0],[0,0,1]));
            ax = gca;
            ax.Color = 'k';
            alpha_mask = ~isnan(ratio_map);
            set(h, 'AlphaData', alpha_mask);
            axis equal;
            grid on;
            ylim([0, size(I, 1)]);
            xlim([0, size(I, 2)]);
            title(sprintf('Nano / Auto ratio map - Slice # %d', z));
            frame = getframe(gcf);
            writeVideo(vidObj_ratio, frame);
            close(gcf);
        end
        close(vidObj_ratio);

    end

    % Clear per-mouse variables to save memory
    clear autoVol_centered dapiVol_centered nanoVol_centered autoVol_registered dapiVol_registered nanoVol_registered selectedVol selectedVolSig slice_data bg_mask_vol correctedVol scaledautoVol
    close all

end