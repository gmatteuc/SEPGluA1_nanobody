function residual_correction(run_settings)
%RESIDUAL_CORRECTION  Scale the autofluorescence onto the nano channel, slice by slice.
%   RESIDUAL_CORRECTION(run_settings) does the work of
%   run_residual_correction, which sets the fields of run_settings (paths,
%   groups_to_process, mice_to_process, atlas_key, doPlotBkg, savePlotBkg,
%   saveRatioMap) and says what each one does.

% settings of run_residual_correction, under the names the code below uses
paths = run_settings.paths;
groups_to_process = run_settings.groups_to_process;
mice_to_process = run_settings.mice_to_process;
atlas_key = run_settings.atlas_key;
doPlotBkg = run_settings.doPlotBkg;
savePlotBkg = run_settings.savePlotBkg;
saveRatioMap = run_settings.saveRatioMap;

%% Add paths

% the toolboxes are on the path from sep_setup_paths; only the atlas folder is
% added here
atlas = get_atlas(atlas_key);
atlas_dir = atlas.dir;
addpath(atlas_dir)

%% Resolve cohort

% check that the registry still lists the adults in their legacy order
get_cohort('verify');

% the mice named, or else every mouse of the groups
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end
fprintf('run_residual_correction: %d mouse/mice selected.\n', numel(cohort));

%% Correct each mouse

for mouse_idx = 1:numel(cohort)

    % correct the mouse, with its figures and videos
    correct_mouse(cohort, mouse_idx, paths, doPlotBkg, savePlotBkg, saveRatioMap);

end

end

% ===== Local functions: per-slice fit =====

function correct_mouse(cohort, mouse_idx, paths, doPlotBkg, savePlotBkg, saveRatioMap)
% One mouse: the nano fitted on the autofluorescence slice by slice, both
% corrections saved, with their figures and videos.

% the mouse's name and group
mouse_name = cohort(mouse_idx).name;
mouse_type = cohort(mouse_idx).group;

% its output folders, and its centred autofluorescence and nano volumes
base_dir = fullfile(paths.data, mouse_type);
output_dir = fullfile(base_dir, mouse_name, 'lightsuite', 'correction_output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
plotDir = fullfile(output_dir, 'diagnostic_plots');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end
autoPath_centered = fullfile(base_dir, mouse_name, 'lightsuite', 'volume_centered', ...
    'chan03_Cy3.tiff');
nanoPath_centered = fullfile(base_dir, mouse_name, 'lightsuite', 'volume_centered', ...
    'chan02_Cy5.tiff');

%% Load the centred volumes

% autofluorescence, then nano, timed
tic
autoVol_centered = loadVolume({autoPath_centered}, 1);
nanoVol_centered = loadVolume({nanoPath_centered}, 1);
toc

%% Fit nano on autofluorescence, slice by slice

% the autofluorescence is the base (I), the nano the signal (J)
selectedVol = autoVol_centered;
selectedVolSig = nanoVol_centered;

% each slice's fit on its reference pixels, and the background of every slice
[H, W, Z] = size(selectedVol);
[slice_data, bg_mask_vol] = fit_reference_pixels(selectedVol, selectedVolSig, ...
    H, W, Z, mouse_name, plotDir, doPlotBkg, savePlotBkg);

%% Summary figure of the fits

% ratio, slope and intercept over the slices; the mean fit for the global way
[average_slope, average_intercept] = plot_regression_summary(slice_data, Z, plotDir);

%% Correct both ways

% 'slicewise' with each slice's fit, 'global' with the mean fit
correction_types = {'slicewise', 'global'};
for ct = 1:numel(correction_types)

    correction_type = correction_types{ct};
    use_per_slice = strcmp(correction_type, 'slicewise');

    correct_and_save(selectedVol, selectedVolSig, slice_data, bg_mask_vol, ...
        average_slope, average_intercept, correction_type, use_per_slice, ...
        mouse_name, output_dir, H, W, Z);

    write_difference_video(selectedVol, selectedVolSig, slice_data, bg_mask_vol, ...
        average_slope, average_intercept, correction_type, use_per_slice, ...
        output_dir, Z);

end

%% Ratio video

% nano / autofluorescence (J / I), slice by slice
if saveRatioMap

    write_ratio_video(selectedVol, selectedVolSig, bg_mask_vol, Z, output_dir);

end

% free the mouse's volumes and close its figures
clear autoVol_centered dapiVol_centered nanoVol_centered autoVol_registered ...
    dapiVol_registered nanoVol_registered selectedVol selectedVolSig slice_data ...
    bg_mask_vol correctedVol scaledautoVol
close all

end

function [slice_data, bg_mask_vol] = fit_reference_pixels(selectedVol, selectedVolSig, ...
    H, W, Z, mouse_name, plotDir, doPlotBkg, savePlotBkg)
% Fit nano on autofluorescence over each slice's reference pixels, with one
% figure per slice saved; the fits, and the background mask of the volume.

% the background of every slice, filled in the loop
bg_mask_vol = false(H, W, Z);
fprintf('Starting within-slice correction analysis on %d slices for %s...\n', Z, ...
    mouse_name);
t0 = tic;

for z = 1:Z

    % the slice of each channel, as single
    I = selectedVol(:, :, z);
    I = im2single(I);
    J = selectedVolSig(:, :, z);
    J = im2single(J);

    % reference pixels and background on the nano slice (select_reference_pixels):
    % the knee searched between the 15th percentile and the 75th in the first nine
    % slices, the 50th after (reason not recorded), the background dilated by a
    % 60-pixel disk, the reference pixels 20% of the way up from the knee's dip
    range_frac = 0.20;
    if z<10
        rangewinmax = 75;
    else
        rangewinmax = 50;
    end
    [ref_pix_mask_J, ~, bg_mask, ~, used_clim, h_diag_J] = select_reference_pixels( ...
        J, 15, rangewinmax, 60, range_frac, doPlotBkg);

    % its figure, saved when drawn, then closed, so one figure per slice does not
    % stay open until the end of the mouse
    if doPlotBkg && savePlotBkg
        saveas(h_diag_J, fullfile(plotDir, ...
            sprintf('reference_pix_selection_slice_J_%03d.png', z)));
        close(h_diag_J);
    end

    % the reference pixels and the background found on the nano serve both channels
    ref_pix_mask = ref_pix_mask_J;
    bg_mask_vol(:, :, z) = bg_mask;

    % the reference pixels' values, and their mean nano / autofluorescence ratio
    basepix = I(ref_pix_mask);
    sigpix = J(ref_pix_mask);
    base_sig_ratio = nanmean(sigpix ./ basepix);

    % scatter of the reference pixels, with the fit when there are two or more
    figure('Name', sprintf('Analysis Slice %03d', z), 'units', 'normalized', ...
        'outerposition', [0 0 1 1]);
    subplot(1, 3, 1);
    scatter(basepix, sigpix, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.1);
    hold on;
    if numel(basepix) > 1

        % nano on autofluorescence, a robust (bisquare) line, so that stray bright
        % pixels weigh less than in least squares
        [b, slope, intercept, delta] = draw_robust_fit(basepix, sigpix);

    end

    % the identity line
    xlim_vals = get(gca, 'XLim');
    ylim_vals = get(gca, 'YLim');
    diag_x = [min(xlim_vals), max(xlim_vals)];
    diag_y = diag_x;
    plot(diag_x, diag_y, '--', 'LineWidth', 1, 'Color', [1, 0, 1]);
    axis square
    xlabel('Base intensity (I)');
    ylabel('Signal intensity (J)');
    title(sprintf('Scatter plot - Slice %03d - Reference pix # %03d - Ratio %0.3f', z, ...
        numel(sigpix), base_sig_ratio));
    grid on;
    axis square;

    % the two images with the reference pixels, then save the figure; a slice with
    % no fit is drawn with the previous slice's slope (ROADMAP, question 22)
    plot_reference_overlays(I, J, ref_pix_mask, used_clim, slope);

    sgtitle('Reference pixels regression')
    saveas(gcf, fullfile(plotDir, sprintf('reference_pix_analysis_slice_%03d.png', z)));
    close(gcf);

    % store the pixels and the fit (NaN when there are fewer than two pixels);
    % slice_data does not exist before the first slice
    if ~exist('slice_data', 'var')
        slice_data = struct();
    end
    slice_data(z).basepix = basepix;
    slice_data(z).sigpix = sigpix;
    if numel(basepix) > 1

        % b is [intercept, slope]
        slice_data(z).p = b;
        slice_data(z).delta = delta;
        slice_data(z).ratio = base_sig_ratio;
        slice_data(z).slope = slope;
        slice_data(z).intercept = intercept;
        slice_data(z).used_clim = used_clim;
    else

        % no fit: NaN, which the mean fit of the global correction leaves out
        slice_data(z).p = NaN(2, 1);
        slice_data(z).delta = NaN(size(basepix));
        slice_data(z).ratio = NaN;
        slice_data(z).slope = NaN;
        slice_data(z).intercept = NaN;
        slice_data(z).used_clim = [NaN, NaN];
    end
end

fprintf('Analysis done in %.1f s\n', toc(t0));

end

function [b, slope, intercept, delta] = draw_robust_fit(basepix, sigpix)
% Robust (bisquare) fit of sigpix on basepix, drawn on the current axes with its
% 95% band (half-width delta) and its equation; returns the fit.

% the robust fit; b is [intercept, slope]
[b, stats] = robustfit(basepix, sigpix, 'bisquare');
slope = b(2);
intercept = b(1);
y_fit = slope * basepix + intercept;

% approximate 95% prediction band, from the residuals' standard deviation
resid = sigpix - y_fit;
resid_std = sqrt(sum(resid.^2) / (length(basepix) - 2));
df = length(basepix) - 2;
t_crit = tinv(0.975, df);
mean_x = mean(basepix);

% half-width at each pixel: t times the residual SD, wider away from the mean
delta = t_crit * resid_std * sqrt(1 + 1/length(basepix) + (basepix - mean_x).^2 / ...
    sum((basepix - mean_x).^2));

% sort, for a smooth band
[basepix_sorted, idx] = sort(basepix);
y_fit_sorted = y_fit(idx);
delta_sorted = delta(idx);

% the band
x_patch = [basepix_sorted; flipud(basepix_sorted); basepix_sorted(1)];
y_patch = [y_fit_sorted + delta_sorted; flipud(y_fit_sorted - delta_sorted); ...
    y_fit_sorted(1) + delta_sorted(1)];
patch(x_patch, y_patch, 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
    'DisplayName', '95% CI');

% the fit line
plot(basepix_sorted, y_fit_sorted, 'r-', 'LineWidth', 2, 'DisplayName', 'Fit');

% axes 10% wider than the data, and the equation at the bottom
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

function plot_reference_overlays(I, J, ref_pix_mask, used_clim, slope)
% Panels 2 and 3 of a slice's figure: the base and signal images with the
% reference pixels in red; the base image's limits are used_clim / slope.

% the base image
subplot(1, 3, 2);
imagesc(I);
axis image off;
colormap(sep_palette('anatomy'));
title('Base image (I)');

% the nano's limits divided by the slope, so the autofluorescence shows on the
% scale of the nano it is fitted to
clim(used_clim*1/slope);
hold on;
draw_reference_squares(ref_pix_mask);
hold off;
colorbar;

% the signal image
subplot(1, 3, 3);
imagesc(J);
axis image off;
colormap(sep_palette('anatomy'));
title('Signal image (J)');
clim(used_clim);
hold on;
draw_reference_squares(ref_pix_mask);
hold off;
colorbar;

end

function draw_reference_squares(ref_pix_mask)
% One red square per reference pixel, on the current axes.

[rows, cols] = find(ref_pix_mask);
if ~isempty(rows)

    % the four corners of each pixel, all squares in one patch object
    x = [cols-0.5, cols+0.5, cols+0.5, cols-0.5]';
    y = [rows-0.5, rows-0.5, rows+0.5, rows+0.5]';
    faces = reshape(1:numel(cols)*4, 4, [])';
    patch('Faces', faces, 'Vertices', [x(:), y(:)], ...
        'FaceColor', 'r', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
end

end

% ===== Local functions: summary figure =====

function [average_slope, average_intercept] = plot_regression_summary(slice_data, Z, ...
    plotDir)
% Ratio, slope and intercept of every slice in one figure, saved; returns the mean
% slope and intercept over the slices, which the global correction uses.

figure('Name', 'Jittered Scatter Plots and Curves of Regression Metrics', ...
    'units', 'normalized', 'outerposition', [0 0 1 1]);

% left: the three metrics of the slices with a fit, jittered
subplot(1, 2, 1);
categories = {'Ratio', 'Slope', 'Intercept'};

% one row per slice: ratio, slope, intercept
values = nan(Z, 3);
for z = 1:Z
    if isfield(slice_data(z), 'ratio')
        values(z, 1) = slice_data(z).ratio;
        values(z, 2) = slice_data(z).slope;
        values(z, 3) = slice_data(z).intercept;
    end
end

% only the slices with a fit (two reference pixels or more)
valid_idx = ~isnan(values(:, 1));
values = values(valid_idx, :);
Z_valid = sum(valid_idx);
slice_indices = find(valid_idx);

plot_metric_scatter(values, Z_valid, categories);

% right: the metrics against the slice
subplot(1, 2, 2);
plot_metric_curves(values, slice_indices);

% mean slope and intercept, NaN slices left out
average_slope = nanmean([slice_data.slope]); %#ok<*NANMEAN>
average_intercept = nanmean([slice_data.intercept]);

% save
saveas(gcf, fullfile(plotDir, 'jittered_regression_metrics_with_curves.png'));
close(gcf);

end

function plot_metric_scatter(values, Z_valid, categories)
% Jittered metrics with their median and inter-quartile range, on the current
% axes; draws the jitter with randn, once per call.

% jitter the dots sideways, by the same amount for the three metrics of a slice
jitter = 0.1 * randn(Z_valid, 1);
x_jittered = repmat(1:3, Z_valid, 1) + [jitter, jitter, jitter];

% the slices from black to dark grey, in order
colors = zeros(Z_valid, 3);
for i = 1:3
    colors(:, i) = linspace(0, 0.5, Z_valid)';
end
colors_expanded = repmat(colors, 3, 1);

% scatter
scatter(x_jittered(:), values(:), 50, colors_expanded, 'filled', 'MarkerFaceAlpha', 1);
hold on;

% median and inter-quartile range of each metric, the intercept on the right axis
medians = median(values, 1, 'omitnan');
iqr_vals = iqr(values, 1);
for i = 1:3

    % the intercept on the right axis (0 to 10), the ratio and slope on the left
    % (0 to 2), both in black
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

    % the median, and a band one inter-quartile range wide centred on it
    y = medians(i);
    y_min = y - iqr_vals(i) / 2;
    y_max = y + iqr_vals(i) / 2;

    % the ratio in magenta, the slope in red, the intercept in blue
    if i == 3
        colline = [0, 0, 1];
    elseif i == 2
        colline = [1, 0, 0];
    elseif i == 1
        colline = [1, 0, 1];
    end

    % the median as a bar, the band where it is finite
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

end

function plot_metric_curves(values, slice_indices)
% Ratio, slope and intercept against the slice index, on the current axes.

% the ratio (magenta) and the slope (red) on the left axis, 0 to 2
hold on;
yyaxis left
plot(slice_indices, values(:, 1), '-', 'LineWidth', 1.5, 'DisplayName', 'Ratio', ...
    'Color', [1, 0, 1]);
plot(slice_indices, values(:, 2), '-', 'LineWidth', 1.5, 'DisplayName', 'Slope', ...
    'Color', [1, 0, 0]);
ylabel('Slope / Ratio');
ylim([0, 2])
ax = gca;
ax.YColor = [0 0 0];

% the intercept (blue) on the right axis, 0 to 10
yyaxis right
plot(slice_indices, values(:, 3), '-', 'LineWidth', 1.5, 'DisplayName', 'Intercept', ...
    'Color', [0, 0, 1]);
ylabel('Intercept');
ylim([0, 10])
ax = gca;
ax.YColor = [0 0 0];
xlabel('Slice Index');
title('Regression Metrics Across Slices');
grid on;
axis square;
set(gca, 'fontsize', 12);

end

% ===== Local functions: correction and videos =====

function correct_and_save(selectedVol, selectedVolSig, slice_data, bg_mask_vol, ...
    average_slope, average_intercept, correction_type, use_per_slice, ...
    mouse_name, output_dir, H, W, Z)
% Nano minus the scaled autofluorescence, with each slice's fit (use_per_slice)
% or the mean fit; writes corrected_volume_<type>.mat and
% scaled_auto_volume_<type>.mat in output_dir.

% the corrected nano, the scaled autofluorescence and the nano, filled slice by slice
correctedVol = zeros(H, W, Z, 'single');
scaledautoVol = zeros(H, W, Z, 'single');
nanoVol = zeros(H, W, Z, 'single');

fprintf('Applying %s correction for %s...\n', correction_type, mouse_name);
t0 = tic;
for z = 1:Z

    % the slice of each channel, as single
    I = im2single(selectedVol(:, :, z));
    J = im2single(selectedVolSig(:, :, z));

    % the slice's own fit, or the mean fit over the slices
    if use_per_slice
        slope = slice_data(z).slope;
        intercept = slice_data(z).intercept;
    else
        slope = average_slope;
        intercept = average_intercept;
    end

    % scale the autofluorescence and subtract it, negative values set to 0; a slice
    % with no fit (NaN slope) comes out all NaN in the slicewise correction
    scaled_I = slope * I + intercept;
    corrected = J - scaled_I;
    corrected(corrected < 0) = 0;
    correctedVol(:, :, z) = corrected;
    scaledautoVol(:, :, z) = scaled_I;
    nanoVol(:, :, z) = J;

end
fprintf('Correction done in %.1f s\n', toc(t0));

% save the corrected volume, then the scaled autofluorescence with the nano
matfile_name = fullfile(output_dir, sprintf('corrected_volume_%s.mat', correction_type));
save(matfile_name, 'correctedVol', 'bg_mask_vol', 'slice_data', 'average_slope', ...
    'average_intercept', 'correction_type', '-v7.3');
matfile_name_bis = fullfile(output_dir, sprintf('scaled_auto_volume_%s.mat', ...
    correction_type));
save(matfile_name_bis, 'scaledautoVol', 'nanoVol', 'bg_mask_vol', 'slice_data', ...
    'average_slope', 'average_intercept', 'correction_type', '-v7.3');

end

function write_difference_video(selectedVol, selectedVolSig, slice_data, bg_mask_vol, ...
    average_slope, average_intercept, correction_type, use_per_slice, ...
    output_dir, Z)
% Video of (nano - scaled auto) / scaled auto, one frame per slice, the background
% in black; writes scaled_difference_video_<type>.mp4 in output_dir.

% one frame per slice, a slice a second
videoFile_diff = fullfile(output_dir, sprintf('scaled_difference_video_%s.mp4', ...
    correction_type));
vidObj_diff = VideoWriter(videoFile_diff, 'MPEG-4');
vidObj_diff.FrameRate = 1;
vidObj_diff.Quality = 100;
open(vidObj_diff);

% colour limits of the relative difference, the same for every slice
clim_min_diff = -1;
clim_max_diff = 1;

for z = 1:Z

    % the slice of each channel, and the fit of this correction
    I = im2single(selectedVol(:, :, z));
    J = im2single(selectedVolSig(:, :, z));

    if use_per_slice
        slope = slice_data(z).slope;
        intercept = slice_data(z).intercept;
    else
        slope = average_slope;
        intercept = average_intercept;
    end

    % relative difference, the background left out
    scaled_I = slope * I + intercept;
    diff_map = (J - scaled_I) ./ scaled_I;
    diff_map(bg_mask_vol(:, :, z)) = NaN;

    % draw it on black, blue to red
    figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
    h = imagesc(diff_map);
    clim(gca, [clim_min_diff, clim_max_diff]);
    colorbar;
    colormap(sep_palette('difference'));

    % the background (NaN) transparent over the black axes
    ax = gca;
    ax.Color = 'k';
    alpha_mask = ~isnan(diff_map);
    set(h, 'AlphaData', alpha_mask);
    axis equal;
    grid on;
    ylim([0, size(I, 1)]);
    xlim([0, size(I, 2)]);

    % the slice and the slope used in the title
    if use_per_slice
        title(sprintf( ...
            'Scaled difference map (relative) - Slice # %d (slicewise %0.2f)', z, slope));
    else
        title(sprintf( ...
            'Scaled difference map (relative) - Slice # %d (global %0.2f)', ...
            z, slope)); %#ok<*UNRCH>
    end

    % add the frame
    frame = getframe(gcf);
    writeVideo(vidObj_diff, frame);
    close(gcf);
end
close(vidObj_diff);

end

function write_ratio_video(selectedVol, selectedVolSig, bg_mask_vol, Z, output_dir)
% Video of nano / autofluorescence, one frame per slice, the background in black;
% writes ratio_map_video.mp4 in output_dir.

% one frame per slice, a slice a second
videoFile_ratio = fullfile(output_dir, 'ratio_map_video.mp4');
vidObj_ratio = VideoWriter(videoFile_ratio, 'MPEG-4');
vidObj_ratio.FrameRate = 1;
vidObj_ratio.Quality = 100;
open(vidObj_ratio);

% colour limits of the ratio, the same for every slice
clim_min_ratio = 0;
clim_max_ratio = 2;

for z = 1:Z

    % the slice of each channel, as single
    I = im2single(selectedVol(:, :, z));
    J = im2single(selectedVolSig(:, :, z));

    % the ratio, the background left out
    I_mask = I;
    I_mask(bg_mask_vol(:, :, z)) = NaN;
    J_mask = J;
    J_mask(bg_mask_vol(:, :, z)) = NaN;
    ratio_map = J_mask ./ I_mask;

    % draw it on black in the palette's difference map: blue below 1 (nano
    % dimmer than autofluorescence), white at 1, red above
    figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
    h = imagesc(ratio_map);
    clim(gca, [clim_min_ratio, clim_max_ratio]);
    colorbar;
    colormap(sep_palette('difference'));

    % the background (NaN) transparent over the black axes
    ax = gca;
    ax.Color = 'k';
    alpha_mask = ~isnan(ratio_map);
    set(h, 'AlphaData', alpha_mask);
    axis equal;
    grid on;
    ylim([0, size(I, 1)]);
    xlim([0, size(I, 2)]);
    title(sprintf('Nano / Auto ratio map - Slice # %d', z));

    % add the frame
    frame = getframe(gcf);
    writeVideo(vidObj_ratio, frame);
    close(gcf);
end
close(vidObj_ratio);

end
