function bg_mask = select_background_pixels(I, p_min, p_max, plot_flag)
%SELECT_BACKGROUND_PIXELS  Mask of the background pixels of one slice image.
%   bg_mask = SELECT_BACKGROUND_PIXELS(I, p_min, p_max, plot_flag) returns a
%   logical mask, the size of I, of the pixels below a threshold found at the
%   knee of the image's percentile curve.
%
%   The pixel values above the most frequent one (usually the zero padding)
%   are taken at the percentiles 1 to 100. On that curve, lightly smoothed,
%   the knee is where background turns into tissue: a peak of the second
%   derivative between the percentiles p_min and p_max, the one with the
%   steepest first derivative within 7 percentiles of it. The peaks and
%   troughs are found at a prominence of 10; at 5, then 2, when there is no
%   peak or no trough, and at 5 when the knee lies above its nearest trough.
%   With no peak at all, the threshold falls back to the percentile p_max,
%   with a warning. A slice that is empty (at most 10% of its pixels above
%   the most frequent value), or whose centre of mass lies outside the middle
%   third of the image, is all background.
%
%   Inputs:
%     I          one slice image
%     p_min      lower end of the percentile window of the knee (default 20)
%     p_max      upper end of that window (default 65)
%     plot_flag  draw a diagnostic figure (default false)

if nargin < 2 || isempty(p_min)
    p_min = 20;
end
if nargin < 3 || isempty(p_max)
    p_max = 65;
end

% diagnostics off by default: run_nano_equalisation calls this with three
% arguments, at four call sites; run_normalise_groups passes the flag
if nargin < 4 || isempty(plot_flag)
    plot_flag = false;
end

% the image as single
I_single = im2single(I);

% the pixel values for the percentile curve, without the most frequent value
% (usually 0, the padding)
p = 1:100;
pix_vals = I_single(:);
pix_vals(pix_vals <= mode(pix_vals)) = [];

% an empty or nearly empty image
bool_empty = or(isempty(pix_vals), numel(pix_vals) <= 0.1 * numel(I(:)));

% a centre of mass implausibly close to the border
[rows, cols] = size(I_single);
[x_grid, y_grid] = meshgrid(1:cols, 1:rows);
total_mass = sum(I_single(:));
if total_mass == 0
    x_center = NaN;
    y_center = NaN;
else
    x_center = sum(sum(x_grid .* I_single)) / total_mass;
    y_center = sum(sum(y_grid .* I_single)) / total_mass;
end
x_out_bool = or(x_center < 0.33 * size(I_single, 2), x_center > 0.66 * size(I_single, 2));
y_out_bool = or(y_center < 0.33 * size(I_single, 1), y_center > 0.66 * size(I_single, 1));
bool_out = or(x_out_bool, y_out_bool);

% all background if the image is empty or off centre
if or(bool_empty, bool_out)
    bg_mask = true(size(I));
    return;
end

% the percentile curve and its knee
[vals, d2, d1, idx_max_bis] = find_knee(pix_vals, p, p_min, p_max);

% no knee at any prominence (it happens on faint channels such as the
% autofluorescence): the percentile p_max is the threshold
if isempty(idx_max_bis)
    idx_max_bis = min(p_max, length(vals));
    warning('select_background_pixels:noKnee', ...
        'No knee detected in percentile curve; falling back to p_max=%d as threshold.', ...
        idx_max_bis);
end

% the threshold at the knee, and the background mask (not dilated)
val_max_bis = vals(idx_max_bis);
bg_mask = I_single < val_max_bis;

% diagnostic figure: the slice, the mask, and the percentile curve with its knee
if plot_flag
    plot_background_diagnostics(I_single, val_max_bis, bg_mask, p, vals, idx_max_bis, ...
        d1, d2);
end

end

% ===== Local functions =====

function [vals, d2, d1, idx_max_bis] = find_knee(pix_vals, p, p_min, p_max)
% The percentile curve, smoothed, its derivatives, and the knee of its second
% derivative in the window (empty if there is none at any prominence).

vals = prctile(pix_vals, p);

% smooth the percentile curve (5-point Gaussian), then its second derivative
w = 5;
x = -(w - 1) / 2 : (w - 1) / 2;
g = exp(-0.5 * (x / (0.3 * (w - 1))).^2);
g = g / sum(g);
vals_smooth = conv(vals, g, 'same');
d2 = [0 diff(diff(vals_smooth)) 0];

% its first derivative, for the knee search
d1 = [0, diff(vals_smooth)];

% peaks and troughs of the second derivative in the percentile window
win = (p >= p_min) & (p <= p_max);
[~, locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 10);
[~, locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 10);

if not(isempty(locs_max)) && not(isempty(locs_min))

    % the knee and its nearest trough, the derivative searched within the window
    [idx_max_bis, idx_min] = knee_and_trough(d1, locs_max, locs_min, p_min + 1, p_max);

    % a knee above its trough is inconsistent: redo with a lower prominence
    if idx_max_bis > idx_min
        [~, locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 5);
        [~, locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 5);

        % the knee again, the derivative searched over every percentile
        idx_max_bis = knee_and_trough(d1, locs_max, locs_min, 1, length(d1));
    end
else

    % no clear peak or trough: redo with a lower prominence, 5, then 2
    [~, locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 5);
    [~, locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 5);
    if not(not(isempty(locs_max)) && not(isempty(locs_min)))
        [~, locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 2);
        [~, locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 2);
    end

    % the knee, the derivative searched over every percentile (its trough is not used)
    idx_max_bis = knee_and_trough(d1, locs_max, locs_min, 1, length(d1));
end

end

function plot_background_diagnostics(I_single, val_max_bis, bg_mask, p, vals, ...
    idx_max_bis, d1, d2)
% The diagnostic figure: the slice, the mask, and the percentile curve with its knee.

val_max = max(I_single(:));
h_diag = figure('name', 'Background mask diagnostics', ...
    'units', 'normalized', 'outerposition', [-0.05 -0.05 0.9 0.9], ...
    'Color', 'w'); %#ok<NASGU>
t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
imagesc(I_single);
axis image off;
colormap(gca, sep_palette('anatomy'));
title('Input slice with ref pixels', 'FontSize', 12);
used_clim = [val_max_bis, val_max];
if used_clim(2) > used_clim(1)
    clim(used_clim);
end
nexttile;
imagesc(bg_mask);
axis image off;
colormap(gca, sep_palette('anatomy'));
title('Background mask', 'FontSize', 12);
nexttile;
yyaxis left

% the percentile curve in blue
h_int = plot(p, vals, '-', 'LineWidth', 2, 'Color', [0 0 0.8]);
hold on;
h_knee = plot(idx_max_bis, val_max_bis, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0 0 0.8], 'MarkerEdgeColor', 'w');
ylabel('Intensity', 'FontSize', 11);
set(gca, 'YColor', [0 0 0.8]);
ylim([min(vals(:)), max(vals(:))*1.05]);
grid on;
yyaxis right

% its first and second derivatives in purple
h_d1 = plot(p, d1, '-', 'LineWidth', 1.5, 'Color', [0.8 0 0.8]);
hold on;
h_d2 = plot(p, d2, '--', 'LineWidth', 1.5, 'Color', [0.8 0 0.8]);
ylabel('Derivatives (1st & 2nd)', 'FontSize', 11);
set(gca, 'YColor', [0.8 0 0.8]);
xlabel('Percentile', 'FontSize', 11);
xlim([0 100]);
axis square;
title(sprintf('Knee Detection (p=%d, thr=%.3g)', idx_max_bis, val_max_bis), ...
    'FontSize', 11, 'FontWeight', 'normal');
legend([h_int, h_d1, h_d2, h_knee], ...
    {'Intensity', '1st Deriv', '2nd Deriv', 'Knee Point'}, ...
    'Location', 'best', 'FontSize', 9);
hold off;
title(t, 'Reference Pixels Estimation Diagnostics', 'FontSize', 14, ...
    'FontWeight', 'bold');

end

function [idx_max_bis, idx_min] = knee_and_trough(d1, locs_max, locs_min, idx_low, ...
    idx_high)
% The knee, the peak whose neighbourhood (7 percentiles each way, kept within
% idx_low to idx_high) has the steepest first derivative d1, and its nearest trough.

% off: the first derivative at the peak itself (reason not recorded)
% [~, chosen_max_idx] = max(max_vals_in_range);
% d1 = [0, diff(vals_smooth)];
% [~,chosen_max_idx] = max(d1(locs_max));
% idx_max_bis = locs_max(chosen_max_idx);
tol_range = 7;
max_vals_in_range = zeros(size(locs_max));
for k = 1:length(locs_max)
    idx_start = max(idx_low, locs_max(k) - tol_range);
    idx_end   = min(idx_high, locs_max(k) + tol_range);
    max_vals_in_range(k) = max(d1(idx_start:idx_end));
end
[~, chosen_max_idx] = max(max_vals_in_range);

% that peak is the knee, which sets the background threshold
idx_max_bis = locs_max(chosen_max_idx);

% the trough nearest to it
% off: the trough of the same rank, or 5 percentiles below the knee (reason
% not recorded)
% if numel(locs_min) >= chosen_max_idx
%     idx_min = locs_min(chosen_max_idx);
% else
%     idx_min = max(idx_max_bis-5, 1);
% end
[~, idx_idx_min] = min(abs(locs_min - idx_max_bis));
idx_min = locs_min(idx_idx_min);

end
