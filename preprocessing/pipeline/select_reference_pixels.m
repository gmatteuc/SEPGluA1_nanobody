function [ref_pix_mask, range_pix, bg_mask, bg_mask_dilated, used_clim, h_diag] = ...
    select_reference_pixels(I, p_min, p_max, disk_px, range_frac, plot_flag)
%SELECT_REFERENCE_PIXELS  Reference pixels and background of a slice.
%   [ref_pix_mask, range_pix, bg_mask, bg_mask_dilated, used_clim, h_diag] =
%       SELECT_REFERENCE_PIXELS(I, p_min, p_max, disk_px, range_frac, plot_flag)
%
%   The percentiles 1 to 100 of the slice I (its most common value left
%   out) are smoothed with a 5-point Gaussian. The knee
%   where the tissue starts is the peak of their second derivative, between
%   percentiles p_min and p_max, where the curve rises fastest (p_max when no
%   peak is found). Pixels below the knee are the background (bg_mask),
%   dilated by a disk of disk_px pixels (bg_mask_dilated). The reference
%   pixels (ref_pix_mask) lie outside the dilated background, from the value
%   at the dip of the second derivative paired with the knee to range_frac of
%   the way up to the upper end, the first percentile past the knee where the
%   second derivative rises above the knee's peak again (100 if none);
%   range_pix is their upper bound.
%
%   Defaults: p_min 20, p_max 65, disk_px 15, range_frac 0.20, plot_flag
%   false. used_clim (display limits, the knee to 1.5 times the upper end)
%   and h_diag (the diagnostic figure) are assigned only when plot_flag is
%   true, so a caller that asks for them otherwise stops.

if nargin < 2 || isempty(p_min)
    p_min = 20;
end
if nargin < 3 || isempty(p_max)
    p_max = 65;
end
if nargin < 4 || isempty(disk_px)
    disk_px = 15;
end
if nargin < 5 || isempty(range_frac)
    range_frac = 0.20;
end
if nargin < 6 || isempty(plot_flag)
    plot_flag = false;
end

% get the input as single
I_single = im2single(I);

% percentile curve, without the most common value
p = 1:100;
pix_vals = I_single(:);
pix_vals(pix_vals==mode(pix_vals)) = [];
vals = prctile(pix_vals, p);

% smooth it (5-point Gaussian) and find the peaks and dips of its second
% derivative between p_min and p_max
w = 5;
x = -(w-1)/2 : (w-1)/2;
g = exp(-0.5*(x/(0.3*(w-1))).^2);
g = g / sum(g);
vals_smooth = conv(vals, g, 'same');
d2 = [0 diff(diff(vals_smooth)) 0];
win = (p >= p_min) & (p <= p_max);
[~, locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 10);
[~, locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 10);

% the knee: the peak where the curve rises fastest, and the dip of the same rank;
% found again with prominence 5 when that dip comes first, 5 then 2 when none is found
if not(isempty(locs_max)) && not(isempty(locs_min))

    % the knee and its dip
    [idx_max_bis, idx_min] = knee_and_dip(vals_smooth, locs_max, locs_min);
    if idx_max_bis>idx_min
        [~, locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 5);
        [~, locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 5);

        % the knee and its dip again
        [idx_max_bis, idx_min] = knee_and_dip(vals_smooth, locs_max, locs_min);
    end
else
    [~, locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 5);
    [~, locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 5);
    if not(not(isempty(locs_max)) && not(isempty(locs_min)))
        [~, locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 2);
        [~, locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 2);
    end

    % the knee and its dip
    [idx_max_bis, idx_min] = knee_and_dip(vals_smooth, locs_max, locs_min);
end

% no knee found (every findpeaks above came back empty): the knee falls back to
% p_max, as in select_background_pixels, rather than stopping below
if isempty(idx_max_bis)
    idx_max_bis = min(p_max, numel(vals));
    warning('select_reference_pixels:noKnee', ...
        'No knee detected in percentile curve; falling back to p_max=%d as threshold.', ...
        idx_max_bis);
end
if isempty(idx_min)
    idx_min = max(idx_max_bis - 5, 1);
end

% the upper end: past the knee's peak of the second derivative (its highest
% value within 5 percentiles), the first percentile where it rises above it again
if idx_max_bis<=90
    half_window = 5;
    start_idx = max(1, idx_max_bis - half_window);
    end_idx = min(length(d2), idx_max_bis + half_window);
    [local_max_val, local_rel_idx] = max(d2(start_idx:end_idx));
    local_max_idx = start_idx + local_rel_idx - 1;
    first_d2_bump_idx = local_max_idx;
    first_d2_bump_val = local_max_val;
else
    first_d2_bump_idx = idx_max_bis+1;
    first_d2_bump_val = d2(idx_max_bis+1);
end
d2_masked = d2;
d2_masked = d2_masked.*((1:100)>first_d2_bump_idx);
idx_max = find(d2_masked>first_d2_bump_val, 1, 'first');
if isempty(idx_max)
    idx_max = 100;
end

% intensities at the upper end, the knee and the dip
val_max = vals(idx_max);
val_max_bis = vals(idx_max_bis);
val_min = vals(idx_min);

% off: a plot of the curve and its derivatives, for checking by eye
% d1 = [0 diff(vals) 0];
% figure; plot(d2); hold on; plot(vals); plot(d1);

% reference pixels: from the dip, range_frac of the way up to the upper end
vals_range_start = val_min;
vals_range_end = range_frac*(val_max-vals_range_start)+vals_range_start;
range_pix = vals_range_end;
ref_pix_range = [vals_range_start, vals_range_end];
ref_pix_mask = and(I_single>ref_pix_range(1), I_single<ref_pix_range(2));

% background below the knee, dilated; the reference pixels kept off it
m = I_single < val_max_bis;
bg_mask = logical(m);
se = strel('disk', disk_px);
bg_mask_dilated = imdilate(bg_mask, se);
ref_pix_mask = and(ref_pix_mask, not(bg_mask_dilated));

% diagnostic figure, on request
if plot_flag
    h_diag = figure('name', 'Background mask diagnostics', 'units', 'normalized', ...
        'outerposition', [0 0 1 1]);
    t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    % the slice with the reference pixels in red
    nexttile;
    imagesc(I_single);
    axis image off;
    colormap(sep_palette('anatomy'));
    title('Input slice with ref pixels');
    used_clim = [val_max_bis, 1.5*val_max];
    clim(used_clim)
    hold on;
    [rows, cols] = find(ref_pix_mask);
    if ~isempty(rows)

        % one square patch per pixel
        x = [cols-0.5, cols+0.5, cols+0.5, cols-0.5]';
        y = [rows-0.5, rows-0.5, rows+0.5, rows+0.5]';
        faces = reshape(1:numel(cols)*4, 4, [])';
        patch('Faces', faces, 'Vertices', [x(:), y(:)], ...
            'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
    hold off;

    % the background mask
    nexttile;
    imagesc(bg_mask);
    axis image off;
    colormap(sep_palette('anatomy'));
    title('Background mask');

    % the percentile curve with the knee (blue), the dip (red), the upper end
    % (magenta) and the range of the reference pixels
    nexttile;
    plot(p, vals, 'LineWidth', 1.5);
    grid on;
    hold on;
    plot(idx_max_bis, val_max_bis, 'o', 'MarkerFaceColor', [0, 0, 1], ...
        'MarkerEdgeColor', [0, 0, 1]);
    plot(idx_min, val_min, 'o', 'MarkerFaceColor', [1, 0, 0], ...
        'MarkerEdgeColor', [1, 0, 0]);
    plot(idx_max, val_max, 'o', 'MarkerFaceColor', [1, 0, 1], ...
        'MarkerEdgeColor', [1, 0, 1]);
    plot([0, 100], [ref_pix_range(1), ref_pix_range(1)], '--', 'Color', [1, 0, 0])
    plot([0, 100], [ref_pix_range(2), ref_pix_range(2)], '--', 'Color', [1, 0, 0])
    xlabel('Percentile');
    ylabel('intensity');
    axis square
    title(sprintf('Percentiles (knee at p=%d, thr=%.3g)', idx_max_bis, val_max_bis));
    hold off;
    title(t, 'Reference pixels estimation diagnostics');
end

end

% ===== Local functions =====

function [idx_max_bis, idx_min] = knee_and_dip(vals_smooth, locs_max, locs_min)
% The knee, the peak of the second derivative where the smoothed curve rises
% fastest, and the dip of the same rank (5 percentiles below the knee if none).

% off: the first peak and dip instead (reason not recorded)
% idx_max_bis=locs_max(1);
% idx_min=locs_min(1);
d1 = diff(vals_smooth);
[~, chosen_max_idx] = max(d1(locs_max));
idx_max_bis = locs_max(chosen_max_idx);
if numel(locs_min)>=chosen_max_idx
    idx_min = locs_min(chosen_max_idx);
else
    idx_min = max(idx_max_bis-5, 1);
end

end
