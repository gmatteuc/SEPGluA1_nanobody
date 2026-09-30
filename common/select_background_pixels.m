function bg_mask = select_background_pixels(I, p_min, p_max, plot_flag)

if nargin < 2 || isempty(p_min), p_min = 20; end
if nargin < 3 || isempty(p_max), p_max = 65; end
% plot_flag had no default, so calling this with three arguments (which is what
% P2bis does, at four different call sites) left it undefined and the function
% errored at "if plot_flag" near the end. Diagnostics off by default; P6bis
% passes the flag explicitly and is unaffected.
if nargin < 4 || isempty(plot_flag), plot_flag = false; end

% get input and cast to single
I_single = im2single(I);

% compute percentile curve
p = 1:100;
pix_vals = I_single(:);
pix_vals(pix_vals<=mode(pix_vals)) = []; % Remove most frequent value (usually 0/padding)

% Check if image is empty or quasi-empty
bool_empty = or(isempty(pix_vals),numel(pix_vals)<=0.1*numel(I(:)));
% Check if center of mass is unrealistically close to border
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
x_out_bool = or(x_center<0.33*size(I_single,2),x_center>0.66*size(I_single,2));
y_out_bool = or(y_center<0.33*size(I_single,1),y_center>0.66*size(I_single,1));
bool_out = or(x_out_bool,y_out_bool);
% Return a full background mask if imege is empty or out
if or(bool_empty,bool_out)
    bg_mask = true(size(I));
    return;
end
% ------------------------------

vals = prctile(pix_vals, p);

% detect smoothed max and min of d2
w = 5;
x = -(w-1)/2 : (w-1)/2;
g = exp(-0.5*(x/(0.3*(w-1))).^2);
g = g / sum(g);
vals_smooth = conv(vals, g, 'same');
d2 = [0 diff(diff(vals_smooth)) 0];

% Find minima and maxima of second derivative in acceptability window
win = (p >= p_min) & (p <= p_max);
[~,locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 10);
[~,locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 10);

if not(isempty(locs_max)) && not(isempty(locs_min))

    % Find max of first derivative in the points detected as peaks of second
    % [~, choosen_max_idx] = max(max_vals_in_range);
    % d1 = [0, diff(vals_smooth)];
    % [~,choosen_max_idx] = max(d1(locs_max));
    % idx_max_bis = locs_max(choosen_max_idx);
    d1 = [0, diff(vals_smooth)];
    tol_range = 7;
    max_vals_in_range = zeros(size(locs_max));
    for k = 1:length(locs_max)
        idx_start = max(p_min+1, locs_max(k) - tol_range);
        idx_end   = min(p_max, locs_max(k) + tol_range);
        max_vals_in_range(k) = max(d1(idx_start:idx_end));
    end
    [~, choosen_max_idx] = max(max_vals_in_range);
    % Pick that max as maximum (that will dtermine background)
    idx_max_bis = locs_max(choosen_max_idx);
    % Pick close minimum
    % if numel(locs_min) >= choosen_max_idx
    %     idx_min = locs_min(choosen_max_idx);
    % else
    %     idx_min = max(idx_max_bis-5, 1);
    % end
    [~, idx_idx_min] = min(abs(locs_min-idx_max_bis));
    idx_min = locs_min(idx_idx_min);
    % Redo with lower prominence threshold if result inconsitent
    if idx_max_bis > idx_min
        [~,locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 5);
        [~,locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 5);
        % [~, choosen_max_idx] = max(max_vals_in_range);
        % d1 = [0, diff(vals_smooth)];
        % [~,choosen_max_idx] = max(d1(locs_max));
        % idx_max_bis = locs_max(choosen_max_idx);
        d1 = [0, diff(vals_smooth)];
        tol_range = 7;
        max_vals_in_range = zeros(size(locs_max));
        for k = 1:length(locs_max)
            idx_start = max(1, locs_max(k) - tol_range);
            idx_end   = min(length(d1), locs_max(k) + tol_range);
            max_vals_in_range(k) = max(d1(idx_start:idx_end));
        end
        [~, choosen_max_idx] = max(max_vals_in_range);
        idx_max_bis = locs_max(choosen_max_idx);
        % if numel(locs_min) >= choosen_max_idx
        %     idx_min = locs_min(choosen_max_idx);
        % else
        %     idx_min = max(idx_max_bis-5, 1);
        % end
        [~, idx_idx_min] = min(abs(locs_min-idx_max_bis));
        idx_min = locs_min(idx_idx_min);
    end
else
    % Redo with lower prominence threshold if result inconsitent
    [~,locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 5);
    [~,locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 5);
    if not(not(isempty(locs_max)) && not(isempty(locs_min)))
        [~,locs_max] = findpeaks(d2 .* win, 'MinPeakProminence', 2);
        [~,locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence', 2);
    end
    % [~, choosen_max_idx] = max(max_vals_in_range);
    % d1 = [0, diff(vals_smooth)];
    % [~,choosen_max_idx] = max(d1(locs_max));
    % idx_max_bis = locs_max(choosen_max_idx);
    d1 = [0, diff(vals_smooth)];
    tol_range = 7;
    max_vals_in_range = zeros(size(locs_max));
    for k = 1:length(locs_max)
        idx_start = max(1, locs_max(k) - tol_range);
        idx_end   = min(length(d1), locs_max(k) + tol_range);
        max_vals_in_range(k) = max(d1(idx_start:idx_end));
    end
    [~, choosen_max_idx] = max(max_vals_in_range);
    idx_max_bis = locs_max(choosen_max_idx);
    % if numel(locs_min) >= choosen_max_idx
    %     idx_min = locs_min(choosen_max_idx);
    % else
    %     idx_min = max(idx_max_bis-5, 1);
    % end
    [~, idx_idx_min] = min(abs(locs_min-idx_max_bis));
    idx_min = locs_min(idx_idx_min);
end

% Fallback: if no knee was detectable (all findpeaks attempts returned empty),
% use p_max as the threshold percentile. This can happen on faint channels
% (e.g. autofluo) where the smoothed-percentile curve has no clear knee.
if isempty(idx_max_bis)
    idx_max_bis = min(p_max, length(vals));
    warning('select_background_pixels:noKnee', ...
        'No knee detected in percentile curve; falling back to p_max=%d as threshold.', ...
        idx_max_bis);
end

% Extract threshold value from the found index
val_max_bis = vals(idx_max_bis);
% Create Background Mask (Undilated)
bg_mask = I_single < val_max_bis;

% Optional diagnostics
if plot_flag
    val_max = max(I_single(:));
    h_diag = figure('name','Background mask diagnostics',...
        'units','normalized','outerposition',[-0.05 -0.05 0.9 0.9], ...
        'Color', 'w'); %#ok<NASGU>
    t = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
    nexttile;
    imagesc(I_single);
    axis image off;
    colormap(gca, gray);
    title('Input slice with ref pixels', 'FontSize', 12);
    used_clim = [val_max_bis, val_max];
    if used_clim(2) > used_clim(1)
        clim(used_clim);
    end
    nexttile;
    imagesc(bg_mask);
    axis image off;
    colormap(gca, gray);
    title('Background mask', 'FontSize', 12);
    nexttile;
    yyaxis left
    h_int = plot(p, vals, '-', 'LineWidth', 2, 'Color', [0 0 0.8]); % Blue
    hold on;
    h_knee = plot(idx_max_bis, val_max_bis, 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', [0 0 0.8], 'MarkerEdgeColor', 'w');
    ylabel('Intensity', 'FontSize', 11);
    set(gca, 'YColor', [0 0 0.8]);
    ylim([min(vals(:)), max(vals(:))*1.05]);
    grid on;
    yyaxis right
    h_d1 = plot(p, d1, '-', 'LineWidth', 1.5, 'Color', [0.8 0 0.8]); % Gold
    hold on;
    h_d2 = plot(p, d2, '--', 'LineWidth', 1.5, 'Color', [0.8 0 0.8]); % Purple
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
    title(t, 'Reference Pixels Estimation Diagnostics', 'FontSize', 14, 'FontWeight', 'bold');
end

end