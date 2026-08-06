function [ref_pix_mask, range_pix, bg_mask, bg_mask_dilated, used_clim, h_diag] = select_reference_pixels(I, p_min, p_max, disk_px, range_frac, plot_flag)

if nargin < 2 || isempty(p_min), p_min = 20; end
if nargin < 3 || isempty(p_max), p_max = 65; end
if nargin < 4 || isempty(disk_px), disk_px = 15; end
if nargin < 5 || isempty(plot_flag), plot_flag = true; end

% get imput and cast to single
I_single = im2single(I);

% compute percentile curve
p = 1:100;
pix_vals=I_single(:);
pix_vals(pix_vals==mode(pix_vals))=[];
vals = prctile(pix_vals, p);

% detect smoothed max and min of d2
w = 5;
x = -(w-1)/2 : (w-1)/2;
g = exp(-0.5*(x/(0.3*(w-1))).^2);
g = g / sum(g);
vals_smooth = conv(vals, g, 'same');
d2 = [0 diff(diff(vals_smooth)) 0];
win = (p >= p_min) & (p <= p_max);
[~,locs_max] = findpeaks(d2 .* win, 'MinPeakProminence',10);
[~,locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence',10);
if not(isempty(locs_max)) && not(isempty(locs_min))
    % idx_max_bis=locs_max(1);
    % idx_min=locs_min(1);
    d1=diff(vals_smooth);
    [~,choosen_max_idx] = max(d1(locs_max));
    idx_max_bis=locs_max(choosen_max_idx);
    if numel(locs_min)>=choosen_max_idx
        idx_min=locs_min(choosen_max_idx);
    else
        idx_min=max(idx_max_bis-5,1);
    end
    if idx_max_bis>idx_min
        [~,locs_max] = findpeaks(d2 .* win, 'MinPeakProminence',5);
        [~,locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence',5);
        % idx_max_bis=locs_max(1);
        % idx_min=locs_min(1);
        d1=diff(vals_smooth);
        [~,choosen_max_idx] = max(d1(locs_max));
        idx_max_bis=locs_max(choosen_max_idx);
        if numel(locs_min)>=choosen_max_idx
            idx_min=locs_min(choosen_max_idx);
        else
            idx_min=max(idx_max_bis-5,1);
        end
    end
else
    [~,locs_max] = findpeaks(d2 .* win, 'MinPeakProminence',5);
    [~,locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence',5);
    if not(not(isempty(locs_max)) && not(isempty(locs_min)))
        [~,locs_max] = findpeaks(d2 .* win, 'MinPeakProminence',2);
        [~,locs_min] = findpeaks(-(d2 .* win), 'MinPeakProminence',2);
    end
    % idx_max_bis=locs_max(1);
    % idx_min=locs_min(1);
    d1=diff(vals_smooth);
    [~,choosen_max_idx] = max(d1(locs_max));
    idx_max_bis=locs_max(choosen_max_idx);
    if numel(locs_min)>=choosen_max_idx
        idx_min=locs_min(choosen_max_idx);
    else
        idx_min=max(idx_max_bis-5,1);
    end
end
% detect d2 final divergence as threshold crossing of d2 above original max bump value
if idx_max_bis<=90
    win = 5;
    start_idx = max(1, idx_max_bis - win);
    end_idx = min(length(d2), idx_max_bis + win);
    [local_max_val, local_rel_idx] = max(d2(start_idx:end_idx));
    local_max_idx = start_idx + local_rel_idx - 1;
    first_d2_bump_idx = local_max_idx;
    first_d2_bump_val = local_max_val;
else
    first_d2_bump_idx=idx_max_bis+1;
    first_d2_bump_val=d2(idx_max_bis+1);
end
d2_masked=d2;
d2_masked=d2_masked.*((1:100)>first_d2_bump_idx);
idx_max=find(d2_masked>first_d2_bump_val,1,'first');
if isempty(idx_max)
    idx_max = 100;
end

% select only last min and last two max
val_max = vals(idx_max);
val_max_bis = vals(idx_max_bis);
val_min = vals(idx_min);

% d1 = [0 diff(vals) 0];
% figure; plot(d2); hold on; plot(vals); plot(d1); 

% select needed pixels
vals_range_start = val_min;
vals_range_end = range_frac*(val_max-vals_range_start)+vals_range_start;
range_pix = vals_range_end;
ref_pix_range = [vals_range_start,vals_range_end];
ref_pix_mask = and(I_single>ref_pix_range(1),I_single<ref_pix_range(2));

% get mask from first max
m = I_single < val_max_bis;
bg_mask = logical(m);
se = strel('disk', disk_px);
bg_mask_dilated = imdilate(bg_mask,se);
ref_pix_mask = and(ref_pix_mask,not(bg_mask_dilated));

% optional diagnostics
if plot_flag
    h_diag = figure('name','Background mask diagnostics','units','normalized','outerposition',[0 0 1 1]);
    t = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
    % input with ref_pix_mask overlay
    nexttile;
    imagesc(I_single); axis image off; colormap gray; title('Input slice with ref pixels');
    used_clim=[val_max_bis,1.5*val_max];
    clim(used_clim)
    hold on;
    [rows, cols] = find(ref_pix_mask);
    if ~isempty(rows)
        % Create patch vertices (expand each pixel to a small square for shading)
        x = [cols-0.5, cols+0.5, cols+0.5, cols-0.5]';
        y = [rows-0.5, rows-0.5, rows+0.5, rows+0.5]';
        faces = reshape(1:numel(cols)*4, 4, [])';
        patch('Faces', faces, 'Vertices', [x(:), y(:)], ...
            'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
    hold off;
    % mask
    nexttile;
    imagesc(bg_mask); axis image off; colormap(gray); title('Background mask');
    % percentile curve + knee
    nexttile;
    plot(p, vals, 'LineWidth', 1.5); grid on; hold on;
    plot(idx_max_bis, val_max_bis, 'o', 'MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1]);
    plot(idx_min, val_min, 'o', 'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0]);
    plot(idx_max, val_max, 'o', 'MarkerFaceColor',[1,0,1],'MarkerEdgeColor',[1,0,1]);
    plot([0,100],[ref_pix_range(1),ref_pix_range(1)],'--','Color',[1,0,0])
    plot([0,100],[ref_pix_range(2),ref_pix_range(2)],'--','Color',[1,0,0])
    xlabel('Percentile'); ylabel('intensity');
    axis square
    title(sprintf('Percentiles (knee at p=%d, thr=%.3g)', idx_max_bis, val_max_bis));
    hold off;
    title(t, 'Reference pixels estimation diagnostics');
end

end