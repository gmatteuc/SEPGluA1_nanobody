function write_lr_video_surpmask_rolling(lr_diff_vol, lr_sum_vol, atlas_vol, brain_mask, ...
    save_dir, video_filename, clim_values, ...
    group_name, label_string_diff, label_string_sum, ...
    surp_diff_vol, surp_sum_vol, surp_thresh, slab_range)

% Ensure directory exists
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

full_video_path = fullfile(save_dir, video_filename);
vidObj = VideoWriter(full_video_path, 'MPEG-4');
vidObj.FrameRate = 15;
vidObj.Quality = 95;
open(vidObj);

[n_slices, ~, n_width] = size(lr_diff_vol);

fprintf('Writing rolling slab video: %s (Slab +/- %d)\n', video_filename, slab_range);

% Initialize figure
fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'k');
set(fh, 'InvertHardcopy', 'off');

for j = 1:n_slices

    % Define slab indices
    z_start = max(1, j - slab_range);
    z_end   = min(n_slices, j + slab_range);
    z_indices = z_start:z_end;

    % Optimization: skip empty slices
    if sum(sum(brain_mask(j,:,:))) == 0
        continue;
    end

    % Get slab mask
    mask_slab_3d = logical(brain_mask(z_indices, :, :));
    slab_mask_2d = squeeze(max(mask_slab_3d, [], 1));

    % Get diff data
    raw_diff = lr_diff_vol(z_indices, :, :);
    raw_diff(~mask_slab_3d) = NaN;
    slab_diff = squeeze(nanmedian(raw_diff, 1)); %#ok<*NANMEDIAN>

    % Get sum data
    raw_sum = lr_sum_vol(z_indices, :, :);
    raw_sum(~mask_slab_3d) = NaN;
    slab_sum = squeeze(nanmedian(raw_sum, 1));

    % Get surprise data
    raw_surp_d = surp_diff_vol(z_indices, :, :);
    raw_surp_d(~mask_slab_3d) = NaN;
    slab_surp_diff = squeeze(nanmedian(raw_surp_d, 1));
    raw_surp_s = surp_sum_vol(z_indices, :, :);
    raw_surp_s(~mask_slab_3d) = NaN;
    slab_surp_sum = squeeze(nanmedian(raw_surp_s, 1));

    % Compute alpha masks
    calc_alpha = @(vol) min(1, max(0, vol ./ surp_thresh));
    alpha_diff = calc_alpha(slab_surp_diff);
    alpha_diff(isnan(alpha_diff)) = 0;
    alpha_mask_diff = alpha_diff .* double(slab_mask_2d);
    alpha_sum  = calc_alpha(slab_surp_sum);
    alpha_sum(isnan(alpha_sum)) = 0;
    alpha_mask_sum  = alpha_sum  .* double(slab_mask_2d);

    % Atlas boundaries
    atlasim = squeeze(atlas_vol(j, :, :));
    atlasim = single(atlasim);
    av_warp_boundaries = gradient(atlasim) ~= 0 & (atlasim > 1);
    [row, col] = ind2sub(size(atlasim), find(av_warp_boundaries));

    % Plotting
    subplot(1, 2, 1);
    h1 = imagesc(slab_diff);
    clim(clim_values);
    set(h1, 'AlphaData', alpha_mask_diff);
    if abs(clim_values(1)) == abs(clim_values(2))
        try colormap(gca, get_color2color_colormap([0, 0, 1], [1, 0, 0])); catch, colormap(gca, jet); end
    else
        colormap(gca, 'hot');
    end
    ax1 = gca; ax1.Color = 'k'; axis equal; axis off; hold on;
    line(col, row, 'Marker', '.', 'LineStyle', 'none', 'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);
    xlim([0, n_width]);
    title([group_name ' - ', label_string_diff], 'Color', 'w', 'FontSize', 12);
    cb1 = colorbar; cb1.Color = 'w'; cb1.Label.String = label_string_diff;

    subplot(1, 2, 2);
    h2 = imagesc(slab_sum);
    clim(clim_values);
    set(h2, 'AlphaData', alpha_mask_sum);
    if abs(clim_values(1)) == abs(clim_values(2))
        try colormap(gca, get_color2color_colormap([0, 0, 1], [1, 0, 0])); catch, colormap(gca, jet); end
    else
        colormap(gca, 'hot');
    end
    ax2 = gca; ax2.Color = 'k'; axis equal; axis off; hold on;
    line(col, row, 'Marker', '.', 'LineStyle', 'none', 'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);
    xlim([0, n_width]);

    title([group_name ' - ', label_string_sum], 'Color', 'w', 'FontSize', 12);
    cb2 = colorbar; cb2.Color = 'w'; cb2.Label.String = label_string_sum;

    sgtitle(['Slice # ' num2str(j) ' (Slab \pm' num2str(slab_range) ')'], 'Color', 'w', 'FontSize', 14);

    frame = getframe(fh);
    writeVideo(vidObj, frame);
    clf(fh);

    if mod(j, 50) == 0
        fprintf('  Frame %d written...\n', j);
    end
end

close(vidObj);
close(fh);
fprintf('Rolling slab video saved: %s\n', full_video_path);
end