function write_lr_video_surpmask_rolling(lr_diff_vol, lr_sum_vol, atlas_vol, ...
    brain_mask, save_dir, video_filename, clim_values, ...
    group_name, label_string_diff, label_string_sum, ...
    surp_diff_vol, surp_sum_vol, surp_thresh, slab_range)
%WRITE_LR_VIDEO_SURPMASK_ROLLING  Surprise-masked left-right video, median over a slab.
%   WRITE_LR_VIDEO_SURPMASK_ROLLING(lr_diff_vol, lr_sum_vol, atlas_vol,
%   brain_mask, save_dir, video_filename, clim_values, group_name,
%   label_string_diff, label_string_sum, surp_diff_vol, surp_sum_vol,
%   surp_thresh, slab_range) writes save_dir\video_filename (MPEG-4, 15 frames
%   per second) as WRITE_LR_VIDEO_SURPMASK does, but each frame shows the
%   median over the planes j - slab_range to j + slab_range, inside brain_mask,
%   of the maps and of their surprise, and each panel takes its opacity from
%   its own surprise volume (surp_diff_vol on the left, surp_sum_vol on the
%   right). The median smooths the maps along AP before the threshold.
%
%   Run by group_differences.
%
%   See also WRITE_LR_VIDEO_SURPMASK.

% create the folder if needed
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% open the video
full_video_path = fullfile(save_dir, video_filename);
vidObj = VideoWriter(full_video_path, 'MPEG-4');
vidObj.FrameRate = 15;
vidObj.Quality = 95;
open(vidObj);

[n_slices, ~, n_width] = size(lr_diff_vol);

fprintf('Writing rolling slab video: %s (Slab +/- %d)\n', video_filename, slab_range);

% one figure, cleared after each frame
fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], ...
    'Color', 'k');
set(fh, 'InvertHardcopy', 'off');

for j = 1:n_slices

    % the planes of the slab
    z_start = max(1, j - slab_range);
    z_end = min(n_slices, j + slab_range);
    z_indices = z_start:z_end;

    % skip the planes with no brain
    if sum(sum(brain_mask(j, :, :))) == 0
        continue;
    end

    % the slab's mask, and its projection: a voxel shown if any plane has brain
    mask_slab_3d = logical(brain_mask(z_indices, :, :));
    slab_mask_2d = squeeze(max(mask_slab_3d, [], 1));

    % median of the difference over the slab
    raw_diff = lr_diff_vol(z_indices, :, :);
    raw_diff(~mask_slab_3d) = NaN;
    slab_diff = squeeze(nanmedian(raw_diff, 1)); %#ok<*NANMEDIAN>

    % median of the sum
    raw_sum = lr_sum_vol(z_indices, :, :);
    raw_sum(~mask_slab_3d) = NaN;
    slab_sum = squeeze(nanmedian(raw_sum, 1));

    % median of the two surprise volumes
    raw_surp_d = surp_diff_vol(z_indices, :, :);
    raw_surp_d(~mask_slab_3d) = NaN;
    slab_surp_diff = squeeze(nanmedian(raw_surp_d, 1));
    raw_surp_s = surp_sum_vol(z_indices, :, :);
    raw_surp_s(~mask_slab_3d) = NaN;
    slab_surp_sum = squeeze(nanmedian(raw_surp_s, 1));

    % opacity: surprise over the threshold, clipped to 0-1, inside the slab's brain
    calc_alpha = @(vol) min(1, max(0, vol ./ surp_thresh));
    alpha_diff = calc_alpha(slab_surp_diff);
    alpha_diff(isnan(alpha_diff)) = 0;
    alpha_mask_diff = alpha_diff .* double(slab_mask_2d);
    alpha_sum = calc_alpha(slab_surp_sum);
    alpha_sum(isnan(alpha_sum)) = 0;
    alpha_mask_sum = alpha_sum .* double(slab_mask_2d);

    % atlas boundaries of the central plane: where the annotation changes along ML
    atlasim = squeeze(atlas_vol(j, :, :));
    atlasim = single(atlasim);
    av_warp_boundaries = gradient(atlasim) ~= 0 & (atlasim > 1);
    [row, col] = ind2sub(size(atlasim), find(av_warp_boundaries));

    % left: the difference; blue-red for symmetric limits (jet if the colormap
    % function is missing)
    subplot(1, 2, 1);
    h1 = imagesc(slab_diff);
    clim(clim_values);
    set(h1, 'AlphaData', alpha_mask_diff);
    if abs(clim_values(1)) == abs(clim_values(2))
        try
            colormap(gca, sep_palette('difference'));
        catch
            colormap(gca, jet);
        end
    else
        colormap(gca, sep_palette('intensity'));
    end
    ax1 = gca;
    ax1.Color = 'k';
    axis equal;
    axis off;
    hold on;
    line(col, row, 'Marker', '.', 'LineStyle', 'none', 'Color', [0.66 0.66 0.66], ...
        'MarkerSize', 0.5);
    xlim([0, n_width]);
    title([group_name ' - ', label_string_diff], 'Color', 'w', 'FontSize', 12);
    cb1 = colorbar;
    cb1.Color = 'w';
    cb1.Label.String = label_string_diff;

    % right: the sum
    subplot(1, 2, 2);
    h2 = imagesc(slab_sum);
    clim(clim_values);
    set(h2, 'AlphaData', alpha_mask_sum);
    if abs(clim_values(1)) == abs(clim_values(2))
        try
            colormap(gca, sep_palette('difference'));
        catch
            colormap(gca, jet);
        end
    else
        colormap(gca, sep_palette('intensity'));
    end
    ax2 = gca;
    ax2.Color = 'k';
    axis equal;
    axis off;
    hold on;
    line(col, row, 'Marker', '.', 'LineStyle', 'none', 'Color', [0.66 0.66 0.66], ...
        'MarkerSize', 0.5);
    xlim([0, n_width]);

    title([group_name ' - ', label_string_sum], 'Color', 'w', 'FontSize', 12);
    cb2 = colorbar;
    cb2.Color = 'w';
    cb2.Label.String = label_string_sum;

    sgtitle(['Slice # ' num2str(j) ' (Slab \pm' num2str(slab_range) ')'], ...
        'Color', 'w', 'FontSize', 14);

    % write the frame, and clear the figure for the next one
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
