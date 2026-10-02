function write_lr_indiv_rolling_video(diff_4d, sum_4d, bg_mask_4d, atlas_vol, ...
    save_dir, video_filename, group_name, mouse_names, slab_range, clim_diff, clim_sum)
%WRITE_LR_INDIV_ROLLING_VIDEO  Every mouse's left-right maps, median over a rolling slab.
%   WRITE_LR_INDIV_ROLLING_VIDEO(diff_4d, sum_4d, bg_mask_4d, atlas_vol,
%   save_dir, video_filename, group_name, mouse_names, slab_range, clim_diff,
%   clim_sum) writes save_dir\video_filename as WRITE_LR_INDIV_VIDEO does, but
%   each frame shows, for every mouse, the median of its maps over the planes
%   j - slab_range to j + slab_range, its background left out, with the atlas
%   boundaries of the central plane drawn over them. A voxel is shown if it is
%   in the atlas and tissue in at least one plane of the slab.
%
%   atlas_vol   the annotation, already cropped to the maps' ML width
%
%   Run by group_differences.
%
%   See also WRITE_LR_INDIV_VIDEO.

% open the video, in a folder made if needed
[vidObj, full_video_path] = open_lr_video(save_dir, video_filename);

[n_slices, ~, ~, n_mice] = size(diff_4d);

fprintf('Writing Rolling Video: %s (Slab +/- %d)\n', video_filename, slab_range);

% one figure, cleared after each frame
fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], ...
    'Color', 'k');
set(fh, 'InvertHardcopy', 'off');

for j = 1:n_slices

    % the planes of the slab
    z_start = max(1, j - slab_range);
    z_end = min(n_slices, j + slab_range);
    z_indices = z_start:z_end;

    % skip the planes with no atlas voxel
    atlas_slice = squeeze(atlas_vol(j, :, :));
    if sum(atlas_slice(:) > 0) == 0
        if mod(j, 100) == 0
            fprintf('  Skipping slice %d\n', j);
        end
        continue;
    end

    % atlas boundaries of the central plane
    [gy, gx] = gradient(single(atlas_slice));
    boundaries = (abs(gx) + abs(gy)) > 0 & (atlas_slice > 0);
    [b_row, b_col] = find(boundaries);

    for k = 1:n_mice
        mouse_name = strrep(mouse_names{k}, '_', ' ');

        % the mouse's slab
        raw_slab_d = diff_4d(z_indices, :, :, k);
        raw_slab_s = sum_4d(z_indices, :, :, k);
        mask_slab = bg_mask_4d(z_indices, :, :, k);

        % background to NaN
        raw_slab_d(logical(mask_slab)) = NaN;
        raw_slab_s(logical(mask_slab)) = NaN;

        % median over the slab
        slab_diff = squeeze(nanmedian(raw_slab_d, 1)); %#ok<NANMEDIAN>
        slab_sum = squeeze(nanmedian(raw_slab_s, 1));

        % shown inside the atlas where any plane of the slab is tissue
        slab_bg = squeeze(min(mask_slab, [], 1));
        valid_pixels = (atlas_slice > 0) & (~slab_bg);
        alpha_data = double(valid_pixels);

        % top row: the difference
        subplot(2, n_mice, k);
        imagesc(slab_diff);
        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        clim(clim_diff);
        colormap(gca, sep_palette('intensity'));
        axis image;
        axis off;
        set(gca, 'Color', 'k');
        hold on;
        plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);
        title(mouse_name, 'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');
        ylabel('LR Diff (Slab)', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
        cb = colorbar('Location', 'westoutside');
        cb.Label.String = '|L - R|';
        cb.Color = 'w';

        % bottom row: the sum
        subplot(2, n_mice, k + n_mice);
        imagesc(slab_sum);
        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        clim(clim_sum);
        colormap(gca, sep_palette('intensity'));
        axis image;
        axis off;
        set(gca, 'Color', 'k');
        hold on;
        plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);
        ylabel('LR Sum (Slab)', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
        cb = colorbar('Location', 'westoutside');
        cb.Label.String = 'L + R';
        cb.Color = 'w';
    end

    sgtitle(['Slice ' num2str(j) ' (Slab \pm' num2str(slab_range) ') - ' group_name], ...
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
fprintf('Video saved: %s\n', full_video_path);
end
