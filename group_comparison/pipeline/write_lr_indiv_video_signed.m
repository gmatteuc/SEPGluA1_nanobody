function write_lr_indiv_video_signed(lr_diff_4d, lr_sum_4d, mask_bg_4d, atlas_vol, ...
    save_dir, video_filename, clim_diff, clim_sum, ...
    group_name, mouse_names)
%WRITE_LR_INDIV_VIDEO_SIGNED  Video of every mouse's signed left-right difference.
%   WRITE_LR_INDIV_VIDEO_SIGNED(lr_diff_4d, lr_sum_4d, mask_bg_4d, atlas_vol,
%   save_dir, video_filename, clim_diff, clim_sum, group_name, mouse_names)
%   writes save_dir\video_filename as WRITE_LR_INDIV_VIDEO does, with the
%   signed difference L - R on top, on a red-blue colormap (red for R > L, blue
%   for L > R), and the sum below in hot. A mouse beyond the end of mouse_names
%   is titled by its number.
%
%   Run by group_differences.
%
%   See also WRITE_LR_INDIV_VIDEO.

% create the folder if needed
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% open the video
full_video_path = fullfile(save_dir, video_filename);
vidObj = VideoWriter(full_video_path, 'MPEG-4');
vidObj.FrameRate = 15;
vidObj.Quality = 95;
open(vidObj)

[n_slices, ~, ~, n_mice] = size(lr_diff_4d);

% the atlas over the maps' ML width
n_depth_data = size(lr_diff_4d, 3);
atlas_subset = double(atlas_vol(:, :, 1:n_depth_data));

% red (negative) to blue (positive)
crwb = get_color2color_colormap([1 0 0], [0 0 1]);

fprintf('Writing signed individual video: %s\n', video_filename);

% one figure, cleared after each frame
fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], ...
    'Color', 'k');
set(fh, 'InvertHardcopy', 'off');

for j = 1:n_slices

    % skip the planes with no atlas voxel
    atlas_slice = squeeze(atlas_subset(j, :, :));
    if sum(atlas_slice(:) > 0) == 0 %#ok<LOGSUM>
        if mod(j, 50) == 0
            fprintf('  Skipping slice %d\n', j);
        end
        continue;
    end

    % atlas boundaries, for the overlay (off)
    [gy, gx] = gradient(atlas_slice);
    boundaries = (abs(gx) + abs(gy)) > 0 & (atlas_slice > 0);
    [b_row, b_col] = find(boundaries); %#ok<ASGLU>

    for k = 1:n_mice
        if k > length(mouse_names)
            m_name = sprintf('Mouse %d', k);
        else
            m_name = strrep(mouse_names{k}, '_', ' ');
        end

        % the mouse's plane
        data_diff = squeeze(lr_diff_4d(j, :, :, k));
        data_sum = squeeze(lr_sum_4d(j, :, :, k));
        mask_bg = squeeze(mask_bg_4d(j, :, :, k));

        % shown inside the atlas, outside the mouse's background
        valid_pixels = (atlas_slice > 0) & (~mask_bg);
        alpha_data = double(valid_pixels);

        % top row: the signed left-right difference
        subplot(2, n_mice, k);
        imagesc(data_diff);
        clim(clim_diff);
        colormap(gca, crwb);

        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        axis image;
        axis off;
        set(gca, 'Color', 'k');
        hold on;

        % atlas overlay, off (reason not recorded)
        % plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);

        title(m_name, 'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');

        ylabel('L - R (Signed)', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
        cb = colorbar('Location', 'eastoutside');
        cb.Label.String = 'R > L  |  L > R';
        cb.Color = 'w';
        cb.Label.Color = 'w';

        % bottom row: the left-right sum, an intensity, in hot
        subplot(2, n_mice, k + n_mice);
        imagesc(data_sum);
        clim(clim_sum);
        colormap(gca, hot);

        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        axis image;
        axis off;
        set(gca, 'Color', 'k');
        hold on;

        % atlas overlay, off (reason not recorded)
        % plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);

        ylabel('L + R (Sum)', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
        cb = colorbar('Location', 'eastoutside');
        cb.Label.String = 'Intensity';
        cb.Color = 'w';
        cb.Label.Color = 'w';
    end

    sgtitle(['Slice # ' num2str(j) ' - Directional Asymmetry (' group_name, ')'], ...
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
