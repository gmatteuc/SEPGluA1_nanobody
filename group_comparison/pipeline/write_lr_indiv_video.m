function write_lr_indiv_video(lr_diff_4d, lr_sum_4d, mask_bg_4d, atlas_vol, ...
    save_dir, video_filename, clim_diff, clim_sum, ...
    group_name, mouse_names)
%WRITE_LR_INDIV_VIDEO  Video of every mouse's left-right difference and sum maps.
%   WRITE_LR_INDIV_VIDEO(lr_diff_4d, lr_sum_4d, mask_bg_4d, atlas_vol,
%   save_dir, video_filename, clim_diff, clim_sum, group_name, mouse_names)
%   writes save_dir\video_filename (MPEG-4, 15 frames per second), one frame
%   per AP plane with atlas voxels: one column per mouse, the difference on top
%   (limits clim_diff) and the sum below (limits clim_sum), both in hot, each
%   shown inside the atlas and outside that mouse's background.
%
%   lr_diff_4d, lr_sum_4d   AP x DV x ML x mouse, one hemisphere
%   mask_bg_4d              AP x DV x ML x mouse, true for background
%   atlas_vol               the annotation, cropped here to the maps' ML width
%   mouse_names             the column titles, in the order of the fourth dimension
%
%   Run by group_differences.

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

[n_slices, ~, ~, n_mice] = size(lr_diff_4d);

% the atlas over the maps' ML width
n_depth_data = size(lr_diff_4d, 3);
atlas_subset = double(atlas_vol(:, :, 1:n_depth_data));

fprintf('Writing individual video: %s\n', video_filename);

% one figure, cleared after each frame
fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], ...
    'Color', 'k');
set(fh, 'InvertHardcopy', 'off');

for j = 1:n_slices

    % skip the planes with no atlas voxel
    atlas_slice = squeeze(atlas_subset(j, :, :));
    if sum(atlas_slice(:) > 0) == 0
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
        mouse_name = strrep(mouse_names{k}, '_', ' ');

        % the mouse's plane
        data_diff = squeeze(lr_diff_4d(j, :, :, k));
        data_sum = squeeze(lr_sum_4d(j, :, :, k));
        mask_bg = squeeze(mask_bg_4d(j, :, :, k));

        % shown inside the atlas, outside the mouse's background
        valid_pixels = (atlas_slice > 0) & (~mask_bg);
        alpha_data = double(valid_pixels);

        % top row: the left-right difference
        subplot(2, n_mice, k);
        imagesc(data_diff);
        clim(clim_diff);
        colormap(gca, hot);
        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        axis image;
        axis off;
        set(gca, 'Color', 'k');
        hold on;

        % atlas overlay, off (reason not recorded)
        % plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);

        title(mouse_name, 'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');

        ylabel('LR Diff', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
        cb = colorbar('Location', 'westoutside');
        cb.Label.String = '|L - R|';
        cb.Color = 'w';
        cb.Label.Color = 'w';

        % bottom row: the left-right sum
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

        ylabel('LR Sum', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
        cb = colorbar('Location', 'westoutside');
        cb.Label.String = 'L + R';
        cb.Color = 'w';
        cb.Label.Color = 'w';
    end

    sgtitle(['Slice # ' num2str(j) ' - Left Hemishpere asymmetry (' group_name, ')'], ...
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
